#include "physicslab/mechanics/Friction.hpp"

#include <cmath>
#include <limits>

namespace pl {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// phi(tau) = (1 - e^{-k tau}) / k  et  psi(tau) = (tau - phi) / k : les deux fonctions de la solution d'une phase
// (limites tau et tau^2/2 quand k -> 0). expm1 et un développement limité évitent la perte de chiffres pour k tau petit.
double phi(double k, double tau) {
    const double x = k * tau;
    return x == 0.0 ? tau : tau * (-std::expm1(-x) / x);
}

double psi(double k, double tau) {
    const double x = k * tau;
    if (x < 1e-2) return tau * tau * (0.5 - x / 6.0 + x * x / 24.0 - x * x * x / 120.0 + x * x * x * x / 720.0);
    return (tau - phi(k, tau)) / k;
}

int sgn(double v) { return (v > 0.0) - (v < 0.0); }

}  // namespace

bool InclineProblem::holds() const { return std::tan(angle) <= muStatic; }

double InclineProblem::stageAcceleration(int dir) const {
    return -gravity * (std::sin(angle) + dir * muKinetic * std::cos(angle));
}

OdeFunction InclineProblem::rhs(InclineModel model) const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const double gs = gravity * std::sin(angle), gc = gravity * std::cos(angle), mu = muKinetic, k = drag, eps = regularization;
    const bool smooth = model == InclineModel::Regularized;
    return [gs, gc, mu, k, eps, smooth](double, const State& y, State& d) {
        const double sign = smooth ? std::tanh(y[1] / eps) : static_cast<double>(sgn(y[1]));
        d[0] = y[1];
        d[1] = -gs - mu * gc * sign - k * y[1];
    };
}

OdeFunction InclineProblem::rhsStage(int dir) const {
    const double A = stageAcceleration(dir), k = drag;
    return [A, k](double, const State& y, State& d) {
        d[0] = y[1];
        d[1] = A - k * y[1];
    };
}

double InclineProblem::energy(double s, double v) const { return 0.5 * v * v + gravity * s * std::sin(angle); }

namespace {

// Durée avant que la vitesse v ne s'annule dans une phase d'accélération A (l'infini si elle ne s'annule pas).
double stopDuration(double v, double A, double k) {
    if (!(v * A < 0.0)) return kInf;
    return k == 0.0 ? -v / A : std::log1p(-k * v / A) / k;
}

}  // namespace

double InclineProblem::firstStopTime() const {
    if (v0 == 0.0) return holds() ? 0.0 : kInf;  // déjà arrêté et collé, ou glisse vers le bas sans jamais s'arrêter
    return stopDuration(v0, stageAcceleration(v0 > 0.0 ? 1 : -1), drag);
}

InclineState InclineProblem::exact(double t) const {
    InclineState st{s0, v0, 0.0, false};
    if (t <= 0.0) return st;

    int dir = sgn(v0);
    if (dir == 0) {
        if (holds()) {
            st.stuck = true;
            return st;
        }
        dir = -1;  // départ au repos sur une pente trop raide : glissement vers le bas
    }

    double elapsed = 0.0;
    for (int stage = 0; stage < 3; ++stage) {
        const double A = stageAcceleration(dir);
        const double tStop = stopDuration(st.v, A, drag);
        const bool ends = t - elapsed > tStop;
        const double tau = ends ? tStop : t - elapsed;

        const double sNew = st.s + st.v * phi(drag, tau) + A * psi(drag, tau);
        st.path += std::abs(sNew - st.s);  // la vitesse ne change pas de signe dans une phase
        st.s = sNew;
        st.v = ends ? 0.0 : st.v * std::exp(-drag * tau) + A * phi(drag, tau);
        if (!ends) return st;

        elapsed += tStop;
        if (holds()) {  // arrêt : adhérence définitive
            st.stuck = true;
            return st;
        }
        dir = -1;       // sinon repart vers le bas
    }
    return st;
}

InclineRun::InclineRun(const InclineProblem& problem) : problem_(problem) {
    state_ = {problem.s0, problem.v0, 0.0, false};
    dir_ = sgn(problem.v0);
    if (dir_ == 0) {
        if (problem.holds()) state_.stuck = true;
        dir_ = -1;
    }
}

void InclineRun::advance(Solver& solver, double dt) {
    double remaining = dt;
    for (int guard = 0; guard < 8 && remaining > 1e-15 * dt; ++guard) {
        if (state_.stuck) break;

        // Une phase : EDO linéaire de sens fixé ; l'événement est le passage de dir * v par zéro (arrêt).
        const OdeFunction f = problem_.rhsStage(dir_);
        const int dir = dir_;
        const EventFunction crossing = [dir](double, const State& y) { return dir * y[1]; };
        State y{state_.s, state_.v};
        const EventStep r = advanceToEvent(solver, f, crossing, time_, y, remaining);

        state_.path += std::abs(y[0] - state_.s);
        state_.s = y[0];
        state_.v = y[1];
        time_ += r.elapsed;
        remaining -= r.elapsed;

        if (r.event) {
            state_.v = 0.0;
            if (stops_++ == 0) firstStop_ = time_;
            if (problem_.holds()) state_.stuck = true;  // adhérence définitive
            else dir_ = -1;                              // repart vers le bas
        }
    }
    if (state_.stuck) time_ += remaining > 0.0 ? remaining : 0.0;  // le temps passe, rien ne bouge
}

double inclineError(const InclineProblem& problem, InclineModel model, Solver& solver, int steps, double tEnd) {
    const double dt = tEnd / steps;
    double s = 0.0;
    if (model == InclineModel::EventDriven) {
        InclineRun run(problem);
        for (int i = 0; i < steps; ++i) run.advance(solver, dt);
        s = run.state().s;
    } else {
        const OdeFunction f = problem.rhs(model);
        State y = problem.initialState();
        double t = 0.0;
        for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);
        s = y[0];
    }
    return std::abs(s - problem.exact(tEnd).s);
}

}  // namespace pl
