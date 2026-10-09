#include "physicslab/mechanics/Pendulum.hpp"

#include <algorithm>
#include <cmath>

namespace pl {
namespace {

constexpr double kPi = constants::pi;

// Moyenne arithmético-géométrique de a et b (a, b > 0).
double agm(double a, double b) {
    for (int i = 0; i < 60 && std::abs(a - b) > 1e-16 * a; ++i) {
        const double next = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = next;
    }
    return 0.5 * (a + b);
}

}  // namespace

double ellipticK(double k) {
    const double kp = std::sqrt(std::max(0.0, 1.0 - k * k));  // module complémentaire
    return kPi / (2.0 * agm(1.0, kp));
}

double PendulumProblem::omega0() const { return std::sqrt(gravity / length); }
double PendulumProblem::smallAnglePeriod() const { return 2.0 * kPi / omega0(); }
double PendulumProblem::gamma() const { return damping / (mass * length * length); }

bool PendulumProblem::hasExactSolution() const {
    return damping == 0.0 && angularVelocity0 == 0.0 && std::abs(theta0) < kPi - 1e-3;
}

double PendulumProblem::period() const {
    if (!hasExactSolution()) return smallAnglePeriod();
    const double a = std::abs(theta0);
    if (a < 1e-9) return smallAnglePeriod();
    return 4.0 * ellipticK(std::sin(0.5 * a)) / omega0();
}

double PendulumProblem::energy(double theta, double omega) const {
    return 0.5 * mass * length * length * omega * omega + mass * gravity * length * (1.0 - std::cos(theta));
}

OdeFunction PendulumProblem::rhs() const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const double w02 = gravity / length, g = gamma();
    return [w02, g](double, const State& y, State& dydt) {
        dydt[0] = y[1];
        dydt[1] = -w02 * std::sin(y[0]) - g * y[1];
    };
}

void PendulumProblem::exact(double t, double& theta, double& omega) const {
    const double sign = theta0 < 0.0 ? -1.0 : 1.0;
    const double a = std::abs(theta0);
    const double w0 = omega0();

    if (a < 1e-6) {  // le terme en theta^3 est sous l'arrondi : oscillateur harmonique
        theta = theta0 * std::cos(w0 * t);
        omega = -theta0 * w0 * std::sin(w0 * t);
        return;
    }

    const double k = std::sin(0.5 * a), kp = std::cos(0.5 * a);  // module et module complémentaire
    const double K = kPi / (2.0 * agm(1.0, kp));
    const double Kp = kPi / (2.0 * agm(1.0, k));
    const double q = std::exp(-kPi * Kp / K);  // nome
    const double x = kPi * w0 * t / (2.0 * K);

    // s(t) = k sn(K + w0 t) et sa dérivée par rapport à x.
    double s = 0.0, ds = 0.0;
    double qn = std::sqrt(q);  // q^{n + 1/2}
    for (int n = 0; n < 2000 && qn > 1e-18; ++n) {
        const double m = 2.0 * n + 1.0;
        const double coef = (n % 2 == 0 ? 1.0 : -1.0) * qn / (1.0 - qn * qn);
        s += coef * std::cos(m * x);
        ds += -coef * m * std::sin(m * x);
        qn *= q;
    }
    const double scale = 2.0 * kPi / K;
    s *= scale;
    ds *= scale * kPi * w0 / (2.0 * K);  // d/dt = (pi w0 / 2K) d/dx

    theta = sign * 2.0 * std::asin(std::clamp(s, -1.0, 1.0));
    omega = sign * 2.0 * ds / std::sqrt(std::max(1.0 - s * s, 1e-300));
}

void PendulumProblem::reference(double t, double& theta, double& omega) const {
    if (hasExactSolution()) {
        exact(t, theta, omega);
        return;
    }
    RK45 solver;
    solver.relTol = 1e-13;
    solver.absTol = 1e-15;
    State y = initialState();
    advance(solver, rhs(), 0.0, y, t);
    theta = y[0];
    omega = y[1];
}

double pendulumError(const PendulumProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);

    double thetaRef, omegaRef;
    problem.reference(t, thetaRef, omegaRef);
    const double dTheta = y[0] - thetaRef, dOmega = (y[1] - omegaRef) / problem.omega0();
    return std::sqrt(dTheta * dTheta + dOmega * dOmega);
}

}  // namespace pl
