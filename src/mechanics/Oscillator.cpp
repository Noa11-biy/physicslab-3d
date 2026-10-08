#include "physicslab/mechanics/Oscillator.hpp"

#include <cmath>
#include <limits>

#include "physicslab/core/Constants.hpp"

namespace pl {

double OscillatorProblem::omega0() const { return std::sqrt(stiffness / mass); }
double OscillatorProblem::period() const { return 2.0 * constants::pi / omega0(); }
double OscillatorProblem::gamma() const { return damping / (2.0 * mass); }
double OscillatorProblem::zeta() const { return damping / (2.0 * std::sqrt(stiffness * mass)); }

DampingRegime OscillatorProblem::regime() const {
    const double z = zeta();
    if (std::abs(z - 1.0) < 1e-9) return DampingRegime::Critical;
    return z < 1.0 ? DampingRegime::Underdamped : DampingRegime::Overdamped;
}

double OscillatorProblem::dampedOmega() const {
    const double d = omega0() * omega0() - gamma() * gamma();
    return d > 0.0 ? std::sqrt(d) : 0.0;
}

double OscillatorProblem::steadyStateAmplitude(double w) const {
    const double a = forceAmplitude / mass;
    const double w0 = omega0(), g = gamma();
    const double d0 = w0 * w0 - w * w, d1 = 2.0 * g * w;
    const double D = d0 * d0 + d1 * d1;
    return D > 0.0 ? a / std::sqrt(D) : std::numeric_limits<double>::infinity();
}

double OscillatorProblem::resonanceOmega() const {
    const double d = omega0() * omega0() - 2.0 * gamma() * gamma();
    return d > 0.0 ? std::sqrt(d) : 0.0;
}

OdeFunction OscillatorProblem::rhs() const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const double m = mass, k = stiffness, c = damping, f0 = forceAmplitude, w = forceFrequency;
    return [m, k, c, f0, w](double t, const State& y, State& dydt) {
        dydt[0] = y[1];
        dydt[1] = (-k * y[0] - c * y[1] + f0 * std::cos(w * t)) / m;
    };
}

void OscillatorProblem::exact(double t, double& x, double& v) const {
    const double w0 = omega0(), g = gamma(), w = forceFrequency;

    // --- solution particulière (régime forcé) et sa valeur initiale ---
    double xp = 0.0, vp = 0.0, xp0 = 0.0, vp0 = 0.0;
    if (forceAmplitude != 0.0) {
        const double a = forceAmplitude / mass;
        const double d0 = w0 * w0 - w * w, d1 = 2.0 * g * w;
        const double D = d0 * d0 + d1 * d1;
        if (D < 1e-26 * w0 * w0 * w0 * w0) {  // résonance exacte sans frottement : croissance séculaire
            xp = a * t * std::sin(w0 * t) / (2.0 * w0);
            vp = a * (std::sin(w0 * t) + w0 * t * std::cos(w0 * t)) / (2.0 * w0);
        } else {
            const double A = a * d0 / D, B = a * d1 / D;
            xp = A * std::cos(w * t) + B * std::sin(w * t);
            vp = w * (-A * std::sin(w * t) + B * std::cos(w * t));
            xp0 = A;
            vp0 = w * B;
        }
    }

    // --- solution homogène avec les conditions initiales corrigées ---
    const double xh0 = x0 - xp0, vh0 = v0 - vp0;
    double xh, vh;
    switch (regime()) {
        case DampingRegime::Underdamped: {
            const double wd = dampedOmega();
            const double e = std::exp(-g * t), c = std::cos(wd * t), s = std::sin(wd * t);
            const double q = (vh0 + g * xh0) / wd;
            xh = e * (xh0 * c + q * s);
            vh = e * (vh0 * c - ((g * vh0 + w0 * w0 * xh0) / wd) * s);
            break;
        }
        case DampingRegime::Critical: {
            const double e = std::exp(-g * t), b = vh0 + g * xh0;
            xh = e * (xh0 + b * t);
            vh = e * (vh0 - g * b * t);
            break;
        }
        default: {  // Overdamped
            const double s = std::sqrt(g * g - w0 * w0);
            const double r1 = -g + s, r2 = -g - s;
            const double c1 = (vh0 - r2 * xh0) / (r1 - r2), c2 = (r1 * xh0 - vh0) / (r1 - r2);
            const double e1 = std::exp(r1 * t), e2 = std::exp(r2 * t);
            xh = c1 * e1 + c2 * e2;
            vh = r1 * c1 * e1 + r2 * c2 * e2;
            break;
        }
    }

    x = xh + xp;
    v = vh + vp;
}

double oscillatorError(const OscillatorProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);
    double xe, ve;
    problem.exact(t, xe, ve);
    const double dx = y[0] - xe, dv = (y[1] - ve) / problem.omega0();
    return std::sqrt(dx * dx + dv * dv);
}

}  // namespace pl
