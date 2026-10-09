#include "physicslab/mechanics/DoublePendulum.hpp"

#include <cmath>
#include <limits>

namespace pl {

OdeFunction DoublePendulumProblem::rhs() const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const double m1 = this->m1, m2 = this->m2, l1 = this->l1, l2 = this->l2, g = gravity;
    return [m1, m2, l1, l2, g](double, const State& y, State& d) {
        const double th1 = y[0], th2 = y[1], w1 = y[2], w2 = y[3];
        const double delta = th1 - th2;
        const double s = std::sin(delta), c = std::cos(delta);
        const double den = 2.0 * m1 + m2 - m2 * std::cos(2.0 * delta);

        d[0] = w1;
        d[1] = w2;
        d[2] = (-g * (2.0 * m1 + m2) * std::sin(th1) - m2 * g * std::sin(th1 - 2.0 * th2) -
                2.0 * s * m2 * (w2 * w2 * l2 + w1 * w1 * l1 * c)) /
               (l1 * den);
        d[3] = (2.0 * s * (w1 * w1 * l1 * (m1 + m2) + g * (m1 + m2) * std::cos(th1) + w2 * w2 * l2 * m2 * c)) /
               (l2 * den);
    };
}

double DoublePendulumProblem::energy(const State& y) const {
    const double th1 = y[0], th2 = y[1], w1 = y[2], w2 = y[3];
    const double kinetic = 0.5 * (m1 + m2) * l1 * l1 * w1 * w1 + 0.5 * m2 * l2 * l2 * w2 * w2 +
                           m2 * l1 * l2 * w1 * w2 * std::cos(th1 - th2);
    const double potential = -(m1 + m2) * gravity * l1 * std::cos(th1) - m2 * gravity * l2 * std::cos(th2);
    return kinetic + potential;
}

void DoublePendulumProblem::positions(const State& y, double& x1, double& y1, double& x2, double& y2) const {
    x1 = l1 * std::sin(y[0]);
    y1 = -l1 * std::cos(y[0]);
    x2 = x1 + l2 * std::sin(y[1]);
    y2 = y1 - l2 * std::cos(y[1]);
}

State DoublePendulumProblem::reference(double t) const {
    RK45 solver;
    solver.relTol = 1e-13;
    solver.absTol = 1e-15;
    State y = initialState();
    advance(solver, rhs(), 0.0, y, t);
    return y;
}

double DoublePendulumProblem::distance(const State& a, const State& b) const {
    const double wScale = std::sqrt(gravity / l1);
    const double d0 = a[0] - b[0], d1 = a[1] - b[1], d2 = (a[2] - b[2]) / wScale, d3 = (a[3] - b[3]) / wScale;
    return std::sqrt(d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3);
}

double doublePendulumError(const DoublePendulumProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);
    return problem.distance(y, problem.reference(t));
}

double lyapunovExponent(const std::vector<double>& t, const std::vector<double>& d, double dMin, double dMax,
                        int* used, double* intercept) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (std::size_t i = 0; i < t.size() && i < d.size(); ++i) {
        if (!(d[i] >= dMin && d[i] <= dMax)) continue;  // écarte aussi les NaN
        const double y = std::log(d[i]);
        sx += t[i]; sy += y; sxx += t[i] * t[i]; sxy += t[i] * y;
        ++n;
    }
    if (used) *used = n;
    if (n < 8) return std::numeric_limits<double>::quiet_NaN();
    const double dn = static_cast<double>(n);
    const double denom = dn * sxx - sx * sx;
    if (!(std::abs(denom) > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double slope = (dn * sxy - sx * sy) / denom;
    if (intercept) *intercept = (sy - slope * sx) / dn;
    return slope;
}

}  // namespace pl
