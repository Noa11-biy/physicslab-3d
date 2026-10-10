#include "physicslab/waves/String.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "physicslab/core/Constants.hpp"

namespace pl::waves {

using constants::pi;

// ------------------------------------------------------- StringProblem -----

double StringProblem::speed() const { return std::sqrt(tension / density); }
double StringProblem::spacing() const { return length / (beads + 1); }
double StringProblem::mass() const { return density * spacing(); }
double StringProblem::stiffness() const { return tension / spacing(); }
double StringProblem::fundamental() const { return speed() / (2.0 * length); }
double StringProblem::continuumOmega(int n) const { return n * pi * speed() / length; }
double StringProblem::maxOmega() const { return 2.0 * std::sqrt(stiffness() / mass()); }
double StringProblem::verletDtLimit() const { return 2.0 / maxOmega(); }
double StringProblem::rk4DtLimit() const { return 2.0 * std::sqrt(2.0) / maxOmega(); }

double StringProblem::chainOmega(int n) const {
    return maxOmega() * std::sin(n * pi / (2.0 * (beads + 1)));
}

std::vector<double> StringProblem::pluck(double x0, double height) const {
    // x0 est tenu à l'intérieur de la corde : un pincement exactement sur un bout n'a pas de sens (triangle de largeur nulle)
    x0 = std::clamp(x0, 1e-9 * length, (1.0 - 1e-9) * length);
    std::vector<double> u(beads);
    for (int j = 1; j <= beads; ++j) {
        const double x = j * spacing();
        u[j - 1] = x < x0 ? height * x / x0 : height * (length - x) / (length - x0);
    }
    return u;
}

std::vector<double> StringProblem::modeShape(int n, double amplitude) const {
    std::vector<double> u(beads);
    for (int j = 1; j <= beads; ++j) u[j - 1] = amplitude * std::sin(n * pi * j / (beads + 1));
    return u;
}

State StringProblem::state(const std::vector<double>& u, const std::vector<double>& v) {
    assert(u.size() == v.size());
    State y(u);
    y.insert(y.end(), v.begin(), v.end());
    return y;
}

OdeFunction StringProblem::rhs() const {
    const int n = beads;
    const double ratio = stiffness() / mass();  // k / m = c² / a²
    return [n, ratio](double, const State& y, State& dydt) {
        for (int j = 0; j < n; ++j) {
            dydt[j] = y[n + j];
            const double left = j > 0 ? y[j - 1] : 0.0;        // u_0 = 0 : bout fixe
            const double right = j < n - 1 ? y[j + 1] : 0.0;   // u_{N+1} = 0
            dydt[n + j] = ratio * (left - 2.0 * y[j] + right);
        }
    };
}

double StringProblem::kineticEnergy(const State& y) const {
    double sum = 0.0;
    for (int j = 0; j < beads; ++j) sum += y[beads + j] * y[beads + j];
    return 0.5 * mass() * sum;
}

double StringProblem::potentialEnergy(const State& y) const {
    double sum = 0.0;
    for (int j = 0; j <= beads; ++j) {
        const double left = j > 0 ? y[j - 1] : 0.0, right = j < beads ? y[j] : 0.0;  // ressort entre la masse j et la masse j + 1
        sum += (right - left) * (right - left);
    }
    return 0.5 * stiffness() * sum;
}

double StringProblem::pluckCoefficient(int n, double x0, double height) const {
    const double l = length;
    return 2.0 * height * l * l / (n * n * pi * pi * x0 * (l - x0)) * std::sin(n * pi * x0 / l);
}

double StringProblem::continuumPluck(double x, double t, double x0, double height, int terms) const {
    double sum = 0.0;
    for (int n = 1; n <= terms; ++n)
        sum += pluckCoefficient(n, x0, height) * std::sin(n * pi * x / length) * std::cos(continuumOmega(n) * t);
    return sum;
}

// --------------------------------------------------------- StringModes -----

StringModes::StringModes(const StringProblem& problem) : p_(problem), n_(problem.beads) {
    omega_.resize(n_);
    sine_.resize(static_cast<std::size_t>(n_) * n_);
    for (int n = 1; n <= n_; ++n) {
        omega_[n - 1] = p_.chainOmega(n);
        for (int j = 1; j <= n_; ++j) sine_[static_cast<std::size_t>(n - 1) * n_ + (j - 1)] = std::sin(n * pi * j / (n_ + 1));
    }
}

std::vector<double> StringModes::amplitudes(const std::vector<double>& f) const {
    assert(static_cast<int>(f.size()) == n_);
    std::vector<double> q(n_, 0.0);
    const double scale = 2.0 / (n_ + 1);
    for (int n = 0; n < n_; ++n) {
        const double* row = &sine_[static_cast<std::size_t>(n) * n_];
        double sum = 0.0;
        for (int j = 0; j < n_; ++j) sum += f[j] * row[j];
        q[n] = scale * sum;
    }
    return q;
}

State StringModes::exact(const std::vector<double>& u0, const std::vector<double>& v0, double t) const {
    const std::vector<double> a = amplitudes(u0), b = amplitudes(v0);
    std::vector<double> q(n_), qdot(n_);
    for (int n = 0; n < n_; ++n) {
        const double w = omega_[n], c = std::cos(w * t), s = std::sin(w * t);
        q[n] = a[n] * c + b[n] / w * s;
        qdot[n] = -a[n] * w * s + b[n] * c;
    }
    State y(2 * static_cast<std::size_t>(n_), 0.0);
    for (int n = 0; n < n_; ++n) {
        const double* row = &sine_[static_cast<std::size_t>(n) * n_];
        for (int j = 0; j < n_; ++j) {
            y[j] += q[n] * row[j];
            y[n_ + j] += qdot[n] * row[j];
        }
    }
    return y;
}

std::vector<double> StringModes::modeEnergies(const State& y) const {
    const std::vector<double> u(y.begin(), y.begin() + n_), v(y.begin() + n_, y.end());
    const std::vector<double> q = amplitudes(u), p = amplitudes(v);
    std::vector<double> e(n_);
    const double weight = p_.mass() * (n_ + 1) / 4.0;
    for (int n = 0; n < n_; ++n) e[n] = weight * (p[n] * p[n] + omega_[n] * omega_[n] * q[n] * q[n]);
    return e;
}

// ----------------------------------------------------------- erreur --------

double stringError(const StringProblem& problem, const std::vector<double>& u0, Solver& solver, int steps, double tEnd) {
    const std::vector<double> v0(problem.beads, 0.0);
    const OdeFunction f = problem.rhs();
    State y = StringProblem::state(u0, v0);
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) {
        advance(solver, f, t, y, dt);
        t += dt;
    }
    const State exact = StringModes(problem).exact(u0, v0, tEnd);
    double worst = 0.0;
    for (int j = 0; j < problem.beads; ++j) {
        const double e = std::abs(y[j] - exact[j]);
        if (!(e <= worst)) worst = e;  // un NaN se propage
    }
    return worst;
}

}  // namespace pl::waves
