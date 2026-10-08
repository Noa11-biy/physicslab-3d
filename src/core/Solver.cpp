#include "physicslab/core/Solver.hpp"

#include <algorithm>
#include <cmath>

namespace pl {

// --- Euler explicite ---------------------------------------------------------------

double ExplicitEuler::step(const OdeFunction& f, double t, State& y, double dt) {
    k_.resize(y.size());
    f(t, y, k_);
    for (std::size_t i = 0; i < y.size(); ++i) y[i] += dt * k_[i];
    return dt;
}

// --- Euler symplectique ------------------------------------------------------------

double SymplecticEuler::step(const OdeFunction& f, double t, State& y, double dt) {
    const std::size_t n = y.size() / 2;
    k_.resize(y.size());
    f(t, y, k_);
    for (std::size_t i = 0; i < n; ++i) y[n + i] += dt * k_[n + i];  // v_{n+1} = v_n + dt a_n
    for (std::size_t i = 0; i < n; ++i) y[i] += dt * y[n + i];       // x_{n+1} = x_n + dt v_{n+1}
    return dt;
}

// --- Verlet des vitesses -----------------------------------------------------------

double VelocityVerlet::step(const OdeFunction& f, double t, State& y, double dt) {
    const std::size_t n = y.size() / 2;
    k1_.resize(y.size());
    k2_.resize(y.size());
    tmp_.resize(y.size());

    f(t, y, k1_);  // k1_[n + i] = a_n
    for (std::size_t i = 0; i < n; ++i) {
        tmp_[i] = y[i] + dt * y[n + i] + 0.5 * dt * dt * k1_[n + i];  // x_{n+1}
        tmp_[n + i] = y[n + i] + dt * k1_[n + i];                    // prédicteur de v_{n+1}
    }
    f(t + dt, tmp_, k2_);  // k2_[n + i] = a_{n+1}
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = tmp_[i];
        y[n + i] += 0.5 * dt * (k1_[n + i] + k2_[n + i]);
    }
    return dt;
}

// --- RK4 ---------------------------------------------------------------------------

double RK4::step(const OdeFunction& f, double t, State& y, double dt) {
    const std::size_t n = y.size();
    k1_.resize(n); k2_.resize(n); k3_.resize(n); k4_.resize(n); tmp_.resize(n);

    f(t, y, k1_);
    for (std::size_t i = 0; i < n; ++i) tmp_[i] = y[i] + 0.5 * dt * k1_[i];
    f(t + 0.5 * dt, tmp_, k2_);
    for (std::size_t i = 0; i < n; ++i) tmp_[i] = y[i] + 0.5 * dt * k2_[i];
    f(t + 0.5 * dt, tmp_, k3_);
    for (std::size_t i = 0; i < n; ++i) tmp_[i] = y[i] + dt * k3_[i];
    f(t + dt, tmp_, k4_);
    for (std::size_t i = 0; i < n; ++i) y[i] += dt / 6.0 * (k1_[i] + 2.0 * k2_[i] + 2.0 * k3_[i] + k4_[i]);
    return dt;
}

// --- RK45 (Dormand-Prince) ---------------------------------------------------------

namespace {

// Tableau de Butcher de Dormand-Prince. La ligne 6 (7e étage) est identique aux poids de l'ordre 5.
constexpr double kC[7] = {0.0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0};
constexpr double kA[7][6] = {
    {},
    {1.0 / 5},
    {3.0 / 40, 9.0 / 40},
    {44.0 / 45, -56.0 / 15, 32.0 / 9},
    {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
    {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
    {35.0 / 384, 0.0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84},
};
// Coefficients de l'estimation d'erreur : (poids ordre 5) - (poids ordre 4).
constexpr double kE[7] = {71.0 / 57600, 0.0, -71.0 / 16695, 71.0 / 1920, -17253.0 / 339200, 22.0 / 525, -1.0 / 40};

}  // namespace

double RK45::step(const OdeFunction& f, double t, State& y, double dt) {
    const std::size_t n = y.size();
    for (State& k : k_) k.resize(n);
    tmp_.resize(n);

    const bool callerLimits = hNext_ > 0.0 && dt < hNext_;  // le plafond de l'appelant coupe le pas
    double h = hNext_ > 0.0 ? std::min(dt, hNext_) : dt;
    const double hMin = 1e-14 * std::max(1.0, std::abs(t));

    for (int attempt = 0;; ++attempt) {
        f(t, y, k_[0]);
        for (int s = 1; s < 7; ++s) {
            for (std::size_t i = 0; i < n; ++i) {
                double acc = 0.0;
                for (int j = 0; j < s; ++j) acc += kA[s][j] * k_[j][i];
                tmp_[i] = y[i] + h * acc;
            }
            f(t + kC[s] * h, tmp_, k_[s]);  // à s = 6, tmp_ contient la solution d'ordre 5
        }
        evaluations_ += 7;

        // Erreur locale estimée, norme RMS pondérée par la tolérance.
        double sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double e = 0.0;
            for (int j = 0; j < 7; ++j) e += kE[j] * k_[j][i];
            e *= h;
            const double scale = absTol + relTol * std::max(std::abs(y[i]), std::abs(tmp_[i]));
            sum += (e / scale) * (e / scale);
        }
        const double err = n > 0 ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
        const double factor = err > 0.0 ? std::clamp(0.9 * std::pow(err, -0.2), 0.2, 5.0) : 5.0;

        if (err <= 1.0 || h <= hMin || attempt >= 100) {
            y = tmp_;
            ++accepted_;
            // Un pas raccourci par l'appelant ne doit pas faire rétrécir le pas courant.
            hNext_ = callerLimits ? std::max(hNext_, h * factor) : h * factor;
            return h;
        }
        ++rejected_;
        h *= factor;
    }
}

double advance(Solver& solver, const OdeFunction& f, double t, State& y, double dt) {
    double elapsed = 0.0;
    while (dt - elapsed > 1e-12 * dt) {
        const double h = solver.step(f, t + elapsed, y, dt - elapsed);
        if (!(h > 0.0)) break;  // garde-fou : un solveur défaillant ne doit pas boucler
        elapsed += h;
    }
    return elapsed;
}

}  // namespace pl
