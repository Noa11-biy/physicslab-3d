#include "physicslab/mechanics/Kepler.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pl {

double kepler::solveEccentricAnomaly(double meanAnomaly, double e, int* iterations) {
    if (iterations) *iterations = 0;
    if (e <= 0.0) return meanAnomaly;  // cercle : E = M

    const double twoPi = 2.0 * constants::pi;
    // Réduction à M dans [0, pi] par périodicité (2 pi) et symétrie (M -> -M donne E -> -E).
    const double turns = std::round(meanAnomaly / twoPi);
    double m = meanAnomaly - turns * twoPi;
    const double sign = m < 0.0 ? -1.0 : 1.0;
    m = std::abs(m);
    if (m == 0.0) return turns * twoPi;  // périastre : E = 0 (la racine est au bord du crochet)

    // f(E) = E - e sin E - M est strictement croissante (f' = 1 - e cos E > 0) et change de signe sur [M, min(M + e, pi)] :
    // ce crochet permet de ne jamais diverger, même pour e proche de 1 où f' est presque nulle près du périastre.
    double lo = m, hi = std::min(m + e, constants::pi);
    // Départ de Danby : proche de la racine pour toute excentricité.
    double E = std::clamp(m + 0.85 * e, lo, hi);

    int n = 0;
    for (; n < 50; ++n) {
        const double f = E - e * std::sin(E) - m;
        if (f == 0.0) break;
        (f < 0.0 ? lo : hi) = E;
        double next = E - f / (1.0 - e * std::cos(E));  // Newton
        if (!(next >= lo && next <= hi)) next = 0.5 * (lo + hi);  // sort du crochet : dichotomie
        const double step = std::abs(next - E);
        E = next;
        if (step <= 4.0 * 2.220446049250313e-16 * std::max(1.0, E)) { ++n; break; }
    }
    if (iterations) *iterations = n;
    return sign * E + turns * twoPi;
}

double KeplerProblem::period() const { return 2.0 * constants::pi * std::sqrt(a * a * a / mu); }
double KeplerProblem::meanMotion() const { return std::sqrt(mu / (a * a * a)); }

State KeplerProblem::initialState() const {
    const double rp = a * (1.0 - e);
    const double vp = std::sqrt(mu * (1.0 + e) / rp);  // vis-viva au périastre
    return {rp, 0.0, 0.0, 0.0, vp, 0.0};
}

OdeFunction KeplerProblem::rhs() const {
    // Copie du paramètre : la fonction reste valide même si le problème est modifié ensuite.
    const double gm = mu;
    return [gm](double, const State& y, State& d) {
        const double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
        const double k = -gm / (r2 * std::sqrt(r2));  // -mu / r^3
        d[0] = y[3];
        d[1] = y[4];
        d[2] = y[5];
        d[3] = k * y[0];
        d[4] = k * y[1];
        d[5] = k * y[2];
    };
}

State KeplerProblem::exact(double t) const {
    const double n = meanMotion();
    const double E = kepler::solveEccentricAnomaly(n * t, e);
    const double c = std::cos(E), s = std::sin(E);
    const double b = a * std::sqrt(1.0 - e * e);   // demi-petit axe
    const double eDot = n / (1.0 - e * c);         // dE/dt, de M' = n = E' (1 - e cos E)
    return {a * (c - e), b * s, 0.0, -a * s * eDot, b * c * eDot, 0.0};
}

double KeplerProblem::energy(const State& y) const {
    const double r = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    const double v2 = y[3] * y[3] + y[4] * y[4] + y[5] * y[5];
    return 0.5 * v2 - mu / r;
}

Vec3 KeplerProblem::angularMomentum(const State& y) const {
    return cross(Vec3{y[0], y[1], y[2]}, Vec3{y[3], y[4], y[5]});
}

Vec3 KeplerProblem::rungeLenz(const State& y) const {
    const Vec3 r{y[0], y[1], y[2]}, v{y[3], y[4], y[5]};
    return cross(v, cross(r, v)) - mu * r.normalized();
}

double KeplerProblem::exactAngularMomentum() const { return std::sqrt(mu * a * (1.0 - e * e)); }

double KeplerProblem::periapsisAngle(const State& y) const {
    const Vec3 A = rungeLenz(y);
    return std::atan2(A.y, A.x);
}

double KeplerProblem::verletPrecessionPerOrbit(double dt) const {
    const double eta2 = 1.0 - e * e;  // (G/L)^2, G et L étant les variables de Delaunay
    return -(constants::pi / 8.0) * (mu * dt * dt / (a * a * a)) * (4.0 + e * e) / (eta2 * eta2 * eta2);
}

double KeplerProblem::distance(const State& s1, const State& s2) const {
    const double vScale = std::sqrt(mu / a);
    double sum = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double dp = (s1[k] - s2[k]) / a, dv = (s1[3 + k] - s2[3 + k]) / vScale;
        sum += dp * dp + dv * dv;
    }
    return std::sqrt(sum);
}

double keplerError(const KeplerProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);
    return problem.distance(y, problem.exact(t));
}

bool PeriapsisTracker::update(double t, const State& y) {
    bool crossed = false;
    if (havePrev_) {
        // s = r.v = r r' : négatif quand l'astre se rapproche, positif quand il s'éloigne.
        const double s0 = prev_[0] * prev_[3] + prev_[1] * prev_[4] + prev_[2] * prev_[5];
        const double s1 = y[0] * y[3] + y[1] * y[4] + y[2] * y[5];
        if (s0 < 0.0 && s1 >= 0.0) {
            const double h = t - tPrev_;
            const double u = s0 / (s0 - s1);  // instant du passage, interpolé linéairement entre les deux pas
            const double u2 = u * u, u3 = u2 * u;
            const double h00 = 2.0 * u3 - 3.0 * u2 + 1.0, h10 = u3 - 2.0 * u2 + u;
            const double h01 = -2.0 * u3 + 3.0 * u2, h11 = u3 - u2;
            const double px = h00 * prev_[0] + h10 * h * prev_[3] + h01 * y[0] + h11 * h * y[3];
            const double py = h00 * prev_[1] + h10 * h * prev_[4] + h01 * y[1] + h11 * h * y[4];
            const double angle = std::atan2(py, px);
            if (count_ == 0) {
                unwrapped_ = 0.0;
            } else {
                double d = angle - lastAngle_;
                d -= 2.0 * constants::pi * std::round(d / (2.0 * constants::pi));  // saut entre -pi et pi : replié
                unwrapped_ += d;
            }
            lastAngle_ = angle;
            ++count_;
            crossed = true;
        }
    }
    prev_ = y;
    tPrev_ = t;
    havePrev_ = true;
    return crossed;
}

double PeriapsisTracker::precessionPerOrbit() const {
    if (count_ < 2) return std::numeric_limits<double>::quiet_NaN();
    return unwrapped_ / (count_ - 1);
}

}  // namespace pl
