#include "physicslab/waves/Dispersion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pl::waves::dispersion {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// w dt = 2 asin(C S) avec S = sqrt(somme des sin²) ; NaN si C S > 1.
double phaseFromS(double cfl, double s, double kdx) {
    if (kdx == 0.0) return 1.0;
    if (cfl == 0.0) return s / (0.5 * kdx);  // limite C -> 0 : espace seul (1D) ; en 2D S est la norme du vecteur des sin
    const double arg = cfl * s;
    if (arg > 1.0) return kNaN;
    return 2.0 * std::asin(arg) / (cfl * kdx);
}

double growthFromS2(double cfl, double s2) {
    const double b = 1.0 - 2.0 * cfl * cfl * s2;
    if (std::abs(b) <= 1.0) return 1.0;
    return std::abs(b) + std::sqrt(b * b - 1.0);
}

}  // namespace

double phaseRatio1D(double x, double cfl) { return phaseFromS(cfl, std::sin(0.5 * x), x); }

double groupRatio1D(double x, double cfl) {
    const double s = cfl * std::sin(0.5 * x);
    if (s > 1.0) return kNaN;
    return std::cos(0.5 * x) / std::sqrt(1.0 - s * s);
}

double phaseRatio2D(double x, double theta, double cfl) {
    const double sx = std::sin(0.5 * x * std::cos(theta)), sy = std::sin(0.5 * x * std::sin(theta));
    return phaseFromS(cfl, std::sqrt(sx * sx + sy * sy), x);
}

double omegaDt2D(double kxdx, double kydx, double cfl) {
    const double sx = std::sin(0.5 * kxdx), sy = std::sin(0.5 * kydx);
    const double arg = cfl * std::sqrt(sx * sx + sy * sy);
    return arg > 1.0 ? kNaN : 2.0 * std::asin(arg);
}

void groupVelocity2D(double kxdx, double kydx, double cfl, double& vx, double& vy) {
    const double sx = std::sin(0.5 * kxdx), sy = std::sin(0.5 * kydx);
    const double s = std::sqrt(sx * sx + sy * sy), arg = cfl * s;
    if (s == 0.0 || arg >= 1.0) {
        vx = vy = arg >= 1.0 ? kNaN : 0.0;
        return;
    }
    const double common = 1.0 / (s * std::sqrt(1.0 - arg * arg));
    vx = sx * std::cos(0.5 * kxdx) * common;
    vy = sy * std::cos(0.5 * kydx) * common;
}

double growthPerStep1D(double x, double cfl) {
    const double s = std::sin(0.5 * x);
    return growthFromS2(cfl, s * s);
}

double growthPerStep2D(double kxdx, double kydx, double cfl) {
    const double sx = std::sin(0.5 * kxdx), sy = std::sin(0.5 * kydx);
    return growthFromS2(cfl, sx * sx + sy * sy);
}

double cflLimit(int dimension) { return dimension == 1 ? 1.0 : 0.70710678118654752440; }

RayPhase2D::RayPhase2D(double omegaDt, double cfl, double dx, int directions) {
    directions = std::max(directions, 2);
    constexpr double halfPi = 1.57079632679489661923;
    // |k| dx tel que w dt (kx dx, ky dx) = omegaDt dans la direction alpha du vecteur d'onde (w croît avec |k| tant que C S <= 1)
    auto wavevector = [&](double alpha, double& kxdx, double& kydx) {
        double lo = 0.0, hi = 3.1;
        for (int it = 0; it < 45; ++it) {
            const double m = 0.5 * (lo + hi), w = omegaDt2D(m * std::cos(alpha), m * std::sin(alpha), cfl);
            if (std::isnan(w) || w > omegaDt) hi = m; else lo = m;
        }
        kxdx = lo * std::cos(alpha);
        kydx = lo * std::sin(alpha);
    };
    f_.resize(directions);
    for (int d = 0; d < directions; ++d) {
        const double phi = halfPi * d / (directions - 1);  // direction du rayon (de la vitesse de groupe)
        double lo = 0.0, hi = halfPi, kx = 0.0, ky = 0.0;
        for (int it = 0; it < 45; ++it) {  // alpha -> direction de la vitesse de groupe : croissante de 0 à 90°
            const double m = 0.5 * (lo + hi);
            wavevector(m, kx, ky);
            double vx, vy;
            groupVelocity2D(kx, ky, cfl, vx, vy);
            if (std::atan2(vy, vx) < phi) lo = m; else hi = m;
        }
        wavevector(0.5 * (lo + hi), kx, ky);
        f_[d] = (kx * std::cos(phi) + ky * std::sin(phi)) / dx;
    }
}

double RayPhase2D::phase(double x, double y) const {
    const double r = std::hypot(x, y);
    if (r == 0.0) return 0.0;
    const double t = std::abs(std::atan2(y, x)) / 1.57079632679489661923 * (static_cast<double>(f_.size()) - 1.0);
    const std::size_t i = std::min(static_cast<std::size_t>(t), f_.size() - 2);
    return r * (f_[i] + (t - static_cast<double>(i)) * (f_[i + 1] - f_[i]));
}

}  // namespace pl::waves::dispersion
