#include "physicslab/waves/Dispersion.hpp"

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

double growthPerStep1D(double x, double cfl) {
    const double s = std::sin(0.5 * x);
    return growthFromS2(cfl, s * s);
}

double growthPerStep2D(double kxdx, double kydx, double cfl) {
    const double sx = std::sin(0.5 * kxdx), sy = std::sin(0.5 * kydx);
    return growthFromS2(cfl, sx * sx + sy * sy);
}

double cflLimit(int dimension) { return dimension == 1 ? 1.0 : 0.70710678118654752440; }

}  // namespace pl::waves::dispersion
