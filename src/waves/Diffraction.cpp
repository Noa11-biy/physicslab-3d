#include "physicslab/waves/Diffraction.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

#include "physicslab/waves/Fft.hpp"

namespace pl::waves::diffraction {

double slitIntensity(double b, double a, int slits, double lambda, double s) {
    const double x = kPi * b * s / lambda;
    const double sinc = std::abs(x) < 1e-12 ? 1.0 : std::sin(x) / x;
    double grating = 1.0;
    if (slits > 1) {
        const double phi = kPi * a * s / lambda, sp = std::sin(phi);
        if (std::abs(sp) > 1e-12) {
            const double r = std::sin(slits * phi) / (slits * sp);
            grating = r * r;
        }
    }
    return sinc * sinc * grating;
}

double gridWavelength(double omega, double dt, double dx, double cfl) {
    const double kdx = 2.0 * std::asin(std::sin(0.5 * omega * dt) / cfl);
    return 2.0 * kPi * dx / kdx;
}

FarField fraunhofer(const std::vector<double>& aperture, double dy, double lambda, int padding) {
    assert(!aperture.empty() && padding >= 1);
    const std::size_t n = nextPowerOfTwo(aperture.size() * static_cast<std::size_t>(padding));
    std::vector<Complex> a(n, Complex(0.0, 0.0));
    for (std::size_t k = 0; k < aperture.size(); ++k) a[k] = aperture[k];
    fft(a);
    FarField out;
    double peak = 0.0;
    for (long m = -static_cast<long>(n / 2); m < static_cast<long>(n / 2); ++m) {
        const double s = lambda * static_cast<double>(m) / (static_cast<double>(n) * dy);
        if (std::abs(s) > 1.0) continue;
        const double v = std::norm(a[static_cast<std::size_t>((m + static_cast<long>(n)) % static_cast<long>(n))]);
        out.sinTheta.push_back(s);
        out.intensity.push_back(v);
        peak = std::max(peak, v);
    }
    if (peak > 0.0)
        for (double& v : out.intensity) v /= peak;
    return out;
}

double apertureAmplitude(const std::vector<double>& aperture, double dy, double k, double xScreen, double px, double py) {
    constexpr int sub = 4;
    std::complex<double> sum(0.0, 0.0);
    for (std::size_t j = 0; j < aperture.size(); ++j) {
        if (aperture[j] <= 0.0) continue;
        for (int s = 0; s < sub; ++s) {
            const double y = (static_cast<double>(j) - 0.5 + (s + 0.5) / sub) * dy;
            const double r = std::hypot(px - xScreen, py - y);
            sum += std::complex<double>(std::cyl_bessel_j(0.0, k * r), std::cyl_neumann(0.0, k * r)) * (dy / sub);
        }
    }
    return 0.5 * k * std::abs(sum);
}

double apertureAmplitudeGrid(const std::vector<double>& aperture, double dy, const dispersion::RayPhase2D& ray, double xScreen, double px, double py) {
    constexpr int sub = 4;
    std::complex<double> sum(0.0, 0.0);
    for (std::size_t j = 0; j < aperture.size(); ++j) {
        if (aperture[j] <= 0.0) continue;
        for (int s = 0; s < sub; ++s) {
            const double y = (static_cast<double>(j) - 0.5 + (s + 0.5) / sub) * dy;
            const double x = px - xScreen, v = py - y;
            sum += std::polar(1.0 / std::sqrt(std::hypot(x, v)), ray.phase(x, v)) * (dy / sub);
        }
    }
    return std::sqrt(ray.axisWavenumber() / (2.0 * kPi)) * std::abs(sum);
}

std::vector<double> findMinima(const std::vector<double>& x, const std::vector<double>& y) {
    std::vector<double> out;
    for (std::size_t i = 1; i + 1 < y.size(); ++i) {
        if (!(y[i] < y[i - 1] && y[i] <= y[i + 1])) continue;
        const double d = y[i - 1] - 2.0 * y[i] + y[i + 1];
        const double h = 0.5 * (x[i + 1] - x[i - 1]);
        out.push_back(d > 0.0 ? x[i] + h * (y[i - 1] - y[i + 1]) / (2.0 * d) : x[i]);
    }
    return out;
}

SlitScreen makeSlitScreen(const SlitScreenSpec& s) {
    SlitScreen out;
    Wave2DParams& p = out.params;
    p.dx = s.dx;
    p.cfl = s.cfl;
    p.speed = 1.0;
    p.cellsX = static_cast<int>(std::lround(s.lengthX / s.dx));
    p.cellsY = static_cast<int>(std::lround(s.lengthY / s.dx));
    p.left = p.right = p.bottom = p.top = Edge::Absorbing;
    p.spongeCells = static_cast<int>(std::lround(s.sponge / s.dx));

    const int nx = p.cellsX + 1, ny = p.cellsY + 1;
    out.screenColumn = static_cast<int>(std::lround(s.screenX / s.dx));
    out.screenColumns = std::max(1, static_cast<int>(std::lround(s.thickness / s.dx)));
    p.wall.assign(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny), 0);
    for (int j = 0; j < ny; ++j)
        for (int c = 0; c < out.screenColumns; ++c) p.wall[static_cast<std::size_t>(j) * nx + out.screenColumn + c] = 1;

    // fentes : W noeuds ouverts chacune, centres espacés de P noeuds, groupe centré sur le noeud du milieu du domaine
    const int W = std::max(1, static_cast<int>(std::lround(s.slitWidth / s.dx)));
    const int P = std::max(W + 1, static_cast<int>(std::lround(s.pitch / s.dx)));
    double centerSum = 0.0;
    const int first0 = static_cast<int>(std::lround(p.cellsY / 2.0 - 0.5 * (s.slits - 1) * P - 0.5 * (W - 1)));
    for (int m = 0; m < s.slits; ++m) {
        const int first = first0 + m * P;  // pas exact en nombre entier de cases
        for (int j = first; j < first + W; ++j)
            for (int c = 0; c < out.screenColumns; ++c) p.wall[static_cast<std::size_t>(j) * nx + out.screenColumn + c] = 0;
        centerSum += 0.5 * (2 * first + W - 1);
    }
    out.centerY = (s.slits > 0 ? centerSum / s.slits : 0.5 * p.cellsY) * s.dx;
    out.centerX = (out.screenColumn + out.screenColumns - 0.5) * s.dx;
    out.slitWidth = W * s.dx;
    out.pitch = P * s.dx;
    out.aperture.resize(static_cast<std::size_t>(ny));
    for (int j = 0; j < ny; ++j) out.aperture[static_cast<std::size_t>(j)] = p.wall[static_cast<std::size_t>(j) * nx + out.screenColumn] ? 0.0 : 1.0;

    out.omega = 2.0 * kPi * p.speed / s.lambda;
    out.source.x = s.sourceX;
    out.source.y0 = 0.0;
    out.source.y1 = s.lengthY;
    out.source.amplitude = s.amplitude;
    out.source.omega = out.omega;
    out.source.rampPeriods = 4.0;
    return out;
}

LockIn::LockIn(std::size_t nodes, double omega, double dt, long long steps) : omega_(omega), dt_(dt), steps_(steps), re_(nodes, 0.0), im_(nodes, 0.0) {}

std::vector<double> LockIn::amplitude() const {
    std::vector<double> a(re_.size(), 0.0);
    if (weights_ <= 0.0) return a;
    for (std::size_t k = 0; k < a.size(); ++k) a[k] = 2.0 * std::hypot(re_[k], im_[k]) / weights_;
    return a;
}

double sampleField(const std::vector<double>& f, int cellsX, int cellsY, double dx, double x, double y) {
    const double fx = std::clamp(x / dx, 0.0, static_cast<double>(cellsX)), fy = std::clamp(y / dx, 0.0, static_cast<double>(cellsY));
    const int i = std::min(static_cast<int>(fx), cellsX - 1), j = std::min(static_cast<int>(fy), cellsY - 1);
    const double tx = fx - i, ty = fy - j;
    const std::size_t stride = static_cast<std::size_t>(cellsX) + 1;
    const std::size_t k = static_cast<std::size_t>(j) * stride + static_cast<std::size_t>(i);
    return (1 - ty) * ((1 - tx) * f[k] + tx * f[k + 1]) + ty * ((1 - tx) * f[k + stride] + tx * f[k + stride + 1]);
}

std::vector<double> sampleArc(const std::vector<double>& field, int cellsX, int cellsY, double dx, double cx, double cy, double radius,
                              const std::vector<double>& thetas) {
    std::vector<double> out;
    out.reserve(thetas.size());
    for (double th : thetas) out.push_back(sampleField(field, cellsX, cellsY, dx, cx + radius * std::cos(th), cy + radius * std::sin(th)));
    return out;
}

}  // namespace pl::waves::diffraction
