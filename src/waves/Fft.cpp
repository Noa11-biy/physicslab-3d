#include "physicslab/waves/Fft.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "physicslab/core/Constants.hpp"

namespace pl::waves {
namespace {

// Cooley-Tukey itératif : permutation par inversion de bits, puis log2(n) étages de « papillons ».
// Les racines de l'unité sont calculées directement (cos, sin) à chaque étage, jamais par multiplications répétées
// (qui accumuleraient l'erreur d'arrondi).
void transform(std::vector<Complex>& a, bool inverse) {
    const std::size_t n = a.size();
    assert(isPowerOfTwo(n));
    if (n < 2) return;

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    std::vector<Complex> w;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2;
        const double angle = (inverse ? 2.0 : -2.0) * constants::pi / static_cast<double>(len);
        w.resize(half);
        for (std::size_t k = 0; k < half; ++k) w[k] = Complex(std::cos(angle * k), std::sin(angle * k));
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t k = 0; k < half; ++k) {
                const Complex u = a[i + k];
                const Complex v = a[i + k + half] * w[k];
                a[i + k] = u + v;
                a[i + k + half] = u - v;
            }
        }
    }

    if (inverse) {
        const double inv = 1.0 / static_cast<double>(n);
        for (Complex& z : a) z *= inv;
    }
}

}  // namespace

bool isPowerOfTwo(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

std::size_t nextPowerOfTwo(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void fft(std::vector<Complex>& a) { transform(a, false); }
void ifft(std::vector<Complex>& a) { transform(a, true); }

std::vector<Complex> dft(const std::vector<Complex>& x) {
    const std::size_t n = x.size();
    std::vector<Complex> out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex sum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            // l'indice j*k est réduit modulo n avant de former l'angle : pas de perte de précision pour les grands j*k
            const double angle = -2.0 * constants::pi * static_cast<double>((j * k) % n) / static_cast<double>(n);
            sum += x[j] * Complex(std::cos(angle), std::sin(angle));
        }
        out[k] = sum;
    }
    return out;
}

std::vector<double> amplitudeSpectrum(const std::vector<double>& x) {
    const std::size_t n = nextPowerOfTwo(x.size());
    std::vector<Complex> a(n, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < x.size(); ++i) a[i] = x[i];
    fft(a);

    const double inv = 1.0 / static_cast<double>(n);
    std::vector<double> spectrum(n / 2 + 1);
    for (std::size_t k = 0; k <= n / 2; ++k) {
        const bool edge = (k == 0 || k == n / 2);
        spectrum[k] = std::abs(a[k]) * inv * (edge ? 1.0 : 2.0);
    }
    return spectrum;
}

std::vector<double> windowedSpectrum(const std::vector<double>& samples, double sampleRate, int padding, double* binWidth) {
    const std::size_t count = samples.size();
    const std::size_t n = nextPowerOfTwo(count * static_cast<std::size_t>(std::max(padding, 1)));
    double mean = 0.0;
    for (double s : samples) mean += s;
    if (count > 0) mean /= static_cast<double>(count);

    std::vector<Complex> a(n, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < count; ++i) {
        const double hann = count > 1 ? 0.5 * (1.0 - std::cos(2.0 * constants::pi * static_cast<double>(i) / static_cast<double>(count - 1))) : 1.0;
        a[i] = (samples[i] - mean) * hann;
    }
    fft(a);

    std::vector<double> spectrum(n / 2 + 1);
    for (std::size_t k = 0; k <= n / 2; ++k) spectrum[k] = std::abs(a[k]);
    if (binWidth) *binWidth = sampleRate / static_cast<double>(n);
    return spectrum;
}

double peakFrequency(const std::vector<double>& samples, double sampleRate, double guess, double searchFraction) {
    if (samples.size() < 8) return 0.0;
    double df = 0.0;
    const std::vector<double> s = windowedSpectrum(samples, sampleRate, 4, &df);
    const double lo = guess * (1.0 - searchFraction), hi = guess * (1.0 + searchFraction);
    std::size_t best = 0;
    double bestValue = 0.0;
    for (std::size_t k = 1; k + 1 < s.size(); ++k) {
        const double f = static_cast<double>(k) * df;
        if (f < lo || f > hi) continue;
        if (s[k] > bestValue) { bestValue = s[k]; best = k; }
    }
    if (best == 0 || bestValue <= 0.0) return 0.0;

    // parabole sur ln|X| en k - 1, k, k + 1 : le sommet est décalé de delta cases
    const double alpha = std::log(std::max(s[best - 1], 1e-300)), beta = std::log(bestValue), gamma = std::log(std::max(s[best + 1], 1e-300));
    const double denom = alpha - 2.0 * beta + gamma;
    const double delta = denom != 0.0 ? 0.5 * (alpha - gamma) / denom : 0.0;
    return (static_cast<double>(best) + delta) * df;
}

}  // namespace pl::waves
