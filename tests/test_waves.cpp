// Validation du domaine Ondes (O0) : FFT contre la transformée directe, équation d'onde 1D et 2D (solutions exactes,
// invariants, ordre de convergence, bords). Même esprit que test_core.cpp : pas de framework, un CHECK minimal ; les valeurs
// mesurées sont affichées pour qu'on puisse les relire (ne jamais se contenter d'un test qui passe).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/waves/Dispersion.hpp"
#include "physicslab/waves/Fft.hpp"
#include "physicslab/waves/String.hpp"
#include "physicslab/waves/Wave.hpp"

namespace {

int g_failures = 0;

void check(bool ok, const char* expr, int line) {
    if (!ok) {
        std::fprintf(stderr, "test_waves.cpp:%d : ECHEC : %s\n", line, expr);
        ++g_failures;
    }
}
#define CHECK(cond) check((cond), #cond, __LINE__)
#define CHECK_NEAR(a, b, tol) check(std::abs((a) - (b)) <= (tol), #a " ~= " #b, __LINE__)

using namespace pl::waves;
using pl::constants::pi;

// Générateur congruentiel : signaux « quelconques » mais identiques à chaque exécution.
struct Lcg {
    unsigned long long s = 12345;
    double next() {  // dans [-1, 1[
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) / 4503599627370496.0 - 1.0;
    }
};

double gaussian(double x, double center, double width) {
    const double z = (x - center) / width;
    return std::exp(-0.5 * z * z);
}

// ------------------------------------------------------------------ FFT -----

void testFft() {
    Lcg rng;
    for (std::size_t n : {1u, 2u, 4u, 8u, 16u, 64u, 256u}) {
        std::vector<Complex> x(n);
        double norm1 = 0.0;
        for (Complex& z : x) {
            z = Complex(rng.next(), rng.next());
            norm1 += std::abs(z);
        }
        std::vector<Complex> a = x;
        fft(a);
        const std::vector<Complex> ref = dft(x);
        double diff = 0.0;
        for (std::size_t k = 0; k < n; ++k) diff = std::max(diff, std::abs(a[k] - ref[k]));
        CHECK(diff <= 1e-13 * std::max(1.0, norm1));

        // Parseval : somme |x|² = (1/n) somme |X|²
        double ex = 0.0, eX = 0.0;
        for (std::size_t k = 0; k < n; ++k) { ex += std::norm(x[k]); eX += std::norm(a[k]); }
        CHECK_NEAR(ex, eX / static_cast<double>(n), 1e-13 * std::max(1.0, ex));

        // aller-retour
        ifft(a);
        double back = 0.0;
        for (std::size_t k = 0; k < n; ++k) back = std::max(back, std::abs(a[k] - x[k]));
        CHECK(back <= 1e-14);
        if (n == 256) std::printf("  FFT n = 256 : ecart a la DFT %.1e, aller-retour %.1e\n", diff, back);
    }

    // une impulsion a un spectre plat ; un signal décalé a le même module
    {
        std::vector<Complex> a(32, 0.0);
        a[0] = 1.0;
        fft(a);
        for (const Complex& z : a) CHECK_NEAR(std::abs(z - Complex(1.0, 0.0)), 0.0, 1e-15);
    }

    // sinusoïde d'amplitude 3 tombant exactement sur la fréquence 5 : une seule raie, de hauteur 3
    {
        const std::size_t n = 64;
        std::vector<double> x(n);
        for (std::size_t j = 0; j < n; ++j) x[j] = 3.0 * std::cos(2.0 * pi * 5.0 * j / n + 0.7) + 0.25;
        const std::vector<double> s = amplitudeSpectrum(x);
        CHECK(s.size() == n / 2 + 1);
        CHECK_NEAR(s[5], 3.0, 1e-13);
        CHECK_NEAR(s[0], 0.25, 1e-13);  // la moyenne
        double other = 0.0;
        for (std::size_t k = 1; k < s.size(); ++k)
            if (k != 5) other = std::max(other, s[k]);
        CHECK(other <= 1e-13);
        std::printf("  raie a k = 5 : %.15f, bruit hors raie %.1e\n", s[5], other);
    }

    // complétion par des zéros : 100 échantillons -> 128 -> 65 valeurs
    CHECK(amplitudeSpectrum(std::vector<double>(100, 1.0)).size() == 65);

    // gaussienne : la FFT échantillonnée reproduit la transformée continue sqrt(2 pi) w exp(-k² w² / 2) (convergence spectrale)
    {
        const std::size_t n = 512;
        const double dx = 0.05, w = 0.4;
        std::vector<Complex> a(n);
        for (std::size_t j = 0; j < n; ++j) {
            const double xj = (j < n / 2 ? static_cast<double>(j) : static_cast<double>(j) - static_cast<double>(n)) * dx;
            a[j] = std::exp(-0.5 * xj * xj / (w * w));
        }
        fft(a);
        double worst = 0.0;
        for (std::size_t k = 0; k < n / 2; ++k) {
            const double omega = 2.0 * pi * static_cast<double>(k) / (static_cast<double>(n) * dx);
            const double analytic = std::sqrt(2.0 * pi) * w * std::exp(-0.5 * omega * omega * w * w);
            worst = std::max(worst, std::abs(a[k].real() * dx - analytic));
            CHECK(std::abs(a[k].imag()) <= 1e-12);
        }
        CHECK(worst <= 1e-12);
        std::printf("  gaussienne : ecart a la transformee continue %.1e\n", worst);
    }
}

// ------------------------------------------------------------------ 1D ------

// Solution exacte d'une impulsion f lâchée au repos dans un milieu infini : moitié vers la droite, moitié vers la gauche.
double dAlembert(const std::function<double(double)>& f, double x, double ct) { return 0.5 * (f(x - ct) + f(x + ct)); }

double maxError1D(const Wave1D& w, const std::function<double(double, double)>& exact) {
    double e = 0.0;
    for (int i = 0; i < w.points(); ++i) e = std::max(e, std::abs(w.u()[i] - exact(w.x(i), w.time())));
    return e;
}

void testWave1DEnergy() {
    struct Case { Edge left, right; double cfl; const char* name; };
    const Case cases[] = {{Edge::Fixed, Edge::Fixed, 0.8, "fixe-fixe  C=0.8"},
                          {Edge::Free, Edge::Free, 0.8, "libre-libre C=0.8"},
                          {Edge::Fixed, Edge::Free, 0.3, "fixe-libre C=0.3"},
                          {Edge::Free, Edge::Fixed, 1.0, "libre-fixe C=1.0"}};
    for (const Case& c : cases) {
        Wave1DParams p;
        p.cells = 300;
        p.cfl = c.cfl;
        p.left = c.left;
        p.right = c.right;
        Wave1D w(p);
        Lcg rng;
        std::vector<double> noise(p.cells + 1), vel(p.cells + 1);
        for (double& v : noise) v = rng.next();
        for (double& v : vel) v = rng.next();
        w.setInitial([&](double x) { return noise[static_cast<int>(std::lround(x / w.dx()))]; },
                     [&](double x) { return vel[static_cast<int>(std::lround(x / w.dx()))]; });
        const double e0 = w.energy();
        double drift = 0.0;
        for (int k = 0; k < 5000; ++k) {
            w.step();
            drift = std::max(drift, std::abs(w.energy() - e0) / e0);
        }
        CHECK(drift <= 1e-12);
        std::printf("  energie 1D %-18s : derive relative max %.1e sur 5000 pas\n", c.name, drift);
    }
}

// Un mode propre de la grille est une solution exacte du schéma : sin²(w dt / 2) = C² sin²(k dx / 2).
void testWave1DModes() {
    struct Case { Edge left, right; double k; bool sine; const char* name; };
    const double length = 2.0;
    const Case cases[] = {{Edge::Fixed, Edge::Fixed, 3.0 * pi / length, true, "fixe-fixe n=3"},
                          {Edge::Free, Edge::Free, 4.0 * pi / length, false, "libre-libre n=4"},
                          {Edge::Fixed, Edge::Free, 2.5 * pi / length, true, "fixe-libre n=2.5"}};
    for (const Case& c : cases) {
        Wave1DParams p;
        p.cells = 200;
        p.length = length;
        p.cfl = 0.8;
        p.left = c.left;
        p.right = c.right;
        Wave1D w(p);
        const double k = c.k;
        w.setInitial([&](double x) { return c.sine ? std::sin(k * x) : std::cos(k * x); });
        w.advance(500);
        const double omega = 2.0 / w.dt() * std::asin(p.cfl * std::sin(0.5 * k * w.dx()));
        const double err = maxError1D(w, [&](double x, double t) { return (c.sine ? std::sin(k * x) : std::cos(k * x)) * std::cos(omega * t); });
        CHECK(err <= 1e-12);
        // la grille ralentit un peu l'onde : w < c k (dispersion numérique, détaillée en O2)
        CHECK(omega < p.speed * k);
        std::printf("  mode %-18s : ecart au mode discret %.1e (w/(ck) = %.8f)\n", c.name, err, omega / (p.speed * k));
    }
}

void testWave1DAlembert() {
    const auto f = [](double x) { return gaussian(x, 2.0, 0.1); };
    // à C = 1 le schéma est exact (propagation d'une case par pas) : rien à mesurer comme ordre
    {
        Wave1DParams p;
        p.cells = 400;
        p.length = 4.0;
        p.cfl = 1.0;
        Wave1D w(p);
        w.setInitial(f);
        w.advance(80);  // dt = dx = 0.01 : t = 0.8
        const double err = maxError1D(w, [&](double x, double t) { return dAlembert(f, x, p.speed * t); });
        CHECK(err <= 1e-13);
        std::printf("  d'Alembert C=1 : ecart %.1e (schema exact)\n", err);
    }
    // à C = 0.5 : ordre 2 (erreur divisée par 4 quand dx est divisé par 2), mesuré à t = 0.8 hors des bords
    {
        double errors[3];
        const int cells[3] = {200, 400, 800};
        for (int r = 0; r < 3; ++r) {
            Wave1DParams p;
            p.cells = cells[r];
            p.length = 4.0;
            p.cfl = 0.5;
            Wave1D w(p);
            w.setInitial(f);
            w.advance(static_cast<long long>(std::lround(0.8 / w.dt())));
            errors[r] = maxError1D(w, [&](double x, double t) { return dAlembert(f, x, p.speed * t); });
        }
        const double s1 = std::log2(errors[0] / errors[1]), s2 = std::log2(errors[1] / errors[2]);
        CHECK(s1 > 1.85 && s1 < 2.15);
        CHECK(s2 > 1.9 && s2 < 2.1);
        std::printf("  d'Alembert C=0.5 : erreurs %.3e, %.3e, %.3e ; ordres %.3f, %.3f\n", errors[0], errors[1], errors[2], s1, s2);
    }
}

// Réflexions : méthode des images. Mur fixe = prolongement impair 2L-périodique, bout libre = prolongement pair.
void testWave1DReflection() {
    const double length = 2.0, x0 = 0.5, width = 0.05;
    const auto f = [&](double x) { return gaussian(x, x0, width); };
    for (int kind = 0; kind < 3; ++kind) {
        const Edge left = (kind == 1) ? Edge::Free : Edge::Fixed;
        const Edge right = (kind == 0) ? Edge::Fixed : (kind == 1 ? Edge::Free : Edge::Free);
        Wave1DParams p;
        p.cells = 400;
        p.length = length;
        p.cfl = 1.0;
        p.left = left;
        p.right = right;
        Wave1D w(p);
        w.setInitial(f);
        // prolongement : réflexion en 0 (signe sI selon le bord gauche), en L (signe sR selon le bord droit)
        const double sL = left == Edge::Fixed ? -1.0 : 1.0, sR = right == Edge::Fixed ? -1.0 : 1.0;
        const auto extended = [&](double x) {
            double sum = 0.0;
            // images : f(x - 2mL...) construit par réflexions successives ; sommer les 4 familles suffit pour t < 5 L / c
            for (int m = -4; m <= 4; ++m) {
                const double shift = 2.0 * length * m;
                const double sgn = std::pow(sL * sR, m);  // chaque aller-retour complet applique sL sR
                sum += sgn * f(x - shift);                              // image directe
                sum += sgn * sL * f(-x - shift);                        // réfléchie sur le bord gauche
            }
            return sum;
        };
        double worst = 0.0;
        for (int stage = 1; stage <= 4; ++stage) {
            w.advance(170);
            const double ct = p.speed * w.time();
            worst = std::max(worst, maxError1D(w, [&](double x, double) { return 0.5 * (extended(x - ct) + extended(x + ct)); }));
        }
        CHECK(worst <= 1e-9);
        std::printf("  reflexions %s : ecart a la methode des images %.1e (t jusqu'a %.2f)\n",
                    kind == 0 ? "fixe-fixe  " : (kind == 1 ? "libre-libre" : "fixe-libre "), worst, w.time());
    }

    // signe de l'écho : un mur fixe renvoie une impulsion inversée, un bout libre la renvoie de même signe
    for (int free = 0; free < 2; ++free) {
        Wave1DParams p;
        p.cells = 400;
        p.length = length;
        p.cfl = 1.0;
        p.left = free ? Edge::Free : Edge::Fixed;
        Wave1D w(p);
        w.setInitial([&](double x) { return gaussian(x, 0.5, width); });
        w.advance(150);  // dt = 0.005 : t = 0.75, l'onde gauche (amplitude ½) s'est réfléchie en x = 0 et repasse en x = 0.25
        const double near = w.u()[static_cast<int>(std::lround(0.25 / w.dx()))];
        CHECK(free ? near > 0.4 : near < -0.4);
    }
}

void testWave1DSponge() {
    // Bosse gaussienne de largeur 0.1 au centre d'un domaine de 8 m ; elle atteint les bords à t = 4 et son écho (s'il y en a un)
    // serait revenu au centre à t = 8. On mesure l'énergie restante dans le domaine à t = 8.
    const auto f = [](double x) { return gaussian(x, 4.0, 0.1); };
    auto residual = [&](int layer, double strength, int order, double* identicalBefore) {
        Wave1DParams p;
        p.cells = 800;
        p.length = 8.0;
        p.cfl = 0.9;
        p.left = p.right = Edge::Absorbing;
        p.spongeCells = layer;
        p.spongeReflection = strength;
        p.spongeOrder = order;
        Wave1D sponge(p);
        sponge.setInitial(f);
        Wave1DParams q = p;
        q.left = q.right = Edge::Fixed;
        Wave1D wall(q);
        wall.setInitial(f);
        const double e0 = sponge.energy();

        // avant l'arrivée de l'onde dans la couche (t = 2 : impulsions en x = 2 et 6, couches pour x < 0.6 et x > 7.4) les deux calculs coïncident
        const long long early = std::lround(2.0 / sponge.dt());
        sponge.advance(early);
        wall.advance(early);
        if (identicalBefore) {
            *identicalBefore = 0.0;
            for (int i = 0; i < sponge.points(); ++i) *identicalBefore = std::max(*identicalBefore, std::abs(sponge.u()[i] - wall.u()[i]));
        }
        const long long end = std::lround(8.0 / sponge.dt());
        sponge.advance(end - early);
        wall.advance(end - early);
        CHECK_NEAR(wall.energy() / e0, 1.0, 1e-12);  // un mur garde tout
        return sponge.energy() / e0;
    };

    double same = 0.0;
    const double kept = residual(60, 0.03, 1, &same);         // défauts : couche de 0.6 m
    const double thick = residual(120, 0.03, 1, nullptr);     // couche deux fois plus épaisse
    const double hard = residual(60, 1e-4, 1, nullptr);       // éponge « plus dure » : sigma_max x 4
    const double smooth = residual(60, 0.03, 3, nullptr);     // profil plus lisse
    CHECK(same <= 1e-14);
    CHECK(kept < 2e-2);            // l'éponge absorbe plus de 98 % de l'énergie (un mur n'en absorbe aucune)
    CHECK(thick < kept / 3.0);     // le résidu chute vite avec l'épaisseur (mesuré : /4)
    CHECK(hard > kept);            // plus dure = pire (zone suramortie pour les grandes longueurs d'onde)
    CHECK(smooth > kept);          // plus lisse = pire aussi pour cette bosse
    std::printf("  eponge 1D (bosse de largeur 0.1, couche 0.6 m, R = 0.03) : energie restante %.2e ; couche 1.2 m %.2e ; R = 1e-4 %.2e ; profil s^3 %.2e ; "
                "avant l'arrivee, ecart %.1e\n", kept, thick, hard, smooth, same);
}

void testWave1DStability() {
    for (double cfl : {0.99, 1.0, 1.02}) {
        Wave1DParams p;
        p.cells = 200;
        p.cfl = cfl;
        Wave1D w(p);
        Lcg rng;
        std::vector<double> noise(p.cells + 1);
        for (double& v : noise) v = 1e-3 * rng.next();
        w.setInitial([&](double x) { return noise[static_cast<int>(std::lround(x / w.dx()))]; });
        w.advance(1500);
        const double m = w.maxAbs();
        std::printf("  stabilite 1D C = %.2f : max |u| apres 1500 pas = %.3e\n", cfl, m);
        if (cfl <= 1.0) CHECK(m < 1.0);
        else CHECK(m > 1e3);
        CHECK(w.stable() == (cfl <= 1.0));
    }
}

// ------------------------------------------------------------------ 2D ------

void testWave2DEnergy() {
    struct Case { Edge l, r, b, t; const char* name; };
    const Case cases[] = {{Edge::Fixed, Edge::Fixed, Edge::Fixed, Edge::Fixed, "4 murs"},
                          {Edge::Free, Edge::Free, Edge::Free, Edge::Free, "4 bords libres"},
                          {Edge::Free, Edge::Fixed, Edge::Fixed, Edge::Free, "mixte"}};
    for (const Case& c : cases) {
        Wave2DParams p;
        p.cellsX = 48;
        p.cellsY = 36;
        p.dx = 0.05;
        p.cfl = 0.6;
        p.left = c.l; p.right = c.r; p.bottom = c.b; p.top = c.t;
        Wave2D w(p);
        Lcg rng;
        std::vector<double> noise(w.u().size()), vel(w.u().size());
        for (double& v : noise) v = rng.next();
        for (double& v : vel) v = rng.next();
        auto at = [&](const std::vector<double>& a, double x, double y) {
            return a[w.index(static_cast<int>(std::lround(x / p.dx)), static_cast<int>(std::lround(y / p.dx)))];
        };
        w.setInitial([&](double x, double y) { return at(noise, x, y); }, [&](double x, double y) { return at(vel, x, y); });
        const double e0 = w.energy();
        double drift = 0.0;
        for (int k = 0; k < 1500; ++k) {
            w.step();
            drift = std::max(drift, std::abs(w.energy() - e0) / e0);
        }
        CHECK(drift <= 1e-12);
        std::printf("  energie 2D %-15s : derive relative max %.1e sur 1500 pas\n", c.name, drift);
    }
}

// Mode propre d'un rectangle fixe : sin²(w dt / 2) = C² (sin²(kx dx / 2) + sin²(ky dx / 2)), exact pour le schéma.
void testWave2DMode() {
    Wave2DParams p;
    p.cellsX = 60;
    p.cellsY = 40;
    p.dx = 0.05;
    p.cfl = 0.65;
    Wave2D w(p);
    const double lx = p.cellsX * p.dx, ly = p.cellsY * p.dx;
    const double kx = 2.0 * pi / lx, ky = 3.0 * pi / ly;
    w.setInitial([&](double x, double y) { return std::sin(kx * x) * std::sin(ky * y); });
    w.advance(300);
    const double s = std::sin(0.5 * kx * p.dx), t = std::sin(0.5 * ky * p.dx);
    const double omega = 2.0 / w.dt() * std::asin(p.cfl * std::sqrt(s * s + t * t));
    double err = 0.0;
    for (int j = 0; j <= p.cellsY; ++j)
        for (int i = 0; i <= p.cellsX; ++i)
            err = std::max(err, std::abs(w.u()[w.index(i, j)] - std::sin(kx * i * p.dx) * std::sin(ky * j * p.dx) * std::cos(omega * w.time())));
    CHECK(err <= 1e-12);
    const double exactOmega = p.speed * std::sqrt(kx * kx + ky * ky);
    CHECK(omega < exactOmega);
    std::printf("  mode 2D (2,3) : ecart au mode discret %.1e (w/(c|k|) = %.6f)\n", err, omega / exactOmega);
}

// Un carré avec une impulsion centrée sur la diagonale est symétrique par échange de x et y ; deux impulsions
// symétriques l'une de l'autre donnent des champs symétriques.
void testWave2DSymmetry() {
    Wave2DParams p;
    p.cellsX = p.cellsY = 64;
    p.dx = 0.05;
    p.cfl = 0.6;
    Wave2D a(p), b(p);
    a.setInitial([](double x, double y) { return gaussian(x, 1.0, 0.15) * gaussian(y, 2.3, 0.15); });
    b.setInitial([](double x, double y) { return gaussian(x, 2.3, 0.15) * gaussian(y, 1.0, 0.15); });
    a.advance(180);
    b.advance(180);
    double worst = 0.0;
    for (int j = 0; j <= p.cellsY; ++j)
        for (int i = 0; i <= p.cellsX; ++i) worst = std::max(worst, std::abs(a.u()[a.index(i, j)] - b.u()[b.index(j, i)]));
    CHECK(worst <= 1e-13);
    std::printf("  symetrie x<->y : ecart %.1e apres 180 pas\n", worst);
}

void testWave2DSponge() {
    Wave2DParams p;
    p.cellsX = p.cellsY = 160;
    p.dx = 0.0125;  // carré de 2 m
    p.cfl = 0.6;
    p.left = p.right = p.bottom = p.top = Edge::Absorbing;
    p.spongeCells = 30;  // couche de 0.375 m, bosse de largeur 0.08 : défauts de la bibliothèque (R = 0.03, profil linéaire)
    const auto bump = [](double x, double y) {
        const double r2 = (x - 1.0) * (x - 1.0) + (y - 1.0) * (y - 1.0);
        return std::exp(-0.5 * r2 / (0.08 * 0.08));
    };
    Wave2D sponge(p);
    sponge.setInitial(bump);
    Wave2DParams q = p;
    q.left = q.right = q.bottom = q.top = Edge::Fixed;
    Wave2D wall(q);
    wall.setInitial(bump);
    const double e0 = sponge.energy();
    const long long steps = std::lround(3.5 / sponge.dt());
    sponge.advance(steps);
    wall.advance(steps);
    const double kept = sponge.energy() / e0, wallKept = wall.energy() / e0;
    CHECK_NEAR(wallKept, 1.0, 1e-12);
    CHECK(kept < 5e-3);
    std::printf("  eponge 2D (bosse de largeur 0.08, couche 0.375 m) : energie restante %.2e apres 3.5 s (mur : %.12f)\n", kept, wallKept);
}

void testWave2DStability() {
    for (double cfl : {0.70, 0.72}) {
        Wave2DParams p;
        p.cellsX = p.cellsY = 60;
        p.dx = 0.05;
        p.cfl = cfl;
        Wave2D w(p);
        Lcg rng;
        std::vector<double> noise(w.u().size());
        for (double& v : noise) v = 1e-3 * rng.next();
        w.setInitial([&](double x, double y) { return noise[w.index(static_cast<int>(std::lround(x / p.dx)), static_cast<int>(std::lround(y / p.dx)))]; });
        w.advance(1500);
        const double m = w.maxAbs();
        std::printf("  stabilite 2D C = %.2f (limite %.4f) : max |u| apres 1500 pas = %.3e\n", cfl, Wave2D::kCflLimit, m);
        if (cfl < Wave2D::kCflLimit) CHECK(m < 1.0);
        else CHECK(m > 1e3);
        CHECK(w.stable() == (cfl < Wave2D::kCflLimit));
    }
}

// ------------------------------------------------------- corde (O1) --------

StringProblem makeString(int beads) {
    StringProblem p;
    p.beads = beads;
    p.length = 1.0;
    p.tension = 4.0;
    p.density = 1.0;  // c = 2 m/s, f1 = 1 Hz
    return p;
}

double maxAbsPositions(const pl::State& y, int beads) {
    double m = 0.0;
    for (int j = 0; j < beads; ++j) {
        const double a = std::abs(y[j]);
        if (!(a <= m)) m = a;
    }
    return m;
}

void testStringFormulas() {
    const StringProblem p = makeString(50);
    CHECK_NEAR(p.speed(), 2.0, 1e-15);
    CHECK_NEAR(p.spacing(), 1.0 / 51.0, 1e-16);
    CHECK_NEAR(p.maxOmega(), 4.0 * 51.0, 1e-9);                       // 2 c / a
    CHECK_NEAR(p.verletDtLimit(), p.spacing() / p.speed(), 1e-16);    // dt_max = a / c : C = 1
    CHECK_NEAR(p.fundamental(), 1.0, 1e-15);
    // dispersion de la chaîne : w_n / (n pi c / L) = sin(t) / t avec t = n pi / (2 (N + 1))
    for (int n : {1, 5, 25, 50}) {
        const double t = n * pi / (2.0 * 51.0);
        CHECK_NEAR(p.chainOmega(n) / p.continuumOmega(n), std::sin(t) / t, 1e-14);
    }
    std::printf("  chaine N = 50 : w_1/(pi c/L) = %.8f, w_50/(50 pi c/L) = %.6f, coupure %.1f rad/s, dt max Verlet %.6f s, RK4 %.6f s\n",
                p.chainOmega(1) / p.continuumOmega(1), p.chainOmega(50) / p.continuumOmega(50), p.maxOmega(), p.verletDtLimit(), p.rk4DtLimit());
    // la limite continue : N grand, les fréquences deviennent celles de la corde
    const StringProblem big = makeString(3000);
    CHECK(std::abs(big.chainOmega(3) / big.continuumOmega(3) - 1.0) < 1e-6);
}

void testStringModes() {
    const StringProblem p = makeString(40);
    const pl::OdeFunction f = p.rhs();
    // un mode est un vecteur propre : l'accélération vaut - w_n² u
    for (int n : {1, 2, 7, 40}) {
        const std::vector<double> u = p.modeShape(n, 0.3);
        pl::State y = StringProblem::state(u, std::vector<double>(p.beads, 0.0)), dy(y.size());
        f(0.0, y, dy);
        const double w2 = p.chainOmega(n) * p.chainOmega(n);
        double worst = 0.0;
        for (int j = 0; j < p.beads; ++j) worst = std::max(worst, std::abs(dy[p.beads + j] + w2 * u[j]));
        CHECK(worst <= 1e-12 * w2);
    }

    // l'énergie se range exactement par mode : somme des E_n = E
    const StringModes modes(p);
    Lcg rng;
    std::vector<double> u(p.beads), v(p.beads);
    for (double& x : u) x = rng.next();
    for (double& x : v) x = rng.next();
    const pl::State y = StringProblem::state(u, v);
    double sum = 0.0;
    for (double e : modes.modeEnergies(y)) sum += e;
    CHECK_NEAR(sum, p.energy(y), 1e-12 * p.energy(y));

    // un mode pur garde toute son énergie dans ce mode
    const pl::State pure = StringProblem::state(p.modeShape(3, 0.01), std::vector<double>(p.beads, 0.0));
    const std::vector<double> e = modes.modeEnergies(pure);
    double others = 0.0;
    for (int n = 1; n <= p.beads; ++n)
        if (n != 3) others = std::max(others, e[n - 1]);
    CHECK_NEAR(e[2], p.energy(pure), 1e-12 * p.energy(pure));
    CHECK(others <= 1e-25 * e[2]);
    std::printf("  energie par mode : somme / E - 1 = %.1e ; mode pur 3, autres modes <= %.1e x E_3\n", sum / p.energy(y) - 1.0, others / e[2]);
}

void testStringExact() {
    const StringProblem p = makeString(30);
    const StringModes modes(p);
    const std::vector<double> u0 = p.pluck(0.3, 0.1), v0(p.beads, 0.0);

    const pl::State y0 = modes.exact(u0, v0, 0.0);
    double init = 0.0;
    for (int j = 0; j < p.beads; ++j) init = std::max(init, std::abs(y0[j] - u0[j]));
    CHECK(init <= 1e-15);

    // contre une référence indépendante : RK45 serré
    pl::RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-15;
    pl::State y = StringProblem::state(u0, v0);
    pl::advance(rk, p.rhs(), 0.0, y, 1.7);
    const pl::State ex = modes.exact(u0, v0, 1.7);
    double err = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i) err = std::max(err, std::abs(y[i] - ex[i]));
    CHECK(err <= 1e-9);
    // l'énergie de la solution exacte ne dépend pas du temps
    CHECK_NEAR(p.energy(ex), p.energy(y0), 1e-12 * p.energy(y0));
    std::printf("  solution exacte de la chaine (N = 30, t = 1.7) : ecart a RK45 serre %.1e, E(t)/E(0) - 1 = %.1e\n", err, p.energy(ex) / p.energy(y0) - 1.0);
}

void testStringContinuum() {
    // série de Fourier du pincement : à t = 0 elle redonne le triangle ; les harmoniques pairs d'un pincement au milieu sont absents
    const StringProblem p = makeString(400);
    const double x0 = 0.3, h = 0.1;
    double worst0 = 0.0;
    for (double x : {0.1, 0.2, 0.3, 0.55, 0.8}) {
        const double tri = x < x0 ? h * x / x0 : h * (1.0 - x) / (1.0 - x0);
        worst0 = std::max(worst0, std::abs(p.continuumPluck(x, 0.0, x0, h, 4000) - tri));
    }
    CHECK(worst0 <= 2e-4);
    CHECK(std::abs(p.pluckCoefficient(2, 0.5, h)) <= 1e-15);
    CHECK(std::abs(p.pluckCoefficient(4, 0.5, h)) <= 1e-15);
    CHECK(std::abs(p.pluckCoefficient(1, 0.5, h)) > 0.01);

    // la chaîne de 400 masses suit la corde continue à 1 % de la hauteur près : l'écart (mesuré 5,6e-4) n'est pas une erreur de schéma
    // mais la DISPERSION de la chaîne : les modes n >~ 30, trop graves, se déphasent au bout de 0,37 s (b_n décroît seulement en 1 / n²)
    const StringModes modes(p);
    const std::vector<double> u0 = p.pluck(x0, h), v0(p.beads, 0.0);
    const double t = 0.37;
    const pl::State y = modes.exact(u0, v0, t);
    double worst = 0.0;
    for (int j = 1; j <= p.beads; ++j) worst = std::max(worst, std::abs(y[j - 1] - p.continuumPluck(j * p.spacing(), t, x0, h, 4000)));
    CHECK(worst <= 1e-3);
    std::printf("  serie de Fourier a t = 0 : ecart au triangle %.1e (4000 termes) ; chaine N = 400 contre corde continue a t = 0.37 : %.1e (h = 0.1)\n", worst0, worst);
}

void testStringHarmonics() {
    // N + 1 = 48 : un pincement à L/3 ou L/2 tombe exactement sur une masse, les modes multiples de 3 (resp. pairs) sont absents
    const StringProblem p = makeString(47);
    const StringModes modes(p);
    const std::vector<double> third = modes.amplitudes(p.pluck(1.0 / 3.0, 0.1));
    const std::vector<double> half = modes.amplitudes(p.pluck(0.5, 0.1));
    double missing3 = 0.0, missing2 = 0.0;
    for (int n = 3; n <= 47; n += 3) missing3 = std::max(missing3, std::abs(third[n - 1]));
    for (int n = 2; n <= 46; n += 2) missing2 = std::max(missing2, std::abs(half[n - 1]));
    CHECK(missing3 <= 1e-15);
    CHECK(missing2 <= 1e-15);
    CHECK(std::abs(third[0]) > 1e-3 && std::abs(third[1]) > 1e-4 && std::abs(half[0]) > 1e-3 && std::abs(half[2]) > 1e-5);
    std::printf("  harmoniques absents : pincement a L/3 -> modes 3, 6, 9... <= %.1e (mode 1 : %.3e) ; a L/2 -> modes pairs <= %.1e (mode 3 : %.3e)\n",
                missing3, third[0], missing2, half[2]);
}

void testStringSolvers() {
    const StringProblem p = makeString(20);
    const std::vector<double> u0 = p.pluck(0.3, 0.1);
    const double tEnd = 2.7;  // pas un multiple de la période (1 s) : voir oscillatorError
    struct Row { const char* name; int first; double expected; double tolerance; };
    const Row rows[] = {{"Euler symplectique", 0, 1.0, 0.25}, {"Verlet", 1, 2.0, 0.1}, {"RK4", 3, 4.0, 0.3}};
    for (const Row& row : rows) {
        // pré-asymptotique à 400-1600 pas (w_max dt jusqu'à 0.57 : ordres apparents 1.4 / 1.5) ; l'ordre se lit sur les pas fins
        double e[4];
        const int steps[4] = {800, 1600, 3200, 6400};
        for (int r = 0; r < 4; ++r) {
            std::unique_ptr<pl::Solver> solver;
            if (row.first == 0) solver = std::make_unique<pl::SymplecticEuler>();
            else if (row.first == 1) solver = std::make_unique<pl::VelocityVerlet>();
            else solver = std::make_unique<pl::RK4>();
            e[r] = stringError(p, u0, *solver, steps[r], tEnd);
        }
        const double s1 = std::log2(e[0] / e[1]), s2 = std::log2(e[1] / e[2]), s3 = std::log2(e[2] / e[3]);
        CHECK(std::abs(s3 - row.expected) < row.tolerance);
        std::printf("  chaine N = 20, t = 2.7 : %-18s erreurs %.3e, %.3e, %.3e, %.3e ; ordres %.3f, %.3f, %.3f\n", row.name, e[0], e[1], e[2], e[3], s1, s2, s3);
    }
}

void testStringStability() {
    // la pulsation de coupure fixe le pas maximal : Verlet pour w_max dt < 2, RK4 pour w_max dt < 2 sqrt(2)
    const StringProblem p = makeString(40);
    Lcg rng;
    std::vector<double> noise(p.beads);
    for (double& x : noise) x = 1e-3 * rng.next();
    const pl::OdeFunction f = p.rhs();
    struct Case { const char* name; double limit; bool verlet; };
    const Case cases[] = {{"Verlet", p.verletDtLimit(), true}, {"RK4", p.rk4DtLimit(), false}};
    for (const Case& c : cases) {
        for (double factor : {0.98, 1.1}) {
            pl::State y = StringProblem::state(noise, std::vector<double>(p.beads, 0.0));
            std::unique_ptr<pl::Solver> solver;
            if (c.verlet) solver = std::make_unique<pl::VelocityVerlet>();
            else solver = std::make_unique<pl::RK4>();
            const double dt = factor * c.limit;
            double t = 0.0;
            for (int i = 0; i < 400; ++i) { solver->step(f, t, y, dt); t += dt; }
            const double m = maxAbsPositions(y, p.beads);
            std::printf("  stabilite %-6s dt = %.2f x limite : max |u| apres 400 pas = %.3e\n", c.name, factor, m);
            if (factor < 1.0) CHECK(m < 1.0);
            else CHECK(m > 1e3);
        }
    }
}

// Le saute-mouton de l'équation d'onde (O0) est EXACTEMENT le Verlet de la chaîne de masses : mêmes nombres à l'arrondi près.
void testStringIsLeapfrog() {
    StringProblem p;
    p.beads = 99;
    p.length = 1.0;
    p.tension = 1.0;
    p.density = 1.0;  // c = 1, a = 0.01
    const double cfl = 0.7, dt = cfl * p.spacing() / p.speed();
    const std::vector<double> u0 = p.pluck(0.3, 0.1);

    Wave1DParams wp;
    wp.cells = 100;
    wp.length = 1.0;
    wp.speed = 1.0;
    wp.cfl = cfl;
    wp.left = wp.right = Edge::Fixed;
    Wave1D wave(wp);
    wave.setInitial([&](double x) { return x < 0.3 ? 0.1 * x / 0.3 : 0.1 * (1.0 - x) / 0.7; });

    pl::State y = StringProblem::state(u0, std::vector<double>(p.beads, 0.0));
    const pl::OdeFunction f = p.rhs();
    pl::VelocityVerlet verlet;
    double t = 0.0, worst = 0.0;
    for (int i = 0; i < 300; ++i) {
        verlet.step(f, t, y, dt);
        t += dt;
        wave.step();
        for (int j = 1; j <= p.beads; ++j) worst = std::max(worst, std::abs(y[j - 1] - wave.u()[j]));
    }
    CHECK(worst <= 1e-13);
    std::printf("  Verlet de la chaine = saute-mouton de l'equation d'onde : ecart max %.1e sur 300 pas (C = 0.7)\n", worst);
}

void testStringFrequencies() {
    // les fréquences lues dans le SPECTRE du mouvement d'une masse sont celles des modes de la chaîne
    const StringProblem p = makeString(60);
    const StringModes modes(p);
    const std::vector<double> u0 = p.pluck(0.37, 0.1), v0(p.beads, 0.0);
    const int probe = 17;        // x = 0.295 L : ne tombe sur aucun noeud des modes 1 à 6
    const double rate = 200.0;   // échantillons par seconde (Nyquist 100 Hz, coupure de la chaîne 38.8 Hz)
    const int count = 16384;
    const double dt = 1.0 / rate;

    std::vector<double> exactSeries(count), rk4Series(count);
    {
        // le même calcul avec RK4 à 5 ms (stable : limite 11.5 ms), un échantillon par pas
        pl::RK4 rk4;
        pl::State y = StringProblem::state(u0, v0);
        const pl::OdeFunction f = p.rhs();
        double t = 0.0;
        for (int i = 0; i < count; ++i) {
            exactSeries[i] = modes.exact(u0, v0, i * dt)[probe];
            rk4Series[i] = y[probe];
            rk4.step(f, t, y, dt);
            t += dt;
        }
    }
    double worstExact = 0.0, worstRk4 = 0.0;
    std::printf("  frequences lues dans le spectre (masse %d, 82 s, Hann + parabole) :\n    n   chaine exacte   spectre exact   spectre RK4    corde continue\n", probe);
    for (int n = 1; n <= 6; ++n) {
        const double chain = p.chainOmega(n) / (2.0 * pi), continuum = n * p.fundamental();
        // fenêtre de ±10 % : les modes 4 et 5 ne sont séparés que de 25 % (une fenêtre de ±20 % avait pris le pic du mode 4 pour le 5)
        const double fe = peakFrequency(exactSeries, rate, chain, 0.1), fr = peakFrequency(rk4Series, rate, chain, 0.1);
        worstExact = std::max(worstExact, std::abs(fe - chain) / chain);
        worstRk4 = std::max(worstRk4, std::abs(fr - chain) / chain);
        std::printf("    %d   %.6f Hz    %.6f Hz   %.6f Hz   %.6f Hz\n", n, chain, fe, fr, continuum);
    }
    CHECK(worstExact <= 5e-4);
    CHECK(worstRk4 <= 5e-4);
    std::printf("  ecart relatif max a la chaine : spectre exact %.1e, spectre RK4 %.1e\n", worstExact, worstRk4);
}

// ------------------------------------------- dispersion et CFL (O2) --------

namespace dsp = pl::waves::dispersion;

void testDispersionFormulas() {
    // C = 1 : la vitesse de phase et celle de groupe valent exactement c (erreurs d'espace et de temps compensées)
    for (double x : {0.1, 0.5, 1.0, 2.0, 3.0}) {
        CHECK_NEAR(dsp::phaseRatio1D(x, 1.0), 1.0, 1e-14);
        CHECK_NEAR(dsp::groupRatio1D(x, 1.0), 1.0, 1e-14);
    }
    // petite phase : 1 - (1 - C²) x² / 24 + O(x^4)
    {
        const double x = 0.1, c = 0.5;
        CHECK_NEAR(dsp::phaseRatio1D(x, c), 1.0 - (1.0 - c * c) * x * x / 24.0, 1e-7);
    }
    // la vitesse de groupe est la dérivée de w : w dt / ... = x * phaseRatio ; différences finies centrées
    for (double c : {0.3, 0.9}) {
        for (double x : {0.4, 1.0, 2.0}) {
            const double h = 1e-5;
            const double fd = ((x + h) * dsp::phaseRatio1D(x + h, c) - (x - h) * dsp::phaseRatio1D(x - h, c)) / (2.0 * h);
            CHECK_NEAR(fd, dsp::groupRatio1D(x, c), 1e-8);
        }
    }
    // C -> 0 : espace seul, c'est la chaîne de masses d'O1
    for (double x : {0.5, 1.5, 3.0}) {
        CHECK_NEAR(dsp::phaseRatio1D(x, 0.0), std::sin(0.5 * x) / (0.5 * x), 1e-15);
        CHECK_NEAR(dsp::phaseRatio1D(x, 1e-9), dsp::phaseRatio1D(x, 0.0), 1e-12);
        CHECK_NEAR(dsp::groupRatio1D(x, 0.0), std::cos(0.5 * x), 1e-15);
    }
    // onde non propagée : NaN
    CHECK(std::isnan(dsp::phaseRatio1D(3.0, 1.2)));
    // 2D : sur l'axe on retrouve le 1D ; l'écart au continu vaut (x² / 24)(cos⁴ + sin⁴ - C²) : plus faible en diagonale
    for (double x : {0.3, 1.0}) CHECK_NEAR(dsp::phaseRatio2D(x, 0.0, 0.6), dsp::phaseRatio1D(x, 0.6), 1e-15);
    {
        const double x = 0.1, c = 0.6, a = 0.3;
        const double axis = dsp::phaseRatio2D(x, a * 0.0, c), diag = dsp::phaseRatio2D(x, dsp::kPi / 4.0, c);
        CHECK_NEAR(axis, 1.0 - x * x / 24.0 * (1.0 - c * c), 1e-7);
        CHECK_NEAR(diag, 1.0 - x * x / 24.0 * (0.5 - c * c), 1e-7);
        // à C = 1/sqrt(2) la grille 2D est exacte (à l'ordre dominant) EN DIAGONALE
        CHECK_NEAR(dsp::phaseRatio2D(x, dsp::kPi / 4.0, dsp::cflLimit(2)), 1.0, 1e-7);
    }
    // croissance : 1 si stable, sinon (C + sqrt(C² - 1))² pour le mode de Nyquist en 1D
    CHECK_NEAR(dsp::growthPerStep1D(2.0, 1.0), 1.0, 0.0);
    {
        const double c = 1.05, x = dsp::kPi;
        const double expected = (c + std::sqrt(c * c - 1.0)) * (c + std::sqrt(c * c - 1.0));
        CHECK_NEAR(dsp::growthPerStep1D(x, c), expected, 1e-13);
    }
    std::printf("  C = 0.5, x = k dx = 1 (6.3 cases par longueur d'onde) : v_phase / c = %.6f, v_group / c = %.6f ; chaine (C -> 0) : %.6f, %.6f\n",
                dsp::phaseRatio1D(1.0, 0.5), dsp::groupRatio1D(1.0, 0.5), dsp::phaseRatio1D(1.0, 0.0), dsp::groupRatio1D(1.0, 0.0));
    std::printf("  2D, C = 0.6, x = 1 : v_phase / c = %.6f sur l'axe, %.6f en diagonale ; C = 1/sqrt(2) : diagonale %.6f (x = 0.1 : %.8f)\n",
                dsp::phaseRatio2D(1.0, 0.0, 0.6), dsp::phaseRatio2D(1.0, dsp::kPi / 4.0, 0.6), dsp::phaseRatio2D(1.0, dsp::kPi / 4.0, dsp::cflLimit(2)),
                dsp::phaseRatio2D(0.1, dsp::kPi / 4.0, dsp::cflLimit(2)));
}

// Un paquet d'ondes (porteuse k0 sous une enveloppe gaussienne large) avance à la vitesse de GROUPE de la grille, pas à c.
void testGroupVelocity() {
    const double k0dx = 1.0, x0 = 8.0, width = 0.5, duration = 6.0;
    std::printf("  paquet d'ondes k0 dx = %.1f (%.1f cases par longueur d'onde), enveloppe 0.5 m, t = 6 s :\n    C     v_groupe/c (formule)   mesure   ecart relatif\n", k0dx, dsp::pointsPerWavelength(k0dx));
    for (double cfl : {1.0, 0.9, 0.5, 0.2}) {
        Wave1DParams p;
        p.cells = 2000;
        p.length = 20.0;
        p.speed = 1.0;
        p.cfl = cfl;
        Wave1D w(p);
        const double k0 = k0dx / w.dx();
        const auto env = [&](double x) { return std::exp(-0.5 * (x - x0) * (x - x0) / (width * width)); };
        // Paquet qui va PUREMENT vers la droite, avec la relation de dispersion de la grille : u^n_i = env(x_i - v_g t_n) cos(k0 x_i - w t_n).
        // (Avec setInitial, la vitesse initiale -c f' est celle du continu : à 6 cases par longueur d'onde elle laisse 8 % d'onde de retour,
        //  et le centre mesuré à C = 1 était 1,5 % trop lent alors que la grille y est exacte.)
        const double vg = p.speed * dsp::groupRatio1D(k0dx, cfl);
        const double omegaDt = 2.0 * std::asin(cfl * std::sin(0.5 * k0dx));
        w.setStates([&](double x) { return env(x + vg * w.dt()) * std::cos(k0 * (x - x0) + omegaDt); },
                    [&](double x) { return env(x) * std::cos(k0 * (x - x0)); });
        w.advance(std::lround(duration / w.dt()));
        double sum = 0.0, sumx = 0.0;
        for (int i = 0; i < w.points(); ++i) {
            const double u2 = w.u()[i] * w.u()[i];
            sum += u2;
            sumx += w.x(i) * u2;
        }
        const double measured = (sumx / sum - x0) / (w.time() * p.speed);
        const double formula = dsp::groupRatio1D(k0dx, cfl);
        CHECK(std::abs(measured / formula - 1.0) <= 3e-3);
        std::printf("    %.1f   %.6f               %.6f   %+.1e\n", cfl, formula, measured, measured / formula - 1.0);
    }
}

// À C = 1 un paquet cohérent avec la grille est simplement décalé d'une case par pas : solution exacte du schéma, à l'arrondi près.
void testSetStatesShift() {
    Wave1DParams p;
    p.cells = 600;
    p.length = 6.0;
    p.cfl = 1.0;
    Wave1D w(p);
    const double dx = w.dx(), k0 = 1.0 / dx;
    // enveloppe étroite et loin des murs : à 8 largeurs du bord la queue gaussienne vaut 1e-14 (à 5 largeurs elle valait 3,7e-6 et le mur fixe
    // imposé en x = 0 donnait exactement cet écart)
    const auto f = [&](double x) { return gaussian(x, 2.0, 0.25) * std::cos(k0 * (x - 2.0)); };
    w.setStates([&](double x) { return f(x + dx); }, f);
    w.advance(150);
    double err = 0.0;
    for (int i = 0; i < w.points(); ++i) err = std::max(err, std::abs(w.u()[i] - f(w.x(i) - 150.0 * dx)));
    CHECK(err <= 1e-12);
    std::printf("  setStates, C = 1 : paquet decale de 150 cases, ecart a f(x - ct) %.1e\n", err);
}

// Fréquence d'un mode propre lue dans le spectre d'une simulation (et non plus par la formule seule) : confirme la vitesse de phase.
void testPhaseVelocityFromSpectrum() {
    const int cells = 200;
    const double length = 2.0, cfl = 0.5;
    std::printf("  vitesse de phase lue dans le spectre (1D, %d cases, C = %.1f) :\n    n   x = k dx   v_phase/c formule   spectre   ecart\n", cells, cfl);
    for (int n : {10, 40, 80, 120}) {
        Wave1DParams p;
        p.cells = cells;
        p.length = length;
        p.cfl = cfl;
        Wave1D w(p);
        const double k = n * dsp::kPi / length;
        w.setInitial([&](double x) { return std::sin(k * x); });
        const int probe = 63;  // x = 0.63 : loin des noeuds
        std::vector<double> series;
        for (int i = 0; i < 8192; ++i) {
            series.push_back(w.u()[probe]);
            w.step();
        }
        const double x = k * w.dx();
        const double predicted = dsp::phaseRatio1D(x, cfl) * p.speed * k / (2.0 * dsp::kPi);  // Hz
        const double read = peakFrequency(series, 1.0 / w.dt(), predicted, 0.05);
        const double ratio = read * 2.0 * dsp::kPi / (p.speed * k);
        CHECK(std::abs(read / predicted - 1.0) <= 2e-5);
        std::printf("    %-3d %.4f    %.6f           %.6f  %+.1e\n", n, x, dsp::phaseRatio1D(x, cfl), ratio, read / predicted - 1.0);
    }
}

// La grille 2D est anisotrope : à |k| égal, une onde en diagonale va plus vite qu'une onde le long d'un axe.
void testAnisotropy2D() {
    Wave2DParams p;
    p.cellsX = p.cellsY = 60;
    p.dx = 0.05;
    p.cfl = 0.65;
    const double lx = p.cellsX * p.dx;
    std::printf("  anisotropie 2D (60 x 60 cases, C = %.2f) : modes (m, n) de |k| presque egal\n    (m, n)   |k| dx   angle   v_phase/c formule   spectre   ecart\n", p.cfl);
    double axisRatio = 0.0, diagRatio = 0.0;
    for (const auto& mn : std::vector<std::pair<int, int>>{{24, 1}, {17, 17}, {12, 22}}) {
        Wave2D w(p);
        const double kx = mn.first * dsp::kPi / lx, ky = mn.second * dsp::kPi / lx;
        w.setInitial([&](double x, double y) { return std::sin(kx * x) * std::sin(ky * y); });
        const int pi_ = 13, pj = 19;
        std::vector<double> series;
        double amp = 0.0;
        for (int i = 0; i < 8192; ++i) {
            series.push_back(w.u()[w.index(pi_, pj)]);
            amp = std::max(amp, std::abs(series.back()));
            w.step();
        }
        CHECK(amp > 1e-3);
        const double kk = std::sqrt(kx * kx + ky * ky), theta = std::atan2(ky, kx), x = kk * p.dx;
        const double formula = dsp::phaseRatio2D(x, theta, p.cfl);
        const double predicted = formula * p.speed * kk / (2.0 * dsp::kPi);
        const double read = peakFrequency(series, 1.0 / w.dt(), predicted, 0.05);
        const double ratio = read * 2.0 * dsp::kPi / (p.speed * kk);
        CHECK(std::abs(read / predicted - 1.0) <= 5e-5);
        if (mn.first == 24) axisRatio = ratio;
        if (mn.first == 17) diagRatio = ratio;
        std::printf("    (%d, %d)  %.4f  %5.1f   %.6f           %.6f  %+.1e\n", mn.first, mn.second, x, theta * 180.0 / dsp::kPi, formula, ratio, read / predicted - 1.0);
    }
    CHECK(diagRatio > axisRatio + 0.02);  // la diagonale est nettement moins dispersive à ce C
}

// Au-delà de la limite CFL, un mode est amplifié à chaque pas du facteur prévu par la théorie (mesuré sur la simulation).
void testGrowthRates() {
    auto geometricGrowth = [](const std::vector<double>& history, int last) {
        const int n = static_cast<int>(history.size());
        return std::pow(history[n - 1] / history[n - 1 - last], 1.0 / last);
    };
    std::printf("  croissance par pas, mode de la grille proche de Nyquist (mesure : moyenne geometrique sur 50 pas, apres 200 pas) :\n    cas           C      formule      mesure\n");
    for (double cfl : {1.02, 1.2}) {  // 1D, 200 cases, mode n = 190
        Wave1DParams p;
        p.cells = 200;
        p.cfl = cfl;
        Wave1D w(p);
        const double k = 190.0 * dsp::kPi / p.length;
        w.setInitial([&](double x) { return std::sin(k * x); });
        std::vector<double> history;
        for (int i = 0; i < 250; ++i) {
            w.step();
            history.push_back(w.maxAbs());
        }
        const double g = geometricGrowth(history, 50), predicted = dsp::growthPerStep1D(k * w.dx(), cfl);
        CHECK(std::abs(g / predicted - 1.0) <= 1e-6);
        std::printf("    1D n = 190   %.2f   %.6f   %.6f\n", cfl, predicted, g);
    }
    {
        // 2D : C = 0.78 > 1/sqrt(2) ; mode (59, 58) sur 60 cases ; et C = 0.70 (stable) : aucune croissance
        for (double cfl : {0.78, 0.70}) {
            Wave2DParams p;
            p.cellsX = p.cellsY = 60;
            p.dx = 0.05;
            p.cfl = cfl;
            Wave2D w(p);
            const double l = p.cellsX * p.dx;
            const double kx = 59.0 * dsp::kPi / l, ky = 58.0 * dsp::kPi / l;
            w.setInitial([&](double x, double y) { return std::sin(kx * x) * std::sin(ky * y); });
            std::vector<double> history;
            const double e0 = w.energy();
            for (int i = 0; i < 250; ++i) {
                w.step();
                history.push_back(w.maxAbs());
            }
            const double predicted = dsp::growthPerStep2D(kx * p.dx, ky * p.dx, cfl);
            if (cfl > dsp::cflLimit(2)) {
                const double g = geometricGrowth(history, 50);
                CHECK(std::abs(g / predicted - 1.0) <= 1e-6);
                std::printf("    2D (59, 58)  %.2f   %.6f   %.6f\n", cfl, predicted, g);
            } else {
                // stable : le mode oscille sans croître. Le maximum INSTANTANÉ oscille (maxima sur 50 pas : 1,0005 d'écart, trop grossier) ;
                // la bonne grandeur est l'énergie, conservée exactement.
                const double drift = std::abs(w.energy() / e0 - 1.0);
                CHECK(drift <= 1e-12 && predicted == 1.0);
                std::printf("    2D (59, 58)  %.2f   %.6f   (energie finale / initiale - 1 = %.1e)\n", cfl, predicted, drift);
            }
        }
    }
}

// L'erreur de dispersion d'une impulsion est proportionnelle à (1 - C²) : plus faible quand C approche de 1, nulle à C = 1.
void testDispersionErrorConstant() {
    const auto f = [](double x) { return gaussian(x, 2.0, 0.1); };
    auto run = [&](double cfl) {
        Wave1DParams p;
        p.cells = 400;
        p.length = 4.0;
        p.cfl = cfl;
        Wave1D w(p);
        w.setInitial(f);
        w.advance(std::lround(0.8 / w.dt()));
        return maxError1D(w, [&](double x, double t) { return dAlembert(f, x, p.speed * t); });
    };
    const double e02 = run(0.2), e05 = run(0.5), e09 = run(0.9), e099 = run(0.99), e1 = run(1.0);
    // rapport à C = 0.5 prévu par le terme dominant : (1 - C²) / (1 - 0.25)
    const double r02 = e02 / e05, r09 = e09 / e05, r099 = e099 / e05;
    CHECK(std::abs(r02 / (0.96 / 0.75) - 1.0) <= 0.1);
    CHECK(std::abs(r09 / (0.19 / 0.75) - 1.0) <= 0.1);
    CHECK(std::abs(r099 / (0.0199 / 0.75) - 1.0) <= 0.15);
    CHECK(e1 <= 1e-13);
    std::printf("  erreur d'une impulsion (400 cases, t = 0.8) : C = 0.2 : %.3e ; 0.5 : %.3e ; 0.9 : %.3e ; 0.99 : %.3e ; 1.0 : %.1e\n"
                "  rapports a C = 0.5 : %.3f (prevu %.3f), %.3f (%.3f), %.3f (%.3f)\n",
                e02, e05, e09, e099, e1, r02, 0.96 / 0.75, r09, 0.19 / 0.75, r099, 0.0199 / 0.75);
}

}  // namespace

int main() {
    std::printf("FFT\n");
    testFft();
    std::printf("Ondes 1D\n");
    testWave1DEnergy();
    testWave1DModes();
    testWave1DAlembert();
    testWave1DReflection();
    testWave1DSponge();
    testWave1DStability();
    std::printf("Ondes 2D\n");
    testWave2DEnergy();
    testWave2DMode();
    testWave2DSymmetry();
    testWave2DSponge();
    testWave2DStability();
    std::printf("Corde (O1)\n");
    testStringFormulas();
    testStringModes();
    testStringExact();
    testStringContinuum();
    testStringHarmonics();
    testStringSolvers();
    testStringStability();
    testStringIsLeapfrog();
    testStringFrequencies();
    std::printf("Dispersion numerique et CFL (O2)\n");
    testDispersionFormulas();
    testSetStatesShift();
    testGroupVelocity();
    testPhaseVelocityFromSpectrum();
    testAnisotropy2D();
    testGrowthRates();
    testDispersionErrorConstant();

    if (g_failures == 0) {
        std::printf("test_waves : OK\n");
        return 0;
    }
    std::fprintf(stderr, "test_waves : %d echec(s)\n", g_failures);
    return 1;
}
