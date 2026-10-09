// Validation du domaine Ondes (O0) : FFT contre la transformée directe, équation d'onde 1D et 2D (solutions exactes,
// invariants, ordre de convergence, bords). Même esprit que test_core.cpp : pas de framework, un CHECK minimal ; les valeurs
// mesurées sont affichées pour qu'on puisse les relire (ne jamais se contenter d'un test qui passe).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/waves/Fft.hpp"
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

    if (g_failures == 0) {
        std::printf("test_waves : OK\n");
        return 0;
    }
    std::fprintf(stderr, "test_waves : %d echec(s)\n", g_failures);
    return 1;
}
