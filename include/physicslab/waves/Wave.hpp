// O0 : équation d'onde u_tt = c² Δu sur une grille régulière, schéma « saute-mouton » (leapfrog, ordre 2 en espace et en temps).
//
//   u^{n+1} = 2 u^n - u^{n-1} + C² (Laplacien discret de u^n),  C = c dt / dx  (nombre de Courant, « CFL »).
//
// C'est le Verlet des orbites appliqué à chaque case de la grille. Stable si C <= 1 (1D) ou C <= 1/sqrt(2) (2D).
// Bords : mur fixe (u = 0, l'onde revient inversée), bout libre (du/dn = 0, l'onde revient de même signe) ou couche
// absorbante (« éponge » : amortissement progressif u_tt + 2 sigma u_t = c² Δu devant un mur fixe).
// Référence en `double` sur CPU ; le GPU `float` viendra avec la cuve à ondes 2D (O3).
#pragma once

#include <functional>
#include <vector>

namespace pl::waves {

enum class Edge { Fixed, Free, Absorbing };

// Amortissement maximal de l'éponge pour qu'une onde de haute fréquence qui traverse la couche (épaisseur `depth`), se réfléchit
// sur le mur et ressort soit réduite à la fraction `reflection` de son amplitude (profil sigma = sigma_max s^order, s = profondeur
// relative) : exp(-2 sigma_max depth / ((order + 1) c)) = reflection. ATTENTION : ce n'est valable que pour les ondes courtes devant
// l'épaisseur de la couche. Pour une impulsion à grandes longueurs d'onde (une bosse d'eau, de la largeur de la couche ou plus), la
// zone devient « suramortie » (u_t domine, elle se comporte comme un mur) : mesuré en test_waves.cpp et PASSATION, une éponge
// PLUS DURE réfléchit PLUS, et un profil plus lisse (order 2, 3) n'aide pas. Les défauts (R = 0,03, order 1) sortent d'un balayage.
double spongeSigmaMax(double speed, double depth, double reflection, int order);

// ------------------------------------------------------------------ 1D ------

struct Wave1DParams {
    int cells = 400;              // nombre d'intervalles : cells + 1 points x_i = i dx
    double length = 2.0;          // [m]
    double speed = 1.0;           // c [m/s]
    double cfl = 0.9;             // C = c dt / dx
    Edge left = Edge::Fixed;
    Edge right = Edge::Fixed;
    int spongeCells = 40;         // épaisseur de la couche absorbante, en intervalles
    double spongeReflection = 0.03;
    int spongeOrder = 1;          // sigma = sigma_max s^order
};

class Wave1D {
public:
    explicit Wave1D(const Wave1DParams& params);

    const Wave1DParams& params() const { return p_; }
    double dx() const { return dx_; }
    double dt() const { return dt_; }
    int points() const { return p_.cells + 1; }
    double x(int i) const { return i * dx_; }
    long long steps() const { return steps_; }
    double time() const { return steps_ * dt_; }
    bool stable() const { return p_.cfl <= 1.0 + 1e-12; }

    // u^n (instant courant) et u^{n-1}.
    const std::vector<double>& u() const { return cur_; }
    const std::vector<double>& uPrevious() const { return prev_; }
    // Amortissement sigma_i de l'éponge au point i (0 hors des couches absorbantes).
    double sigma(int i) const { return sigma_[i]; }

    // Position u0(x) et vitesse v0(x) à t = 0 (v0 vide = au repos). Les extrémités fixes ou absorbantes sont remises à 0.
    // L'état d'avant, u^{-1}, vient d'un développement de Taylor à l'ordre 2 : le schéma démarre sans erreur d'ordre 1.
    void setInitial(const std::function<double(double)>& u0, const std::function<double(double)>& v0 = {});
    // Impose directement les deux niveaux de temps u^{-1} (previous) et u^0 (current) : pour un état COHÉRENT AVEC LA GRILLE (par exemple un paquet
    // d'ondes qui va purement vers la droite, avec sa vraie relation de dispersion). Avec setInitial, la vitesse initiale est prise au sens
    // continu et un paquet de quelques cases par longueur d'onde reçoit une onde de retour parasite. Bords fixes ou absorbants remis à 0.
    void setStates(const std::function<double(double)>& previous, const std::function<double(double)>& current);

    void step();
    void advance(long long n) { for (long long k = 0; k < n; ++k) step(); }

    // Énergie ½ ∫ (u_t² + c² u_x²) dx (par unité de masse linéique) sous sa forme discrète CONSERVÉE exactement par le schéma
    // sans éponge : ½ Σ w_i dx ((u^n - u^{n-1})/dt)² + ½ c² Σ dx (D u^n)(D u^{n-1}), D = différence avant / dx, w_i = ½ aux extrémités.
    // Avec une éponge elle décroît (c'est voulu).
    double energy() const;
    double maxAbs() const;

private:
    Wave1DParams p_;
    double dx_, dt_;
    long long steps_ = 0;
    std::vector<double> prev_, cur_, next_, sigma_;
};

// ------------------------------------------------------------------ 2D ------

struct Wave2DParams {
    int cellsX = 200;
    int cellsY = 200;
    double dx = 0.01;             // côté d'une case [m] (cases carrées)
    double speed = 1.0;           // c [m/s]
    double cfl = 0.6;             // C = c dt / dx, stable si C <= 1/sqrt(2)
    Edge left = Edge::Fixed, right = Edge::Fixed, bottom = Edge::Fixed, top = Edge::Fixed;
    int spongeCells = 30;
    double spongeReflection = 0.03;
    int spongeOrder = 1;
};

class Wave2D {
public:
    explicit Wave2D(const Wave2DParams& params);

    const Wave2DParams& params() const { return p_; }
    double dx() const { return p_.dx; }
    double dt() const { return dt_; }
    int pointsX() const { return p_.cellsX + 1; }
    int pointsY() const { return p_.cellsY + 1; }
    int index(int i, int j) const { return j * (p_.cellsX + 1) + i; }
    long long steps() const { return steps_; }
    double time() const { return steps_ * dt_; }
    static constexpr double kCflLimit = 0.70710678118654752440;  // 1/sqrt(2)
    bool stable() const { return p_.cfl <= kCflLimit + 1e-12; }

    const std::vector<double>& u() const { return cur_; }
    const std::vector<double>& uPrevious() const { return prev_; }
    double sigma(int i, int j) const { return sigma_[index(i, j)]; }

    void setInitial(const std::function<double(double, double)>& u0, const std::function<double(double, double)>& v0 = {});

    void step();
    void advance(long long n) { for (long long k = 0; k < n; ++k) step(); }

    // Même énergie discrète conservée qu'en 1D (poids ½ sur les bords, ¼ aux coins).
    double energy() const;
    double maxAbs() const;

private:
    // Laplacien discret (sans 1/dx²) de f au point (i, j) : 5 points, avec point fantôme symétrique sur un bord libre.
    double laplacian(const std::vector<double>& f, int i, int j) const;
    void clampBoundaries(std::vector<double>& f) const;

    Wave2DParams p_;
    double dt_;
    long long steps_ = 0;
    std::vector<double> prev_, cur_, next_, sigma_;
};

}  // namespace pl::waves
