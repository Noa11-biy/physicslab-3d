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
    // Vitesse c(i, j) en chaque noeud (taille (cellsX + 1)(cellsY + 1), rangée comme les champs) : vide = vitesse uniforme `speed`. Un milieu
    // à vitesse variable est un milieu d'indice n = c0 / c(x, y) (cuve à ondes de profondeur variable) : réfraction de Snell-Descartes.
    std::vector<double> speedField;
};

// Source ponctuelle sinusoïdale : u_tt = ... + A sin(w t) ramp(t) delta(x - xs), « douce » (montée en sin² sur `rampPeriods` périodes) pour ne pas
// lancer de transitoire large bande. Numériquement : u^{n+1}[noeud] += (dt² / dx²) A sin(w t^n) ramp.
struct PointSource {
    double x = 0.0, y = 0.0;     // position [m] (arrondie au noeud le plus proche)
    double amplitude = 1.0;      // A
    double omega = 6.283185307179586;  // pulsation [rad/s]
    double rampPeriods = 2.0;
};

class Wave2D {
public:
    explicit Wave2D(const Wave2DParams& params);

    // Vitesse c(x, y) en chaque noeud à partir d'une fonction (taille (cellsX + 1)(cellsY + 1), à passer dans Wave2DParams::speedField).
    static std::vector<double> makeSpeedField(const Wave2DParams& geometry, const std::function<double(double, double)>& c);

    const Wave2DParams& params() const { return p_; }
    double dx() const { return p_.dx; }
    double dt() const { return dt_; }
    int pointsX() const { return p_.cellsX + 1; }
    int pointsY() const { return p_.cellsY + 1; }
    int index(int i, int j) const { return j * (p_.cellsX + 1) + i; }
    long long steps() const { return steps_; }
    double time() const { return steps_ * dt_; }
    static constexpr double kCflLimit = 0.70710678118654752440;  // 1/sqrt(2)
    // Le nombre de Courant est rapporté à la vitesse MAXIMALE du milieu (c'est elle qui fixe dt).
    double speedMax() const { return cMax_; }
    bool stable() const { return p_.cfl <= kCflLimit + 1e-12; }

    const std::vector<double>& u() const { return cur_; }
    const std::vector<double>& uPrevious() const { return prev_; }
    double sigma(int i, int j) const { return sigma_[index(i, j)]; }
    // Coefficients du schéma, repris tels quels par le calcul GPU : kx[k] = (dt² / dx²) c²_face entre les noeuds k et k + 1 (de même ky entre k et
    // k + (cellsX + 1)), avec c²_face = moyenne arithmétique des c² des deux noeuds ; sigmaField() : amortissement de l'éponge en chaque noeud.
    const std::vector<double>& coefficientX() const { return kx_; }
    const std::vector<double>& coefficientY() const { return ky_; }
    const std::vector<double>& sigmaField() const { return sigma_; }

    void setInitial(const std::function<double(double, double)>& u0, const std::function<double(double, double)>& v0 = {});
    // Deux niveaux de temps imposés (voir Wave1D::setStates) : un paquet purement progressif cohérent avec la grille.
    void setStates(const std::function<double(double, double)>& previous, const std::function<double(double, double)>& current);

    void clearSources() { sources_.clear(); }
    void addSource(const PointSource& s) { sources_.push_back(s); }
    const std::vector<PointSource>& sources() const { return sources_; }

    void step();
    void advance(long long n) { for (long long k = 0; k < n; ++k) step(); }

    // Même énergie discrète conservée qu'en 1D (poids ½ sur les bords, ¼ aux coins), pour l'équation u_tt = div(c² grad u) : le potentiel est
    // ½ sum c²_face (Du^n)(Du^{n-1}) avec la même c²_face que le schéma, donc conservée exactement aussi quand c varie (sans éponge ni source).
    double energy() const;
    double maxAbs() const;

private:
    // Divergence discrète (c² grad u) dt² / dx² au noeud (i, j), avec point fantôme symétrique sur un bord libre : remplace C² * laplacien.
    double weightedLaplacian(const std::vector<double>& f, int i, int j) const;
    void clampBoundaries(std::vector<double>& f) const;

    Wave2DParams p_;
    double dt_ = 0.0, cMax_ = 0.0, ratio_ = 0.0;  // ratio_ = dt² / dx²
    long long steps_ = 0;
    std::vector<double> prev_, cur_, next_, sigma_, kx_, ky_;
    std::vector<PointSource> sources_;
};

}  // namespace pl::waves
