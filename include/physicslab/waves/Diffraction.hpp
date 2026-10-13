// O4 : fentes de Young et diffraction. Une onde plane (source en ligne) arrive sur un écran percé de n fentes de largeur b, de pas a ; on mesure l'amplitude
// derrière l'écran et on la compare à la théorie de FRAUNHOFER (champ lointain, Fresnel b²/(lambda R) << 1).
//
// Théorie. Le champ lointain est la transformée de Fourier de l'ouverture : A(theta) ~ F[ouverture](sin(theta) / lambda). Pour n fentes de largeur b et de pas a :
//     I(theta) / I(0) = sinc²(pi b s / lambda) * [ sin(n phi) / (n sin(phi)) ]²,   s = sin(theta),   phi = pi a s / lambda,   sinc(x) = sin(x)/x.
//   - le facteur de FENTE sinc² (diffraction par une fente) s'annule pour s = m lambda / b (m = ±1, ±2...) : premier minimum sin(theta) = lambda / b ;
//   - le facteur de RÉSEAU s'annule pour phi = m pi / n sauf aux multiples de pi (pics principaux s = m lambda / a). Deux fentes (n = 2) : cos²(phi), minima en
//     s = (m + 1/2) lambda / a, interfrange angulaire lambda / a, soit i = lambda D / a sur un écran à la distance D (petits angles).
// Mesure : la source est sinusoïdale ; on démodule le champ à la pulsation w sur une fenêtre de Hann (détection synchrone) : amplitude A(x, y) en chaque noeud, sans
// dépendre du nombre de pas par période. Les angles se lisent sur un ARC de cercle centré sur l'écran (le champ lointain y est à distance constante).
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "physicslab/waves/Dispersion.hpp"
#include "physicslab/waves/Wave.hpp"

namespace pl::waves::diffraction {

constexpr double kPi = 3.14159265358979323846;

// I(theta) / I(0) ci-dessus (n = 1 : fente simple ; b = largeur, a = pas entre centres, lambda = longueur d'onde). a et b dans la même unité que lambda.
double slitIntensity(double b, double a, int slits, double lambda, double sinTheta);

// Longueur d'onde RÉELLE sur la grille (axe x ou y) pour la pulsation omega : la grille dispersive donne k dx = 2 asin(sin(w dt / 2) / C) (voir Dispersion.hpp), donc
// lambda_grille = 2 pi dx / (k dx) un peu plus grand que 2 pi c / w. Les prédictions de la mesure l'utilisent.
double gridWavelength(double omega, double dt, double dx, double cfl);

// Champ lointain par FFT de l'ouverture (1 = ouvert, 0 = fermé, échantillonnée tous les dy) complétée par des zéros (padding fois). Renvoie I(sin theta) normalisée à
// son maximum pour |sin theta| <= 1, triée par sin theta croissant.
struct FarField {
    std::vector<double> sinTheta;
    std::vector<double> intensity;
};
FarField fraunhofer(const std::vector<double>& aperture, double dy, double lambda, int padding = 32);

// Champ à DISTANCE FINIE (proche comme lointain) derrière une ouverture, par l'intégrale de Rayleigh-Sommerfeld (deuxième espèce, écran rigide) en 2D :
//     u(P) = -(k / 2) * integrale_ouverture H0^(1)(k r) dy',   r = distance du point P de l'ouverture à P,   pour une onde incidente exp(i k x) d'amplitude 1
// (H0^(1) = J0 + i Y0 : fonction de Hankel ; le facteur vient de ∫ H0 dy' = 2 exp(i k x) / k sur tout un plan : sans écran on retrouve l'amplitude 1). La limite
// r -> infini redonne Fraunhofer. `aperture[j]` (1 ouvert, 0 fermé) est rangé tous les dy en y = j dy, la face arrière de l'écran est en x = xScreen ; k = 2 pi / lambda
// de la grille. Renvoie l'amplitude |u(px, py)| (4 points par case, règle du point milieu).
double apertureAmplitude(const std::vector<double>& aperture, double dy, double k, double xScreen, double px, double py);

// Même intégrale avec la phase de la GRILLE : le noyau de Hankel H0(k r) ~ sqrt(2 / (pi k r)) exp(i (k r - pi/4)) devient exp(i Phi(dx, dy)) / sqrt(r), Phi = RayPhase2D
// (phase de l'onde de la grille dont la vitesse de groupe va de l'élément d'ouverture au point), soit |u| = sqrt(k / (2 pi)) |somme exp(i Phi) / sqrt(r) dy'|. Valable
// pour k r >> 1 (au-delà de quelques longueurs d'onde de l'écran). C'est ce que le calcul sur grille donne vraiment ; l'écart avec apertureAmplitude (continu) est
// l'anisotropie de la grille, qui croît avec le nombre de longueurs d'onde parcourues.
double apertureAmplitudeGrid(const std::vector<double>& aperture, double dy, const dispersion::RayPhase2D& ray, double xScreen, double px, double py);

// Minima locaux STRICTS d'une courbe échantillonnée y(x) (x croissant), positions affinées par une parabole sur les trois points autour du minimum.
std::vector<double> findMinima(const std::vector<double>& x, const std::vector<double>& y);

// ---- l'expérience -------------------------------------------------------------------------------------------------------------------------------------------------

struct SlitScreenSpec {
    double lambda = 1.0;          // longueur d'onde dans le continu : w = 2 pi c / lambda [m]
    double dx = 0.1;              // côté d'une case [m]
    double cfl = 0.625;           // C = c dt / dx (limite 1/sqrt(2) = 0,7071)
    double lengthX = 40.0, lengthY = 45.0;  // domaine [m]
    double sourceX = 4.0;         // abscisse de la source en ligne
    double screenX = 10.0;        // abscisse de l'écran (face avant)
    double thickness = 0.2;       // épaisseur de l'écran [m] (au moins une case)
    int slits = 2;                // nombre de fentes (0 : écran plein, pour un test d'opacité)
    double slitWidth = 2.0;       // b [m], arrondie à un nombre entier de cases
    double pitch = 4.0;           // a [m], entre centres, arrondi à un nombre entier de cases
    double sponge = 3.0;          // épaisseur des couches absorbantes [m] (gauche, droite, haut, bas)
    double amplitude = 1.0;       // source en ligne
};

struct SlitScreen {
    Wave2DParams params;
    LineSource source;
    double omega = 0.0;
    int screenColumn = 0;         // premier noeud de l'écran en x
    int screenColumns = 1;        // nombre de colonnes de noeuds de l'écran
    double centerX = 0.0, centerY = 0.0;  // centre de l'écran (milieu de la face arrière) : centre des arcs de mesure
    double slitWidth = 0.0, pitch = 0.0;  // largeur et pas EFFECTIFS (nombres entiers de cases)
    std::vector<double> aperture;         // 1 = ouvert le long de l'écran, noeud par noeud en y (cellsY + 1 valeurs)
};

SlitScreen makeSlitScreen(const SlitScreenSpec& spec);

// Détection synchrone sur `steps` pas : somme pondérée (Hann) de u(t) exp(-i w t) ; amplitude(x, y) = 2 |somme| / somme des poids (l'amplitude d'une sinusoïde).
class LockIn {
public:
    LockIn(std::size_t nodes, double omega, double dt, long long steps);
    template <class T>
    void add(const std::vector<T>& u, long long step);  // step = numéro absolu du pas (le temps est step * dt), appelé une fois par pas
    long long count() const { return added_; }
    std::vector<double> amplitude() const;

private:
    double omega_, dt_;
    long long steps_, added_ = 0;
    double weights_ = 0.0;
    std::vector<double> re_, im_;
};

template <class T>
void LockIn::add(const std::vector<T>& u, long long step) {
    const double w = 0.5 * (1.0 - std::cos(2.0 * kPi * (static_cast<double>(added_) + 0.5) / static_cast<double>(steps_)));
    const double c = w * std::cos(omega_ * step * dt_), s = w * std::sin(omega_ * step * dt_);
    for (std::size_t k = 0; k < re_.size(); ++k) {
        re_[k] += c * static_cast<double>(u[k]);
        im_[k] += s * static_cast<double>(u[k]);
    }
    weights_ += w;
    ++added_;
}

// Valeur interpolée (bilinéaire) d'un champ rangé comme Wave2D (j * (cellsX + 1) + i) au point (x, y) [m].
double sampleField(const std::vector<double>& field, int cellsX, int cellsY, double dx, double x, double y);

// Amplitude le long d'un arc de rayon R centré en (cx, cy), pour les angles donnés (theta = 0 : vers +x).
std::vector<double> sampleArc(const std::vector<double>& field, int cellsX, int cellsY, double dx, double cx, double cy, double radius, const std::vector<double>& thetas);

}  // namespace pl::waves::diffraction
