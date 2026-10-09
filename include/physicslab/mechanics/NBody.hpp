// M4b : gravitation à N corps (CPU, double). Référence de validation du futur calcul GPU (M7).
//
// N masses ponctuelles m_i en interaction newtonienne. L'accélération du corps i (adoucissement de Plummer, eps >= 0) :
//     a_i = G sum_{j != i} m_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
// dérive de l'énergie potentielle (adoucie) U = -G sum_{i<j} m_i m_j / sqrt(|r_i - r_j|^2 + eps^2) : avec eps > 0, l'énergie
// totale E = T + U reste exactement conservée (sans adoucissement, une collision rapprochée donnerait une force infinie).
// La force sur i par j est exactement l'opposée de celle sur j par i (3e loi de Newton), d'où les invariants :
//     impulsion  P = sum m_i v_i      (symétrie par translation),   centre de masse en mouvement rectiligne uniforme
//     moment cinétique  L = sum m_i r_i x v_i   (symétrie par rotation)
//     énergie  E                      (symétrie par translation dans le temps)
// Le calcul est en O(N^2) : on boucle sur les paires i < j et chaque force est appliquée aux deux corps (impulsion
// conservée à l'arrondi près par construction).
//
// Cas de référence, avec solution connue :
//   * 2 corps : le mouvement relatif est le problème de Kepler de GM = G (m1 + m2) (voir Kepler.hpp).
//   * triangle de Lagrange : 3 masses égales aux sommets d'un triangle équilatéral de côté s tournent rigidement autour de
//     leur centre de masse avec  omega^2 = 3 G m / s^3  (force centripète = résultante des deux attractions : sqrt(3) G m / s^2).
//   * « huit » de Chenciner-Montgomery (2000) : 3 masses égales suivent la même courbe en huit, décalées d'un tiers de période.
//     Conditions initiales numériques de Simó : période T ~ 6,32591398 (G = m = 1), E ~ -1,28714, P = L = 0.
//
// Convention des solveurs : y = [x0 y0 z0 x1 ... | vx0 vy0 vz0 vx1 ...] (positions | vitesses, n = 3N) : le tableau d'accélérations
// est exactement la seconde moitié de f(t, y).
#pragma once

#include <vector>

#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

namespace nbody {

// « Force sur chaque particule », découplée de l'intégrateur (interface reprise par le calcul GPU).
// positions : 3N doubles (x0 y0 z0 x1 ...), masses : N, acc : 3N doubles écrasés.
void accelerations(const double* positions, const double* masses, int n, double G, double softening, double* acc);

// Énergie potentielle (adoucie) du système.
double potentialEnergy(const double* positions, const double* masses, int n, double G, double softening);

inline constexpr double kFigureEightPeriod = 6.32591398;  // G = m = 1

}  // namespace nbody

struct NBodyProblem {
    double G = 1.0;                  // constante gravitationnelle (1 = unités normalisées)
    double softening = 0.0;          // longueur d'adoucissement eps
    std::vector<double> mass;        // N masses
    std::vector<Vec3> position, velocity;  // conditions initiales

    int count() const { return static_cast<int>(mass.size()); }
    State initialState() const;
    OdeFunction rhs() const;

    // Grandeurs conservées d'un état quelconque.
    double kineticEnergy(const State& y) const;
    double potentialEnergy(const State& y) const;
    double energy(const State& y) const { return kineticEnergy(y) + potentialEnergy(y); }
    Vec3 momentum(const State& y) const;
    Vec3 angularMomentum(const State& y) const;   // par rapport à l'origine
    Vec3 centerOfMass(const State& y) const;
    double totalMass() const;

    // État à l'instant t, par RK45 haute précision (fiable tant que le chaos n'a pas amplifié l'erreur).
    State reference(double t) const;

    // Distance dans l'espace des phases : positions en unités de `length`, vitesses en unités de sqrt(G M / length)
    // (length = distance moyenne au centre de masse à t = 0, M = masse totale).
    double distance(const State& a, const State& b) const;

    // --- cas de référence ---
    static NBodyProblem twoBody(double m1, double m2, double semiMajor, double eccentricity, double G = 1.0);  // périastre, repère du centre de masse
    static NBodyProblem lagrangeTriangle(double side = 1.0, double mass = 1.0, double G = 1.0);
    static NBodyProblem figureEight();
    // N corps de masse totale 1 dans une boule de rayon `radius`, vitesses aléatoires à l'équilibre du viriel (2T = -U),
    // centre de masse au repos à l'origine. Déterministe pour un `seed` donné.
    static NBodyProblem randomCluster(int n, unsigned seed, double radius = 1.0, double softening = 0.05);
    // Sphère de Plummer à l'équilibre (M7) : N corps égaux, masse totale 1, G = 1, rayon d'échelle `scaleRadius` a. Densité
    // rho(r) ~ (1 + r²/a²)^(-5/2) ; positions et vitesses tirées par la méthode d'Aarseth, Hénon et Wielen (1974), qui échantillonne
    // la fonction de distribution exacte (équilibre sans mise à l'échelle). Rayons limités à 10 a (1,5 % de la masse écartée).
    // Centre de masse immobile à l'origine. Théorie (N grand, sans adoucissement) : E = -3 pi / 64 G M² / a, 2T = -U, rayon de
    // demi-masse 1,3048 a. Déterministe pour un `seed` donné, et en O(N) (contrairement à randomCluster).
    static NBodyProblem plummer(int n, unsigned seed, double scaleRadius = 1.0, double softening = 0.05);
    // Deux sphères de Plummer de n/2 corps et de masse 1/2 chacune, centres en (±separation/2, 0, 0) qui se rapprochent à la vitesse
    // relative `relativeSpeed` (repère du centre de masse). Chaque sphère garde l'équilibre qu'elle a seule.
    static NBodyProblem plummerCollision(int n, unsigned seed, double separation = 6.0, double relativeSpeed = 0.4,
                                         double scaleRadius = 0.6, double softening = 0.05);
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie la distance finale à la référence RK45.
double nbodyError(const NBodyProblem& problem, Solver& solver, int steps, double tEnd);

}  // namespace pl
