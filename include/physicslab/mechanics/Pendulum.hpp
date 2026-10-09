// M3 : pendule simple (tige rigide sans masse de longueur L, masse m au bout).
//
// Théorie. Moment cinétique par rapport au pivot : m L^2 theta'' = -m g L sin(theta) - b theta'.
//     theta'' = -w0^2 sin(theta) - gamma theta'     avec  w0 = sqrt(g / L),  gamma = b / (m L^2).
// Énergie (sans frottement conservée) : E = 1/2 m L^2 theta'^2 + m g L (1 - cos theta).
//
// Petits angles : sin(theta) ~ theta  =>  oscillateur harmonique de période T0 = 2 pi / w0, indépendante de
// l'amplitude. Le vrai pendule est plus lent quand l'amplitude croît.
//
// Solution exacte (sans frottement, lâché sans vitesse depuis theta0 < pi). Avec k = sin(theta0 / 2) :
//     sin(theta / 2) = k sn(K + w0 t, k)       (sn : sinus elliptique de Jacobi, K = K(k))
//     T = 4 K(k) / w0
// sn(K + x) est paire et périodique : son développement de Fourier (DLMF 22.11.1) donne, avec la nome
// q = exp(-pi K' / K) et x = pi w0 t / (2 K),
//     k sn(K + w0 t) = (2 pi / K) Sum_{n>=0} (-1)^n q^{n+1/2} / (1 - q^{2n+1}) cos((2n+1) x)
// série qui converge géométriquement (raison q^2 < 1).
// Dans tous les autres cas (frottement, vitesse initiale, rotation) il n'y a pas de solution élémentaire :
// la référence est alors un RK45 à tolérance 1e-13.
#pragma once

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"

namespace pl {

// Intégrale elliptique complète de première espèce K(k) = Int_0^{pi/2} dphi / sqrt(1 - k^2 sin^2 phi),
// pour un module 0 <= k < 1 (moyenne arithmético-géométrique).
double ellipticK(double k);

struct PendulumProblem {
    double length = 1.0;            // L [m]
    double mass = 1.0;              // m [kg]
    double gravity = constants::g0; // g [m/s^2]
    double damping = 0.0;           // b [N.m.s] : couple de frottement -b theta'
    double theta0 = 1.0;            // angle initial par rapport à la verticale basse [rad]
    double angularVelocity0 = 0.0;  // vitesse angulaire initiale [rad/s]

    double omega0() const;          // pulsation des petits angles sqrt(g / L)
    double smallAnglePeriod() const;  // T0 = 2 pi / w0
    double gamma() const;           // b / (m L^2)

    // La solution exacte existe : pas de frottement, lâché sans vitesse, |theta0| < pi.
    bool hasExactSolution() const;
    // Période réelle 4 K(sin(theta0 / 2)) / w0 si hasExactSolution(), sinon T0.
    double period() const;

    double energy(double theta, double omega) const;  // origine : position basse au repos

    // y = [theta | omega] (convention des solveurs, n = 1)
    State initialState() const { return {theta0, angularVelocity0}; }
    OdeFunction rhs() const;

    // Solution exacte ; précondition hasExactSolution().
    void exact(double t, double& theta, double& omega) const;
    // Référence : exacte si possible, sinon RK45 haute précision.
    void reference(double t, double& theta, double& omega) const;
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie l'erreur finale dans l'espace des phases,
// sqrt(dtheta^2 + (domega / w0)^2), par rapport à la référence. Choisir tEnd non multiple de la période
// (voir oscillatorError : sur-convergence apparente sinon).
double pendulumError(const PendulumProblem& problem, Solver& solver, int steps, double tEnd);

}  // namespace pl
