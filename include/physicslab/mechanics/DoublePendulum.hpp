// M3 (suite) : pendule double, le système chaotique de référence.
//
// Deux tiges sans masse (longueurs l1, l2), masses m1 et m2, angles theta1, theta2 mesurés depuis la verticale
// basse. Lagrangien L = T - V avec
//     T = 1/2 (m1 + m2) l1^2 w1^2 + 1/2 m2 l2^2 w2^2 + m2 l1 l2 w1 w2 cos(theta1 - theta2)
//     V = -(m1 + m2) g l1 cos(theta1) - m2 g l2 cos(theta2)
// Les équations d'Euler-Lagrange, résolues en w1' et w2', donnent (delta = theta1 - theta2,
// D = 2 m1 + m2 - m2 cos(2 delta)) :
//     w1' = [-g (2 m1 + m2) sin(theta1) - m2 g sin(theta1 - 2 theta2) - 2 sin(delta) m2 (w2^2 l2 + w1^2 l1 cos(delta))]
//           / (l1 D)
//     w2' = [2 sin(delta) (w1^2 l1 (m1 + m2) + g (m1 + m2) cos(theta1) + w2^2 l2 m2 cos(delta))] / (l2 D)
// E = T + V est conservée. Pas de solution analytique : la référence est un RK45 à tolérance 1e-13, fiable
// tant que l'erreur n'a pas été amplifiée par le chaos (temps de Lyapunov, voir lyapunovExponent).
//
// Convention des solveurs : y = [theta1, theta2 | w1, w2] (positions | vitesses, n = 2).
// Le hamiltonien n'est pas séparable (la cinétique dépend des angles) : Euler symplectique et Verlet restent
// d'ordre 1 et 2 mais ne sont plus symplectiques ici.
#pragma once

#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"

namespace pl {

struct DoublePendulumProblem {
    double m1 = 1.0, m2 = 1.0;      // [kg]
    double l1 = 1.0, l2 = 1.0;      // [m]
    double gravity = constants::g0;
    double theta1 = 2.0, theta2 = 2.0;   // angles initiaux [rad]
    double omega1 = 0.0, omega2 = 0.0;   // vitesses angulaires initiales [rad/s]

    State initialState() const { return {theta1, theta2, omega1, omega2}; }
    OdeFunction rhs() const;
    double energy(const State& y) const;

    // Positions des deux masses, pivot à l'origine, y vers le haut.
    void positions(const State& y, double& x1, double& y1, double& x2, double& y2) const;

    // État à l'instant t, par RK45 haute précision.
    State reference(double t) const;

    // Distance dans l'espace des phases entre deux états : angles en rad, vitesses divisées par sqrt(g / l1).
    double distance(const State& a, const State& b) const;
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie la distance finale à la référence.
// Valable à tEnd modéré : le chaos amplifie ensuite toute erreur exponentiellement.
double doublePendulumError(const DoublePendulumProblem& problem, Solver& solver, int steps, double tEnd);

// Exposant de Lyapunov : ajuste ln d(t) = ln d0 + lambda t par moindres carrés sur les points où dMin <= d <= dMax
// (phase de croissance exponentielle, avant saturation). Renvoie lambda [1/s], ou NaN s'il y a moins de 8 points.
// `used` reçoit le nombre de points retenus, `intercept` l'ordonnée à l'origine ln d0 de la droite ajustée.
double lyapunovExponent(const std::vector<double>& t, const std::vector<double>& d, double dMin, double dMax,
                        int* used = nullptr, double* intercept = nullptr);

}  // namespace pl
