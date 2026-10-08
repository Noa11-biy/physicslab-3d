// M1 : chute libre et projectile avec frottement visqueux linéaire.
//
// Théorie (par particule de masse m, pesanteur g, frottement F = -b v, k = b/m) :
//     dv/dt = g - k v     =>     v(t) = v0 e^{-kt} + g phi(t)
//                                r(t) = r0 + v0 phi(t) + g psi(t)
//   avec phi(t) = (1 - e^{-kt}) / k   et   psi(t) = (t - phi(t)) / k.
// Quand k -> 0 on retrouve la parabole : phi -> t, psi -> t^2 / 2.
// Cette solution exacte sert de référence pour valider tous les solveurs.
#pragma once

#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"
#include "physicslab/core/World.hpp"

namespace pl {

struct ProjectileProblem {
    Vec3 position0{0.0, 2.0, 0.0};   // [m]
    Vec3 velocity0{8.0, 11.0, 0.0};  // [m/s]
    double mass = 1.0;               // [kg]
    Vec3 gravity{0.0, -constants::g0, 0.0};
    double linearDrag = 0.0;         // b [kg/s]

    // Monde prêt à être intégré, avec une particule dans l'état initial.
    World makeWorld() const;

    // Solution exacte.
    Vec3 position(double t) const;
    Vec3 velocity(double t) const;
    double energy(double t) const;  // Ec + Ep (Ep = -m g.r) évaluée sur la solution exacte

    // Premier instant où y(t) = 0 en retombant (0 si on part du sol sans monter ;
    // l'infini si la pesanteur ne rabaisse jamais la particule).
    double landingTime() const;
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie |r_numérique - r_exact| à l'arrivée.
// Sert aux tests d'ordre de convergence et à l'étude de convergence de l'interface.
double integrationError(const ProjectileProblem& problem, Solver& solver, int steps, double tEnd);

}  // namespace pl
