// Monde physique : un ensemble de particules ponctuelles dans un champ de pesanteur uniforme.
// Les modules suivants y brancheront d'autres forces ; l'interface Solver reste la même.
#pragma once

#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

struct Particle {
    Vec3 position;
    Vec3 velocity;
    double mass = 1.0;  // [kg]
};

// Grandeurs censées se conserver (ou dont on mesure la dérive) : validation en direct.
struct Invariants {
    double kinetic = 0.0;     // énergie cinétique [J]
    double potential = 0.0;   // énergie potentielle de pesanteur [J]
    Vec3 momentum;            // quantité de mouvement [kg.m/s]
    Vec3 angularMomentum;     // moment cinétique par rapport à l'origine [kg.m^2/s]

    double total() const { return kinetic + potential; }
};

class World {
public:
    double time = 0.0;                                  // [s]
    Vec3 gravity{0.0, -constants::g0, 0.0};             // accélération de pesanteur [m/s^2]
    std::vector<Particle> particles;

    // Avance le monde d'un pas avec le solveur donné ; renvoie la durée avancée.
    double step(Solver& solver, double dt);

    Invariants invariants() const;

private:
    State state_;  // [positions (3N) | vitesses (3N)], tampon réutilisé
};

}  // namespace pl
