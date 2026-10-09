// Détection d'événements pendant l'intégration d'une EDO : instant exact où une grandeur g(t, y) s'annule
// (contact avec le sol, vitesse nulle d'un bloc qui frotte...).
//
// Pourquoi : la plupart des phénomènes de M5 (frottement sec, chocs) ont une dynamique DISCONTINUE : la force ou la vitesse
// change brutalement quand g(t, y) change de signe. Un pas fixe franchit la discontinuité au milieu d'un pas : l'erreur est alors
// d'ordre 1 quel que soit le schéma (RK4 compris). La parade est de s'arrêter pile sur l'événement, de traiter le changement
// de régime (inversion de la vitesse, adhérence...), puis de repartir.
//
// Méthode : on avance d'un pas complet ; si g change de signe, on cherche l'instant du changement par bissection sur la durée du
// pas, en repartant chaque fois de l'état initial du pas (le solveur n'a pas besoin de « sortie dense »). 100 itérations ramènent
// la durée à la précision de la machine. Coût négligeable pour de petits systèmes.
#pragma once

#include <functional>

#include "physicslab/core/Solver.hpp"

namespace pl {

// g(t, y) : un événement a lieu quand elle change de signe (valeur exactement nulle au départ = pas d'événement).
using EventFunction = std::function<double(double t, const State& y)>;

struct EventStep {
    double elapsed = 0.0;   // durée réellement avancée (<= dt)
    bool event = false;     // vrai si g a changé de signe : y est alors l'état à l'instant de l'événement
};

// Avance y de t à t + dt (en plusieurs pas internes si le solveur est adaptatif), mais s'arrête à l'instant où g s'annule,
// s'il a lieu pendant ce pas. `y` contient alors l'état à cet instant, au plus près (g(y) est de l'ordre de l'erreur du solveur).
EventStep advanceToEvent(Solver& solver, const OdeFunction& f, const EventFunction& g, double t, State& y, double dt);

}  // namespace pl
