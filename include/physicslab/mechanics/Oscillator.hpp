// M2 : oscillateur harmonique à un degré de liberté (ressort-masse), libre, amorti ou forcé.
//
// Théorie. Masse m, raideur k, frottement visqueux c, force extérieure F0 cos(w t) :
//     m x'' = -k x - c x' + F0 cos(w t)
//     x'' + 2 g x' + w0^2 x = a cos(w t)     avec  w0 = sqrt(k/m),  g = c / (2 m),  a = F0 / m,
//                                            zeta = g / w0.
// Solution = homogène (régime transitoire) + particulière (régime forcé permanent).
//   Homogène, selon zeta :
//     zeta < 1  (pseudo-périodique)  x_h = e^{-g t} (P cos wd t + Q sin wd t),   wd = sqrt(w0^2 - g^2)
//     zeta = 1  (critique)           x_h = e^{-g t} (P + (v_h0 + g P) t)
//     zeta > 1  (apériodique)        x_h = C1 e^{r1 t} + C2 e^{r2 t},            r = -g +- sqrt(g^2 - w0^2)
//   Particulière : x_p = A cos w t + B sin w t, avec  D = (w0^2 - w^2)^2 + (2 g w)^2,
//                  A = a (w0^2 - w^2) / D,  B = a (2 g w) / D   (amplitude X = a / sqrt(D)).
//   Résonance sans frottement (w = w0, D = 0) : x_p = a t sin(w0 t) / (2 w0), amplitude croissante.
// Les conditions initiales portent sur la somme : x_h(0) = x0 - x_p(0), v_h(0) = v0 - v_p(0).
#pragma once

#include "physicslab/core/Solver.hpp"

namespace pl {

enum class DampingRegime { Underdamped, Critical, Overdamped };

struct OscillatorProblem {
    double mass = 1.0;             // m [kg]
    double stiffness = 10.0;       // k [N/m]
    double damping = 0.0;          // c [kg/s]
    double forceAmplitude = 0.0;   // F0 [N]
    double forceFrequency = 0.0;   // w [rad/s]
    double x0 = 1.0;               // élongation initiale [m]
    double v0 = 0.0;               // vitesse initiale [m/s]

    // --- grandeurs caractéristiques ---
    double omega0() const;         // pulsation propre sqrt(k/m) [rad/s]
    double period() const;         // période propre 2 pi / w0 [s]
    double gamma() const;          // taux d'amortissement c / (2 m) [1/s]
    double zeta() const;           // taux d'amortissement réduit c / (2 sqrt(k m))
    DampingRegime regime() const;
    double dampedOmega() const;    // pseudo-pulsation sqrt(w0^2 - g^2) (0 si zeta >= 1)
    double steadyStateAmplitude(double omega) const;  // X(w) du régime forcé permanent
    double resonanceOmega() const;                    // w du maximum de X(w) (0 s'il n'y en a pas)

    // --- énergie ---
    double energy(double x, double v) const { return 0.5 * mass * v * v + 0.5 * stiffness * x * x; }

    // --- intégration numérique : état y = [x | v] (convention des solveurs, n = 1) ---
    State initialState() const { return {x0, v0}; }
    OdeFunction rhs() const;

    // --- solution exacte ---
    void exact(double t, double& x, double& v) const;
    double position(double t) const { double x, v; exact(t, x, v); return x; }
    double velocity(double t) const { double x, v; exact(t, x, v); return v; }
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie l'erreur finale dans l'espace des phases,
// sqrt(dx^2 + (dv / w0)^2). Pour mesurer un ordre de convergence, choisir tEnd qui n'est PAS un multiple de la
// période : à t = n T un oscillateur repasse par son état initial et l'erreur de phase d'ordre 1 s'annule
// (sur-convergence apparente : Verlet ordre 4 au lieu de 2, Euler symplectique ordre 2 au lieu de 1).
double oscillatorError(const OscillatorProblem& problem, Solver& solver, int steps, double tEnd);

}  // namespace pl
