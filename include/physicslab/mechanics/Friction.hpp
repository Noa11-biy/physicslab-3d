// M5a : frottement sec de Coulomb, bloc glissant sur un plan incliné.
//
// Axe s le long de la pente, positif vers le HAUT (le bloc est lancé vers le haut avec v0 > 0). Par unité de masse, avec la pente
// theta, la réaction normale N/m = g cos(theta), le frottement dynamique mu_d, statique mu_s >= mu_d et une résistance linéaire
// k = b/m (optionnelle) :
//     glissement :  s'' = -g sin(theta) - mu_d g cos(theta) sgn(s') - k s'
//     à l'arrêt :   le frottement statique s'adapte pour retenir le bloc tant que  g sin(theta) <= mu_s g cos(theta),
//                   soit  tan(theta) <= mu_s  ; sinon le bloc repart en glissant vers le bas.
// La loi de Coulomb est DISCONTINUE : la force change de signe quand la vitesse s'annule (sgn) et l'adhérence est une inégalité.
//
// Dans une phase de glissement, le sens de la vitesse dir (+1 montée, -1 descente) est fixé : l'EDO devient linéaire
//     v' = A - k v,     A = -g (sin(theta) + dir mu_d cos(theta))
// et se résout exactement : v(t) = A/k + (v0 - A/k) e^{-kt},  s(t) = s0 + A t/k + (v0 - A/k)(1 - e^{-kt})/k
// (k -> 0 : mouvement uniformément accéléré). La phase s'arrête quand v s'annule, à t = ln(1 - k v0/A)/k (-v0/A si k = 0), possible
// seulement si A v0 < 0. Ensuite : adhérence définitive si tan(theta) <= mu_s, sinon glissement vers le bas (dir = -1).
// La solution exacte est donc DÉFINIE PAR MORCEAUX, au plus trois morceaux (montée, descente, repos).
//
// Bilan d'énergie sans résistance (k = 0) : E = v^2/2 + g s sin(theta) diminue du travail de frottement mu_d g cos(theta) x (chemin parcouru).
//
// Trois façons de l'intégrer numériquement (voir InclineModel) : naïve (sgn dans l'EDO), régularisée (sgn lissé) et « événement +
// adhérence » (arrêt exact à v = 0, voir Events.hpp). Convention des solveurs : y = [s | v] (positions | vitesses, n = 1).
#pragma once

#include <limits>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Events.hpp"
#include "physicslab/core/Solver.hpp"

namespace pl {

enum class InclineModel {
    Naive,         // s'' = ... - mu_d g cos(theta) sgn(s'), sgn(0) = 0 : ignore mu_s et ne s'arrête jamais exactement
    Regularized,   // sgn(v) remplacé par tanh(v / eps) : continu mais raide, le bloc « rampe » au lieu d'adhérer
    EventDriven    // arrêt exact à v = 0 (bissection) puis adhérence ou changement de sens : modèle correct
};

struct InclineState {
    double s = 0.0, v = 0.0;
    double path = 0.0;     // chemin parcouru (intégrale de |v|) [m]
    bool stuck = false;    // adhérence (au repos définitivement)
};

struct InclineProblem {
    double angle = 0.5;            // pente theta [rad], 0 <= theta < pi/2
    double gravity = constants::g0;
    double muStatic = 0.5;         // mu_s >= mu_d
    double muKinetic = 0.4;        // mu_d
    double drag = 0.0;             // k = b/m [1/s]
    double s0 = 0.0, v0 = 6.0;     // position et vitesse initiales le long de la pente (+ = vers le haut) [m, m/s]

    State initialState() const { return {s0, v0}; }
    bool holds() const;                       // tan(theta) <= mu_s : le bloc au repos reste collé
    double stageAcceleration(int dir) const;  // A = -g (sin(theta) + dir mu_d cos(theta))
    double regularization = 0.05;             // eps [m/s] de tanh(v/eps) pour InclineModel::Regularized

    OdeFunction rhs(InclineModel model) const;     // Naive ou Regularized (EventDriven : voir rhsStage)
    OdeFunction rhsStage(int dir) const;           // frottement de sens fixé : EDO linéaire d'une phase

    InclineState exact(double t) const;
    double energy(double s, double v) const;       // v^2/2 + g s sin(theta) par unité de masse [J/kg]
    // Durée jusqu'au premier arrêt (v = 0) si le bloc s'arrête ; l'infini sinon.
    double firstStopTime() const;
};

// Intègre le bloc avec « événement + adhérence » : le pas est interrompu pile à v = 0, puis on repart (ou on colle).
class InclineRun {
public:
    explicit InclineRun(const InclineProblem& problem);

    void advance(Solver& solver, double dt);   // avance de dt exactement

    double time() const { return time_; }
    const InclineState& state() const { return state_; }
    int stops() const { return stops_; }       // nombre d'arrêts détectés (0, 1 ou 2)
    double stopTime() const { return firstStop_; }  // instant du premier arrêt (NaN s'il n'a pas eu lieu)

private:
    InclineProblem problem_;
    InclineState state_;
    double time_ = 0.0;
    int dir_ = 1;
    int stops_ = 0;
    double firstStop_ = std::numeric_limits<double>::quiet_NaN();
};

// Intègre jusqu'à tEnd en `steps` pas égaux avec le modèle choisi ; renvoie |s_numérique - s_exact| à l'arrivée.
double inclineError(const InclineProblem& problem, InclineModel model, Solver& solver, int steps, double tEnd);

}  // namespace pl
