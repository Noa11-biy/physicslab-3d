// M5b : chocs entre deux corps lisses (sans frottement), loi de restitution de Newton.
//
// Deux corps de masses m1, m2, vitesses v1, v2, et n le vecteur unitaire normal au contact, dirigé de 1 vers 2. La vitesse
// d'approche normale est  v_n = (v1 - v2) . n  (positive si les corps se rapprochent). Un choc dure un temps négligeable : seules
// les vitesses changent, par une IMPULSION J le long de n, égale et opposée sur les deux corps (3e loi de Newton) :
//     v1' = v1 - (J/m1) n,      v2' = v2 + (J/m2) n.
// Loi de restitution de Newton : la vitesse relative normale est inversée et réduite par le coefficient e, 0 <= e <= 1 :
//     (v1' - v2') . n = -e (v1 - v2) . n         =>     J = (1 + e) mu_r v_n,     mu_r = m1 m2 / (m1 + m2).
// Conséquences : l'impulsion totale m1 v1 + m2 v2 est conservée ; la composante tangentielle de chaque vitesse est inchangée ;
// l'énergie cinétique perdue vaut  ½ (1 - e^2) mu_r v_n^2  (e = 1 : choc élastique, rien de perdu ; e = 0 : choc « mou », les
// deux corps ont la même vitesse normale après). Deux billes de même masse, l'une au repos, en choc élastique : elles repartent à
// 90 degrés ; si e < 1, à moins de 90 degrés (le produit scalaire des vitesses finales vaut (1 - e^2)/4 de v_n^2 par unité de masse).
//
// Détection du contact : le contact a lieu à l'instant où la distance entre les surfaces s'annule. Un pas fixe ne le voit qu'APRÈS (les
// corps se sont interpénétrés d'une fraction de pas) : le choc est appliqué trop tard, la trajectoire suivante est décalée de
// O(dt) quel que soit le solveur. L'événement (core/Events.hpp) trouve l'instant exact. Voir ContactModel.
#pragma once

#include <limits>

#include "physicslab/core/Events.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

namespace collision {

struct Result1D {
    double v1 = 0.0, v2 = 0.0;   // vitesses après le choc
    double energyLoss = 0.0;     // énergie cinétique perdue (>= 0)
};

// Choc frontal de deux masses sur une droite orientée : le corps 1 est à gauche du corps 2, le choc a lieu si v1 > v2.
// Si v1 <= v2 les corps s'éloignent : rien ne change.
Result1D collide1D(double m1, double v1, double m2, double v2, double e);

struct Result3D {
    Vec3 v1, v2;
    double impulse = 0.0;        // J, le long de la normale
    double energyLoss = 0.0;
};

// Choc de deux sphères lisses, normale unitaire n dirigée de 1 vers 2. Si elles s'éloignent (v_n <= 0) : rien ne change.
Result3D collideSpheres(double m1, const Vec3& v1, double m2, const Vec3& v2, const Vec3& n, double e);

}  // namespace collision

// Comment l'instant du contact est trouvé (voir en-tête).
enum class ContactModel {
    Naive,         // pas fixe, contact détecté APRÈS le pas (interpénétration), choc appliqué sur l'état courant
    EventDriven    // arrêt exact à l'instant du contact (bissection), choc, puis reprise
};

// Deux billes (disques) en mouvement libre dans le plan xy, qui se heurtent au plus une fois. Convention des solveurs :
// y = [x1 y1 x2 y2 | vx1 vy1 vx2 vy2] (positions | vitesses, n = 4). Hors choc le mouvement est rectiligne uniforme : TOUS les solveurs
// sont exacts entre deux chocs, et toute l'erreur d'un schéma vient de l'instant du contact.
struct TwoBallProblem {
    double m1 = 1.0, m2 = 1.0;                  // masses [kg]
    double r1 = 0.5, r2 = 0.5;                  // rayons [m]
    Vec3 p1{-4.0, 0.3, 0.0}, p2{0.0, 0.0, 0.0}; // positions initiales [m]
    Vec3 v1{3.0, 0.0, 0.0}, v2{0.0, 0.0, 0.0};  // vitesses initiales [m/s]
    double restitution = 1.0;                   // e

    State initialState() const;
    OdeFunction rhs() const;                    // mouvement libre
    double gap(const State& y) const;           // distance entre les surfaces (<= 0 : contact ou interpénétration)
    double kineticEnergy(const State& y) const;
    Vec3 momentum(const State& y) const;

    // Instant du contact (l'infini si les billes ne se touchent jamais) et solution exacte.
    double collisionTime() const;
    State exact(double t) const;
    // Angle [rad] entre les vitesses finales des deux billes (exact, après le choc).
    double outgoingAngle() const;
};

// Pilote de TwoBallProblem avec l'un des deux modèles de détection.
class TwoBallRun {
public:
    TwoBallRun(const TwoBallProblem& problem, ContactModel model);

    void advance(Solver& solver, double dt);   // avance de dt (un seul pas pour le modèle naïf)

    double time() const { return time_; }
    const State& state() const { return y_; }
    int collisions() const { return collisions_; }
    double collisionTime() const { return firstCollision_; }  // instant détecté du premier choc (NaN s'il n'a pas eu lieu)

private:
    void collide();

    TwoBallProblem problem_;
    ContactModel model_;
    State y_;
    double time_ = 0.0;
    int collisions_ = 0;
    double firstCollision_ = std::numeric_limits<double>::quiet_NaN();
};

// Intègre jusqu'à tEnd en `steps` pas égaux ; renvoie la plus grande erreur de position des deux billes par rapport à la solution exacte.
double twoBallError(const TwoBallProblem& problem, ContactModel model, Solver& solver, int steps, double tEnd);

}  // namespace pl
