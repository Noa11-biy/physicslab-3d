// M5b : balle rebondissante sur un sol plat (y = 0), restitution de Newton e, résistance linéaire k = b/m optionnelle.
//
// Pendant un vol (plan xy, y vers le haut, pesanteur g) : x'' = -k x',  y'' = -g - k y'  : solution exacte de M1 (projectile avec frottement
// linéaire). À l'impact (y = 0, vy < 0) : vy -> -e vy, vx inchangée (sol lisse). Sans résistance (k = 0) et lâchée sans vitesse de la hauteur
// h0, avec t0 = sqrt(2 h0 / g) la durée de la première chute :
//     vitesse après le n-ième rebond  v_n = e^n g t0,     hauteur  h_n = e^(2n) h0,     durée du n-ième vol  2 e^n t0.
// La somme de ces durées est FINIE : la balle s'arrête de rebondir au bout de  t0 (1 + e) / (1 - e)  (suite géométrique : une infinité de
// rebonds en temps fini, « accumulation de Zénon »). Un calcul qui ne le sait pas continue de vibrer sur le sol. Ici la balle est déclarée
// au repos dès que sa vitesse de rebond tombe sous `restSpeed` : elle glisse alors sans frottement sur le sol (la composante horizontale
// continue, freinée par la résistance de l'air), comme l'EDO du vol le prévoit à y = 0, vy = 0.
//
// Convention des solveurs : y = [x, y | vx, vy] (positions | vitesses, n = 2).
#pragma once

#include <limits>
#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/mechanics/Collision.hpp"

namespace pl {

struct BounceState {
    double x = 0.0, y = 0.0, vx = 0.0, vy = 0.0;
    int bounces = 0;          // nombre de rebonds déjà effectués
    bool resting = false;     // la balle a cessé de rebondir (elle glisse sur le sol)
};

struct BounceProblem {
    double restitution = 0.8;                  // e, 0 <= e <= 1
    double gravity = constants::g0;
    double drag = 0.0;                         // k = b/m [1/s]
    double x0 = 0.0, y0 = 2.0;                 // position initiale [m] (y0 >= 0)
    double vx0 = 0.0, vy0 = 0.0;               // vitesse initiale [m/s]
    // Vitesse de rebond en dessous de laquelle la balle est au repos [m/s]. En pratique >= 1e-8 : en dessous la durée du vol
    // (2 v / g) n'est plus distinguable de 0 pour la recherche du contact.
    double restSpeed = 1e-4;

    State initialState() const { return {x0, y0, vx0, vy0}; }
    OdeFunction rhsFlight() const;             // vol : pesanteur et résistance
    OdeFunction rhsRest() const;               // au sol : seule la composante horizontale évolue (x' = vx, vx' = -k vx)

    BounceState exact(double t) const;
    double energy(double y, double vx, double vy) const { return 0.5 * (vx * vx + vy * vy) + gravity * y; }  // par unité de masse

    // Instants des `maxCount` premiers impacts (moins s'il y en a moins avant le repos).
    std::vector<double> impactTimes(int maxCount) const;
    // Instant où la balle cesse de rebondir (l'infini si e = 1 : elle ne s'arrête jamais).
    double restTime() const;
};

// Pilote de BounceProblem : modèle naïf (pas fixe, rebond appliqué après le pas : y est ramené à 0 et vy inversée) ou événement.
class BounceRun {
public:
    BounceRun(const BounceProblem& problem, ContactModel model);

    void advance(Solver& solver, double dt);   // avance de dt (un seul pas pour le modèle naïf)

    double time() const { return time_; }
    const State& state() const { return y_; }
    int bounces() const { return bounces_; }
    bool resting() const { return resting_; }          // seulement pour l'événement : le naïf ne s'arrête jamais
    double restTime() const { return restTime_; }      // instant où le repos a été décidé (NaN avant)

private:
    BounceProblem problem_;
    ContactModel model_;
    State y_;
    double time_ = 0.0;
    int bounces_ = 0;
    bool resting_ = false;
    bool launch_ = false;      // la balle vient de quitter le sol : décollage par un pas minuscule avant de guetter le prochain impact
    double restTime_ = std::numeric_limits<double>::quiet_NaN();
};

// Intègre jusqu'à tEnd en `steps` pas égaux ; renvoie la distance entre la position calculée et la position exacte.
double bounceError(const BounceProblem& problem, ContactModel model, Solver& solver, int steps, double tEnd);

}  // namespace pl
