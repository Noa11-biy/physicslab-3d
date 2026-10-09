// M4 : gravitation, problème de Kepler à deux corps.
//
// Deux corps m1, m2 qui s'attirent par F = G m1 m2 / r^2 se réduisent à UN corps fictif de masse réduite
// mu_r = m1 m2 / (m1 + m2) placé à la position relative r = r2 - r1 :
//     r'' = -GM r / |r|^3          avec GM = G (m1 + m2)
// (le centre de masse, lui, est en mouvement rectiligne uniforme). On travaille ici directement avec ce problème
// relatif, par unité de masse réduite : `mu` désigne GM [m^3/s^2], ce qui permet les unités normalisées (mu = 1).
//
// Grandeurs conservées (par unité de masse réduite), avec L = r x v :
//     énergie             E = v^2/2 - mu/r
//     moment cinétique    L = r x v                (l'orbite reste dans un plan, et la loi des aires en découle)
//     Laplace-Runge-Lenz  A = v x L - mu r/|r|     (pointe vers le périastre, |A| = mu e : l'orbite est FERMÉE)
// Pour une ellipse de demi-grand axe a et d'excentricité e : E = -mu/(2a), |L| = sqrt(mu a (1 - e^2)), |A| = mu e.
//
// Solution exacte (e < 1) par l'anomalie excentrique E_a :
//     équation de Kepler   M = E_a - e sin(E_a),   M = n t  (anomalie moyenne), n = sqrt(mu/a^3)
//     x = a (cos E_a - e),  y = a sqrt(1 - e^2) sin E_a,  r = a (1 - e cos E_a)
//     3e loi : T = 2 pi sqrt(a^3/mu)
// L'équation de Kepler n'a pas de solution en formules : on la résout par Newton (kepler::solveEccentricAnomaly).
//
// Convention des solveurs : y = [x, y, z | vx, vy, vz] (positions | vitesses, n = 3). L'orbite est dans le plan xy,
// périastre sur +x, sens trigonométrique. Le hamiltonien v^2/2 - mu/r est séparable : Euler symplectique et
// Verlet sont symplectiques ici.
#pragma once

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

namespace kepler {

// Résout M = E_a - e sin(E_a) pour E_a (0 <= e < 1, M quelconque). `iterations` reçoit le nombre de pas de Newton.
double solveEccentricAnomaly(double meanAnomaly, double e, int* iterations = nullptr);

}  // namespace kepler

struct KeplerProblem {
    double mu = 1.0;   // GM [m^3/s^2] ; 1 = unités normalisées
    double a = 1.0;    // demi-grand axe [m]
    double e = 0.5;    // excentricité, 0 <= e < 1

    double period() const;       // T = 2 pi sqrt(a^3/mu)
    double meanMotion() const;   // n = 2 pi / T

    // Départ au périastre : r = a (1 - e) sur +x, vitesse tangentielle maximale sur +y.
    State initialState() const;
    OdeFunction rhs() const;

    // Solution exacte à l'instant t (anomalie moyenne nulle au périastre).
    State exact(double t) const;

    // Grandeurs conservées d'un état quelconque.
    double energy(const State& y) const;                 // v^2/2 - mu/r [J/kg]
    Vec3 angularMomentum(const State& y) const;          // r x v
    Vec3 rungeLenz(const State& y) const;                // v x L - mu r/|r|
    double exactEnergy() const { return -mu / (2.0 * a); }
    double exactAngularMomentum() const;                 // sqrt(mu a (1 - e^2))

    // Angle du périastre dans le plan (direction de A), en radians : sa dérive mesure la précession numérique.
    double periapsisAngle(const State& y) const;

    // Précession par orbite [rad] que Verlet (et Euler symplectique) introduisent à pas fixe dt, prédite par le
    // hamiltonien modifié H~ = H + dt^2 (1/12 p.H''.p - 1/24 |grad V|^2) moyenné sur l'orbite (orbite à e fixé, dt petit
    // devant la période) : -(pi/8) (mu dt^2 / a^3) (4 + e^2) / (1 - e^2)^3. Négative = rétrograde.
    double verletPrecessionPerOrbit(double dt) const;

    // Distance dans l'espace des phases : positions divisées par a, vitesses par sqrt(mu/a).
    double distance(const State& s1, const State& s2) const;
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie la distance finale à la solution exacte.
double keplerError(const KeplerProblem& problem, Solver& solver, int steps, double tEnd);

// Mesure de la précession : repère les passages au périastre d'une trajectoire échantillonnée et suit la rotation
// de la direction du périastre. Un périastre est le point où r.v passe de négatif à positif (r minimal) ; sa position
// est interpolée entre deux pas (Hermite cubique) puis son angle polaire est « déroulé » (sans saut de 2 pi).
// Pour la solution exacte la précession est nulle (orbite fermée) ; un schéma numérique la fait apparaître.
class PeriapsisTracker {
public:
    // À appeler après chaque pas avec le temps et l'état (plan xy). Renvoie true si un périastre vient d'être franchi.
    bool update(double t, const State& y);

    int count() const { return count_; }                 // nombre de périastres détectés
    double precession() const { return unwrapped_; }    // rotation cumulée depuis le 1er périastre détecté [rad]
    // Rotation moyenne par orbite [rad] ; NaN tant que deux périastres n'ont pas été détectés.
    double precessionPerOrbit() const;
    void reset() { *this = PeriapsisTracker{}; }

private:
    bool havePrev_ = false;
    double tPrev_ = 0.0;
    State prev_;
    int count_ = 0;
    double lastAngle_ = 0.0;   // angle polaire du dernier périastre, dans (-pi, pi]
    double unwrapped_ = 0.0;
};

}  // namespace pl
