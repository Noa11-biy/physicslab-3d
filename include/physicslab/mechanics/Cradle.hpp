// M5c : berceau de Newton. Une chaîne de billes identiques qui se touchent sur un rail horizontal (1D), les `launched` premières
// lancées vers les autres. Deux modèles de choc à comparer : le contact de Hertz (ici) et les impulsions séquentielles (plus bas).
//
// Contact de Hertz. Deux sphères élastiques de rayon R se touchent sur une petite surface qui grandit avec l'écrasement : avec
// delta = 2R - (distance des centres) > 0 l'interpénétration (écrasement total des deux billes), la théorie de Hertz (1881) donne
//     F = k delta^(3/2)         (répulsive, nulle si delta <= 0 : une bille ne tire jamais sa voisine),
//     U = (2/5) k delta^(5/2)   (l'énergie élastique stockée, dérivée de F).
// La raideur d'un ressort, dF/d(delta) = (3/2) k sqrt(delta), CROÎT avec l'écrasement : un ressort non linéaire. Chaque bille ne
// subit que ses deux voisines, et la force d'un contact est opposée de chaque côté (3e loi) : l'impulsion P = sum m v est conservée, et,
// sans amortissement, l'énergie E = sum (1/2) m v^2 + sum U(delta_i) aussi.
//
// Deux billes identiques (masse m, masse réduite mu = m/2), vitesse relative d'approche v. En coordonnée relative,
//     mu delta'' = -k delta^(3/2),   delta(0) = 0,   delta'(0) = v,
// d'énergie 1/2 mu delta'^2 + (2/5) k delta^(5/2) = 1/2 mu v^2. L'écrasement est maximal quand delta' = 0 :
//     delta_max = (5 mu v^2 / (4 k))^(2/5).
// La durée du contact (aller et retour) vaut, avec le changement de variable s = (delta/delta_max)^(5/2) qui ramène à une fonction bêta,
//     T = 2 (delta_max / v) int_0^1 du / sqrt(1 - u^(5/2)) = (4/5) B(2/5, 1/2) delta_max / v = 2,9433 delta_max / v.
// La sortie est élastique (énergie conservée par le contact) : les deux billes échangent leurs vitesses, comme dans un choc parfait.
// Lois d'échelle : delta_max ~ v^(4/5) k^(-2/5), T ~ v^(-1/5) k^(-2/5) : un contact plus raide est plus court, mais le RÉSULTAT d'une
// chaîne ne dépend ni de k ni de v (une seule échelle de longueur (m v^2/k)^(2/5) apparaît dans les équations).
//
// Amortissement optionnel (Hunt-Crossley, 1975) : F = k delta^(3/2) (1 + (3/2) alpha delta'), delta' = vitesse d'écrasement, ramenée à 0
// si le facteur devient négatif (la force ne devient jamais attractive). Dissipation pendant le contact seulement. À la main (trajectoire
// conservative à l'ordre dominant), l'énergie relative perdue en un choc vaut  (4/5) alpha k v delta_max^(5/2) = alpha mu v^3, soit
// 2 alpha v de l'énergie d'approche : restitution e = 1 - alpha v au premier ordre. L'impulsion reste conservée ; l'énergie E de
// energy() ne l'est plus (elle décroît) dès que alpha > 0.
//
// Convention des solveurs : y = [x_0 ... x_{N-1} | v_0 ... v_{N-1}] (positions | vitesses, n = N). Les billes sont numérotées de
// gauche à droite ; au départ elles se touchent (centres espacés de 2R, aucune compression) et les `launched` premières ont la vitesse
// `speed` : leur contact avec la bille suivante démarre à t = 0, ou à t = gap / speed si le groupe lancé est reculé de `gap`.
#pragma once

#include <array>
#include <vector>

#include "physicslab/core/Events.hpp"
#include "physicslab/core/Solver.hpp"

namespace pl {

namespace hertz {

// (4/5) Gamma(2/5) Gamma(1/2) / Gamma(9/10) = 2,9433 : durée du contact en unités de delta_max / v.
double contactConstant();

// Force de Hertz k delta^(3/2) et énergie élastique (2/5) k delta^(5/2) ; nulles si delta <= 0.
double force(double stiffness, double delta);
double potentialEnergy(double stiffness, double delta);

// Choc de deux billes identiques, masse réduite mu, vitesse d'approche v (références exactes, voir l'en-tête).
double maxCompression(double reducedMass, double approachSpeed, double stiffness);
double contactDuration(double reducedMass, double approachSpeed, double stiffness);

}  // namespace hertz

namespace cradle {

// Une bille de vitesse v frappe deux billes identiques au repos et qui se touchent. L'impulsion (v1 + v2 + v3 = v) et l'énergie
// (v1^2 + v2^2 + v3^2 = v^2) NE DÉTERMINENT PAS le résultat : on en tire v1 v2 + v1 v3 + v2 v3 = 0, c'est-à-dire un cercle de solutions
// dans l'espace des vitesses. Avec les vitesses rangées (v1 <= v2 <= v3, plus aucun choc ensuite) et v1 = -x v :
//     v2, v3 racines de  t^2 - (1 + x) v t + x (1 + x) v^2 = 0,   de discriminant (1 + x)(1 - 3x) v^2 >= 0  =>  0 <= x <= 1/3.
// x = 0 : « une entre, une sort » (0, 0, v). x = 1/3 : (-v/3, 2v/3, 2v/3), les deux dernières billes partent ensemble. Seule la
// DYNAMIQUE du contact choisit le point de la courbe. `x` hors de [0, 1/3] est ramené à la borne la plus proche.
std::array<double, 3> threeBallFamily(double x, double speed = 1.0);

// Modèle 1 : impulsions séquentielles. Les chocs sont instantanés et BINAIRES : on applique à la fois un seul choc de deux billes
// voisines qui se rapprochent (restitution de Newton, collision::collide1D), puis on recommence jusqu'à ce qu'aucune paire ne se
// rapproche. Avec e = 1 et des masses égales chaque choc échange les vitesses des deux billes : le résultat est le tri des vitesses,
// donc « n entrent, n sortent », QUEL QUE SOIT l'ordre. Avec e < 1 l'ordre de résolution peut changer le résultat : le modèle
// n'a aucune dynamique qui le départage. Ce qu'il ne sait jamais faire, c'est un choc à trois billes simultané (la courbe de
// threeBallFamily lui est inaccessible, hors de son point x = 0 quand e = 1).
enum class ResolveOrder {
    LeftToRight,   // passes de la gauche vers la droite
    RightToLeft    // passes de la droite vers la gauche
};

struct ImpulseResult {
    std::vector<double> velocities;   // vitesses après tous les chocs
    int collisions = 0;               // chocs binaires appliqués
    int sweeps = 0;                   // passes complètes sur la chaîne (la dernière n'applique rien)
    bool converged = true;            // faux si la limite de passes est atteinte (e proche de 0 : convergence géométrique lente)
    double energyLoss = 0.0;          // énergie cinétique perdue, pour des billes de masse 1
};

// Résout les chocs entre billes identiques voisines (même masse) à partir des vitesses données. Une paire se rapproche si
// v_gauche - v_droite dépasse 1e-12 fois la plus grande vitesse (seuil d'arrondi pour e < 1).
ImpulseResult sequentialImpulses(std::vector<double> velocities, double restitution, ResolveOrder order = ResolveOrder::LeftToRight);

}  // namespace cradle

struct CradleOutcome {
    std::vector<double> velocities;   // vitesses de chaque bille à la fin de la collision
    double time = 0.0;                // instant (au pas près) où plus aucune bille ne peut en toucher une autre
    double maxCompression = 0.0;      // plus forte interpénétration observée (échantillonnée aux pas)
    bool finished = false;            // faux si tMax ou le budget de pas ont été atteints avant
    int steps = 0;
};

struct CradleProblem {
    int balls = 5;                    // N, nombre de billes (>= 2)
    int launched = 1;                 // nombre de billes lancées (les `launched` premières), 1 <= launched < balls
    double mass = 1.0;                // masse de chaque bille
    double radius = 0.5;              // rayon R (le rail est idéal : pas de pesanteur ni de frottement)
    double stiffness = 1.0e4;         // k de Hertz [force / longueur^(3/2)] : unités normalisées, contact court devant le mouvement
    double speed = 1.0;               // vitesse des billes lancées
    double damping = 0.0;             // alpha de Hunt-Crossley [temps/longueur] : 0 = contact conservatif (voir ci-dessus)
    double gap = 0.0;                 // distance de vol libre avant le contact : le groupe lancé est reculé de `gap` (>= 0)

    State initialState() const;
    OdeFunction rhs() const;

    // Écrasement du contact entre la bille `pair` et la suivante (peut être négatif : les billes sont séparées).
    double compression(const State& y, int pair) const;

    double kineticEnergy(const State& y) const;
    double potentialEnergy(const State& y) const;        // somme des énergies de Hertz
    double energy(const State& y) const { return kineticEnergy(y) + potentialEnergy(y); }
    double momentum(const State& y) const;

    // Vrai quand plus aucune bille ne peut en toucher une autre : pour chaque paire voisine, plus de compression et vitesses
    // rangées (la gauche n'est pas plus rapide que la droite). L'état est alors définitif : le mouvement est libre.
    bool collisionOver(const State& y) const;

    // Durée d'un contact isolé (deux billes, vitesse `speed`) divisée par 100 : un pas raisonnable pour un schéma d'ordre 4.
    double suggestedStep() const;

    // Intègre depuis l'état initial jusqu'à la fin de la collision (ou tMax), par pas de dt (le solveur peut en faire plus s'il
    // est adaptatif).
    CradleOutcome run(Solver& solver, double dt, double tMax) const;

    // État à l'instant t depuis l'état initial, par RK45 serré (relTol 1e-13) : référence pour mesurer l'ordre des schémas.
    State reference(double t) const;
};

// Intègre jusqu'à tEnd en `steps` pas égaux et renvoie la distance finale à la référence RK45. Distance dans l'espace des phases :
// positions en unités de  speed * (durée d'un contact isolé),  vitesses en unités de `speed`.
double cradleError(const CradleProblem& problem, Solver& solver, int steps, double tEnd);

}  // namespace pl
