// M6 : corps rigide. Un solide indéformable tourne autour de son centre de masse (ou d'un point fixe). Son orientation est un
// quaternion unitaire q (rotation repère du corps -> repère fixe, v_fixe = q v_corps q*), sa vitesse angulaire w est exprimée dans le
// repère du CORPS, où le tenseur d'inertie est constant et diagonal I = diag(I1, I2, I3) (axes principaux).
//
// Équations du mouvement. Le moment cinétique dans le repère du corps est L_c = I w, dans le repère fixe L = q L_c q*. La loi
// dL/dt = couple, écrite dans le repère tournant, donne les équations d'EULER :
//     I w' + w x (I w) = tau_corps       soit      w1' = [(I2 - I3) w2 w3 + tau1] / I1   et permutations circulaires.
// L'orientation évolue par la cinématique des quaternions :
//     q' = 1/2 q (0, w)     (produit de Hamilton, w dans le repère du corps).
// Sans couple, trois grandeurs sont conservées : le VECTEUR L dans le repère fixe (symétrie par rotation), l'énergie
// E = 1/2 w . I w (symétrie par translation dans le temps) et la norme |q| = 1 (contrainte, que l'intégration numérique ne respecte
// pas toute seule). En repère du corps L_c tourne, mais sa norme est constante.
//
// Cas de référence avec solution connue (voir FreeBodyProblem) : solide symétrique libre (I1 = I2), solide asymétrique libre
// (fonctions elliptiques de Jacobi pour w(t)) ; la toupie pesante de Lagrange est dans HeavyTopProblem.
//
// Convention des solveurs : y = [q.w q.x q.y q.z | w1 w2 w3] (7 nombres). Ce n'est PAS la convention « positions | vitesses » :
// seuls les solveurs génériques (Euler explicite, RK4, RK45) s'appliquent ; Euler symplectique et Verlet sont exclus. Pour respecter
// la contrainte |q| = 1 on dispose des intégrateurs d'orientation (RotationIntegrator, plus bas).
#pragma once

#include <functional>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Mat3.hpp"
#include "physicslab/core/Quaternion.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

// Moments d'inertie principaux de solides homogènes de masse m (axes principaux x, y, z passant par le centre de masse).
namespace inertia {
Vec3 sphere(double m, double r);                       // pleine : 2/5 m r^2 (les trois)
Vec3 box(double m, double a, double b, double c);      // côtés a, b, c selon x, y, z : I_x = m (b^2 + c^2) / 12, ...
Vec3 cylinder(double m, double r, double h);           // axe z : I_z = m r^2 / 2, I_x = I_y = m (3 r^2 + h^2) / 12
}  // namespace inertia

struct RotationState {
    Quaternion q;
    Vec3 omega;   // vitesse angulaire, repère du corps
};

State packRotation(const RotationState& s);            // [q.w q.x q.y q.z | w1 w2 w3]
RotationState unpackRotation(const State& y);

// Fonctions elliptiques de Jacobi sn, cn, dn de paramètre m = k^2 (0 <= m < 1), par leurs séries de Fourier en fonction du nome
// q = exp(-pi K'/K) (DLMF 22.11.1 à 22.11.3) : sn(u) = (2 pi / (K k)) sum_{n>=0} q^(n+1/2) sin((2n+1) zeta) / (1 - q^(2n+1)), avec
// zeta = pi u / (2 K), et de même cn (cosinus, 1 + q^(2n+1)) et dn. Elles sont définies par sn' = cn dn, cn' = -sn dn, dn' = -m sn cn
// et (sn, cn, dn)(0) = (0, 1, 1) ; période 4 K en u.
void jacobiSnCnDn(double u, double m, double& sn, double& cn, double& dn);

// Couple (repère du corps) en fonction de l'orientation, par exemple la pesanteur d'une toupie ; vide = pas de couple.
using TorqueFunction = std::function<Vec3(const Quaternion&)>;

// Accélération angulaire I^-1 [(I w) x w + couple] et second membre complet de y' = f(y).
Vec3 rotationOmegaDot(const Vec3& inertia, const Vec3& omega, const Vec3& torque = {});
OdeFunction rotationRhs(const Vec3& inertia, const TorqueFunction& torque = {});

// Solide libre (sans couple) : trois moments d'inertie principaux, orientation et vitesse angulaire initiales.
struct FreeBodyProblem {
    Vec3 inertia{1.0, 2.0, 3.0};
    Quaternion q0;
    Vec3 omega0{0.5, 1.0, 0.3};

    Vec3 omegaDot(const Vec3& omega) const { return rotationOmegaDot(inertia, omega); }
    State initialState() const { return packRotation({q0, omega0}); }
    OdeFunction rhs() const { return rotationRhs(inertia); }

    double energy(const RotationState& s) const;
    Vec3 angularMomentumBody(const RotationState& s) const;    // I w
    Vec3 angularMomentumSpace(const RotationState& s) const;   // q (I w) q* : constant au cours du temps

    // Solution EXACTE du solide symétrique (I1 = I2, sinon assertion). Avec Omega = (I3 - I1) w3 / I1 : w3 constante et (w1, w2) tourne
    // de l'angle Omega t dans le repère du corps. Comme L = I1 w + (I3 - I1) w3 e3, la vitesse angulaire (repère fixe) vaut
    // (L / I1) L_chapeau - Omega e3 : une PRÉCESSION autour du vecteur fixe L à la vitesse L / I1 composée à une ROTATION PROPRE autour
    // de l'axe e3 du corps à la vitesse -Omega, d'où  q(t) = q_L(L t / I1) q0 q_3(-Omega t)  (rotations autour de L fixe, à gauche, et de
    // l'axe 3 du corps, à droite). L'axe e3 décrit un cône d'axe L, de demi-angle constant (cos = I3 w3 / L).
    RotationState exactSymmetric(double t) const;

    // Solide ASYMÉTRIQUE (I1 < I2 < I3), départ w0 = (a, 0, c) avec a, c > 0 (sinon assertion), hors séparatrice 2 E I2 = L^2.
    // Solution de Landau-Lifchitz (Mécanique, § 37) : avec w = (a cn(lambda t), b sn(lambda t), c dn(lambda t)) si 2 E I2 < L^2 (rotation
    // autour de l'axe 3, w3 presque constante), ou w = (a dn, b sn, c cn) si 2 E I2 > L^2 (autour de l'axe 1) ; b = (I3 - I1) a c / (I2 lambda)
    // et, selon le cas, lambda = c sqrt((I3 - I2)(I3 - I1) / (I1 I2)) ou a sqrt((I3 - I1)(I2 - I1) / (I2 I3)). Le paramètre est
    //     m = (I2 - I1)(2 E I3 - L^2) / ((I3 - I2)(L^2 - 2 E I1))     (ou son inverse dans le second cas),   m < 1  <=>  hors séparatrice.
    // L'orientation q(t) n'a pas de forme fermée élémentaire (fonctions theta) : on la calcule avec `reference`.
    Vec3 exactAsymmetricOmega(double t) const;
    double asymmetricParameter() const;    // m
    double asymmetricPeriod() const;       // période de w(t) : 4 K(k) / lambda

    // Instabilité de l'axe intermédiaire (I1 < I2 < I3) : autour de l'axe 2 à la vitesse w2, une perturbation (w1, w3) croît comme
    // exp(lambda t) avec lambda = |w2| sqrt((I3 - I2)(I2 - I1) / (I1 I3)) (de w1'' = lambda^2 w1, en linéarisant les équations d'Euler).
    double intermediateAxisGrowthRate() const;

    // État à l'instant t par RK45 très serré (relTol 1e-13), sans renormalisation : référence pour l'asymétrique.
    RotationState reference(double t) const;
};

// ---------------------------------------------------------------------------------------------------------------------------------
// Intégrateurs d'orientation. La contrainte |q| = 1 fait de l'orientation un point d'une SPHÈRE (S^3), pas d'un espace vectoriel : un
// schéma de Runge-Kutta ordinaire en sort. Quatre approches à comparer (couple optionnel, fonction de l'orientation en repère du corps) :
//   * EulerRotation   : Euler explicite sur (q, w), sans renormalisation. Ordre 1. |q| croît de sqrt(1 + (h w)^2 / 4) à chaque pas.
//   * RK4Rotation     : RK4 sur le vecteur de 7 nombres, puis renormalisation de q (projection sur la sphère) ; sans renormalisation si
//                       renormalize = false (la dérive de |q| est alors d'ordre élevé mais cumulative).
//   * LieHeunRotation : schéma de GROUPE DE LIE d'ordre 2 : l'orientation est toujours mise à jour par un produit avec l'exponentielle
//                       d'une rotation, q <- q exp(h <w>), donc |q| = 1 à l'arrondi près ; <w> = moyenne (Heun) de w sur le pas.
//   * SplittingRotation : découpage symplectique de Dullweber-Leimkuhler-McLachlan (1997). Le hamiltonien libre H = sum_k pi_k^2 / (2 I_k)
//                       (pi = I w, moment cinétique dans le repère du corps) est la somme de trois rotateurs H_k = pi_k^2 / (2 I_k), chacun
//                       exactement soluble : autour de l'axe k, w_k = pi_k / I_k est constant, pi tourne de l'angle -theta autour de e_k et q
//                       est multiplié à droite par la rotation de theta = w_k tau autour de e_k. Le pas de Strang
//                       R1(h/2) R2(h/2) R3(h) R2(h/2) R1(h/2) est d'ordre 2, symplectique, réversible. Comme Q pi (le moment cinétique dans le
//                       repère fixe) est inchangé par chaque sous-pas (Q tourne de +theta, pi de -theta), L_fixe est conservé à l'ARRONDI PRÈS
//                       quel que soit le pas, et |q| = 1. Un couple se traite par des demi-coups (kick) de part et d'autre : pi += (h/2) tau(q).
class RotationIntegrator {
public:
    virtual ~RotationIntegrator() = default;
    virtual const char* name() const = 0;
    virtual void step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) = 0;
};

class EulerRotation final : public RotationIntegrator {
public:
    const char* name() const override { return "Euler explicite"; }
    void step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) override;
};

class RK4Rotation final : public RotationIntegrator {
public:
    explicit RK4Rotation(bool renormalize = true) : renormalize_(renormalize) {}
    const char* name() const override { return renormalize_ ? "RK4 + renormalisation" : "RK4 sans renormalisation"; }
    void step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) override;

private:
    bool renormalize_;
    RK4 rk4_;
};

class LieHeunRotation final : public RotationIntegrator {
public:
    const char* name() const override { return "Groupe de Lie (Heun)"; }
    void step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) override;
};

class SplittingRotation final : public RotationIntegrator {
public:
    const char* name() const override { return "Découpage symplectique"; }
    void step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) override;
};

// Quaternion unitaire de la rotation d'angle |v| autour de v (exponentielle de l'algèbre so(3)) ; l'identité pour v = 0.
Quaternion quaternionExp(const Vec3& rotationVector);

// Avance de tEnd en `steps` pas égaux.
RotationState integrateRotation(RotationIntegrator& integrator, const Vec3& inertia, const TorqueFunction& torque, RotationState s, double tEnd,
                                int steps);

// Distance entre deux états : sqrt(|R_a - R_b|_F^2 + |w_a - w_b|^2), orientations normalisées (insensible au signe de q).
double rotationDistance(const RotationState& a, const RotationState& b);

// ---------------------------------------------------------------------------------------------------------------------------------
// Toupie pesante de Lagrange : corps de révolution (I1 = I2, I3 axial, moments AUTOUR DU PIVOT) fixé en un point, dans la pesanteur. Le
// centre de masse est sur l'axe e3, à la distance `lever` du pivot. Avec theta l'angle de l'axe e3 avec la verticale ascendante (u = cos theta
// = e3_fixe . z), le couple de la pesanteur, écrit dans le repère du corps, est  tau = m l e3 x (R^T g)  (g = (0, 0, -g)) : il est
// perpendiculaire à e3. Trois intégrales premières (Lagrange) :
//     E = 1/2 w . I w + m g l cos theta,     L_z (vertical, pas de couple autour de la verticale),     L_3 = I3 w3 (couple perpendiculaire à e3).
// En angles d'Euler (phi précession, theta nutation, psi rotation propre) ils réduisent le mouvement à une quadrature en u = cos theta :
//     u'^2 = f(u) = (2/I1)(E' - m g l u)(1 - u^2) - (L_z - L_3 u)^2 / I1^2,        E' = E - L_3^2 / (2 I3),
// polynôme cubique de racines u1 < u2 < 1 < u3 : la nutation oscille entre u1 et u2, avec la période T = 4 K(k) / sqrt(beta (u3 - u1)),
// beta = 2 m g l / I1, k^2 = (u2 - u1)/(u3 - u1). Cas particuliers exacts :
//   * précession régulière (theta constant) : I1 phi'^2 cos theta - I3 w3 phi' + m g l = 0 (branches lente et rapide), w3 = psi' + phi' cos theta ;
//   * toupie « endormie » (theta = 0) : stable si I3^2 w3^2 > 4 I1 m g l ; sinon une inclinaison croît au taux
//     gamma = sqrt(4 I1 m g l - I3^2 w3^2) / (2 I1).
struct SteadyPrecession {
    bool valid = false;      // faux si le discriminant est négatif (spin trop faible : pas de précession régulière) ou si la branche n'existe pas
    double theta = 0.0, spin = 0.0;    // inclinaison et composante w3 de la vitesse angulaire (repère du corps)
    double phiDot = 0.0, psiDot = 0.0; // vitesses de précession et de rotation propre
};

struct NutationRange {
    bool valid = false;      // faux sans nutation (trois racines réelles distinctes introuvables)
    double u1 = 0.0, u2 = 0.0, u3 = 0.0;
};

struct HeavyTopProblem {
    double mass = 1.0;
    Vec3 inertia{1.2, 1.2, 0.4};       // autour du pivot ; I1 = I2 requis par les formules fermées
    double lever = 0.5;                // pivot -> centre de masse, le long de e3
    double gravity = constants::g0;
    Quaternion q0;
    Vec3 omega0{0.0, 0.0, 10.0};

    TorqueFunction torque() const;
    State initialState() const { return packRotation({q0, omega0}); }
    OdeFunction rhs() const { return rotationRhs(inertia, torque()); }

    double cosTheta(const Quaternion& q) const { return q.rotate({0.0, 0.0, 1.0}).z; }
    double energy(const RotationState& s) const;
    Vec3 angularMomentumSpace(const RotationState& s) const;
    double verticalMomentum(const RotationState& s) const { return angularMomentumSpace(s).z; }
    double axialMomentum(const RotationState& s) const { return inertia.z * s.omega.z; }

    // État à l'instant t par RK45 très serré (relTol 1e-13).
    RotationState reference(double t) const;

    // Précession régulière d'inclinaison theta et de spin w3 (> 0) ; branche lente (phi' petite) ou rapide.
    SteadyPrecession steadyPrecession(double theta, double spin, bool slow) const;
    void startSteady(const SteadyPrecession& sp);                               // q0 = q_x(theta), w0 = (0, phi' sin theta, w3)
    RotationState exactSteady(const SteadyPrecession& sp, double t) const;     // q(t) = q_z(phi' t) q_x(theta) q_z(psi' t)

    // Racines u1 < u2 < u3 du polynôme cubique de la nutation, à partir de l'état initial, et période de nutation.
    NutationRange nutation() const;
    double nutationPeriod() const;

    // Spin critique de la toupie endormie : I3 w3 = 2 sqrt(I1 m g l).
    double sleepingCriticalSpin() const;
};

}  // namespace pl
