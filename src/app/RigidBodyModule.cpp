#include "RigidBodyModule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <implot.h>

#include "physicslab/core/Constants.hpp"

namespace pl {
namespace {

constexpr double kPi = constants::pi;
constexpr double kDeg = 180.0 / kPi;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kFloor = 1e-17;           // plancher des courbes en échelle logarithmique
constexpr std::size_t kTrailLength = 3000;

// Couleur de chaque intégrateur (mêmes couleurs que SolverSet quand le schéma existe là-bas).
const float kMethodColors[5][3] = {{1.00f, 0.55f, 0.15f},    // Euler : orange
                                   {0.35f, 0.85f, 0.45f},    // RK4 + renormalisation : vert
                                   {0.30f, 0.80f, 0.80f},    // RK4 sans renormalisation : turquoise
                                   {0.85f, 0.45f, 0.95f},    // groupe de Lie : violet
                                   {0.80f, 0.85f, 0.25f}};   // découpage symplectique : citron
const float kLColor[3] = {1.0f, 0.82f, 0.30f};               // moment cinétique
const float kWColor[3] = {0.35f, 1.0f, 0.45f};               // vitesse angulaire
const float kAxisColors[3][3] = {{0.90f, 0.25f, 0.25f}, {0.30f, 0.80f, 0.30f}, {0.30f, 0.50f, 0.95f}};

// Repère physique (z vers le haut, pesanteur -z) -> repère de la scène (y vers le haut) : (x, y, z) -> (x, z, -y), une rotation propre.
Vec3 toScene(const Vec3& p) { return {p.x, p.z, -p.y}; }

// q v q* pour un quaternion QUELCONQUE (pas forcément unitaire) : un quaternion qui a dérivé dilate le solide de |q|^2, ce qu'on veut voir.
Vec3 sandwich(const Quaternion& q, const Vec3& v) {
    const Vec3 u{q.x, q.y, q.z};
    return (q.w * q.w - u.norm2()) * v + 2.0 * dot(u, v) * u + 2.0 * q.w * cross(u, v);
}

Vertex vertexAt(const Vec3& scenePos, const float* c, float dim = 1.0f) {
    return {static_cast<float>(scenePos.x), static_cast<float>(scenePos.y), static_cast<float>(scenePos.z), c[0] * dim, c[1] * dim, c[2] * dim};
}

// Géométrie du corps en repère du corps : suite de couples de points (segments).
void addLine(std::vector<Vec3>& out, const Vec3& a, const Vec3& b) {
    out.push_back(a);
    out.push_back(b);
}

void addCircle(std::vector<Vec3>& out, double z, double r, int n = 28) {
    for (int i = 0; i < n; ++i) {
        const double a0 = 2.0 * kPi * i / n, a1 = 2.0 * kPi * (i + 1) / n;
        addLine(out, {r * std::cos(a0), r * std::sin(a0), z}, {r * std::cos(a1), r * std::sin(a1), z});
    }
}

std::vector<Vec3> cylinderLines(double r, double h) {
    std::vector<Vec3> v;
    addCircle(v, 0.5 * h, r);
    addCircle(v, -0.5 * h, r);
    for (int k = 0; k < 4; ++k) {
        const double a = 0.5 * kPi * k;
        addLine(v, {r * std::cos(a), r * std::sin(a), -0.5 * h}, {r * std::cos(a), r * std::sin(a), 0.5 * h});
    }
    addLine(v, {0.0, 0.0, 0.5 * h}, {r, 0.0, 0.5 * h});             // repère sur le couvercle : on voit tourner le corps sur lui-même
    addLine(v, {0.0, 0.0, -0.5 * h - 0.5}, {0.0, 0.0, 0.5 * h + 0.5});   // axe de symétrie
    return v;
}

std::vector<Vec3> boxLines(double hx, double hy, double hz) {
    std::vector<Vec3> v;
    for (int i = 0; i < 4; ++i) {   // arêtes parallèles à chaque axe
        const double sy = (i & 1) ? hy : -hy, sz = (i & 2) ? hz : -hz;
        addLine(v, {-hx, sy, sz}, {hx, sy, sz});
        const double sx = (i & 1) ? hx : -hx;
        addLine(v, {sx, -hy, sz}, {sx, hy, sz});
        const double sy2 = (i & 1) ? hy : -hy;
        addLine(v, {sx, sy2, -hz}, {sx, sy2, hz});
    }
    return v;
}

std::vector<Vec3> topLines(double rod, double discRadius, double discPos) {
    std::vector<Vec3> v;
    addLine(v, {0.0, 0.0, 0.0}, {0.0, 0.0, rod});                        // tige de l'axe
    addCircle(v, discPos, discRadius);                                     // volant
    for (int k = 0; k < 4; ++k) {
        const double a = 0.5 * kPi * k;
        addLine(v, {0.0, 0.0, discPos}, {discRadius * std::cos(a), discRadius * std::sin(a), discPos});
    }
    addLine(v, {0.0, 0.0, discPos}, {discRadius * 1.15, 0.0, discPos});   // repère de rotation propre
    return v;
}

// Quaternion qui amène le vecteur (en repère du corps) lBody sur la verticale +z.
Quaternion alignWithVertical(const Vec3& lBody) {
    const Vec3 u = lBody.normalized(), z{0.0, 0.0, 1.0};
    const Vec3 axis = cross(u, z);
    const double s = axis.norm(), c = dot(u, z);
    if (s < 1e-12) return c > 0.0 ? Quaternion{} : Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, kPi);
    return Quaternion::fromAxisAngle(axis, std::atan2(s, c));
}

// Limites d'un axe pour une courbe presque constante : au moins `minSpan` d'étendue autour de la plage observée, sinon les graduations
// se répètent (« 0,825336 » partout) et la courbe plate occupe tout le cadre.
void paddedLimits(const Series& s, double minSpan, double& lo, double& hi) {
    lo = std::numeric_limits<double>::max();
    hi = -std::numeric_limits<double>::max();
    for (double v : s.y) {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    if (s.size() == 0) { lo = 0.0; hi = 1.0; }
    const double span = hi - lo;
    if (span < minSpan) {
        const double mid = 0.5 * (hi + lo);
        lo = mid - 0.5 * minSpan;
        hi = mid + 0.5 * minSpan;
    } else {
        lo -= 0.05 * span;
        hi += 0.05 * span;
    }
}

const char* const kScenarioNames[3] = {"Solide symétrique libre", "Solide asymétrique libre", "Toupie pesante"};

}  // namespace

RigidBodyModule::RigidBodyModule() {
    ref_.relTol = 1e-12;
    ref_.absTol = 1e-14;
    reset();
}

// ------------------------------- textes --------------------------------

const char* RigidBodyModule::explanation(Level level) const {
    static const char* const kText[3][6] = {
        // ---- solide symétrique libre
        {"Lancez un frisbee ou une pièce en la faisant tourner : si le lancer est bien droit l'objet tourne tranquillement, s'il est un peu de "
         "travers il « vacille ». Dans l'espace rien ne le pousse : il garde son « élan de rotation », une flèche jaune qui reste fixe, et c'est "
         "l'axe de l'objet qui tourne autour d'elle en dessinant un cône (la courbe tracée).\n\nLe solide bleu est la vraie trajectoire, le vert "
         "est le calcul de l'ordinateur : ils se superposent.",
         "Un objet qui tourne garde son moment cinétique : un vecteur (flèche jaune) dont la direction et la taille ne changent pas tant que rien "
         "ne le pousse. La rotation instantanée (flèche verte) n'est pas alignée avec lui dès que l'objet n'est pas une sphère : voilà pourquoi "
         "l'axe vacille.\n\nUn disque très plat vacille une fois par tour de spin (comme un frisbee mal lancé), une sphère ne vacille pas, un "
         "cigare très fin vacille aussi une fois par tour, mais dans l'autre sens. Curseurs : forme, spin, écart de l'axe.",
         "Le moment cinétique L = I × ω joue pour la rotation le rôle de p = m × v : I est le moment d'inertie (la « masse » de rotation, différente "
         "pour chaque axe). Sans force, L (direction et valeur) et l'énergie ½ I ω² sont constants.\n\nPour un solide symétrique, deux moments "
         "d'inertie sont égaux (I1) et le troisième (I3) diffère. L'axe tourne autour de L à la vitesse L / I1 (rad/s). Le tableau compare le "
         "calcul à la formule exacte : Euler (« méthode simple ») fait gonfler le solide.",
         "Équations d'Euler (repère du corps) : I1 ω1' = (I2 − I3) ω2 ω3, et permutations circulaires. Pour I1 = I2 : ω3 est constante et (ω1, ω2) "
         "tourne à Ω = (I3 − I1) ω3 / I1 dans le corps. Dans l'espace, l'axe décrit un cône autour de L à la vitesse L / I1, avec "
         "cos θ = I3 ω3 / L.\n\nL'énergie E = ½ Σ I_i ω_i² et le vecteur L (repère fixe) sont conservés : voyez le tableau. Le calcul naïf (Euler) "
         "les viole.",
         "Orientation par un quaternion unitaire q : q' = ½ q (0, ω). Solution exacte : q(t) = q_L(L t / I1) q0 q_3(−Ω t) : précession autour du vecteur "
         "fixe L à gauche, rotation propre autour de l'axe e3 du corps à droite.\n\nCinq intégrateurs : Euler (|q| croît de √(1 + (h ω)²/4) par pas), "
         "RK4 avec ou sans renormalisation (projection sur la sphère |q| = 1), groupe de Lie (q ← q exp(h <ω>)) et découpage symplectique (rotations "
         "exactes autour des axes principaux : L et |q| exacts, E borné). « Analyse » : dérives et ordres de convergence mesurés.",
         "Intégrer sur S³ : RK4 sort de la contrainte |q| = 1 (la renormalisation la rétablit mais pas L) ; un schéma de groupe de Lie la respecte "
         "(|q| = 1 à l'arrondi) sans conserver L ni E (Heun, mesuré sur un corps I = (1, 2, 3) proche de l'axe instable : dérive d'énergie de 1,4·10⁻² "
         "sur 1000 s à h = 0,05, divergence à h = 0,2). Découpage de Dullweber-Leimkuhler-McLachlan : H = Σ π_k²/(2 I_k) est la somme de trois rotateurs "
         "exactement solubles ; pas de Strang d'ordre 2, symplectique, réversible ; Q π (L fixe) est inchangé à chaque sous-pas : conservation à "
         "l'arrondi pour tout pas, E borné en h² sans dérive (mesuré).\n\nPour I1 = I2 la solution exacte est fermée : l'écart affiché est celui du schéma seul."},
        // ---- solide asymétrique libre
        {"Lancez un livre ou une raquette de tennis en les faisant tourner. Autour de leur grande ou de leur petite direction, l'objet tourne "
         "proprement. Autour de la direction « du milieu », il bascule tout seul d'un coup, puis revient, encore et encore : c'est l'effet de la "
         "raquette de tennis. Choisissez l'axe et regardez : bleu = la vraie trajectoire, vert = le calcul de l'ordinateur.",
         "Un objet a trois axes principaux. Autour de l'axe le plus long et de l'axe le plus court, une petite perturbation reste petite : la rotation "
         "est stable. Autour de l'axe du milieu, la perturbation grandit de plus en plus vite (de façon exponentielle) et l'objet se retourne "
         "régulièrement. Pourtant son énergie et son moment cinétique (flèche jaune) ne changent pas.",
         "Les trois moments d'inertie de la boîte 3 × 2 × 1 (masse 1) sont I1 ≈ 0,42 < I2 ≈ 0,83 < I3 ≈ 1,08 : le plus petit pour l'axe long, le "
         "plus grand pour l'axe court. L = I × ω et E = ½ I ω² sont constants. Le tableau compare les calculs à la référence : près de l'axe du "
         "milieu, une toute petite erreur grandit vite, c'est un cas difficile pour un ordinateur.",
         "Équations d'Euler : I1 ω1' = (I2 − I3) ω2 ω3, etc. Autour de l'axe 2 à la vitesse ω2 : ω1'' = λ² ω1 avec "
         "λ = ω2 √((I3 − I2)(I2 − I1)/(I1 I3)) : croissance en e^(λ t), instable ; autour des axes 1 et 3 les deux facteurs ont des signes opposés : "
         "oscillation, stable.\n\nL'énergie et le vecteur L (repère fixe) sont conservés ; le tableau suit leurs dérives pour chaque méthode.",
         "Solution exacte de ω(t) par les fonctions de Jacobi pour un départ (a, 0, c) : ω = (a cn, b sn, c dn) ou (a dn, b sn, c cn) selon le côté de la "
         "séparatrice 2 E I2 = L², période 4 K(k)/λ. L'orientation se calcule numériquement (référence : RK45 serré).\n\nCourbe « polhode » (Analyse) : "
         "trajectoire de ω dans le corps, intersection de l'ellipsoïde d'énergie et de la sphère |L| = constante : des courbes fermées autour des "
         "axes stables ; l'axe du milieu est un point selle.",
         "Les mêmes cinq intégrateurs sur un corps asymétrique. Mesuré sur 1000 s (corps I = (1, 2, 3), h = 0,05, départ (0,1 ; 2 ; 0,1) proche de l'axe "
         "instable) : découpage L à 6·10⁻¹⁴ et E borné (7·10⁻⁵) ; RK4 E −3·10⁻⁶ et L 2·10⁻⁶ ; Heun E +1,4·10⁻² ; Euler diverge. Ordres mesurés : "
         "Euler 1, RK4 4, Heun 2, découpage 2.\n\nQuand le départ est (a, 0, c), la solution de Jacobi donne ω(t) à 10⁻¹² près (mesuré)."},
        // ---- toupie pesante
        {"Une toupie qui tourne vite ne tombe pas : au lieu de se coucher, son axe tourne lentement autour de la verticale, c'est la précession. Si "
         "elle tourne trop lentement, elle s'incline franchement. Dressée et « endormie », elle reste droite tant qu'elle tourne assez vite. La "
         "courbe tracée est le sommet de la toupie. Bleu : le vrai mouvement, vert : le calcul de l'ordinateur.",
         "La pesanteur tire le centre de la toupie vers le bas et crée un couple : au lieu de la coucher, il fait tourner son axe autour de la "
         "verticale (précession gyroscopique). Plus la toupie tourne vite, plus la précession est lente.\n\nQuatre cas : précession lente régulière, "
         "précession rapide (rare), nutation (la toupie « ondule » en précessant) et toupie endormie (droite, en train de tourner).",
         "Pour une toupie qui tourne vite, la vitesse de précession vaut environ m g l / (I ω) : elle diminue quand ω augmente. Conservation : "
         "l'énergie (mouvement + hauteur), le moment cinétique autour de la verticale et celui autour de l'axe de la toupie ne changent pas. Le "
         "tableau suit l'énergie et le moment vertical de chaque méthode.",
         "Couple de la pesanteur τ = r × F = m l e3 × (Rᵀ g), perpendiculaire à l'axe. Précession régulière à inclinaison θ constante : "
         "m g l = φ' (I3 ω3 − I1 φ' cos θ), deux solutions φ' (lente et rapide) si I3² ω3² ≥ 4 I1 m g l cos θ. Si le spin est trop faible il n'y a "
         "pas de précession régulière : la toupie s'incline et oscille.",
         "Intégrales premières : E, L_z, L_3. Avec u = cos θ : u'² = f(u) = (2/I1)(E' − m g l u)(1 − u²) − (L_z − L_3 u)²/I1², cubique de racines "
         "u1 < u2 < 1 < u3 ; la nutation oscille entre u1 et u2 avec la période T = 4 K(k)/√(β (u3 − u1)), β = 2 m g l / I1, k² = (u2 − u1)/(u3 − u1). "
         "Toupie endormie stable si I3² ω3² > 4 I1 m g l ; sinon l'inclinaison croît en e^(γ t), γ = √(4 I1 m g l − I3² ω3²)/(2 I1).",
         "Intégrateurs avec couple : le découpage ajoute des demi-coups de couple τ(q) autour des rotateurs exacts ; le moment vertical L_z (le couple n'a "
         "aucune composante verticale) est conservé à l'arrondi pour tout pas (mesuré 3·10⁻¹³ sur 100 s à h = 0,005, contre 5·10⁻⁶ pour RK4 et "
         "7,5·10⁻³ pour Heun). Ordres mesurés avec couple : RK4 4, Heun 2, découpage 2.\n\nSous le spin critique la toupie ne « tombe » pas : fixée "
         "au pivot et sans perte, elle effectue une grande nutation entre la verticale et environ 1,25 rad (mesuré à 0,7 fois le spin critique) puis remonte."}};
    return kText[scenario_][static_cast<int>(level)];
}

// ---------------------------- intégrateurs -----------------------------

bool RigidBodyModule::isShown(Level level, int m) const {
    if (!atLeast(level, Level::College)) return m == kRK4;
    if (!atLeast(level, Level::Etudiant)) return m == kEuler || m == kRK4;
    return show_[m];
}

std::string RigidBodyModule::label(Level level, int m) const {
    if (!atLeast(level, Level::College) && m == kRK4) return "Ordinateur";
    if (!atLeast(level, Level::Etudiant)) {
        if (m == kEuler) return "Méthode simple (Euler)";
        if (m == kRK4) return "Méthode précise (RK4)";
    }
    return runs_[m].integrator->name();
}

std::string RigidBodyModule::shortLabel(int m) const {
    static const char* names[kMethodCount] = {"Euler", "RK4 renorm.", "RK4 brut", "Lie-Heun", "Découpage"};
    return names[m];
}

const float* RigidBodyModule::color(int m) const { return kMethodColors[m]; }

// ------------------------------ physique -------------------------------

double RigidBodyModule::energyOf(const RotationState& s) const {
    RotationState n = s;
    n.q = s.q.normalized();   // l'énergie de la toupie dépend de l'inclinaison : on la lit sur l'orientation normalisée
    return scenario_ == kTop ? top_.energy(n) : free_.energy(n);
}

double RigidBodyModule::momentumDrift(const RotationState& s) const {
    RotationState n = s;
    n.q = s.q.normalized();
    if (scenario_ == kTop) return std::abs(top_.verticalMomentum(n) - momentum0_) / std::abs(momentum0_);
    return (free_.angularMomentumSpace(n) - lSpace0_).norm() / momentum0_;
}

Vec3 RigidBodyModule::axisTip(int scenario) const {
    if (scenario == kSymmetric) return {0.0, 0.0, 1.7};
    if (scenario == kTop) return {0.0, 0.0, 3.0};
    return axis_ == 0 ? Vec3{2.0, 0.0, 0.0} : (axis_ == 1 ? Vec3{0.0, 2.0, 0.0} : Vec3{0.0, 0.0, 2.0});
}

RotationState RigidBodyModule::referenceState() const {
    if (scenario_ == kSymmetric) return free_.exactSymmetric(clock_.time);
    return unpackRotation(yRef_);
}

void RigidBodyModule::reset() {
    steadyFallback_ = false;
    steady_ = SteadyPrecession{};
    bodyLines_.clear();

    if (scenario_ == kSymmetric) {
        free_ = FreeBodyProblem{};
        free_.inertia = {1.0, 1.0, shapeRatio_};
        free_.omega0 = {wobble_, 0.0, spin_};
        free_.q0 = alignWithVertical({free_.inertia.x * free_.omega0.x, free_.inertia.y * free_.omega0.y, free_.inertia.z * free_.omega0.z});
        inertia_ = free_.inertia;
        torque_ = {};
        initial_ = {free_.q0, free_.omega0};
        // Cylindre de rayon 1 dont le rapport I3 / I1 = 6 / (3 + h²) vaut shapeRatio_ (h = 0 : disque infiniment plat).
        const double h = std::sqrt(std::max(6.0 / shapeRatio_ - 3.0, 0.0));
        const double scale = 2.4 / std::max(2.0, h);
        bodyLines_ = cylinderLines(scale, std::max(scale * h, 0.06));
        duration_ = 20.0;
    } else if (scenario_ == kAsymmetric) {
        free_ = FreeBodyProblem{};
        free_.inertia = inertia::box(1.0, 3.0, 2.0, 1.0);     // boîte 3 x 2 x 1 : I1 < I2 < I3 (axe long, du milieu, court)
        const double w = 2.0, e = perturb_;
        if (axis_ == 0) free_.omega0 = {w, 0.0, e};
        else if (axis_ == 1) free_.omega0 = {e, w, 0.5 * e};
        else free_.omega0 = {e, 0.0, w};
        free_.q0 = alignWithVertical({free_.inertia.x * free_.omega0.x, free_.inertia.y * free_.omega0.y, free_.inertia.z * free_.omega0.z});
        inertia_ = free_.inertia;
        torque_ = {};
        initial_ = {free_.q0, free_.omega0};
        bodyLines_ = boxLines(1.5 * 0.8, 1.0 * 0.8, 0.5 * 0.8);
        axisLines_[0] = {Vec3{0.0, 0.0, 0.0}, Vec3{1.9, 0.0, 0.0}};   // les trois axes du corps, en couleur
        axisLines_[1] = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 1.3, 0.0}};
        axisLines_[2] = {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.9}};
        duration_ = 20.0;
    } else {
        top_ = HeavyTopProblem{};
        const double spin = topMode_ == kSleeping ? sleepFactor_ * top_.sleepingCriticalSpin() : topSpin_;
        const double tilt = topMode_ == kSleeping ? 0.03 : topTilt_;
        top_.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, tilt);
        top_.omega0 = {0.0, 0.0, spin};
        if (topMode_ == kSteadySlow || topMode_ == kSteadyFast) {
            steady_ = top_.steadyPrecession(tilt, spin, topMode_ == kSteadySlow);
            if (steady_.valid) top_.startSteady(steady_);
            else steadyFallback_ = true;                        // spin trop faible : la toupie part en nutation
        }
        inertia_ = top_.inertia;
        torque_ = top_.torque();
        initial_ = {top_.q0, top_.omega0};
        bodyLines_ = topLines(3.0, 1.2, 2.0);
        duration_ = (topMode_ == kSteadySlow) ? 14.0 : 10.0;
    }
    rhs_ = rotationRhs(inertia_, torque_);
    yRef_ = packRotation(initial_);
    energy0_ = energyOf(initial_);
    lSpace0_ = scenario_ == kTop ? Vec3{} : free_.angularMomentumSpace(initial_);
    momentum0_ = scenario_ == kTop ? top_.verticalMomentum(initial_) : lSpace0_.norm();
    omegaScale_ = std::max(initial_.omega.norm(), 1e-9);
    convergenceTime_ = scenario_ == kTop ? 2.0 : 6.3;

    for (int m = 0; m < kMethodCount; ++m) {
        Run& r = runs_[m];
        switch (m) {
            case kEuler: r.integrator = std::make_unique<EulerRotation>(); break;
            case kRK4: r.integrator = std::make_unique<RK4Rotation>(true); break;
            case kRK4Raw: r.integrator = std::make_unique<RK4Rotation>(false); break;
            case kLie: r.integrator = std::make_unique<LieHeunRotation>(); break;
            default: r.integrator = std::make_unique<SplittingRotation>(); break;
        }
        r.state = initial_;
        r.diverged = false;
        r.error.clear();
        r.energy.clear();
        r.momentum.clear();
        r.norm.clear();
        r.trail.clear();
    }
    for (Series& s : refW_) s.clear();
    refTheta_.clear();
    refCos_.clear();
    refTrail_.clear();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void RigidBodyModule::step(double h) {
    for (Run& r : runs_) {
        if (r.diverged) continue;
        r.integrator->step(inertia_, torque_, r.state, h);
        if (!std::isfinite(r.state.q.w) || !std::isfinite(r.state.omega.x) || r.state.omega.norm() > 1e6 || r.state.q.norm() > 1e6) r.diverged = true;
    }
    if (scenario_ != kSymmetric) advance(ref_, rhs_, clock_.time, yRef_, h, 1000);
}

void RigidBodyModule::sample() {
    const double t = clock_.time;
    const RotationState ref = referenceState();
    if (scenario_ == kTop) {
        const double u = top_.cosTheta(ref.q.normalized());
        refCos_.add(t, u);
        refTheta_.add(t, std::acos(std::clamp(u, -1.0, 1.0)) * kDeg);
    } else {
        refW_[0].add(t, ref.omega.x);
        refW_[1].add(t, ref.omega.y);
        refW_[2].add(t, ref.omega.z);
    }
    const Vec3 tip = axisTip(scenario_);
    refTrail_.push_back(vertexAt(toScene(sandwich(ref.q, tip)), kBlue, 0.6f));
    if (refTrail_.size() > kTrailLength) refTrail_.erase(refTrail_.begin(), refTrail_.begin() + 300);

    for (int m = 0; m < kMethodCount; ++m) {
        Run& r = runs_[m];
        if (r.diverged) continue;
        r.error.add(t, std::max(rotationDistance(r.state, ref), kFloor));
        r.energy.add(t, std::max(std::abs(energyOf(r.state) - energy0_) / std::abs(energy0_), kFloor));
        r.momentum.add(t, std::max(momentumDrift(r.state), kFloor));
        r.norm.add(t, std::max(std::abs(r.state.q.norm() - 1.0), kFloor));
        r.trail.push_back(vertexAt(toScene(sandwich(r.state.q, tip)), kMethodColors[m], 0.6f));
        if (r.trail.size() > kTrailLength) r.trail.erase(r.trail.begin(), r.trail.begin() + 300);
    }
}

void RigidBodyModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;
    clock_.advance(frameSeconds, timeScale_, dt_, duration_, [&](double h) { step(h); }, [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur des intégrateurs à pas fixe à `convergenceTime_`, pour différents pas, contre la référence (formule exacte ou RK45 serré).
void RigidBodyModule::computeConvergence() {
    static const int kStepCounts[] = {50, 100, 200, 400, 800, 1600, 3200};
    static const int kMethods[kFixedMethods] = {kEuler, kRK4, kLie, kSplitting};
    RotationState ref;
    if (scenario_ == kSymmetric) ref = free_.exactSymmetric(convergenceTime_);
    else if (scenario_ == kAsymmetric) ref = free_.reference(convergenceTime_);
    else ref = top_.reference(convergenceTime_);

    for (int i = 0; i < kFixedMethods; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        std::unique_ptr<RotationIntegrator> integrator;
        switch (kMethods[i]) {
            case kEuler: integrator = std::make_unique<EulerRotation>(); break;
            case kRK4: integrator = std::make_unique<RK4Rotation>(true); break;
            case kLie: integrator = std::make_unique<LieHeunRotation>(); break;
            default: integrator = std::make_unique<SplittingRotation>(); break;
        }
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            const double e = rotationDistance(integrateRotation(*integrator, inertia_, torque_, initial_, convergenceTime_, n), ref);
            if (!std::isfinite(e) || e <= 1e-13) continue;
            const double dt = convergenceTime_ / n;
            c.x.push_back(dt);
            c.y.push_back(e);
            if (e < 0.3) {   // la pente ne se lit que dans le régime asymptotique (erreur petite)
                const double lx = std::log10(dt), ly = std::log10(e);
                sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
                ++m;
            }
        }
        c.slope = kNaN;
        if (m >= 3) {
            const double dm = static_cast<double>(m);
            c.slope = (dm * sxy - sx * sy) / (dm * sxx - sx * sx);
        }
    }
    convergenceDirty_ = false;
}

// Ce que la théorie prévoit, en mots.
std::string RigidBodyModule::verdict() const {
    if (scenario_ == kSymmetric) {
        const double i1 = inertia_.x, i3 = inertia_.z;
        const Vec3 lb{i1 * free_.omega0.x, i1 * free_.omega0.y, i3 * free_.omega0.z};
        const double l = lb.norm(), omegaRel = (i3 - i1) * free_.omega0.z / i1, cosT = i3 * free_.omega0.z / l;
        return strf("L'axe tourne autour de L à %.2f rad/s (un tour en %.2f s), à %.1f° de L ; (ω1, ω2) tourne dans le corps à Ω = %.2f rad/s.", l / i1,
                    2.0 * kPi * i1 / l, std::acos(std::clamp(cosT, -1.0, 1.0)) * kDeg, omegaRel);
    }
    if (scenario_ == kAsymmetric) {
        if (axis_ != 1) return strf("Rotation STABLE autour de l'axe %s : la perturbation de %.3f reste petite.", axis_ == 0 ? "long" : "court", perturb_);
        const double lambda = free_.omega0.y * std::sqrt((inertia_.z - inertia_.y) * (inertia_.y - inertia_.x) / (inertia_.x * inertia_.z));
        return strf("Rotation INSTABLE autour de l'axe du milieu : la perturbation croît comme e^(λ t) avec λ = %.3f /s (×%.1f à chaque seconde).", lambda, std::exp(lambda));
    }
    const double critical = top_.sleepingCriticalSpin();
    if (topMode_ == kSleeping) {
        const double spin = sleepFactor_ * critical;
        if (spin >= critical) return strf("Toupie endormie STABLE : spin %.1f rad/s au-dessus du spin critique %.1f rad/s.", spin, critical);
        const double gamma = std::sqrt(4.0 * inertia_.x * top_.mass * top_.gravity * top_.lever - inertia_.z * inertia_.z * spin * spin) / (2.0 * inertia_.x);
        return strf("Toupie endormie INSTABLE : spin %.1f rad/s sous le spin critique %.1f rad/s ; l'inclinaison croît comme e^(γ t), γ = %.2f /s.", spin, critical, gamma);
    }
    if (steadyFallback_) return "Spin trop faible pour une précession régulière à cette inclinaison : la toupie s'incline et oscille (nutation).";
    if (steady_.valid) return strf("Précession régulière à φ' = %.3f rad/s (un tour de l'axe en %.1f s), inclinaison constante %.0f°.", steady_.phiDot, 2.0 * kPi / steady_.phiDot, steady_.theta * kDeg);
    const NutationRange n = top_.nutation();
    if (!n.valid) return "Mouvement sans nutation périodique simple.";
    return strf("Nutation : cos θ oscille entre %.3f et %.3f (θ de %.1f° à %.1f°) avec la période T = %.3f s.", n.u2, n.u1, std::acos(std::clamp(n.u2, -1.0, 1.0)) * kDeg,
                std::acos(std::clamp(n.u1, -1.0, 1.0)) * kDeg, top_.nutationPeriod());
}

// --------------------------------- UI ----------------------------------

void RigidBodyModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool etudiant = atLeast(level, Level::Etudiant);
    bool changed = false;

    ImGui::SeparatorText("Scénario");
    int s = scenario_;
    for (int i = 0; i < 3; ++i) ImGui::RadioButton(kScenarioNames[i], &s, i);
    if (s != scenario_) {
        scenario_ = static_cast<Scenario>(s);
        dt_ = scenario_ == kTop ? 0.004 : 0.01;
        reset();
        return;
    }

    if (interesse && scenario_ != kTop) {   // pour la toupie les « cas types » sont les quatre mouvements, plus bas
        ImGui::TextDisabled("Cas types");
        if (scenario_ == kSymmetric) {
            if (ImGui::SmallButton("Cigare")) { shapeRatio_ = 0.3; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Sphère")) { shapeRatio_ = 1.0; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Disque")) { shapeRatio_ = 1.9; changed = true; }
        } else if (scenario_ == kAsymmetric) {
            if (ImGui::SmallButton("Axe long")) { axis_ = 0; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Axe du milieu (instable)")) { axis_ = 1; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Axe court")) { axis_ = 2; changed = true; }
        }
    }

    ImGui::SeparatorText("Paramètres");
    pushSliderWidth();
    if (scenario_ == kSymmetric) {
        if (interesse) changed |= sliderD("Forme (cigare ↔ disque)", &shapeRatio_, 0.2, 1.95, "%.2f");
        changed |= sliderD(interesse ? "Spin ω3 (rad/s)" : "Rotation sur lui-même", &spin_, 0.5, 6.0, interesse ? "%.2f" : "");
        changed |= sliderD(interesse ? "Écart de l'axe (rad/s)" : "Lancer de travers", &wobble_, 0.0, 2.0, interesse ? "%.2f" : "");
    } else if (scenario_ == kAsymmetric) {
        ImGui::TextDisabled("Axe de rotation");
        changed |= ImGui::RadioButton("Axe long", &axis_, 0);
        changed |= ImGui::RadioButton("Axe du milieu", &axis_, 1);
        changed |= ImGui::RadioButton("Axe court", &axis_, 2);
        changed |= sliderD(interesse ? "Perturbation (rad/s)" : "Petit écart au départ", &perturb_, 1e-3, 0.5, interesse ? "%.3f" : "", ImGuiSliderFlags_Logarithmic);
    } else {
        if (interesse) {
            ImGui::TextDisabled("Mouvement");
            int mode = topMode_;
            changed |= ImGui::RadioButton("Précession lente", &mode, kSteadySlow);
            changed |= ImGui::RadioButton("Précession rapide", &mode, kSteadyFast);
            changed |= ImGui::RadioButton("Nutation", &mode, kNutation);
            changed |= ImGui::RadioButton("Toupie endormie", &mode, kSleeping);
            topMode_ = static_cast<TopMode>(mode);
        }
        if (topMode_ == kSleeping) {
            changed |= sliderD(interesse ? "Spin / spin critique" : "Vitesse de rotation", &sleepFactor_, 0.3, 2.0, interesse ? "%.2f" : "");
        } else {
            changed |= sliderD(interesse ? "Spin ω3 (rad/s)" : "Vitesse de rotation", &topSpin_, 6.0, 40.0, interesse ? "%.1f" : "");
            changed |= sliderD(interesse ? "Inclinaison θ (rad)" : "Inclinaison", &topTilt_, 0.2, 1.4, interesse ? "%.2f" : "");
        }
    }
    popSliderWidth();

    if (college) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
        popSliderWidth();
    }
    if (etudiant) {
        ImGui::TextDisabled("Méthodes affichées");
        for (int m = 0; m < kMethodCount; ++m) {
            ImGui::PushStyleColor(ImGuiCol_CheckMark, toImVec4(kMethodColors[m]));
            ImGui::Checkbox(label(level, m).c_str(), &show_[m]);
            ImGui::PopStyleColor();
        }
    }
    if (changed) reset();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !clock_.finished ? "Pause" : "Lecture")) {
        if (clock_.finished) reset();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) reset();
    pushSliderWidth();
    if (interesse) sliderD("Vitesse du temps", &timeScale_, 0.05, 4.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();
    if (college) ImGui::Text("t = %.2f / %.1f s", clock_.time, duration_);
    for (int m = 0; m < kMethodCount; ++m) {
        if (isShown(level, m) && runs_[m].diverged) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.40f, 1.0f));
            ImGui::TextWrapped("%s : le calcul a divergé (pas trop grand pour ce schéma).", label(level, m).c_str());
            ImGui::PopStyleColor();
        }
    }
    if (interesse) wrapped("Théorie : " + verdict());
    if (interesse) wrapped("Bleu : trajectoire exacte. Flèche jaune : moment cinétique. Flèche verte : rotation instantanée.", true);
}

void RigidBodyModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped(scenario_ == kTop ? "Même penchée, une toupie qui tourne vite garde son énergie et son « élan de rotation » vertical : c'est ce qui "
                                                "la fait tourner autour de la verticale au lieu de tomber."
                                              : "Quand rien ne pousse l'objet, son énergie de rotation et son « élan de rotation » (une flèche fixe dans l'espace) ne "
                                                "changent jamais. Pourtant l'objet, lui, ne cesse de tourner dans tous les sens.");
        return;
    }
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    const bool top = scenario_ == kTop;

    ImGui::SeparatorText("Théorie");
    if (scenario_ == kSymmetric) {
        const double l = std::sqrt(inertia_.x * inertia_.x * (free_.omega0.x * free_.omega0.x + free_.omega0.y * free_.omega0.y) +
                                   inertia_.z * inertia_.z * free_.omega0.z * free_.omega0.z);
        wrapped(strf("I1 = I2 = %.2f, I3 = %.2f ; |L| = %.3f ; E = %.3f", inertia_.x, inertia_.z, l, energy0_));
    } else if (scenario_ == kAsymmetric) {
        wrapped(strf("Boîte 3 × 2 × 1 : I1 = %.3f < I2 = %.3f < I3 = %.3f ; |L| = %.3f ; E = %.3f", inertia_.x, inertia_.y, inertia_.z, momentum0_, energy0_));
        if (etudiant && axis_ != 1) {
            wrapped(strf("Solution de Jacobi : m = k² = %.4f, période de ω(t) = %.4f s.", free_.asymmetricParameter(), free_.asymmetricPeriod()));
        } else if (etudiant) {
            wrapped("Départ (ε, ω2, ε/2) : pas de formule de Jacobi directe (il faudrait l'intégrale elliptique incomplète) : référence RK45 serré.", true);
        }
    } else {
        wrapped(strf("I1 = %.2f, I3 = %.2f, m g l = %.3f ; spin critique de la toupie endormie : %.2f rad/s", inertia_.x, inertia_.z, top_.mass * top_.gravity * top_.lever,
                     top_.sleepingCriticalSpin()));
        if (etudiant) {
            const NutationRange n = top_.nutation();
            if (n.valid) wrapped(strf("Cubique de nutation : u1 = %.4f, u2 = %.4f, u3 = %.3f ; période T = %.4f s", n.u1, n.u2, n.u3, top_.nutationPeriod()));
        }
    }
    wrapped("→ " + verdict());

    const RotationState ref = referenceState();
    ImGui::SeparatorText("Comparaison à l'instant courant");
    std::vector<TableRow> rowsA, rowsB;
    rowsA.push_back({"Référence", kBlue, {"-", strf("%+.1e", (energyOf(ref) - energy0_) / std::abs(energy0_)), strf("%.0e", momentumDrift(ref))}});
    rowsB.push_back({"Référence", kBlue, {strf("%.0e", std::abs(ref.q.norm() - 1.0)), "-", "-"}});
    for (int m = 0; m < kMethodCount; ++m) {
        if (!isShown(level, m)) continue;
        const Run& r = runs_[m];
        const std::string name = etudiant ? shortLabel(m) : label(level, m);
        if (r.diverged) {
            rowsA.push_back({name, kMethodColors[m], {"diverge", "-", "-"}});
            rowsB.push_back({name, kMethodColors[m], {"diverge", "-", "-"}});
            continue;
        }
        const double omegaErr = (r.state.omega - ref.omega).norm();
        const double orientErr = std::sqrt(std::max(rotationDistance(r.state, ref) * rotationDistance(r.state, ref) - omegaErr * omegaErr, 0.0));
        rowsA.push_back({name, kMethodColors[m],
                         {strf("%.1e", rotationDistance(r.state, ref)), strf("%+.1e", (energyOf(r.state) - energy0_) / std::abs(energy0_)), strf("%.1e", momentumDrift(r.state))}});
        rowsB.push_back({name, kMethodColors[m], {strf("%.1e", std::abs(r.state.q.norm() - 1.0)), strf("%.1e", omegaErr), strf("%.1e", orientErr)}});
    }
    drawResultTable("comparaison", {"écart", "ΔE/E", top ? "ΔLz/Lz" : "ΔL/L"}, rowsA);
    if (etudiant) {
        ImGui::SeparatorText("Norme du quaternion et détail de l'écart");
        drawResultTable("detail", {"||q|−1|", "écart ω", "écart R"}, rowsB);
    }
    wrapped(top ? "écart : distance dans l'espace des phases (matrice de rotation et ω). ΔLz/Lz : dérive du moment cinétique VERTICAL (conservé : le couple de la pesanteur n'a pas de composante verticale)."
                : "écart : distance dans l'espace des phases (matrice de rotation et ω). ΔL/L : dérive du vecteur moment cinétique dans le repère FIXE (constant).",
            true);
    if (lycee && !etudiant) wrapped("Euler (méthode simple) fait gonfler le solide : sa norme |q| croît à chaque pas.", true);
    if (chercheur) wrapped("Le découpage symplectique conserve L (ou L_z) à l'arrondi pour tout pas ; la renormalisation de RK4 rétablit |q| = 1 mais pas L.", true);
}

void RigidBodyModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showError = atLeast(level, Level::College);
    const bool showEnergy = atLeast(level, Level::Lycee);
    const int cols = 1 + (showError ? 1 : 0) + (showEnergy ? 1 : 0);

    auto plotRuns = [&](Series Run::*member) {
        for (int m = 0; m < kMethodCount; ++m) {
            if (!isShown(level, m) || (runs_[m].*member).size() == 0) continue;
            const Series& s = runs_[m].*member;
            ImPlot::PlotLine(plotLabel(level, m).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(kMethodColors[m], s.offset));
        }
    };

    if (!ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) return;
    if (scenario_ == kTop) {
        if (ImPlot::BeginPlot("Inclinaison de la toupie θ(t)")) {
            ImPlot::SetupAxes("t (s)", "θ (degrés)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
            double lo, hi;
            paddedLimits(refTheta_, 4.0, lo, hi);
            ImPlot::SetupAxisLimits(ImAxis_Y1, lo, hi, ImPlotCond_Always);
            if (refTheta_.size() > 0) ImPlot::PlotLine("Référence", refTheta_.x.data(), refTheta_.y.data(), refTheta_.size(), lineSpec(kBlue, refTheta_.offset));
            ImPlot::EndPlot();
        }
    } else if (ImPlot::BeginPlot("Vitesse de rotation autour des axes du corps")) {
        ImPlot::SetupAxes("t (s)", "ω (rad/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        static const char* names[3] = {"ω1 (axe long / 1)", "ω2 (axe 2)", "ω3 (axe 3)"};
        for (int k = 0; k < 3; ++k)
            if (refW_[k].size() > 0) ImPlot::PlotLine(names[k], refW_[k].x.data(), refW_[k].y.data(), refW_[k].size(), lineSpec(kAxisColors[k], refW_[k].offset));
        ImPlot::EndPlot();
    }
    if (showError && ImPlot::BeginPlot("Écart à la référence")) {
        ImPlot::SetupAxes("t (s)", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        plotRuns(&Run::error);
        ImPlot::EndPlot();
    }
    if (showEnergy && ImPlot::BeginPlot("Énergie : |ΔE/E|")) {
        ImPlot::SetupAxes("t (s)", "écart relatif", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        plotRuns(&Run::energy);
        ImPlot::EndPlot();
    }
    ImPlot::EndSubplots();
}

void RigidBodyModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    if (convergenceDirty_) computeConvergence();
    const bool top = scenario_ == kTop;
    const int cols = 4;

    auto markerSpec = [](const float* c) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f);
    };
    auto plotRuns = [&](Series Run::*member) {
        for (int m = 0; m < kMethodCount; ++m) {
            if (!isShown(level, m) || (runs_[m].*member).size() == 0) continue;
            const Series& s = runs_[m].*member;
            ImPlot::PlotLine(plotLabel(level, m).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(kMethodColors[m], s.offset));
        }
    };

    if (!ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) return;
    if (ImPlot::BeginPlot(top ? "Moment vertical : |ΔLz/Lz|" : "Moment cinétique : |ΔL/L|")) {
        ImPlot::SetupAxes("t (s)", "écart relatif", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        plotRuns(&Run::momentum);
        ImPlot::EndPlot();
    }
    if (ImPlot::BeginPlot("Norme du quaternion : ||q| − 1|")) {
        ImPlot::SetupAxes("t (s)", "écart", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        plotRuns(&Run::norm);
        ImPlot::EndPlot();
    }
    if (ImPlot::BeginPlot(strf("Convergence : erreur à t = %.1f s", convergenceTime_).c_str())) {
        ImPlot::SetupAxes("dt (s)", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        static const int kMethods[kFixedMethods] = {kEuler, kRK4, kLie, kSplitting};
        for (int i = 0; i < kFixedMethods; ++i) {
            const Curve& c = convergence_[i];
            if (c.x.size() < 2) continue;
            const std::string base = shortLabel(kMethods[i]);
            const std::string name = std::isfinite(c.slope) ? strf("%s (pente %.2f)", base.c_str(), c.slope) : base;
            ImPlot::PlotLine(name.c_str(), c.x.data(), c.y.data(), static_cast<int>(c.x.size()), markerSpec(kMethodColors[kMethods[i]]));
        }
        ImPlot::EndPlot();
    }
    if (top) {
        if (ImPlot::BeginPlot("Nutation : cos θ entre u1 et u2")) {
            ImPlot::SetupAxes("t (s)", "cos θ", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
            double lo, hi;
            paddedLimits(refCos_, 0.04, lo, hi);
            ImPlot::SetupAxisLimits(ImAxis_Y1, lo, hi, ImPlotCond_Always);
            if (refCos_.size() > 0) ImPlot::PlotLine("Référence", refCos_.x.data(), refCos_.y.data(), refCos_.size(), lineSpec(kBlue, refCos_.offset));
            const NutationRange n = top_.nutation();
            if (n.valid) {
                const double lo = n.u1, hi = n.u2;
                ImPlot::PlotInfLines("u1", &lo, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kLColor), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
                ImPlot::PlotInfLines("u2", &hi, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kLColor), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            }
            ImPlot::EndPlot();
        }
    } else {
        // Polhode : trajectoire de ω dans le repère du corps. Solide symétrique : un cercle dans le plan (ω1, ω2) (ω3 est constante) ; asymétrique :
        // une courbe fermée dans le plan (ω1, ω3).
        const bool sym = scenario_ == kSymmetric;
        const int yIndex = sym ? 1 : 2;
        if (ImPlot::BeginPlot(sym ? "Polhode : ω dans le corps (ω1, ω2)" : "Polhode : ω dans le corps (ω1, ω3)", ImVec2(-1.0f, 0.0f), sym ? ImPlotFlags_Equal : 0)) {
            ImPlot::SetupAxes("ω1 (rad/s)", sym ? "ω2 (rad/s)" : "ω3 (rad/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (refW_[0].size() > 0) ImPlot::PlotLine("Référence", refW_[0].y.data(), refW_[yIndex].y.data(), refW_[0].size(), lineSpec(kBlue, refW_[0].offset));
            ImPlot::EndPlot();
        }
    }
    ImPlot::EndSubplots();
}

// ------------------------------ rendu 3D --------------------------------

void RigidBodyModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = scenario_ == kTop ? 1.4f : 0.8f;
    camera.target[2] = 0.0f;
    camera.distance = scenario_ == kTop ? 6.2f : 9.5f;
    camera.yaw = 0.6f;
    camera.pitch = 0.35f;
}

void RigidBodyModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 1.0f);
    renderer.draw(Primitive::Lines, grid);
    const Level level = ctx.level;

    auto drawWire = [&](const RotationState& s, const float* c, float dim, double scale) {
        std::vector<Vertex> v;
        v.reserve(bodyLines_.size());
        for (const Vec3& p : bodyLines_) v.push_back(vertexAt(toScene(sandwich(s.q, scale * p)), c, dim));
        renderer.draw(Primitive::Lines, v);
    };

    const RotationState ref = referenceState();
    // Référence : un peu plus grande, pour rester visible derrière les solides calculés qui la recouvrent.
    drawWire(ref, kBlue, 1.0f, 1.04);
    if (scenario_ == kAsymmetric) {   // axes du corps en couleur (référence)
        for (int k = 0; k < 3; ++k) {
            std::vector<Vertex> v;
            for (const Vec3& p : axisLines_[k]) v.push_back(vertexAt(toScene(sandwich(ref.q, p)), kAxisColors[k]));
            renderer.draw(Primitive::Lines, v);
        }
    }
    for (int m = 0; m < kMethodCount; ++m) {
        if (!isShown(level, m)) continue;
        const Run& r = runs_[m];
        if (r.diverged) continue;
        if (!r.trail.empty()) renderer.draw(Primitive::LineStrip, r.trail);
        drawWire(r.state, kMethodColors[m], 1.0f, 1.0);
    }
    if (!refTrail_.empty()) renderer.draw(Primitive::LineStrip, refTrail_);

    const float pivot[3] = {0.8f, 0.8f, 0.85f};
    renderer.draw(Primitive::Points, {vertexAt({0.0, 0.0, 0.0}, pivot)}, 7.0f * ctx.uiScale);

    if (atLeast(level, Level::Interesse)) {
        if (scenario_ != kTop) {   // moment cinétique fixe dans l'espace (jaune) et vitesse angulaire instantanée (verte)
            const Vec3 l = lSpace0_.normalized() * 3.6;
            renderer.draw(Primitive::Lines, {vertexAt(toScene({0.0, 0.0, 0.0}), kLColor), vertexAt(toScene(l), kLColor)});
            renderer.draw(Primitive::Points, {vertexAt(toScene(l), kLColor)}, 8.0f * ctx.uiScale);
        }
        const Vec3 w = sandwich(ref.q, ref.omega) * (2.4 / omegaScale_);
        renderer.draw(Primitive::Lines, {vertexAt(toScene({0.0, 0.0, 0.0}), kWColor), vertexAt(toScene(w), kWColor)});
        renderer.draw(Primitive::Points, {vertexAt(toScene(w), kWColor)}, 6.0f * ctx.uiScale);
    }
}

}  // namespace pl
