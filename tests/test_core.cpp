// Validation du socle : maths de base + comparaison du solveur à des solutions analytiques.
// Pas de framework externe : un CHECK minimal suffit pour M0.
#include <cmath>
#include <cstdio>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Mat3.hpp"
#include "physicslab/core/Quaternion.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"
#include "physicslab/core/World.hpp"
#include "physicslab/mechanics/Oscillator.hpp"
#include "physicslab/mechanics/Projectile.hpp"

namespace {

int g_failures = 0;

void check(bool ok, const char* expr, int line) {
    if (!ok) {
        std::fprintf(stderr, "test_core.cpp:%d : ECHEC : %s\n", line, expr);
        ++g_failures;
    }
}
#define CHECK(cond) check((cond), #cond, __LINE__)
#define CHECK_NEAR(a, b, tol) check(std::abs((a) - (b)) <= (tol), #a " ~= " #b, __LINE__)

using namespace pl;

void testVec3() {
    const Vec3 a{1, 2, 3}, b{4, 5, 6};
    CHECK_NEAR(dot(a, b), 32.0, 1e-12);
    const Vec3 c = cross(a, b);  // (-3, 6, -3)
    CHECK_NEAR(c.x, -3.0, 1e-12);
    CHECK_NEAR(c.y, 6.0, 1e-12);
    CHECK_NEAR(c.z, -3.0, 1e-12);
    CHECK_NEAR(dot(c, a), 0.0, 1e-12);  // a x b est orthogonal à a
    CHECK_NEAR(Vec3(3, 4, 0).norm(), 5.0, 1e-12);
    CHECK_NEAR(Vec3(0, 0, 0).normalized().norm(), 0.0, 0.0);  // pas de NaN
}

void testMat3() {
    Mat3 a;
    a.m[0][0] = 2; a.m[0][1] = 1; a.m[0][2] = 0;
    a.m[1][0] = 0; a.m[1][1] = 3; a.m[1][2] = 1;
    a.m[2][0] = 1; a.m[2][1] = 0; a.m[2][2] = 4;
    const Mat3 p = a * a.inverse();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK_NEAR(p.m[i][j], i == j ? 1.0 : 0.0, 1e-12);

    const Vec3 v{1, 2, 3}, w{-2, 0.5, 4};
    const Vec3 s = Mat3::skew(v) * w;  // [v]x w = v x w
    const Vec3 c = cross(v, w);
    CHECK_NEAR(s.x, c.x, 1e-12);
    CHECK_NEAR(s.y, c.y, 1e-12);
    CHECK_NEAR(s.z, c.z, 1e-12);
}

void testQuaternion() {
    const double pi = constants::pi;
    const Quaternion q = Quaternion::fromAxisAngle({0, 0, 1}, 0.5 * pi);  // 90 deg autour de z
    const Vec3 r = q.rotate({1, 0, 0});
    CHECK_NEAR(r.x, 0.0, 1e-12);
    CHECK_NEAR(r.y, 1.0, 1e-12);
    CHECK_NEAR(r.z, 0.0, 1e-12);

    // La matrice de rotation et la rotation par quaternion doivent coïncider.
    const Quaternion q2 = Quaternion::fromAxisAngle({1, 2, 3}, 0.7);
    const Vec3 v{0.3, -1.2, 2.5};
    const Vec3 a = q2.rotate(v), b = q2.toMat3() * v;
    CHECK_NEAR(a.x, b.x, 1e-12);
    CHECK_NEAR(a.y, b.y, 1e-12);
    CHECK_NEAR(a.z, b.z, 1e-12);
    CHECK_NEAR(a.norm(), v.norm(), 1e-12);  // une rotation conserve les longueurs

    // q * q^-1 = identité
    const Quaternion id = q2 * q2.conjugate();
    CHECK_NEAR(id.w, 1.0, 1e-12);
    CHECK_NEAR(id.x, 0.0, 1e-12);
}

// y' = -y, y(0) = 1  =>  y(1) = e^-1. Euler est d'ordre 1 : diviser dt par 2 divise l'erreur par ~2.
double eulerErrorDecay(double dt) {
    ExplicitEuler solver;
    State y{1.0};
    const OdeFunction f = [](double, const State& s, State& d) { d[0] = -s[0]; };
    const int steps = static_cast<int>(std::lround(1.0 / dt));
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += solver.step(f, t, y, dt);
    return std::abs(y[0] - std::exp(-1.0));
}

void testEulerOrder() {
    const double e1 = eulerErrorDecay(0.01);
    const double e2 = eulerErrorDecay(0.005);
    CHECK(e1 > 0.0);
    CHECK_NEAR(e1 / e2, 2.0, 0.05);
}

// Chute libre : solution analytique r(t) = r0 + v0 t + g t^2 / 2.
World makeProjectile(double mass) {
    World w;
    Particle p;
    p.mass = mass;
    p.position = {0.0, 10.0, 0.0};
    p.velocity = {3.0, 8.0, 0.0};
    w.particles.push_back(p);
    return w;
}

double runFreeFall(double dt, double tEnd, double& positionError, double& energyDrift) {
    World w = makeProjectile(2.0);
    const Invariants i0 = w.invariants();
    ExplicitEuler solver;
    const int steps = static_cast<int>(std::lround(tEnd / dt));
    for (int i = 0; i < steps; ++i) w.step(solver, dt);

    const Particle& p = w.particles[0];
    const double t = w.time;
    const Vec3 ref = Vec3{0.0, 10.0, 0.0} + t * Vec3{3.0, 8.0, 0.0} + 0.5 * t * t * w.gravity;
    positionError = (p.position - ref).norm();
    energyDrift = w.invariants().total() - i0.total();
    return t;
}

void testFreeFall() {
    const double m = 2.0, g = constants::g0, tEnd = 1.0;

    double err1 = 0, drift1 = 0, err2 = 0, drift2 = 0;
    const double t1 = runFreeFall(0.001, tEnd, err1, drift1);
    runFreeFall(0.0005, tEnd, err2, drift2);

    // Position : erreur d'ordre 1 en dt (Euler explicite).
    CHECK(err1 < 0.01);
    CHECK_NEAR(err1 / err2, 2.0, 0.05);

    // Dérive d'énergie calculée à la main : dE = (1/2) m g^2 t dt, exactement.
    CHECK_NEAR(drift1, 0.5 * m * g * g * t1 * 0.001, 1e-6);
    CHECK_NEAR(drift1 / drift2, 2.0, 0.01);

    // Quantité de mouvement : p_x constante, p_y = m (v0y - g t) à la précision machine.
    World w = makeProjectile(m);
    ExplicitEuler solver;
    for (int i = 0; i < 1000; ++i) w.step(solver, 0.001);
    const Invariants inv = w.invariants();
    CHECK_NEAR(inv.momentum.x, m * 3.0, 1e-9);
    CHECK_NEAR(inv.momentum.y, m * (8.0 - g * w.time), 1e-9);

    // Moment cinétique L_z = m (x vy - y vx) : varie sous la pesanteur, valeur initiale connue.
    const Invariants i0 = makeProjectile(m).invariants();
    CHECK_NEAR(i0.angularMomentum.z, m * (0.0 * 8.0 - 10.0 * 3.0), 1e-12);
}

// --- M1 : projectile avec frottement ---------------------------------------------------

ProjectileProblem dragProblem() {
    ProjectileProblem p;
    p.position0 = {0.0, 2.0, 0.0};
    p.velocity0 = {8.0, 11.0, 0.0};
    p.mass = 1.5;
    p.linearDrag = 1.2;  // k = b/m = 0.8 1/s
    return p;
}

void testProjectileAnalytic() {
    ProjectileProblem none = dragProblem();
    none.linearDrag = 0.0;

    // Sans frottement : parabole r0 + v0 t + g t^2 / 2.
    const double t = 1.7;
    const Vec3 parabola = none.position0 + t * none.velocity0 + 0.5 * t * t * none.gravity;
    CHECK_NEAR((none.position(t) - parabola).norm(), 0.0, 1e-12);

    // Frottement infime : la forme exacte de phi/psi (0/0) doit rester stable et rejoindre la parabole.
    ProjectileProblem tiny = none;
    tiny.linearDrag = 1e-12;
    CHECK_NEAR((tiny.position(t) - parabola).norm(), 0.0, 1e-9);

    // Les deux branches (Taylor / forme exacte) se raccordent sans saut.
    ProjectileProblem drag = dragProblem();
    const double tSwitch = 1e-3 / (drag.linearDrag / drag.mass);
    CHECK_NEAR((drag.position(tSwitch * 0.999999) - drag.position(tSwitch * 1.000001)).norm(), 0.0, 1e-4);

    // La vitesse est bien la dérivée de la position (différence centrée).
    const double h = 1e-5, t0 = 0.9;
    const Vec3 fd = (drag.position(t0 + h) - drag.position(t0 - h)) / (2.0 * h);
    CHECK_NEAR((fd - drag.velocity(t0)).norm(), 0.0, 1e-7);

    // Vitesse limite : v_y -> -m g / b.
    CHECK_NEAR(drag.velocity(60.0).y, -drag.mass * constants::g0 / drag.linearDrag, 1e-9);

    // Durée de vol sans frottement : (vy + sqrt(vy^2 + 2 g h)) / g.
    const double g = constants::g0, vy = none.velocity0.y, y0 = none.position0.y;
    CHECK_NEAR(none.landingTime(), (vy + std::sqrt(vy * vy + 2.0 * g * y0)) / g, 1e-9);
    CHECK_NEAR(drag.position(drag.landingTime()).y, 0.0, 1e-9);
    CHECK(drag.landingTime() < none.landingTime());  // l'air raccourcit le vol
}

// Rapport d'erreur quand on double le nombre de pas : 2^ordre pour une méthode d'ordre `ordre`.
double errorRatio(Solver& solver, const ProjectileProblem& p, int steps, double tEnd) {
    return integrationError(p, solver, steps, tEnd) / integrationError(p, solver, 2 * steps, tEnd);
}

void testConvergenceOrders() {
    const ProjectileProblem p = dragProblem();
    const double tEnd = 1.5;

    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    CHECK_NEAR(errorRatio(euler, p, 40, tEnd), 2.0, 0.15);       // ordre 1
    CHECK_NEAR(errorRatio(symplectic, p, 40, tEnd), 2.0, 0.15);  // ordre 1
    CHECK_NEAR(errorRatio(verlet, p, 40, tEnd), 4.0, 0.4);       // ordre 2
    CHECK_NEAR(errorRatio(rk4, p, 10, tEnd), 16.0, 2.0);         // ordre 4

    // À pas égal, plus l'ordre est élevé, plus l'erreur est petite.
    const double eE = integrationError(p, euler, 30, tEnd);
    const double eV = integrationError(p, verlet, 30, tEnd);
    const double eR = integrationError(p, rk4, 30, tEnd);
    CHECK(eR < eV);
    CHECK(eV < eE);
}

void testPolynomialExactness() {
    // Accélération constante : la solution est un polynôme de degré 2, Verlet et RK4 l'intègrent exactement.
    ProjectileProblem p = dragProblem();
    p.linearDrag = 0.0;
    VelocityVerlet verlet;
    RK4 rk4;
    SymplecticEuler symplectic;
    CHECK(integrationError(p, verlet, 7, 1.0) < 1e-12);
    CHECK(integrationError(p, rk4, 7, 1.0) < 1e-12);
    CHECK(integrationError(p, symplectic, 7, 1.0) > 1e-3);  // Euler symplectique, lui, ne l'est pas

    // Et l'énergie reste constante pour Verlet.
    World w = p.makeWorld();
    const double e0 = w.invariants().total();
    for (int i = 0; i < 100; ++i) w.step(verlet, 0.05);
    CHECK_NEAR(w.invariants().total(), e0, 1e-10);
}

void testRK45() {
    const ProjectileProblem p = dragProblem();

    auto run = [&](double relTol, double dt, RK45& solver) {
        solver.relTol = relTol;
        solver.absTol = relTol * 1e-2;
        World w = p.makeWorld();
        const double advanced = w.step(solver, dt);
        CHECK_NEAR(advanced, dt, 1e-12);  // World::step couvre bien tout l'intervalle demandé
        return (w.particles[0].position - p.position(w.time)).norm();
    };

    RK45 loose, tight;
    const double eLoose = run(1e-4, 1.5, loose);
    const double eTight = run(1e-10, 1.5, tight);
    CHECK(eTight < eLoose);
    CHECK(eTight < 1e-7);
    CHECK(tight.acceptedSteps() > loose.acceptedSteps());  // tolérance plus fine => plus de pas
    CHECK(tight.evaluations() >= 7 * tight.acceptedSteps());
}

// --- M2 : oscillateur harmonique -------------------------------------------------------

OscillatorProblem makeOscillator(double zeta, double forceAmplitude, double forceOmega) {
    OscillatorProblem p;
    p.mass = 1.3;
    p.stiffness = 12.0;
    p.x0 = 0.7;
    p.v0 = -0.4;
    p.damping = zeta * 2.0 * std::sqrt(p.stiffness * p.mass);
    p.forceAmplitude = forceAmplitude;
    p.forceFrequency = forceOmega;
    return p;
}

// La solution exacte vérifie-t-elle x' = v et m v' = -k x - c v + F0 cos(w t) ? (différences centrées)
void checkOscillatorResidual(const OscillatorProblem& p, double t) {
    const double h = 1e-5;
    double x, v, xa, va, xb, vb;
    p.exact(t, x, v);
    p.exact(t + h, xa, va);
    p.exact(t - h, xb, vb);
    CHECK_NEAR((xa - xb) / (2.0 * h), v, 1e-6);
    const double acceleration = (-p.stiffness * x - p.damping * v + p.forceAmplitude * std::cos(p.forceFrequency * t)) / p.mass;
    CHECK_NEAR((va - vb) / (2.0 * h), acceleration, 1e-6);
}

void testOscillatorExact() {
    const double w0 = makeOscillator(0, 0, 0).omega0();
    const double zetas[] = {0.0, 0.2, 1.0, 2.5};
    for (double z : zetas) {
        for (int forced = 0; forced < 2; ++forced) {
            const OscillatorProblem p = makeOscillator(z, forced ? 3.0 : 0.0, 0.8 * w0);
            double x, v;
            p.exact(0.0, x, v);
            CHECK_NEAR(x, p.x0, 1e-12);  // conditions initiales respectées
            CHECK_NEAR(v, p.v0, 1e-12);
            for (double t : {0.3, 1.1, 2.9}) checkOscillatorResidual(p, t);
        }
    }

    // Résonance exacte sans frottement (solution séculaire) et ses voisines.
    const OscillatorProblem res = makeOscillator(0.0, 2.0, w0);
    for (double t : {0.5, 2.0, 4.0}) checkOscillatorResidual(res, t);
    OscillatorProblem near = res;
    near.forceFrequency = w0 * (1.0 + 1e-7);
    CHECK_NEAR(near.position(3.0), res.position(3.0), 1e-5);  // continuité à travers la résonance

    // Régime critique : continu en zeta (sous-critique, critique, sur-critique).
    const double t = 1.3;
    const double below = makeOscillator(1.0 - 1e-7, 0, 0).position(t);
    const double crit = makeOscillator(1.0, 0, 0).position(t);
    const double above = makeOscillator(1.0 + 1e-7, 0, 0).position(t);
    CHECK_NEAR(below, crit, 1e-5);
    CHECK_NEAR(above, crit, 1e-5);
    CHECK(makeOscillator(1.0 - 1e-7, 0, 0).regime() == DampingRegime::Underdamped);
    CHECK(makeOscillator(1.0, 0, 0).regime() == DampingRegime::Critical);
    CHECK(makeOscillator(1.0 + 1e-7, 0, 0).regime() == DampingRegime::Overdamped);

    // Oscillateur libre : périodique et d'énergie constante.
    const OscillatorProblem freeOsc = makeOscillator(0.0, 0, 0);
    double x, v;
    freeOsc.exact(freeOsc.period(), x, v);
    CHECK_NEAR(x, freeOsc.x0, 1e-12);
    CHECK_NEAR(v, freeOsc.v0, 1e-12);
    freeOsc.exact(7.3, x, v);
    CHECK_NEAR(freeOsc.energy(x, v), freeOsc.energy(freeOsc.x0, freeOsc.v0), 1e-12);

    // Amplitude du régime permanent : le maximum de |x| après le transitoire vaut X(w).
    const OscillatorProblem driven = makeOscillator(0.2, 3.0, 0.8 * w0);
    double peak = 0.0;
    const double tStart = 80.0, span = 2.0 * constants::pi / driven.forceFrequency;
    for (int i = 0; i <= 20000; ++i) peak = std::max(peak, std::abs(driven.position(tStart + span * i / 20000.0)));
    CHECK_NEAR(peak, driven.steadyStateAmplitude(driven.forceFrequency), 1e-6);

    // Le maximum de X(w) est bien en w_res = sqrt(w0^2 - 2 g^2).
    const double wRes = driven.resonanceOmega();
    CHECK(driven.steadyStateAmplitude(wRes) > driven.steadyStateAmplitude(wRes * 0.97));
    CHECK(driven.steadyStateAmplitude(wRes) > driven.steadyStateAmplitude(wRes * 1.03));
}

void testOscillatorConvergence() {
    const OscillatorProblem p = makeOscillator(0.15, 2.0, 2.2);
    const double tEnd = 3.0 * p.period();

    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n) { return oscillatorError(p, s, n, tEnd) / oscillatorError(p, s, 2 * n, tEnd); };
    // Euler n'est asymptotique (rapport 2) qu'à petit pas : à 400 pas le rapport vaut encore 2,17.
    CHECK_NEAR(ratio(euler, 1600), 2.0, 0.15);
    CHECK_NEAR(ratio(symplectic, 1600), 2.0, 0.15);
    CHECK_NEAR(ratio(verlet, 100), 4.0, 0.4);
    CHECK_NEAR(ratio(rk4, 40), 16.0, 2.0);

    // Oscillateur libre observé à 2,7 périodes (instant générique) : Euler symplectique d'ordre 1, Verlet d'ordre 2.
    // Observé à un multiple exact de la période il paraîtrait d'ordre 2 et 4 (sur-convergence), voir oscillatorError.
    OscillatorProblem freeOsc = makeOscillator(0.0, 0.0, 0.0);
    freeOsc.v0 = 0.0;
    auto freeRatio = [&](Solver& s, int n) {
        const double tf = 2.7 * freeOsc.period();
        return oscillatorError(freeOsc, s, n, tf) / oscillatorError(freeOsc, s, 2 * n, tf);
    };
    CHECK_NEAR(freeRatio(verlet, 100), 4.0, 0.4);
    CHECK_NEAR(freeRatio(symplectic, 1600), 2.0, 0.15);
}

// Le coeur de M2 : l'énergie d'un oscillateur idéal selon le solveur (k = 4 pi^2, m = 1 : T = 1 s).
void testOscillatorEnergy() {
    OscillatorProblem p;
    p.mass = 1.0;
    p.stiffness = 4.0 * constants::pi * constants::pi;
    p.x0 = 1.0;
    p.v0 = 0.0;
    const double w0 = p.omega0(), dt = 0.02;  // 50 pas par période
    const int periods = 100, perPeriod = 50, steps = periods * perPeriod;
    const OdeFunction f = p.rhs();
    const double e0 = p.energy(p.x0, p.v0);

    auto run = [&](Solver& solver, double& maxFirst, double& maxLast, double& meanFirst, double& meanLast, double& finalRatio) {
        State y = p.initialState();
        double t = 0.0;
        maxFirst = maxLast = meanFirst = meanLast = 0.0;
        for (int i = 1; i <= steps; ++i) {
            t += advance(solver, f, t, y, dt);
            const double rel = p.energy(y[0], y[1]) / e0;
            if (i <= 10 * perPeriod) maxFirst = std::max(maxFirst, std::abs(rel - 1.0));
            if (i > (periods - 10) * perPeriod) maxLast = std::max(maxLast, std::abs(rel - 1.0));
            if (i <= 10 * perPeriod) meanFirst += rel / (10 * perPeriod);
            if (i > (periods - 10) * perPeriod) meanLast += rel / (10 * perPeriod);
            finalRatio = rel;
        }
    };

    double maxFirst, maxLast, meanFirst, meanLast, finalRatio = 0.0;

    // Euler explicite : chaque pas multiplie l'énergie par exactement (1 + w0^2 dt^2). Instable pour tout dt.
    ExplicitEuler euler;
    run(euler, maxFirst, maxLast, meanFirst, meanLast, finalRatio);
    const double predicted = std::pow(1.0 + w0 * w0 * dt * dt, steps);
    CHECK_NEAR(finalRatio / predicted, 1.0, 1e-9);
    CHECK(finalRatio > 100.0);

    // Euler symplectique : énergie bornée, sans dérive séculaire (l'amplitude des oscillations ne grandit pas).
    SymplecticEuler symplectic;
    run(symplectic, maxFirst, maxLast, meanFirst, meanLast, finalRatio);
    CHECK(maxFirst < 0.15);
    CHECK(maxLast < 1.01 * maxFirst);

    // Verlet : Hamiltonien modifié conservé => énergie bornée (ordre 2), moyenne sur 10 périodes sans dérive
    // (écart très inférieur à l'amplitude des oscillations de l'énergie).
    VelocityVerlet verlet;
    run(verlet, maxFirst, maxLast, meanFirst, meanLast, finalRatio);
    CHECK(maxFirst < 0.02);
    CHECK(maxLast < 1.01 * maxFirst);
    CHECK(std::abs(meanLast - meanFirst) < 0.1 * maxFirst);

    // RK4 : très légère dissipation numérique, de l'ordre de (w0 dt)^6 par pas.
    RK4 rk4;
    run(rk4, maxFirst, maxLast, meanFirst, meanLast, finalRatio);
    CHECK(finalRatio < 1.0);
    CHECK(finalRatio > 0.999);
}

}  // namespace

int main() {
    testVec3();
    testMat3();
    testQuaternion();
    testEulerOrder();
    testFreeFall();
    testProjectileAnalytic();
    testConvergenceOrders();
    testPolynomialExactness();
    testRK45();
    testOscillatorExact();
    testOscillatorConvergence();
    testOscillatorEnergy();

    if (g_failures == 0) {
        std::puts("test_core : OK");
        return 0;
    }
    std::fprintf(stderr, "test_core : %d échec(s)\n", g_failures);
    return 1;
}
