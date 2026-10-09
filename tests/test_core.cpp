// Validation du socle : maths de base + comparaison du solveur à des solutions analytiques.
// Pas de framework externe : un CHECK minimal suffit pour M0.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Events.hpp"
#include "physicslab/core/Mat3.hpp"
#include "physicslab/core/Quaternion.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"
#include "physicslab/core/World.hpp"
#include "physicslab/mechanics/Bounce.hpp"
#include "physicslab/mechanics/Collision.hpp"
#include "physicslab/mechanics/Cradle.hpp"
#include "physicslab/mechanics/DoublePendulum.hpp"
#include "physicslab/mechanics/Friction.hpp"
#include "physicslab/mechanics/Kepler.hpp"
#include "physicslab/mechanics/NBody.hpp"
#include "physicslab/mechanics/Oscillator.hpp"
#include "physicslab/mechanics/Pendulum.hpp"
#include "physicslab/mechanics/Projectile.hpp"
#include "physicslab/mechanics/RigidBody.hpp"

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


// --- M3 : pendule simple et pendule double ----------------------------------------------

PendulumProblem makePendulum(double theta0) {
    PendulumProblem p;
    p.length = 1.7;
    p.mass = 0.8;
    p.gravity = 9.81;
    p.theta0 = theta0;
    return p;
}

void testElliptic() {
    CHECK_NEAR(ellipticK(0.0), constants::pi / 2.0, 1e-15);
    CHECK_NEAR(ellipticK(0.5), 1.685750354812596, 1e-14);                // K(m = k^2 = 1/4)
    CHECK_NEAR(ellipticK(std::sqrt(0.5)), 1.854074677301372, 1e-14);      // K(m = 1/2)
}

void testPendulumExact() {
    for (double th0 : {0.05, 0.6, 1.5, 2.5, 3.05, -1.2}) {
        const PendulumProblem p = makePendulum(th0);
        CHECK(p.hasExactSolution());

        double th, om;
        p.exact(0.0, th, om);
        CHECK_NEAR(th, th0, 1e-12);  // conditions initiales
        CHECK_NEAR(om, 0.0, 1e-12);

        // theta' = omega et omega' = -w0^2 sin(theta) (différences centrées)
        const double h = 1e-5, w02 = p.gravity / p.length;
        for (double t : {0.4, 1.3, 3.1}) {
            double a, aw, b, bw, c, cw;
            p.exact(t, a, aw);
            p.exact(t + h, b, bw);
            p.exact(t - h, c, cw);
            CHECK_NEAR((b - c) / (2.0 * h), aw, 1e-6);
            CHECK_NEAR((bw - cw) / (2.0 * h), -w02 * std::sin(a), 1e-5);
            CHECK_NEAR(p.energy(a, aw), p.energy(th0, 0.0), 1e-10);  // énergie conservée
        }

        // Après une période, retour à l'état initial.
        p.exact(p.period(), th, om);
        CHECK_NEAR(th, th0, 1e-9);
        CHECK_NEAR(om, 0.0, 1e-8);
    }

    // Indépendance de l'implémentation : la série elliptique doit coïncider avec un RK45 à tolérance 1e-13.
    const PendulumProblem big = makePendulum(2.5);
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    State y = big.initialState();
    advance(rk, big.rhs(), 0.0, y, 7.0);
    double th, om;
    big.exact(7.0, th, om);
    CHECK_NEAR(th, y[0], 1e-9);
    CHECK_NEAR(om, y[1], 1e-9);

    // Valeurs connues : à 60 degrés T / T0 = 1,0732 ; aux petits angles T ~ T0 (1 + theta0^2 / 16).
    const PendulumProblem sixty = makePendulum(constants::pi / 3.0);
    CHECK_NEAR(sixty.period() / sixty.smallAnglePeriod(), 1.073182, 1e-5);
    const PendulumProblem tiny = makePendulum(1e-2);
    CHECK_NEAR(tiny.period() / tiny.smallAnglePeriod(), 1.0 + 1e-4 / 16.0, 1e-9);

    PendulumProblem damped = makePendulum(1.0);
    damped.damping = 0.1;
    CHECK(!damped.hasExactSolution());
    PendulumProblem pushed = makePendulum(1.0);
    pushed.angularVelocity0 = 0.5;
    CHECK(!pushed.hasExactSolution());
}

void testPendulumConvergence() {
    const PendulumProblem p = makePendulum(2.0);
    const double tEnd = 2.7 * p.period();

    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n) { return pendulumError(p, s, n, tEnd) / pendulumError(p, s, 2 * n, tEnd); };
    CHECK_NEAR(ratio(euler, 3200), 2.0, 0.15);
    CHECK_NEAR(ratio(symplectic, 3200), 2.0, 0.15);
    CHECK_NEAR(ratio(verlet, 200), 4.0, 0.4);
    // RK4 n'est asymptotique qu'à petit pas sur ce grand angle : 80 pas donnent encore un rapport de 52.
    CHECK_NEAR(ratio(rk4, 2560), 16.0, 2.0);

    // Pendule amorti : pas de solution élémentaire, la référence est le RK45 serré et les ordres sont les mêmes.
    PendulumProblem damped = p;
    damped.damping = 0.3;
    CHECK(!damped.hasExactSolution());
    auto dratio = [&](Solver& s, int n) { return pendulumError(damped, s, n, tEnd) / pendulumError(damped, s, 2 * n, tEnd); };
    CHECK_NEAR(dratio(verlet, 200), 4.0, 0.4);
    CHECK_NEAR(dratio(rk4, 2560), 16.0, 2.5);
}

// Énergie du pendule non linéaire sur 100 périodes : c'est ici que le caractère symplectique de Verlet se voit.
void testPendulumEnergy() {
    const PendulumProblem p = makePendulum(2.5);  // grande amplitude : très non linéaire
    const double T = p.period();
    const int perPeriod = 80, periods = 100, steps = perPeriod * periods;
    const double dt = T / perPeriod;
    const OdeFunction f = p.rhs();
    const double e0 = p.energy(p.theta0, 0.0);

    auto run = [&](Solver& solver, double& maxFirst, double& maxLast, double& finalRatio) {
        State y = p.initialState();
        double t = 0.0;
        maxFirst = maxLast = 0.0;
        for (int i = 1; i <= steps; ++i) {
            t += advance(solver, f, t, y, dt);
            const double rel = p.energy(y[0], y[1]) / e0;
            if (i <= 10 * perPeriod) maxFirst = std::max(maxFirst, std::abs(rel - 1.0));
            if (i > (periods - 10) * perPeriod) maxLast = std::max(maxLast, std::abs(rel - 1.0));
            finalRatio = rel;
        }
    };

    double maxFirst, maxLast, finalRatio = 0.0;

    ExplicitEuler euler;
    run(euler, maxFirst, maxLast, finalRatio);
    CHECK(finalRatio > 10.0);  // l'énergie explose

    SymplecticEuler symplectic;
    run(symplectic, maxFirst, maxLast, finalRatio);
    CHECK(maxLast < 1.05 * maxFirst);  // bornée : l'amplitude des oscillations de E ne grandit pas

    VelocityVerlet verlet;
    run(verlet, maxFirst, maxLast, finalRatio);
    CHECK(maxFirst < 0.05);
    CHECK(maxLast < 1.05 * maxFirst);

    RK4 rk4;
    run(rk4, maxFirst, maxLast, finalRatio);
    CHECK(std::abs(finalRatio - 1.0) < 1e-3);
}

DoublePendulumProblem chaoticProblem() {
    DoublePendulumProblem p;
    p.m1 = 1.0; p.m2 = 1.3; p.l1 = 1.0; p.l2 = 0.8;
    p.theta1 = 2.0; p.theta2 = 2.0;
    return p;
}

void testDoublePendulumEquations() {
    // Cohérence des équations transcrites : E = T + V est constante le long du flot, donc dE/dt = grad(E) . f = 0.
    // (une faute de signe ou de facteur dans rhs() fait échouer ce test)
    const DoublePendulumProblem p = chaoticProblem();
    const OdeFunction f = p.rhs();
    const State samples[] = {{0.7, -0.4, 0.3, 0.2}, {2.0, 2.0, 0.0, 0.0}, {-1.1, 2.6, 1.5, -2.0}, {3.0, 0.1, -0.8, 1.9}};
    for (const State& y : samples) {
        State d(4);
        f(0.0, y, d);
        const double eps = 1e-6;
        State up = y, dn = y;
        for (int i = 0; i < 4; ++i) { up[i] += eps * d[i]; dn[i] -= eps * d[i]; }
        CHECK_NEAR((p.energy(up) - p.energy(dn)) / (2.0 * eps), 0.0, 1e-6);
    }

    // Limite m2 -> 0 : le pendule du haut devient un pendule simple de longueur l1.
    DoublePendulumProblem light = chaoticProblem();
    light.m2 = 1e-12;
    State y = {0.7, -0.4, 0.3, 0.2}, d(4);
    light.rhs()(0.0, y, d);
    CHECK_NEAR(d[2], -(light.gravity / light.l1) * std::sin(0.7), 1e-8);

    // Pendule double au repos en bas : équilibre.
    y = {0.0, 0.0, 0.0, 0.0};
    p.rhs()(0.0, y, d);
    for (double v : d) CHECK_NEAR(v, 0.0, 1e-15);

    // Géométrie : theta = 0 pend verticalement, theta1 = pi/2 horizontalement.
    double x1, y1, x2, y2;
    p.positions({constants::pi / 2.0, 0.0, 0.0, 0.0}, x1, y1, x2, y2);
    CHECK_NEAR(x1, p.l1, 1e-12);
    CHECK_NEAR(y1, 0.0, 1e-12);
    CHECK_NEAR(x2, p.l1, 1e-12);
    CHECK_NEAR(y2, -p.l2, 1e-12);
}

void testDoublePendulumNumerics() {
    const DoublePendulumProblem p = chaoticProblem();

    // La référence RK45 conserve l'énergie sur 20 s de mouvement chaotique.
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    State y = p.initialState();
    const double e0 = p.energy(y);
    double t = 0.0;
    double worst = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += advance(rk, p.rhs(), t, y, 0.01);
        worst = std::max(worst, std::abs(p.energy(y) - e0));
    }
    CHECK(worst < 1e-9 * (p.m1 + p.m2) * p.gravity * (p.l1 + p.l2));

    // Ordres de convergence à t = 2 s (avant que le chaos n'amplifie tout), mesurés sur la distance à la référence.
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n) {
        return doublePendulumError(p, s, n, 2.0) / doublePendulumError(p, s, 2 * n, 2.0);
    };
    CHECK_NEAR(ratio(euler, 3200), 2.0, 0.2);
    CHECK_NEAR(ratio(symplectic, 3200), 2.0, 0.2);
    CHECK_NEAR(ratio(verlet, 200), 4.0, 0.5);
    CHECK_NEAR(ratio(rk4, 100), 16.0, 3.0);
}

// Le chaos en chiffres : deux départs qui diffèrent de 1e-9 rad. Pendule double : écart exponentiel,
// exposant de Lyapunov positif. Pendule simple : écart qui croît lentement (déphasage), pas d'exponentielle.
void testChaosVersusRegular() {
    const double delta = 1e-9, tEnd = 25.0, dt = 0.01;

    // distance entre l'orbite et son jumeau perturbé, échantillonnée tous les dt
    auto twinDistances = [&](const OdeFunction& rhs, State a, State b, auto&& distance, std::vector<double>& ts,
                             std::vector<double>& ds) {
        RK45 s1, s2;
        s1.relTol = s2.relTol = 1e-13;
        s1.absTol = s2.absTol = 1e-15;
        double t = 0.0;
        while (t < tEnd - 1e-9) {
            advance(s1, rhs, t, a, dt);
            advance(s2, rhs, t, b, dt);
            t += dt;
            ts.push_back(t);
            ds.push_back(distance(a, b));
        }
    };

    // pendule double
    {
        const DoublePendulumProblem p = chaoticProblem();
        State a = p.initialState(), b = a;
        b[0] += delta;
        std::vector<double> ts, ds;
        twinDistances(p.rhs(), a, b, [&](const State& u, const State& v) { return p.distance(u, v); }, ts, ds);
        const double growth = *std::max_element(ds.begin(), ds.end()) / delta;
        CHECK(growth > 1e5);  // amplification d'au moins 5 ordres de grandeur
        int used = 0;
        const double lambda = lyapunovExponent(ts, ds, 1e2 * delta, 1e-1, &used);
        CHECK(used >= 50);
        CHECK(lambda > 0.5 && lambda < 5.0);
    }
    // pendule simple (lâché de 2 rad), jumeau décalé de delta
    {
        const PendulumProblem p = makePendulum(2.0);
        State a = p.initialState(), b = a;
        b[0] += delta;
        std::vector<double> ts, ds;
        twinDistances(p.rhs(), a, b,
                      [&](const State& u, const State& v) {
                          const double dTheta = u[0] - v[0], dOmega = (u[1] - v[1]) / p.omega0();
                          return std::sqrt(dTheta * dTheta + dOmega * dOmega);
                      },
                      ts, ds);
        const double growth = *std::max_element(ds.begin(), ds.end()) / delta;
        CHECK(growth < 200.0);  // croissance au plus linéaire en t (déphasage), loin de l'exponentielle
    }
}

void testLyapunovFit() {
    // données synthétiques d(t) = d0 e^{0.7 t}, saturées à 1 : l'ajustement doit retrouver 0.7 sur la fenêtre exponentielle
    std::vector<double> t, d;
    for (int i = 0; i < 4000; ++i) {
        t.push_back(0.01 * i);
        d.push_back(std::min(1e-9 * std::exp(0.7 * 0.01 * i), 1.0));
    }
    int used = 0;
    CHECK_NEAR(lyapunovExponent(t, d, 1e-8, 1e-1, &used), 0.7, 1e-9);
    CHECK(used > 100);
    CHECK(std::isnan(lyapunovExponent(t, d, 2.0, 3.0, &used)));  // aucune donnée dans la fenêtre
    CHECK(used == 0);
}

// --- M4 : gravitation, Kepler à 2 corps ------------------------------------------------

void testKeplerEquation() {
    const double pi = constants::pi;
    // M = E - e sin E sur plusieurs tours (M va de -14,8 à +14,8 rad), jusqu'aux orbites très excentriques.
    for (double e : {0.0, 0.3, 0.9, 0.99, 0.999}) {
        for (int i = -40; i <= 40; ++i) {
            const double M = 0.37 * i;
            int iterations = 0;
            const double E = kepler::solveEccentricAnomaly(M, e, &iterations);
            CHECK_NEAR(E - e * std::sin(E), M, 1e-13);
            CHECK(iterations <= 12);
        }
    }
    CHECK_NEAR(kepler::solveEccentricAnomaly(0.0, 0.9), 0.0, 1e-15);  // périastre
    CHECK_NEAR(kepler::solveEccentricAnomaly(pi, 0.9), pi, 1e-14);    // apoastre
    CHECK_NEAR(kepler::solveEccentricAnomaly(1.234, 0.0), 1.234, 1e-15);  // cercle : E = M
}

void testKeplerExact() {
    const double pi = constants::pi;
    for (double e : {0.0, 0.5, 0.9}) {
        KeplerProblem p;
        p.mu = 1.7;
        p.a = 2.5;
        p.e = e;
        const double T = p.period();
        // 3e loi de Kepler : T^2 mu / a^3 = 4 pi^2
        CHECK_NEAR(T * T * p.mu / (p.a * p.a * p.a), 4.0 * pi * pi, 1e-10);

        // Périastre au départ : r = a (1 - e), vitesse de vis-viva.
        const State s0 = p.initialState();
        CHECK_NEAR(s0[0], p.a * (1.0 - e), 1e-14);
        CHECK_NEAR(s0[1], 0.0, 1e-14);
        CHECK_NEAR(s0[4], std::sqrt(p.mu * (1.0 + e) / (p.a * (1.0 - e))), 1e-13);
        CHECK(p.distance(p.exact(0.0), s0) < 1e-14);

        // Apoastre à T/2 : r = a (1 + e), vitesse opposée, plus lente.
        const State sa = p.exact(0.5 * T);
        CHECK_NEAR(sa[0], -p.a * (1.0 + e), 1e-12);
        CHECK_NEAR(sa[1], 0.0, 1e-12);
        CHECK_NEAR(sa[4], -std::sqrt(p.mu * (1.0 - e) / (p.a * (1.0 + e))), 1e-12);

        // Périodicité : après une période on retrouve l'état initial.
        CHECK(p.distance(p.exact(T), s0) < 1e-12);
        CHECK(p.distance(p.exact(7.0 * T), s0) < 1e-11);

        // L'état exact vérifie bien l'EDO r'' = -mu r/|r|^3 (différences finies centrées).
        const OdeFunction f = p.rhs();
        const double h = 1e-4;
        for (double t : {0.2 * T, 0.45 * T, 0.8 * T}) {
            const State yp = p.exact(t + h), ym = p.exact(t - h), y = p.exact(t);
            State d(6);
            f(t, y, d);
            for (int k = 0; k < 3; ++k) {
                CHECK_NEAR((yp[k] - ym[k]) / (2.0 * h), y[3 + k], 1e-6);       // x' = v
                CHECK_NEAR((yp[3 + k] - ym[3 + k]) / (2.0 * h), d[3 + k], 1e-6);  // v' = a(rhs)
                CHECK_NEAR(d[k], y[3 + k], 1e-15);  // la 1re moitié de rhs est la vitesse
            }
            // Accélération attractive, de module mu/r^2
            const double r = std::sqrt(y[0] * y[0] + y[1] * y[1]);
            CHECK_NEAR(std::sqrt(d[3] * d[3] + d[4] * d[4]), p.mu / (r * r), 1e-13);
            CHECK(d[3] * y[0] + d[4] * y[1] < 0.0);
        }

        // Invariants le long de la solution exacte : E, |L| et le vecteur de Runge-Lenz (|A| = mu e, vers le périastre).
        for (double t : {0.0, 0.13 * T, 0.5 * T, 0.77 * T, 3.3 * T}) {
            const State y = p.exact(t);
            CHECK_NEAR(p.energy(y), -p.mu / (2.0 * p.a), 1e-12);
            const Vec3 L = p.angularMomentum(y);
            CHECK_NEAR(L.z, std::sqrt(p.mu * p.a * (1.0 - e * e)), 1e-12);
            CHECK_NEAR(L.x, 0.0, 1e-13);
            CHECK_NEAR(L.y, 0.0, 1e-13);
            CHECK_NEAR(p.exactAngularMomentum(), L.z, 1e-12);
            const Vec3 A = p.rungeLenz(y);
            CHECK_NEAR(A.x, p.mu * e, 1e-12);
            CHECK_NEAR(A.y, 0.0, 1e-12);
            CHECK_NEAR(A.z, 0.0, 1e-13);
        }
        CHECK_NEAR(p.exactEnergy(), p.energy(s0), 1e-13);
        if (e > 0.0) CHECK_NEAR(p.periapsisAngle(p.exact(0.31 * T)), 0.0, 1e-12);
    }
}

// Ordres de convergence sur l'orbite, mesurés à un instant qui n'est pas un multiple de la période (voir PASSATION).
void testKeplerConvergence() {
    KeplerProblem p;
    p.e = 0.5;
    const double T = p.period();
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n, double tEnd) { return keplerError(p, s, n, tEnd) / keplerError(p, s, 2 * n, tEnd); };
    // Euler est encore pré-asymptotique à 2,7 T (rapport 1,76 à 25600 pas : l'erreur de phase croît trop vite) :
    // on le mesure plus tôt, à 0,35 T. Les autres schémas sont asymptotiques à 2,7 T avec ces nombres de pas.
    CHECK_NEAR(ratio(euler, 800, 0.35 * T), 2.0, 0.15);        // ordre 1
    CHECK_NEAR(ratio(symplectic, 25600, 2.7 * T), 2.0, 0.15);  // ordre 1
    CHECK_NEAR(ratio(verlet, 800, 2.7 * T), 4.0, 0.4);         // ordre 2
    CHECK_NEAR(ratio(rk4, 3200, 2.7 * T), 16.0, 2.0);          // ordre 4

    RK45 rk45;
    rk45.relTol = 1e-10;
    rk45.absTol = 1e-12;
    CHECK(keplerError(p, rk45, 200, 2.7 * T) < 1e-7);
}

// Intègre `periods` orbites à pas fixe T/perPeriod et mesure ce qui se conserve (ou non).
struct KeplerRun {
    double maxFirst = 0.0, maxLast = 0.0;  // |E/E0 - 1| maximal sur les 10 premières / 10 dernières orbites
    double lzDrift = 0.0;                  // |Lz/Lz0 - 1| maximal
    double energyEnd = 0.0;                // énergie finale
    double precession = 0.0;               // précession moyenne par orbite [rad]
};

KeplerRun runKepler(const KeplerProblem& p, Solver& solver, int perPeriod, int periods) {
    KeplerRun run;
    State y = p.initialState();
    const OdeFunction f = p.rhs();
    const double dt = p.period() / perPeriod, e0 = p.energy(y), lz0 = p.angularMomentum(y).z;
    PeriapsisTracker tracker;
    double t = 0.0;
    for (int i = 1; i <= perPeriod * periods; ++i) {
        t += advance(solver, f, t, y, dt);
        tracker.update(t, y);
        const double rel = std::abs(p.energy(y) / e0 - 1.0);
        if (i <= 10 * perPeriod) run.maxFirst = std::max(run.maxFirst, rel);
        if (i > (periods - 10) * perPeriod) run.maxLast = std::max(run.maxLast, rel);
        run.lzDrift = std::max(run.lzDrift, std::abs(p.angularMomentum(y).z / lz0 - 1.0));
    }
    run.energyEnd = p.energy(y);
    run.precession = tracker.precessionPerOrbit();
    return run;
}

void testKeplerEnergyAndPrecession() {
    KeplerProblem p;
    p.e = 0.5;
    const double e0 = p.exactEnergy();
    const int perPeriod = 200, periods = 100;

    // Euler explicite : l'orbite spirale vers l'extérieur. L'énergie et le moment cinétique augmentent à chaque tour.
    {
        ExplicitEuler euler;
        State y = p.initialState();
        const OdeFunction f = p.rhs();
        const double dt = p.period() / 1000.0;
        double t = 0.0, previousE = e0, previousL = p.exactAngularMomentum();
        for (int orbit = 1; orbit <= 10; ++orbit) {
            for (int i = 0; i < 1000; ++i) t += advance(euler, f, t, y, dt);
            CHECK(p.energy(y) > previousE);
            CHECK(p.angularMomentum(y).z > previousL);
            previousE = p.energy(y);
            previousL = p.angularMomentum(y).z;
        }
        CHECK(-p.mu / (2.0 * previousE) > 2.0 * p.a);  // le demi-grand axe a plus que doublé en 10 tours
    }
    {   // et à 200 pas par période l'astre finit par s'échapper (E > 0 : orbite ouverte)
        ExplicitEuler euler;
        CHECK(runKepler(p, euler, perPeriod, periods).energyEnd > 0.0);
    }

    // Euler symplectique et Verlet : énergie bornée sans dérive, moment cinétique conservé à l'arrondi (symétrie de
    // la force centrale), mais précession rétrograde en dt^2 que la théorie du hamiltonien modifié prédit.
    const double dt = p.period() / perPeriod;
    const double predicted = p.verletPrecessionPerOrbit(dt);
    CHECK(predicted < 0.0);
    {
        SymplecticEuler symplectic;
        const KeplerRun r = runKepler(p, symplectic, perPeriod, periods);
        CHECK(r.maxLast < 1.05 * r.maxFirst);
        CHECK(r.maxFirst < 0.1);
        CHECK(r.lzDrift < 1e-9);
        CHECK_NEAR(r.precession / predicted, 1.0, 0.02);
    }
    {
        VelocityVerlet verlet;
        const KeplerRun r = runKepler(p, verlet, perPeriod, periods);
        CHECK(r.maxLast < 1.05 * r.maxFirst);
        CHECK(r.maxFirst < 0.01);
        CHECK(r.lzDrift < 1e-9);
        CHECK_NEAR(r.precession / predicted, 1.0, 0.02);

        // En dt^2 : pas divisé par 2 => précession divisée par 4 ; pas divisé par 4 => par 16.
        VelocityVerlet coarse, fine;
        const double p100 = runKepler(p, coarse, 100, 60).precession, p400 = runKepler(p, fine, 400, 60).precession;
        CHECK_NEAR(p100 / p400, 16.0, 1.0);
    }

    // RK4 : précession ~ 1000 fois plus faible, mais l'énergie dérive lentement (dissipation) au lieu de rester bornée.
    {
        RK4 rk4;
        const KeplerRun r = runKepler(p, rk4, perPeriod, periods);
        CHECK(r.energyEnd < e0);                    // E devient plus négative : l'orbite se resserre
        CHECK(r.maxLast > 2.0 * r.maxFirst);        // la dérive grandit
        CHECK(std::abs(r.precession) < 0.01 * std::abs(predicted));
    }
}

// La formule de précession du hamiltonien modifié est valable pour toute excentricité (mesuré : écart <= 0,2 %
// de e = 0,1 à 0,9) et indépendante des unités (seul compte le nombre de pas par période).
void testKeplerPrecessionTheory() {
    struct Case { double e, a, mu; int perPeriod; };
    for (const Case& c : {Case{0.1, 1.0, 1.0, 200}, Case{0.2, 1.0, 1.0, 200}, Case{0.7, 1.0, 1.0, 800},
                          Case{0.9, 1.0, 1.0, 3200}, Case{0.5, 3.7, 2.3, 400}}) {
        KeplerProblem p;
        p.e = c.e; p.a = c.a; p.mu = c.mu;
        VelocityVerlet verlet;
        const KeplerRun r = runKepler(p, verlet, c.perPeriod, 40);
        CHECK_NEAR(r.precession / p.verletPrecessionPerOrbit(p.period() / c.perPeriod), 1.0, 0.01);
    }
}

// Forte excentricité : les pas fixes ne résolvent pas le périastre, le pas adaptatif si.
void testKeplerAdaptiveStep() {
    KeplerProblem p;
    p.e = 0.9;
    const double T = p.period();

    // RK4 à 200 pas par période : la vitesse au périastre est 4,4 fois la vitesse circulaire, le pas est trop grand.
    RK4 rk4;
    CHECK(runKepler(p, rk4, 200, 2).maxFirst > 0.1);

    // RK45 : on observe lui-même ses pas. Au périastre (r = 0,1) ils sont ~ 300 fois plus courts qu'à l'apoastre (r = 1,9).
    RK45 rk45;
    rk45.relTol = 1e-10;
    rk45.absTol = 1e-12;
    State y = p.initialState();
    const OdeFunction f = p.rhs();
    double t = 0.0, hPeri = 1e9, hApo = 0.0;
    const double tEnd = 2.0 * T, hMax = T / 10.0;
    int steps = 0;
    while (t < tEnd - 1e-14) {
        const double r = std::sqrt(y[0] * y[0] + y[1] * y[1]);
        const double h = rk45.step(f, t, y, std::min(hMax, tEnd - t));
        t += h;
        ++steps;
        if (r < 0.15) hPeri = std::min(hPeri, h);
        if (r > 1.7) hApo = std::max(hApo, h);
    }
    CHECK(hPeri < 0.02 * hApo);
    CHECK(steps < 1500);
    CHECK(p.distance(y, p.exact(t)) < 1e-5);
}

// --- M4b : gravitation à N corps -----------------------------------------------------

void testNBodyAccelerations() {
    // 2 corps, sans adoucissement : a1 = G m2 / r^2 vers le corps 2, a2 = G m1 / r^2 vers le corps 1 (r = 5, direction 3-4-5).
    {
        const double pos[6] = {0, 0, 0, 3, 4, 0};
        const double m[2] = {2.0, 7.0};
        double a[6];
        nbody::accelerations(pos, m, 2, 1.5, 0.0, a);
        const double g1 = 1.5 * 7.0 / 25.0, g2 = 1.5 * 2.0 / 25.0;
        CHECK_NEAR(a[0], g1 * 0.6, 1e-15);
        CHECK_NEAR(a[1], g1 * 0.8, 1e-15);
        CHECK_NEAR(a[2], 0.0, 1e-15);
        CHECK_NEAR(a[3], -g2 * 0.6, 1e-15);
        CHECK_NEAR(a[4], -g2 * 0.8, 1e-15);
        CHECK_NEAR(a[5], 0.0, 1e-15);
        CHECK_NEAR(nbody::potentialEnergy(pos, m, 2, 1.5, 0.0), -1.5 * 2.0 * 7.0 / 5.0, 1e-14);
    }
    // Adoucissement de Plummer : a = G m r / (r^2 + eps^2)^(3/2), U = -G m1 m2 / sqrt(r^2 + eps^2) ; fini même à r = 0.
    {
        const double pos[6] = {0, 0, 0, 3, 0, 0};
        const double m[2] = {1.0, 1.0};
        double a[6];
        nbody::accelerations(pos, m, 2, 1.0, 2.0, a);
        CHECK_NEAR(a[0], 3.0 / std::pow(13.0, 1.5), 1e-15);
        CHECK_NEAR(nbody::potentialEnergy(pos, m, 2, 1.0, 2.0), -1.0 / std::sqrt(13.0), 1e-15);
        const double same[6] = {1, 1, 1, 1, 1, 1};
        nbody::accelerations(same, m, 2, 1.0, 0.1, a);
        for (double v : a) CHECK(v == 0.0);
    }
    // Amas de 7 corps : 3e loi de Newton (somme des forces nulle) et force = -gradient de l'énergie potentielle.
    const NBodyProblem c = NBodyProblem::randomCluster(7, 42, 1.0, 0.1);
    CHECK(c.count() == 7);
    State y = c.initialState();
    std::vector<double> acc(21);
    nbody::accelerations(y.data(), c.mass.data(), 7, c.G, c.softening, acc.data());
    Vec3 sum;
    for (int i = 0; i < 7; ++i) sum += c.mass[i] * Vec3{acc[3 * i], acc[3 * i + 1], acc[3 * i + 2]};
    CHECK_NEAR(sum.norm(), 0.0, 1e-14);

    const double h = 1e-5;
    for (int i = 0; i < 7; ++i) {
        for (int k = 0; k < 3; ++k) {
            State plus = y, minus = y;
            plus[3 * i + k] += h;
            minus[3 * i + k] -= h;
            const double dU = (nbody::potentialEnergy(plus.data(), c.mass.data(), 7, c.G, c.softening) -
                               nbody::potentialEnergy(minus.data(), c.mass.data(), 7, c.G, c.softening)) / (2.0 * h);
            CHECK_NEAR(acc[3 * i + k], -dU / c.mass[i], 1e-7);
        }
    }
}

void testNBodyDiagnostics() {
    // Grandeurs de l'amas : cohérence avec les définitions, centre de masse et impulsion nuls, viriel 2T = -U.
    const NBodyProblem c = NBodyProblem::randomCluster(7, 42, 1.0, 0.1);
    const State y = c.initialState();
    CHECK_NEAR(c.totalMass(), 1.0, 1e-14);
    CHECK_NEAR(c.centerOfMass(y).norm(), 0.0, 1e-14);
    CHECK_NEAR(c.momentum(y).norm(), 0.0, 1e-14);
    CHECK_NEAR(2.0 * c.kineticEnergy(y) / (-c.potentialEnergy(y)), 1.0, 1e-12);
    CHECK_NEAR(c.energy(y), c.kineticEnergy(y) + c.potentialEnergy(y), 0.0);
    CHECK(c.kineticEnergy(y) > 0.0 && c.potentialEnergy(y) < 0.0);
    // même graine : même amas ; graine différente : autre amas
    const NBodyProblem same = NBodyProblem::randomCluster(7, 42, 1.0, 0.1), other = NBodyProblem::randomCluster(7, 43, 1.0, 0.1);
    CHECK(same.initialState() == y);
    CHECK(other.initialState() != y);
    // Moment cinétique d'une particule isolée : m r x v
    NBodyProblem one;
    one.mass = {2.0};
    one.position = {{1, 2, 3}};
    one.velocity = {{-1, 0.5, 2}};
    const State s = one.initialState();
    const Vec3 L = one.angularMomentum(s), expected = 2.0 * cross({1, 2, 3}, {-1, 0.5, 2});
    CHECK_NEAR((L - expected).norm(), 0.0, 1e-15);
    CHECK_NEAR(one.kineticEnergy(s), 0.5 * 2.0 * (1 + 0.25 + 4), 1e-15);
    CHECK_NEAR(one.potentialEnergy(s), 0.0, 0.0);
}

void testLagrangeTriangle() {
    const double s = 2.0, m = 0.7, G = 1.3;
    const NBodyProblem p = NBodyProblem::lagrangeTriangle(s, m, G);
    CHECK(p.count() == 3);
    const double omega = std::sqrt(3.0 * G * m / (s * s * s)), R = s / std::sqrt(3.0), T = 2.0 * constants::pi / omega;
    const State y0 = p.initialState();

    // Géométrie : sommets à R du centre, côtés s, vitesse v = omega z x r.
    for (int i = 0; i < 3; ++i) {
        const Vec3 ri{y0[3 * i], y0[3 * i + 1], y0[3 * i + 2]}, vi{y0[9 + 3 * i], y0[9 + 3 * i + 1], y0[9 + 3 * i + 2]};
        CHECK_NEAR(ri.norm(), R, 1e-14);
        CHECK_NEAR((vi - omega * cross({0, 0, 1}, ri)).norm(), 0.0, 1e-14);
        const int j = (i + 1) % 3;
        CHECK_NEAR((ri - Vec3{y0[3 * j], y0[3 * j + 1], y0[3 * j + 2]}).norm(), s, 1e-14);
    }
    CHECK_NEAR(p.centerOfMass(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.momentum(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.angularMomentum(y0).z, 3.0 * m * R * R * omega, 1e-13);
    CHECK_NEAR(p.energy(y0), -1.5 * G * m * m / s, 1e-13);  // 3 (1/2 m R^2 omega^2) - 3 G m^2 / s

    // Les accélérations valent -omega^2 r (force centripète) : la rotation rigide est bien solution.
    State f(18);
    p.rhs()(0.0, y0, f);
    for (int i = 0; i < 9; ++i) {
        CHECK_NEAR(f[9 + i], -omega * omega * y0[i], 1e-13);
        CHECK_NEAR(f[i], y0[9 + i], 0.0);
    }

    // Solution exacte : rotation d'angle omega t de tout l'état. RK4 à 800 pas/période la retrouve.
    RK4 rk4;
    const OdeFunction rhs = p.rhs();
    State y = y0;
    const double tEnd = 1.37 * T, dt = T / 800.0;
    double t = 0.0;
    while (t < tEnd - 1e-12) t += advance(rk4, rhs, t, y, std::min(dt, tEnd - t));
    const double c = std::cos(omega * t), sn = std::sin(omega * t);
    State exact = y0;
    for (int block = 0; block < 2; ++block) {  // bloc 0 : positions, bloc 1 : vitesses
        for (int body = 0; body < 3; ++body) {
            const int k = 9 * block + 3 * body;
            exact[k] = c * y0[k] - sn * y0[k + 1];
            exact[k + 1] = sn * y0[k] + c * y0[k + 1];
        }
    }
    CHECK(p.distance(y, exact) < 1e-8);  // mesuré : 1,4e-9 à 800 pas/période

    // Ordre 4 : doubler le nombre de pas divise l'erreur par 16 (mesuré 16,6 entre 800 et 1600 pas par période).
    auto rotationError = [&](int perPeriod) {
        RK4 solver;
        State z = y0;
        const double end = 1.37 * T, h = T / perPeriod;
        double time = 0.0;
        while (time < end - 1e-12) time += advance(solver, rhs, time, z, std::min(h, end - time));
        const double cc = std::cos(omega * time), ss = std::sin(omega * time);
        State ex = y0;
        for (int block = 0; block < 2; ++block)
            for (int body = 0; body < 3; ++body) {
                const int k = 9 * block + 3 * body;
                ex[k] = cc * y0[k] - ss * y0[k + 1];
                ex[k + 1] = ss * y0[k] + cc * y0[k + 1];
            }
        return p.distance(z, ex);
    };
    CHECK_NEAR(rotationError(800) / rotationError(1600), 16.0, 2.0);
}

void testTwoBodyMatchesKepler() {
    const double m1 = 1.0, m2 = 0.3, a = 1.7, e = 0.6, G = 2.0;
    const NBodyProblem p = NBodyProblem::twoBody(m1, m2, a, e, G);
    KeplerProblem k;
    k.mu = G * (m1 + m2);
    k.a = a;
    k.e = e;
    const double reduced = m1 * m2 / (m1 + m2);

    const State y0 = p.initialState();
    CHECK(p.count() == 2);
    CHECK_NEAR(p.centerOfMass(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.momentum(y0).norm(), 0.0, 1e-14);
    // L et E du système = masse réduite x grandeurs par unité de masse du problème de Kepler
    CHECK_NEAR(p.energy(y0), reduced * k.exactEnergy(), 1e-13);
    CHECK_NEAR(p.angularMomentum(y0).z, reduced * k.exactAngularMomentum(), 1e-13);

    // Le mouvement relatif r2 - r1 suit la solution exacte de Kepler de GM = G (m1 + m2) ; le centre de masse reste au repos.
    RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-14;
    State y = y0;
    const OdeFunction rhs = p.rhs();
    double t = 0.0;
    const double tEnd = 1.3 * k.period();
    for (double target : {0.4 * tEnd, 0.7 * tEnd, tEnd}) {
        t += advance(rk, rhs, t, y, target - t);
        const State rel = k.exact(t);
        for (int c = 0; c < 3; ++c) {
            CHECK_NEAR(y[3 + c] - y[c], rel[c], 1e-9);                // position relative (mesuré : 1,5e-11)
            CHECK_NEAR(y[9 + c] - y[6 + c], rel[3 + c], 1e-9);        // vitesse relative
        }
        CHECK_NEAR(p.centerOfMass(y).norm(), 0.0, 1e-9);
    }
}

void testFigureEight() {
    const NBodyProblem p = NBodyProblem::figureEight();
    CHECK(p.count() == 3);
    CHECK(p.G == 1.0 && p.softening == 0.0);
    const State y0 = p.initialState();
    // Données de Simó : x1 = -x2, x3 = 0, v1 = v2 = -v3/2 ; P = 0, L = 0 par symétrie, E = -1,28714 (T = 1,2129 ; U = -2,5)
    CHECK_NEAR(p.momentum(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.angularMomentum(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.centerOfMass(y0).norm(), 0.0, 1e-14);
    CHECK_NEAR(p.potentialEnergy(y0), -2.5, 1e-7);   // |x1| = 1 (r12 = 2, r13 = r23 = 1)
    CHECK_NEAR(p.energy(y0), -1.287142, 1e-6);

    // Une période après, les trois corps sont revenus : même état à 1e-6 près (mesuré : 8e-8, limité par les 9 chiffres
    // de Simó), et T = 6,32591398 est bien la période : à 1e-4 T d'écart, le retour est 20000 fois moins bon.
    const State yT = p.reference(nbody::kFigureEightPeriod);
    CHECK(p.distance(yT, y0) < 1e-6);
    CHECK(p.distance(p.reference(1.0001 * nbody::kFigureEightPeriod), y0) > 1e-3);
    // Chorégraphie : à T/3 les corps ont tourné d'un cran sur la même courbe, l'ensemble des positions est inchangé.
    const State y3 = p.reference(nbody::kFigureEightPeriod / 3.0);
    std::vector<bool> taken(3, false);
    for (int i = 0; i < 3; ++i) {
        int match = -1;
        for (int j = 0; j < 3; ++j) {
            const double d = std::hypot(y3[3 * i] - y0[3 * j], y3[3 * i + 1] - y0[3 * j + 1]);
            if (d < 1e-5 && !taken[j]) match = j;
        }
        CHECK(match >= 0);
        if (match >= 0) taken[match] = true;
    }
}

// Ordres de convergence sur le huit (mêmes remarques que pour Kepler : Euler et RK4 se mesurent à 0,35 T, Verlet à 2,7 T).
void testNBodyConvergence() {
    const NBodyProblem p = NBodyProblem::figureEight();
    const double T = nbody::kFigureEightPeriod;
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n, double tEnd) { return nbodyError(p, s, n, tEnd) / nbodyError(p, s, 2 * n, tEnd); };
    CHECK_NEAR(ratio(euler, 1600, 0.35 * T), 2.0, 0.15);        // ordre 1
    CHECK_NEAR(ratio(symplectic, 800, 0.35 * T), 2.0, 0.15);    // ordre 1
    CHECK_NEAR(ratio(verlet, 800, 2.7 * T), 4.0, 0.4);          // ordre 2
    CHECK_NEAR(ratio(rk4, 400, 0.35 * T), 16.0, 1.5);           // ordre 4
}

// Ce que chaque schéma conserve sur un amas de 6 corps (t = 20) : l'impulsion par construction des forces opposées
// (même pour Euler), le moment cinétique seulement pour les schémas à structure symplectique (la force est centrale).
void testNBodyConservation() {
    const NBodyProblem c = NBodyProblem::randomCluster(6, 42, 1.0, 0.05);
    const OdeFunction f = c.rhs();
    const int steps = 20000;
    const double dt = 1e-3;

    struct Result { double maxEnergy = 0.0, maxMomentum = 0.0, maxAngular = 0.0; };
    auto run = [&](Solver& solver) {
        Result r;
        State y = c.initialState();
        const double e0 = c.energy(y);
        const Vec3 l0 = c.angularMomentum(y);
        double t = 0.0;
        for (int i = 0; i < steps; ++i) {
            t += advance(solver, f, t, y, dt);
            r.maxEnergy = std::max(r.maxEnergy, std::abs(c.energy(y) / e0 - 1.0));
            r.maxMomentum = std::max(r.maxMomentum, c.momentum(y).norm());
            r.maxAngular = std::max(r.maxAngular, (c.angularMomentum(y) - l0).norm());
        }
        return r;
    };

    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    const Result e = run(euler), s = run(symplectic), v = run(verlet), r = run(rk4);
    for (const Result& x : {e, s, v, r}) CHECK(x.maxMomentum < 1e-13);  // impulsion : tous
    CHECK(e.maxAngular > 1e-4);                                          // Euler dérive
    CHECK(s.maxAngular < 1e-13 && v.maxAngular < 1e-13);                 // symplectiques : exact à l'arrondi
    CHECK(r.maxAngular > 1e-12 && r.maxAngular < 1e-6);                  // RK4 : petit mais non nul (invariant quadratique)
    // Énergie : Euler s'effondre (>10 %), les schémas d'ordre supérieur la gardent, dans l'ordre Euler symp. < Verlet < RK4.
    CHECK(e.maxEnergy > 0.1);
    CHECK(s.maxEnergy < 0.1 && s.maxEnergy > v.maxEnergy);
    CHECK(v.maxEnergy < 5e-3 && v.maxEnergy > r.maxEnergy);
    CHECK(r.maxEnergy < 1e-5);
}

// Chaos : un écart de 1e-9 sur un amas est amplifié exponentiellement ; sur le huit il ne croît que linéairement (stable).
void testNBodyChaos() {
    auto growth = [](const NBodyProblem& base) {
        NBodyProblem twin = base;
        twin.position[0].x += 1e-9;
        RK45 a, b;
        a.relTol = b.relTol = 1e-13;
        a.absTol = b.absTol = 1e-15;
        State ya = base.initialState(), yb = twin.initialState();
        advance(a, base.rhs(), 0.0, ya, 2.0);
        advance(b, twin.rhs(), 0.0, yb, 2.0);
        const double early = base.distance(ya, yb);
        advance(a, base.rhs(), 2.0, ya, 14.0);
        advance(b, twin.rhs(), 2.0, yb, 14.0);
        return base.distance(ya, yb) / early;  // amplification entre t = 2 et t = 16
    };
    CHECK(growth(NBodyProblem::figureEight()) < 30.0);                          // mesuré : 7,3 (croissance linéaire)
    CHECK(growth(NBodyProblem::randomCluster(6, 42, 1.0, 0.05)) > 1e3);         // mesuré : 6,6e4
}

// M7 : sphère de Plummer (conditions initiales de l'amas du module GPU). Théorie pour N grand, G = M = a = 1, sans troncature :
// E = -3 pi / 64 = -0,1473, 2T = -U, rayon de demi-masse 1,3048. Mesuré sur 7 graines à N = 4000 : E de -0,147 à -0,159 (moyenne
// -0,152 : retirer les corps au-delà de 10 a, soit 1,5 % de la masse, creuse un peu le puits), 2T/|U| de 0,96 à 1,00, rayon de
// demi-masse de 1,264 à 1,330.
void testPlummer() {
    const auto halfMassRadius = [](const NBodyProblem& p) {
        std::vector<double> r;
        for (const Vec3& x : p.position) r.push_back(x.norm());
        std::sort(r.begin(), r.end());
        return r[r.size() / 2];
    };
    const auto virial = [](const NBodyProblem& p, const State& y) { return 2.0 * p.kineticEnergy(y) / -p.potentialEnergy(y); };

    const NBodyProblem p = NBodyProblem::plummer(4000, 42, 1.0, 0.0);
    const State y = p.initialState();
    CHECK_NEAR(p.totalMass(), 1.0, 1e-12);
    CHECK(p.centerOfMass(y).norm() < 1e-13);   // centre de masse au repos à l'origine
    CHECK(p.momentum(y).norm() < 1e-13);
    CHECK(p.energy(y) < -0.14 && p.energy(y) > -0.17);       // mesuré : -0,1505
    CHECK(virial(p, y) > 0.93 && virial(p, y) < 1.05);       // équilibre du viriel : mesuré 0,995
    CHECK(halfMassRadius(p) > 1.24 && halfMassRadius(p) < 1.37);   // mesuré : 1,307

    // Déterminisme, et changement d'échelle exact : positions x a, vitesses / sqrt(a), donc E / a.
    const NBodyProblem same = NBodyProblem::plummer(4000, 42, 1.0, 0.0), other = NBodyProblem::plummer(4000, 43, 1.0, 0.0);
    CHECK(same.position[0].x == p.position[0].x && same.velocity[3999].z == p.velocity[3999].z);
    CHECK(other.position[0].x != p.position[0].x);
    const NBodyProblem wide = NBodyProblem::plummer(4000, 42, 2.0, 0.0);
    CHECK_NEAR(wide.energy(wide.initialState()) / p.energy(y), 0.5, 1e-9);
    CHECK_NEAR(halfMassRadius(wide) / halfMassRadius(p), 2.0, 1e-12);

    // L'équilibre tient en dynamique : N = 500, 500 pas de 0,01 (t = 5, environ une traversée), rayon de demi-masse à 10 % et
    // énergie conservée (mesuré : +3,4 % et 4e-7).
    const NBodyProblem small = NBodyProblem::plummer(500, 3, 1.0, 0.05);
    State ys = small.initialState();
    VelocityVerlet verlet;
    const OdeFunction f = small.rhs();
    const double e0 = small.energy(ys), r0 = halfMassRadius(small);
    double t = 0.0;
    for (int i = 0; i < 500; ++i) t += verlet.step(f, t, ys, 0.01);
    NBodyProblem moved = small;
    for (int i = 0; i < 500; ++i) moved.position[i] = {ys[3 * i], ys[3 * i + 1], ys[3 * i + 2]};
    CHECK(std::abs(halfMassRadius(moved) / r0 - 1.0) < 0.1);
    CHECK(std::abs(small.energy(ys) / e0 - 1.0) < 1e-4);

    // Collision : deux sphères de masse 1/2 qui se rapprochent, système lié, centre de masse au repos à l'origine.
    const NBodyProblem c = NBodyProblem::plummerCollision(2000, 7, 6.0, 0.4, 0.6, 0.0);
    const State yc = c.initialState();
    CHECK_NEAR(c.totalMass(), 1.0, 1e-12);
    CHECK(c.momentum(yc).norm() < 1e-13 && c.centerOfMass(yc).norm() < 1e-13);
    CHECK(c.energy(yc) < -0.13 && c.energy(yc) > -0.16);     // lié ; mesuré : -0,148 (attendu -0,144 sans troncature)
    double leftMass = 0.0, leftX = 0.0, leftVx = 0.0;
    for (int i = 0; i < 1000; ++i) {
        leftMass += c.mass[i];
        leftX += c.position[i].x / 1000.0;
        leftVx += c.velocity[i].x / 1000.0;
    }
    CHECK_NEAR(leftMass, 0.5, 1e-12);
    CHECK_NEAR(leftX, -3.0, 1e-12);    // centre de la sphère de gauche en -séparation/2
    CHECK_NEAR(leftVx, 0.2, 1e-12);    // elle avance à la moitié de la vitesse relative
}

// --- M5 : événements, frottement sec, chocs ----------------------------------------------

void testEventDetection() {
    // Chute libre y'' = -g depuis y = 1 m : l'événement « y = 0 » a lieu à t* = sqrt(2/g), bien avant la fin d'un pas de 1 s.
    const double g = constants::g0;
    const OdeFunction f = [g](double, const State& y, State& d) { d[0] = y[1]; d[1] = -g; };
    const EventFunction ground = [](double, const State& y) { return y[0]; };
    const double tStar = std::sqrt(2.0 / g);

    RK4 rk4;
    RK45 rk45;
    for (Solver* solver : {static_cast<Solver*>(&rk4), static_cast<Solver*>(&rk45)}) {
        State y{1.0, 0.0};
        const EventStep r = advanceToEvent(*solver, f, ground, 0.0, y, 1.0);
        CHECK(r.event);
        CHECK_NEAR(r.elapsed, tStar, 1e-12);
        CHECK_NEAR(y[0], 0.0, 1e-12);
        CHECK_NEAR(y[1], -g * tStar, 1e-11);
    }

    // Pas sans événement : on avance de la durée demandée, comme advance().
    State y{1.0, 0.0};
    EventStep r = advanceToEvent(rk4, f, ground, 0.0, y, 0.2);
    CHECK(!r.event);
    CHECK_NEAR(r.elapsed, 0.2, 1e-15);
    CHECK_NEAR(y[0], 1.0 - 0.5 * g * 0.04, 1e-13);

    // Départ exactement sur la surface (g = 0) en s'en éloignant : ce n'est pas un événement.
    State up{0.0, 1.0};
    r = advanceToEvent(rk4, f, ground, 0.0, up, 0.1);
    CHECK(!r.event);
    CHECK_NEAR(r.elapsed, 0.1, 1e-15);

    // Événement dépendant du temps : g(t, y) = t - 0.37 s'annule à t = 0.37, quel que soit l'état.
    State any{0.0, 0.0};
    r = advanceToEvent(rk4, f, [](double t, const State&) { return t - 0.37; }, 0.0, any, 1.0);
    CHECK(r.event);
    CHECK_NEAR(r.elapsed, 0.37, 1e-12);

    // Deux événements dans un même pas (g change deux fois de signe) : non détectés, c'est une limite documentée du pas.
    // On vérifie au moins que le comportement reste défini (pas d'événement, avance complète).
    State twice{0.0, 0.0};
    r = advanceToEvent(rk4, f, [](double t, const State&) { return (t - 0.2) * (t - 0.6); }, 0.0, twice, 1.0);
    CHECK(r.elapsed > 0.0);
}

InclineProblem makeIncline(double angle, double muStatic, double muKinetic, double s0, double v0, double drag = 0.0) {
    InclineProblem p;
    p.angle = angle;
    p.muStatic = muStatic;
    p.muKinetic = muKinetic;
    p.s0 = s0;
    p.v0 = v0;
    p.drag = drag;
    return p;
}

// Budget de pas : un solveur adaptatif qui s'effondre sur une dynamique discontinue rend la main au lieu de boucler.
void testAdvanceBudget() {
    // Bloc qui doit rester collé (tan(theta) < mu_d), modèle naïf sgn(v) : après l'arrêt, la solution « glisse » le long de la
    // surface v = 0, où l'erreur locale reste d'ordre h : RK45 réduit son pas sans fin (jusqu'à 1e-14) pour une tolérance 1e-8.
    const InclineProblem p = makeIncline(0.2, 0.6, 0.5, 0.0, -3.0, 0.7);
    const OdeFunction f = p.rhs(InclineModel::Naive);
    RK45 rk;
    rk.relTol = 1e-8;
    rk.absTol = 1e-10;
    State y = p.initialState();
    double t = 0.0;
    bool exhausted = false;
    const double dt = 0.0125;
    for (int i = 0; i < 400 && !exhausted; ++i) {
        rk.resetStats();
        const double elapsed = advance(rk, f, t, y, dt, 300);
        CHECK(rk.acceptedSteps() <= 300);           // le budget est respecté
        CHECK(elapsed > 0.0);                       // et du temps a bien été avancé
        t += elapsed;
        if (elapsed < dt * (1.0 - 1e-9)) exhausted = true;
    }
    CHECK(exhausted);                               // RK45 s'est bloqué
    CHECK(t > p.firstStopTime());                   // après l'arrêt du bloc, pas avant

    // Budget suffisant : même résultat que sans budget. Sans budget (défaut), comportement inchangé.
    const OdeFunction decay = [](double, const State& s, State& d) { d[0] = -s[0]; };
    RK4 a, b;
    State ya{1.0}, yb{1.0};
    const double ea = advance(a, decay, 0.0, ya, 1.0);
    const double eb = advance(b, decay, 0.0, yb, 1.0, 5);
    CHECK_NEAR(ea, 1.0, 1e-15);
    CHECK_NEAR(eb, 1.0, 1e-15);
    CHECK_NEAR(ya[0], yb[0], 0.0);
}

void testInclineExact() {
    const double g = constants::g0;

    // A. tan(0,3) = 0,31 <= mu_s : le bloc monte, s'arrête et RESTE COLLÉ.
    {
        const InclineProblem p = makeIncline(0.3, 0.5, 0.4, 1.0, 6.0);
        CHECK(p.holds());
        const double A = -g * (std::sin(0.3) + 0.4 * std::cos(0.3));
        CHECK_NEAR(p.stageAcceleration(+1), A, 1e-14);
        const double t1 = -p.v0 / A, s1 = p.s0 - p.v0 * p.v0 / (2.0 * A);
        CHECK_NEAR(p.firstStopTime(), t1, 1e-13);

        InclineState st = p.exact(0.0);
        CHECK_NEAR(st.s, 1.0, 1e-15);
        CHECK_NEAR(st.v, 6.0, 1e-15);
        CHECK(!st.stuck);
        st = p.exact(0.5 * t1);
        CHECK_NEAR(st.s, p.s0 + p.v0 * 0.5 * t1 + 0.5 * A * 0.25 * t1 * t1, 1e-13);
        CHECK_NEAR(st.v, p.v0 + A * 0.5 * t1, 1e-13);
        CHECK_NEAR(st.path, st.s - p.s0, 1e-13);
        st = p.exact(t1 + 3.0);
        CHECK_NEAR(st.s, s1, 1e-12);
        CHECK_NEAR(st.v, 0.0, 0.0);
        CHECK_NEAR(st.path, s1 - p.s0, 1e-12);
        CHECK(st.stuck);
    }

    // B. tan(0,6) = 0,68 > mu_s : le bloc monte, s'arrête, puis REDESCEND (glissement uniformément accéléré vers le bas).
    {
        const InclineProblem p = makeIncline(0.6, 0.5, 0.4, 0.0, 5.0);
        CHECK(!p.holds());
        const double Aup = -g * (std::sin(0.6) + 0.4 * std::cos(0.6)), Adown = -g * (std::sin(0.6) - 0.4 * std::cos(0.6));
        CHECK(Adown < 0.0);
        const double t1 = -p.v0 / Aup, s1 = -p.v0 * p.v0 / (2.0 * Aup);
        const double dt = 0.8;
        const InclineState st = p.exact(t1 + dt);
        CHECK_NEAR(st.s, s1 + 0.5 * Adown * dt * dt, 1e-12);
        CHECK_NEAR(st.v, Adown * dt, 1e-12);
        CHECK_NEAR(st.path, s1 + (s1 - st.s), 1e-12);  // montée puis descente
        CHECK(!st.stuck);
    }

    // C. Vers le bas, tan(0,2) < mu_d : le frottement FREINE la descente, le bloc s'arrête puis reste collé.
    {
        const InclineProblem p = makeIncline(0.2, 0.6, 0.5, 0.0, -3.0);
        const double A = -g * (std::sin(0.2) - 0.5 * std::cos(0.2));
        CHECK(A > 0.0);
        const double t1 = -p.v0 / A, s1 = p.v0 * p.v0 / (-2.0 * A);  // s1 < 0
        CHECK_NEAR(p.firstStopTime(), t1, 1e-13);
        const InclineState st = p.exact(10.0 * t1);
        CHECK_NEAR(st.s, s1, 1e-12);
        CHECK(st.stuck);
        CHECK_NEAR(st.path, -s1, 1e-12);
    }

    // D. Départ au repos : adhérence si tan(theta) <= mu_s, sinon glissement vers le bas.
    {
        const InclineProblem hold = makeIncline(0.3, 0.5, 0.4, 2.0, 0.0);
        const InclineState a = hold.exact(7.0);
        CHECK_NEAR(a.s, 2.0, 0.0);
        CHECK(a.stuck);
        CHECK_NEAR(hold.firstStopTime(), 0.0, 0.0);
        const InclineProblem slide = makeIncline(0.6, 0.5, 0.4, 2.0, 0.0);
        const double Adown = -g * (std::sin(0.6) - 0.4 * std::cos(0.6));
        const InclineState b = slide.exact(1.5);
        CHECK_NEAR(b.s, 2.0 + 0.5 * Adown * 2.25, 1e-12);
        CHECK(!b.stuck);
        CHECK(std::isinf(slide.firstStopTime()));
    }

    // Bilan d'énergie sans résistance : E(t) + mu_d g cos(theta) x chemin = E(0), pour toutes les phases.
    for (const InclineProblem& p : {makeIncline(0.3, 0.5, 0.4, 1.0, 6.0), makeIncline(0.6, 0.5, 0.4, 0.0, 5.0),
                                   makeIncline(0.2, 0.6, 0.5, 0.0, -3.0)}) {
        const double e0 = p.energy(p.s0, p.v0);
        for (double t : {0.1, 0.5, 0.9, 1.3, 2.0, 3.5, 8.0}) {
            const InclineState st = p.exact(t);
            CHECK_NEAR(p.energy(st.s, st.v) + p.muKinetic * g * std::cos(p.angle) * st.path, e0, 1e-11);
        }
    }

    // Avec résistance k : v' = A - k v dans chaque phase (différences finies), v continue à l'arrêt, arrêt en ln(1 - k v0 / A)/k.
    {
        const InclineProblem p = makeIncline(0.6, 0.5, 0.4, 0.0, 5.0, 0.7);
        const double Aup = p.stageAcceleration(+1), Adown = p.stageAcceleration(-1);
        const double t1 = std::log(1.0 - p.drag * p.v0 / Aup) / p.drag;
        CHECK_NEAR(p.firstStopTime(), t1, 1e-12);
        CHECK_NEAR(p.exact(t1).v, 0.0, 1e-12);
        const double h = 1e-4;
        for (double t : {0.3 * t1, 0.8 * t1, t1 + 0.4, t1 + 1.5}) {
            const InclineState c = p.exact(t), up = p.exact(t + h), dn = p.exact(t - h);
            CHECK_NEAR((up.s - dn.s) / (2.0 * h), c.v, 1e-6);                            // s' = v
            const double A = t < t1 ? Aup : Adown;
            CHECK_NEAR((up.v - dn.v) / (2.0 * h), A - p.drag * c.v, 1e-6);               // v' = A - k v
        }
        // k -> 0 : on retrouve le mouvement uniformément accéléré (exp évitée par expm1 : pas de perte de chiffres)
        const InclineProblem tiny = makeIncline(0.6, 0.5, 0.4, 0.0, 5.0, 1e-12), none = makeIncline(0.6, 0.5, 0.4, 0.0, 5.0, 0.0);
        for (double t : {0.4, 1.0, 2.5}) CHECK_NEAR(tiny.exact(t).s, none.exact(t).s, 1e-9);
    }
}

// « Événement + adhérence » : arrêt exact à v = 0 puis repos ou demi-tour. Retrouve les ordres des solveurs (ici avec résistance
// k = 0,7 pour que la solution ne soit pas un simple polynôme : sans résistance Verlet et RK4 sont exacts, rien à mesurer).
void testInclineEventDriven() {
    struct Regime { InclineProblem p; double tEnd; bool stuck; };
    const Regime regimes[] = {{makeIncline(0.3, 0.5, 0.4, 1.0, 6.0, 0.7), 5.0, true},    // A : monte, colle
                              {makeIncline(0.6, 0.5, 0.4, 0.0, 5.0, 0.7), 4.0, false},   // B : monte, redescend
                              {makeIncline(0.2, 0.6, 0.5, 0.0, -3.0, 0.7), 5.0, true}}; // C : descente freinée, colle
    for (const Regime& r : regimes) {
        const InclineState ex = r.p.exact(r.tEnd);
        auto run = [&](Solver& solver, double stopTolerance) {
            InclineRun sim(r.p);
            for (int i = 0; i < static_cast<int>(std::lround(r.tEnd / 0.05)); ++i) sim.advance(solver, 0.05);
            CHECK(sim.stops() == 1);                   // un seul arrêt dans ces trois scénarios
            CHECK_NEAR(sim.stopTime(), r.p.firstStopTime(), stopTolerance);  // instant d'arrêt trouvé par bissection
            CHECK(sim.state().stuck == r.stuck);
            CHECK_NEAR(sim.time(), r.tEnd, 1e-12);
            return sim.state();
        };
        RK4 rk4;
        VelocityVerlet verlet;
        RK45 rk45;
        const InclineState a = run(rk4, 1e-6), b = run(verlet, 1e-2), c = run(rk45, 1e-8);  // l'arrêt hérite de l'ordre du solveur
        CHECK_NEAR(a.s, ex.s, 1e-6);       // mesuré : 4e-8 à 9e-8
        CHECK_NEAR(a.v, ex.v, 1e-7);
        CHECK_NEAR(a.path, ex.path, 1e-6);
        CHECK_NEAR(b.s, ex.s, 1e-2);       // mesuré : 1,4e-3 (ordre 2)
        CHECK_NEAR(c.s, ex.s, 1e-8);       // mesuré : 1e-10 (tolérance par défaut)
        CHECK_NEAR(c.path, ex.path, 1e-8);
    }

    // Sans résistance chaque phase est un polynôme de degré 2 : Verlet et RK4 sont exacts à l'arrondi.
    {
        const InclineProblem p = makeIncline(0.6, 0.5, 0.4, 0.0, 5.0, 0.0);
        RK4 rk4;
        VelocityVerlet verlet;
        for (Solver* s : {static_cast<Solver*>(&rk4), static_cast<Solver*>(&verlet)}) {
            InclineRun sim(p);
            for (int i = 0; i < 80; ++i) sim.advance(*s, 0.05);
            CHECK_NEAR(sim.state().s, p.exact(4.0).s, 1e-12);
        }
    }

    // Ordres de convergence de l'erreur de position à t = 4 s (cas B) : ceux des solveurs, l'événement ayant supprimé la discontinuité.
    const InclineProblem b = regimes[1].p;
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    auto ratio = [&](Solver& s, int n) {
        return inclineError(b, InclineModel::EventDriven, s, n, 4.0) / inclineError(b, InclineModel::EventDriven, s, 2 * n, 4.0);
    };
    CHECK_NEAR(ratio(euler, 200), 2.0, 0.15);
    CHECK_NEAR(ratio(symplectic, 200), 2.0, 0.15);
    CHECK_NEAR(ratio(verlet, 100), 4.0, 0.4);
    CHECK_NEAR(ratio(rk4, 100), 16.0, 1.5);

    // Départ au repos : adhérence exacte (rien ne bouge, le temps passe) ou glissement exact.
    {
        const InclineProblem hold = makeIncline(0.3, 0.5, 0.4, 2.0, 0.0);
        InclineRun sim(hold);
        for (int i = 0; i < 10; ++i) sim.advance(rk4, 0.1);
        CHECK(sim.state().stuck);
        CHECK_NEAR(sim.state().s, 2.0, 0.0);
        CHECK_NEAR(sim.state().v, 0.0, 0.0);
        CHECK_NEAR(sim.time(), 1.0, 1e-12);
        CHECK(sim.stops() == 0);
        CHECK(std::isnan(sim.stopTime()));   // pas d'arrêt : rien à signaler

        const InclineProblem slide = makeIncline(0.6, 0.5, 0.4, 2.0, 0.0);
        InclineRun down(slide);
        for (int i = 0; i < 20; ++i) down.advance(rk4, 0.1);
        CHECK(!down.state().stuck);
        CHECK_NEAR(down.state().s, slide.exact(2.0).s, 1e-12);
    }
}

// Modèle naïf (sgn dans l'EDO) : la discontinuité fait tomber TOUS les schémas à l'ordre 1, le bloc ne s'arrête jamais exactement,
// et il ignore mu_s. C'est la raison d'être de l'événement.
void testInclineNaive() {
    const InclineProblem a = makeIncline(0.3, 0.5, 0.4, 1.0, 6.0, 0.7);  // monte puis doit rester collé
    InclineProblem smooth = a;                                            // même code, mais sans frottement sec : EDO lisse
    smooth.muStatic = smooth.muKinetic = 0.0;

    RK4 rk4;
    VelocityVerlet verlet;
    auto ratio = [](const InclineProblem& p, Solver& s, int n) {
        return inclineError(p, InclineModel::Naive, s, n, 5.0) / inclineError(p, InclineModel::Naive, s, 2 * n, 5.0);
    };
    // mesuré : RK4 2,05 et 1,94 (n = 400, 800) contre 16,1 sans frottement sec ; Verlet 2,0 contre 4,0
    CHECK_NEAR(ratio(a, rk4, 400), 2.0, 0.4);
    CHECK_NEAR(ratio(a, rk4, 800), 2.0, 0.4);
    CHECK_NEAR(ratio(smooth, rk4, 200), 16.0, 1.5);
    CHECK_NEAR(ratio(a, verlet, 400), 2.0, 0.3);
    CHECK_NEAR(ratio(smooth, verlet, 400), 4.0, 0.3);

    // À pas égal (n = 400), l'événement est des millions de fois plus précis (mesuré : 6,7e-2 contre 3,3e-10).
    const double naive = inclineError(a, InclineModel::Naive, rk4, 400, 5.0);
    const double evented = inclineError(a, InclineModel::EventDriven, rk4, 400, 5.0);
    CHECK(naive > 1e6 * evented);

    // mu_d < tan(theta) <= mu_s : le bloc au repos doit RESTER COLLÉ (adhérence). Le modèle naïf ne connaît que mu_d : il glisse.
    {
        const InclineProblem p = makeIncline(0.45, 0.5, 0.4, 0.0, 0.0);
        CHECK(p.holds());
        CHECK_NEAR(p.exact(5.0).s, 0.0, 0.0);
        State y = p.initialState();
        const OdeFunction f = p.rhs(InclineModel::Naive);
        double t = 0.0;
        for (int i = 0; i < 500; ++i) t += advance(rk4, f, t, y, 0.01);
        CHECK(y[0] < -5.0);                       // mesuré : -9,2 m en 5 s
        InclineRun sim(p);
        for (int i = 0; i < 500; ++i) sim.advance(rk4, 0.01);
        CHECK_NEAR(sim.state().s, 0.0, 0.0);      // l'événement + adhérence : exactement immobile
    }

    // Vitesse résiduelle : après son arrêt le bloc (qui devrait rester immobile) garde une vitesse proportionnelle à dt (mesuré :
    // 1,18 |A| dt) ; avec l'événement elle est exactement nulle.
    {
        const InclineProblem p = makeIncline(0.2, 0.6, 0.5, 0.0, -3.0, 0.0);
        const double tStop = p.firstStopTime(), A = std::abs(p.stageAcceleration(-1));
        auto residual = [&](double dt) {
            State y = p.initialState();
            const OdeFunction f = p.rhs(InclineModel::Naive);
            double t = 0.0, vmax = 0.0;
            const int n = static_cast<int>(std::lround(5.0 / dt));
            for (int i = 0; i < n; ++i) {
                t += advance(rk4, f, t, y, dt);
                if (t > tStop + 0.3) vmax = std::max(vmax, std::abs(y[1]));
            }
            return vmax;
        };
        const double v1 = residual(0.01), v2 = residual(0.005);
        CHECK(v1 > 0.5 * A * 0.01 && v1 < 2.0 * A * 0.01);
        CHECK_NEAR(v1 / v2, 2.0, 0.2);
        InclineRun sim(p);
        for (int i = 0; i < 500; ++i) sim.advance(rk4, 0.01);
        CHECK_NEAR(sim.state().v, 0.0, 0.0);
    }
}

// Régularisation sgn(v) -> tanh(v / eps) : continue, mais ce n'est pas le bon modèle. Le bloc rampe au lieu d'adhérer, et l'erreur
// vient de eps (modèle), pas du pas : diminuer dt n'y change rien (mesuré : même s(5 s) pour dt = 0,01 et 0,001).
void testInclineRegularized() {
    InclineProblem p = makeIncline(0.2, 0.6, 0.5, 0.0, -3.0, 0.0);  // freine puis doit coller
    const double exact = p.exact(5.0).s;
    auto run = [&](double eps, double dt) {
        p.regularization = eps;
        RK4 rk4;
        State y = p.initialState();
        const OdeFunction f = p.rhs(InclineModel::Regularized);
        double t = 0.0;
        for (int i = 0; i < static_cast<int>(std::lround(5.0 / dt)); ++i) t += advance(rk4, f, t, y, dt);
        return y;
    };
    CHECK_NEAR(run(0.05, 0.01)[0], run(0.05, 0.001)[0], 1e-4);       // indépendant du pas
    const double e02 = std::abs(run(0.2, 0.001)[0] - exact), e005 = std::abs(run(0.05, 0.001)[0] - exact),
                 e001 = std::abs(run(0.01, 0.001)[0] - exact);
    CHECK(e001 < e005 && e005 < e02);                                  // l'écart au vrai modèle suit eps (mesuré : 0,017 ; 0,085 ; 0,35)
    CHECK(e005 > 0.05);
    CHECK(run(0.05, 0.001)[1] != 0.0);                                 // la vitesse n'est jamais nulle : le bloc rampe

    // Là où le bloc doit rester collé (mu_d < tan <= mu_s), il s'éloigne en rampant (mesuré : -9,4 m en 5 s).
    InclineProblem q = makeIncline(0.45, 0.5, 0.4, 0.0, 0.0);
    q.regularization = 0.05;
    RK45 rk45;
    State y = q.initialState();
    const OdeFunction f = q.rhs(InclineModel::Regularized);
    double t = 0.0;
    for (int i = 0; i < 5; ++i) t += advance(rk45, f, t, y, 1.0);
    CHECK(y[0] < -5.0);
}

// --- M5b : chocs et rebonds --------------------------------------------------------------

void testCollide1D() {
    using collision::collide1D;
    // e = 1, masses égales : les vitesses s'échangent, rien n'est perdu.
    collision::Result1D r = collide1D(2.0, 3.0, 2.0, -1.0, 1.0);
    CHECK_NEAR(r.v1, -1.0, 1e-15);
    CHECK_NEAR(r.v2, 3.0, 1e-15);
    CHECK_NEAR(r.energyLoss, 0.0, 1e-14);
    // e = 0 : choc mou, vitesse commune = impulsion / masse totale.
    r = collide1D(1.0, 4.0, 3.0, -2.0, 0.0);
    CHECK_NEAR(r.v1, -0.5, 1e-15);
    CHECK_NEAR(r.v2, -0.5, 1e-15);
    // Mur (masse énorme) : v1' = -e v1, le mur ne bouge pas.
    r = collide1D(1.0, 5.0, 1e12, 0.0, 0.6);
    CHECK_NEAR(r.v1, -3.0, 1e-9);
    CHECK_NEAR(r.v2, 0.0, 1e-11);
    // Cas général : impulsion conservée, vitesse relative inversée et réduite par e, énergie perdue = 1/2 (1 - e^2) mu v_rel^2.
    for (double e : {0.0, 0.3, 0.8, 1.0}) {
        const double m1 = 1.7, v1 = 4.2, m2 = 0.6, v2 = -1.1;
        const collision::Result1D c = collide1D(m1, v1, m2, v2, e);
        CHECK_NEAR(m1 * c.v1 + m2 * c.v2, m1 * v1 + m2 * v2, 1e-13);
        CHECK_NEAR(c.v1 - c.v2, -e * (v1 - v2), 1e-13);
        const double mu = m1 * m2 / (m1 + m2), vRel = v1 - v2;
        CHECK_NEAR(c.energyLoss, 0.5 * (1.0 - e * e) * mu * vRel * vRel, 1e-13);
        CHECK_NEAR(0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2 - 0.5 * m1 * c.v1 * c.v1 - 0.5 * m2 * c.v2 * c.v2, c.energyLoss, 1e-12);
    }
    // Les corps s'éloignent déjà : pas de choc.
    r = collide1D(1.0, 1.0, 1.0, 2.0, 0.5);
    CHECK_NEAR(r.v1, 1.0, 0.0);
    CHECK_NEAR(r.v2, 2.0, 0.0);
    CHECK_NEAR(r.energyLoss, 0.0, 0.0);
}

void testCollideSpheres() {
    using collision::collideSpheres;
    const Vec3 n{0.6, 0.8, 0.0};
    for (double e : {0.0, 0.5, 0.9, 1.0}) {
        const double m1 = 1.3, m2 = 0.4;
        const Vec3 v1{2.0, -1.0, 0.5}, v2{-1.5, 0.7, 0.2};
        const collision::Result3D r = collideSpheres(m1, v1, m2, v2, n, e);
        const double mu = m1 * m2 / (m1 + m2), vn = dot(v1 - v2, n);
        CHECK(vn > 0.0);
        CHECK_NEAR((m1 * r.v1 + m2 * r.v2 - m1 * v1 - m2 * v2).norm(), 0.0, 1e-14);   // impulsion conservée
        CHECK_NEAR(dot(r.v1 - r.v2, n), -e * vn, 1e-14);                                // restitution de Newton
        CHECK_NEAR((r.v1 - dot(r.v1, n) * n - (v1 - dot(v1, n) * n)).norm(), 0.0, 1e-14);  // tangentielles inchangées
        CHECK_NEAR((r.v2 - dot(r.v2, n) * n - (v2 - dot(v2, n) * n)).norm(), 0.0, 1e-14);
        CHECK_NEAR(r.impulse, (1.0 + e) * mu * vn, 1e-14);
        CHECK_NEAR(r.energyLoss, 0.5 * (1.0 - e * e) * mu * vn * vn, 1e-14);
        const double before = 0.5 * m1 * v1.norm2() + 0.5 * m2 * v2.norm2();
        CHECK_NEAR(before - 0.5 * m1 * r.v1.norm2() - 0.5 * m2 * r.v2.norm2(), r.energyLoss, 1e-13);
    }
    // Deux billes égales, l'une au repos : à 90 degrés si e = 1 ; sinon v1'.v2' = (1 - e^2)/4 v_n^2 (angle aigu).
    for (double alpha : {0.2, 0.7, 1.2}) {
        const Vec3 normal{std::cos(alpha), std::sin(alpha), 0.0}, v1{3.0, 0.0, 0.0}, v2;
        const double vn = dot(v1, normal);
        for (double e : {1.0, 0.6}) {
            const collision::Result3D r = collideSpheres(1.0, v1, 1.0, v2, normal, e);
            CHECK_NEAR(dot(r.v1, r.v2), (1.0 - e * e) / 4.0 * vn * vn, 1e-14);
        }
    }
    // Elles s'éloignent (v_n <= 0) : rien ne change.
    const collision::Result3D apart = collideSpheres(1.0, {-1.0, 0.0, 0.0}, 2.0, {0.5, 0.0, 0.0}, {1.0, 0.0, 0.0}, 0.5);
    CHECK_NEAR(apart.v1.x, -1.0, 0.0);
    CHECK_NEAR(apart.v2.x, 0.5, 0.0);
    CHECK_NEAR(apart.energyLoss, 0.0, 0.0);
}

void testTwoBallExact() {
    // Choc frontal de deux billes égales (r = 0,5) : contact à t = 1 s (distance 4 - 1 = 3 m à 3 m/s), elles échangent leurs vitesses.
    {
        TwoBallProblem p;
        p.p1 = {-4.0, 0.0, 0.0};
        CHECK_NEAR(p.gap(p.initialState()), 3.0, 1e-15);
        CHECK_NEAR(p.collisionTime(), 1.0, 1e-14);
        const State at = p.exact(1.0);
        CHECK_NEAR(at[0], -1.0, 1e-14);
        const State after = p.exact(2.0);
        CHECK_NEAR(after[0], -1.0, 1e-13);   // la bille 1 s'est arrêtée
        CHECK_NEAR(after[2], 3.0, 1e-13);    // la bille 2 est partie à 3 m/s
        CHECK_NEAR(after[4], 0.0, 1e-13);
        CHECK_NEAR(after[6], 3.0, 1e-13);
    }
    // Choc oblique, e = 1 : le contact a lieu quand les centres sont à r1 + r2, les billes repartent à 90 degrés ; P et E conservées.
    {
        TwoBallProblem p;
        p.p1 = {-4.0, 0.6, 0.0};
        const double tc = p.collisionTime();
        CHECK_NEAR(tc, (12.0 - 2.4) / 9.0, 1e-14);
        const State at = p.exact(tc);
        CHECK_NEAR(std::hypot(at[0] - at[2], at[1] - at[3]), 1.0, 1e-13);
        CHECK_NEAR(p.gap(at), 0.0, 1e-13);
        CHECK_NEAR(p.outgoingAngle(), constants::pi / 2.0, 1e-13);
        const double e0 = p.kineticEnergy(p.initialState());
        const Vec3 p0 = p.momentum(p.initialState());
        for (double t : {0.5, 1.5, 3.0}) {
            const State y = p.exact(t);
            CHECK_NEAR(p.kineticEnergy(y), e0, 1e-13);
            CHECK_NEAR((p.momentum(y) - p0).norm(), 0.0, 1e-13);
        }
        const State y = p.exact(3.0);
        CHECK_NEAR(y[4] * y[6] + y[5] * y[7], 0.0, 1e-13);   // vitesses finales orthogonales
    }
    // Restitution e = 0,6, masses inégales : l'énergie perdue est celle de la formule ; l'angle n'est plus de 90 degrés.
    {
        TwoBallProblem p;
        p.p1 = {-4.0, 0.6, 0.0};
        p.m2 = 3.0;
        p.restitution = 0.6;
        const State before = p.initialState(), after = p.exact(4.0);
        const Vec3 n = Vec3{p.exact(p.collisionTime())[2] - p.exact(p.collisionTime())[0], p.exact(p.collisionTime())[3] - p.exact(p.collisionTime())[1], 0.0}.normalized();
        const double mu = p.m1 * p.m2 / (p.m1 + p.m2), vn = dot(p.v1 - p.v2, n);
        CHECK_NEAR(p.kineticEnergy(before) - p.kineticEnergy(after), 0.5 * (1.0 - 0.36) * mu * vn * vn, 1e-13);
        CHECK_NEAR((p.momentum(after) - p.momentum(before)).norm(), 0.0, 1e-13);
    }
    // Pas de choc : trajectoires qui se manquent, ou billes qui s'éloignent.
    {
        TwoBallProblem miss;
        miss.p1 = {-4.0, 2.0, 0.0};
        CHECK(std::isinf(miss.collisionTime()));
        CHECK_NEAR(miss.exact(2.0)[0], 2.0, 1e-14);
        TwoBallProblem apart;
        apart.p1 = {-4.0, 0.0, 0.0};
        apart.v1 = {-3.0, 0.0, 0.0};
        CHECK(std::isinf(apart.collisionTime()));
    }
}

void testBounceExact() {
    const double g = constants::g0;
    BounceProblem p;
    p.restitution = 0.8;
    p.y0 = 2.0;
    p.vx0 = 1.5;
    p.restSpeed = 1e-7;
    const double e = p.restitution, h0 = 2.0, t0 = std::sqrt(2.0 * h0 / g);

    // Première chute libre.
    BounceState s = p.exact(0.5 * t0);
    CHECK_NEAR(s.y, h0 - 0.5 * g * 0.25 * t0 * t0, 1e-13);
    CHECK_NEAR(s.vy, -g * 0.5 * t0, 1e-13);
    CHECK(s.bounces == 0 && !s.resting);

    // Instants d'impact : t_n = t0 + 2 t0 e (1 - e^(n-1)) / (1 - e) ; sommet du n-ième rebond à t_n + e^n t0 : hauteur e^(2n) h0.
    const std::vector<double> impacts = p.impactTimes(6);
    CHECK(impacts.size() == 6);
    for (int n = 1; n <= 6 && n <= static_cast<int>(impacts.size()); ++n)
        CHECK_NEAR(impacts[n - 1], t0 + 2.0 * t0 * e * (1.0 - std::pow(e, n - 1)) / (1.0 - e), 1e-9);
    for (int n = 1; n <= 5; ++n) {
        s = p.exact(impacts[n - 1] + std::pow(e, n) * t0);
        CHECK_NEAR(s.y, std::pow(e, 2 * n) * h0, 1e-9);
        CHECK_NEAR(s.vy, 0.0, 1e-8);
        CHECK(s.bounces == n);
        CHECK_NEAR(s.vx, 1.5, 1e-15);   // sol lisse : la vitesse horizontale ne change pas
    }

    // Accumulation : la balle s'arrête de rebondir à t0 (1 + e) / (1 - e) (à l'effet du seuil près), puis glisse à vx constante.
    const double tRest = t0 * (1.0 + e) / (1.0 - e);
    CHECK_NEAR(p.restTime(), tRest, 1e-5);
    s = p.exact(tRest + 3.0);
    CHECK(s.resting);
    CHECK_NEAR(s.y, 0.0, 0.0);
    CHECK_NEAR(s.vy, 0.0, 0.0);
    CHECK_NEAR(s.x, 1.5 * (tRest + 3.0), 1e-5);
    for (double t : {0.3, 1.7, 4.0, tRest + 2.0}) CHECK_NEAR(p.exact(t).x, 1.5 * t, 1e-13);

    // Restitution e = 0 : un seul choc, puis repos. e = 1 : jamais de repos, la balle remonte à la hauteur de départ.
    BounceProblem soft = p;
    soft.restitution = 0.0;
    CHECK_NEAR(soft.restTime(), t0, 1e-12);
    BounceProblem elastic = p;
    elastic.restitution = 1.0;
    CHECK(std::isinf(elastic.restTime()));
    const std::vector<double> ei = elastic.impactTimes(4);
    CHECK_NEAR(elastic.exact(ei[2] + t0).y, h0, 1e-9);
    CHECK(!elastic.exact(50.0).resting);

    // Départ du sol vers le haut, ou au repos sur le sol.
    BounceProblem fromGround;
    fromGround.y0 = 0.0;
    fromGround.vy0 = 5.0;
    fromGround.restitution = 0.5;
    CHECK_NEAR(fromGround.impactTimes(1)[0], 2.0 * 5.0 / g, 1e-9);
    BounceProblem lying;
    lying.y0 = 0.0;
    CHECK(lying.exact(3.0).resting);
    CHECK_NEAR(lying.restTime(), 0.0, 0.0);

    // Avec résistance k = 0,7 : vx(t) = vx0 e^{-kt} (vols ET repos), contact y = 0 aux instants d'impact, EDO du vol vérifiée
    // par différences finies, vitesse de rebond = -e x vitesse d'impact.
    BounceProblem drag = p;
    drag.drag = 0.7;
    const std::vector<double> di = drag.impactTimes(4);
    for (double t : {0.2, 0.9, di[1] + 0.1, drag.restTime() + 1.0}) CHECK_NEAR(drag.exact(t).vx, 1.5 * std::exp(-0.7 * t), 1e-12);
    for (double ti : di) CHECK_NEAR(drag.exact(ti).y, 0.0, 1e-12);
    const double h = 1e-5;
    for (double t : {0.3, di[0] + 0.15, di[1] + 0.1}) {
        const BounceState c = drag.exact(t), up = drag.exact(t + h), dn = drag.exact(t - h);
        CHECK_NEAR((up.y - dn.y) / (2.0 * h), c.vy, 1e-6);
        CHECK_NEAR((up.vy - dn.vy) / (2.0 * h), -g - 0.7 * c.vy, 1e-5);
        CHECK_NEAR((up.vx - dn.vx) / (2.0 * h), -0.7 * c.vx, 1e-5);
    }
    for (int n = 0; n < 3; ++n) {
        const double vBefore = drag.exact(di[n] - 1e-9).vy, vAfter = drag.exact(di[n] + 1e-9).vy;
        CHECK_NEAR(vAfter, -e * vBefore, 1e-7);
    }
}

double fitSlope(const std::vector<double>& x, const std::vector<double>& y) {  // pente log-log par moindres carrés
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = static_cast<double>(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double lx = std::log10(x[i]), ly = std::log10(y[i]);
        sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// Balle rebondissante, modèle « événement » : impact localisé exactement, repos décidé au même instant que la solution exacte.
void testBounceEventDriven() {
    const double g = constants::g0;
    BounceProblem p;   // e = 0,8, chute de 2 m, restSpeed 1e-4
    p.vx0 = 1.5;

    // Sans résistance chaque vol est un polynôme : Verlet et RK4 sont exacts, même à travers l'accumulation de Zénon et le repos.
    RK4 rk4;
    VelocityVerlet verlet;
    for (Solver* s : {static_cast<Solver*>(&rk4), static_cast<Solver*>(&verlet)}) {
        BounceRun run(p, ContactModel::EventDriven);
        double t = 0.0;
        for (int i = 0; i < 800; ++i) {
            run.advance(*s, 0.01);
            t += 0.01;
            if (i == 109 || i == 289 || i == 509 || i == 799) {   // t = 1,1 ; 2,9 ; 5,1 ; 8
                const BounceState ex = p.exact(t);
                CHECK_NEAR(run.state()[0], ex.x, 1e-11);
                CHECK_NEAR(run.state()[1], ex.y, 1e-11);
                CHECK_NEAR(run.state()[2], ex.vx, 1e-11);
                CHECK_NEAR(run.state()[3], ex.vy, 1e-11);
            }
        }
        CHECK(run.resting());
        CHECK_NEAR(run.restTime(), p.restTime(), 1e-8);
        CHECK(run.bounces() == p.exact(8.0).bounces);   // 50 rebonds avant le repos
        CHECK_NEAR(run.state()[1], 0.0, 0.0);
        CHECK_NEAR(run.state()[3], 0.0, 0.0);
    }

    // Avec résistance k = 0,7 : RK4 à dt = 0,02 reste à 1e-7 de la solution exacte pendant 6 s (rebonds, accumulation et repos compris).
    BounceProblem drag = p;
    drag.drag = 0.7;
    {
        BounceRun run(drag, ContactModel::EventDriven);
        double worst = 0.0;
        for (int i = 0; i < 300; ++i) {
            run.advance(rk4, 0.02);
            const BounceState ex = drag.exact(0.02 * (i + 1));
            worst = std::max(worst, std::hypot(run.state()[0] - ex.x, run.state()[1] - ex.y));
        }
        CHECK(worst < 1e-7);   // mesuré : 3e-8 à dt = 0,05
        CHECK(run.resting());
    }

    // Ordres des solveurs (erreur de position à t = 2,37 s, après 3 rebonds, dans un vol) : ceux des schémas, la discontinuité ayant disparu.
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    auto ratio = [&](Solver& s, int n) {
        return bounceError(drag, ContactModel::EventDriven, s, n, 2.37) / bounceError(drag, ContactModel::EventDriven, s, 2 * n, 2.37);
    };
    CHECK_NEAR(ratio(euler, 1600), 2.0, 0.2);
    CHECK_NEAR(ratio(symplectic, 1600), 2.0, 0.2);
    CHECK_NEAR(ratio(verlet, 200), 4.0, 0.4);
    CHECK_NEAR(ratio(rk4, 100), 16.0, 1.5);

    // Régressions : la balle ne doit jamais traverser le sol, même quand les rebonds deviennent plus courts que le pas de calcul
    // (cas d'un pas de 0,05 s alors que les derniers vols durent quelques millisecondes).
    for (double dt : {0.05, 0.02, 0.007}) {
        BounceRun run(drag, ContactModel::EventDriven);
        double minY = 0.0;
        for (int i = 0; i < static_cast<int>(std::lround(8.0 / dt)); ++i) {
            run.advance(rk4, dt);
            minY = std::min(minY, run.state()[1]);
        }
        CHECK(minY >= -1e-12);
        CHECK(run.resting());
        CHECK(run.bounces() == drag.exact(8.0).bounces);                 // 46 rebonds, comme la solution exacte
        CHECK_NEAR(run.restTime(), drag.restTime(), 1e-5);               // mesuré : 2e-6 (dt = 0,05), 8e-8 (0,02), 2e-9 (0,007) : ordre 4
    }

    // Lancée du sol vers le haut, puis posée au sol dès le départ, puis e = 0.
    BounceProblem launched;
    launched.y0 = 0.0;
    launched.vy0 = 5.0;
    launched.restitution = 0.5;
    {
        BounceRun run(launched, ContactModel::EventDriven);
        for (int i = 0; i < 120; ++i) run.advance(rk4, 0.01);   // t = 1,2 s : le premier impact (2 v0 / g = 1,02 s) a eu lieu
        CHECK_NEAR(run.state()[1], launched.exact(1.2).y, 1e-11);
        CHECK_NEAR(run.state()[3], launched.exact(1.2).vy, 1e-11);
        CHECK(run.bounces() == 1 && launched.exact(1.2).bounces == 1);
    }
    BounceProblem lying;
    lying.y0 = 0.0;
    BounceRun still(lying, ContactModel::EventDriven);
    still.advance(rk4, 0.5);
    CHECK(still.resting());
    CHECK_NEAR(still.state()[1], 0.0, 0.0);
    BounceProblem soft;
    soft.restitution = 0.0;
    BounceRun mud(soft, ContactModel::EventDriven);
    for (int i = 0; i < 100; ++i) mud.advance(rk4, 0.01);
    CHECK(mud.resting());
    CHECK_NEAR(mud.restTime(), std::sqrt(4.0 / g), 1e-9);
    CHECK(mud.bounces() == 1);
    (void)g;
}

// Modèle naïf (pas fixe, rebond appliqué après le pas) : le contact est vu trop tard, tous les schémas tombent à l'ordre 1,
// et la balle ne s'arrête jamais. Mesures : voir PASSATION.
void testBounceNaive() {
    BounceProblem drag;   // e = 0,8, chute de 2 m
    drag.vx0 = 1.5;
    drag.drag = 0.7;
    RK4 rk4;

    // Pente log-log de l'erreur de position à t = 2,37 s en fonction du pas : ~ 1 pour le naïf (RK4 compris), ~ 4 pour l'événement.
    const std::vector<int> counts = {50, 100, 200, 400, 800, 1600};
    std::vector<double> dts, naive, evented;
    for (int n : counts) {
        dts.push_back(2.37 / n);
        naive.push_back(bounceError(drag, ContactModel::Naive, rk4, n, 2.37));
        evented.push_back(bounceError(drag, ContactModel::EventDriven, rk4, n, 2.37));
    }
    const double slopeNaive = fitSlope(dts, naive), slopeEvent = fitSlope(dts, evented);
    CHECK(slopeNaive > 0.6 && slopeNaive < 1.4);     // mesuré : 0,94
    CHECK(slopeEvent > 3.6 && slopeEvent < 4.4);     // mesuré : 4,0
    // À pas égal (400) l'événement est plus de 100000 fois plus précis (mesuré : 1,5e-2 contre 7e-11).
    CHECK(naive[3] > 1e5 * evented[3]);

    // Le premier impact est vu au pas suivant (jamais avant), donc en retard de moins d'un pas.
    {
        BounceProblem p;
        p.vx0 = 0.0;
        const double tImpact = std::sqrt(4.0 / constants::g0), dt = 0.001;
        BounceRun run(p, ContactModel::Naive);
        double detected = -1.0;
        for (int i = 0; i < 1000 && detected < 0.0; ++i) {
            run.advance(rk4, dt);
            if (run.bounces() == 1) detected = run.time();
        }
        CHECK(detected >= tImpact && detected < tImpact + dt);
    }

    // La balle naïve ne s'arrête jamais : bien après l'instant de repos exact elle rebondit encore, avec une vitesse de l'ordre de g dt,
    // d'autant plus petite que le pas est petit ; la balle événement, elle, est posée à l'instant exact avec vy = 0.
    BounceProblem lie;
    lie.vx0 = 0.0;   // e = 0,8, restSpeed 1e-4
    const double tRest = lie.restTime();
    double previous = 1e9;
    for (double dt : {0.02, 0.01, 0.005}) {
        BounceRun naiveRun(lie, ContactModel::Naive), eventRun(lie, ContactModel::EventDriven);
        double vMax = 0.0;
        for (int i = 0; i < static_cast<int>(std::lround(8.0 / dt)); ++i) {
            naiveRun.advance(rk4, dt);
            eventRun.advance(rk4, dt);
            if (naiveRun.time() > tRest + 0.5) vMax = std::max(vMax, std::abs(naiveRun.state()[3]));
        }
        CHECK(vMax > constants::g0 * dt);   // mesuré : 2,3 à 4 g dt
        CHECK(vMax < previous);
        previous = vMax;
        CHECK(!naiveRun.resting());
        CHECK(eventRun.resting());
        CHECK_NEAR(eventRun.restTime(), tRest, 1e-8);
    }
}

// Deux billes : hors choc le mouvement est rectiligne, donc TOUS les solveurs sont exacts ; l'erreur ne vient que de l'instant du contact.
void testTwoBallContact() {
    TwoBallProblem p;
    p.p1 = {-4.0, 0.6, 0.0};
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    std::vector<Solver*> solvers = {&euler, &symplectic, &verlet, &rk4};

    for (Solver* s : solvers) {
        for (int n : {30, 120, 480}) CHECK(twoBallError(p, ContactModel::EventDriven, *s, n, 3.0) < 1e-12);  // mesuré : <= 8e-14
    }
    {
        TwoBallRun run(p, ContactModel::EventDriven);
        for (int i = 0; i < 300; ++i) run.advance(rk4, 0.01);
        CHECK(run.collisions() == 1);
        CHECK_NEAR(run.collisionTime(), p.collisionTime(), 1e-12);
        CHECK_NEAR(run.time(), 3.0, 1e-12);
    }

    // Naïf : même erreur pour les quatre solveurs (ordre 1 pour tous, RK4 compris) : proportionnelle au pas (x 4 pour un pas 4 fois plus grand).
    for (Solver* s : solvers) {
        const double coarse = twoBallError(p, ContactModel::Naive, *s, 60, 3.0), fine = twoBallError(p, ContactModel::Naive, *s, 240, 3.0);
        CHECK(coarse / fine > 3.5 && coarse / fine < 4.8);          // mesuré : 4,2
        CHECK_NEAR(coarse, twoBallError(p, ContactModel::Naive, euler, 60, 3.0), 1e-12);
    }
    CHECK(twoBallError(p, ContactModel::Naive, rk4, 480, 3.0) > 1e-3);   // mesuré : 2,5e-2

    // Le contact naïf est détecté au pas suivant, en retard de moins d'un pas ; l'impulsion totale est tout de même conservée.
    {
        TwoBallRun run(p, ContactModel::Naive);
        for (int i = 0; i < 300; ++i) run.advance(rk4, 0.01);
        CHECK(run.collisions() == 1);
        CHECK(run.collisionTime() >= p.collisionTime() && run.collisionTime() < p.collisionTime() + 0.01 + 1e-12);
        const State y = run.state();
        CHECK_NEAR((p.momentum(y) - p.momentum(p.initialState())).norm(), 0.0, 1e-13);
    }
}

// ---- M5c : berceau de Newton, contact de Hertz ----

// Deux billes identiques de masse 1, raideur 1e4, vitesse d'approche 1 (unités normalisées) : une bille lancée sur une bille au repos.
CradleProblem makeCradle(int balls, int launched, double stiffness = 1e4, double speed = 1.0) {
    CradleProblem p;
    p.balls = balls;
    p.launched = launched;
    p.stiffness = stiffness;
    p.speed = speed;
    return p;
}

void testHertzReference() {
    // Constante du temps de contact (4/5) Gamma(2/5) Gamma(1/2) / Gamma(9/10) : 2,9433 (valeur classique 2,943).
    CHECK_NEAR(hertz::contactConstant(), 2.9433, 1e-4);
    CHECK_NEAR(hertz::force(1e4, 0.02), 1e4 * std::pow(0.02, 1.5), 1e-12);
    CHECK_NEAR(hertz::force(1e4, 0.0), 0.0, 0.0);
    CHECK_NEAR(hertz::force(1e4, -0.01), 0.0, 0.0);               // pas de contact : pas de force (et jamais de force attractive)
    CHECK_NEAR(hertz::potentialEnergy(1e4, 0.02), 0.4 * 1e4 * std::pow(0.02, 2.5), 1e-12);
    CHECK_NEAR(hertz::potentialEnergy(1e4, -0.01), 0.0, 0.0);

    // delta_max : l'énergie cinétique relative 1/2 mu v^2 devient (2/5) k delta^(5/2) (bilan d'énergie à l'arrêt relatif).
    const double mu = 0.5, v = 1.0, k = 1e4;
    const double dMax = hertz::maxCompression(mu, v, k);
    CHECK_NEAR(hertz::potentialEnergy(k, dMax), 0.5 * mu * v * v, 1e-14);
    // Lois d'échelle : delta_max ~ v^(4/5), T ~ delta_max / v ~ v^(-1/5), delta_max ~ k^(-2/5).
    CHECK_NEAR(hertz::maxCompression(mu, 2.0 * v, k) / dMax, std::pow(2.0, 0.8), 1e-12);
    CHECK_NEAR(hertz::contactDuration(mu, 2.0 * v, k) / hertz::contactDuration(mu, v, k), std::pow(2.0, -0.2), 1e-12);
    CHECK_NEAR(hertz::maxCompression(mu, v, 2.0 * k) / dMax, std::pow(2.0, -0.4), 1e-12);
    CHECK_NEAR(hertz::contactDuration(mu, v, k), hertz::contactConstant() * dMax / v, 1e-15);
}

// Choc de deux billes intégré par RK45 serré, comparé aux références exactes (l'événement et l'EDO ne partagent aucun code
// avec la fonction bêta de la formule).
void testHertzTwoBalls() {
    CradleProblem p = makeCradle(2, 1);
    const double mu = 0.5 * p.mass;
    const double dMax = hertz::maxCompression(mu, p.speed, p.stiffness);
    const double T = hertz::contactDuration(mu, p.speed, p.stiffness);

    RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-14;
    const OdeFunction f = p.rhs();
    State y = p.initialState();
    const double E0 = p.energy(y);
    CHECK_NEAR(E0, 0.5 * p.mass * p.speed * p.speed, 1e-15);
    CHECK_NEAR(p.compression(y, 0), 0.0, 1e-15);                // les billes se touchent sans compression au départ

    // Compression maximale : vitesse relative nulle (v0 = v1), à la mi-durée du contact.
    const EventFunction relativeSpeed = [](double, const State& s) { return s[2] - s[3]; };
    EventStep r = advanceToEvent(rk, f, relativeSpeed, 0.0, y, 2.0 * T);
    CHECK(r.event);
    CHECK_NEAR(p.compression(y, 0), dMax, 1e-8 * dMax);
    CHECK_NEAR(r.elapsed, 0.5 * T, 1e-8 * T);
    CHECK_NEAR(p.energy(y), E0, 1e-10);
    CHECK_NEAR(p.potentialEnergy(y), 0.5 * mu * p.speed * p.speed, 1e-9);   // toute l'énergie relative est stockée dans le contact

    // Fin du contact : la compression revient à 0 à t = T ; sortie élastique = échange des vitesses.
    const EventFunction compression = [&p](double, const State& s) { return p.compression(s, 0); };
    const double tHalf = r.elapsed;
    r = advanceToEvent(rk, f, compression, tHalf, y, 2.0 * T);
    CHECK(r.event);
    CHECK_NEAR(tHalf + r.elapsed, T, 1e-8 * T);
    CHECK_NEAR(y[2], 0.0, 1e-9);
    CHECK_NEAR(y[3], p.speed, 1e-9);
    CHECK_NEAR(p.energy(y), E0, 1e-10);
    CHECK_NEAR(p.momentum(y), p.mass * p.speed, 1e-12);
}

// Chaîne de 5 billes (une lancée) intégrée par RK4 : invariants, fin de collision, vitesses ordonnées.
void testCradleInvariants() {
    CradleProblem p = makeCradle(5, 1);
    RK4 rk4;
    const double dt = p.suggestedStep();
    CHECK(dt > 0.0 && dt < hertz::contactDuration(0.5 * p.mass, p.speed, p.stiffness));

    const State y0 = p.initialState();
    CHECK(!p.collisionOver(y0));                                // la bille lancée touche déjà la première bille au repos
    CHECK_NEAR(p.momentum(y0), p.mass * p.speed, 1e-15);
    CHECK_NEAR(p.kineticEnergy(y0), 0.5 * p.mass * p.speed * p.speed, 1e-15);
    CHECK_NEAR(p.potentialEnergy(y0), 0.0, 1e-15);

    const CradleOutcome out = p.run(rk4, dt, 5.0);
    CHECK(out.finished);
    CHECK(out.time > 0.0 && out.time < 1.0);                    // l'onde traverse 4 contacts en bien moins d'une unité de temps
    CHECK(out.maxCompression > 0.0);
    CHECK((int)out.velocities.size() == 5);

    double momentum = 0.0, kinetic = 0.0;
    for (double v : out.velocities) { momentum += p.mass * v; kinetic += 0.5 * p.mass * v * v; }
    CHECK_NEAR(momentum, p.mass * p.speed, 1e-12);              // impulsion conservée (forces égales et opposées)
    CHECK_NEAR(kinetic, 0.5 * p.mass * p.speed * p.speed, 1e-5);   // énergie : plus aucun contact, tout est cinétique
    for (int i = 0; i + 1 < 5; ++i) CHECK(out.velocities[i] <= out.velocities[i + 1] + 1e-12);   // plus aucun choc ne peut avoir lieu
    CHECK(out.velocities[4] > 0.5);                             // la dernière bille part
}

// Issue de la collision de N billes (RK45 serré : l'erreur d'intégration est très inférieure aux écarts testés).
CradleOutcome runCradleTight(int balls, int launched, double stiffness = 1e4, double speed = 1.0) {
    const CradleProblem p = makeCradle(balls, launched, stiffness, speed);
    RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-14;
    return p.run(rk, p.suggestedStep(), 50.0);
}

// Une bille de vitesse v frappe deux billes au repos : impulsion et énergie ne fixent PAS le résultat. Ce n'est pas un point
// mais une courbe : v1 = -x v, et v2, v3 racines de t^2 - (1 + x) t + x (1 + x) = 0, avec 0 <= x <= 1/3.
void testThreeBallFamily() {
    for (double x : {0.0, 0.05, 0.071, 0.1, 0.25, 1.0 / 3.0}) {
        const std::array<double, 3> v = cradle::threeBallFamily(x, 2.0);
        CHECK_NEAR(v[0] + v[1] + v[2], 2.0, 1e-14);                                  // impulsion
        CHECK_NEAR(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], 4.0, 1e-13);              // énergie
        CHECK_NEAR(v[0], -2.0 * x, 1e-15);
        CHECK(v[0] <= v[1] + 1e-15 && v[1] <= v[2] + 1e-15);                          // plus aucun choc possible
    }
    // Les deux extrémités : « une entre, une sort » (x = 0), et les deux dernières billes ensemble (x = 1/3).
    const std::array<double, 3> a = cradle::threeBallFamily(0.0, 1.0), b = cradle::threeBallFamily(1.0 / 3.0, 1.0);
    CHECK_NEAR(a[0], 0.0, 1e-15);
    CHECK_NEAR(a[1], 0.0, 1e-15);
    CHECK_NEAR(a[2], 1.0, 1e-15);
    CHECK_NEAR(b[0], -1.0 / 3.0, 1e-15);
    CHECK_NEAR(b[1], 2.0 / 3.0, 1e-12);       // racine double (discriminant nul)
    CHECK_NEAR(b[2], 2.0 / 3.0, 1e-12);
    // Exemple à la main : x = 0,1 donne v2, v3 = (1,1 -/+ sqrt(0,77)) / 2.
    const std::array<double, 3> c = cradle::threeBallFamily(0.1, 1.0);
    CHECK_NEAR(c[1], 0.5 * (1.1 - std::sqrt(0.77)), 1e-14);
    CHECK_NEAR(c[2], 0.5 * (1.1 + std::sqrt(0.77)), 1e-14);
    // Hors de l'intervalle : ramené aux bornes (pas de racine d'un nombre négatif).
    const std::array<double, 3> d = cradle::threeBallFamily(0.5, 1.0), e = cradle::threeBallFamily(-0.2, 1.0);
    CHECK_NEAR(d[0], -1.0 / 3.0, 1e-15);
    CHECK_NEAR(e[2], 1.0, 1e-15);
}

// La dynamique du contact tranche : le résultat de Hertz ne dépend ni de la raideur ni de la vitesse (une seule échelle de longueur),
// et il tombe sur la courbe des solutions, strictement entre ses deux extrémités.
void testHertzThreeBalls() {
    const CradleOutcome ref = runCradleTight(3, 1);
    CHECK(ref.finished);
    const double x = -ref.velocities[0];
    CHECK_NEAR(x, 0.070952, 2e-6);                          // mesuré : 0,070952 ; ni 0 (impulsions) ni 1/3
    CHECK_NEAR(ref.velocities[1], 0.076403, 2e-6);
    CHECK_NEAR(ref.velocities[2], 0.994549, 2e-6);
    const std::array<double, 3> onCurve = cradle::threeBallFamily(x, 1.0);
    CHECK_NEAR(ref.velocities[1], onCurve[1], 1e-9);        // le point mesuré est bien sur la famille (impulsion et énergie conservées)
    CHECK_NEAR(ref.velocities[2], onCurve[2], 1e-9);

    // Mêmes vitesses relatives quelles que soient la raideur et la vitesse (rapportées à la vitesse initiale).
    const CradleOutcome stiff = runCradleTight(3, 1, 1e6, 1.0), fast = runCradleTight(3, 1, 1e4, 3.0), slow = runCradleTight(3, 1, 2e3, 0.2);
    for (int i = 0; i < 3; ++i) {
        CHECK_NEAR(stiff.velocities[i], ref.velocities[i], 1e-7);
        CHECK_NEAR(fast.velocities[i] / 3.0, ref.velocities[i], 1e-7);
        CHECK_NEAR(slow.velocities[i] / 0.2, ref.velocities[i], 1e-7);
    }
    CHECK(stiff.time < 0.2 * ref.time);                     // un contact plus raide est plus court (T ~ k^(-2/5) : 100 fois plus raide, 6,3 fois plus court)
}

// Chaîne de N billes, n lancées : symétrie exacte (miroir + changement de repère galiléen) entre n lancées et N - n lancées,
// et résultats de référence (RK45 serré) comme ancrage de non-régression.
void testHertzChain() {
    for (int N : {3, 4, 5, 6}) {
        for (int n = 1; n < N; ++n) {
            const CradleOutcome a = runCradleTight(N, n), b = runCradleTight(N, N - n);
            CHECK(a.finished && b.finished);
            for (int i = 0; i < N; ++i) CHECK_NEAR(a.velocities[i], 1.0 - b.velocities[N - 1 - i], 1e-8);
        }
    }
    const CradleOutcome five = runCradleTight(5, 1);   // une bille sur quatre : la dernière part à 0,989, les autres restent petites
    const double expected1[5] = {-0.071085, -0.030274, -0.014464, 0.127045, 0.988777};
    for (int i = 0; i < 5; ++i) CHECK_NEAR(five.velocities[i], expected1[i], 2e-6);
    const CradleOutcome two = runCradleTight(5, 2);    // deux sur trois : pas tout à fait « deux entrent, deux sortent »
    const double expected2[5] = {-0.112615, -0.041960, 0.214486, 0.800367, 1.139722};
    for (int i = 0; i < 5; ++i) CHECK_NEAR(two.velocities[i], expected2[i], 2e-6);
    for (const CradleOutcome* o : {&five, &two}) {
        double momentum = 0.0, kinetic = 0.0;
        for (double v : o->velocities) { momentum += v; kinetic += 0.5 * v * v; }
        CHECK_NEAR(momentum, o == &five ? 1.0 : 2.0, 1e-10);
        CHECK_NEAR(kinetic, o == &five ? 0.5 : 1.0, 1e-9);
    }
    // Le résultat de RK4 au pas conseillé est le même à 1e-5 près (contact lisse à l'échelle du pas).
    {
        const CradleProblem p = makeCradle(5, 1);
        RK4 rk4;
        const CradleOutcome o = p.run(rk4, p.suggestedStep(), 50.0);
        for (int i = 0; i < 5; ++i) CHECK_NEAR(o.velocities[i], five.velocities[i], 1e-5);
    }
}

// Vitesses initiales d'une chaîne de N billes dont les n premières sont lancées à la vitesse 1.
std::vector<double> launchVelocities(int N, int n) {
    std::vector<double> v(N, 0.0);
    for (int i = 0; i < n; ++i) v[i] = 1.0;
    return v;
}

// Impulsions séquentielles : un choc binaire à la fois (collideSpheres de M5b), jusqu'à ce qu'aucune paire voisine ne se rapproche.
void testSequentialImpulses() {
    using cradle::ResolveOrder;
    // e = 1 et masses égales : chaque choc ÉCHANGE les vitesses, donc le résultat est le tri des vitesses : n entrent, n sortent,
    // quel que soit l'ordre de résolution.
    for (ResolveOrder order : {ResolveOrder::LeftToRight, ResolveOrder::RightToLeft}) {
        for (int N = 3; N <= 7; ++N) {
            for (int n = 1; n < N; ++n) {
                const cradle::ImpulseResult r = cradle::sequentialImpulses(launchVelocities(N, n), 1.0, order);
                CHECK(r.converged);
                CHECK(r.collisions >= 1);
                CHECK_NEAR(r.energyLoss, 0.0, 1e-15);
                for (int i = 0; i < N; ++i) CHECK_NEAR(r.velocities[i], i < N - n ? 0.0 : 1.0, 1e-15);
            }
        }
    }
    // Chaque bille ne change de vitesse qu'à cause de ses voisines et chaque choc conserve l'impulsion, même avec perte.
    // e = 0 : les billes finissent ensemble (impulsion n v partagée entre N billes) ; tolérance de la convergence 1e-12.
    for (ResolveOrder order : {ResolveOrder::LeftToRight, ResolveOrder::RightToLeft}) {
        const cradle::ImpulseResult r = cradle::sequentialImpulses(launchVelocities(5, 2), 0.0, order);
        CHECK(r.converged);
        for (double v : r.velocities) CHECK_NEAR(v, 2.0 / 5.0, 1e-9);
        CHECK_NEAR(r.energyLoss, 0.5 * 2.0 - 0.5 * 5.0 * 0.4 * 0.4, 1e-9);   // 1 - 0,4 = 0,6 : perte parfaitement inélastique
    }
    // e = 0,9 : l'impulsion est conservée, l'énergie diminue, les vitesses finissent rangées (plus aucun choc possible).
    for (ResolveOrder order : {ResolveOrder::LeftToRight, ResolveOrder::RightToLeft}) {
        for (int N = 3; N <= 7; ++N) {
            const cradle::ImpulseResult r = cradle::sequentialImpulses(launchVelocities(N, 1), 0.9, order);
            CHECK(r.converged);
            double sum = 0.0, kinetic = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += r.velocities[i];
                kinetic += 0.5 * r.velocities[i] * r.velocities[i];
                if (i + 1 < N) CHECK(r.velocities[i] <= r.velocities[i + 1] + 1e-11);
            }
            CHECK_NEAR(sum, 1.0, 1e-12);
            CHECK(r.energyLoss > 0.0 && r.energyLoss < 0.5);
            CHECK_NEAR(kinetic, 0.5 - r.energyLoss, 1e-12);   // pertes cumulées choc par choc = différence des énergies cinétiques
        }
    }
    // Entrée dégénérée : aucune bille qui s'approche d'une autre, rien ne bouge.
    const cradle::ImpulseResult none = cradle::sequentialImpulses({0.0, 1.0, 2.0}, 0.5);
    CHECK(none.collisions == 0 && none.sweeps == 1 && none.converged);
    CHECK_NEAR(none.velocities[1], 1.0, 0.0);
}

// Où l'ordre de résolution compte et où il ne compte pas (mesuré, et contraire à l'intuition).
void testSequentialOrder() {
    using cradle::ResolveOrder;
    // Berceau : seule la paire de tête se rapproche à chaque instant, la causalité impose la suite des chocs (0,1) puis (1,2)...
    // Les deux ordres de passe appliquent donc les mêmes chocs dans le même ordre : résultat identique, même avec perte (e < 1).
    int cases = 0;
    for (int N = 3; N <= 9; ++N) {
        for (int n = 1; n < N; ++n) {
            for (double e : {0.99, 0.9, 0.7, 0.5, 0.2, 0.05, 0.0}) {
                const cradle::ImpulseResult a = cradle::sequentialImpulses(launchVelocities(N, n), e, ResolveOrder::LeftToRight);
                const cradle::ImpulseResult b = cradle::sequentialImpulses(launchVelocities(N, n), e, ResolveOrder::RightToLeft);
                CHECK(a.converged && b.converged);
                for (int i = 0; i < N; ++i) CHECK_NEAR(a.velocities[i], b.velocities[i], 1e-15);
                ++cases;
            }
        }
    }
    CHECK(cases == 245);
    // Bille prise entre deux voisines qui s'approchent toutes deux : deux chocs simultanés, l'ordre devient un choix du modèle.
    // Avec e = 1 (échange de vitesses = tri) il n'y paraît pas ; avec e = 0,5 l'écart est de l'ordre de 0,1 (mesuré : 0,07 sur (1; 0,5; 0; 0)).
    const std::vector<double> squeezed = {1.0, 0.5, 0.0, 0.0};
    {
        const cradle::ImpulseResult a = cradle::sequentialImpulses(squeezed, 1.0, ResolveOrder::LeftToRight);
        const cradle::ImpulseResult b = cradle::sequentialImpulses(squeezed, 1.0, ResolveOrder::RightToLeft);
        for (int i = 0; i < 4; ++i) CHECK_NEAR(a.velocities[i], b.velocities[i], 1e-15);
    }
    {
        const cradle::ImpulseResult a = cradle::sequentialImpulses(squeezed, 0.5, ResolveOrder::LeftToRight);
        const cradle::ImpulseResult b = cradle::sequentialImpulses(squeezed, 0.5, ResolveOrder::RightToLeft);
        double gap = 0.0;
        for (int i = 0; i < 4; ++i) gap = std::max(gap, std::abs(a.velocities[i] - b.velocities[i]));
        CHECK(gap > 0.05 && gap < 0.1);
        double pa = 0.0, pb = 0.0;
        for (int i = 0; i < 4; ++i) { pa += a.velocities[i]; pb += b.velocities[i]; }
        CHECK_NEAR(pa, 1.5, 1e-12);                    // l'impulsion reste conservée dans les deux cas
        CHECK_NEAR(pb, 1.5, 1e-12);
    }
}

// Écart initial entre le groupe lancé et le reste de la chaîne : le contact commence plus tard, sans changer le résultat.
void testCradleGap() {
    CradleProblem p = makeCradle(3, 1);
    p.gap = 0.05;
    const State y0 = p.initialState();
    CHECK_NEAR(y0[0], -0.05, 1e-15);                      // la bille lancée est reculée de `gap`
    CHECK_NEAR(y0[1], 1.0, 1e-15);
    CHECK_NEAR(p.compression(y0, 0), -0.05, 1e-15);      // séparées : compression négative
    CHECK_NEAR(p.compression(y0, 1), 0.0, 1e-15);
    CHECK(!p.collisionOver(y0));                         // la bille approche : le choc n'a pas encore eu lieu

    const CradleOutcome flush = runCradleTight(3, 1);
    RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-14;
    const CradleOutcome spaced = p.run(rk, p.suggestedStep(), 50.0);
    CHECK(spaced.finished);
    for (int i = 0; i < 3; ++i) CHECK_NEAR(spaced.velocities[i], flush.velocities[i], 1e-8);
    CHECK_NEAR(spaced.time - flush.time, p.gap / p.speed, 2.0 * p.suggestedStep());   // le choc a lieu après 0,05 s de vol libre
}

// Référence RK45 serré et erreur des solveurs à pas fixe sur le contact de Hertz.
void testCradleConvergence() {
    CradleProblem p = makeCradle(3, 1);
    p.gap = 0.0123;                                     // le contact démarre à un instant qui n'est pas sur la grille des pas
    const double tEnd = 0.2;                            // bien après la fin de la collision (0,094 + 0,0123 s)

    // La référence conserve l'impulsion et l'énergie et reproduit le résultat de run() après la collision.
    const State ref = p.reference(tEnd);
    CHECK_NEAR(p.momentum(ref), p.mass * p.speed, 1e-12);
    CHECK_NEAR(p.energy(ref), 0.5 * p.mass * p.speed * p.speed, 1e-10);
    CHECK(p.collisionOver(ref));
    const CradleOutcome out = runCradleTight(3, 1);
    for (int i = 0; i < 3; ++i) CHECK_NEAR(ref[3 + i], out.velocities[i], 1e-8);

    // L'erreur diminue quand le pas diminue, et un schéma d'ordre supérieur fait mieux à pas égal.
    ExplicitEuler euler;
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    RK4 rk4;
    for (Solver* s : {static_cast<Solver*>(&euler), static_cast<Solver*>(&symplectic), static_cast<Solver*>(&verlet), static_cast<Solver*>(&rk4)}) {
        const double coarse = cradleError(p, *s, 100, tEnd), fine = cradleError(p, *s, 400, tEnd);
        CHECK(fine < coarse);
        CHECK(coarse > 0.0);
    }
    CHECK(cradleError(p, rk4, 200, tEnd) < cradleError(p, verlet, 200, tEnd));
    CHECK(cradleError(p, verlet, 200, tEnd) < cradleError(p, euler, 200, tEnd));

    // Ordres effectifs mesurés (rapport des erreurs pour des pas divisés par 4, de 800 à 3200 pas ou de 400 à 1600).
    // Euler : ordre 1 (4,07 mesuré) ; Verlet : ordre 2 (15,4 mesuré, soit 3,9 par doublement).
    const double eulerRatio = cradleError(p, euler, 800, tEnd) / cradleError(p, euler, 3200, tEnd);
    CHECK(eulerRatio > 3.8 && eulerRatio < 4.4);
    const double verletRatio = cradleError(p, verlet, 400, tEnd) / cradleError(p, verlet, 1600, tEnd);
    CHECK(verletRatio > 14.0 && verletRatio < 17.0);
    // Euler symplectique donne EXACTEMENT les erreurs de Verlet (à l'arrondi : écart mesuré <= 5e-13) : les deux schémas sont
    // conjugués par un demi-pas de vitesse, qui est l'identité hors contact (état de départ et état final sans force).
    for (int steps : {50, 100, 400, 1600}) CHECK_NEAR(cradleError(p, symplectic, steps, tEnd), cradleError(p, verlet, steps, tEnd), 1e-11);
    // RK4 n'atteint PAS l'ordre 4 : la force k delta^(3/2) n'est pas lisse en delta = 0 (début et fin de chaque contact), ce qui
    // limite l'ordre à 2,5 environ ; l'ordre apparent dépend du décalage du début du contact sur la grille (mesuré 2,29 à 2,75).
    for (double gap : {0.0, 0.0123, 0.0231, 0.0377}) {
        CradleProblem q = p;
        q.gap = gap;
        const double order = std::log2(cradleError(q, rk4, 200, tEnd) / cradleError(q, rk4, 3200, tEnd)) / 4.0;
        CHECK(order > 2.0 && order < 3.0);
    }
    // RK45 à tolérance par défaut (1e-8) : erreur de l'ordre de 3e-7 (mesuré) ; 88 pas acceptés pour traverser tout [0 ; 0,2] d'un coup.
    RK45 rk45;
    CHECK(cradleError(p, rk45, 20, tEnd) < 1e-6);
    CHECK(rk45.acceptedSteps() > 0 && rk45.acceptedSteps() < 400);
}

// Amortissement de Hunt-Crossley F = k delta^(3/2) (1 + (3/2) alpha delta') : à la main, l'énergie relative perdue vaut
// 2 alpha v de l'énergie d'approche (premier ordre), donc la restitution est e = 1 - alpha v.
void testHertzDamping() {
    RK45 rk;
    rk.relTol = 1e-12;
    rk.absTol = 1e-14;
    for (double alpha : {0.01, 0.05, 0.1}) {
        CradleProblem p = makeCradle(2, 1);
        p.damping = alpha;
        const CradleOutcome out = p.run(rk, p.suggestedStep(), 5.0);
        CHECK(out.finished);
        const double restitution = out.velocities[1] - out.velocities[0];   // vitesse de séparation / vitesse d'approche (1)
        CHECK_NEAR(out.velocities[0] + out.velocities[1], 1.0, 1e-10);       // la force d'amortissement est aussi égale et opposée
        CHECK_NEAR(restitution, 1.0 - alpha, 1.5 * alpha * alpha);           // premier ordre : e = 1/(1 + alpha v) = 1 - alpha v + (alpha v)^2 (mesuré)
        CHECK(restitution < 1.0 && restitution > 0.8);
        CHECK(out.velocities[0] > 0.0);                                      // pas de recul : l'amortissement ne rend jamais la force attractive
    }
    // Sans amortissement (valeur par défaut) la restitution est exactement 1 : contact conservatif de Hertz.
    {
        const CradleProblem p = makeCradle(2, 1);
        CHECK_NEAR(p.damping, 0.0, 0.0);
        const CradleOutcome out = p.run(rk, p.suggestedStep(), 5.0);
        CHECK_NEAR(out.velocities[1] - out.velocities[0], 1.0, 1e-9);
    }
    // Chaîne amortie : impulsion conservée, énergie perdue (strictement), vitesses rangées. Le recul des premières billes diminue.
    {
        CradleProblem p = makeCradle(5, 1);
        p.damping = 0.05;
        const CradleOutcome out = p.run(rk, p.suggestedStep(), 50.0);
        CHECK(out.finished);
        double momentum = 0.0, kinetic = 0.0;
        for (double v : out.velocities) { momentum += v; kinetic += 0.5 * v * v; }
        CHECK_NEAR(momentum, 1.0, 1e-10);
        CHECK(kinetic < 0.5 - 0.01);
        for (int i = 0; i + 1 < 5; ++i) CHECK(out.velocities[i] <= out.velocities[i + 1] + 1e-12);
        CHECK(out.velocities[0] > runCradleTight(5, 1).velocities[0]);       // moins de recul qu'en contact conservatif (-0,071)
    }
}

// ---- M6 : corps rigide ----

// Moments d'inertie par intégration numérique directe (somme de Riemann au point milieu sur une grille) : contrôle indépendant des
// formules fermées. `inside(x, y, z)` décrit le solide ; la densité est uniforme, de masse totale m.
Vec3 gridInertia(double m, double halfX, double halfY, double halfZ, int n, bool (*inside)(double, double, double)) {
    double volume = 0.0, ixx = 0.0, iyy = 0.0, izz = 0.0;
    const double dx = 2.0 * halfX / n, dy = 2.0 * halfY / n, dz = 2.0 * halfZ / n, dv = dx * dy * dz;
    for (int i = 0; i < n; ++i) {
        const double x = -halfX + (i + 0.5) * dx;
        for (int j = 0; j < n; ++j) {
            const double y = -halfY + (j + 0.5) * dy;
            for (int k = 0; k < n; ++k) {
                const double z = -halfZ + (k + 0.5) * dz;
                if (!inside(x, y, z)) continue;
                volume += dv;
                ixx += (y * y + z * z) * dv;
                iyy += (x * x + z * z) * dv;
                izz += (x * x + y * y) * dv;
            }
        }
    }
    const double rho = m / volume;   // masse volumique déduite du volume de la grille (même erreur de bord sur les deux)
    return {rho * ixx, rho * iyy, rho * izz};
}

void testInertiaTensors() {
    // Valeurs à la main.
    const Vec3 s = inertia::sphere(2.0, 0.5);                   // 2/5 m r^2 = 0,2
    CHECK_NEAR(s.x, 0.2, 1e-15);
    CHECK_NEAR(s.y, 0.2, 1e-15);
    CHECK_NEAR(s.z, 0.2, 1e-15);
    const Vec3 b = inertia::box(3.0, 1.0, 2.0, 3.0);            // m/12 (b^2 + c^2), ... = 3,25 ; 2,5 ; 1,25
    CHECK_NEAR(b.x, 3.25, 1e-15);
    CHECK_NEAR(b.y, 2.5, 1e-15);
    CHECK_NEAR(b.z, 1.25, 1e-15);
    const Vec3 c = inertia::cylinder(2.0, 0.5, 2.0);            // axe z : 1/2 m r^2 = 0,25 ; m/12 (3 r^2 + h^2) = 0,79167
    CHECK_NEAR(c.z, 0.25, 1e-15);
    CHECK_NEAR(c.x, 2.0 / 12.0 * (0.75 + 4.0), 1e-15);
    CHECK_NEAR(c.y, c.x, 0.0);

    // Contrôle par intégration directe sur une grille (la sphère et le cylindre ont une erreur de bord de l'ordre de 1/n).
    const Vec3 gs = gridInertia(2.0, 0.5, 0.5, 0.5, 100, [](double x, double y, double z) { return x * x + y * y + z * z <= 0.25; });
    const Vec3 gb = gridInertia(3.0, 0.5, 1.0, 1.5, 60, [](double, double, double) { return true; });
    const Vec3 gc = gridInertia(2.0, 0.5, 0.5, 1.0, 100, [](double x, double y, double) { return x * x + y * y <= 0.25; });
    CHECK_NEAR(gs.x / s.x, 1.0, 5e-3);
    CHECK_NEAR(gs.z / s.z, 1.0, 5e-3);
    CHECK_NEAR(gb.x / b.x, 1.0, 1e-3);
    CHECK_NEAR(gb.y / b.y, 1.0, 1e-3);
    CHECK_NEAR(gb.z / b.z, 1.0, 1e-3);
    CHECK_NEAR(gc.z / c.z, 1.0, 5e-3);
    CHECK_NEAR(gc.x / c.x, 1.0, 5e-3);

    // Un solide réel vérifie l'inégalité triangulaire I_i + I_j >= I_k (égalité pour une figure plane).
    for (const Vec3& v : {s, b, c}) {
        CHECK(v.x + v.y >= v.z - 1e-12 && v.y + v.z >= v.x - 1e-12 && v.x + v.z >= v.y - 1e-12);
    }
}

// Équations d'Euler I w' = (I w) x w, orientation q' = 1/2 q (0, w).
void testEulerEquations() {
    FreeBodyProblem p;
    p.inertia = {1.0, 2.0, 3.0};
    // À la main : w = (1, 2, 3) -> w1' = (I2 - I3)/I1 w2 w3 = -6 ; w2' = (I3 - I1)/I2 w3 w1 = 3 ; w3' = (I1 - I2)/I3 w1 w2 = -2/3.
    const Vec3 d = p.omegaDot({1.0, 2.0, 3.0});
    CHECK_NEAR(d.x, -6.0, 1e-15);
    CHECK_NEAR(d.y, 3.0, 1e-15);
    CHECK_NEAR(d.z, -2.0 / 3.0, 1e-15);
    // Rotation autour d'un axe principal : état stationnaire (aucun couple en repère du corps).
    for (const Vec3& w : {Vec3{5.0, 0.0, 0.0}, Vec3{0.0, -2.0, 0.0}, Vec3{0.0, 0.0, 7.0}}) CHECK_NEAR(p.omegaDot(w).norm(), 0.0, 1e-15);
    // Solide sphérique (I1 = I2 = I3) : jamais de précession, w constant quel que soit w.
    p.inertia = {0.4, 0.4, 0.4};
    CHECK_NEAR(p.omegaDot({1.0, -2.0, 3.0}).norm(), 0.0, 1e-15);

    // Cinématique : q = identité, w = (0, 0, 2) -> q' = 1/2 (0, w) = (0, 0, 0, 1).
    p.inertia = {1.0, 2.0, 3.0};
    p.q0 = Quaternion{};
    p.omega0 = {0.0, 0.0, 2.0};
    State dy(7);
    p.rhs()(0.0, p.initialState(), dy);
    CHECK_NEAR(dy[0], 0.0, 1e-15);
    CHECK_NEAR(dy[1], 0.0, 1e-15);
    CHECK_NEAR(dy[2], 0.0, 1e-15);
    CHECK_NEAR(dy[3], 1.0, 1e-15);
    CHECK_NEAR(dy[6], 0.0, 1e-15);
    // Aller-retour de l'état : (q, w) <-> vecteur de 7 nombres.
    const RotationState rs = unpackRotation(p.initialState());
    CHECK_NEAR(rs.q.w, 1.0, 0.0);
    CHECK_NEAR(rs.omega.z, 2.0, 0.0);
    CHECK((packRotation(rs) == p.initialState()));
}

// Invariants : énergie, moment cinétique (vecteur, dans le repère FIXE), norme du quaternion le long d'une trajectoire de référence.
void testFreeBodyInvariants() {
    FreeBodyProblem p;
    p.inertia = {1.0, 2.0, 3.0};
    p.q0 = Quaternion::fromAxisAngle({1.0, 2.0, 3.0}, 0.8);
    p.omega0 = {0.7, 1.3, -0.4};
    const RotationState s0{p.q0, p.omega0};

    // E = 1/2 sum I_i w_i^2 et L_corps = I w, à la main.
    CHECK_NEAR(p.energy(s0), 0.5 * (1.0 * 0.49 + 2.0 * 1.69 + 3.0 * 0.16), 1e-15);
    const Vec3 lb = p.angularMomentumBody(s0);
    CHECK_NEAR(lb.x, 0.7, 1e-15);
    CHECK_NEAR(lb.y, 2.6, 1e-15);
    CHECK_NEAR(lb.z, -1.2, 1e-15);
    // Le moment cinétique dans le repère fixe est R L_corps : même norme.
    const Vec3 ls = p.angularMomentumSpace(s0);
    CHECK_NEAR(ls.norm(), lb.norm(), 1e-14);
    const Vec3 manual = p.q0.toMat3() * lb;
    CHECK_NEAR((ls - manual).norm(), 0.0, 1e-14);

    // Le long de la trajectoire exacte (RK45 serré) : tout est conservé, y compris le VECTEUR L dans le repère fixe.
    for (double t : {1.0, 5.0, 15.0, 40.0}) {
        const RotationState s = p.reference(t);
        CHECK_NEAR(p.energy(s), p.energy(s0), 1e-11);
        CHECK_NEAR((p.angularMomentumSpace(s) - ls).norm(), 0.0, 1e-11);
        CHECK_NEAR(s.q.norm(), 1.0, 1e-11);
        // Le corps a bien tourné (le test ne doit pas passer sur un état figé).
        CHECK((s.omega - p.omega0).norm() > 1e-3);
    }
}

// Distance entre deux orientations : norme de Frobenius de la différence des matrices de rotation (insensible au signe q <-> -q).
double orientationDistance(const Quaternion& a, const Quaternion& b) {
    const Mat3 ra = a.toMat3(), rb = b.toMat3();
    double sum = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) sum += (ra.m[i][j] - rb.m[i][j]) * (ra.m[i][j] - rb.m[i][j]);
    return std::sqrt(sum);
}

// Parcourt la trajectoire de référence (RK45 serré) en UNE seule intégration, en appelant visit(t, état) tous les dt : appeler reference(t)
// à chaque instant ré-intègre depuis 0 chaque fois (coût quadratique : 76 s mesurés pour ce seul test en Debug).
template <class Visit>
void walkReference(const OdeFunction& f, State y, double dt, int count, Visit&& visit) {
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    for (int i = 1; i <= count; ++i) {
        advance(rk, f, (i - 1) * dt, y, dt, 100000);
        visit(i * dt, unpackRotation(y));
    }
}

// Solide symétrique libre (I1 = I2) : w3 constante, (w1, w2) tourne à Omega = (I3 - I1) w3 / I1 dans le repère du corps, l'axe du corps
// décrit un cône autour du vecteur fixe L à la vitesse L / I1. Formule fermée de l'orientation : q(t) = q_L(L t / I1) q0 q_3(-Omega t).
void testSymmetricTopExact() {
    struct Case { Vec3 inertia; Vec3 omega; };
    const Case cases[] = {{{2.0, 2.0, 1.0}, {0.9, -0.4, 3.0}},     // allongé (I3 < I1) : le « cigare »
                          {{1.0, 1.0, 2.0}, {0.5, 0.8, 4.0}},      // aplati (I3 > I1) : le disque
                          {{1.5, 1.5, 1.5}, {1.0, -2.0, 0.7}}};    // sphérique : Omega = 0, rotation uniforme
    for (const Case& c : cases) {
        FreeBodyProblem p;
        p.inertia = c.inertia;
        p.q0 = Quaternion::fromAxisAngle({0.3, -1.0, 0.6}, 1.1);
        p.omega0 = c.omega;
        const RotationState s0{p.q0, p.omega0};
        const double L = p.angularMomentumBody(s0).norm(), Omega = (c.inertia.z - c.inertia.x) * c.omega.z / c.inertia.x;
        const Vec3 Ls = p.angularMomentumSpace(s0);

        // t = 0 : on retrouve les conditions initiales.
        const RotationState e0 = p.exactSymmetric(0.0);
        CHECK_NEAR(orientationDistance(e0.q, p.q0), 0.0, 1e-14);
        CHECK_NEAR((e0.omega - p.omega0).norm(), 0.0, 1e-15);

        for (double t : {0.7, 3.0, 10.0, 25.0}) {
            const RotationState e = p.exactSymmetric(t), r = p.reference(t);
            CHECK_NEAR((e.omega - r.omega).norm(), 0.0, 1e-10);                       // w(t) = formule : accord avec l'intégration serrée (mesuré 2e-12)
            CHECK_NEAR(orientationDistance(e.q, r.q), 0.0, 1e-10);                    // q(t) aussi : orientation entière (mesuré 4e-12)
            CHECK_NEAR(e.q.norm(), 1.0, 1e-14);                                       // la formule fermée reste exactement unitaire
            CHECK_NEAR(e.omega.z, c.omega.z, 1e-15);                                  // w3 constante
            CHECK_NEAR(std::hypot(e.omega.x, e.omega.y), std::hypot(c.omega.x, c.omega.y), 1e-14);   // module de (w1, w2) constant
            CHECK_NEAR((p.angularMomentumSpace(e) - Ls).norm(), 0.0, 1e-13);          // L fixe dans le repère fixe
            CHECK_NEAR(p.energy(e), p.energy(s0), 1e-13);
        }
        // Le cône : l'angle entre l'axe 3 du corps (dans le repère fixe) et L est constant, cos(theta) = I3 w3 / L.
        const Vec3 Lhat = Ls / L;
        for (double t : {0.0, 2.0, 7.5}) {
            const Vec3 axis3 = p.exactSymmetric(t).q.rotate({0.0, 0.0, 1.0});
            CHECK_NEAR(dot(axis3, Lhat), c.inertia.z * c.omega.z / L, 1e-13);
        }
        // Précession à la vitesse L / I1 : après une période 2 pi I1 / L l'axe du corps est revenu dans sa direction initiale.
        const double period = 2.0 * constants::pi * c.inertia.x / L;
        const Vec3 a0 = p.q0.rotate({0.0, 0.0, 1.0}), a1 = p.exactSymmetric(period).q.rotate({0.0, 0.0, 1.0});
        CHECK_NEAR((a1 - a0).norm(), 0.0, 1e-13);
        (void)Omega;
    }
    // Rotation pure autour de l'axe de symétrie : L est parallèle à l'axe 3, le corps tourne simplement de w3 t autour de lui.
    FreeBodyProblem spin;
    spin.inertia = {2.0, 2.0, 1.0};
    spin.q0 = Quaternion::fromAxisAngle({1.0, 1.0, 0.0}, 0.5);
    spin.omega0 = {0.0, 0.0, 3.0};
    const Quaternion expected = spin.q0 * Quaternion::fromAxisAngle({0.0, 0.0, 1.0}, 3.0 * 4.0);
    CHECK_NEAR(orientationDistance(spin.exactSymmetric(4.0).q, expected), 0.0, 1e-13);
}

// Fonctions elliptiques de Jacobi sn, cn, dn (paramètre m = k^2) : les trois sont l'UNIQUE solution du système sn' = cn dn, cn' = -sn dn,
// dn' = -m sn cn avec (sn, cn, dn)(0) = (0, 1, 1). On compare donc la série de Fourier à une intégration RK4 de ce système (aucun code commun).
void testJacobiFunctions() {
    for (double m : {1e-14, 0.01, 0.3, 0.5, 0.75, 0.95, 0.999}) {
        // Intégration de référence jusqu'à u = 5,3 (plus d'une période pour m petit : 4K ~ 6,3).
        State y{0.0, 1.0, 1.0};
        const OdeFunction f = [m](double, const State& s, State& d) { d[0] = s[1] * s[2]; d[1] = -s[0] * s[2]; d[2] = -m * s[0] * s[1]; };
        RK4 rk4;
        double u = 0.0;
        const double du = 1e-3;
        for (int i = 0; i < 5300; ++i) {
            if (i % 530 == 0 && i > 0) {
                double sn, cn, dn;
                jacobiSnCnDn(u, m, sn, cn, dn);
                CHECK_NEAR(sn, y[0], 1e-10);
                CHECK_NEAR(cn, y[1], 1e-10);
                CHECK_NEAR(dn, y[2], 1e-10);
                CHECK_NEAR(sn * sn + cn * cn, 1.0, 5e-14);         // identités algébriques (mesuré : 4e-15, et 1,4e-14 pour m = 1e-14 : préfacteur 1/k)
                CHECK_NEAR(dn * dn + m * sn * sn, 1.0, 1e-14);
            }
            u += rk4.step(f, u, y, du);
        }
    }
    // Valeurs particulières : sn(K) = 1, cn(K) = 0, dn(K) = sqrt(1 - m), sn(0) = 0, sn(-u) = -sn(u), période 4K.
    const double m = 0.64, k = 0.8, kk = ellipticK(k);
    double sn, cn, dn, sn2, cn2, dn2;
    jacobiSnCnDn(kk, m, sn, cn, dn);
    CHECK_NEAR(sn, 1.0, 1e-14);
    CHECK_NEAR(cn, 0.0, 1e-14);
    CHECK_NEAR(dn, std::sqrt(1.0 - m), 1e-14);
    jacobiSnCnDn(0.0, m, sn, cn, dn);
    CHECK_NEAR(sn, 0.0, 1e-15);
    CHECK_NEAR(cn, 1.0, 1e-15);
    CHECK_NEAR(dn, 1.0, 1e-15);
    jacobiSnCnDn(0.9, m, sn, cn, dn);
    jacobiSnCnDn(-0.9, m, sn2, cn2, dn2);
    CHECK_NEAR(sn2, -sn, 1e-15);
    CHECK_NEAR(cn2, cn, 1e-15);
    jacobiSnCnDn(0.9 + 4.0 * kk, m, sn2, cn2, dn2);
    CHECK_NEAR(sn2, sn, 1e-12);
    CHECK_NEAR(cn2, cn, 1e-12);
    // m = 0 : fonctions trigonométriques.
    jacobiSnCnDn(1.3, 0.0, sn, cn, dn);
    CHECK_NEAR(sn, std::sin(1.3), 1e-15);
    CHECK_NEAR(cn, std::cos(1.3), 1e-15);
    CHECK_NEAR(dn, 1.0, 1e-15);
}

// Solide asymétrique libre (I1 < I2 < I3), départ w = (a, 0, c) : w(t) en fonctions de Jacobi (voir RigidBody.hpp), comparé à l'intégration
// serrée des équations d'Euler. Deux régimes séparés par la séparatrice 2 E I2 = L^2 (rotation autour de l'axe 3 ou de l'axe 1).
void testAsymmetricExact() {
    struct Case { Vec3 omega; double m; bool nearAxis3; };
    const Case cases[] = {{{1.5, 0.0, 1.0}, 0.75, true},      // 2 E I2 = 10,5 < L^2 = 11,25 : w3 en dn, w1 en cn
                          {{2.0, 0.0, 0.5}, 0.1875, false},   // 2 E I2 = 9,5 > L^2 = 6,25 : w1 en dn, w3 en cn
                          {{1.5, 0.0, 0.8}, 0.853, false}};   // près de la séparatrice : 8,34 contre 8,01
    for (const Case& c : cases) {
        FreeBodyProblem p;
        p.inertia = {1.0, 2.0, 3.0};
        p.q0 = Quaternion{};
        p.omega0 = c.omega;
        const RotationState s0{p.q0, p.omega0};
        CHECK_NEAR(p.asymmetricParameter(), c.m, 5e-4);        // paramètre m = k^2 (à la main : voir les commentaires)
        const double period = p.asymmetricPeriod();
        CHECK(period > 0.0);
        for (double t : {0.0, 0.4, 2.0, 7.0, 20.0, 40.0}) {
            const Vec3 w = p.exactAsymmetricOmega(t);
            const RotationState r = p.reference(t);
            CHECK_NEAR((w - r.omega).norm(), 0.0, 1e-11);                                               // mesuré : <= 6e-13
            CHECK_NEAR(p.energy({p.q0, w}), p.energy(s0), 1e-13);                                       // E et |L| sont des invariants de la formule
            CHECK_NEAR(p.angularMomentumBody({p.q0, w}).norm(), p.angularMomentumBody(s0).norm(), 1e-13);
        }
        // Période de w(t) : 4 K(k) / lambda. Contrôle par l'intégration serrée, qui ne connaît pas la formule.
        for (double t : {0.0, 1.3, 5.0}) {
            const Vec3 a = p.reference(t).omega, b = p.reference(t + period).omega;
            CHECK_NEAR((a - b).norm(), 0.0, 1e-11);                                                     // mesuré : <= 1e-13
        }
        // Pas de période plus courte : à une demi-période w1 ou w3 a changé de signe (cn) : la vitesse n'est pas revenue.
        CHECK((p.reference(0.5 * period).omega - p.omega0).norm() > 0.5);
    }
}

// Instabilité de l'axe intermédiaire (la « raquette de tennis ») : autour de l'axe 2, une petite perturbation (w1, w3) croît comme
// exp(lambda t), avec lambda = w2 sqrt((I3 - I2)(I2 - I1) / (I1 I3)). Les axes extrêmes sont stables (simple oscillation).
void testIntermediateAxisInstability() {
    FreeBodyProblem p;
    p.inertia = {1.0, 2.0, 3.0};
    p.q0 = Quaternion{};
    const double lambda = std::sqrt((3.0 - 2.0) * (2.0 - 1.0) / (1.0 * 3.0)) * 1.0;           // w2 = 1
    p.omega0 = {0.0, 0.0, 0.0};
    CHECK_NEAR(p.intermediateAxisGrowthRate(), 0.0, 1e-15);                                     // w2 = 0 : pas de rotation, pas de croissance
    p.omega0 = {1e-9, 1.0, 1e-9};
    CHECK_NEAR(p.intermediateAxisGrowthRate(), lambda, 1e-15);

    // Pente de ln|w1| entre t = 8 et 14 (amplitude de 1e-9 e^(0,58 t) < 3e-6 : toujours linéaire), loin du mode décroissant.
    const double w1a = std::abs(p.reference(8.0).omega.x), w1b = std::abs(p.reference(14.0).omega.x);
    const double measured = std::log(w1b / w1a) / 6.0;
    CHECK_NEAR(measured, lambda, 5e-4);                         // mesuré : 0,577411 pour 0,577350 prévu (écart relatif 1e-4)
    CHECK(w1b > 100.0 * 1e-9);                                  // la perturbation a grossi d'au moins deux ordres de grandeur

    // Axes extrêmes : une perturbation de 1e-3 reste de l'ordre de 1e-3 pendant toute la durée (stable).
    for (const Vec3& spin : {Vec3{1.0, 1e-3, 1e-3}, Vec3{1e-3, 1e-3, 1.0}}) {
        FreeBodyProblem q = p;
        q.omega0 = spin;
        double worst = 0.0;
        walkReference(q.rhs(), q.initialState(), 0.5, 200, [&](double, const RotationState& r) {
            const Vec3& w = r.omega;
            worst = std::max(worst, (spin.x > 0.5) ? std::hypot(w.y, w.z) : std::hypot(w.x, w.y));
        });
        CHECK(worst < 5e-3);
    }
}

// Intégrateurs d'orientation : cas où la réponse est connue à l'avance.
void testRotationIntegratorBasics() {
    // Distance entre deux états : insensible au signe de q (q et -q sont la même rotation), et la partie vitesse compte.
    const RotationState a{Quaternion::fromAxisAngle({1, 2, 3}, 0.7), {0.1, 0.2, 0.3}};
    CHECK_NEAR(rotationDistance(a, a), 0.0, 0.0);
    CHECK_NEAR(rotationDistance(a, {Quaternion{-a.q.w, -a.q.x, -a.q.y, -a.q.z}, a.omega}), 0.0, 1e-15);
    CHECK_NEAR(rotationDistance(a, {a.q, {0.1, 0.5, 0.3}}), 0.3, 1e-15);

    // Rotation pure autour d'un axe principal : w constante, q(t) = q0 q_axe(w t). Réponse exacte pour tout pas.
    const Vec3 inertia{1.0, 2.0, 3.0};
    const RotationState s0{Quaternion::fromAxisAngle({1, 2, 3}, 0.7), {0.0, 2.0, 0.0}};
    const Quaternion exact = s0.q * Quaternion::fromAxisAngle({0.0, 1.0, 0.0}, 2.0 * 1.0);   // t = 1
    LieHeunRotation lie;
    SplittingRotation splitting;
    for (int steps : {1, 4, 50}) {                       // même avec UN SEUL pas de 1 s : l'exponentielle est exacte pour w constante
        const RotationState l = integrateRotation(lie, inertia, {}, s0, 1.0, steps);
        const RotationState d = integrateRotation(splitting, inertia, {}, s0, 1.0, steps);
        CHECK_NEAR(orientationDistance(l.q, exact), 0.0, 1e-14);
        CHECK_NEAR(orientationDistance(d.q, exact), 0.0, 1e-14);
        CHECK_NEAR(l.q.norm(), 1.0, 1e-15);
        CHECK_NEAR(d.q.norm(), 1.0, 1e-15);
        CHECK_NEAR(l.omega.y, 2.0, 1e-15);
        CHECK_NEAR(d.omega.y, 2.0, 1e-15);
    }
    RK4Rotation rk4;
    const RotationState r = integrateRotation(rk4, inertia, {}, s0, 1.0, 100);
    CHECK_NEAR(orientationDistance(r.q, exact), 0.0, 1e-9);             // ordre 4 : (h w)^5 / 120 par pas, très petit pour h w = 0,02
    CHECK_NEAR(r.q.norm(), 1.0, 1e-15);                                 // avec renormalisation

    // Euler explicite : chaque pas multiplie la norme de q par sqrt(1 + (h w)^2 / 4) EXACTEMENT (l'incrément q (0, w) est orthogonal à q),
    // donc |q| = (1 + h^2 w^2 / 4)^(N/2) : la contrainte |q| = 1 n'est pas respectée toute seule.
    EulerRotation euler;
    const double h = 0.01;
    const RotationState e = integrateRotation(euler, inertia, {}, s0, 1.0, 100);
    CHECK_NEAR(e.q.norm(), std::pow(1.0 + h * h * 2.0 * 2.0 / 4.0, 50.0), 1e-12);
    CHECK(e.q.norm() > 1.004);
    CHECK_NEAR(e.omega.y, 2.0, 1e-15);

    // RK4 sans renormalisation : la dérive de la norme est d'ordre élevé mais non nulle (elle s'accumule d'un pas à l'autre).
    RK4Rotation raw(false);
    const RotationState u = integrateRotation(raw, inertia, {}, s0, 1.0, 100);
    CHECK(std::abs(u.q.norm() - 1.0) < 1e-10);
}

// Découpage symplectique (Dullweber-Leimkuhler-McLachlan) : chaque sous-pas est une rotation EXACTE autour d'un axe principal, donc
// |q| = 1 et le vecteur L dans le repère fixe sont conservés à l'arrondi près, quel que soit le pas ; l'énergie reste bornée.
void testSplittingConservation() {
    const Vec3 inertia{1.0, 2.0, 3.0};
    const RotationState s0{Quaternion::fromAxisAngle({0.3, -1.0, 0.6}, 1.1), {0.1, 2.0, 0.1}};   // près de l'axe instable : mouvement violent
    FreeBodyProblem p;
    p.inertia = inertia;
    const Vec3 l0 = p.angularMomentumSpace(s0);
    const double e0 = p.energy(s0);

    SplittingRotation splitting;
    RotationState s = s0;
    double worstL = 0.0, worstQ = 0.0, worstE = 0.0;
    const double h = 0.1;                                   // h |w| = 0,2 : pas grossier
    for (int i = 0; i < 5000; ++i) {                        // 500 s
        splitting.step(inertia, {}, s, h);
        worstL = std::max(worstL, (p.angularMomentumSpace(s) - l0).norm());
        worstQ = std::max(worstQ, std::abs(s.q.norm() - 1.0));
        worstE = std::max(worstE, std::abs(p.energy(s) - e0) / e0);
    }
    CHECK(worstL < 1e-12);
    CHECK(worstQ < 1e-13);
    CHECK(worstE < 3.5e-4);                                 // mesuré : 2,74e-4 (borne de l'erreur d'énergie en h^2, sans dérive séculaire)

    // L'erreur d'énergie est en h^2 et ne dérive pas : même maximum sur chaque tiers du calcul (mouvement périodique).
    auto energyBound = [&](double step, int steps) {
        RotationState t = s0;
        double worst = 0.0;
        for (int i = 0; i < steps; ++i) {
            splitting.step(inertia, {}, t, step);
            worst = std::max(worst, std::abs(p.energy(t) - e0) / e0);
        }
        return worst;
    };
    const double fine = energyBound(0.05, 10000), coarse = energyBound(0.1, 5000);
    CHECK(coarse / fine > 3.7 && coarse / fine < 4.3);       // mesuré : 3,98
    const double lastThird = [&] {
        RotationState t = s0;
        double worst = 0.0;
        for (int i = 0; i < 5000; ++i) {
            splitting.step(inertia, {}, t, 0.1);
            if (i >= 3334) worst = std::max(worst, std::abs(p.energy(t) - e0) / e0);
        }
        return worst;
    }();
    CHECK(lastThird > 0.98 * coarse);                        // le dernier tiers atteint déjà le maximum : aucune dérive
}

// Mesures de convergence et de dérive sur une trajectoire de solide asymétrique (corps proche de l'axe intermédiaire : mouvement violent).
struct DriftStats {
    double maxEnergy = 0.0, maxNorm = 0.0, maxMomentum = 0.0, finalEnergy = 0.0;
    int divergedAt = -1;   // pas où l'état explose (|w| > 1e6 ou non fini), -1 sinon
};

DriftStats rotationDrift(RotationIntegrator& integrator, double h, int steps) {
    FreeBodyProblem p;
    p.inertia = {1.0, 2.0, 3.0};
    const RotationState s0{Quaternion::fromAxisAngle({0.3, -1.0, 0.6}, 1.1), {0.1, 2.0, 0.1}};
    const Vec3 l0 = p.angularMomentumSpace(s0);
    const double e0 = p.energy(s0);
    RotationState s = s0;
    DriftStats stats;
    for (int i = 0; i < steps; ++i) {
        integrator.step(p.inertia, {}, s, h);
        if (!std::isfinite(s.q.w) || !std::isfinite(s.omega.x) || s.omega.norm() > 1e6) {
            stats.divergedAt = i;
            break;
        }
        stats.finalEnergy = (p.energy(s) - e0) / e0;
        stats.maxEnergy = std::max(stats.maxEnergy, std::abs(stats.finalEnergy));
        stats.maxNorm = std::max(stats.maxNorm, std::abs(s.q.norm() - 1.0));
        stats.maxMomentum = std::max(stats.maxMomentum, (p.angularMomentumSpace(s) - l0).norm() / l0.norm());
    }
    return stats;
}

// Ordres de convergence mesurés (erreur à t = 6,3 contre RK45 serré) : le rapport pour des pas divisés par 4 vaut 4^ordre.
void testRotationOrders() {
    FreeBodyProblem p;
    p.inertia = {1.0, 2.0, 3.0};
    p.q0 = Quaternion::fromAxisAngle({0.3, -1.0, 0.6}, 1.1);
    p.omega0 = {0.9, 0.5, 1.1};
    const RotationState s0{p.q0, p.omega0};
    const double tEnd = 6.3;
    const RotationState ref = p.reference(tEnd);
    auto error = [&](RotationIntegrator& in, int n) { return rotationDistance(integrateRotation(in, p.inertia, {}, s0, tEnd, n), ref); };

    EulerRotation euler;
    RK4Rotation rk4;
    LieHeunRotation lie;
    SplittingRotation splitting;
    // De 100 à 400 pas : Euler ordre 1 (4,29 mesuré), RK4 ordre 4 (263 mesuré, 256 théorique), Heun et découpage ordre 2 (15,8 et 15,9).
    const double eulerRatio = error(euler, 100) / error(euler, 400);
    CHECK(eulerRatio > 4.0 && eulerRatio < 4.6);
    const double rk4Ratio = error(rk4, 100) / error(rk4, 400);
    CHECK(rk4Ratio > 240.0 && rk4Ratio < 290.0);
    const double lieRatio = error(lie, 100) / error(lie, 400);
    CHECK(lieRatio > 15.0 && lieRatio < 16.6);
    const double splitRatio = error(splitting, 100) / error(splitting, 400);
    CHECK(splitRatio > 15.0 && splitRatio < 16.6);
    // À pas égal le découpage est environ 3 fois plus précis que Heun (constante d'erreur plus faible) ; RK4 les écrase tous deux aux pas fins.
    CHECK(error(splitting, 100) < 0.5 * error(lie, 100));
    CHECK(error(rk4, 400) < 1e-3 * error(splitting, 400));
    CHECK(error(rk4, 1600) < 1e-10);                          // mesuré : 3,5e-11
}

// Dérives sur 1000 s (20000 pas de 0,05) : ce que la renormalisation corrige et ce qu'elle ne corrige pas.
void testRotationDrift() {
    EulerRotation euler;
    RK4Rotation rk4, raw(false);
    LieHeunRotation lie;
    SplittingRotation splitting;

    const DriftStats d = rotationDrift(splitting, 0.05, 20000);
    CHECK(d.divergedAt < 0);
    CHECK(d.maxMomentum < 1e-12);                             // mesuré : 6e-14 : L conservé à l'arrondi près
    CHECK(d.maxNorm < 1e-14);                                 // mesuré : 3e-16
    CHECK(d.maxEnergy < 1e-4);                                // mesuré : 6,9e-5, borné
    CHECK(std::abs(d.finalEnergy) < 1e-4);

    const DriftStats r = rotationDrift(rk4, 0.05, 20000);
    CHECK(r.divergedAt < 0);
    CHECK(r.maxNorm < 1e-15);                                 // la renormalisation tient |q| = 1
    CHECK(r.finalEnergy < -1e-6 && r.finalEnergy > -1e-5);    // mesuré : -3,1e-6 : RK4 perd lentement de l'énergie
    CHECK(r.maxMomentum > 1e-7 && r.maxMomentum < 1e-5);      // mesuré : 1,8e-6 : L dérive (la renormalisation ne rétablit pas L)
    CHECK(d.maxMomentum < 1e-6 * r.maxMomentum);              // le découpage conserve L au moins un million de fois mieux

    const DriftStats u = rotationDrift(raw, 0.05, 20000);
    CHECK(u.maxNorm > 1e-7 && u.maxNorm < 1e-4);              // mesuré : 4e-6 : la norme dérive sans renormalisation
    CHECK(u.maxMomentum > r.maxMomentum);                     // et L aussi, d'un ordre de grandeur de plus (1,7e-5 contre 1,8e-6)

    const DriftStats l = rotationDrift(lie, 0.05, 20000);
    CHECK(l.divergedAt < 0);
    CHECK(l.maxNorm < 1e-14);                                 // le groupe de Lie tient |q| = 1 ...
    CHECK(l.maxEnergy > 1e-3 && l.maxEnergy < 0.1);           // ... mais PAS l'énergie : mesuré 1,4e-2 (dérive séculaire) ...
    CHECK(l.maxMomentum > 1e-3);                              // ... ni L : mesuré 7e-3. Garder |q| = 1 ne suffit pas.
    CHECK(rotationDrift(lie, 0.2, 5000).divergedAt > 0);      // à pas 0,2 : diverge (mesuré au pas 2691)

    CHECK(rotationDrift(euler, 0.05, 20000).divergedAt > 0 && rotationDrift(euler, 0.05, 20000).divergedAt < 5000);   // mesuré : 2251
}

// ---- M6 : toupie pesante de Lagrange (corps symétrique I1 = I2, point fixe, pesanteur) ----

HeavyTopProblem makeTop() {
    HeavyTopProblem t;
    t.mass = 1.0;
    t.inertia = {1.2, 1.2, 0.4};   // moments autour du PIVOT
    t.lever = 0.5;                 // pivot -> centre de masse
    return t;
}

// Intégrales premières : l'énergie E = 1/2 w.Iw + m g l cos(theta), le moment cinétique VERTICAL L_z (la pesanteur n'a pas de couple autour de
// la verticale) et le moment cinétique AXIAL L_3 = I3 w3 (le couple m g l e3 x (...) est perpendiculaire à e3).
void testHeavyTopInvariants() {
    HeavyTopProblem top = makeTop();
    top.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 0.7);
    top.omega0 = {0.8, 1.5, 9.0};
    const RotationState s0{top.q0, top.omega0};
    CHECK_NEAR(top.cosTheta(top.q0), std::cos(0.7), 1e-15);
    // E initiale à la main : 1/2 (1,2 (0,64 + 2,25) + 0,4 * 81) + 1 * g * 0,5 * cos 0,7.
    CHECK_NEAR(top.energy(s0), 0.5 * (1.2 * (0.64 + 2.25) + 0.4 * 81.0) + 0.5 * constants::g0 * std::cos(0.7), 1e-12);
    // Couple à la main : theta = 0 (axe vertical) : aucun couple ; theta = 90 degrés : m g l.
    CHECK_NEAR(top.torque()(Quaternion{}).norm(), 0.0, 1e-15);
    const Vec3 horizontal = top.torque()(Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 0.5 * constants::pi));
    CHECK_NEAR(horizontal.norm(), top.mass * constants::g0 * top.lever, 1e-12);
    CHECK_NEAR(horizontal.z, 0.0, 1e-12);                       // le couple est perpendiculaire à l'axe du corps

    for (double t : {0.5, 2.0, 7.0, 15.0}) {
        const RotationState s = top.reference(t);
        CHECK_NEAR(top.energy(s), top.energy(s0), 1e-10);
        CHECK_NEAR(top.verticalMomentum(s), top.verticalMomentum(s0), 1e-10);
        CHECK_NEAR(top.axialMomentum(s), top.axialMomentum(s0), 1e-10);
        CHECK_NEAR(s.q.norm(), 1.0, 1e-10);
    }
    // La toupie bouge vraiment (nutation) : l'inclinaison n'est pas constante.
    double lo = 1.0, hi = -1.0;
    for (int i = 0; i <= 60; ++i) {
        const double u = top.cosTheta(top.reference(0.05 * i).q);
        lo = std::min(lo, u);
        hi = std::max(hi, u);
    }
    CHECK(hi - lo > 1e-3);
}

// Précession régulière : theta constant, phi' = [I3 w3 -/+ sqrt(I3^2 w3^2 - 4 I1 m g l cos(theta))] / (2 I1 cos(theta)) (branche lente / rapide),
// d'où une solution exacte q(t) = q_z(phi' t) q_x(theta) q_z(psi' t) avec psi' = w3 - phi' cos(theta).
void testHeavyTopSteadyPrecession() {
    HeavyTopProblem top = makeTop();
    for (double theta : {0.6, 1.2}) {
        for (bool slow : {true, false}) {
            const SteadyPrecession sp = top.steadyPrecession(theta, 25.0, slow);
            CHECK(sp.valid);
            // Équation de l'équilibre de theta : I1 phi'^2 cos(theta) - I3 w3 phi' + m g l = 0.
            const double cosT = std::cos(theta);
            CHECK_NEAR(top.inertia.x * sp.phiDot * sp.phiDot * cosT - top.inertia.z * 25.0 * sp.phiDot + top.mass * constants::g0 * top.lever, 0.0, 1e-9);
            CHECK(slow ? sp.phiDot < 2.0 : sp.phiDot > 8.0);   // lente : environ m g l / (I3 w3) = 0,49 ; rapide : environ I3 w3 / (I1 cos) = 8,3 / cos
            top.startSteady(sp);

            for (double t : {0.0, 0.3, 2.0, 9.0}) {
                const RotationState e = top.exactSteady(sp, t), r = top.reference(t);
                CHECK_NEAR(orientationDistance(e.q, r.q), 0.0, 1e-8);
                CHECK_NEAR((e.omega - r.omega).norm(), 0.0, 1e-8);
                CHECK_NEAR(top.cosTheta(r.q), cosT, 1e-8);      // l'inclinaison ne bouge pas (pas de nutation)
            }
        }
    }
    // Axe horizontal (cos theta = 0) : une seule précession, phi' = m g l / (I3 w3), sans branche rapide finie.
    const SteadyPrecession flat = top.steadyPrecession(0.5 * constants::pi, 25.0, true);
    CHECK(flat.valid);
    CHECK_NEAR(flat.phiDot, top.mass * constants::g0 * top.lever / (top.inertia.z * 25.0), 1e-12);
    // Pas de précession régulière si le discriminant est négatif : spin trop faible à theta donné (la toupie tombe).
    CHECK(!top.steadyPrecession(0.3, 2.0, true).valid);
}

// Nutation de la toupie : u = cos(theta) oscille entre deux racines u1 < u2 du polynôme cubique f(u) = u'^2 (Landau-Lifchitz, § 35) :
//   f(u) = (2/I1)(E' - m g l u)(1 - u^2) - (L_z - L_3 u)^2 / I1^2,   E' = E - L_3^2 / (2 I3),
// de racine supérieure u3 > 1. La période de nutation vaut T = 4 K(k) / sqrt(beta (u3 - u1)), beta = 2 m g l / I1, k^2 = (u2 - u1) / (u3 - u1).
void testNutationPeriod() {
    HeavyTopProblem top = makeTop();
    top.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 0.9);
    top.omega0 = {0.0, 2.0, 20.0};                               // theta' = 0 au départ : on part d'un point de rebroussement
    const NutationRange n = top.nutation();
    CHECK(n.valid);
    CHECK(n.u1 < n.u2 && n.u2 < 1.0 && n.u3 > 1.0);
    const double u0 = top.cosTheta(top.q0);
    CHECK(std::abs(u0 - n.u1) < 1e-10 || std::abs(u0 - n.u2) < 1e-10);   // le départ (theta' = 0) est une des deux racines

    const double period = top.nutationPeriod();
    CHECK(period > 0.05 && period < 5.0);
    // À une demi-période l'autre point de rebroussement, à une période retour au départ.
    const double other = std::abs(u0 - n.u1) < 1e-10 ? n.u2 : n.u1;
    CHECK_NEAR(top.cosTheta(top.reference(0.5 * period).q), other, 1e-11);      // mesuré : 2e-14 (T = 1,069185)
    CHECK_NEAR(top.cosTheta(top.reference(period).q), u0, 1e-11);
    CHECK_NEAR(top.cosTheta(top.reference(3.0 * period).q), u0, 1e-10);
    // Ni plus tôt : au tiers de période la toupie est ailleurs.
    CHECK(std::abs(top.cosTheta(top.reference(period / 3.0).q) - u0) > 1e-3);
    // Entre les deux, u reste dans [u1, u2].
    for (double f : {0.1, 0.27, 0.5, 0.71, 0.93}) {
        const double u = top.cosTheta(top.reference(f * period).q);
        CHECK(u >= n.u1 - 1e-9 && u <= n.u2 + 1e-9);
    }
}

// Toupie endormie (theta = 0) : stable si I3^2 w3^2 > 4 I1 m g l. En dessous, une petite inclinaison croît comme exp(gamma t),
// gamma = sqrt(4 I1 m g l - I3^2 w3^2) / (2 I1) (de I1 zeta'' + i I3 w3 zeta' - m g l zeta = 0, zeta = theta e^(i phi), linéarisée).
void testSleepingTop() {
    HeavyTopProblem top = makeTop();
    const double critical = top.sleepingCriticalSpin();
    CHECK_NEAR(critical, 2.0 * std::sqrt(top.inertia.x * top.mass * constants::g0 * top.lever) / top.inertia.z, 1e-12);
    CHECK(critical > 11.0 && critical < 13.0);

    // Au-dessus du seuil : une inclinaison de 0,01 rad reste petite (oscillation bornée).
    const double bounds[] = {0.025, 0.015, 0.012};              // mesuré : 0,0203 ; 0,0134 ; 0,0106 (plus le spin est grand, moins la toupie s'écarte)
    int index = 0;
    for (double factor : {1.15, 1.5, 3.0}) {
        top.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 0.01);
        top.omega0 = {0.0, 0.0, factor * critical};
        double worst = 0.0;
        walkReference(top.rhs(), top.initialState(), 0.05, 400, [&](double, const RotationState& s) { worst = std::max(worst, std::acos(std::min(1.0, top.cosTheta(s.q)))); });
        CHECK(worst < bounds[index++]);
    }
    // En dessous : croissance exponentielle au taux gamma (mesuré par la pente de ln theta), puis chute.
    const double spin = 0.7 * critical;
    const double gamma = std::sqrt(4.0 * top.inertia.x * top.mass * constants::g0 * top.lever - top.inertia.z * top.inertia.z * spin * spin) / (2.0 * top.inertia.x);
    top.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 1e-6);
    top.omega0 = {0.0, 0.0, spin};
    const double a = std::acos(top.cosTheta(top.reference(3.0).q)), b = std::acos(top.cosTheta(top.reference(5.0).q));
    CHECK_NEAR(std::log(b / a) / 2.0, gamma, 1e-3 * gamma);      // mesuré : 1,44357 pour 1,44358 prévu
    // Puis la grande nutation : la toupie s'incline jusqu'à environ 1,25 rad (mesuré) mais, fixée au pivot et sans perte, elle REMONTE vers la
    // verticale (theta(20 s) ~ 1e-5) : elle ne « tombe » pas, l'énergie est conservée.
    double maxTheta = 0.0, thetaAtEnd = 0.0;
    walkReference(top.rhs(), top.initialState(), 0.5, 40, [&](double, const RotationState& s) {
        thetaAtEnd = std::acos(std::min(1.0, top.cosTheta(s.q)));
        maxTheta = std::max(maxTheta, thetaAtEnd);
    });
    CHECK(maxTheta > 1.0 && maxTheta < 1.6);
    CHECK(thetaAtEnd < 0.01);
}

// Les intégrateurs avec un couple (la toupie) : mêmes ordres que sans couple, et le découpage symplectique garde le moment vertical L_z à
// l'arrondi près (le couple de la pesanteur n'a pas de composante verticale : les « coups » ne le changent pas).
void testRotationIntegratorsWithTorque() {
    HeavyTopProblem top = makeTop();
    top.q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, 0.9);
    top.omega0 = {0.0, 2.0, 20.0};
    const RotationState s0{top.q0, top.omega0};
    const TorqueFunction torque = top.torque();
    const double tEnd = 2.0;                                    // près de deux périodes de nutation (1,069)
    const RotationState ref = top.reference(tEnd);
    auto error = [&](RotationIntegrator& in, int n) { return rotationDistance(integrateRotation(in, top.inertia, torque, s0, tEnd, n), ref); };

    RK4Rotation rk4;
    LieHeunRotation lie;
    SplittingRotation splitting;
    // De 200 à 800 pas (pré-asymptotique en dessous) : ordre 4 pour RK4 (263 mesuré), ordre 2 pour Heun (16,3) et le découpage (16,0).
    const double rk4Ratio = error(rk4, 200) / error(rk4, 800);
    CHECK(rk4Ratio > 245.0 && rk4Ratio < 285.0);
    const double lieRatio = error(lie, 200) / error(lie, 800);
    CHECK(lieRatio > 15.0 && lieRatio < 17.5);
    const double splitRatio = error(splitting, 200) / error(splitting, 800);
    CHECK(splitRatio > 15.0 && splitRatio < 16.8);
    CHECK(error(rk4, 800) < 5e-3 * error(splitting, 800));      // RK4 est nettement le plus précis sur un horizon court (mesuré : rapport 1,3e-3)

    // Euler explicite explose dès 50 pas (h w3 = 0,8) ; les autres restent finis.
    EulerRotation euler;
    const RotationState bad = integrateRotation(euler, top.inertia, torque, s0, tEnd, 50);
    CHECK(!std::isfinite(bad.q.w) || bad.omega.norm() > 1e3);

    // 100 s de nutation à h = 0,005 : L_z (et E) selon le schéma.
    const double lz0 = top.verticalMomentum(s0), e0 = top.energy(s0);
    auto drift = [&](RotationIntegrator& in, double& maxLz, double& maxE, double& maxQ) {
        RotationState s = s0;
        maxLz = maxE = maxQ = 0.0;
        for (int i = 0; i < 20000; ++i) {
            in.step(top.inertia, torque, s, 0.005);
            maxLz = std::max(maxLz, std::abs(top.verticalMomentum(s) - lz0) / std::abs(lz0));
            maxE = std::max(maxE, std::abs(top.energy(s) - e0) / std::abs(e0));
            maxQ = std::max(maxQ, std::abs(s.q.norm() - 1.0));
        }
    };
    double lzSplit, eSplit, qSplit, lzRk4, eRk4, qRk4, lzLie, eLie, qLie;
    drift(splitting, lzSplit, eSplit, qSplit);
    drift(rk4, lzRk4, eRk4, qRk4);
    drift(lie, lzLie, eLie, qLie);
    CHECK(lzSplit < 1e-11);                                     // mesuré : 2,8e-13
    CHECK(eSplit < 1e-4);                                       // mesuré : 3,9e-5 (borné)
    CHECK(qSplit < 1e-14);
    CHECK(lzRk4 > 1e-6 && lzRk4 < 1e-4);                        // mesuré : 5,1e-6 : RK4 laisse dériver L_z (1e7 fois plus que le découpage)
    CHECK(eRk4 < 5e-6);                                         // mesuré : 7,5e-7 : à ce pas RK4 conserve mieux l'énergie sur 100 s
    CHECK(lzLie > 1e-3);                                        // mesuré : 7,5e-3 : le groupe de Lie d'ordre 2 dérive en L_z ...
    CHECK(eLie > 1e-4 && eLie < 0.05);                          // ... et en énergie (mesuré 2,3e-3)
    CHECK(qLie < 1e-14 && qRk4 < 1e-14);

    // Précession régulière sur 10 s en 2000 pas : theta reste constant au moins à 1e-5 près pour RK4 et le découpage ; Heun s'en écarte de 1e-3.
    HeavyTopProblem st = makeTop();
    const SteadyPrecession sp = st.steadyPrecession(0.8, 25.0, true);
    st.startSteady(sp);
    const RotationState exact = st.exactSteady(sp, 10.0);
    const RotationState a = integrateRotation(rk4, st.inertia, st.torque(), {st.q0, st.omega0}, 10.0, 2000);
    const RotationState b = integrateRotation(splitting, st.inertia, st.torque(), {st.q0, st.omega0}, 10.0, 2000);
    const RotationState c = integrateRotation(lie, st.inertia, st.torque(), {st.q0, st.omega0}, 10.0, 2000);
    CHECK_NEAR(rotationDistance(a, exact), 0.0, 1e-4);          // mesuré : 5e-5
    CHECK_NEAR(rotationDistance(b, exact), 0.0, 5e-4);          // mesuré : 2e-4
    CHECK(rotationDistance(c, exact) > 1e-3);                   // mesuré : 8e-3
    CHECK_NEAR(st.cosTheta(a.q), std::cos(0.8), 1e-5);          // mesuré : 1,3e-6
    CHECK_NEAR(st.cosTheta(b.q), std::cos(0.8), 5e-5);          // mesuré : 8e-6
    CHECK(std::abs(st.cosTheta(c.q) - std::cos(0.8)) > 1e-4);   // mesuré : 8,8e-4
}

void testPeriapsisTracker() {
    // Solution exacte échantillonnée grossièrement (50 points par période, décalés pour ne pas tomber pile sur un
    // périastre) : exactement 5 périastres sur 5,3 périodes, tous à l'angle 0 (orbite fermée, pas de précession).
    KeplerProblem p;
    p.e = 0.5;
    const double T = p.period();
    PeriapsisTracker tracker;
    CHECK(std::isnan(tracker.precessionPerOrbit()));
    int detections = 0;
    for (int i = 0; i < 265; ++i) {
        const double t = (i + 0.3) * T / 50.0;
        if (tracker.update(t, p.exact(t))) ++detections;
    }
    CHECK(detections == 5);
    CHECK(tracker.count() == 5);
    CHECK_NEAR(tracker.precession(), 0.0, 1e-6);
    CHECK_NEAR(tracker.precessionPerOrbit(), 0.0, 1e-6);

    // Une orbite qui tourne de delta par tour : on fait tourner l'état exact de k*delta à chaque période.
    // La précession mesurée doit valoir delta par orbite, y compris au-delà de pi (angle déroulé).
    const double delta = 0.9;
    PeriapsisTracker turning;
    for (int i = 0; i < 400; ++i) {
        const double t = (i + 0.3) * T / 50.0;
        State y = p.exact(t);
        const double rot = delta * t / T;  // rotation progressive du plan de l'orbite
        const double c = std::cos(rot), s = std::sin(rot);
        State z = y;
        z[0] = c * y[0] - s * y[1]; z[1] = s * y[0] + c * y[1];
        // vitesse vraie = vitesse tournée + vitesse d'entraînement (omega x r), avec omega = delta / T
        const double w = delta / T;
        z[3] = c * y[3] - s * y[4] - w * z[1]; z[4] = s * y[3] + c * y[4] + w * z[0];
        turning.update(t, z);
    }
    CHECK(turning.count() == 7);
    CHECK_NEAR(turning.precessionPerOrbit(), delta, 1e-5);
    CHECK(turning.precession() > constants::pi);  // 6 orbites * 0.9 = 5,4 rad : l'angle n'est pas replié dans (-pi, pi]

    turning.reset();
    CHECK(turning.count() == 0);
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
    testElliptic();
    testPendulumExact();
    testPendulumConvergence();
    testPendulumEnergy();
    testDoublePendulumEquations();
    testDoublePendulumNumerics();
    testChaosVersusRegular();
    testLyapunovFit();
    testKeplerEquation();
    testKeplerExact();
    testPeriapsisTracker();
    testKeplerConvergence();
    testKeplerEnergyAndPrecession();
    testKeplerPrecessionTheory();
    testKeplerAdaptiveStep();
    testNBodyAccelerations();
    testNBodyDiagnostics();
    testLagrangeTriangle();
    testTwoBodyMatchesKepler();
    testFigureEight();
    testNBodyConvergence();
    testNBodyConservation();
    testNBodyChaos();
    testPlummer();
    testEventDetection();
    testAdvanceBudget();
    testInclineExact();
    testInclineEventDriven();
    testInclineNaive();
    testInclineRegularized();
    testCollide1D();
    testCollideSpheres();
    testTwoBallExact();
    testBounceExact();
    testBounceEventDriven();
    testBounceNaive();
    testTwoBallContact();
    testHertzReference();
    testHertzTwoBalls();
    testCradleInvariants();
    testThreeBallFamily();
    testHertzThreeBalls();
    testHertzChain();
    testSequentialImpulses();
    testSequentialOrder();
    testCradleGap();
    testCradleConvergence();
    testHertzDamping();
    testInertiaTensors();
    testEulerEquations();
    testFreeBodyInvariants();
    testSymmetricTopExact();
    testJacobiFunctions();
    testAsymmetricExact();
    testIntermediateAxisInstability();
    testRotationIntegratorBasics();
    testSplittingConservation();
    testRotationOrders();
    testRotationDrift();
    testHeavyTopInvariants();
    testHeavyTopSteadyPrecession();
    testNutationPeriod();
    testSleepingTop();
    testRotationIntegratorsWithTorque();

    if (g_failures == 0) {
        std::puts("test_core : OK");
        return 0;
    }
    std::fprintf(stderr, "test_core : %d échec(s)\n", g_failures);
    return 1;
}
