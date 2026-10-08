// Validation du socle : maths de base + comparaison du solveur à des solutions analytiques.
// Pas de framework externe : un CHECK minimal suffit pour M0.
#include <algorithm>
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
#include "physicslab/mechanics/DoublePendulum.hpp"
#include "physicslab/mechanics/Friction.hpp"
#include "physicslab/mechanics/Kepler.hpp"
#include "physicslab/mechanics/NBody.hpp"
#include "physicslab/mechanics/Oscillator.hpp"
#include "physicslab/mechanics/Pendulum.hpp"
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
    testEventDetection();
    testAdvanceBudget();
    testInclineExact();
    testInclineEventDriven();
    testInclineNaive();
    testInclineRegularized();

    if (g_failures == 0) {
        std::puts("test_core : OK");
        return 0;
    }
    std::fprintf(stderr, "test_core : %d échec(s)\n", g_failures);
    return 1;
}
