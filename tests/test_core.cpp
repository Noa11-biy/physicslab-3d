// Validation du socle : maths de base + comparaison du solveur à des solutions analytiques.
// Pas de framework externe : un CHECK minimal suffit pour M0.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Mat3.hpp"
#include "physicslab/core/Quaternion.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"
#include "physicslab/core/World.hpp"
#include "physicslab/mechanics/DoublePendulum.hpp"
#include "physicslab/mechanics/Kepler.hpp"
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

    if (g_failures == 0) {
        std::puts("test_core : OK");
        return 0;
    }
    std::fprintf(stderr, "test_core : %d échec(s)\n", g_failures);
    return 1;
}
