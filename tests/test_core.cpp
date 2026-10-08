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

}  // namespace

int main() {
    testVec3();
    testMat3();
    testQuaternion();
    testEulerOrder();
    testFreeFall();

    if (g_failures == 0) {
        std::puts("test_core : OK");
        return 0;
    }
    std::fprintf(stderr, "test_core : %d échec(s)\n", g_failures);
    return 1;
}
