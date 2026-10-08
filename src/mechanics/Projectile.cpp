#include "physicslab/mechanics/Projectile.hpp"

#include <cmath>
#include <limits>

namespace pl {
namespace {

// phi = (1 - e^{-kt}) / k  et  psi = (t - phi) / k.
// Pour k t petit, la forme exacte perd tous ses chiffres significatifs (0/0) : on utilise
// alors le développement de Taylor (erreur relative < 1e-14 pour |k t| < 1e-3).
void dragFunctions(double k, double t, double& phi, double& psi) {
    const double x = k * t;
    if (std::abs(x) < 1e-3) {
        phi = t * (1.0 - x / 2.0 + x * x / 6.0 - x * x * x / 24.0);
        psi = t * t * (0.5 - x / 6.0 + x * x / 24.0 - x * x * x / 120.0);
    } else {
        phi = -std::expm1(-x) / k;
        psi = (t - phi) / k;
    }
}

}  // namespace

World ProjectileProblem::makeWorld() const {
    World w;
    w.gravity = gravity;
    w.linearDrag = linearDrag;
    Particle p;
    p.mass = mass;
    p.position = position0;
    p.velocity = velocity0;
    w.particles.push_back(p);
    return w;
}

Vec3 ProjectileProblem::position(double t) const {
    double phi, psi;
    dragFunctions(linearDrag / mass, t, phi, psi);
    return position0 + phi * velocity0 + psi * gravity;
}

Vec3 ProjectileProblem::velocity(double t) const {
    const double k = linearDrag / mass;
    double phi, psi;
    dragFunctions(k, t, phi, psi);
    return std::exp(-k * t) * velocity0 + phi * gravity;
}

double ProjectileProblem::energy(double t) const {
    return 0.5 * mass * velocity(t).norm2() - mass * dot(gravity, position(t));
}

double ProjectileProblem::landingTime() const {
    if (gravity.y >= 0.0) return std::numeric_limits<double>::infinity();

    auto height = [this](double t) { return position(t).y; };

    // On part juste après 0 pour ne pas confondre "départ depuis le sol" et "retombée".
    double lo = 1e-9;
    if (height(lo) <= 0.0) return 0.0;

    double hi = 1.0;
    for (int i = 0; i < 200 && height(hi) > 0.0; ++i) hi *= 2.0;

    for (int i = 0; i < 200; ++i) {  // dichotomie
        const double mid = 0.5 * (lo + hi);
        (height(mid) > 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double integrationError(const ProjectileProblem& problem, Solver& solver, int steps, double tEnd) {
    World w = problem.makeWorld();
    const double dt = tEnd / steps;
    for (int i = 0; i < steps; ++i) w.step(solver, dt);
    return (w.particles[0].position - problem.position(w.time)).norm();
}

}  // namespace pl
