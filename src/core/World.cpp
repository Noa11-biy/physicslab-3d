#include "physicslab/core/World.hpp"

namespace pl {

double World::step(Solver& solver, double dt) {
    if (!(dt > 0.0)) return 0.0;

    const std::size_t n = particles.size();
    const std::size_t split = 3 * n;  // state_ = [positions | vitesses]

    state_.resize(6 * n);
    for (std::size_t i = 0; i < n; ++i) {
        const Particle& p = particles[i];
        state_[3 * i + 0] = p.position.x;
        state_[3 * i + 1] = p.position.y;
        state_[3 * i + 2] = p.position.z;
        state_[split + 3 * i + 0] = p.velocity.x;
        state_[split + 3 * i + 1] = p.velocity.y;
        state_[split + 3 * i + 2] = p.velocity.z;
    }

    // dx/dt = v ; dv/dt = g - (b/m) v
    const Vec3 g = gravity;
    const double drag = linearDrag;
    const std::vector<Particle>* ps = &particles;  // seules les masses sont lues, elles ne changent pas pendant le pas
    const OdeFunction f = [split, n, g, drag, ps](double, const State& y, State& dydt) {
        for (std::size_t i = 0; i < split; ++i) dydt[i] = y[split + i];
        for (std::size_t i = 0; i < n; ++i) {
            const double k = drag / (*ps)[i].mass;
            dydt[split + 3 * i + 0] = g.x - k * y[split + 3 * i + 0];
            dydt[split + 3 * i + 1] = g.y - k * y[split + 3 * i + 1];
            dydt[split + 3 * i + 2] = g.z - k * y[split + 3 * i + 2];
        }
    };

    // Un solveur adaptatif peut avancer de moins que demandé : on recommence jusqu'à couvrir dt.
    double elapsed = 0.0;
    while (dt - elapsed > 1e-12 * dt) {
        const double advanced = solver.step(f, time + elapsed, state_, dt - elapsed);
        if (!(advanced > 0.0)) break;  // garde-fou : un solveur défaillant ne doit pas boucler
        elapsed += advanced;
    }

    for (std::size_t i = 0; i < n; ++i) {
        Particle& p = particles[i];
        p.position = {state_[3 * i + 0], state_[3 * i + 1], state_[3 * i + 2]};
        p.velocity = {state_[split + 3 * i + 0], state_[split + 3 * i + 1], state_[split + 3 * i + 2]};
    }
    time += elapsed;
    return elapsed;
}

Invariants World::invariants() const {
    Invariants inv;
    for (const Particle& p : particles) {
        inv.kinetic += 0.5 * p.mass * p.velocity.norm2();
        inv.potential -= p.mass * dot(gravity, p.position);  // U = -m g.r
        inv.momentum += p.mass * p.velocity;
        inv.angularMomentum += p.mass * cross(p.position, p.velocity);
    }
    return inv;
}

}  // namespace pl
