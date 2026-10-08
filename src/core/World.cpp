#include "physicslab/core/World.hpp"

namespace pl {

double World::step(Solver& solver, double dt) {
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

    // dx/dt = v ; dv/dt = g (même accélération pour toutes les particules).
    const Vec3 g = gravity;
    const OdeFunction f = [split, n, g](double, const State& y, State& dydt) {
        for (std::size_t i = 0; i < split; ++i) dydt[i] = y[split + i];
        for (std::size_t i = 0; i < n; ++i) {
            dydt[split + 3 * i + 0] = g.x;
            dydt[split + 3 * i + 1] = g.y;
            dydt[split + 3 * i + 2] = g.z;
        }
    };

    const double advanced = solver.step(f, time, state_, dt);

    for (std::size_t i = 0; i < n; ++i) {
        Particle& p = particles[i];
        p.position = {state_[3 * i + 0], state_[3 * i + 1], state_[3 * i + 2]};
        p.velocity = {state_[split + 3 * i + 0], state_[split + 3 * i + 1], state_[split + 3 * i + 2]};
    }
    time += advanced;
    return advanced;
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
