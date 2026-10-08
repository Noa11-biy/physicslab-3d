#include "physicslab/mechanics/Cradle.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "physicslab/mechanics/Collision.hpp"

namespace pl {

namespace hertz {

double contactConstant() {
    static const double c = 0.8 * std::tgamma(0.4) * std::tgamma(0.5) / std::tgamma(0.9);
    return c;
}

double force(double stiffness, double delta) { return delta > 0.0 ? stiffness * delta * std::sqrt(delta) : 0.0; }

double potentialEnergy(double stiffness, double delta) {
    return delta > 0.0 ? 0.4 * stiffness * delta * delta * std::sqrt(delta) : 0.0;
}

double maxCompression(double reducedMass, double approachSpeed, double stiffness) {
    return std::pow(5.0 * reducedMass * approachSpeed * approachSpeed / (4.0 * stiffness), 0.4);
}

double contactDuration(double reducedMass, double approachSpeed, double stiffness) {
    return contactConstant() * maxCompression(reducedMass, approachSpeed, stiffness) / approachSpeed;
}

}  // namespace hertz

namespace cradle {

std::array<double, 3> threeBallFamily(double x, double speed) {
    x = std::clamp(x, 0.0, 1.0 / 3.0);
    const double sum = 1.0 + x;
    const double root = std::sqrt(std::max(0.0, sum * (1.0 - 3.0 * x)));   // discriminant (1 + x)(1 - 3x)
    return {-x * speed, 0.5 * (sum - root) * speed, 0.5 * (sum + root) * speed};
}

ImpulseResult sequentialImpulses(std::vector<double> v, double restitution, ResolveOrder order) {
    constexpr int kMaxSweeps = 100000;
    ImpulseResult r;
    double scale = 0.0;
    for (double x : v) scale = std::max(scale, std::abs(x));
    const double tol = 1e-12 * scale;
    const int n = static_cast<int>(v.size());

    bool moved = true;
    while (moved && r.sweeps < kMaxSweeps) {
        moved = false;
        ++r.sweeps;
        for (int k = 0; k + 1 < n; ++k) {
            const int i = order == ResolveOrder::LeftToRight ? k : n - 2 - k;
            if (v[i] - v[i + 1] <= tol) continue;   // la paire ne se rapproche pas
            const collision::Result1D c = collision::collide1D(1.0, v[i], 1.0, v[i + 1], restitution);
            v[i] = c.v1;
            v[i + 1] = c.v2;
            r.energyLoss += c.energyLoss;
            ++r.collisions;
            moved = true;
        }
    }
    r.converged = !moved;
    r.velocities = std::move(v);
    return r;
}

}  // namespace cradle

State CradleProblem::initialState() const {
    State y(2 * balls, 0.0);
    for (int i = 0; i < balls; ++i) {
        y[i] = 2.0 * radius * i;
        if (i < launched) {
            y[i] -= gap;
            y[balls + i] = speed;
        }
    }
    return y;
}

double CradleProblem::compression(const State& y, int pair) const { return 2.0 * radius - (y[pair + 1] - y[pair]); }

OdeFunction CradleProblem::rhs() const {
    const CradleProblem p = *this;   // copie : la fonction survit au problème
    return [p](double, const State& y, State& dydt) {
        const int n = p.balls;
        for (int i = 0; i < n; ++i) {
            dydt[i] = y[n + i];
            dydt[n + i] = 0.0;
        }
        for (int i = 0; i + 1 < n; ++i) {
            // Force du contact i : repousse la bille i vers la gauche et la bille i+1 vers la droite (3e loi).
            const double delta = p.compression(y, i);
            double f = hertz::force(p.stiffness, delta) / p.mass;
            if (p.damping > 0.0 && delta > 0.0) {
                const double squeezeSpeed = y[n + i] - y[n + i + 1];   // delta' : les billes se rapprochent si positive
                f *= std::max(0.0, 1.0 + 1.5 * p.damping * squeezeSpeed);
            }
            dydt[n + i] -= f;
            dydt[n + i + 1] += f;
        }
    };
}

double CradleProblem::kineticEnergy(const State& y) const {
    double sum = 0.0;
    for (int i = 0; i < balls; ++i) sum += y[balls + i] * y[balls + i];
    return 0.5 * mass * sum;
}

double CradleProblem::potentialEnergy(const State& y) const {
    double sum = 0.0;
    for (int i = 0; i + 1 < balls; ++i) sum += hertz::potentialEnergy(stiffness, compression(y, i));
    return sum;
}

double CradleProblem::momentum(const State& y) const {
    double sum = 0.0;
    for (int i = 0; i < balls; ++i) sum += y[balls + i];
    return mass * sum;
}

bool CradleProblem::collisionOver(const State& y) const {
    const double lengthTol = 1e-12 * 2.0 * radius, speedTol = 1e-12 * std::abs(speed);   // bruit d'arrondi
    for (int i = 0; i + 1 < balls; ++i) {
        if (compression(y, i) > lengthTol) return false;
        if (y[balls + i] > y[balls + i + 1] + speedTol) return false;
    }
    return true;
}

double CradleProblem::suggestedStep() const { return hertz::contactDuration(0.5 * mass, speed, stiffness) / 100.0; }

CradleOutcome CradleProblem::run(Solver& solver, double dt, double tMax) const {
    const OdeFunction f = rhs();
    State y = initialState();
    CradleOutcome out;
    double t = 0.0;
    out.finished = collisionOver(y);
    while (!out.finished && t < tMax) {
        const double h = std::min(dt, tMax - t);
        const double advanced = advance(solver, f, t, y, h, 1000);   // budget : on ne boucle jamais indéfiniment
        t += advanced;
        ++out.steps;
        for (int i = 0; i + 1 < balls; ++i) out.maxCompression = std::max(out.maxCompression, compression(y, i));
        if (advanced < h) break;
        out.finished = collisionOver(y);
    }
    out.time = t;
    out.velocities.assign(y.begin() + balls, y.end());
    return out;
}

State CradleProblem::reference(double t) const {
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    State y = initialState();
    advance(rk, rhs(), 0.0, y, t, 1000000);
    return y;
}

double cradleError(const CradleProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    for (int i = 0; i < steps; ++i) advance(solver, f, i * dt, y, dt, 1000);
    const State ref = problem.reference(tEnd);
    const double length = problem.speed * hertz::contactDuration(0.5 * problem.mass, problem.speed, problem.stiffness);
    double sum = 0.0;
    const int n = problem.balls;
    for (int i = 0; i < n; ++i) {
        const double dx = (y[i] - ref[i]) / length, dv = (y[n + i] - ref[n + i]) / problem.speed;
        sum += dx * dx + dv * dv;
    }
    return std::sqrt(sum);
}

}  // namespace pl
