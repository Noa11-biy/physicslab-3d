#include "physicslab/mechanics/Collision.hpp"

#include <algorithm>
#include <cmath>

namespace pl {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

collision::Result1D collision::collide1D(double m1, double v1, double m2, double v2, double e) {
    const double vRel = v1 - v2;
    if (!(vRel > 0.0)) return {v1, v2, 0.0};  // ils s'éloignent déjà
    const double mu = m1 * m2 / (m1 + m2);
    const double J = (1.0 + e) * mu * vRel;   // impulsion échangée
    return {v1 - J / m1, v2 + J / m2, 0.5 * (1.0 - e * e) * mu * vRel * vRel};
}

collision::Result3D collision::collideSpheres(double m1, const Vec3& v1, double m2, const Vec3& v2, const Vec3& n, double e) {
    const double vn = dot(v1 - v2, n);
    if (!(vn > 0.0)) return {v1, v2, 0.0, 0.0};
    const double mu = m1 * m2 / (m1 + m2);
    const double J = (1.0 + e) * mu * vn;
    return {v1 - (J / m1) * n, v2 + (J / m2) * n, J, 0.5 * (1.0 - e * e) * mu * vn * vn};
}

// ------------------------------- deux billes ---------------------------------------

State TwoBallProblem::initialState() const {
    return {p1.x, p1.y, p2.x, p2.y, v1.x, v1.y, v2.x, v2.y};
}

OdeFunction TwoBallProblem::rhs() const {
    return [](double, const State& y, State& d) {
        for (int i = 0; i < 4; ++i) {
            d[i] = y[4 + i];  // x' = v
            d[4 + i] = 0.0;   // pas de force : mouvement libre
        }
    };
}

double TwoBallProblem::gap(const State& y) const { return std::hypot(y[0] - y[2], y[1] - y[3]) - (r1 + r2); }

double TwoBallProblem::kineticEnergy(const State& y) const {
    return 0.5 * m1 * (y[4] * y[4] + y[5] * y[5]) + 0.5 * m2 * (y[6] * y[6] + y[7] * y[7]);
}

Vec3 TwoBallProblem::momentum(const State& y) const {
    return {m1 * y[4] + m2 * y[6], m1 * y[5] + m2 * y[7], 0.0};
}

double TwoBallProblem::collisionTime() const {
    // |d + w t| = R avec d = p1 - p2, w = v1 - v2 : (w.w) t^2 + 2 (d.w) t + (d.d - R^2) = 0.
    const double dx = p1.x - p2.x, dy = p1.y - p2.y, wx = v1.x - v2.x, wy = v1.y - v2.y, R = r1 + r2;
    const double a = wx * wx + wy * wy, b = dx * wx + dy * wy, c = dx * dx + dy * dy - R * R;
    if (!(a > 0.0) || !(b < 0.0)) return c <= 0.0 && b < 0.0 ? 0.0 : kInf;   // immobiles l'une par rapport à l'autre, ou elles s'éloignent
    const double disc = b * b - a * c;
    if (disc < 0.0) return kInf;                                               // elles se manquent
    if (c <= 0.0) return 0.0;                                                  // déjà en contact, en approche
    return c / (-b + std::sqrt(disc));  // plus petite racine, sans soustraction de deux nombres proches
}

namespace {

// Vitesses après le choc de problème `p` (contact à l'instant tc), normale de 1 vers 2.
collision::Result3D afterCollision(const TwoBallProblem& p, double tc) {
    const double x1 = p.p1.x + p.v1.x * tc, y1 = p.p1.y + p.v1.y * tc, x2 = p.p2.x + p.v2.x * tc, y2 = p.p2.y + p.v2.y * tc;
    const Vec3 n = Vec3{x2 - x1, y2 - y1, 0.0}.normalized();
    return collision::collideSpheres(p.m1, p.v1, p.m2, p.v2, n, p.restitution);
}

}  // namespace

State TwoBallProblem::exact(double t) const {
    const double tc = collisionTime();
    if (!(t > tc)) {
        return {p1.x + v1.x * t, p1.y + v1.y * t, p2.x + v2.x * t, p2.y + v2.y * t, v1.x, v1.y, v2.x, v2.y};
    }
    const collision::Result3D r = afterCollision(*this, tc);
    const double dt = t - tc;
    return {p1.x + v1.x * tc + r.v1.x * dt, p1.y + v1.y * tc + r.v1.y * dt, p2.x + v2.x * tc + r.v2.x * dt,
            p2.y + v2.y * tc + r.v2.y * dt, r.v1.x, r.v1.y, r.v2.x, r.v2.y};
}

double TwoBallProblem::outgoingAngle() const {
    const double tc = collisionTime();
    if (!std::isfinite(tc)) return kNaN;
    const collision::Result3D r = afterCollision(*this, tc);
    const double n1 = r.v1.norm(), n2 = r.v2.norm();
    if (!(n1 > 0.0) || !(n2 > 0.0)) return kNaN;  // une bille est arrêtée : angle indéfini
    return std::acos(std::clamp(dot(r.v1, r.v2) / (n1 * n2), -1.0, 1.0));
}

TwoBallRun::TwoBallRun(const TwoBallProblem& problem, ContactModel model)
    : problem_(problem), model_(model), y_(problem.initialState()) {}

void TwoBallRun::collide() {
    const Vec3 n = Vec3{y_[2] - y_[0], y_[3] - y_[1], 0.0}.normalized();
    const collision::Result3D r =
        collision::collideSpheres(problem_.m1, {y_[4], y_[5], 0.0}, problem_.m2, {y_[6], y_[7], 0.0}, n, problem_.restitution);
    y_[4] = r.v1.x; y_[5] = r.v1.y; y_[6] = r.v2.x; y_[7] = r.v2.y;
    if (collisions_++ == 0) firstCollision_ = time_;
}

void TwoBallRun::advance(Solver& solver, double dt) {
    const OdeFunction f = problem_.rhs();
    if (model_ == ContactModel::Naive) {
        // Un pas, puis on regarde : si les billes se sont interpénétrées en s'approchant, on applique le choc sur l'état courant.
        pl::advance(solver, f, time_, y_, dt);
        time_ += dt;
        if (problem_.gap(y_) < 0.0) {
            const Vec3 n = Vec3{y_[2] - y_[0], y_[3] - y_[1], 0.0}.normalized();
            if (dot(Vec3{y_[4] - y_[6], y_[5] - y_[7], 0.0}, n) > 0.0) collide();
        }
        return;
    }

    // Événement : g = max(distance entre les surfaces, 0) change de signe exactement au contact ; après le choc g reste nulle (elles s'éloignent).
    const EventFunction touching = [this](double, const State& y) { return std::max(problem_.gap(y), 0.0); };
    double remaining = dt;
    for (int guard = 0; guard < 8 && remaining > 1e-15 * dt; ++guard) {
        const EventStep r = advanceToEvent(solver, f, touching, time_, y_, remaining);
        time_ += r.elapsed;
        remaining -= r.elapsed;
        if (r.event) collide();
    }
}

double twoBallError(const TwoBallProblem& problem, ContactModel model, Solver& solver, int steps, double tEnd) {
    TwoBallRun run(problem, model);
    const double dt = tEnd / steps;
    for (int i = 0; i < steps; ++i) run.advance(solver, dt);
    const State ex = problem.exact(tEnd);
    double err = 0.0;
    for (int i = 0; i < 4; ++i) err = std::max(err, std::abs(run.state()[i] - ex[i]));
    return err;
}

}  // namespace pl
