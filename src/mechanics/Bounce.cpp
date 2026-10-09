#include "physicslab/mechanics/Bounce.hpp"

#include <algorithm>
#include <cmath>

#include "physicslab/mechanics/Projectile.hpp"

namespace pl {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// phi(tau) = (1 - e^{-k tau}) / k : intégrale de e^{-k s}, limite tau quand k -> 0.
double phi(double k, double tau) {
    const double x = k * tau;
    return x == 0.0 ? tau : tau * (-std::expm1(-x) / x);
}

// Parcourt la suite exacte des vols et des rebonds. S'arrête au temps tTarget (état renvoyé), ou quand `maxImpacts` instants d'impact
// ont été collectés, ou au repos. restStart reçoit l'instant du repos (l'infini s'il n'a pas eu lieu).
BounceState trace(const BounceProblem& p, double tTarget, std::vector<double>* impacts, std::size_t maxImpacts, double* restStart) {
    const double k = p.drag, e = p.restitution;
    double tCur = 0.0, x = p.x0, y = p.y0, vx = p.vx0, vy = p.vy0;
    int bounces = 0;
    if (restStart) *restStart = kInf;

    for (int guard = 0; guard < 200000; ++guard) {
        // Au sol, vitesse de rebond sous le seuil : la balle cesse de rebondir et glisse (x' = vx, vx' = -k vx).
        if (y <= 0.0 && vy <= p.restSpeed) {
            if (restStart) *restStart = tCur;
            const double tau = std::isfinite(tTarget) ? std::max(tTarget - tCur, 0.0) : 0.0;
            return {x + vx * phi(k, tau), 0.0, vx * std::exp(-k * tau), 0.0, bounces, true};
        }

        ProjectileProblem fp;
        fp.position0 = {x, y, 0.0};
        fp.velocity0 = {vx, vy, 0.0};
        fp.mass = 1.0;
        fp.gravity = {0.0, -p.gravity, 0.0};
        fp.linearDrag = k;
        const double flight = fp.landingTime();

        if (tTarget - tCur <= flight) {
            const double tau = tTarget - tCur;
            const Vec3 pos = fp.position(tau), vel = fp.velocity(tau);
            return {pos.x, pos.y, vel.x, vel.y, bounces, false};
        }
        if (!(flight > 0.0)) {  // vol de durée nulle (vitesse de rebond infime) : traité comme un repos
            if (restStart) *restStart = tCur;
            return {x, 0.0, vx * std::exp(-k * std::max(tTarget - tCur, 0.0)), 0.0, bounces, true};
        }

        // Impact : la composante verticale est inversée et réduite par e, l'horizontale est inchangée (sol lisse).
        const Vec3 pImpact = fp.position(flight), vImpact = fp.velocity(flight);
        tCur += flight;
        x = pImpact.x;
        y = 0.0;
        vx = vImpact.x;
        vy = -e * vImpact.y;
        ++bounces;
        if (impacts) {
            impacts->push_back(tCur);
            if (impacts->size() >= maxImpacts) return {x, 0.0, vx, vy, bounces, false};
        }
    }
    return {kNaN, kNaN, kNaN, kNaN, bounces, false};  // garde-fou : suite anormalement longue
}

}  // namespace

OdeFunction BounceProblem::rhsFlight() const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const double g = gravity, k = drag;
    return [g, k](double, const State& y, State& d) {
        d[0] = y[2];
        d[1] = y[3];
        d[2] = -k * y[2];
        d[3] = -g - k * y[3];
    };
}

OdeFunction BounceProblem::rhsRest() const {
    const double k = drag;
    return [k](double, const State& y, State& d) {
        d[0] = y[2];
        d[1] = 0.0;
        d[2] = -k * y[2];
        d[3] = 0.0;
    };
}

BounceState BounceProblem::exact(double t) const { return trace(*this, t, nullptr, 0, nullptr); }

std::vector<double> BounceProblem::impactTimes(int maxCount) const {
    std::vector<double> times;
    if (maxCount > 0) trace(*this, kInf, &times, static_cast<std::size_t>(maxCount), nullptr);
    return times;
}

double BounceProblem::restTime() const {
    if (restitution >= 1.0 && drag == 0.0) return kInf;  // la balle ne perd jamais d'énergie : elle rebondit indéfiniment
    double rest = kInf;
    trace(*this, kInf, nullptr, 0, &rest);
    return rest;
}

// ---------------------------------- pilote ------------------------------------------

BounceRun::BounceRun(const BounceProblem& problem, ContactModel model)
    : problem_(problem), model_(model), y_(problem.initialState()) {
    if (model == ContactModel::EventDriven && y_[1] <= 0.0 && y_[3] <= problem.restSpeed) {  // posée au sol dès le départ
        y_[1] = 0.0;
        y_[3] = 0.0;
        resting_ = true;
        restTime_ = 0.0;
    } else if (model == ContactModel::EventDriven && y_[1] <= 0.0) {
        launch_ = true;  // lancée depuis le sol vers le haut
    }
}

void BounceRun::advance(Solver& solver, double dt) {
    const OdeFunction flight = problem_.rhsFlight();
    if (model_ == ContactModel::Naive) {
        // Un pas, puis on regarde : si la balle est passée sous le sol en descendant, on la ramène à y = 0 et on inverse vy.
        // Le contact a eu lieu quelque part dans le pas : le temps passé « sous le sol » est perdu, et au repos la balle vibre sans fin.
        pl::advance(solver, flight, time_, y_, dt);
        time_ += dt;
        if (y_[1] < 0.0 && y_[3] < 0.0) {
            y_[1] = 0.0;
            y_[3] = -problem_.restitution * y_[3];
            ++bounces_;
        }
        return;
    }

    // Événement : l'impact est localisé exactement (g = max(y, 0) s'annule), la vitesse est inversée, puis on repart ou on se pose.
    const EventFunction ground = [](double, const State& y) { return std::max(y[1], 0.0); };
    const OdeFunction rest = problem_.rhsRest();
    double remaining = dt;
    for (int guard = 0; guard < 100000 && remaining > 1e-15 * dt; ++guard) {
        if (resting_) {
            pl::advance(solver, rest, time_, y_, remaining);
            time_ += remaining;
            return;
        }
        // Départ du sol (juste après un rebond) : g vaut exactement 0, et si le vol entier tient dans le pas, aucun changement de signe
        // n'est visible : la balle traverserait le sol (cela arrive dès que les rebonds deviennent plus courts que le pas, près de
        // l'accumulation de Zénon). On décolle donc d'abord par un pas minuscule (1/1000 du temps de montée), après quoi g > 0.
        if (launch_) {
            launch_ = false;
            const double kick = std::min(remaining, 1e-3 * std::max(y_[3], 0.0) / problem_.gravity);
            if (kick > 0.0) {
                pl::advance(solver, flight, time_, y_, kick);
                time_ += kick;
                remaining -= kick;
            }
            continue;
        }
        const EventStep r = advanceToEvent(solver, flight, ground, time_, y_, remaining);
        time_ += r.elapsed;
        remaining -= r.elapsed;
        if (r.event) {
            y_[1] = 0.0;
            y_[3] = -problem_.restitution * y_[3];
            ++bounces_;
            if (y_[3] <= problem_.restSpeed) {  // le rebond serait inférieur au seuil : la balle s'arrête de rebondir
                y_[3] = 0.0;
                resting_ = true;
                restTime_ = time_;
            } else {
                launch_ = true;
            }
        }
    }
}

double bounceError(const BounceProblem& problem, ContactModel model, Solver& solver, int steps, double tEnd) {
    BounceRun run(problem, model);
    const double dt = tEnd / steps;
    for (int i = 0; i < steps; ++i) run.advance(solver, dt);
    const BounceState ex = problem.exact(tEnd);
    return std::hypot(run.state()[0] - ex.x, run.state()[1] - ex.y);
}

}  // namespace pl
