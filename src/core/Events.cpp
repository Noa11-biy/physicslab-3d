#include "physicslab/core/Events.hpp"

#include <algorithm>

namespace pl {

EventStep advanceToEvent(Solver& solver, const OdeFunction& f, const EventFunction& g, double t, State& y, double dt) {
    const State start = y;
    const double g0 = g(t, start);
    const double elapsed = advance(solver, f, t, y, dt);
    const double g1 = g(t + elapsed, y);

    // Changement de signe de g sur le pas. Un départ exactement sur la surface (g0 = 0) n'est pas un événement.
    const bool crossed = (g0 < 0.0 && g1 >= 0.0) || (g0 > 0.0 && g1 <= 0.0);
    if (!crossed) return {elapsed, false};

    // Bissection sur la durée du pas : [lo, hi] encadre l'événement, g(hi) est de l'autre côté de g0 (ou nul).
    double lo = 0.0, hi = elapsed;
    State atHi = y;
    for (int i = 0; i < 100 && hi - lo > 8.9e-16 * hi; ++i) {
        const double mid = 0.5 * (lo + hi);
        State s = start;
        advance(solver, f, t, s, mid);
        const double gm = g(t + mid, s);
        if ((g0 < 0.0 && gm < 0.0) || (g0 > 0.0 && gm > 0.0)) {
            lo = mid;
        } else {
            hi = mid;
            atHi = s;
        }
    }
    y = atHi;
    return {hi, true};
}

}  // namespace pl
