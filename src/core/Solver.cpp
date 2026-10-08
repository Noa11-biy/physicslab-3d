#include "physicslab/core/Solver.hpp"

namespace pl {

double ExplicitEuler::step(const OdeFunction& f, double t, State& y, double dt) {
    k_.resize(y.size());
    f(t, y, k_);
    for (std::size_t i = 0; i < y.size(); ++i) y[i] += dt * k_[i];
    return dt;
}

}  // namespace pl
