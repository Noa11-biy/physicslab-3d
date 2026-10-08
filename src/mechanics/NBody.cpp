#include "physicslab/mechanics/NBody.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "physicslab/mechanics/Kepler.hpp"

namespace pl {

void nbody::accelerations(const double* positions, const double* masses, int n, double G, double softening, double* acc) {
    std::fill(acc, acc + 3 * n, 0.0);
    const double eps2 = softening * softening;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = positions[3 * j] - positions[3 * i];
            const double dy = positions[3 * j + 1] - positions[3 * i + 1];
            const double dz = positions[3 * j + 2] - positions[3 * i + 2];
            const double r2 = dx * dx + dy * dy + dz * dz + eps2;
            const double inv = 1.0 / (r2 * std::sqrt(r2));  // 1 / (r^2 + eps^2)^(3/2)
            const double gi = G * masses[j] * inv, gj = G * masses[i] * inv;
            // La même force est appliquée aux deux corps, de signes opposés : impulsion conservée à l'arrondi près.
            acc[3 * i] += gi * dx;
            acc[3 * i + 1] += gi * dy;
            acc[3 * i + 2] += gi * dz;
            acc[3 * j] -= gj * dx;
            acc[3 * j + 1] -= gj * dy;
            acc[3 * j + 2] -= gj * dz;
        }
    }
}

double nbody::potentialEnergy(const double* positions, const double* masses, int n, double G, double softening) {
    const double eps2 = softening * softening;
    double u = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = positions[3 * j] - positions[3 * i];
            const double dy = positions[3 * j + 1] - positions[3 * i + 1];
            const double dz = positions[3 * j + 2] - positions[3 * i + 2];
            u -= G * masses[i] * masses[j] / std::sqrt(dx * dx + dy * dy + dz * dz + eps2);
        }
    }
    return u;
}

State NBodyProblem::initialState() const {
    const int n = count();
    State y(6 * n);
    for (int i = 0; i < n; ++i) {
        y[3 * i] = position[i].x;
        y[3 * i + 1] = position[i].y;
        y[3 * i + 2] = position[i].z;
        y[3 * n + 3 * i] = velocity[i].x;
        y[3 * n + 3 * i + 1] = velocity[i].y;
        y[3 * n + 3 * i + 2] = velocity[i].z;
    }
    return y;
}

OdeFunction NBodyProblem::rhs() const {
    // Copie des paramètres : la fonction reste valide même si le problème est modifié ensuite.
    const int n = count();
    return [masses = mass, G = G, eps = softening, n](double, const State& y, State& d) {
        std::copy(y.begin() + 3 * n, y.begin() + 6 * n, d.begin());          // dx/dt = v
        nbody::accelerations(y.data(), masses.data(), n, G, eps, d.data() + 3 * n);  // dv/dt = a
    };
}

double NBodyProblem::kineticEnergy(const State& y) const {
    const int n = count();
    double t = 0.0;
    for (int i = 0; i < n; ++i) {
        const double* v = &y[3 * n + 3 * i];
        t += 0.5 * mass[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    return t;
}

double NBodyProblem::potentialEnergy(const State& y) const {
    return nbody::potentialEnergy(y.data(), mass.data(), count(), G, softening);
}

Vec3 NBodyProblem::momentum(const State& y) const {
    const int n = count();
    Vec3 p;
    for (int i = 0; i < n; ++i) p += mass[i] * Vec3{y[3 * n + 3 * i], y[3 * n + 3 * i + 1], y[3 * n + 3 * i + 2]};
    return p;
}

Vec3 NBodyProblem::angularMomentum(const State& y) const {
    const int n = count();
    Vec3 l;
    for (int i = 0; i < n; ++i)
        l += mass[i] * cross({y[3 * i], y[3 * i + 1], y[3 * i + 2]}, {y[3 * n + 3 * i], y[3 * n + 3 * i + 1], y[3 * n + 3 * i + 2]});
    return l;
}

double NBodyProblem::totalMass() const {
    double m = 0.0;
    for (double mi : mass) m += mi;
    return m;
}

Vec3 NBodyProblem::centerOfMass(const State& y) const {
    const int n = count();
    Vec3 c;
    for (int i = 0; i < n; ++i) c += mass[i] * Vec3{y[3 * i], y[3 * i + 1], y[3 * i + 2]};
    return c / totalMass();
}

State NBodyProblem::reference(double t) const {
    RK45 solver;
    solver.relTol = 1e-13;
    solver.absTol = 1e-15;
    State y = initialState();
    advance(solver, rhs(), 0.0, y, t);
    return y;
}

double NBodyProblem::distance(const State& a, const State& b) const {
    // Échelles de l'état initial : longueur L0 = distance moyenne au centre de masse, vitesse sqrt(G M / L0).
    const int n = count();
    const State y0 = initialState();
    const Vec3 c0 = centerOfMass(y0);
    double length = 0.0;
    for (int i = 0; i < n; ++i) length += (Vec3{y0[3 * i], y0[3 * i + 1], y0[3 * i + 2]} - c0).norm();
    length = n > 0 && length > 0.0 ? length / n : 1.0;
    const double vScale = std::sqrt(G * totalMass() / length);

    double sum = 0.0;
    for (int k = 0; k < 3 * n; ++k) {
        const double dp = (a[k] - b[k]) / length, dv = (a[3 * n + k] - b[3 * n + k]) / vScale;
        sum += dp * dp + dv * dv;
    }
    return std::sqrt(sum);
}

NBodyProblem NBodyProblem::twoBody(double m1, double m2, double semiMajor, double eccentricity, double G) {
    KeplerProblem k;
    k.mu = G * (m1 + m2);
    k.a = semiMajor;
    k.e = eccentricity;
    const State rel = k.initialState();  // position et vitesse relatives r2 - r1, v2 - v1
    const Vec3 r{rel[0], rel[1], rel[2]}, v{rel[3], rel[4], rel[5]};
    const double total = m1 + m2;

    NBodyProblem p;
    p.G = G;
    p.mass = {m1, m2};
    p.position = {-(m2 / total) * r, (m1 / total) * r};  // centre de masse à l'origine
    p.velocity = {-(m2 / total) * v, (m1 / total) * v};
    return p;
}

NBodyProblem NBodyProblem::lagrangeTriangle(double side, double massEach, double G) {
    const double radius = side / std::sqrt(3.0);
    const double omega = std::sqrt(3.0 * G * massEach / (side * side * side));
    NBodyProblem p;
    p.G = G;
    p.mass = {massEach, massEach, massEach};
    for (int i = 0; i < 3; ++i) {
        const double phi = constants::pi * (0.5 + 2.0 * i / 3.0);  // 90°, 210°, 330°
        p.position.push_back({radius * std::cos(phi), radius * std::sin(phi), 0.0});
        p.velocity.push_back({-omega * radius * std::sin(phi), omega * radius * std::cos(phi), 0.0});
    }
    return p;
}

NBodyProblem NBodyProblem::figureEight() {
    // Chenciner-Montgomery, conditions initiales numériques de C. Simó : x1 = -x2, x3 = 0, v1 = v2 = -v3/2.
    NBodyProblem p;
    p.mass = {1.0, 1.0, 1.0};
    p.position = {{0.97000436, -0.24308753, 0.0}, {-0.97000436, 0.24308753, 0.0}, {0.0, 0.0, 0.0}};
    p.velocity = {{0.466203685, 0.43236573, 0.0}, {0.466203685, 0.43236573, 0.0}, {-0.93240737, -0.86473146, 0.0}};
    return p;
}

NBodyProblem NBodyProblem::randomCluster(int n, unsigned seed, double radius, double softening) {
    std::mt19937 rng(seed);
    auto uniform = [&] { return (rng() + 0.5) / 4294967296.0; };  // dans ]0, 1[
    auto gaussian = [&] { return std::sqrt(-2.0 * std::log(uniform())) * std::cos(2.0 * constants::pi * uniform()); };

    NBodyProblem p;
    p.softening = softening;
    p.mass.assign(n, 1.0 / n);
    for (int i = 0; i < n; ++i) {
        Vec3 r;
        do {
            r = {2.0 * uniform() - 1.0, 2.0 * uniform() - 1.0, 2.0 * uniform() - 1.0};
        } while (r.norm2() > 1.0);  // tirage uniforme dans la boule unité
        p.position.push_back(radius * r);
        p.velocity.push_back({gaussian(), gaussian(), gaussian()});
    }

    // Centre de masse au repos à l'origine (masses égales : moyenne simple).
    Vec3 meanPos, meanVel;
    for (int i = 0; i < n; ++i) { meanPos += p.position[i]; meanVel += p.velocity[i]; }
    meanPos /= n;
    meanVel /= n;
    for (int i = 0; i < n; ++i) { p.position[i] -= meanPos; p.velocity[i] -= meanVel; }

    // Échelle des vitesses : équilibre du viriel, 2T = -U.
    State y = p.initialState();
    const double kinetic = p.kineticEnergy(y), potential = p.potentialEnergy(y);
    if (kinetic > 0.0 && potential < 0.0) {
        const double scale = std::sqrt(-potential / (2.0 * kinetic));
        for (Vec3& v : p.velocity) v *= scale;
    }
    return p;
}

double nbodyError(const NBodyProblem& problem, Solver& solver, int steps, double tEnd) {
    const OdeFunction f = problem.rhs();
    State y = problem.initialState();
    const double dt = tEnd / steps;
    double t = 0.0;
    for (int i = 0; i < steps; ++i) t += advance(solver, f, t, y, dt);
    return problem.distance(y, problem.reference(t));
}

}  // namespace pl
