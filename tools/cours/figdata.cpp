// Données des figures du cours : exécute les vrais solveurs du projet (physicslab_core) et écrit des CSV.
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "physicslab/core/Solver.hpp"
#include "physicslab/mechanics/DoublePendulum.hpp"
#include "physicslab/mechanics/Oscillator.hpp"

using namespace pl;

static std::unique_ptr<Solver> make(int i) {
    switch (i) {
        case 0: return std::make_unique<ExplicitEuler>();
        case 1: return std::make_unique<SymplecticEuler>();
        case 2: return std::make_unique<VelocityVerlet>();
        default: return std::make_unique<RK4>();
    }
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : ".";
    const double pi = std::acos(-1.0);

    // 1. Convergence de l'oscillateur (erreur dans l'espace des phases à t = 2,7 T)
    {
        OscillatorProblem p;  // m = 1, k = 10, x0 = 1
        const double tEnd = 2.7 * p.period();
        std::FILE* f = std::fopen((dir + "/convergence.csv").c_str(), "w");
        std::fprintf(f, "steps,euler,symplectic,verlet,rk4\n");
        for (int steps = 50; steps <= 6400; steps *= 2) {
            std::fprintf(f, "%d", steps);
            for (int i = 0; i < 4; ++i) {
                auto s = make(i);
                std::fprintf(f, ",%.6e", oscillatorError(p, *s, steps, tEnd));
            }
            std::fprintf(f, "\n");
        }
        std::fclose(f);
    }

    // 2. Énergie de l'oscillateur au fil de 30 périodes, 20 pas par période
    {
        OscillatorProblem p;
        const double T = p.period();
        const int perPeriod = 20, periods = 30;
        const double dt = T / perPeriod;
        const OdeFunction rhs = p.rhs();
        std::vector<State> y(4, p.initialState());
        std::vector<std::unique_ptr<Solver>> s;
        for (int i = 0; i < 4; ++i) s.push_back(make(i));
        const double e0 = p.energy(y[0][0], y[0][1]);
        std::FILE* f = std::fopen((dir + "/energy.csv").c_str(), "w");
        std::fprintf(f, "t,euler,symplectic,verlet,rk4\n");
        double t = 0.0;
        std::fprintf(f, "0,1,1,1,1\n");
        for (int n = 1; n <= perPeriod * periods; ++n) {
            for (int i = 0; i < 4; ++i) s[i]->step(rhs, t, y[i], dt);
            t += dt;
            std::fprintf(f, "%.6f", t / T);
            for (int i = 0; i < 4; ++i) std::fprintf(f, ",%.8e", p.energy(y[i][0], y[i][1]) / e0);
            std::fprintf(f, "\n");
        }
        std::fclose(f);
    }

    // 3. Orbites de Kepler (e = 0,5, GM = a = 1), 200 pas par orbite, 10 orbites
    {
        const OdeFunction kepler = [](double, const State& y, State& d) {
            const double r2 = y[0] * y[0] + y[1] * y[1], r3 = r2 * std::sqrt(r2);
            d[0] = y[2];
            d[1] = y[3];
            d[2] = -y[0] / r3;
            d[3] = -y[1] / r3;
        };
        const double e = 0.5, T = 2.0 * pi;
        const State y0 = {1.0 - e, 0.0, 0.0, std::sqrt((1.0 + e) / (1.0 - e))};
        const char* names[3] = {"euler", "verlet", "rk4"};
        const int idx[3] = {0, 2, 3};
        const int stepsPerOrbit = 200, orbits = 10;
        const double dt = T / stepsPerOrbit;
        for (int k = 0; k < 3; ++k) {
            auto s = make(idx[k]);
            State y = y0;
            std::FILE* f = std::fopen((dir + "/kepler_" + names[k] + ".csv").c_str(), "w");
            std::fprintf(f, "x,y\n%.8f,%.8f\n", y[0], y[1]);
            double t = 0.0;
            for (int n = 0; n < stepsPerOrbit * orbits; ++n) {
                s->step(kepler, t, y, dt);
                t += dt;
                std::fprintf(f, "%.8f,%.8f\n", y[0], y[1]);
            }
            std::fclose(f);
            const double energy = 0.5 * (y[2] * y[2] + y[3] * y[3]) - 1.0 / std::sqrt(y[0] * y[0] + y[1] * y[1]);
            std::printf("kepler %s : E finale = %.6f (exacte %.6f)\n", names[k], energy, -0.5);
        }
    }

    // 4. Pendule double : écart entre deux trajectoires voisines (1e-9 sur theta1), RK45 serré
    {
        DoublePendulumProblem a;
        a.m1 = a.m2 = a.l1 = a.l2 = 1.0;
        a.theta1 = 120.0 * pi / 180.0;
        a.theta2 = -10.0 * pi / 180.0;
        DoublePendulumProblem b = a;
        b.theta1 += 1e-9;
        RK45 sa, sb;
        sa.relTol = sb.relTol = 1e-13;
        sa.absTol = sb.absTol = 1e-15;
        State ya = a.initialState(), yb = b.initialState();
        const OdeFunction fa = a.rhs(), fb = b.rhs();
        std::FILE* f = std::fopen((dir + "/double_pendulum.csv").c_str(), "w");
        std::fprintf(f, "t,distance\n");
        double t = 0.0;
        const double dt = 0.05;
        std::fprintf(f, "0,%.6e\n", a.distance(ya, yb));
        for (int n = 1; n <= 400; ++n) {
            advance(sa, fa, t, ya, dt);
            advance(sb, fb, t, yb, dt);
            t += dt;
            std::fprintf(f, "%.4f,%.6e\n", t, a.distance(ya, yb));
        }
        std::fclose(f);
    }
    return 0;
}
