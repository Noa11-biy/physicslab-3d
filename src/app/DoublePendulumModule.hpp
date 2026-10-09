// M3 (suite) : pendule double, système chaotique. Le même pendule est intégré par 5 solveurs dans le même plan, face
// à une référence RK45 haute précision et à un "jumeau" dont le départ diffère de 1e-9 rad : l'écart exponentiel entre
// la référence et son jumeau donne l'exposant de Lyapunov, et l'horizon de prédictibilité de chaque méthode.
#pragma once

#include <array>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/DoublePendulum.hpp"

namespace pl {

class DoublePendulumModule final : public SimulationModule {
public:
    DoublePendulumModule();

    const char* title() const override { return "M3b - Pendule double (chaos)"; }
    const char* explanation(Level level) const override;

    void reset() override;
    void update(double frameSeconds) override;

    void drawControls(const UiContext& ctx) override;
    void drawInvariants(const UiContext& ctx) override;
    void drawGraphs(const UiContext& ctx) override;
    bool hasAnalysis(const UiContext& ctx) const override { return atLeast(ctx.level, Level::Lycee); }
    void drawAnalysis(const UiContext& ctx) override;

    void drawScene(Renderer& renderer, const UiContext& ctx) override;
    void frameCamera(Camera& camera) const override;

private:
    struct Run {
        State y;                           // [theta1, theta2 | omega1, omega2]
        Series theta2, energy, error;      // energy : E(t) - E(0)
        std::vector<Vertex> trail;         // trajectoire de la seconde masse
        bool diverged = false;
        double horizon = -1.0;             // premier instant où l'erreur dépasse le seuil (< 0 : pas encore)
    };
    struct ConvergenceCurve {
        std::vector<double> dt, err;
        double slope = 0.0;
    };

    void sample();
    void computeConvergence();

    // paramètres (unités SI)
    double m1_ = 1.0, m2_ = 1.0;        // [kg]
    double l1_ = 1.0, l2_ = 1.0;        // [m]
    double gravity_ = constants::g0;
    double theta1Deg_ = 120.0, theta2Deg_ = -10.0;  // angles initiaux [deg]
    double omega1_ = 0.0, omega2_ = 0.0;            // vitesses angulaires initiales [rad/s]
    double perturbation_ = 1e-9;        // écart initial du jumeau sur theta1 [rad]
    double durationSec_ = 20.0;         // [s]
    double dt_ = 1.0 / 60.0;            // pas de calcul [s]
    double relTol_ = 1e-8;              // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    DoublePendulumProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    RK45 refSolver_, twinSolver_;       // référence et jumeau : RK45 à tolérance 1e-13
    State refY_, twinY_;
    Series theta2Ref_;
    std::vector<Vertex> refTrail_, twinTrail_;
    std::vector<double> twinT_, twinD_;  // distance référence - jumeau en fonction du temps (pour l'exposant de Lyapunov)
    Series configX_;                     // espace des configurations (theta1, theta2) de la référence
    double e0_ = 0.0;
    std::array<ConvergenceCurve, SolverSet::kFixedStep> convergence_;
    bool convergenceDirty_ = true;
    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
