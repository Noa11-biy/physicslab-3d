// M3 : pendule simple. Le même pendule est intégré par 5 solveurs, chacun dans son couloir 3D, et comparé à la
// solution exacte (fonctions elliptiques) quand elle existe, sinon à un RK45 haute précision.
#pragma once

#include <array>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Pendulum.hpp"

namespace pl {

class PendulumModule final : public SimulationModule {
public:
    PendulumModule();

    const char* title() const override { return "M3 - Pendule simple"; }
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
        State y;                           // [theta | omega]
        Series theta, phase, energy, error;  // phase : theta en abscisse, omega en ordonnée
        bool diverged = false;
    };
    struct ConvergenceCurve {
        std::vector<double> dt, err;
        double slope = 0.0;
    };

    void sample();
    void computeConvergence();
    void referenceState(double& theta, double& omega) const;

    // paramètres (unités SI)
    double length_ = 1.0;           // [m]
    double mass_ = 1.0;             // [kg]
    double gravity_ = constants::g0;
    double damping_ = 0.0;          // b [N.m.s]
    double thetaDeg_ = 60.0;        // angle initial [deg]
    double omegaInit_ = 0.0;        // vitesse angulaire initiale [rad/s]
    double durationPeriods_ = 10.0; // en périodes (de la solution exacte si elle existe)
    double dt_ = 1.0 / 30.0;        // pas de calcul [s] (grossier à dessein : les erreurs se voient)
    double relTol_ = 1e-8;          // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    PendulumProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    RK45 refSolver_;                // référence numérique quand il n'y a pas de solution exacte
    State refY_;
    Series thetaRef_, phaseRef_, energyRef_, harmonic_;  // référence et approximation des petits angles
    std::array<ConvergenceCurve, SolverSet::kFixedStep> convergence_;
    bool convergenceDirty_ = true;
    StepClock clock_;
    double endTime_ = 1.0;
    bool running_ = true;
};

}  // namespace pl
