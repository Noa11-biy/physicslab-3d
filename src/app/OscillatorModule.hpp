// M2 : ressort-masse (oscillateur harmonique) libre, amorti ou forcé. Le même oscillateur est intégré
// par 5 solveurs, chacun dans son couloir 3D, et comparé à la solution exacte (régime transitoire + permanent).
#pragma once

#include <array>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Oscillator.hpp"

namespace pl {

class OscillatorModule final : public SimulationModule {
public:
    OscillatorModule();

    const char* title() const override { return "M2 - Ressort-masse (oscillateur harmonique)"; }
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
    // L'état [x | v] d'un solveur et ses courbes.
    struct Run {
        State y;
        Series x, phase, energy, error;  // phase : x en abscisse, v en ordonnée
        bool diverged = false;           // arrêté si |x| devient absurde (Euler explicite est instable)
    };
    struct ConvergenceCurve {
        std::vector<double> dt, err;
        double slope = 0.0;
    };

    void sample();
    void computeConvergence();

    // paramètres (unités SI)
    double mass_ = 1.0;                  // [kg]
    double stiffness_ = 10.0;            // [N/m]
    double damping_ = 0.0;               // c [kg/s]
    double x0_ = 1.0;                    // [m]
    double v0_ = 0.0;                    // [m/s]
    bool forced_ = false;
    double forceAmp_ = 3.0;              // F0 [N]
    double forceOmega_ = 3.1622776601683795;  // [rad/s], par défaut la pulsation propre de (1 kg, 10 N/m)
    double durationPeriods_ = 10.0;      // durée de la simulation en périodes propres
    double dt_ = 1.0 / 30.0;             // pas de calcul [s] (grossier à dessein : les erreurs se voient)
    double relTol_ = 1e-8;               // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    OscillatorProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    Series xExact_, phaseExact_, energyExact_;
    std::array<ConvergenceCurve, SolverSet::kFixedStep> convergence_;
    bool convergenceDirty_ = true;       // recalculée seulement quand le panneau Analyse la demande
    double time_ = 0.0;
    double endTime_ = 1.0;
    double accumulator_ = 0.0;
    double lastSampleTime_ = 0.0;
    bool running_ = true;
    bool finished_ = false;
};

}  // namespace pl
