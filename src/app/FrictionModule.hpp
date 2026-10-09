// M5a : frottement sec de Coulomb sur un plan incliné. Le même bloc est intégré par 5 solveurs avec l'un des trois modèles
// (naïf, régularisé, événement + adhérence) face à la solution exacte par morceaux. Ce module montre pourquoi une force
// DISCONTINUE fait tomber tous les schémas à l'ordre 1, pourquoi le bloc « naïf » ne s'arrête jamais, et comment un événement
// (arrêt exact à v = 0) puis un état d'adhérence règlent le problème.
#pragma once

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Friction.hpp"

namespace pl {

class FrictionModule final : public SimulationModule {
public:
    FrictionModule();

    const char* title() const override { return "M5a - Frottement sec (Coulomb) sur plan incliné"; }
    const char* explanation(Level level) const override;

    void reset() override;
    void update(double frameSeconds) override;

    void drawControls(const UiContext& ctx) override;
    void drawInvariants(const UiContext& ctx) override;
    void drawGraphs(const UiContext& ctx) override;
    bool hasAnalysis(const UiContext& ctx) const override { return atLeast(ctx.level, Level::Etudiant); }
    void drawAnalysis(const UiContext& ctx) override;

    void drawScene(Renderer& renderer, const UiContext& ctx) override;
    void frameCamera(Camera& camera) const override;

private:
    struct Run {
        State y;                                    // [s, v]
        std::unique_ptr<InclineRun> event;          // seulement pour le modèle « événement + adhérence »
        Series position, velocity, error, energy;
        bool stalled = false;                       // solveur adaptatif bloqué (pas devenu trop petit)
        double stallTime = 0.0;
        double stopTime = std::numeric_limits<double>::quiet_NaN();  // premier arrêt détecté
    };
    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    void sample();
    void applyPreset(double angleDeg, double muS, double muK, double v0);
    void computeConvergence();
    std::string verdict() const;
    Vertex onSlope(double s, double normal, float z, const float* color, float dim = 1.0f) const;

    // paramètres
    double angleDeg_ = 26.0;     // pente [deg]
    double muStatic_ = 0.5;      // adhérence mu_s
    double muKinetic_ = 0.4;     // glissement mu_d
    double launch_ = 6.0;        // vitesse de lancement vers le haut [m/s]
    double drag_ = 0.0;          // résistance k [1/s]
    double epsilon_ = 0.05;      // régularisation [m/s]
    double durationSec_ = 5.0;
    double dt_ = 0.01;
    double relTol_ = 1e-8;
    double timeScale_ = 1.0;
    InclineModel model_ = InclineModel::Naive;

    // simulation
    InclineProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    Series exactPosition_, exactVelocity_, exactEnergy_, exactFriction_;
    double lo_ = -10.0, hi_ = 5.0;       // extrémités de la rampe le long de la pente [m]
    float scale_ = 1.0f;                 // unités de scène par mètre
    std::array<Curve, SolverSet::kFixedStep> convergence_;
    bool convergenceDirty_ = true;
    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
