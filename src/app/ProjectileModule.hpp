// M1 : projectile avec frottement. Un même lancer est intégré en parallèle par 5 solveurs
// et comparé en direct à la solution analytique.
#pragma once

#include <array>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/core/Constants.hpp"
#include "physicslab/core/World.hpp"
#include "physicslab/mechanics/Projectile.hpp"

namespace pl {

class ProjectileModule final : public SimulationModule {
public:
    ProjectileModule();

    const char* title() const override { return "M1 - Projectile avec frottement"; }
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
    // Un solveur et son monde : positions, courbes et trajectoire 3D.
    struct Run {
        World world;
        Series y, posError, energyError;
        std::vector<Vertex> trail;
    };
    // Erreur finale en fonction du pas, pour un solveur à pas fixe.
    struct ConvergenceCurve {
        std::vector<double> dt, err;
        double slope = 0.0;  // ordre mesuré (pente en log-log)
    };

    void sample();
    void computeConvergence();

    // paramètres du lancer (unités SI) et du calcul
    double speed_ = 14.0;               // [m/s]
    double angleDeg_ = 55.0;            // [deg]
    double height_ = 2.0;               // [m]
    double gravity_ = constants::g0;    // [m/s^2]
    double mass_ = 1.0;                 // [kg]
    double drag_ = 0.0;                 // b [kg/s]
    double dt_ = 1.0 / 30.0;            // pas de calcul [s] (grossier à dessein : les erreurs se voient)
    double relTol_ = 1e-8;              // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    ProjectileProblem problem_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    Series refY_;                       // hauteur exacte
    std::vector<Vertex> refPath_;       // trajectoire exacte complète
    std::array<ConvergenceCurve, SolverSet::kFixedStep> convergence_;
    double convergenceEnd_ = 1.0;       // instant final de l'étude de convergence
    bool convergenceForcedDrag_ = false;  // l'étude impose un frottement (b nul => solution polynomiale)
    double landingTime_ = 0.0;
    double accumulator_ = 0.0;
    double lastSampleTime_ = 0.0;
    bool running_ = true;
    bool landed_ = false;
};

}  // namespace pl
