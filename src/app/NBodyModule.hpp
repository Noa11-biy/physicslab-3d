// M4b : gravitation à N corps. Le même système est intégré par 5 solveurs face à une référence RK45 haute précision et à un
// « jumeau » dont le départ diffère de 1e-9 : l'écart référence-jumeau dit si le système est chaotique (amas : croissance
// exponentielle) ou stable (huit de Chenciner-Montgomery : croissance linéaire). Unités normalisées : G = 1.
#pragma once

#include <array>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/NBody.hpp"

namespace pl {

class NBodyModule final : public SimulationModule {
public:
    NBodyModule();

    const char* title() const override { return "M4b - Gravitation : problème à N corps"; }
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
    enum Scenario { kEight = 0, kLagrange, kCluster, kScenarioCount };

    struct Run {
        State y;
        Series energy, error, angular;             // dE/|E0|, distance à la référence, |L - L0|
        std::vector<std::vector<Vertex>> trails;   // une trace par corps
        bool diverged = false;
        double horizon = -1.0;                     // premier instant où l'écart dépasse le seuil (< 0 : pas encore)
    };
    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    void selectScenario(Scenario s);
    void buildProblem();
    void sample();
    void computeConvergence();

    // paramètres
    Scenario scenario_ = kEight;
    double clusterCount_ = 6.0;      // nombre de corps de l'amas
    unsigned seed_ = 42;             // tirage de l'amas
    double mass3Factor_ = 1.0;       // masse du 3e corps du huit, en unités de la masse des deux autres
    double softening_ = 0.05;        // adoucissement de l'amas
    double durationSec_ = 19.0;      // durée simulée (unités normalisées)
    double dt_ = 0.005;              // pas de calcul
    double relTol_ = 1e-8;           // tolérance relative de RK45
    double perturbation_ = 1e-9;     // écart initial du jumeau sur la position du corps 1
    double timeScale_ = 1.0;

    // simulation
    NBodyProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    RK45 refSolver_, twinSolver_;
    State refY_, twinY_;
    std::vector<std::vector<Vertex>> refTrails_;
    std::vector<Series> bodyX_, pairDistance_;   // x(t) de chaque corps ; distances des premières paires
    Series minDistance_;                         // plus petite distance entre deux corps
    std::vector<double> twinT_, twinD_;
    double e0_ = -1.0;
    Vec3 l0_;
    double convergenceTime_ = 2.0;
    std::array<Curve, SolverSet::kFixedStep> convergence_;
    bool convergenceDirty_ = true;
    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
