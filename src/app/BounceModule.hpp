// M5b : chocs et rebonds. Deux scénarios : une balle qui rebondit sur le sol (accumulation de Zénon : une infinité de rebonds en temps
// fini) et le choc de deux billes (frontal ou oblique). Le même système est intégré par 5 solveurs avec l'un de deux modèles de
// détection du contact (naïf : vu après le pas ; événement : instant exact), face à la solution exacte.
#pragma once

#include <array>
#include <memory>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Bounce.hpp"
#include "physicslab/mechanics/Collision.hpp"

namespace pl {

class BounceModule final : public SimulationModule {
public:
    BounceModule();

    const char* title() const override { return "M5b - Chocs et rebonds"; }
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
    enum Scenario { kBall = 0, kBalls };

    struct Run {
        std::unique_ptr<BounceRun> ball;       // scénario de la balle
        std::unique_ptr<TwoBallRun> balls;     // scénario des deux billes
        Series primary, secondary, error, energy;
        std::vector<Vertex> trail, trail2;
    };
    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    void sample();
    void computeConvergence();
    std::string verdict() const;

    // paramètres
    Scenario scenario_ = kBall;
    double restitution_ = 0.8;       // balle
    double height_ = 2.0;            // hauteur de chute [m]
    double horizontal_ = 1.5;        // vitesse horizontale [m/s]
    double drag_ = 0.0;              // résistance k [1/s]
    double restSpeed_ = 1e-4;        // seuil d'arrêt [m/s]
    double restitutionBalls_ = 1.0;  // billes
    double impact_ = 0.6;            // décalage vertical de la bille 1 (paramètre d'impact) [m]
    double approach_ = 3.0;          // vitesse de la bille 1 [m/s]
    double massRatio_ = 1.0;         // m2 / m1
    double dt_ = 0.02;
    double relTol_ = 1e-8;
    double timeScale_ = 1.0;
    ContactModel model_ = ContactModel::Naive;

    // simulation
    BounceProblem ball_;
    TwoBallProblem balls_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    Series exactPrimary_, exactSecondary_, exactEnergy_;
    std::vector<double> impactX_, impactT_;       // impacts exacts de la balle (jusqu'à 40)
    double duration_ = 8.0;
    double scale_ = 1.0, offsetX_ = 0.0;
    std::array<Curve, SolverSet::kFixedStep> convergence_;
    double convergenceTime_ = 2.0;
    bool convergenceDirty_ = true;
    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
