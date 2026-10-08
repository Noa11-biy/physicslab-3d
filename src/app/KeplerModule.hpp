// M4a : gravitation, orbite de Kepler à deux corps (problème réduit). La même orbite est intégrée par 5 solveurs face à la
// solution exacte (équation de Kepler). En interne tout est normalisé (GM = 1, demi-grand axe = 1, période = 2 pi) : l'orbite
// ne dépend que de l'excentricité. Le demi-grand axe (UA) et la masse de l'étoile ne servent qu'à l'affichage en unités réelles.
#pragma once

#include <array>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Kepler.hpp"

namespace pl {

class KeplerModule final : public SimulationModule {
public:
    KeplerModule();

    const char* title() const override { return "M4a - Gravitation : orbite de Kepler"; }
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
        State y;                                         // [x, y, z | vx, vy, vz], unités normalisées
        Series radius, energy, momentum, stepSize, precession;  // r [UA], dE/|E0|, dL/L0, pas RK45, précession cumulée [deg]
        std::vector<Vertex> trail;
        PeriapsisTracker tracker;
        bool diverged = false;
        double lastStep = 0.0;                           // dernier pas accepté (RK45)
    };
    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    void sample();
    void advanceRun(int i, Run& run, double t, double h);
    void buildGeometry();
    void computeConvergence();
    void computePrecessionScan();

    // conversions normalisé -> réel
    double lengthUnit() const;   // [m]
    double timeUnit() const;     // [s]

    // paramètres
    double ecc_ = 0.5;             // excentricité
    double semiMajorAU_ = 1.0;     // demi-grand axe [UA] (affichage)
    double massSun_ = 1.0;         // masse de l'étoile [masses solaires] (affichage)
    double stepsPerOrbit_ = 200.0; // pas de calcul par période
    double orbits_ = 10.0;         // durée simulée [orbites]
    double relTol_ = 1e-8;         // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    KeplerProblem problem_;
    OdeFunction rhs_;
    SolverSet solvers_;
    std::array<Run, SolverSet::kCount> runs_;
    Series exactRadius_, exactKinetic_, exactPotential_, exactTotal_;
    std::vector<Vertex> ellipse_, sectors_;              // orbite exacte et rayons à intervalles de temps égaux
    double e0_ = -0.5;
    std::array<Curve, SolverSet::kFixedStep> convergence_;
    Curve scanSymplectic_, scanVerlet_, scanTheory_;      // précession par orbite selon le pas
    bool convergenceDirty_ = true, scanDirty_ = true;
    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
