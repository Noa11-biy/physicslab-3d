// M5c : berceau de Newton. Une rangée de billes identiques sur un rail : les premières sont lancées sur les autres. Deux modèles de choc
// tournent côte à côte (au premier plan le contact de Hertz, qui résout la dynamique ; derrière les impulsions instantanées binaires) et
// donnent des résultats différents, ce que ni la conservation de l'impulsion ni celle de l'énergie ne permettent de trancher.
#pragma once

#include <array>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/Cradle.hpp"

namespace pl {

class CradleModule final : public SimulationModule {
public:
    CradleModule();

    const char* title() const override { return "M5c - Berceau de Newton"; }
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
    static constexpr int kMaxBalls = 7;
    enum Display { kBoth = 0, kHertz, kImpulses };

    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    void step(double h);
    void sample();
    bool slowMotion() const;
    void computeConvergence();

    double contactForce(const State& y, int pair) const;            // force de Hertz (avec amortissement éventuel)
    double phaseDistance(const State& a, const State& b) const;     // même métrique que cradleError
    double impulseVelocity(int ball, double t) const;               // modèle des impulsions : exact, par morceaux
    double impulsePosition(int ball, double t) const;

    // paramètres
    int balls_ = 5;
    int launched_ = 1;
    double restitution_ = 1.0;     // e des chocs instantanés
    double stiffness_ = 1.0e4;     // k de Hertz
    double damping_ = 0.0;         // alpha de Hunt-Crossley
    double speed_ = 1.0;           // vitesse de lancement
    double gap_ = 1.0;             // distance avant le contact
    double stepsPerContact_ = 100.0;
    double relTol_ = 1e-8;
    double timeScale_ = 1.0;
    int solverIndex_ = SolverSet::kRK4;
    Display display_ = kBoth;

    // simulation
    CradleProblem prob_;
    OdeFunction rhs_;
    State y_, yRef_;
    RK45 ref_;                     // référence en pas de temps commun (tolérance 1e-12)
    SolverSet solvers_;
    std::vector<double> v0_;       // vitesses initiales (impulsions)
    cradle::ImpulseResult impulses_;
    CradleOutcome outcome_;        // issue de référence du modèle de Hertz (RK45 serré)
    double contactTime_ = 0.0;     // durée d'un contact isolé
    double dt_ = 0.0;
    double duration_ = 4.0;
    double viewMin_ = -2.0, viewMax_ = 8.0;
    double forceScale_ = 1.0;      // force maximale d'un choc de deux billes (pour les barres)
    double energy0_ = 1.0;         // énergie initiale (normalisation des courbes)

    // mesures
    double maxCompression_ = 0.0;
    double contactStart_ = -1.0, contactEnd_ = -1.0;   // premier contact de la bille de tête du groupe lancé

    // courbes
    std::array<Series, kMaxBalls> vel_, impVel_;
    std::array<Series, kMaxBalls - 1> force_;
    Series kinetic_, potential_, total_, error_;

    std::array<Curve, SolverSet::kFixedStep> convergence_;
    double convergenceTime_ = 0.2;
    bool convergenceDirty_ = true;

    StepClock clock_;
    bool running_ = true;
    bool diverged_ = false;       // le schéma choisi a explosé (pas trop grand) : la simulation s'arrête
};

}  // namespace pl
