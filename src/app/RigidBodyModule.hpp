// M6 : corps rigide. Trois scénarios : un solide symétrique libre (cône de précession, formule exacte), un solide asymétrique libre (la
// « raquette de tennis » : l'axe du milieu est instable), la toupie pesante de Lagrange (précession, nutation, toupie endormie). Le même
// mouvement est calculé par cinq intégrateurs d'orientation (Euler, RK4 avec et sans renormalisation, groupe de Lie, découpage symplectique)
// face à une référence (formule exacte ou RK45 serré) : ce qu'on voit, c'est ce que chaque schéma fait de l'énergie, du moment cinétique et de |q|.
#pragma once

#include <array>
#include <memory>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/RigidBody.hpp"

namespace pl {

class RigidBodyModule final : public SimulationModule {
public:
    RigidBodyModule();

    const char* title() const override { return "M6 - Corps rigide"; }
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
    enum Scenario { kSymmetric = 0, kAsymmetric, kTop };
    enum TopMode { kSteadySlow = 0, kSteadyFast, kNutation, kSleeping };
    enum Method { kEuler = 0, kRK4, kRK4Raw, kLie, kSplitting, kMethodCount };
    static constexpr int kFixedMethods = 4;   // pour l'étude de convergence : Euler, RK4, Lie, découpage (RK4 brut = RK4 au pas près)

    struct Run {
        std::unique_ptr<RotationIntegrator> integrator;
        RotationState state;
        bool diverged = false;
        Series error, energy, momentum, norm;
        std::vector<Vertex> trail;
    };
    struct Curve {
        std::vector<double> x, y;
        double slope = 0.0;
    };

    // intégrateurs et affichage
    bool isShown(Level level, int method) const;
    std::string label(Level level, int method) const;
    std::string shortLabel(int method) const;
    std::string plotLabel(Level level, int method) const { return atLeast(level, Level::Etudiant) ? shortLabel(method) : label(level, method); }
    const float* color(int method) const;

    // physique commune aux trois scénarios
    double energyOf(const RotationState& s) const;
    double momentumDrift(const RotationState& s) const;      // dérive relative du moment cinétique conservé (L fixe, ou L_z pour la toupie)
    Vec3 axisTip(int scenario) const;                         // point du corps dont on trace la trajectoire
    RotationState referenceState() const;                     // formule exacte (symétrique) ou RK45 serré en pas de temps commun

    void step(double h);
    void sample();
    void computeConvergence();
    std::string verdict() const;

    // paramètres
    Scenario scenario_ = kSymmetric;
    double shapeRatio_ = 0.5;          // I3 / I1 du solide symétrique (0,2 : cigare, 2 : disque)
    double spin_ = 3.0;                // composante w3 du solide symétrique
    double wobble_ = 0.8;              // composante transverse (écart de l'axe de rotation)
    int axis_ = 1;                     // axe de rotation du solide asymétrique : 0 long, 1 du milieu, 2 court
    double perturb_ = 0.05;            // petite perturbation
    TopMode topMode_ = kSteadySlow;
    double topSpin_ = 25.0;            // w3 de la toupie (ou, en mode « endormie », rapporté au spin critique : voir sleepFactor_)
    double sleepFactor_ = 0.7;         // spin / spin critique (toupie endormie)
    double topTilt_ = 0.6;             // inclinaison initiale theta
    double dt_ = 0.01;
    double timeScale_ = 1.0;
    std::array<bool, kMethodCount> show_{true, true, false, true, true};
    bool steadyFallback_ = false;      // précession régulière impossible (spin trop faible) : la toupie part en nutation

    // problème courant
    FreeBodyProblem free_;
    HeavyTopProblem top_;
    Vec3 inertia_{1.0, 1.0, 1.0};
    TorqueFunction torque_;
    OdeFunction rhs_;
    RotationState initial_;
    State yRef_;
    RK45 ref_;
    double duration_ = 20.0;
    double energy0_ = 1.0, momentum0_ = 1.0;
    Vec3 lSpace0_;                     // moment cinétique initial (repère fixe), pour les flèches
    double omegaScale_ = 1.0;          // |w| de départ (échelle des flèches)
    SteadyPrecession steady_;          // précession régulière de la toupie, si elle existe
    std::vector<Vec3> bodyLines_;      // segments du corps (repère du corps), par paires
    std::array<std::vector<Vec3>, 3> axisLines_;   // les trois axes du solide asymétrique
    double convergenceTime_ = 6.3;

    // courbes
    std::array<Run, kMethodCount> runs_;
    std::array<Series, 3> refW_;       // w1, w2, w3 de la référence (solides libres)
    Series refTheta_;                  // inclinaison de la toupie (degrés)
    Series refCos_;
    std::vector<Vertex> refTrail_;

    std::array<Curve, kFixedMethods> convergence_;
    bool convergenceDirty_ = true;

    StepClock clock_;
    bool running_ = true;
};

}  // namespace pl
