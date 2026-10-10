// O1 : corde vibrante. Une corde fixée aux deux bouts, pincée ou lancée dans un seul mode, est vue comme une chaîne de N masses reliées par des
// ressorts (le ressort-masse de M2, N fois) et intégrée par les solveurs de la Mécanique ; elle est comparée à la solution exacte par modes
// propres. On y lit les harmoniques f_n = n c / (2 L), les harmoniques absents selon le point de pincement, la dispersion de la chaîne
// et, par FFT du mouvement d'une masse, les fréquences des modes.
#pragma once

#include <memory>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/waves/String.hpp"

namespace pl {

class StringModule final : public SimulationModule {
public:
    StringModule();

    const char* domain() const override { return "Ondes"; }
    const char* title() const override { return "O1 - Corde vibrante (modes propres et harmoniques)"; }
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
    static constexpr std::size_t kMaxRecord = 65536;  // échantillons gardés pour la FFT

    int probeBead() const;                 // masse observée (indice 0 .. N-1), vers x = 0.29 L : ne tombe sur aucun noeud des modes 1 à 6
    double maxAbsPositions() const;
    double stabilityFactorLimit() const;   // limite de dt / (a / c) du solveur choisi (0 : instable à tout pas)
    void sample();

    // paramètres (unités SI)
    double length_ = 1.0;                  // L [m]
    double tension_ = 4.0;                 // T [N]
    double density_ = 1.0;                 // mu [kg/m]
    int beads_ = 59;                       // N (N + 1 = 60 : divisible par 2, 3, 4, 5, 6 : les harmoniques absents sont exactement nuls)
    double pluckPos_ = 0.37;               // x0 / L
    double height_ = 0.1;                  // hauteur du pincement [m] (ou amplitude du mode pur)
    int startMode_ = 0;                    // 0 : pincement ; n >= 1 : la corde part dans le mode pur n
    int solverIndex_ = SolverSet::kRK4;
    double dtFactor_ = 0.2;                // dt = dtFactor a / c (limite de Verlet : 1) ; = kDefaultDtFactor du .cpp
    double timeScale_ = 1.0;
    bool running_ = true;

    // simulation
    waves::StringProblem problem_;
    std::unique_ptr<waves::StringModes> modes_;
    std::vector<double> u0_, v0_;
    State y_, exact_;                      // état calculé et solution exacte à l'instant courant
    std::unique_ptr<Solver> solver_;
    OdeFunction rhs_;
    double dt_ = 0.0, time_ = 0.0, accumulator_ = 0.0, lastSample_ = 0.0, endTime_ = 1.0, energy0_ = 1.0;
    bool finished_ = false, diverged_ = false;

    // courbes et enregistrement
    Series probe_, probeExact_, drift_;
    std::vector<double> record_;           // u de la masse observée, un échantillon par pas
    std::vector<double> spectrum_;         // cache du spectre (Hann, remplissage x4), rafraîchi tous les 8 affichages
    double binWidth_ = 0.0;
    std::vector<double> readFreq_;         // fréquences des modes 1 à 6 lues dans ce spectre (0 : pas de pic)
    int spectrumAge_ = 1000;

    void refreshSpectrum();                // met à jour spectrum_ et readFreq_ si l'enregistrement a avancé
};

}  // namespace pl
