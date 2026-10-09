// O0 : une impulsion sur une grille. Une bosse (gaussienne) lâchée au repos sur une corde (1D) ou une surface d'eau (2D) se
// sépare en deux ondes qui s'éloignent ; aux bords elles disparaissent (éponge), rebondissent sur un mur ou sur un bout libre.
// Le calcul (CPU `double`, saute-mouton) est comparé à la solution exacte de d'Alembert (méthode des images) en 1D.
#pragma once

#include <memory>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/waves/Wave.hpp"

namespace pl {

class WaveModule final : public SimulationModule {
public:
    WaveModule();

    const char* domain() const override { return "Ondes"; }
    const char* title() const override { return "O0 - Une impulsion sur une grille (corde et surface)"; }
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
    enum EdgeKind { kAbsorbing = 0, kWall = 1, kFree = 2 };

    double length() const;                       // côté du domaine [m]
    double cflNow() const { return dimension_ == 1 ? cfl1D_ : cfl2D_; }
    double cflLimit() const;
    double dtNow() const;
    double bump1D(double x) const;
    double exact1D(double x, double t) const;    // d'Alembert + images : l'onde exacte dans un domaine à bords idéaux
    std::vector<double> cut() const;             // u le long de l'axe x, au milieu du domaine
    double energyNow() const;
    double maxAbsNow() const;
    double timeNow() const;
    double dxNow() const;
    void sample();

    // paramètres (unités SI)
    int dimension_ = 2;               // 1 = corde, 2 = surface
    double speed_ = 1.0;              // c [m/s]
    double width_ = 0.2;              // largeur (écart-type) de la bosse [m]
    int edge_ = kAbsorbing;
    int cells1D_ = 400;
    int cells2D_ = 120;
    double cfl1D_ = 0.9;
    double cfl2D_ = 0.6;
    double spongeFraction_ = 0.15;    // épaisseur de la couche absorbante, part du côté
    double spongeStrength_ = 0.03;    // « réflexion visée » aux hautes fréquences (voir spongeSigmaMax)
    int spongeOrder_ = 1;
    double timeScale_ = 1.0;
    bool running_ = true;

    // simulation
    std::unique_ptr<waves::Wave1D> wave1_;
    std::unique_ptr<waves::Wave2D> wave2_;
    double initialEnergy_ = 1.0;
    double accumulator_ = 0.0;
    double lastSample_ = 0.0;
    double endTime_ = 1.0;
    bool finished_ = false;
    bool diverged_ = false;
    Series energy_;                   // E / E0 en fonction du temps
};

}  // namespace pl
