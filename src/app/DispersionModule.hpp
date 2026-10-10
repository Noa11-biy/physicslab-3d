// O2 : dispersion numérique et condition CFL. Une bosse (large bande) ou un paquet d'ondes (une longueur d'onde) se propage sur la grille du saute-mouton
// d'O0 : la grille ralentit les ondes courtes (vitesse de phase et de groupe < c), déforme la bosse (ondulations derrière elle), rend la surface 2D
// anisotrope (la diagonale ne va pas à la même vitesse qu'un axe) et explose si le nombre de Courant C = c dt / dx dépasse la limite (CFL).
// Les formules (waves/Dispersion) sont confrontées au calcul : paquet contre vitesse de groupe, étude de convergence, taux de croissance de l'instabilité.
#pragma once

#include <memory>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/waves/Wave.hpp"

namespace pl {

class DispersionModule final : public SimulationModule {
public:
    DispersionModule();

    const char* domain() const override { return "Ondes"; }
    const char* title() const override { return "O2 - Dispersion numérique et condition CFL"; }
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
    enum Shape { kPulse = 0, kPacket = 1 };

    double cflNow() const { return dimension_ == 1 ? cfl1D_ : cfl2D_; }
    double cflLimit() const;
    double dxNow() const { return dimension_ == 1 ? wave1_->dx() : wave2_->dx(); }
    double dtNow() const { return dimension_ == 1 ? wave1_->dt() : wave2_->dt(); }
    double timeNow() const { return dimension_ == 1 ? wave1_->time() : wave2_->time(); }
    double maxAbsNow() const { return dimension_ == 1 ? wave1_->maxAbs() : wave2_->maxAbs(); }
    double carrierPhase() const;                 // k0 dx de la porteuse (ou 2 pi / cases par longueur d'onde)
    double exact1D(double x, double t) const;    // onde exacte du continu (d'Alembert), sans erreur de grille
    double lag1D() const;                        // retard [m] du centre d'énergie du calcul sur le centre exact (moitié droite)
    void frontRadii2D(double& axis, double& diagonal) const;  // rayons du front de l'anneau le long d'un axe et de la diagonale [m]
    void computeConvergence();
    double pulseError(int cellsPerWavelength, double cfl, double time) const;  // étude 1D : erreur max d'une bosse contre d'Alembert

    // paramètres
    int dimension_ = 1;
    int shape_ = kPulse;
    int ppw_ = 5;                         // cases par longueur d'onde (0,8 m) : assez grossier pour que la déformation se voie
    double cfl1D_ = 0.5;
    double cfl2D_ = 0.5;
    double timeScale_ = 1.0;
    bool running_ = true;

    // simulation
    std::unique_ptr<waves::Wave1D> wave1_;
    std::unique_ptr<waves::Wave2D> wave2_;
    double accumulator_ = 0.0;
    double endTime_ = 1.0;
    bool finished_ = false;
    bool diverged_ = false;
    std::vector<double> peaks_;           // max |u| à chaque pas (pour mesurer la croissance d'une instabilité)

    // étude de convergence (Analyse)
    struct Convergence {
        std::vector<double> dx, err;      // erreur d'une bosse selon la taille des cases, au C actuel
        std::vector<double> cfl, errC;    // erreur à 24 cases par longueur d'onde selon C
        double slope = 0.0;
        bool exactScheme = false;         // C = 1 : erreur d'arrondi, pas d'ordre à mesurer
    } conv_;
    bool convDirty_ = true;
    double convCfl_ = -1.0;
};

}  // namespace pl
