// O4 : fentes de Young et diffraction. Une onde plane (source en ligne à gauche) arrive sur un écran percé de 1 à 5 fentes ; on regarde ce qui sort : la tache de diffraction
// d'une fente (largeur b), les franges d'interférence de deux fentes (pas a), puis le réseau. L'intensité est lue sur un écran d'observation à la distance D = 12 m derrière
// les fentes et comparée à la théorie de Fraunhofer (champ lointain), à la formule des franges (interfrange i = lambda D / a) et, au niveau 6, à l'intégrale de
// Rayleigh-Sommerfeld avec la phase de la grille. Le calcul se fait comme en O3 : CPU en double, GPU en float, ou les deux (le CPU vérifie le GPU).
#pragma once

#include <chrono>
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/render/GpuWave2D.hpp"
#include "physicslab/waves/Diffraction.hpp"
#include "physicslab/waves/Wave.hpp"

namespace pl {

class YoungModule final : public SimulationModule {
public:
    explicit YoungModule(std::string shaderDir);
    ~YoungModule() override;

    const char* domain() const override { return "Ondes"; }
    const char* title() const override { return "O4 - Fentes de Young et diffraction"; }
    const char* explanation(Level level) const override;

    void reset() override;
    void update(double frameSeconds) override;

    void drawControls(const UiContext& ctx) override;
    void drawInvariants(const UiContext& ctx) override;
    void drawGraphs(const UiContext& ctx) override;
    bool hasAnalysis(const UiContext&) const override { return false; }
    void drawAnalysis(const UiContext&) override {}

    void drawScene(Renderer& renderer, const UiContext& ctx) override;
    void frameCamera(Camera& camera) const override;

private:
    enum Backend { kCpu = 0, kGpu, kBoth };

    // géométrie [m] : cuve de 26 x 24, source en ligne en x = 3, écran en x = 7, écran d'observation à D derrière l'écran, éponge de 2,5 m
    static constexpr double kWidth = 26.0, kHeight = 24.0, kSourceX = 3.0, kScreenX = 7.0, kObsDistance = 12.0, kSponge = 2.5;
    static constexpr double kCfl = 0.625, kDx = 0.1;   // 10 cases par mètre ; la période dure P pas, P multiple de 4 : lambda = 0.0625 P (m)
    static constexpr double kViewHalf = 9.0;            // l'écran d'observation est tracé pour |y - yc| <= 9 m

    double lambda() const { return 0.0625 * 4 * quarter_; }
    double timeNow() const { return backend_ == kGpu ? gpu_.time() : cpu_->time(); }
    int nodesX() const { return cpu_->pointsX(); }
    int nodesY() const { return cpu_->pointsY(); }
    const std::vector<double>& field() const { return *shown_; }
    double obsX() const { return screen_.centerX + kObsDistance; }
    double readyTime() const { return (kObsDistance + kScreenX - kSourceX) / 1.0 + 8.0; }  // front d'onde + montée en 4 périodes + marge
    bool steady() const { return timeNow() >= readyTime() && snapshots_.size() == 4; }

    void refreshField(bool compare);   // relit le champ GPU ; compare = mesurer l'écart GPU / CPU (mode « les deux »)
    void takeSnapshot();
    void computeProfile();        // amplitude sur l'écran d'observation, intensité normalisée, minima mesurés
    void ensureGridTheory();      // Rayleigh-Sommerfeld avec la phase de la grille (niveau 6), calculé à la demande

    // paramètres
    int slits_ = 2;
    double width_ = 2.0;          // b [m]
    double pitch_ = 4.0;          // a [m]
    int quarter_ = 4;             // quart de période en pas : lambda = 1 m par défaut (P = 16)
    int backend_ = kGpu;
    double timeScale_ = 2.0;
    bool running_ = true;
    bool showFraunhofer_ = true;

    // simulation
    std::string shaderDir_;
    GpuWave2D gpu_;
    bool gpuOk_ = false;
    waves::diffraction::SlitScreenSpec spec_;
    waves::diffraction::SlitScreen screen_;
    std::unique_ptr<waves::Wave2D> cpu_;
    std::vector<float> gpuField_;
    std::vector<double> gpuFieldD_;
    std::vector<double> display_;
    const std::vector<double>* shown_ = nullptr;
    double accumulator_ = 0.0;
    double endTime_ = 120.0;
    bool finished_ = false;
    int maxStepsPerFrame_ = 24;
    double cpuMsPerStep_ = 0.0, gpuMsPerStep_ = 0.0;
    double gpuVsCpuMax_ = 0.0, gpuVsCpuL2_ = 0.0;

    // mesures sur l'écran d'observation (colonne de noeuds en obsX) : quatre instantanés espacés d'un quart de période donnent l'amplitude A = hypot(u0 - u2, u1 - u3) / 2
    std::deque<std::vector<double>> snapshots_;
    int obsColumn_ = 0;
    std::vector<double> profileY_;        // y - yc [m] (|y| <= kViewHalf)
    std::vector<double> profileI_;        // intensité mesurée, normalisée au maximum
    std::vector<double> theoryY_;         // abscisses (y - yc) des courbes théoriques : 361 points entre -9 et 9 m
    std::vector<double> fraunhoferI_;     // théorie de Fraunhofer sur l'écran (sin theta = y / sqrt(D² + y²)), maximum 1
    std::vector<double> gridI_;           // Rayleigh-Sommerfeld avec la phase de la grille, normalisée
    bool gridValid_ = false;
    std::unique_ptr<waves::dispersion::RayPhase2D> ray_;
    std::vector<double> minimaY_;         // minima mesurés (y - yc) du profil
    double measuredFringe_ = 0.0;         // interfrange mesuré (écart entre les deux minima autour du centre), 0 : pas encore
};

}  // namespace pl
