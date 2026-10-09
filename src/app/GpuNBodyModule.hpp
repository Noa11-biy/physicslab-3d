// M7 : gravitation à N corps sur carte graphique. Un amas de Plummer (ou deux amas qui se rencontrent) de N étoiles, intégré par
// kick-drift-kick (Verlet des vitesses) soit par le GPU en float (compute shader, GpuNBody), soit par le CPU en double (la
// référence), soit par les deux côte à côte : « le CPU vérifie le GPU ». Unités normalisées : G = 1, masse totale 1, a = 1.
#pragma once

#include <chrono>
#include <string>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/mechanics/NBody.hpp"
#include "physicslab/render/GpuNBody.hpp"

namespace pl {

class GpuNBodyModule final : public SimulationModule {
public:
    explicit GpuNBodyModule(std::string shaderDir);
    ~GpuNBodyModule() override;

    const char* title() const override { return "M7 - Gravitation : N corps sur carte graphique"; }
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
    enum Scenario { kCluster = 0, kCollision };
    enum Backend { kGpu = 0, kCpu, kBoth };

    // Un calcul (GPU ou CPU) : son état et l'historique de ses grandeurs conservées.
    struct Track {
        State y;                           // [positions | vitesses] en double (relu du GPU pour la piste GPU)
        Series energy, momentum, angular;  // (E - E0)/|E0| ; |P - P0| / sum m|v| ; |L - L0| / sum m|r x v|
        Series stepMs;                     // temps de calcul par pas (ms)
        double stepMillis = 0.0;           // moyenne glissante de stepMs (temps vécu : attente du GPU comprise)
        double kernelMillis = 0.0;         // piste GPU : moyenne glissante du temps des seuls calculs GPU par pas
        double lastEnergy = 0.0, lastMomentum = 0.0, lastAngular = 0.0;   // dernières valeurs relatives
        bool diverged = false;
        void clear();
    };

    // Écart entre les accélérations du shader et celles du CPU en double, au même état.
    struct Validation {
        bool valid = false;
        int count = 0;
        double time = 0.0;                       // instant simulé de la mesure
        double maxRel = 0.0, medianRel = 0.0, p99Rel = 0.0, globalRel = 0.0;
        double potentialCpu = 0.0, potentialGpu = 0.0;
        double cpuMillis = 0.0;
        std::vector<double> centers, counts;     // histogramme de log10 de l'erreur relative par étoile
    };

    // Une mesure du banc d'essai (bouton du niveau 6).
    struct BenchRow {
        int n = 0;
        double gpuGints = 0.0;                   // milliards d'interactions par seconde
        double cpuGints = -1.0;                  // < 0 : non mesuré
        double error = -1.0;                     // écart global des accélérations GPU float / CPU double
    };

    int bodyCount() const { return problem_.count(); }
    bool useGpu() const { return gpuOk_ && backend_ != kCpu; }
    bool useCpu() const { return backend_ != kGpu || !gpuOk_; }
    int maxBodies() const;

    void buildProblem();
    void stepGpu(int steps);
    void stepCpu(int steps);
    void recordConserved(Track& t);
    void recordEnergy(Track& t, double potential);
    void runValidation();
    void runBenchmark();
    Vertex starVertex(const State& y, int i, const float* color, float speedScale) const;

    // paramètres
    std::string shaderDir_;
    Scenario scenario_ = kCluster;
    Backend backend_ = kGpu;
    double bodies_ = 2000.0;      // nombre d'étoiles
    unsigned seed_ = 7;
    double softening_ = 0.05;
    double dt_ = 0.004;
    double timeScale_ = 1.0;
    bool running_ = true;
    bool pendingReset_ = false;       // un paramètre a changé : on recrée l'amas dès qu'aucun curseur n'est tenu
    Level level_ = Level::College;

    // simulation
    GpuNBody gpu_;
    bool gpuOk_ = false;
    NBodyProblem problem_;
    Track gpuTrack_, cpuTrack_;
    std::vector<double> cpuAcc_;
    bool cpuAccValid_ = false;
    Series gap_;                      // distance dans l'espace des phases entre GPU et CPU (mode « les deux »)
    double time_ = 0.0, accumulator_ = 0.0;
    double achieved_ = 1.0;           // part des pas demandés réellement calculés (1 = temps réel)
    double e0_ = -1.0, pScale_ = 1.0, lScale_ = 1.0, speedScale_ = 1.0;
    Vec3 p0_, l0_;
    bool exactEnergy_ = true;         // E0 calculée en double (sinon estimée par le shader)
    std::chrono::steady_clock::time_point lastCpuEnergy_;
    std::chrono::steady_clock::time_point lastValidation_;

    // analyse
    Validation validation_;
    std::vector<BenchRow> bench_;
};

}  // namespace pl
