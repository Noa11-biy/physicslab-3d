// O3 : cuve à ondes. Une surface d'eau 2D dont la vitesse des vagues c(x, y) varie (zone moins profonde = plus lente) : une interface plane ou une lentille.
// Une source ponctuelle (rides circulaires, amplitude en 1/sqrt(r)) ou un paquet d'ondes plan (réflexion, réfraction de Snell-Descartes, réflexion totale).
// Le calcul se fait au choix sur le CPU en double (la référence), sur le GPU en float (compute shader), ou sur les deux côte à côte : « le CPU vérifie le GPU ».
#pragma once

#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include "SimulationModule.hpp"
#include "UiCommon.hpp"
#include "physicslab/render/GpuWave2D.hpp"
#include "physicslab/waves/Wave.hpp"

namespace pl {

class WaveTankModule final : public SimulationModule {
public:
    explicit WaveTankModule(std::string shaderDir);
    ~WaveTankModule() override;

    const char* domain() const override { return "Ondes"; }
    const char* title() const override { return "O3 - Cuve à ondes : réflexion et réfraction"; }
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
    enum Medium { kUniform = 0, kInterface, kLens };
    enum Source { kPoint = 0, kPacket };
    enum Backend { kCpu = 0, kGpu, kBoth };

    // géométrie [m] : cuve de 13 x 10 ; interface verticale en x = 7 ; lentille de rayon 1,5 centrée en (7,5 ; 5)
    static constexpr double kWidth = 13.0, kHeight = 10.0, kInterfaceX = 7.0, kLensX = 7.5, kLensRadius = 1.5;
    static constexpr double kPacketSigma = 0.4, kPacketDistance = 1.8;  // enveloppe du paquet ; distance de son centre à l'interface au départ

    double cLeft() const;            // vitesse à gauche de l'interface (ou partout, sans interface)
    double cRight() const;
    double speedAt(double x, double y) const;
    double wavelength() const { return cLeft() / freq_; }   // longueur d'onde du paquet dans le milieu d'arrivée
    double sourceX() const { return medium_ == kUniform ? 6.5 : 3.0; }  // position de la source ponctuelle
    double sourceY() const { return 5.0; }
    double timeNow() const { return backend_ == kGpu ? gpu_.time() : cpu_->time(); }
    bool packetCrossesInterface() const { return source_ == kPacket && medium_ == kInterface; }
    const std::vector<double>& field() const { return *shown_; }
    int nodesX() const { return 26 * ppw_ + 1; }
    int nodesY() const { return 20 * ppw_ + 1; }

    void refreshField();
    void regionStats(double xmin, double xmax, double& energy, double& cx, double& cy) const;
    void sampleSeries();
    void measureAngles();
    void updateAmplitudes();

    // paramètres
    int medium_ = kInterface;
    int source_ = kPoint;
    int backend_ = kCpu;
    int ppw_ = 12;                  // cases par 0,5 m (finesse) : 312 x 240 noeuds par défaut
    double freq_ = 2.0;             // fréquence [Hz]
    double ratio_ = 0.6;            // vitesse du milieu lent / vitesse du milieu rapide
    bool swapMedia_ = false;        // milieu lent à gauche (le paquet arrive de la zone lente : réflexion totale possible)
    double angleDeg_ = 30.0;        // angle d'incidence du paquet
    double timeScale_ = 1.0;
    bool running_ = true;

    // simulation
    std::string shaderDir_;
    GpuWave2D gpu_;
    bool gpuOk_ = false;
    std::unique_ptr<waves::Wave2D> cpu_;
    std::vector<float> gpuField_;
    std::vector<double> gpuFieldD_;
    std::vector<double> display_;   // champ écrêté pour le dessin
    const std::vector<double>* shown_ = nullptr;
    double accumulator_ = 0.0;
    double endTime_ = 60.0;
    bool finished_ = false;
    int maxStepsPerFrame_ = 20;     // adapté pour garder l'interface fluide (le CPU double est lent en build Debug)
    double cpuMsPerStep_ = 0.0, gpuMsPerStep_ = 0.0;
    double gpuVsCpuMax_ = 0.0, gpuVsCpuL2_ = 0.0;   // écarts relatifs en mode « les deux »

    // mesures
    double totalEnergy0_ = 1.0;
    double hitTime_ = 0.0;           // instant où le centre du paquet atteint l'interface
    Series leftEnergy_, rightEnergy_;
    double lastSeriesTime_ = -1.0;
    bool haveT_ = false, haveR_ = false;
    double txA_ = 0, tyA_ = 0, rxA_ = 0, ryA_ = 0, tA_ = 0, rA_ = 0;   // centres mémorisés (transmis, réfléchi) et leurs instants
    double measuredT_ = 0.0, measuredR_ = 0.0;                           // angles mesurés [degrés par rapport à la normale], 0 : pas encore
    double measuredTEnergy_ = 0.0;
    // amplitude de la source ponctuelle le long de +x : maximum de |u| sur les deux dernières périodes
    std::vector<double> ampRadii_;
    std::vector<double> amplitude_;
    std::vector<std::pair<double, std::vector<float>>> ampHistory_;
};

}  // namespace pl
