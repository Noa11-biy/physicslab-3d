// M7 : gravitation à N corps calculée par des compute shaders OpenGL 4.5 (shaders/nbody.comp et nbody_step.comp).
//
// Même contrat que `nbody::accelerations` (M4b, CPU double) : positions, masses, n, G, adoucissement de Plummer -> accélérations.
// Le calcul reste en O(N^2) : un thread par corps cible, boucle sur les sources par tuiles chargées en mémoire partagée.
// Précision `Float` (usage normal, gros volumes) ou `Double` (si le pilote l'accepte ; émulé donc très lent sur certaines cartes :
// sert à séparer l'erreur d'arrondi de l'erreur de parallélisme). Les entrées sont des `double` : le centrage sur le barycentre se
// fait en double avant la conversion. En plus des accélérations, la classe intègre le mouvement sur le GPU par kick-drift-kick
// (Verlet des vitesses, ordre 2, symplectique) avec positions, vitesses et accélérations qui restent dans des SSBO.
// Un contexte OpenGL 4.5 doit être actif ; aucune destruction automatique (appeler shutdown() tant que le contexte existe).
#pragma once

#include <string>
#include <utility>
#include <vector>

namespace pl {

enum class GpuPrecision { Float, Double };

// Capacités du pilote utiles au calcul (valeurs lues à l'exécution, rien n'est supposé).
struct GpuLimits {
    std::string renderer, version;
    int maxGroupInvocations = 0;     // threads par groupe de travail
    int maxGroupSize[3] = {0, 0, 0};
    int maxSharedBytes = 0;          // mémoire partagée par groupe
    long long maxStorageBlockBytes = 0;
    bool doublePrecision = false;    // extension GL_ARB_gpu_shader_fp64 annoncée (cœur depuis GL 4.0)
};
GpuLimits queryGpuLimits();

// Positions de `n` corps (3n doubles) translatées pour que le barycentre soit à l'origine ; sert au GPU et aux tests.
// Si `center` n'est pas nul, il reçoit le barycentre retranché (3 doubles).
void centerOnBarycenter(const double* positions, const double* masses, int n, double* out, double* center = nullptr);

class GpuNBody {
public:
    // Compile `shaderDir`/nbody.comp et nbody_step.comp avec un groupe de `workGroupSize` threads. Faux si le pilote refuse
    // (message sur stderr). 256 threads par groupe : le plus rapide mesuré sur Intel UHD (N = 16384 : 8,3 ms contre 9,2 ms à 128
    // et 17,7 ms à 32).
    bool init(const std::string& shaderDir, GpuPrecision precision = GpuPrecision::Float, int workGroupSize = 256);
    void shutdown();
    bool ready() const { return program_ != 0 && stepProgram_ != 0; }

    // Centrer sur le barycentre avant la conversion en float (oui par défaut). Les accélérations ne dépendent pas de l'origine :
    // le centrage ne sert qu'à limiter l'erreur d'arrondi sur les positions et sur les différences r_j - r_i.
    void setCentering(bool on) { center_ = on; }

    // --- accélérations seules (à positions données) ---
    // Envoie positions (3n doubles) et masses (n doubles) au GPU.
    void setBodies(const double* positions, const double* masses, int n);
    // Lance le calcul (G, adoucissement) et attend sa fin ; le temps de calcul des shaders est mesuré par requête GL_TIME_ELAPSED.
    void compute(double G, double softening);
    // Relit les accélérations (3n doubles) du dernier compute().
    void readAccelerations(double* acc) const;
    // Équivalent GPU de nbody::accelerations (envoi + calcul + relecture).
    void accelerations(const double* positions, const double* masses, int n, double G, double softening, double* acc);
    // Énergie potentielle (adoucie) des positions pour lesquelles les accélérations viennent d'être calculées (compute(), ou fin de
    // step()) : U = 1/2 sum m_i phi_i, où le shader range phi_i dans la 4e composante des accélérations. Chaque phi_i est sommé en
    // float (ou double) sur le GPU, la somme sur les corps est faite ici en double. Renvoie 0 si les accélérations ne sont pas à jour.
    double potentialEnergy() const;

    // --- mouvement ---
    // Envoie l'état initial : positions (3n), vitesses (3n), masses (n). Les positions sont relues dans le repère d'origine.
    void setState(const double* positions, const double* velocities, const double* masses, int n);
    // Avance de `steps` pas de `dt` (kick-drift-kick) : par pas, 1 envoi de calcul des accélérations et 1 de mise à jour, tout
    // sur le GPU. Au retour l'état est synchrone (vitesses et positions au même instant). Les accélérations du dernier point
    // restent en mémoire : un nouvel appel avec les mêmes G et adoucissement ne les recalcule pas.
    void step(double G, double softening, double dt, int steps);
    // Relit positions et vitesses (3n doubles chacune).
    void readState(double* positions, double* velocities) const;

    int count() const { return n_; }
    // Temps GPU du dernier compute() ou step() : calcul des accélérations seulement (la mise à jour des positions, élément par
    // élément, est négligeable devant le O(N^2)), hors transferts.
    double lastKernelSeconds() const { return kernelSeconds_; }
    int workGroupSize() const { return wg_; }
    GpuPrecision precision() const { return precision_; }

private:
    void uploadBodies(const double* positions, const double* masses, int n);
    void setReal(unsigned int program, int location, double value) const;
    void dispatchStepKernel(double kick, double drift);
    // Envoie le calcul des accélérations, en envois courts. `blocking` : attend chaque envoi (mesure du temps immédiate) ; sinon les
    // envois s'empilent sans synchronisation et leurs temps sont relus plus tard par collectTimings().
    void dispatchAccel(double G, double softening, bool blocking);
    void collectTimings();   // attend les envois en cours, ajoute leur temps à kernelSeconds_ et remet à jour le débit mesuré

    unsigned int program_ = 0;      // nbody.comp : accélérations
    unsigned int stepProgram_ = 0;  // nbody_step.comp : coup de pied et dérive
    unsigned int bodies_ = 0;       // SSBO 0 : (x, y, z, masse) par corps
    unsigned int accel_ = 0;        // SSBO 1 : (ax, ay, az, 0) par corps
    unsigned int velocities_ = 0;   // SSBO 2 : (vx, vy, vz, 0) par corps
    std::vector<unsigned int> queryPool_;                       // requêtes GL_TIME_ELAPSED réutilisées
    std::vector<std::pair<unsigned int, double>> pending_;      // (requête, interactions) des envois pas encore relus
    double queuedSeconds_ = 0.0;                                // durée estimée du travail empilé
    int uCount_ = -1, uFirst_ = -1, uG_ = -1, uEps2_ = -1;
    int uStepCount_ = -1, uDtKick_ = -1, uDtDrift_ = -1;
    int n_ = 0;
    int capacity_ = 0;          // corps pour lesquels les tampons sont dimensionnés
    int wg_ = 256;
    bool center_ = true;
    bool stateReady_ = false;   // des vitesses ont été envoyées (setState)
    bool accelValid_ = false;   // accel_ correspond aux positions actuelles pour (accelG_, accelEps_)
    double accelG_ = 0.0, accelEps_ = 0.0;
    double offset_[3] = {0.0, 0.0, 0.0};  // barycentre retranché à l'envoi
    GpuPrecision precision_ = GpuPrecision::Float;
    double kernelSeconds_ = 0.0;
    // Une seule commande GPU ne doit pas durer plus d'environ 2 s (Windows relance le pilote graphique au-delà) : le calcul est
    // découpé en envois d'environ `targetDispatchSeconds_`. Le débit (interactions/s) est mesuré à chaque envoi, en commençant
    // par un seul groupe de travail, car il varie de plusieurs ordres de grandeur d'une carte et d'une précision à l'autre.
    double targetDispatchSeconds_ = 0.05;
    double rate_ = 0.0;
};

}  // namespace pl
