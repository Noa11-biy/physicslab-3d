// O3 : équation d'onde 2D en milieu variable (waves::Wave2D) calculée par un compute shader OpenGL 4.5 (shaders/wave2d.comp), en float.
//
// Le CPU `double` (Wave2D) reste la référence : la classe en reprend les coefficients tels quels (faces kx, ky, éponge, bords, sources) et l'état courant
// (u^{n-1}, u^n), puis avance par pas identiques sur le GPU, avec trois tampons qui tournent (précédent, courant, suivant) sans copie. On relit le champ
// (float) pour le dessiner ou pour le comparer au CPU. Un contexte OpenGL 4.5 doit être actif ; appeler shutdown() tant qu'il existe.
#pragma once

#include <string>
#include <vector>

#include "physicslab/waves/Wave.hpp"

namespace pl {

class GpuWave2D {
public:
    // Compile `shaderDir`/wave2d.comp. Faux si le pilote refuse (message sur stderr).
    bool init(const std::string& shaderDir);
    void shutdown();
    bool ready() const { return program_ != 0; }

    // Envoie la géométrie, les coefficients, l'éponge, les bords, les sources et l'état (u^{n-1}, u^n) de `wave`, et repart de son instant.
    void upload(const waves::Wave2D& wave);

    // Avance de `count` pas (le calcul est mis en file : relire le champ attend sa fin). Chaque pas est une petite commande (< 1 ms pour 512 x 512 sur Intel UHD).
    void step(int count);

    // Champ u^n courant (nx * ny floats), après attente des pas en file.
    void readField(std::vector<float>& out) const;
    // Attend la fin des pas en file ; renvoie la durée écoulée du GPU pour ces pas (secondes), mesurée par requête GL_TIME_ELAPSED autour d'un lot de pas.
    double benchmark(int steps);

    int pointsX() const { return nx_; }
    int pointsY() const { return ny_; }
    long long steps() const { return steps_; }
    double time() const { return steps_ * dt_; }

private:
    void bindBuffers() const;
    void sourceValues(int* nodes, float* values, int& count) const;

    unsigned int program_ = 0;
    unsigned int field_[3] = {0, 0, 0};     // u en rotation : précédent, courant, suivant
    unsigned int coefX_ = 0, coefY_ = 0, sigma_ = 0;
    int slot_ = 1;                          // indice du tampon courant dans field_ (précédent : slot_ - 1, suivant : slot_ + 1 modulo 3)
    int nx_ = 0, ny_ = 0;
    double dt_ = 0.0, ratio_ = 0.0;          // ratio_ = dt² / dx²
    long long steps_ = 0;
    int freeEdge_[4] = {0, 0, 0, 0};         // gauche, droite, bas, haut : 1 = libre
    std::vector<waves::PointSource> sources_;
    int uNx_ = -1, uNy_ = -1, uDt_ = -1, uFree_ = -1, uSourceCount_ = -1, uSourceNode_ = -1, uSourceValue_ = -1;
    int cellsX_ = 0;
    double dx_ = 0.0;
};

}  // namespace pl
