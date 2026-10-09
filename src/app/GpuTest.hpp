// Validation du calcul GPU (M7) : `physicslab --gpu-test`. Compare les accélérations du compute shader à `nbody::accelerations`
// (CPU double) sur des états fixés, mesure les temps et le débit. Le contexte OpenGL 4.5 est ouvert par l'application.
#pragma once

#include <string>

namespace pl {

// Renvoie 0 si toutes les vérifications passent. `maxN` limite la taille la plus grande des mesures de temps.
int runGpuTest(const std::string& shaderDir, int maxN);

}  // namespace pl
