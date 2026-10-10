// Compilation d'un compute shader OpenGL 4.5 depuis un fichier, partagée par GpuNBody (M7) et GpuWave2D (O3).
#pragma once

#include <string>

namespace pl {

// Lit `path`, insère `defines` juste après la ligne #version (le fichier doit commencer par #version), compile et lie un programme de calcul.
// Renvoie 0 en cas d'échec (message sur stderr, préfixé par `tag`). Un contexte OpenGL 4.5 doit être actif.
unsigned int buildComputeProgram(const std::string& path, const std::string& defines, const char* tag);

}  // namespace pl
