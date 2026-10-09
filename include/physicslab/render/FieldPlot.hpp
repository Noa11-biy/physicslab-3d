// Dessin d'un champ scalaire (hauteur d'une onde) : des sommets pour le Renderer, sans OpenGL ici.
// Palette divergente : bleu pour les valeurs négatives (creux), gris sombre à zéro, rouge pour les positives (bosses).
#pragma once

#include <vector>

#include "physicslab/render/Renderer.hpp"

namespace pl {

// Couleur de la palette divergente pour t dans [-1, 1] (au-delà : saturée).
void divergingColor(double t, float rgb[3]);

// Profil d'un champ 1D : y = heightScale u(x), en LineStrip, centré sur x = 0 (le point i est en (i - (n-1)/2) dx). La couleur suit u / colorScale.
std::vector<Vertex> makeProfile(const std::vector<double>& u, double dx, float heightScale, float colorScale, float z = 0.0f);

// Fil de fer coloré d'un champ 2D (nx × ny points, rangés ligne par ligne) : une ligne du maillage sur `stride`, en Lines.
// Le champ est centré sur l'origine du plan (x, z) ; la hauteur y = heightScale u ; la couleur suit u / colorScale.
std::vector<Vertex> makeRelief(const std::vector<double>& u, int nx, int ny, double dx, float heightScale, float colorScale, int stride);

}  // namespace pl
