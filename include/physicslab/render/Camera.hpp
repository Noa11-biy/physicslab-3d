// Caméra orbitale (axe Y vers le haut). Calculs en float : c'est du rendu, pas de la physique.
#pragma once

namespace pl {

class Camera {
public:
    float target[3] = {10.0f, 3.0f, 0.0f};  // point regardé
    float distance = 38.0f;                 // distance à la cible
    float yaw = 0.55f;                      // rotation autour de Y [rad]
    float pitch = 0.30f;                    // élévation [rad]
    float fovY = 0.8f;                      // champ de vision vertical [rad]

    void orbit(float dxPixels, float dyPixels);  // clic gauche
    void pan(float dxPixels, float dyPixels);    // clic droit / milieu
    void zoom(float wheelTicks);                 // molette

    // Matrice projection * vue, colonnes d'abord (format OpenGL).
    void viewProjection(float aspect, float out[16]) const;
};

}  // namespace pl
