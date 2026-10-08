// Rendu OpenGL 4.5 minimal : lignes et points colorés. Ne connaît rien de la physique.
#pragma once

#include <string>
#include <vector>

#include "physicslab/render/Camera.hpp"

namespace pl {

struct Vertex {
    float x, y, z;  // position
    float r, g, b;  // couleur
};

// Rectangle en pixels du framebuffer, origine en haut à gauche de la fenêtre.
struct ViewRect {
    int x, y, w, h;
};

enum class Primitive { Lines, LineStrip, Points };

class Renderer {
public:
    // Charge et compile les shaders depuis `shaderDir`. Un contexte OpenGL 4.5 doit être actif.
    bool init(const std::string& shaderDir);
    void shutdown();

    // Efface toute la fenêtre puis restreint le dessin à `view`.
    void beginFrame(int framebufferWidth, int framebufferHeight, const ViewRect& view,
                    const Camera& camera, const float clearColor[3]);

    void draw(Primitive primitive, const std::vector<Vertex>& vertices, float pointSize = 1.0f);

private:
    unsigned int vao_ = 0;
    unsigned int vbo_ = 0;
    unsigned int program_ = 0;
    int uViewProj_ = -1;
    int uPointSize_ = -1;
    int uRound_ = -1;
};

// Grille dans le plan y = 0 : (2 * halfCells) cellules de côté cellSize.
std::vector<Vertex> makeGrid(int halfCells, float cellSize);

// Axes X (rouge), Y (vert), Z (bleu) depuis l'origine.
std::vector<Vertex> makeAxes(float length);

}  // namespace pl
