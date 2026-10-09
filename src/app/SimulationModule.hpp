// Interface commune des simulations : chaque module de physique la réalise et l'application
// n'a plus qu'à héberger le module actif (fenêtre, panneaux, caméra, niveau pédagogique).
#pragma once

#include "physicslab/core/Level.hpp"
#include "physicslab/render/Camera.hpp"
#include "physicslab/render/Renderer.hpp"

namespace pl {

struct UiContext {
    Level level;
    float uiScale;
};

class SimulationModule {
public:
    virtual ~SimulationModule() = default;

    virtual const char* domain() const { return "Mécanique"; }  // titre de section dans le menu Simulation
    virtual const char* title() const = 0;                // entrée de menu et en-tête du panneau Explication
    virtual const char* explanation(Level level) const = 0;

    virtual void reset() = 0;
    virtual void update(double frameSeconds) = 0;         // avance la simulation de la durée d'une image

    // Contenu des panneaux : l'application ouvre la fenêtre ImGui, le module la remplit.
    virtual void drawControls(const UiContext& ctx) = 0;    // "Simulation"
    virtual void drawInvariants(const UiContext& ctx) = 0;  // "Invariants"
    virtual void drawGraphs(const UiContext& ctx) = 0;      // "Graphes"
    virtual bool hasAnalysis(const UiContext&) const { return false; }
    virtual void drawAnalysis(const UiContext&) {}          // "Analyse" (affiché si hasAnalysis)

    // Dessin 3D, appelé après Renderer::beginFrame.
    virtual void drawScene(Renderer& renderer, const UiContext& ctx) = 0;
    virtual void frameCamera(Camera& camera) const = 0;     // cadrage initial de la scène
};

}  // namespace pl
