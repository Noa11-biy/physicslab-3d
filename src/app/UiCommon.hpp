// Éléments d'interface partagés par les modules : courbes temps réel, curseurs, et l'ensemble
// des 5 solveurs comparés (couleurs, visibilité et libellés selon le niveau pédagogique).
#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <imgui.h>
#include <implot.h>

#include "physicslab/core/Level.hpp"
#include "physicslab/core/Solver.hpp"

namespace pl {

inline constexpr float kBlue[3] = {0.35f, 0.65f, 1.0f};  // couleur de la solution exacte

inline ImVec4 toImVec4(const float c[3], float alpha = 1.0f) { return ImVec4(c[0], c[1], c[2], alpha); }

// Texte formaté à la printf.
std::string strf(const char* format, ...);

// Style de courbe ImPlot : couleur, épaisseur 2, décalage du tampon circulaire.
inline ImPlotSpec lineSpec(const float c[3], int offset = 0) {
    return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Offset, offset);
}

// Curseur sur un double (ImGui::SliderFloat ne prend que des float).
bool sliderD(const char* label, double* v, double lo, double hi, const char* fmt, ImGuiSliderFlags flags = 0);

// Largeur des curseurs : colonne de libellés calée sur le plus long, curseurs alignés et textes jamais coupés.
void pushSliderWidth();
void popSliderWidth();

// Courbe temps réel dans un tampon circulaire (compatible avec l'argument `offset` d'ImPlot).
struct Series {
    static constexpr int kCapacity = 4000;
    std::vector<double> x, y;
    int offset = 0;

    int size() const { return static_cast<int>(x.size()); }
    void clear() { x.clear(); y.clear(); offset = 0; }
    void add(double xx, double yy) {
        if (size() < kCapacity) {
            x.push_back(xx);
            y.push_back(yy);
        } else {
            x[offset] = xx;
            y[offset] = yy;
            offset = (offset + 1) % kCapacity;
        }
    }
};

// Horloge à pas fixe pour les simulations : accumule le temps des images, avance par pas de dt, échantillonne
// au plus tous les 1/240 s, et raccourcit le dernier pas pour s'arrêter pile à endTime.
struct StepClock {
    double time = 0.0;
    double accumulator = 0.0;
    double lastSample = 0.0;
    bool finished = false;

    void reset() { time = accumulator = lastSample = 0.0; finished = false; }

    // step(h) fait avancer tous les solveurs de h ; sample() enregistre les courbes.
    template <class StepFn, class SampleFn>
    void advance(double frameSeconds, double timeScale, double dt, double endTime, StepFn&& step, SampleFn&& sample) {
        if (finished) return;
        accumulator += std::min(frameSeconds, 0.1) * timeScale;  // borne anti "spirale de la mort"
        const double sampleInterval = std::max(dt, 1.0 / 240.0);
        int guard = 0;
        while (!finished && guard++ < 5000) {
            const double remaining = endTime - time;
            const double h = std::min(dt, remaining);
            if (accumulator < h) break;
            step(h);
            time += h;
            accumulator -= h;
            finished = remaining <= dt * (1.0 + 1e-9);
            if (finished || time - lastSample >= sampleInterval - 1e-12) {
                lastSample = time;
                sample();
            }
        }
    }
};

// Tableau de résultats : une ligne par méthode, première colonne "Méthode" (colorée), puis des colonnes numériques.
// Une cellule "-" s'affiche en grisé.
struct TableRow {
    std::string name;
    const float* color = nullptr;
    std::vector<std::string> cells;  // une par en-tête
};
void drawResultTable(const char* id, const std::vector<std::string>& headers, const std::vector<TableRow>& rows);

// Euler, Euler symplectique, Verlet, RK4, RK45 : instances, couleurs et choix d'affichage.
class SolverSet {
public:
    enum Index { kEuler = 0, kSymplectic, kVerlet, kRK4, kRK45, kCount };
    static constexpr int kFixedStep = 4;  // les 4 premiers sont à pas fixe (étude de convergence)

    SolverSet();

    Solver& solver(int i) { return *entries_[i].solver; }
    RK45& rk45() { return *rk45_; }
    const float* color(int i) const { return entries_[i].color; }
    bool& show(int i) { return entries_[i].show; }

    // Niveaux 1-2 : le seul "ordinateur" (RK4). Niveaux 3-4 : Euler contre RK4. Niveau 5+ : au choix.
    bool isShown(Level level, int i) const;
    std::string label(Level level, int i) const;
    // Version abrégée pour les tableaux étroits.
    std::string shortLabel(Level level, int i) const;

    // Cases à cocher des méthodes (niveau 5+).
    void drawToggles(Level level);

    static std::unique_ptr<Solver> makeFixedStep(int i);

private:
    struct Entry {
        std::unique_ptr<Solver> solver;
        float color[3] = {1, 1, 1};
        bool show = true;
    };
    std::array<Entry, kCount> entries_;
    RK45* rk45_ = nullptr;
};

}  // namespace pl
