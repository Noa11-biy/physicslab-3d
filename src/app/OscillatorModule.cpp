#include "OscillatorModule.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

#include <implot.h>

#include "physicslab/core/Constants.hpp"

namespace pl {
namespace {

constexpr double kTwoPi = 2.0 * constants::pi;
constexpr float kLaneSpacing = 2.6f;   // distance entre deux couloirs [m]
constexpr float kWallX = -6.0f;        // position du mur

float laneZ(int lane) { return kLaneSpacing * static_cast<float>(lane); }

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

// Arêtes d'un cube de côté 1 centré en (cx, cy, cz).
void addCube(std::vector<Vertex>& v, float cx, float cy, float cz, const float* c) {
    const float h = 0.5f;
    const float xs[2] = {cx - h, cx + h}, ys[2] = {cy - h, cy + h}, zs[2] = {cz - h, cz + h};
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            addSegment(v, xs[0], ys[a], zs[b], xs[1], ys[a], zs[b], c);  // arêtes parallèles à X
            addSegment(v, xs[a], ys[0], zs[b], xs[a], ys[1], zs[b], c);  // parallèles à Y
            addSegment(v, xs[a], ys[b], zs[0], xs[a], ys[b], zs[1], c);  // parallèles à Z
        }
}

// Un couloir : mur, ressort en hélice, masse cubique et repère de la position d'équilibre.
void drawLane(Renderer& renderer, int lane, double position, const float* color) {
    const float z = laneZ(lane);
    const float x = static_cast<float>(std::clamp(position, -12.0, 12.0));  // une méthode qui diverge reste dans le cadre
    const float y = 0.5f;

    std::vector<Vertex> lines;
    const float dim[3] = {color[0] * 0.55f, color[1] * 0.55f, color[2] * 0.55f};
    // mur
    addSegment(lines, kWallX, 0.0f, z - 0.9f, kWallX, 1.6f, z - 0.9f, dim);
    addSegment(lines, kWallX, 1.6f, z - 0.9f, kWallX, 1.6f, z + 0.9f, dim);
    addSegment(lines, kWallX, 1.6f, z + 0.9f, kWallX, 0.0f, z + 0.9f, dim);
    addSegment(lines, kWallX, 0.0f, z + 0.9f, kWallX, 0.0f, z - 0.9f, dim);
    // position d'équilibre
    const float tick[3] = {0.6f, 0.6f, 0.6f};
    addSegment(lines, 0.0f, 0.0f, z - 0.9f, 0.0f, 0.0f, z + 0.9f, tick);
    // masse
    addCube(lines, x, y, z, color);
    renderer.draw(Primitive::Lines, lines);

    // ressort : hélice de 10 spires entre le mur et la face gauche de la masse
    std::vector<Vertex> coil;
    const int turns = 10, perTurn = 14, n = turns * perTurn;
    const float x0 = kWallX, x1 = x - 0.5f;
    for (int i = 0; i <= n; ++i) {
        const float s = static_cast<float>(i) / n;
        const float a = static_cast<float>(kTwoPi) * turns * s;
        coil.push_back({x0 + (x1 - x0) * s, y + 0.32f * std::sin(a), z + 0.32f * std::cos(a), dim[0] + 0.15f, dim[1] + 0.15f, dim[2] + 0.15f});
    }
    renderer.draw(Primitive::LineStrip, coil);
}

}  // namespace

OscillatorModule::OscillatorModule() { reset(); }

const char* OscillatorModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Un ressort tire la masse vers sa position de repos : plus on l'étire, plus il tire fort. "
                   "Lâchée, la masse va et vient, comme un enfant sur une balançoire. Un ressort plus raide fait osciller "
                   "plus vite ; avec du frottement (air, huile), le mouvement s'éteint peu à peu.\n\n"
                   "La masse bleue est ce que fait vraiment la nature, la verte est le calcul de l'ordinateur.";
        case Level::Interesse:
            return "Une masse plus lourde oscille plus lentement ; un ressort plus raide, plus vite.\n\n"
                   "Si l'on pousse le ressort à intervalles réguliers, l'effet dépend du rythme : au rythme propre du "
                   "ressort, chaque poussée s'ajoute aux précédentes et le mouvement grandit sans cesse. C'est la "
                   "résonance (la balançoire, le verre qui se brise sous une note, le pont qui vibre).\n\n"
                   "Essayez « Régler sur le rythme du ressort ».";
        case Level::College:
            return "Loi de Hooke : F = k × x  (k en N/m).\n"
                   "Période : T = 2π √(m / k) ; fréquence : f = 1 / T (en Hz).\n"
                   "Énergies : Ec = ½ m v² ; Ep = ½ k x². Sans frottement, Ec + Ep reste constante.\n\n"
                   "Pourtant la méthode simple (Euler, orange) ajoute un peu d'énergie à chaque pas : l'oscillation grandit "
                   "alors que le ressort est parfait ! La méthode précise (RK4, vert) reste fidèle. "
                   "Diminuez dt pour réduire l'erreur d'Euler.";
        case Level::Lycee:
            return "2e loi de Newton : m a = −k x − c v + F(t)  (rappel du ressort, frottement, force extérieure).\n"
                   "Pulsation propre : ω0 = √(k / m) ; sans frottement x(t) = A cos(ω0 t + φ).\n\n"
                   "Taux d'amortissement ζ = c / (2 √(k m)) :\n"
                   "  ζ < 1 : oscillations amorties (pseudo-périodique)\n"
                   "  ζ = 1 : régime critique, retour le plus rapide à l'équilibre sans osciller\n"
                   "  ζ > 1 : retour lent sans oscillation (apériodique)\n\n"
                   "Forcé à la pulsation ω, l'amplitude est maximale près de ω0 : résonance. Le portrait de phase (x, v) "
                   "d'un oscillateur sans frottement est une ellipse : E = ½ m v² + ½ k x² est constante.";
        case Level::Etudiant:
            return "EDO : x'' + 2γ x' + ω0² x = (F0/m) cos ωt, avec γ = c / 2m. Solution = x_h + x_p.\n"
                   "  x_h = e^(−γt) (A cos ω_d t + B sin ω_d t), ω_d = √(ω0² − γ²)\n"
                   "  x_p = X cos(ωt − φ), X = (F0/m) / √((ω0² − ω²)² + (2γω)²)\n"
                   "Lagrangien L = ½ m x'² − ½ k x² ; Hamiltonien H = p²/2m + ½ k x².\n\n"
                   "Euler explicite multiplie l'énergie par (1 + ω0² dt²) à chaque pas : spirale divergente dans le portrait "
                   "de phase. Euler symplectique et Verlet gardent une énergie bornée (courbes fermées) ; RK4 dissipe très "
                   "légèrement.";
        case Level::Chercheur:
            return "Oscillateur libre : un pas à un pas s'écrit y_{n+1} = G y_n avec y = (x, v).\n"
                   "Euler explicite : G = [[1, dt], [−ω0² dt, 1]], valeurs propres 1 ± iω0 dt de module √(1 + ω0² dt²) > 1 : "
                   "instable pour tout dt.\n"
                   "Euler symplectique et Verlet : det G = 1 (application symplectique), stables si ω0 dt < 2. Verlet "
                   "conserve exactement un hamiltonien modifié H_mod = H + O(dt²) (analyse rétrograde) : énergie bornée, "
                   "sans dérive.\n"
                   "RK4 : |R(iω0 dt)|² = 1 − (ω0 dt)^6 / 72 + … : dissipation numérique d'ordre 6.\n\n"
                   "Résonance : facteur de qualité Q = 1 / (2ζ), pic en ω_res = √(ω0² − 2γ²), largeur à mi-puissance "
                   "Δω ≈ ω0 / Q. Décrément logarithmique δ = 2πζ / √(1 − ζ²).";
    }
    return "";
}

void OscillatorModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.5f;
    camera.target[2] = 3.0f;
    camera.distance = 21.0f;
    camera.yaw = 0.35f;
    camera.pitch = 0.5f;
}

// ------------------------------ simulation -----------------------------

void OscillatorModule::reset() {
    problem_.mass = mass_;
    problem_.stiffness = stiffness_;
    problem_.damping = damping_;
    problem_.forceAmplitude = forced_ ? forceAmp_ : 0.0;
    problem_.forceFrequency = forceOmega_;
    problem_.x0 = x0_;
    problem_.v0 = v0_;
    rhs_ = problem_.rhs();
    endTime_ = durationPeriods_ * problem_.period();

    for (Run& r : runs_) {
        r.y = problem_.initialState();
        r.x.clear();
        r.phase.clear();
        r.energy.clear();
        r.error.clear();
        r.diverged = false;
    }
    xExact_.clear();
    phaseExact_.clear();
    energyExact_.clear();
    solvers_.rk45().resetStats();

    time_ = 0.0;
    accumulator_ = 0.0;
    running_ = true;
    finished_ = false;
    convergenceDirty_ = true;
    sample();
}

void OscillatorModule::sample() {
    lastSampleTime_ = time_;
    double xe, ve;
    problem_.exact(time_, xe, ve);
    xExact_.add(time_, xe);
    phaseExact_.add(xe, ve);
    energyExact_.add(time_, std::max(problem_.energy(xe, ve), 1e-12));  // plancher : l'axe est logarithmique

    for (Run& r : runs_) {
        if (r.diverged) continue;
        r.x.add(time_, r.y[0]);
        r.phase.add(r.y[0], r.y[1]);
        r.energy.add(time_, std::max(problem_.energy(r.y[0], r.y[1]), 1e-12));
        r.error.add(time_, std::max(std::abs(r.y[0] - xe), 1e-12));
    }
}

void OscillatorModule::update(double frameSeconds) {
    if (!running_ || finished_) return;

    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;  // borne anti "spirale de la mort"
    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    const double sampleInterval = std::max(dt_, 1.0 / 240.0);

    int guard = 0;
    while (!finished_ && guard++ < 5000) {
        const double remaining = endTime_ - time_;
        const double h = std::min(dt_, remaining);  // dernier pas raccourci : on s'arrête pile à la fin
        if (accumulator_ < h) break;

        for (int i = 0; i < SolverSet::kCount; ++i) {
            Run& r = runs_[i];
            if (r.diverged) continue;
            advance(solvers_.solver(i), rhs_, time_, r.y, h);
            if (!std::isfinite(r.y[0]) || std::abs(r.y[0]) > 1e6) r.diverged = true;
        }
        time_ += h;
        accumulator_ -= h;
        finished_ = remaining <= dt_ * (1.0 + 1e-9);
        if (finished_ || time_ - lastSampleTime_ >= sampleInterval - 1e-12) sample();
    }
    if (finished_) running_ = false;
}

// Erreur finale de chaque solveur à pas fixe pour 7 pas différents, puis pente log-log (ordre mesuré).
void OscillatorModule::computeConvergence() {
    static const int kStepCounts[] = {200, 400, 800, 1600, 3200, 6400, 12800};  // Euler n'est asymptotique qu'à petit pas
    const double tEnd = 2.7 * problem_.period();  // pas un multiple de la période : voir oscillatorError

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        ConvergenceCurve& c = convergence_[i];
        c.dt.clear();
        c.err.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        for (int n : kStepCounts) {
            const double e = oscillatorError(problem_, *solver, n, tEnd);
            if (e > 1e-12) {  // sous ce seuil on mesure l'arrondi, pas la méthode
                c.dt.push_back(tEnd / n);
                c.err.push_back(e);
            }
        }

        c.slope = 0.0;
        const std::size_t m = c.dt.size();
        if (m >= 2) {
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (std::size_t k = 0; k < m; ++k) {
                const double x = std::log10(c.dt[k]), y = std::log10(c.err[k]);
                sx += x; sy += y; sxx += x * x; sxy += x * y;
            }
            const double dm = static_cast<double>(m);
            c.slope = (dm * sxy - sx * sy) / (dm * sxx - sx * sx);
        }
    }
    convergenceDirty_ = false;
}

// --------------------------------- UI ----------------------------------

void OscillatorModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    bool changed = false;

    pushSliderWidth();
    ImGui::SeparatorText("Ressort");
    changed |= sliderD(college ? "Raideur k (N/m)" : "Raideur du ressort", &stiffness_, 2.0, 60.0, college ? "%.1f" : "");
    if (atLeast(level, Level::Interesse))
        changed |= sliderD("Masse m (kg)", &mass_, 0.2, 10.0, "%.2f", ImGuiSliderFlags_Logarithmic);
    changed |= sliderD(college ? "Frottement c (kg/s)" : "Frottement", &damping_, 0.0, 16.0, college ? "%.2f" : "");
    if (atLeast(level, Level::Interesse))
        changed |= sliderD(college ? "Étirement initial (m)" : "Étirement initial", &x0_, -3.0, 3.0, college ? "%.2f" : "");
    if (atLeast(level, Level::Etudiant)) changed |= sliderD("Vitesse initiale (m/s)", &v0_, -6.0, 6.0, "%.2f");
    popSliderWidth();

    if (lycee) {
        ImGui::TextDisabled("Régime d'amortissement");
        struct Preset { const char* name; double zeta; };
        const Preset presets[] = {{"Libre", 0.0}, {"Sous-amorti", 0.15}, {"Critique", 1.0}, {"Sur-amorti", 2.5}};
        for (int i = 0; i < 4; ++i) {
            if (i > 0) ImGui::SameLine();
            if (ImGui::Button(presets[i].name)) {
                damping_ = presets[i].zeta * 2.0 * std::sqrt(stiffness_ * mass_);
                changed = true;
            }
        }
    }

    if (atLeast(level, Level::Interesse)) {
        ImGui::SeparatorText(college ? "Force extérieure" : "Pousser le ressort");
        changed |= ImGui::Checkbox(college ? "Force périodique F0 cos(ω t)" : "Pousser régulièrement", &forced_);
        if (forced_) {
            pushSliderWidth();
            if (level == Level::Interesse) {
                changed |= sliderD("Rythme des poussées", &forceOmega_, 0.3, 25.0, "", ImGuiSliderFlags_Logarithmic);
            } else if (level == Level::College || level == Level::Lycee) {
                double f = forceOmega_ / kTwoPi;
                if (sliderD("Fréquence f (Hz)", &f, 0.05, 4.0, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                    forceOmega_ = kTwoPi * f;
                    changed = true;
                }
            } else {
                changed |= sliderD("Pulsation ω (rad/s)", &forceOmega_, 0.3, 25.0, "%.3f", ImGuiSliderFlags_Logarithmic);
            }
            if (college) changed |= sliderD("Force F0 (N)", &forceAmp_, 0.1, 20.0, "%.1f", ImGuiSliderFlags_Logarithmic);
            popSliderWidth();
            if (ImGui::Button(college ? "Régler sur ω0 (résonance)" : "Régler sur le rythme du ressort")) {
                forceOmega_ = problem_.omega0();
                changed = true;
            }
        }
    }

    if (college) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Etudiant))
            changed |= sliderD("Durée (périodes)", &durationPeriods_, 2.0, 100.0, "%.0f");
        if (atLeast(level, Level::Chercheur))
            changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
        popSliderWidth();
    }

    solvers_.drawToggles(level);
    if (changed) reset();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !finished_ ? "Pause" : "Lecture")) {
        if (finished_) reset();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) reset();
    pushSliderWidth();
    if (atLeast(level, Level::Interesse))
        sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        double xe, ve;
        problem_.exact(time_, xe, ve);
        ImGui::Text("t = %.2f s / %.2f s", time_, endTime_);
        ImGui::Text("x = %.3f m   v = %.3f m/s", xe, ve);

        ImGui::SeparatorText("Grandeurs");
        const double w0 = problem_.omega0(), zeta = problem_.zeta(), gamma = problem_.gamma();
        ImGui::Text("T = 2π√(m/k) = %.3f s", problem_.period());
        ImGui::Text("f = 1/T = %.3f Hz", 1.0 / problem_.period());
        if (lycee) {
            const char* names[] = {"sous-amorti (pseudo-périodique)", "critique", "sur-amorti (apériodique)"};
            ImGui::Text("ω0 = %.3f rad/s", w0);
            ImGui::Text("ζ = %.3f : %s", zeta, names[static_cast<int>(problem_.regime())]);
            if (forced_) ImGui::Text("Amplitude permanente X = %.3f m", problem_.steadyStateAmplitude(forceOmega_));
        }
        if (atLeast(level, Level::Etudiant)) {
            ImGui::Text("γ = c/2m = %.3f 1/s", gamma);
            if (problem_.regime() == DampingRegime::Underdamped) ImGui::Text("ω_d = %.3f rad/s", problem_.dampedOmega());
        }
        if (atLeast(level, Level::Chercheur)) {
            if (zeta > 0.0) ImGui::Text("Q = 1/(2ζ) = %.2f", 1.0 / (2.0 * zeta));
            if (zeta > 0.0 && zeta < 1.0) ImGui::Text("δ = 2πζ/√(1−ζ²) = %.3f", kTwoPi * zeta / std::sqrt(1.0 - zeta * zeta));
            if (forced_ && problem_.resonanceOmega() > 0.0)
                ImGui::Text("ω_res = %.3f rad/s  X_res = %.3f m", problem_.resonanceOmega(),
                            problem_.steadyStateAmplitude(problem_.resonanceOmega()));
        }
    }
}

void OscillatorModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant le mouvement, l'énergie passe sans cesse du ressort (étiré) à la masse (qui va vite). "
                           "Sans frottement, leur total reste le même.");
        return;
    }

    double xe, ve;
    problem_.exact(time_, xe, ve);
    const double exactEnergy = problem_.energy(xe, ve);
    const double e0 = problem_.energy(x0_, v0_);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool researcher = atLeast(level, Level::Chercheur);

    // Chercheur : dE et dE/E0 remplacent E (redondant) pour que le tableau tienne dans le panneau.
    std::vector<std::string> headers;
    if (!researcher) headers.push_back("E (J)");
    if (lycee) headers.push_back("dE (J)");
    if (researcher) headers.push_back("dE/E0");
    headers.push_back("dx (m)");

    std::vector<TableRow> rows;
    {
        TableRow ref{"Exacte", kBlue, {}};
        if (!researcher) ref.cells.push_back(strf("%.3f", exactEnergy));
        if (lycee) ref.cells.push_back("-");
        if (researcher) ref.cells.push_back("-");
        ref.cells.push_back("-");
        rows.push_back(ref);
    }
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(level, i)) continue;
        const Run& r = runs_[i];
        TableRow row{solvers_.shortLabel(level, i), solvers_.color(i), {}};
        if (r.diverged) {
            row.cells.assign(headers.size(), "-");
            row.cells[0] = "diverge";
        } else {
            const double energy = problem_.energy(r.y[0], r.y[1]);
            if (!researcher) row.cells.push_back(strf("%.3f", energy));
            if (lycee) row.cells.push_back(strf("%+.1e", energy - exactEnergy));
            if (researcher) row.cells.push_back(strf("%+.1e", (energy - exactEnergy) / e0));
            row.cells.push_back(strf("%.1e", std::abs(r.y[0] - xe)));
        }
        rows.push_back(row);
    }
    drawResultTable("erreurs", headers, rows);

    if (lycee) {
        ImGui::SeparatorText("Solution exacte");
        ImGui::Text("Ec = ½ m v² = %.3f J", 0.5 * mass_ * ve * ve);
        ImGui::Text("Ep = ½ k x² = %.3f J", 0.5 * stiffness_ * xe * xe);
        if (damping_ > 0.0 || forced_) ImGui::TextDisabled("Frottement ou force extérieure : E varie.");
        else ImGui::TextDisabled("Sans frottement : E constante.");
    }
    if (researcher) {
        ImGui::SeparatorText("RK45");
        ImGui::Text("%d pas acceptés, %d rejetés", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps());
        ImGui::Text("%d évaluations de f", solvers_.rk45().evaluations());
        if (damping_ == 0.0 && !forced_ && !runs_[SolverSet::kEuler].diverged) {
            const double w0 = problem_.omega0();
            const double steps = std::round(time_ / dt_);
            const double predicted = std::pow(1.0 + w0 * w0 * dt_ * dt_, steps);  // chaque pas multiplie E par 1 + w0² dt²
            const double measured = problem_.energy(runs_[SolverSet::kEuler].y[0], runs_[SolverSet::kEuler].y[1]) / e0;
            ImGui::SeparatorText("Euler explicite");
            ImGui::Text("E/E0 mesuré = %.4f", measured);
            ImGui::TextDisabled("(1 + ω0² dt²)^n = %.4f", predicted);
        }
    }
}

void OscillatorModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showEnergy = atLeast(level, Level::College);
    const bool showError = atLeast(level, Level::Lycee);
    const int cols = 1 + (showEnergy ? 1 : 0) + (showError ? 1 : 0);

    auto spec = [](const float c[3], int offset) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Offset, offset);
    };
    auto plotRuns = [&](Series Run::*member) {
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i) || (runs_[i].*member).size() == 0) continue;
            const Series& s = runs_[i].*member;
            ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                             spec(solvers_.color(i), s.offset));
        }
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Élongation x(t)")) {
            ImPlot::SetupAxes("t (s)", "x (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (xExact_.size() > 0)
                ImPlot::PlotLine("Exacte", xExact_.x.data(), xExact_.y.data(), xExact_.size(), spec(kBlue, xExact_.offset));
            plotRuns(&Run::x);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie E(t)")) {
            ImPlot::SetupAxes("t (s)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (atLeast(level, Level::Lycee)) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (energyExact_.size() > 0)
                ImPlot::PlotLine("Exacte", energyExact_.x.data(), energyExact_.y.data(), energyExact_.size(),
                                 spec(kBlue, energyExact_.offset));
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Erreur |x - x_exact|")) {
            ImPlot::SetupAxes("t (s)", "m", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void OscillatorModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showResonance = forced_;
    const bool showConvergence = atLeast(level, Level::Etudiant);
    if (showConvergence && convergenceDirty_) computeConvergence();
    const int cols = 1 + (showResonance ? 1 : 0) + (showConvergence ? 1 : 0);

    auto spec = [](const float c[3], int offset) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Offset, offset);
    };

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        // Portrait de phase : une ellipse fermée pour un oscillateur idéal bien intégré.
        if (ImPlot::BeginPlot("Portrait de phase (x, v)")) {
            ImPlot::SetupAxes("x (m)", "v (m/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (phaseExact_.size() > 0)
                ImPlot::PlotLine("Exacte", phaseExact_.x.data(), phaseExact_.y.data(), phaseExact_.size(),
                                 spec(kBlue, phaseExact_.offset));
            for (int i = 0; i < SolverSet::kCount; ++i) {
                if (!solvers_.isShown(level, i) || runs_[i].phase.size() == 0) continue;
                const Series& s = runs_[i].phase;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                                 spec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }

        // Courbe de résonance : amplitude du régime permanent rapportée à la déformation statique F0/k.
        if (showResonance && ImPlot::BeginPlot("Résonance : X(ω) / (F0/k)")) {
            ImPlot::SetupAxes("ω (rad/s)", "amplification", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            const double w0 = problem_.omega0();
            const double staticDeflection = std::max(forceAmp_ / stiffness_, 1e-12);
            auto curve = [&](double zeta, std::vector<double>& xs, std::vector<double>& ys) {
                OscillatorProblem p = problem_;
                p.forceAmplitude = forceAmp_;
                p.damping = zeta * 2.0 * std::sqrt(stiffness_ * mass_);
                xs.clear();
                ys.clear();
                for (int i = 0; i <= 300; ++i) {
                    const double w = 2.5 * w0 * i / 300.0;
                    xs.push_back(w);
                    ys.push_back(std::min(p.steadyStateAmplitude(w) / staticDeflection, 50.0));  // plafonné (pic infini)
                }
            };
            std::vector<double> xs, ys;
            const float grey[3] = {0.5f, 0.5f, 0.55f};
            for (double z : {0.1, 0.3, 1.0}) {
                curve(z, xs, ys);
                char name[32];
                std::snprintf(name, sizeof(name), "ζ = %.1f", z);
                ImPlot::PlotLine(name, xs.data(), ys.data(), static_cast<int>(xs.size()), spec(grey, 0));
            }
            curve(problem_.zeta(), xs, ys);
            const float accent[3] = {0.27f, 0.70f, 0.88f};
            char name[32];
            std::snprintf(name, sizeof(name), "ζ = %.2f (actuel)", problem_.zeta());
            ImPlot::PlotLine(name, xs.data(), ys.data(), static_cast<int>(xs.size()), spec(accent, 0));
            const double w = forceOmega_;
            ImPlot::PlotInfLines("ω appliquée", &w, 1, spec(kBlue, 0));
            ImPlot::EndPlot();
        }

        if (showConvergence && ImPlot::BeginPlot("Convergence : erreur à t = 2,7 T")) {
            ImPlot::SetupAxes("dt (s)", "erreur (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kFixedStep; ++i) {
                const ConvergenceCurve& c = convergence_[i];
                if (!solvers_.show(i) || c.dt.size() < 2) continue;
                char name[96];
                std::snprintf(name, sizeof(name), "%s (pente %.2f)", solvers_.solver(i).name(), c.slope);
                const ImPlotSpec markers(ImPlotProp_LineColor, toImVec4(solvers_.color(i)), ImPlotProp_LineWeight, 2.0f,
                                         ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f);
                ImPlot::PlotLine(name, c.dt.data(), c.err.data(), static_cast<int>(c.dt.size()), markers);
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

void OscillatorModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    // Couloir 0 : la solution exacte ; puis un couloir par méthode affichée.
    double xe, ve;
    problem_.exact(time_, xe, ve);
    drawLane(renderer, 0, xe, kBlue);
    int lane = 1;
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        // Une méthode qui a divergé reste figée sur le bord du cadre, en rouge.
        const bool diverged = runs_[i].diverged;
        const float red[3] = {0.9f, 0.2f, 0.2f};
        drawLane(renderer, lane++, diverged ? 12.0 : runs_[i].y[0], diverged ? red : solvers_.color(i));
    }
}

}  // namespace pl
