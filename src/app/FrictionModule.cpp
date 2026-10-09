#include "FrictionModule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <implot.h>

namespace pl {
namespace {

constexpr double kPi = constants::pi;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr int kStepBudget = 400;       // appels au solveur par pas de calcul : au-delà, un solveur adaptatif est déclaré bloqué
constexpr float kBlockHalfLength = 0.65f, kBlockHeight = 0.9f;  // taille du bloc en unités de scène
constexpr float kLaneSpacing = 1.4f;
const float kWeightColor[3] = {0.60f, 0.70f, 1.00f};
const float kNormalColor[3] = {0.40f, 1.00f, 0.50f};
const float kFrictionColor[3] = {1.00f, 0.40f, 0.35f};
const float kGrey[3] = {0.55f, 0.55f, 0.6f};

struct Preset {
    const char* name;
    double angleDeg, muStatic, muKinetic, launch;
};
const Preset kPresets[] = {{"Pente douce : il colle", 17.2, 0.5, 0.4, 6.0},
                           {"Pente forte : il redescend", 34.4, 0.5, 0.4, 5.0},
                           {"Descente freinée", 11.5, 0.6, 0.5, -3.0},
                           {"Piège : adhérence", 26.0, 0.5, 0.4, 6.0}};

void arrow(Renderer& renderer, float x, float y, float z, float dx, float dy, const float* c, float uiScale) {
    renderer.draw(Primitive::Lines, {{x, y, z, c[0], c[1], c[2]}, {x + dx, y + dy, z, c[0], c[1], c[2]}});
    renderer.draw(Primitive::Points, {{x + dx, y + dy, z, c[0], c[1], c[2]}}, 7.0f * uiScale);
}

}  // namespace

FrictionModule::FrictionModule() {
    reset();
}

const char* FrictionModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Un bloc est lancé vers le haut d'une pente. Le frottement le freine : il s'arrête, puis soit il reste "
                   "collé (pente douce), soit il redescend (pente raide).\n\n"
                   "Le bloc bleu suit les vraies lois, le vert est le calcul de l'ordinateur. Regardez l'instant où le bloc "
                   "s'arrête : l'ordinateur a du mal avec les arrêts et les départs, et peut laisser glisser un bloc qui "
                   "devrait rester collé.\n\n"
                   "Cochez « Calcul amélioré » pour voir la différence.";
        case Level::Interesse:
            return "Il y a deux frottements : quand le bloc glisse, c'est le frottement de glissement ; quand il est arrêté, "
                   "c'est l'adhérence, plus forte, qui peut le retenir. Il reste collé tant que la pente n'est pas trop raide.\n\n"
                   "Essayez les quatre cas. Le dernier est un piège : la pente est trop raide pour que le bloc glisse "
                   "« tout seul » en bas, mais l'adhérence le retient. Un calcul qui ne connaît que le frottement de "
                   "glissement le laisse partir.";
        case Level::College:
            return "Force de frottement : F = μ × N, avec N = m g cos θ (réaction du plan). Glissement vers le haut : "
                   "décélération a = g (sin θ + μd cos θ). Le bloc s'arrête après t = v0 / a, à la distance v0² / (2a).\n"
                   "Il reste collé si tan θ ≤ μs (coefficient d'adhérence) ; sinon il redescend avec "
                   "a = g (sin θ − μd cos θ).\n\n"
                   "Dans « Invariants », comparez la position et la vitesse finales : le calcul naïf (orange, vert) n'a jamais "
                   "une vitesse exactement nulle.";
        case Level::Lycee:
            return "Énergie par unité de masse : E = ½ v² + g s sin θ. Elle diminue du travail du frottement, "
                   "μd g cos θ × (chemin parcouru) : la chaleur. Flèches sur le bloc exact : poids (bleu), réaction normale "
                   "(verte), frottement (rouge). À l'arrêt, le frottement d'adhérence égale la composante du poids g sin θ "
                   "(rouge vers le haut), tant qu'elle ne dépasse pas μs g cos θ.\n\n"
                   "Dans « Analyse » : la loi de Coulomb. Le frottement change de signe quand la vitesse s'annule : "
                   "c'est une force discontinue, et c'est ce qui met les calculs en difficulté.";
        case Level::Etudiant:
            return "s'' = −g sin θ − μd g cos θ sgn(s') − k s' : second membre discontinu en s' = 0 (équation à second membre "
                   "discontinu, de Filippov). Dans chaque phase le sens est fixé et l'EDO est linéaire : "
                   "v' = A − k v, A = −g (sin θ + dir μd cos θ), de solution exacte (voir la passation).\n"
                   "Les théorèmes de convergence supposent un second membre régulier : ici RK4 retombe à l'ordre 1 (graphe de "
                   "convergence, modèle naïf), Verlet aussi. Le modèle naïf ignore μs : si μd < tan θ ≤ μs, il converge vers "
                   "une mauvaise solution (le plateau de l'erreur).\n"
                   "RK45 y échoue autrement : sur la surface v = 0 avec adhérence, la solution « glisse » et l'erreur locale reste "
                   "d'ordre h quel que soit h : le pas s'effondre, on le plafonne et on le déclare bloqué.";
        case Level::Chercheur:
            return "Trois traitements. (1) Événement + machine à états (modèle correct ici) : l'arrêt est localisé par bissection "
                   "sur le pas, la vitesse est imposée nulle, puis adhérence si tan θ ≤ μs, sinon changement de sens ; les ordres "
                   "des solveurs sont retrouvés. (2) Régularisation sgn(v) → tanh(v/ε) : continue mais raide ; l'erreur vient "
                   "du modèle (ε), pas du pas, et le bloc rampe (fluage) au lieu d'adhérer. (3) Complémentarité : le frottement de "
                   "Coulomb est une inclusion différentielle, que les schémas « time-stepping » (Moreau-Jean) résolvent par un "
                   "problème de complémentarité à chaque pas : c'est la voie des moteurs multicorps (M5c, M6).\n"
                   "Sans résistance (k = 0) chaque phase est un polynôme de degré 2 : Verlet et RK4 sont exacts, rien à mesurer.";
    }
    return "";
}

void FrictionModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 4.0f;
    camera.target[2] = -3.0f;
    camera.distance = 26.0f;
    camera.yaw = 0.35f;
    camera.pitch = 0.3f;
}

// ------------------------------ simulation -----------------------------

void FrictionModule::applyPreset(double angleDeg, double muS, double muK, double v0) {
    angleDeg_ = angleDeg;
    muStatic_ = muS;
    muKinetic_ = muK;
    launch_ = v0;
    reset();
}

void FrictionModule::reset() {
    problem_.angle = angleDeg_ * kPi / 180.0;
    problem_.muStatic = muStatic_;
    problem_.muKinetic = muKinetic_;
    problem_.drag = drag_;
    problem_.s0 = 0.0;
    problem_.v0 = launch_;
    problem_.regularization = epsilon_;
    rhs_ = problem_.rhs(model_);

    // Rampe : couvre la trajectoire exacte, plus de la place pour qu'un calcul fautif glisse vers le bas.
    double sMin = 0.0, sMax = 0.0;
    for (int i = 0; i <= 200; ++i) {
        const double s = problem_.exact(durationSec_ * i / 200.0).s;
        sMin = std::min(sMin, s);
        sMax = std::max(sMax, s);
    }
    lo_ = sMin - 8.0;
    hi_ = sMax + 3.0;
    scale_ = static_cast<float>(std::clamp(16.0 / (hi_ - lo_), 0.2, 2.5));

    for (Run& r : runs_) {
        r.y = problem_.initialState();
        r.event = model_ == InclineModel::EventDriven ? std::make_unique<InclineRun>(problem_) : nullptr;
        r.position.clear();
        r.velocity.clear();
        r.error.clear();
        r.energy.clear();
        r.stalled = false;
        r.stallTime = 0.0;
        r.stopTime = kNaN;
    }
    exactPosition_.clear();
    exactVelocity_.clear();
    exactEnergy_.clear();
    exactFriction_.clear();
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void FrictionModule::sample() {
    const double t = clock_.time;
    const InclineState ex = problem_.exact(t);
    exactPosition_.add(t, ex.s);
    exactVelocity_.add(t, ex.v);
    exactEnergy_.add(t, problem_.energy(ex.s, ex.v));
    const double gc = problem_.gravity * std::cos(problem_.angle), gs = problem_.gravity * std::sin(problem_.angle);
    exactFriction_.add(t, ex.stuck ? gs : (ex.v > 0.0 ? -problem_.muKinetic * gc : problem_.muKinetic * gc));

    for (Run& r : runs_) {
        if (r.stalled) continue;
        r.position.add(t, r.y[0]);
        r.velocity.add(t, r.y[1]);
        r.error.add(t, std::max(std::abs(r.y[0] - ex.s), 1e-12));
        r.energy.add(t, problem_.energy(r.y[0], r.y[1]));
    }
}

void FrictionModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;

    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    clock_.advance(
        frameSeconds, timeScale_, dt_, durationSec_,
        [&](double h) {
            for (int i = 0; i < SolverSet::kCount; ++i) {
                Run& r = runs_[i];
                if (r.stalled) continue;
                const double before = r.y[1], t = clock_.time;
                if (r.event) {
                    r.event->advance(solvers_.solver(i), h);
                    r.y = {r.event->state().s, r.event->state().v};
                    r.stopTime = r.event->stopTime();
                } else {
                    // Budget de pas : RK45 peut s'effondrer sur la discontinuité (voir testAdvanceBudget).
                    const double done = advance(solvers_.solver(i), rhs_, t, r.y, h, kStepBudget);
                    if (done < h * (1.0 - 1e-9)) {
                        r.stalled = true;
                        r.stallTime = t + done;
                        continue;
                    }
                    // Premier arrêt : la vitesse change de signe (ou s'annule) pendant ce pas ; instant interpolé linéairement.
                    if (std::isnan(r.stopTime) && before != 0.0 && before * r.y[1] <= 0.0)
                        r.stopTime = t + h * std::abs(before) / (std::abs(before) + std::abs(r.y[1]));
                }
            }
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur de position à la fin de la simulation pour différents pas, pour les 4 solveurs à pas fixe et le modèle courant.
void FrictionModule::computeConvergence() {
    static const int kStepCounts[] = {25, 50, 100, 200, 400, 800, 1600};
    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            const double e = inclineError(problem_, model_, *solver, n, durationSec_);
            if (!std::isfinite(e) || e <= 1e-12) continue;  // sous ce seuil : arrondi
            const double dt = durationSec_ / n;
            c.x.push_back(dt);
            c.y.push_back(e);
            const double lx = std::log10(dt), ly = std::log10(e);
            sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
            ++m;
        }
        c.slope = kNaN;
        if (m >= 3) {
            const double dm = static_cast<double>(m);
            c.slope = (dm * sxy - sx * sy) / (dm * sxx - sx * sx);
        }
    }
    convergenceDirty_ = false;
}

// Ce que la théorie prévoit pour le bloc, en mots.
std::string FrictionModule::verdict() const {
    const bool holds = problem_.holds();
    if (problem_.v0 == 0.0) return holds ? "Au repos : il reste collé." : "Au repos : il glisse vers le bas.";
    const double A = problem_.stageAcceleration(problem_.v0 > 0.0 ? 1 : -1);
    const bool stops = problem_.v0 * A < 0.0;
    if (!stops) return "Il descend en accélérant, sans jamais s'arrêter.";
    return holds ? "Il s'arrête, puis reste collé." : "Il s'arrête, puis redescend.";
}

// --------------------------------- UI ----------------------------------

void FrictionModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool etudiant = atLeast(level, Level::Etudiant);
    bool changed = false;

    if (interesse) {
        ImGui::SeparatorText("Cas types");
        for (const Preset& p : kPresets) {
            if (ImGui::Button(p.name)) {
                applyPreset(p.angleDeg, p.muStatic, p.muKinetic, p.launch);
                return;  // état réinitialisé : on redessine au prochain tour
            }
        }
    }

    ImGui::SeparatorText("Pente et frottement");
    pushSliderWidth();
    changed |= sliderD(college ? "Pente θ (deg)" : "Pente", &angleDeg_, 0.0, 60.0, college ? "%.1f" : "");
    bool staticChanged = false, kineticChanged = false;
    if (interesse) {
        staticChanged = sliderD(college ? "Adhérence μs" : "Adhérence (au repos)", &muStatic_, 0.0, 1.2, "%.2f");
        kineticChanged = sliderD(college ? "Glissement μd" : "Glissement", &muKinetic_, 0.0, 1.2, "%.2f");
    }
    changed |= staticChanged || kineticChanged;
    if (kineticChanged && muKinetic_ > muStatic_) muStatic_ = muKinetic_;   // toujours μs >= μd
    if (staticChanged && muStatic_ < muKinetic_) muKinetic_ = muStatic_;
    changed |= sliderD(college ? "Vitesse de lancement (m/s)" : "Vitesse de lancement", &launch_, -8.0, 10.0, college ? "%.1f" : "");
    if (etudiant) changed |= sliderD("Résistance k (1/s)", &drag_, 0.0, 2.0, "%.2f");
    popSliderWidth();

    // Modèle de calcul : au niveau 1-2 une simple case, ensuite les trois modèles.
    ImGui::SeparatorText("Calcul");
    if (!college) {
        bool improved = model_ == InclineModel::EventDriven;
        if (ImGui::Checkbox("Calcul amélioré (détecte l'arrêt)", &improved)) {
            model_ = improved ? InclineModel::EventDriven : InclineModel::Naive;
            changed = true;
        }
    } else {
        int m = static_cast<int>(model_);
        bool modelChanged = ImGui::RadioButton("Naïf : sgn(v) dans l'équation", &m, static_cast<int>(InclineModel::Naive));
        modelChanged |= ImGui::RadioButton("Régularisé : tanh(v/ε)", &m, static_cast<int>(InclineModel::Regularized));
        modelChanged |= ImGui::RadioButton("Événement + adhérence", &m, static_cast<int>(InclineModel::EventDriven));
        if (modelChanged) {
            model_ = static_cast<InclineModel>(m);
            changed = true;
        }
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Lycee)) changed |= sliderD("Durée (s)", &durationSec_, 2.0, 15.0, "%.1f");
        if (etudiant && model_ == InclineModel::Regularized)
            changed |= sliderD("Régularisation ε (m/s)", &epsilon_, 0.002, 0.5, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Chercheur))
            changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-10, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
        popSliderWidth();
    }

    solvers_.drawToggles(level);
    if (changed) reset();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !clock_.finished ? "Pause" : "Lecture")) {
        if (clock_.finished) reset();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) reset();
    pushSliderWidth();
    if (interesse) sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) ImGui::Text("t = %.2f / %.1f s", clock_.time, durationSec_);
    if (interesse) wrapped("Théorie : " + verdict());

    for (int i = 0; i < SolverSet::kCount; ++i)
        if (solvers_.isShown(level, i) && runs_[i].stalled)
            wrapped(strf("%s : BLOQUÉ à t = %.2f s (pas devenu trop petit : la dynamique est discontinue).",
                         solvers_.label(level, i).c_str(), runs_[i].stallTime), true);
}

void FrictionModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Le frottement transforme le mouvement en chaleur : l'énergie du bloc diminue. Quand il s'arrête, "
                           "c'est l'adhérence qui décide : le bloc reste collé ou repart.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const double g = problem_.gravity, th = problem_.angle;

    ImGui::SeparatorText("Théorie");
    wrapped(strf("tan θ = %.3f ; μs = %.2f ; μd = %.2f", std::tan(th), problem_.muStatic, problem_.muKinetic));
    wrapped(strf("Adhérence au repos : %s (tan θ %s μs)", problem_.holds() ? "oui" : "non", problem_.holds() ? "≤" : ">"));
    if (problem_.v0 != 0.0) {
        const double A = problem_.stageAcceleration(problem_.v0 > 0.0 ? 1 : -1);
        wrapped(strf("Accélération pendant le lancer : %.2f m/s²", A));
        if (std::isfinite(problem_.firstStopTime()))
            wrapped(strf("Arrêt à t = %.3f s, à s = %.3f m", problem_.firstStopTime(), problem_.exact(problem_.firstStopTime()).s));
    }
    wrapped("→ " + verdict());
    if (lycee)
        wrapped(strf("Frottement : glissement μd g cos θ = %.2f m/s² ; à retenir g sin θ = %.2f ; maximum d'adhérence μs g cos θ = %.2f",
                     problem_.muKinetic * g * std::cos(th), g * std::sin(th), problem_.muStatic * g * std::cos(th)));

    const double tStopExact = problem_.v0 != 0.0 ? problem_.firstStopTime() : kNaN;
    const InclineState ex = problem_.exact(clock_.time);
    const double exEnergy = problem_.energy(ex.s, ex.v);

    std::vector<std::string> headersA = {"s (m)", "v (m/s)", "écart (m)"};
    std::vector<std::string> headersB = {"t_arrêt", "Δt_arr.", "E−E_ex"};
    std::vector<TableRow> rowsA, rowsB;
    rowsA.push_back({"Exacte", kBlue, {strf("%+.3f", ex.s), strf("%+.3f", ex.v), "-"}});
    rowsB.push_back({"Exacte", kBlue, {std::isfinite(tStopExact) ? strf("%.3f", tStopExact) : "-", "-", "-"}});

    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(level, i)) continue;
        const Run& r = runs_[i];
        TableRow a{solvers_.shortLabel(level, i), solvers_.color(i), {}};
        TableRow b = a;
        if (r.stalled) {
            a.cells = {"bloqué", "-", "-"};
            b.cells = {"-", "-", "-"};
        } else {
            a.cells = {strf("%+.3f", r.y[0]), strf("%+.3f", r.y[1]), strf("%.1e", std::abs(r.y[0] - ex.s))};
            b.cells.push_back(std::isfinite(r.stopTime) ? strf("%.3f", r.stopTime) : "-");
            b.cells.push_back(std::isfinite(r.stopTime) && std::isfinite(tStopExact) ? strf("%.1e", std::abs(r.stopTime - tStopExact)) : "-");
            b.cells.push_back(strf("%+.1e", problem_.energy(r.y[0], r.y[1]) - exEnergy));
        }
        rowsA.push_back(a);
        rowsB.push_back(b);
    }

    ImGui::SeparatorText("Comparaison à l'instant courant");
    drawResultTable("position", headersA, rowsA);
    if (lycee) {
        ImGui::SeparatorText("Arrêt et énergie");
        drawResultTable("arret", headersB, rowsB);
    }
    wrapped("écart : |s − s exacte|. Une vitesse non nulle alors que le bloc devrait être collé est une erreur de calcul.", true);
    if (lycee) wrapped("t_arrêt : premier instant où v s'annule ; E−E_ex : énergie (J/kg) par rapport à l'exacte.", true);
    if (etudiant && model_ == InclineModel::Regularized)
        wrapped(strf("Régularisation ε = %.3f m/s : la vitesse ne vaut jamais exactement 0.", epsilon_), true);
}

void FrictionModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showError = atLeast(level, Level::College);
    const bool showEnergy = atLeast(level, Level::Lycee);
    const int cols = 2 + (showError ? 1 : 0) + (showEnergy ? 1 : 0);

    auto plotRuns = [&](Series Run::*member) {
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i) || (runs_[i].*member).size() == 0) continue;
            const Series& s = runs_[i].*member;
            ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(solvers_.color(i), s.offset));
        }
    };
    auto plotExact = [&](const Series& s) {
        if (s.size() > 0) ImPlot::PlotLine("Exacte", s.x.data(), s.y.data(), s.size(), lineSpec(kBlue, s.offset));
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Position s(t) le long de la pente")) {
            ImPlot::SetupAxes("t (s)", "s (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactPosition_);
            plotRuns(&Run::position);
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Vitesse v(t)")) {
            ImPlot::SetupAxes("t (s)", "v (m/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactVelocity_);
            plotRuns(&Run::velocity);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Écart de position")) {
            ImPlot::SetupAxes("t (s)", "|s − s exacte| (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie ½v² + g s sin θ")) {
            ImPlot::SetupAxes("t (s)", "J/kg", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactEnergy_);
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void FrictionModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool researcher = atLeast(level, Level::Chercheur);
    if (etudiant && convergenceDirty_) computeConvergence();
    const int cols = 2 + (etudiant ? 1 : 0) + (researcher ? 1 : 0);

    auto markerSpec = [](const float* c) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle,
                          ImPlotProp_MarkerSize, 4.0f);
    };

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        // Loi de Coulomb : f(v) = -mu_d g cos(theta) sgn(v), et tout l'intervalle [-mu_s, mu_s] g cos(theta) à v = 0.
        // Le bloc reste collé si la force à retenir g sin(theta) tombe dans cet intervalle.
        if (ImPlot::BeginPlot("Loi de Coulomb : frottement f(v)")) {
            const double gc = problem_.gravity * std::cos(problem_.angle), gs = problem_.gravity * std::sin(problem_.angle);
            const double vm = 6.0, kin = problem_.muKinetic * gc, stat = problem_.muStatic * gc;
            ImPlot::SetupAxes("v (m/s)", "f (m/s²)", ImPlotAxisFlags_None, ImPlotAxisFlags_None);
            ImPlot::SetupAxisLimits(ImAxis_X1, -vm, vm, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -1.5 * std::max(stat, gs), 1.5 * std::max(stat, gs), ImPlotCond_Always);
            const double xn[2] = {-vm, 0.0}, yn[2] = {kin, kin}, xp[2] = {0.0, vm}, yp[2] = {-kin, -kin};
            ImPlot::PlotLine("glissement", xn, yn, 2, lineSpec(kFrictionColor, 0));
            ImPlot::PlotLine("glissement ", xp, yp, 2, lineSpec(kFrictionColor, 0));
            const double xs[2] = {0.0, 0.0}, ys[2] = {-stat, stat};
            ImPlot::PlotLine("adhérence (intervalle)", xs, ys, 2, lineSpec(kNormalColor, 0));
            const double hold = gs;
            ImPlot::PlotInfLines("à retenir : g sin θ", &hold, 1,
                                 ImPlotSpec(ImPlotProp_LineColor, toImVec4(kWeightColor), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            if (etudiant) {  // régularisation : tanh continue, mais jamais l'intervalle d'adhérence
                std::vector<double> xr, yr;
                for (int k = 0; k <= 200; ++k) {
                    const double v = -vm + 2.0 * vm * k / 200.0;
                    xr.push_back(v);
                    yr.push_back(-kin * std::tanh(v / problem_.regularization));
                }
                ImPlot::PlotLine("régularisé", xr.data(), yr.data(), static_cast<int>(xr.size()), lineSpec(kGrey, 0));
            }
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Frottement subi par le bloc exact")) {
            ImPlot::SetupAxes("t (s)", "f (m/s²)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (exactFriction_.size() > 0)
                ImPlot::PlotLine("f(t)", exactFriction_.x.data(), exactFriction_.y.data(), exactFriction_.size(), lineSpec(kFrictionColor, exactFriction_.offset));
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot("Convergence : erreur finale")) {
            ImPlot::SetupAxes("dt (s)", "|s − s exacte| (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kFixedStep; ++i) {
                const Curve& c = convergence_[i];
                if (!solvers_.show(i) || c.x.size() < 2) continue;
                const std::string name = std::isfinite(c.slope) ? strf("%s (pente %.2f)", solvers_.solver(i).name(), c.slope)
                                                                  : strf("%s", solvers_.solver(i).name());
                ImPlot::PlotLine(name.c_str(), c.x.data(), c.y.data(), static_cast<int>(c.x.size()), markerSpec(solvers_.color(i)));
            }
            ImPlot::EndPlot();
        }

        // Vitesse résiduelle : un bloc qui devrait être collé garde une vitesse non nulle (naïf, régularisé) ; avec l'événement, exactement 0.
        if (researcher && ImPlot::BeginPlot("|v| (échelle log)")) {
            ImPlot::SetupAxes("t (s)", "|v| (m/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kCount; ++i) {
                const Series& s = runs_[i].velocity;
                if (!solvers_.isShown(level, i) || s.size() == 0) continue;
                std::vector<double> mag(s.y.size());
                for (std::size_t k = 0; k < mag.size(); ++k) mag[k] = std::max(std::abs(s.y[k]), 1e-18);
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), mag.data(), s.size(), lineSpec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

// Point de la pente à l'abscisse s, décalé de `normal` unités de scène perpendiculairement à la pente.
Vertex FrictionModule::onSlope(double s, double normal, float z, const float* color, float dim) const {
    const double th = problem_.angle, mid = 0.5 * (lo_ + hi_);
    const float x = scale_ * static_cast<float>((s - mid) * std::cos(th)) - static_cast<float>(normal * std::sin(th));
    const float y = scale_ * static_cast<float>((s - lo_) * std::sin(th)) + static_cast<float>(normal * std::cos(th));
    return {x, y, z, color[0] * dim, color[1] * dim, color[2] * dim};
}

void FrictionModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const double th = problem_.angle;
    int shown = 0;
    for (int i = 0; i < SolverSet::kCount; ++i) if (solvers_.isShown(ctx.level, i)) ++shown;
    const float zMax = 1.2f, zMin = -kLaneSpacing * static_cast<float>(shown) - 0.8f;

    // Sol et rampe : deux triangles (avant et arrière) reliés, et un repère tous les mètres.
    const Vertex low = onSlope(lo_, 0.0, 0.0f, kGrey), high = onSlope(hi_, 0.0, 0.0f, kGrey);
    std::vector<Vertex> frame;
    for (float z : {zMin, zMax}) {
        Vertex a = low, b = high, c = {high.x, 0.0f, z, kGrey[0], kGrey[1], kGrey[2]}, d = {low.x, 0.0f, z, kGrey[0], kGrey[1], kGrey[2]};
        a.z = b.z = z;
        frame.insert(frame.end(), {a, b, b, c, c, d, d, a});
    }
    for (const Vertex& v : {low, high, Vertex{high.x, 0.0f, 0.0f, kGrey[0], kGrey[1], kGrey[2]}, Vertex{low.x, 0.0f, 0.0f, kGrey[0], kGrey[1], kGrey[2]}}) {
        Vertex front = v, back = v;
        front.z = zMax;
        back.z = zMin;
        frame.insert(frame.end(), {front, back});
    }
    renderer.draw(Primitive::Lines, frame);

    std::vector<Vertex> ticks;
    const float white[3] = {1.0f, 1.0f, 1.0f};
    for (int k = static_cast<int>(std::ceil(lo_)); k <= static_cast<int>(std::floor(hi_)); ++k) {
        const float len = k == 0 ? 0.55f : 0.25f;
        const float* c = k == 0 ? white : kGrey;
        ticks.push_back(onSlope(k, 0.0, zMax, c));
        ticks.push_back(onSlope(k, len, zMax, c));
    }
    renderer.draw(Primitive::Lines, ticks);

    // Bloc : parallélépipède fil de fer posé sur la pente, aux deux faces z0 et z1.
    auto drawBlock = [&](double s, float z0, float z1, const float* color, float width) {
        const double sc = std::clamp(s, lo_, hi_);
        const Vertex base = onSlope(sc, 0.0, 0.0f, color);
        const float ux = static_cast<float>(std::cos(th)), uy = static_cast<float>(std::sin(th));
        const float nx = -uy, ny = ux;
        auto corner = [&](float along, float up, float z) {
            return Vertex{base.x + along * ux + up * nx, base.y + along * uy + up * ny, z, color[0], color[1], color[2]};
        };
        const float hl = kBlockHalfLength * width, hh = kBlockHeight * width;
        std::vector<Vertex> edges;
        for (float z : {z0, z1}) {
            const Vertex p0 = corner(-hl, 0, z), p1 = corner(hl, 0, z), p2 = corner(hl, hh, z), p3 = corner(-hl, hh, z);
            edges.insert(edges.end(), {p0, p1, p1, p2, p2, p3, p3, p0});
        }
        for (const std::pair<float, float>& pt : {std::pair<float, float>{-hl, 0}, {hl, 0}, {hl, hh}, {-hl, hh}}) {
            edges.push_back(corner(pt.first, pt.second, z0));
            edges.push_back(corner(pt.first, pt.second, z1));
        }
        renderer.draw(Primitive::Lines, edges);
        renderer.draw(Primitive::Points, {corner(0.0f, 0.5f * hh, 0.5f * (z0 + z1))}, 6.0f * ctx.uiScale);
    };

    // Bloc exact au premier plan, puis un couloir par méthode affichée.
    const InclineState ex = problem_.exact(clock_.time);
    drawBlock(ex.s, -0.6f, 0.6f, kBlue, 1.0f);
    int lane = 0;
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        ++lane;
        const Run& r = runs_[i];
        if (r.stalled) continue;
        const float zc = -kLaneSpacing * static_cast<float>(lane);
        drawBlock(r.y[0], zc - 0.5f, zc + 0.5f, solvers_.color(i), 0.9f);
    }

    // Forces sur le bloc exact (niveau 4+) : longueur proportionnelle à la force, le poids m g valant 2 unités.
    if (atLeast(ctx.level, Level::Lycee)) {
        const Vertex c = onSlope(std::clamp(ex.s, lo_, hi_), 0.5f * kBlockHeight, 0.0f, kBlue);
        const float len = 2.0f, ux = static_cast<float>(std::cos(th)), uy = static_cast<float>(std::sin(th));
        const float cosTh = static_cast<float>(std::cos(th));
        arrow(renderer, c.x, c.y, 0.0f, 0.0f, -len, kWeightColor, ctx.uiScale);                               // poids
        arrow(renderer, c.x, c.y, 0.0f, -uy * cosTh * len, ux * cosTh * len, kNormalColor, ctx.uiScale);      // réaction normale (m g cos θ)
        float friction = 0.0f;  // composante le long de la pente, positive vers le haut
        if (ex.stuck) friction = static_cast<float>(std::sin(th));
        else friction = static_cast<float>(problem_.muKinetic * std::cos(th)) * (ex.v > 0.0 ? -1.0f : 1.0f);
        if (friction != 0.0f) arrow(renderer, c.x, c.y, 0.0f, ux * friction * len, uy * friction * len, kFrictionColor, ctx.uiScale);
    }
}

}  // namespace pl
