#include "ProjectileModule.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

#include <implot.h>

namespace pl {

ProjectileModule::ProjectileModule() { reset(); }

const char* ProjectileModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Quand on lance une balle, elle monte, ralentit, puis retombe : la Terre l'attire vers le bas, "
                   "et l'air la freine un peu.\n\n"
                   "La ligne bleue montre ce que fait vraiment la nature. Le point vert est l'ordinateur, qui calcule "
                   "le mouvement pas à pas : il suit presque parfaitement.";
        case Level::Interesse:
            return "Sans air, la balle décrit une parabole. Avec de l'air, elle est freinée : elle monte moins haut et "
                   "retombe plus près. Sur la Lune (pas d'air, pesanteur six fois plus faible), la même balle irait "
                   "beaucoup plus loin.\n\n"
                   "Le mouvement vers l'avant et le mouvement vers le haut sont indépendants : la gravité n'agit "
                   "que vers le bas.";
        case Level::College:
            return "Poids : P = m x g  (g = 9,81 N/kg sur Terre). Vitesse : v = d / t.\n"
                   "Énergie de mouvement : Ec = 1/2 x m x v^2 ; énergie de hauteur : Ep = m x g x h.\n\n"
                   "Un ordinateur ne sait pas suivre un mouvement continu : il avance par petits pas de durée dt.\n"
                   "- Méthode simple (Euler, orange) : elle suppose la vitesse constante pendant dt. Elle se trompe un "
                   "peu à chaque pas et les erreurs s'additionnent.\n"
                   "- Méthode précise (RK4, vert) : elle regarde la vitesse à plusieurs instants du pas.\n\n"
                   "Essayez d'augmenter dt : l'orange s'écarte de la courbe bleue, le vert reste dessus.";
        case Level::Lycee:
            return "2e loi de Newton : m a = m g - b v  (poids + frottement opposé à la vitesse).\n\n"
                   "Sans frottement (b = 0) :  x = x0 + v0x t ;  y = y0 + v0y t - 1/2 g t^2  (parabole).\n"
                   "Avec frottement, k = b / m : la vitesse tend vers la vitesse limite v_lim = g / k "
                   "(verticale : m g / b).\n\n"
                   "Erreur d'Euler : proportionnelle à dt. Diviser dt par 2 divise l'erreur par 2 ; pour RK4 on la "
                   "divise par 16.\n"
                   "Énergie mécanique : Em = Ec + Ep, conservée sans frottement, dissipée avec.";
        case Level::Etudiant:
            return "EDO : dx/dt = v ; dv/dt = g - k v  (k = b/m). Lagrangien avec dissipation de Rayleigh :\n"
                   "d/dt (dL/dq') - dL/dq = - b q'.\n\n"
                   "Solution exacte : v(t) = v0 e^(-kt) + g phi(t) ; r(t) = r0 + v0 phi(t) + g psi(t),\n"
                   "phi = (1 - e^(-kt)) / k ; psi = (t - phi) / k.\n\n"
                   "Ordres des méthodes : Euler 1, Euler symplectique 1, Verlet des vitesses 2, RK4 4. "
                   "Le panneau Analyse trace l'erreur finale en fonction de dt : en log-log, la pente est l'ordre.";
        case Level::Chercheur:
            return "RK45 de Dormand-Prince : paire imbriquée d'ordres 5 et 4. L'erreur locale estimée e = ||y5 - y4|| "
                   "(norme RMS pondérée par absTol + relTol |y|) décide du pas : accepté si e <= 1, puis "
                   "h_nouveau = 0,9 h e^(-1/5) borné dans [0,2 h ; 5 h]. 7 évaluations de f par essai.\n\n"
                   "Euler : erreur locale O(dt^2), globale O(dt) (consistance + stabilité => convergence). "
                   "Sans frottement, dérive d'énergie d'Euler : dE = 1/2 m g^2 t dt.\n\n"
                   "Verlet et RK4 intègrent exactement une accélération constante (solution polynomiale de degré <= 2, "
                   "<= 4). Le caractère symplectique ne se voit qu'avec une force non constante : module M2 (ressort).";
    }
    return "";
}

void ProjectileModule::frameCamera(Camera& camera) const {
    camera.target[0] = 10.0f;
    camera.target[1] = 3.0f;
    camera.target[2] = 0.0f;
    camera.distance = 38.0f;
    camera.yaw = 0.55f;
    camera.pitch = 0.30f;
}

// ------------------------------ simulation -----------------------------

void ProjectileModule::reset() {
    const double a = angleDeg_ * constants::pi / 180.0;
    problem_.position0 = {0.0, height_, 0.0};
    problem_.velocity0 = {speed_ * std::cos(a), speed_ * std::sin(a), 0.0};
    problem_.mass = mass_;
    problem_.gravity = {0.0, -gravity_, 0.0};
    problem_.linearDrag = drag_;
    landingTime_ = problem_.landingTime();

    for (Run& r : runs_) {
        r.world = problem_.makeWorld();
        r.y.clear();
        r.posError.clear();
        r.energyError.clear();
        r.trail.clear();
    }
    refY_.clear();
    solvers_.rk45().resetStats();

    refPath_.clear();
    for (int i = 0; i <= 120; ++i) {
        const Vec3 r = problem_.position(landingTime_ * i / 120.0);
        refPath_.push_back({static_cast<float>(r.x), static_cast<float>(r.y), static_cast<float>(r.z),
                            kBlue[0], kBlue[1], kBlue[2]});
    }

    accumulator_ = 0.0;
    landed_ = false;
    running_ = true;
    sample();
    computeConvergence();
}

void ProjectileModule::sample() {
    const double t = runs_[0].world.time;
    lastSampleTime_ = t;
    const Vec3 ref = problem_.position(t);
    const double exactEnergy = problem_.energy(t);

    refY_.add(t, ref.y);
    for (int i = 0; i < SolverSet::kCount; ++i) {
        Run& r = runs_[i];
        const Particle& p = r.world.particles[0];
        const float* c = solvers_.color(i);
        r.y.add(t, p.position.y);
        r.posError.add(t, std::max((p.position - ref).norm(), 1e-12));  // plancher : l'axe est logarithmique
        r.energyError.add(t, r.world.invariants().total() - exactEnergy);
        r.trail.push_back({static_cast<float>(p.position.x), static_cast<float>(p.position.y),
                           static_cast<float>(p.position.z), c[0], c[1], c[2]});
    }
}

void ProjectileModule::update(double frameSeconds) {
    if (!running_ || landed_) return;

    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;  // borne anti "spirale de la mort"
    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    const double sampleInterval = std::max(dt_, 1.0 / 240.0);

    int guard = 0;
    while (!landed_ && guard++ < 5000) {
        const double remaining = landingTime_ - runs_[0].world.time;
        const double h = std::min(dt_, remaining);  // dernier pas raccourci : tous arrivent pile à l'impact exact
        if (accumulator_ < h) break;

        for (int i = 0; i < SolverSet::kCount; ++i) runs_[i].world.step(solvers_.solver(i), h);
        accumulator_ -= h;
        landed_ = remaining <= dt_ * (1.0 + 1e-9);
        if (landed_ || runs_[0].world.time - lastSampleTime_ >= sampleInterval - 1e-12) sample();
    }
    if (landed_) running_ = false;
}

// Erreur finale de chaque solveur à pas fixe pour 7 pas différents, puis pente log-log (ordre mesuré).
void ProjectileModule::computeConvergence() {
    static const int kStepCounts[] = {5, 10, 20, 40, 80, 160, 320};
    convergenceEnd_ = std::max(0.2, landingTime_);

    // Sans frottement la solution est un polynôme de degré 2 que Verlet et RK4 intègrent exactement :
    // il n'y aurait rien à mesurer. On impose alors un frottement de référence (k = 0,8 /s).
    ProjectileProblem study = problem_;
    convergenceForcedDrag_ = study.linearDrag < 0.3 * study.mass;
    if (convergenceForcedDrag_) study.linearDrag = 0.8 * study.mass;

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        ConvergenceCurve& c = convergence_[i];
        c.dt.clear();
        c.err.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        for (int n : kStepCounts) {
            const double e = integrationError(study, *solver, n, convergenceEnd_);
            if (e > 1e-12) {  // sous ce seuil on mesure l'arrondi, pas la méthode
                c.dt.push_back(convergenceEnd_ / n);
                c.err.push_back(e);
            }
        }

        // Moindres carrés sur (log10 dt, log10 err).
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
}

// --------------------------------- UI ----------------------------------

void ProjectileModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    bool changed = false;

    pushSliderWidth();
    ImGui::SeparatorText("Lancer");
    changed |= sliderD(college ? "Vitesse (m/s)" : "Force du lancer", &speed_, 1.0, 30.0, college ? "%.1f" : "");
    if (atLeast(level, Level::Interesse)) changed |= sliderD("Angle (deg)", &angleDeg_, 5.0, 85.0, "%.0f");

    if (level == Level::Interesse) {
        static const char* planets[] = {"Terre", "Lune", "Mars", "Jupiter"};
        static const double gs[] = {constants::g0, 1.62, 3.71, 24.79};
        static int planet = 0;
        if (ImGui::Combo("Planète", &planet, planets, 4)) {
            gravity_ = gs[planet];
            changed = true;
        }
    }
    if (college) {
        changed |= sliderD("Pesanteur g (m/s²)", &gravity_, 0.5, 30.0, "%.2f");
        changed |= sliderD("Masse (kg)", &mass_, 0.1, 20.0, "%.1f", ImGuiSliderFlags_Logarithmic);
    }
    if (atLeast(level, Level::Lycee)) {
        changed |= sliderD("Hauteur initiale (m)", &height_, 0.0, 20.0, "%.1f");
        const double a = angleDeg_ * constants::pi / 180.0;
        ImGui::TextDisabled("v0 = (%.2f ; %.2f) m/s", speed_ * std::cos(a), speed_ * std::sin(a));
    }
    changed |= sliderD(college ? "Frottement b (kg/s)" : "Résistance de l'air", &drag_, 0.0, 3.0, college ? "%.2f" : "");

    if (college) {
        ImGui::SeparatorText("Calcul");
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
    }
    if (atLeast(level, Level::Chercheur))
        changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    solvers_.drawToggles(level);
    if (changed) reset();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !landed_ ? "Pause" : "Lecture")) {
        if (landed_) reset();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) reset();
    pushSliderWidth();
    if (atLeast(level, Level::Interesse))
        sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        const Run& ref = runs_[SolverSet::kRK4];
        const Particle& p = ref.world.particles[0];
        ImGui::Text("t = %.3f s", ref.world.time);
        ImGui::Text("x = %.2f m   y = %.2f m", p.position.x, p.position.y);
        ImGui::Text("vitesse = %.2f m/s", p.velocity.norm());
    }
    if (landed_) {
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Impact au sol !");
        if (college) ImGui::Text("portée = %.2f m, durée = %.3f s", problem_.position(landingTime_).x, landingTime_);
    }
}

void ProjectileModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant le vol, l'énergie du mouvement et l'énergie de la hauteur s'échangent. "
                           "Sans air, leur total reste le même.");
        return;
    }

    const double t = runs_[0].world.time;
    const Vec3 exactPos = problem_.position(t);
    const double exactEnergy = problem_.energy(t);
    const double e0 = problem_.energy(0.0);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool researcher = atLeast(level, Level::Chercheur);

    // Chercheur : dE et dE/E0 remplacent Em (redondant) pour que le tableau tienne dans le panneau.
    std::vector<std::string> headers;
    if (!researcher) headers.push_back("Em (J)");
    if (lycee) headers.push_back("dE (J)");
    if (researcher) headers.push_back("dE/E0");
    headers.push_back("dr (m)");

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
        const double energy = r.world.invariants().total();
        TableRow row{solvers_.shortLabel(level, i), solvers_.color(i), {}};
        if (!researcher) row.cells.push_back(strf("%.3f", energy));
        if (lycee) row.cells.push_back(strf("%+.1e", energy - exactEnergy));
        if (researcher) row.cells.push_back(strf("%+.1e", (energy - exactEnergy) / e0));
        row.cells.push_back(strf("%.1e", (r.world.particles[0].position - exactPos).norm()));
        rows.push_back(row);
    }
    drawResultTable("erreurs", headers, rows);

    if (lycee) {
        const Vec3 v = problem_.velocity(t);
        const Vec3 L = problem_.mass * cross(exactPos, v);
        ImGui::SeparatorText("Solution exacte");
        ImGui::Text("px = %7.3f kg.m/s", problem_.mass * v.x);
        ImGui::Text("py = %7.3f kg.m/s", problem_.mass * v.y);
        ImGui::Text("Lz = %7.3f kg.m²/s", L.z);
        if (drag_ > 0.0) ImGui::TextDisabled("Le frottement dissipe E et p.");
        else ImGui::TextDisabled("Sans frottement : px et Em constantes.");
    }
    if (researcher) {
        ImGui::SeparatorText("RK45");
        ImGui::Text("%d pas acceptés, %d rejetés", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps());
        ImGui::Text("%d évaluations de f", solvers_.rk45().evaluations());
        if (drag_ == 0.0) {
            const double predicted = 0.5 * mass_ * gravity_ * gravity_ * t * dt_;
            ImGui::TextDisabled("Euler : dE théorique = %+.3e J", predicted);
            ImGui::TextDisabled("(1/2 m g² t dt)");
        }
    }
}

void ProjectileModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showError = atLeast(level, Level::College);
    const bool showEnergy = atLeast(level, Level::Lycee);
    const int cols = 1 + (showError ? 1 : 0) + (showEnergy ? 1 : 0);

    auto spec = [](const float c[3], int offset) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Offset, offset);
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Hauteur y(t)")) {
            ImPlot::SetupAxes("t (s)", "y (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (refY_.size() > 0) ImPlot::PlotLine("Exacte", refY_.x.data(), refY_.y.data(), refY_.size(), spec(kBlue, refY_.offset));
            for (int i = 0; i < SolverSet::kCount; ++i) {
                if (!solvers_.isShown(level, i) || runs_[i].y.size() == 0) continue;
                const Series& s = runs_[i].y;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                                 spec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Erreur de position |r - r_exact|")) {
            ImPlot::SetupAxes("t (s)", "m", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kCount; ++i) {
                if (!solvers_.isShown(level, i) || runs_[i].posError.size() == 0) continue;
                const Series& s = runs_[i].posError;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                                 spec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Erreur d'énergie Em - Em_exacte")) {
            ImPlot::SetupAxes("t (s)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            for (int i = 0; i < SolverSet::kCount; ++i) {
                if (!solvers_.isShown(level, i) || runs_[i].energyError.size() == 0) continue;
                const Series& s = runs_[i].energyError;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                                 spec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void ProjectileModule::drawAnalysis(const UiContext&) {
    ImGui::TextWrapped("Erreur de position à t = %.2f s en fonction du pas dt. En échelle log-log, la pente est "
                       "l'ordre de la méthode (RK45, adaptatif, n'y figure pas).",
                       convergenceEnd_);
    if (convergenceForcedDrag_)
        ImGui::TextDisabled("Frottement k = 0,8 /s imposé : sans frottement, Verlet et RK4 sont exacts (polynôme).");
    if (ImPlot::BeginPlot("##convergence", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("dt (s)", "erreur (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        for (int i = 0; i < SolverSet::kFixedStep; ++i) {
            const ConvergenceCurve& c = convergence_[i];
            if (!solvers_.show(i) || c.dt.size() < 2) continue;
            char name[96];
            std::snprintf(name, sizeof(name), "%s (pente %.2f)", solvers_.solver(i).name(), c.slope);
            const ImPlotSpec spec(ImPlotProp_LineColor, toImVec4(solvers_.color(i)), ImPlotProp_LineWeight, 2.0f,
                                  ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f);
            ImPlot::PlotLine(name, c.dt.data(), c.err.data(), static_cast<int>(c.dt.size()), spec);
        }
        ImPlot::EndPlot();
    }
}

// ------------------------------ rendu 3D --------------------------------

void ProjectileModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const Level level = ctx.level;

    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    renderer.draw(Primitive::LineStrip, refPath_);  // solution exacte
    for (int i = 0; i < SolverSet::kCount; ++i)
        if (solvers_.isShown(level, i)) renderer.draw(Primitive::LineStrip, runs_[i].trail);

    // Billes : l'exacte en bleu (grosse), les numériques par-dessus (plus petites).
    const Vec3 exact = problem_.position(runs_[0].world.time);
    const std::vector<Vertex> exactBall = {{static_cast<float>(exact.x), static_cast<float>(exact.y),
                                            static_cast<float>(exact.z), kBlue[0], kBlue[1], kBlue[2]}};
    renderer.draw(Primitive::Points, exactBall, 18.0f * ctx.uiScale);
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(level, i)) continue;
        const Vec3 p = runs_[i].world.particles[0].position;
        const float* c = solvers_.color(i);
        const std::vector<Vertex> ball = {{static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z),
                                           c[0], c[1], c[2]}};
        renderer.draw(Primitive::Points, ball, 11.0f * ctx.uiScale);
    }
}

}  // namespace pl
