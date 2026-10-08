#include "DoublePendulumModule.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

#include <implot.h>

namespace pl {
namespace {

constexpr double kPi = constants::pi;
constexpr float kPivotHeight = 5.2f;    // hauteur du pivot [m]
constexpr float kScale = 1.3f;          // 1 m de tige = 1,3 unité de scène
constexpr double kHorizonThreshold = 0.1;  // une méthode est "perdue" quand son écart à la référence dépasse 0,1
constexpr std::size_t kTrailLength = 700;
const float kTwinColor[3] = {0.55f, 0.90f, 1.0f};

// Tiges, masses et trace de la seconde masse à l'état y. `layer` décale légèrement en z pour que les pendules
// superposés restent distinguables.
void drawDoublePendulum(Renderer& renderer, const DoublePendulumProblem& p, const State& y, const float* color,
                        float layer, const std::vector<Vertex>& trail, float pointSize) {
    double x1, y1, x2, y2;
    p.positions(y, x1, y1, x2, y2);
    const float ax = static_cast<float>(x1) * kScale, ay = kPivotHeight + static_cast<float>(y1) * kScale;
    const float bx = static_cast<float>(x2) * kScale, by = kPivotHeight + static_cast<float>(y2) * kScale;

    if (!trail.empty()) renderer.draw(Primitive::LineStrip, trail);
    const std::vector<Vertex> rods = {{0.0f, kPivotHeight, layer, color[0], color[1], color[2]},
                                      {ax, ay, layer, color[0], color[1], color[2]},
                                      {ax, ay, layer, color[0], color[1], color[2]},
                                      {bx, by, layer, color[0], color[1], color[2]}};
    renderer.draw(Primitive::Lines, rods);
    renderer.draw(Primitive::Points, {{ax, ay, layer, color[0], color[1], color[2]},
                                      {bx, by, layer, color[0], color[1], color[2]}},
                  pointSize);
}

void appendTrail(std::vector<Vertex>& trail, const DoublePendulumProblem& p, const State& y, const float* color, float layer) {
    double x1, y1, x2, y2;
    p.positions(y, x1, y1, x2, y2);
    const float dim = 0.7f;
    trail.push_back({static_cast<float>(x2) * kScale, kPivotHeight + static_cast<float>(y2) * kScale, layer,
                     color[0] * dim, color[1] * dim, color[2] * dim});
    if (trail.size() > kTrailLength) trail.erase(trail.begin(), trail.begin() + 100);
}

}  // namespace

DoublePendulumModule::DoublePendulumModule() {
    refSolver_.relTol = twinSolver_.relTol = 1e-13;
    refSolver_.absTol = twinSolver_.absTol = 1e-15;
    reset();
}

const char* DoublePendulumModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Deux pendules accrochés bout à bout. Au début le mouvement paraît ordinaire, puis il devient "
                   "imprévisible : c'est le chaos.\n\n"
                   "Le pendule bleu est la référence, le vert est le calcul de l'ordinateur. Regardez : ils se confondent un "
                   "moment, puis se séparent. Même un ordinateur très précis finit par se tromper, parce que le chaos "
                   "agrandit la moindre erreur.\n\n"
                   "Essayez « Départ au hasard ».";
        case Level::Interesse:
            return "Chaos veut dire : un tout petit changement au départ entraîne, plus tard, un résultat complètement "
                   "différent. C'est l'« effet papillon » : c'est pourquoi la météo ne se prévoit pas à plus de quelques jours.\n\n"
                   "Ici l'ordinateur calcule avec de minuscules erreurs d'arrondi. Elles grossissent à chaque seconde et "
                   "le pendule vert finit par quitter le bleu. Regardez après combien de secondes : c'est l'horizon de "
                   "prédictibilité.";
        case Level::College:
            return "Position de chaque masse : x1 = l1 sin θ1 ; y1 = −l1 cos θ1 ; x2 = x1 + l2 sin θ2 ; y2 = y1 − l2 cos θ2.\n"
                   "Énergie de hauteur : Ep = m × g × h, énergie de mouvement : Ec = ½ m v². Sans frottement, le total est "
                   "constant.\n\n"
                   "La méthode simple (Euler, orange) crée de l'énergie : le pendule s'emballe et fait des tours complets. "
                   "La méthode précise (RK4, vert) conserve l'énergie, mais finit tout de même par s'écarter de la "
                   "référence, à cause du chaos.";
        case Level::Lycee:
            return "Écrire la 2e loi de Newton avec les tensions des deux tiges est pénible : on passe par l'énergie. "
                   "T = ½ (m1+m2) l1² θ1'² + ½ m2 l2² θ2'² + m2 l1 l2 θ1' θ2' cos(θ1 − θ2) ; "
                   "V = −(m1+m2) g l1 cos θ1 − m2 g l2 cos θ2.\n\n"
                   "Le jumeau (cyan) part avec un écart de 1e-9 rad. Son écart d à la référence croît comme "
                   "d(t) ≈ d0 e^(λ t) : exponentiellement (échelle logarithmique : une droite). Le temps de Lyapunov "
                   "1/λ dit en combien de temps l'erreur est multipliée par e ≈ 2,7.\n"
                   "Un pendule simple, lui, n'amplifie l'écart que lentement.";
        case Level::Etudiant:
            return "Lagrangien L = T − V ; équations d'Euler-Lagrange, résolues en θ1'' et θ2'' :\n"
                   "θ1'' = [−g(2m1+m2) sin θ1 − m2 g sin(θ1−2θ2) − 2 sin(Δ) m2 (θ2'² l2 + θ1'² l1 cos Δ)] / (l1 D)\n"
                   "θ2'' = [2 sin(Δ)(θ1'² l1 (m1+m2) + g(m1+m2) cos θ1 + θ2'² l2 m2 cos Δ)] / (l2 D)\n"
                   "avec Δ = θ1 − θ2 et D = 2m1 + m2 − m2 cos 2Δ. E = T + V est conservée.\n\n"
                   "Le hamiltonien n'est pas séparable (la cinétique dépend des angles) : Euler symplectique et Verlet "
                   "restent d'ordre 1 et 2 mais ne sont plus symplectiques. L'ordre de convergence se mesure à t = 2 s, "
                   "avant que le chaos n'amplifie tout.";
        case Level::Chercheur:
            return "Exposant de Lyapunov maximal : λ = lim 1/t · ln(d(t) / d0). λ > 0 caractérise le chaos ; ici il est "
                   "estimé par moindres carrés sur la phase exponentielle (20 d0 < d < 0,1), avant saturation.\n"
                   "Horizon de prédictibilité : t_H ≈ (1/λ) ln(Δ / d0). Améliorer la précision d'un facteur 10 ne gagne que "
                   "ln(10)/λ secondes : en chaos, la précision se paie très cher. En simple précision (GPU, d0 ≈ 1e-7) l'horizon "
                   "est plus court qu'en double précision (1e-16) : voir M7.\n"
                   "Pas de solution analytique : la référence est un RK45 de tolérance 1e-13, fiable tant que l'erreur "
                   "amplifiée reste petite (t ≲ (1/λ) ln(1/1e-13)).\n"
                   "Régimes : à basse énergie les orbites sont régulières (tores KAM), à haute énergie la mer chaotique "
                   "envahit l'espace des phases ; sections de Poincaré.";
    }
    return "";
}

void DoublePendulumModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 3.4f;
    camera.target[2] = 0.0f;
    camera.distance = 13.0f;
    camera.yaw = 0.3f;
    camera.pitch = 0.15f;
}

// ------------------------------ simulation -----------------------------

void DoublePendulumModule::reset() {
    problem_.m1 = m1_;
    problem_.m2 = m2_;
    problem_.l1 = l1_;
    problem_.l2 = l2_;
    problem_.gravity = gravity_;
    problem_.theta1 = theta1Deg_ * kPi / 180.0;
    problem_.theta2 = theta2Deg_ * kPi / 180.0;
    problem_.omega1 = omega1_;
    problem_.omega2 = omega2_;
    rhs_ = problem_.rhs();

    for (Run& r : runs_) {
        r.y = problem_.initialState();
        r.theta2.clear();
        r.energy.clear();
        r.error.clear();
        r.trail.clear();
        r.diverged = false;
        r.horizon = -1.0;
    }
    refY_ = problem_.initialState();
    twinY_ = refY_;
    twinY_[0] += perturbation_;
    theta2Ref_.clear();
    configX_.clear();
    refTrail_.clear();
    twinTrail_.clear();
    twinT_.clear();
    twinD_.clear();
    e0_ = problem_.energy(refY_);
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void DoublePendulumModule::sample() {
    const double t = clock_.time;
    theta2Ref_.add(t, refY_[1]);
    configX_.add(refY_[0], refY_[1]);
    appendTrail(refTrail_, problem_, refY_, kBlue, 0.0f);
    appendTrail(twinTrail_, problem_, twinY_, kTwinColor, 0.0f);
    twinT_.push_back(t);
    twinD_.push_back(std::max(problem_.distance(refY_, twinY_), 1e-300));

    for (int i = 0; i < SolverSet::kCount; ++i) {
        Run& r = runs_[i];
        if (r.diverged) continue;
        r.theta2.add(t, r.y[1]);
        r.energy.add(t, problem_.energy(r.y) - e0_);
        const double d = problem_.distance(r.y, refY_);
        r.error.add(t, std::max(d, 1e-12));
        if (r.horizon < 0.0 && d > kHorizonThreshold) r.horizon = t;
        appendTrail(r.trail, problem_, r.y, solvers_.color(i), 0.04f * static_cast<float>(i + 1));
    }
}

void DoublePendulumModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;

    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    clock_.advance(
        frameSeconds, timeScale_, dt_, durationSec_,
        [&](double h) {
            for (int i = 0; i < SolverSet::kCount; ++i) {
                Run& r = runs_[i];
                if (r.diverged) continue;
                advance(solvers_.solver(i), rhs_, clock_.time, r.y, h);
                if (!std::isfinite(r.y[0]) || !std::isfinite(r.y[1]) || !std::isfinite(r.y[2])) r.diverged = true;
            }
            advance(refSolver_, rhs_, clock_.time, refY_, h);
            advance(twinSolver_, rhs_, clock_.time, twinY_, h);
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur finale à t = 2 s de chaque solveur à pas fixe : avant que le chaos n'amplifie tout, l'ordre se mesure.
void DoublePendulumModule::computeConvergence() {
    static const int kStepCounts[] = {400, 800, 1600, 3200, 6400, 12800, 25600};

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        ConvergenceCurve& c = convergence_[i];
        c.dt.clear();
        c.err.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        for (int n : kStepCounts) {
            const double e = doublePendulumError(problem_, *solver, n, 2.0);
            if (e > 1e-10) {  // sous ce seuil on mesure la précision de la référence
                c.dt.push_back(2.0 / n);
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

void DoublePendulumModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    bool changed = false;

    pushSliderWidth();
    ImGui::SeparatorText("Départ");
    changed |= sliderD(college ? "Angle du haut θ1 (deg)" : "Écart du haut", &theta1Deg_, -179.0, 179.0, college ? "%.1f" : "");
    changed |= sliderD(college ? "Angle du bas θ2 (deg)" : "Écart du bas", &theta2Deg_, -179.0, 179.0, college ? "%.1f" : "");
    if (atLeast(level, Level::Etudiant)) {
        changed |= sliderD("Vitesse θ1' (rad/s)", &omega1_, -8.0, 8.0, "%.2f");
        changed |= sliderD("Vitesse θ2' (rad/s)", &omega2_, -8.0, 8.0, "%.2f");
    }
    popSliderWidth();

    if (atLeast(level, Level::Interesse) && ImGui::Button("Départ au hasard")) {
        static std::mt19937 rng{std::random_device{}()};
        std::uniform_real_distribution<double> angle(40.0, 170.0), sign(0.0, 1.0);
        theta1Deg_ = angle(rng) * (sign(rng) < 0.5 ? -1.0 : 1.0);
        theta2Deg_ = angle(rng) * (sign(rng) < 0.5 ? -1.0 : 1.0);
        changed = true;
    }

    if (college) {
        ImGui::SeparatorText("Pendule");
        pushSliderWidth();
        changed |= sliderD("Masse du haut m1 (kg)", &m1_, 0.2, 5.0, "%.2f", ImGuiSliderFlags_Logarithmic);
        changed |= sliderD("Masse du bas m2 (kg)", &m2_, 0.2, 5.0, "%.2f", ImGuiSliderFlags_Logarithmic);
        changed |= sliderD("Tige du haut l1 (m)", &l1_, 0.3, 2.0, "%.2f");
        changed |= sliderD("Tige du bas l2 (m)", &l2_, 0.3, 2.0, "%.2f");
        changed |= sliderD("Pesanteur g (m/s²)", &gravity_, 0.5, 30.0, "%.2f");
        popSliderWidth();

        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Etudiant)) changed |= sliderD("Durée (s)", &durationSec_, 5.0, 60.0, "%.0f");
        if (atLeast(level, Level::Chercheur)) {
            changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
            changed |= sliderD("Écart du jumeau (rad)", &perturbation_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
        }
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
    if (atLeast(level, Level::Interesse))
        sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        ImGui::Text("t = %.2f s / %.0f s", clock_.time, durationSec_);
        ImGui::Text("E = %.3f J", problem_.energy(refY_));
    }
    if (atLeast(level, Level::Lycee)) {
        int used = 0;
        const double lambda = lyapunovExponent(twinT_, twinD_, 20.0 * perturbation_, 0.1, &used);
        ImGui::SeparatorText("Chaos");
        if (std::isfinite(lambda)) {
            ImGui::Text("λ = %.3f 1/s   (%d points)", lambda, used);
            ImGui::Text("Temps de Lyapunov 1/λ = %.2f s", 1.0 / lambda);
            if (atLeast(level, Level::Etudiant) && lambda > 0.0)
                ImGui::Text("Horizon de départ ≈ %.1f s", std::log(0.1 / perturbation_) / lambda);
        } else {
            ImGui::TextDisabled("λ : trop peu de points exploitables pour l'instant");
        }
    }
}

void DoublePendulumModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant le mouvement, l'énergie se partage entre les deux masses et passe de la hauteur au "
                           "mouvement. Sans frottement, leur total reste le même.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool researcher = atLeast(level, Level::Chercheur);
    const double t = clock_.time;

    std::vector<std::string> headers;
    headers.push_back(researcher ? "dE/E0" : "dE (J)");
    if (lycee) headers.push_back("dist");
    headers.push_back("t_H (s)");

    std::vector<TableRow> rows;
    {
        TableRow ref{"Référence", kBlue, {}};
        ref.cells.push_back(strf(researcher ? "%+.1e" : "%+.1e", (problem_.energy(refY_) - e0_) / (researcher ? e0_ : 1.0)));
        if (lycee) ref.cells.push_back("-");
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
            const double dE = problem_.energy(r.y) - e0_;
            row.cells.push_back(strf("%+.1e", researcher ? dE / e0_ : dE));
            if (lycee) row.cells.push_back(strf("%.1e", problem_.distance(r.y, refY_)));
            row.cells.push_back(r.horizon >= 0.0 ? strf("%.2f", r.horizon) : "-");
        }
        rows.push_back(row);
    }
    drawResultTable("erreurs", headers, rows);
    ImGui::TextDisabled("t_H : instant où l'écart à la référence dépasse %.1f.", kHorizonThreshold);
    ImGui::TextDisabled("« - » : la méthode suit encore la référence.");

    if (lycee) {
        ImGui::SeparatorText("Jumeau");
        ImGui::Text("écart initial = %.0e rad", perturbation_);
        ImGui::Text("écart à t = %.1f s : %.2e", t, twinD_.empty() ? 0.0 : twinD_.back());
        if (!twinD_.empty() && twinD_.front() > 0.0 && twinD_.back() > 0.0)
            ImGui::Text("amplification : x %.1e", twinD_.back() / std::max(twinD_.front(), 1e-300));
    }
    if (researcher) {
        ImGui::SeparatorText("RK45");
        ImGui::Text("%d pas acceptés, %d rejetés", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps());
        ImGui::Text("%d évaluations de f", solvers_.rk45().evaluations());
    }
}

void DoublePendulumModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showEnergy = atLeast(level, Level::College);
    const bool showError = atLeast(level, Level::Lycee);
    const int cols = 1 + (showEnergy ? 1 : 0) + (showError ? 1 : 0);

    auto plotRuns = [&](Series Run::*member) {
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i) || (runs_[i].*member).size() == 0) continue;
            const Series& s = runs_[i].*member;
            ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                             lineSpec(solvers_.color(i), s.offset));
        }
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Angle du bas θ2(t)")) {
            ImPlot::SetupAxes("t (s)", "θ2 (rad)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (theta2Ref_.size() > 0)
                ImPlot::PlotLine("Référence", theta2Ref_.x.data(), theta2Ref_.y.data(), theta2Ref_.size(),
                                 lineSpec(kBlue, theta2Ref_.offset));
            plotRuns(&Run::theta2);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie E(t) − E(0)")) {
            ImPlot::SetupAxes("t (s)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Écart à la référence")) {
            ImPlot::SetupAxes("t (s)", "distance (espace des phases)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            const double threshold = kHorizonThreshold;
            const float grey[3] = {0.55f, 0.55f, 0.6f};
            ImPlot::PlotInfLines("seuil", &threshold, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(grey), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void DoublePendulumModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool etudiant = atLeast(level, Level::Etudiant);
    if (etudiant && convergenceDirty_) computeConvergence();
    const int cols = 1 + (etudiant ? 2 : 0);

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        // Sensibilité aux conditions initiales : d(t) croît exponentiellement, la pente de ln d est lambda.
        if (ImPlot::BeginPlot("Sensibilité : écart référence - jumeau")) {
            ImPlot::SetupAxes("t (s)", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (!twinT_.empty())
                ImPlot::PlotLine("Écart", twinT_.data(), twinD_.data(), static_cast<int>(twinT_.size()), lineSpec(kTwinColor, 0));

            int used = 0;
            double intercept = 0.0;
            const double lambda = lyapunovExponent(twinT_, twinD_, 20.0 * perturbation_, 0.1, &used, &intercept);
            if (std::isfinite(lambda)) {
                // droite d0 e^(lambda t) sur la durée de la phase exponentielle
                std::vector<double> xs, ys;
                const double tEnd = std::max(clock_.time, 1.0);
                for (int i = 0; i <= 50; ++i) {
                    const double t = tEnd * i / 50.0;
                    xs.push_back(t);
                    ys.push_back(std::min(std::exp(intercept + lambda * t), 10.0));
                }
                const float orange[3] = {1.0f, 0.55f, 0.15f};
                const std::string name = strf("e^(λt), λ = %.2f /s", lambda);
                ImPlot::PlotLine(name.c_str(), xs.data(), ys.data(), static_cast<int>(xs.size()), lineSpec(orange, 0));
            }
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot("Espace des configurations (θ1, θ2)")) {
            ImPlot::SetupAxes("θ1 (rad)", "θ2 (rad)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (configX_.size() > 0)
                ImPlot::PlotLine("Référence", configX_.x.data(), configX_.y.data(), configX_.size(), lineSpec(kBlue, configX_.offset));
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot("Convergence : erreur à t = 2 s")) {
            ImPlot::SetupAxes("dt (s)", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kFixedStep; ++i) {
                const ConvergenceCurve& c = convergence_[i];
                if (!solvers_.show(i) || c.dt.size() < 2) continue;
                const std::string name = strf("%s (pente %.2f)", solvers_.solver(i).name(), c.slope);
                const ImPlotSpec markers(ImPlotProp_LineColor, toImVec4(solvers_.color(i)), ImPlotProp_LineWeight, 2.0f,
                                         ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f);
                ImPlot::PlotLine(name.c_str(), c.dt.data(), c.err.data(), static_cast<int>(c.dt.size()), markers);
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

void DoublePendulumModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    // Support.
    const std::vector<Vertex> support = {{-0.8f, kPivotHeight, 0.0f, 0.5f, 0.5f, 0.55f}, {0.8f, kPivotHeight, 0.0f, 0.5f, 0.5f, 0.55f}};
    renderer.draw(Primitive::Lines, support);

    const std::vector<Vertex> none;
    // Jumeau (niveau 4+), puis méthodes, puis la référence par-dessus.
    if (atLeast(ctx.level, Level::Lycee))
        drawDoublePendulum(renderer, problem_, twinY_, kTwinColor, -0.04f, twinTrail_, 9.0f * ctx.uiScale);
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        const Run& r = runs_[i];
        if (r.diverged) continue;
        drawDoublePendulum(renderer, problem_, r.y, solvers_.color(i), 0.04f * static_cast<float>(i + 1), r.trail,
                           11.0f * ctx.uiScale);
    }
    drawDoublePendulum(renderer, problem_, refY_, kBlue, 0.0f, refTrail_, 15.0f * ctx.uiScale);
}

}  // namespace pl
