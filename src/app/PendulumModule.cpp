#include "PendulumModule.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

#include <implot.h>

#include "physicslab/mechanics/Oscillator.hpp"

namespace pl {
namespace {

constexpr double kPi = constants::pi;
constexpr float kLaneSpacing = 2.6f;  // distance entre deux couloirs [m]
constexpr float kPivotHeight = 5.8f;  // hauteur du pivot [m]
constexpr float kLengthScale = 1.5f;  // 1 m de fil = 1,5 unité de scène

float laneZ(int lane) { return kLaneSpacing * static_cast<float>(lane); }

// Un couloir : pivot, tige et masse à l'angle theta.
void drawPendulumLane(Renderer& renderer, int lane, double theta, double length, const float* color, float uiScale) {
    const float z = laneZ(lane);
    const float r = static_cast<float>(length) * kLengthScale;
    const float bx = r * static_cast<float>(std::sin(theta));
    const float by = kPivotHeight - r * static_cast<float>(std::cos(theta));
    const float dim[3] = {color[0] * 0.55f, color[1] * 0.55f, color[2] * 0.55f};

    const std::vector<Vertex> lines = {
        {-0.7f, kPivotHeight, z, dim[0], dim[1], dim[2]}, {0.7f, kPivotHeight, z, dim[0], dim[1], dim[2]},  // support
        {0.0f, kPivotHeight, z, color[0], color[1], color[2]}, {bx, by, z, color[0], color[1], color[2]},    // tige
        {0.0f, kPivotHeight, z, dim[0], dim[1], dim[2]}, {0.0f, 0.0f, z, dim[0] * 0.5f, dim[1] * 0.5f, dim[2] * 0.5f},  // verticale
    };
    renderer.draw(Primitive::Lines, lines);
    renderer.draw(Primitive::Points, {{bx, by, z, color[0], color[1], color[2]}}, 20.0f * uiScale);
}

}  // namespace

PendulumModule::PendulumModule() {
    refSolver_.relTol = 1e-13;
    refSolver_.absTol = 1e-15;
    reset();
}

const char* PendulumModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Un pendule, c'est une masse au bout d'un fil. Écarté puis lâché, il va et vient à un rythme régulier : "
                   "c'est ce qui faisait marcher les horloges anciennes. Plus le fil est long, plus il est lent.\n\n"
                   "Le pendule bleu est ce que fait vraiment la nature, le vert est le calcul de l'ordinateur.";
        case Level::Interesse:
            return "La masse ne change rien au rythme (Galilée l'a observé dans la cathédrale de Pise) : seuls comptent la "
                   "longueur du fil et la pesanteur. Essayez de changer la masse. Sur la Lune, où la gravité est six fois plus "
                   "faible, le même pendule bat bien plus lentement.\n\n"
                   "Pour de petits écarts, le rythme ne dépend pas de l'écart ; pour un grand écart, le pendule est un peu "
                   "plus lent.";
        case Level::College:
            return "Période (petits écarts) : T = 2π √(L / g). Exemple : L = 1 m, g = 9,81 m/s² donne T ≈ 2,01 s. "
                   "Fréquence : f = 1 / T.\n"
                   "Énergie de hauteur : Ep = m × g × h, avec h = L (1 − cos θ). Énergie de mouvement : Ec = ½ m v².\n\n"
                   "La méthode simple (Euler, orange) ajoute de l'énergie à chaque pas : le pendule monte de plus en plus "
                   "haut, voire fait des tours complets, ce qui est impossible dans la réalité. La méthode précise "
                   "(RK4, vert) reste fidèle.";
        case Level::Lycee:
            return "Moment cinétique : m L² θ'' = −m g L sin θ − b θ', soit θ'' = −(g / L) sin θ sans frottement.\n"
                   "Petits angles : sin θ ≈ θ, donc θ'' = −ω0² θ, l'oscillateur harmonique de M2 avec ω0 = √(g / L). "
                   "Pour un grand écart l'approximation est fausse : le vrai pendule est plus lent (la courbe grise, "
                   "petits angles, se décale de la bleue).\n"
                   "Correction d'amplitude : T ≈ T0 (1 + θ0² / 16).\n"
                   "Énergie : E = ½ m L² θ'² + m g L (1 − cos θ).";
        case Level::Etudiant:
            return "Lagrangien L = ½ m L² θ'² + m g L cos θ ; Hamiltonien H = p² / (2 m L²) − m g L cos θ, p = m L² θ'.\n"
                   "Portrait de phase (θ, θ') : orbites fermées (libration) à l'intérieur de la séparatrice "
                   "θ' = ±2 ω0 cos(θ/2) (E = 2 m g L), orbites ouvertes (rotation) à l'extérieur. Essayez une grande vitesse "
                   "initiale.\n"
                   "Solution exacte (lâché sans vitesse) : sin(θ/2) = k sn(K + ω0 t, k), k = sin(θ0/2), T = 4 K(k) / ω0 "
                   "(K : intégrale elliptique complète de première espèce).\n"
                   "H est séparable : Verlet et Euler symplectique gardent une énergie bornée.";
        case Level::Chercheur:
            return "Série de Fourier de sn, nome q = e^(−π K'/K) : k sn(K + ω0 t) = (2π / K) Σ (−1)ⁿ q^(n+½) / (1 − q^(2n+1)) "
                   "cos((2n+1) π ω0 t / 2K), convergence géométrique de raison q².\n"
                   "T(θ0) = T0 · 2 K(sin(θ0/2)) / π ; T → ∞ quand θ0 → π (séparatrice : point selle instable).\n"
                   "Avec frottement ou vitesse initiale il n'y a pas de solution élémentaire : la référence est un RK45 de "
                   "tolérance 1e-13.\n"
                   "H = p²/2mL² + V(θ) est séparable : Verlet est symplectique et conserve exactement un hamiltonien modifié "
                   "H_mod = H + O(dt²) (énergie bornée sur des temps très longs), tandis que RK4 dissipe lentement.";
    }
    return "";
}

void PendulumModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 3.0f;
    camera.target[2] = 3.0f;
    camera.distance = 19.0f;
    camera.yaw = 0.5f;
    camera.pitch = 0.3f;
}

// ------------------------------ simulation -----------------------------

void PendulumModule::reset() {
    problem_.length = length_;
    problem_.mass = mass_;
    problem_.gravity = gravity_;
    problem_.damping = damping_;
    problem_.theta0 = thetaDeg_ * kPi / 180.0;
    problem_.angularVelocity0 = omegaInit_;
    rhs_ = problem_.rhs();
    endTime_ = durationPeriods_ * problem_.period();

    for (Run& r : runs_) {
        r.y = problem_.initialState();
        r.theta.clear();
        r.phase.clear();
        r.energy.clear();
        r.error.clear();
        r.diverged = false;
    }
    refY_ = problem_.initialState();
    thetaRef_.clear();
    phaseRef_.clear();
    energyRef_.clear();
    harmonic_.clear();
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void PendulumModule::referenceState(double& theta, double& omega) const {
    if (problem_.hasExactSolution()) {
        problem_.exact(clock_.time, theta, omega);
    } else {
        theta = refY_[0];
        omega = refY_[1];
    }
}

void PendulumModule::sample() {
    const double t = clock_.time;
    double thRef, omRef;
    referenceState(thRef, omRef);
    thetaRef_.add(t, thRef);
    phaseRef_.add(thRef, omRef);
    energyRef_.add(t, std::max(problem_.energy(thRef, omRef), 1e-12));  // plancher : l'axe est logarithmique

    // Approximation des petits angles : l'oscillateur harmonique linéarisé (réutilise M2), avec frottement éventuel.
    OscillatorProblem linear;
    linear.mass = 1.0;
    linear.stiffness = problem_.omega0() * problem_.omega0();
    linear.damping = problem_.gamma();
    linear.x0 = problem_.theta0;
    linear.v0 = problem_.angularVelocity0;
    harmonic_.add(t, linear.position(t));

    const double w0 = problem_.omega0();
    for (Run& r : runs_) {
        if (r.diverged) continue;
        r.theta.add(t, r.y[0]);
        r.phase.add(r.y[0], r.y[1]);
        r.energy.add(t, std::max(problem_.energy(r.y[0], r.y[1]), 1e-12));
        const double dTheta = r.y[0] - thRef, dOmega = (r.y[1] - omRef) / w0;
        r.error.add(t, std::max(std::sqrt(dTheta * dTheta + dOmega * dOmega), 1e-12));
    }
}

void PendulumModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;

    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    clock_.advance(
        frameSeconds, timeScale_, dt_, endTime_,
        [&](double h) {
            for (int i = 0; i < SolverSet::kCount; ++i) {
                Run& r = runs_[i];
                if (r.diverged) continue;
                advance(solvers_.solver(i), rhs_, clock_.time, r.y, h);
                if (!std::isfinite(r.y[0]) || std::abs(r.y[0]) > 1e6) r.diverged = true;
            }
            if (!problem_.hasExactSolution()) advance(refSolver_, rhs_, clock_.time, refY_, h);
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur finale de chaque solveur à pas fixe pour 7 pas différents, puis pente log-log (ordre mesuré).
void PendulumModule::computeConvergence() {
    static const int kStepCounts[] = {400, 800, 1600, 3200, 6400, 12800, 25600};  // le grand angle exige de petits pas
    const double tEnd = 2.7 * problem_.period();  // pas un multiple de la période : voir oscillatorError

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        ConvergenceCurve& c = convergence_[i];
        c.dt.clear();
        c.err.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        for (int n : kStepCounts) {
            const double e = pendulumError(problem_, *solver, n, tEnd);
            if (e > 1e-11) {  // sous ce seuil on mesure l'arrondi et la précision de la référence
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

void PendulumModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    bool changed = false;

    pushSliderWidth();
    ImGui::SeparatorText("Pendule");
    changed |= sliderD(college ? "Longueur L (m)" : "Longueur du fil", &length_, 0.2, 3.0, college ? "%.2f" : "");
    changed |= sliderD(college ? "Écart initial (deg)" : "Écart initial", &thetaDeg_, 2.0, 175.0, college ? "%.0f" : "");
    if (atLeast(level, Level::Interesse))
        changed |= sliderD(college ? "Masse m (kg)" : "Masse", &mass_, 0.1, 10.0, college ? "%.2f" : "", ImGuiSliderFlags_Logarithmic);

    if (level == Level::Interesse) {
        static const char* planets[] = {"Terre", "Lune", "Mars", "Jupiter"};
        static const double gs[] = {constants::g0, 1.62, 3.71, 24.79};
        static int planet = 0;
        if (ImGui::Combo("Planète", &planet, planets, 4)) {
            gravity_ = gs[planet];
            changed = true;
        }
    }
    if (college) changed |= sliderD("Pesanteur g (m/s²)", &gravity_, 0.5, 30.0, "%.2f");
    if (lycee) changed |= sliderD("Frottement b (N·m·s)", &damping_, 0.0, 1.5, "%.3f");
    if (atLeast(level, Level::Etudiant)) changed |= sliderD("Vitesse initiale (rad/s)", &omegaInit_, -10.0, 10.0, "%.2f");
    popSliderWidth();

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
        double th, om;
        referenceState(th, om);
        ImGui::Text("t = %.2f s / %.2f s", clock_.time, endTime_);
        ImGui::Text("θ = %.1f°   θ' = %.3f rad/s", th * 180.0 / kPi, om);

        ImGui::SeparatorText("Grandeurs");
        ImGui::Text("T0 = 2π√(L/g) = %.3f s", problem_.smallAnglePeriod());
        if (problem_.hasExactSolution()) {
            ImGui::Text("T exacte = %.3f s  (x %.3f)", problem_.period(), problem_.period() / problem_.smallAnglePeriod());
        }
        if (lycee) ImGui::Text("ω0 = √(g/L) = %.3f rad/s", problem_.omega0());
        if (atLeast(level, Level::Etudiant)) {
            const double e0 = problem_.energy(problem_.theta0, problem_.angularVelocity0);
            const double eSeparatrix = 2.0 * problem_.mass * problem_.gravity * problem_.length;
            ImGui::Text("E0 = %.3f J   séparatrice : %.3f J", e0, eSeparatrix);
            ImGui::Text("Mouvement : %s", e0 < eSeparatrix ? "libration (va et vient)" : "rotation (tours complets)");
        }
        if (atLeast(level, Level::Chercheur) && problem_.hasExactSolution()) {
            const double k = std::sin(0.5 * std::abs(problem_.theta0));
            ImGui::Text("k = sin(θ0/2) = %.4f   K(k) = %.4f", k, ellipticK(k));
        }
        ImGui::TextDisabled("Référence : %s", problem_.hasExactSolution() ? "solution exacte (elliptique)" : "RK45, tolérance 1e-13");
    }
}

void PendulumModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant le mouvement, l'énergie passe sans cesse de la hauteur (en haut) au mouvement (en bas). "
                           "Sans frottement, leur total reste le même.");
        return;
    }

    double thRef, omRef;
    referenceState(thRef, omRef);
    const double exactEnergy = problem_.energy(thRef, omRef);
    const double e0 = problem_.energy(problem_.theta0, problem_.angularVelocity0);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool researcher = atLeast(level, Level::Chercheur);
    const double w0 = problem_.omega0();

    std::vector<std::string> headers;
    if (!researcher) headers.push_back("E (J)");
    if (lycee) headers.push_back("dE (J)");
    if (researcher) headers.push_back("dE/E0");
    headers.push_back("dθ (rad)");

    std::vector<TableRow> rows;
    {
        TableRow ref{problem_.hasExactSolution() ? "Exacte" : "Référence", kBlue, {}};
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
            const double dTheta = r.y[0] - thRef, dOmega = (r.y[1] - omRef) / w0;
            if (!researcher) row.cells.push_back(strf("%.3f", energy));
            if (lycee) row.cells.push_back(strf("%+.1e", energy - exactEnergy));
            if (researcher) row.cells.push_back(strf("%+.1e", (energy - exactEnergy) / e0));
            row.cells.push_back(strf("%.1e", std::sqrt(dTheta * dTheta + dOmega * dOmega)));
        }
        rows.push_back(row);
    }
    drawResultTable("erreurs", headers, rows);

    if (lycee) {
        ImGui::SeparatorText("Référence");
        const double kinetic = 0.5 * problem_.mass * problem_.length * problem_.length * omRef * omRef;
        ImGui::Text("Ec = ½ m L² θ'² = %.3f J", kinetic);
        ImGui::Text("Ep = m g L (1 − cos θ) = %.3f J", exactEnergy - kinetic);
        ImGui::TextDisabled(damping_ > 0.0 ? "Le frottement dissipe E." : "Sans frottement : E constante.");
    }
    if (researcher) {
        ImGui::SeparatorText("RK45");
        ImGui::Text("%d pas acceptés, %d rejetés", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps());
        ImGui::Text("%d évaluations de f", solvers_.rk45().evaluations());
    }
}

void PendulumModule::drawGraphs(const UiContext& ctx) {
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
        if (ImPlot::BeginPlot("Angle θ(t)")) {
            ImPlot::SetupAxes("t (s)", "θ (rad)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            const char* refName = problem_.hasExactSolution() ? "Exacte" : "Référence";
            if (thetaRef_.size() > 0)
                ImPlot::PlotLine(refName, thetaRef_.x.data(), thetaRef_.y.data(), thetaRef_.size(), lineSpec(kBlue, thetaRef_.offset));
            if (atLeast(level, Level::Lycee) && harmonic_.size() > 0) {
                const float grey[3] = {0.6f, 0.6f, 0.65f};
                ImPlot::PlotLine("Petits angles", harmonic_.x.data(), harmonic_.y.data(), harmonic_.size(),
                                 lineSpec(grey, harmonic_.offset));
            }
            plotRuns(&Run::theta);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie E(t)")) {
            ImPlot::SetupAxes("t (s)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (atLeast(level, Level::Lycee)) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (energyRef_.size() > 0)
                ImPlot::PlotLine("Référence", energyRef_.x.data(), energyRef_.y.data(), energyRef_.size(),
                                 lineSpec(kBlue, energyRef_.offset));
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Erreur (θ, θ'/ω0)")) {
            ImPlot::SetupAxes("t (s)", "rad", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void PendulumModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showConvergence = atLeast(level, Level::Etudiant);
    if (showConvergence && convergenceDirty_) computeConvergence();
    const int cols = 2 + (showConvergence ? 1 : 0);

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        // Portrait de phase, avec la séparatrice qui sépare les va-et-vient des tours complets.
        if (ImPlot::BeginPlot("Portrait de phase (θ, θ')")) {
            ImPlot::SetupAxes("θ (rad)", "θ' (rad/s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (atLeast(level, Level::Etudiant)) {
                std::vector<double> xs, up, down;
                for (int i = 0; i <= 200; ++i) {
                    const double th = -kPi + 2.0 * kPi * i / 200.0;
                    xs.push_back(th);
                    up.push_back(2.0 * problem_.omega0() * std::cos(0.5 * th));
                    down.push_back(-up.back());
                }
                const float grey[3] = {0.55f, 0.55f, 0.6f};
                ImPlot::PlotLine("Séparatrice", xs.data(), up.data(), static_cast<int>(xs.size()), lineSpec(grey, 0));
                ImPlot::PlotLine("##sep2", xs.data(), down.data(), static_cast<int>(xs.size()), lineSpec(grey, 0));
            }
            if (phaseRef_.size() > 0)
                ImPlot::PlotLine(problem_.hasExactSolution() ? "Exacte" : "Référence", phaseRef_.x.data(), phaseRef_.y.data(),
                                 phaseRef_.size(), lineSpec(kBlue, phaseRef_.offset));
            for (int i = 0; i < SolverSet::kCount; ++i) {
                if (!solvers_.isShown(level, i) || runs_[i].phase.size() == 0) continue;
                const Series& s = runs_[i].phase;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(),
                                 lineSpec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }

        // La période dépend de l'amplitude : T / T0 = 2 K(sin(theta0/2)) / pi, croissante, infinie à theta0 = pi.
        if (ImPlot::BeginPlot("Période selon l'amplitude")) {
            ImPlot::SetupAxes("θ0 (deg)", "T / T0", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            std::vector<double> xs, ys;
            for (int deg = 1; deg <= 178; ++deg) {
                xs.push_back(deg);
                ys.push_back(2.0 * ellipticK(std::sin(0.5 * deg * kPi / 180.0)) / kPi);
            }
            ImPlot::PlotLine("2K(k)/π", xs.data(), ys.data(), static_cast<int>(xs.size()), lineSpec(kBlue, 0));
            const double a = std::min(std::abs(thetaDeg_), 178.0);
            const double xNow = a, yNow = 2.0 * ellipticK(std::sin(0.5 * a * kPi / 180.0)) / kPi;
            const float orange[3] = {1.0f, 0.55f, 0.15f};
            const ImPlotSpec marker(ImPlotProp_LineColor, toImVec4(orange), ImPlotProp_Marker, ImPlotMarker_Circle,
                                    ImPlotProp_MarkerSize, 7.0f);
            ImPlot::PlotScatter("Actuel", &xNow, &yNow, 1, marker);
            ImPlot::EndPlot();
        }

        if (showConvergence && ImPlot::BeginPlot("Convergence : erreur à t = 2,7 T")) {
            ImPlot::SetupAxes("dt (s)", "erreur", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
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

void PendulumModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    // Couloir 0 : la référence ; puis un couloir par méthode affichée.
    double thRef, omRef;
    referenceState(thRef, omRef);
    drawPendulumLane(renderer, 0, thRef, length_, kBlue, ctx.uiScale);
    int lane = 1;
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        const Run& r = runs_[i];
        const float red[3] = {0.9f, 0.2f, 0.2f};
        drawPendulumLane(renderer, lane++, r.diverged ? 0.0 : r.y[0], length_, r.diverged ? red : solvers_.color(i), ctx.uiScale);
    }
}

}  // namespace pl
