#include "NBodyModule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include <implot.h>

#include "physicslab/mechanics/DoublePendulum.hpp"  // lyapunovExponent (ajustement exponentiel générique)

namespace pl {
namespace {

constexpr float kScale = 4.0f;               // 1 unité de longueur = 4 unités de scène
constexpr double kHorizonThreshold = 0.1;    // une méthode est « perdue » quand son écart à la référence dépasse 0,1
constexpr double kChaosAmplification = 1000.0;  // amplification du jumeau au-delà de laquelle on parle de chaos
constexpr std::size_t kTrailLength = 900;
constexpr int kMaxPairs = 6;
const float kTwinColor[3] = {0.55f, 0.90f, 1.0f};
const float kGrey[3] = {0.55f, 0.55f, 0.6f};

// Une couleur par corps (la référence), distinctes de celles des solveurs.
const float kBodyColors[10][3] = {{1.00f, 0.45f, 0.35f}, {0.40f, 0.85f, 0.55f}, {0.45f, 0.62f, 1.00f}, {1.00f, 0.85f, 0.30f},
                                  {0.85f, 0.55f, 0.95f}, {0.30f, 0.90f, 0.90f}, {1.00f, 0.65f, 0.25f}, {0.70f, 0.90f, 0.30f},
                                  {0.95f, 0.45f, 0.70f}, {0.75f, 0.75f, 0.80f}};

// Position (x, y, z) du corps i -> scène : le plan du mouvement (x, y) devient le plan horizontal, z monte.
Vertex pt(const State& y, int i, const float* c, float dim = 1.0f) {
    return {static_cast<float>(y[3 * i]) * kScale, static_cast<float>(y[3 * i + 2]) * kScale,
            -static_cast<float>(y[3 * i + 1]) * kScale, c[0] * dim, c[1] * dim, c[2] * dim};
}

bool isFinite(const State& y) {
    for (double v : y) if (!std::isfinite(v)) return false;
    return true;
}

}  // namespace

NBodyModule::NBodyModule() {
    refSolver_.relTol = twinSolver_.relTol = 1e-13;
    refSolver_.absTol = twinSolver_.absTol = 1e-15;
    selectScenario(kEight);
}

const char* NBodyModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Deux astres qui s'attirent, Newton sait prévoir leur danse pour toujours. À trois, plus personne ne sait "
                   "écrire la réponse : il faut calculer pas à pas.\n\n"
                   "Mais il existe des cas magiques. Le « huit » : trois étoiles identiques se poursuivent sur la même courbe "
                   "en forme de 8, sans jamais se heurter. Trouvé par ordinateur en 1993.\n\n"
                   "Essayez « Amas au hasard » : le mouvement devient imprévisible.";
        case Level::Interesse:
            return "Chaque astre est attiré par tous les autres : ce n'est plus une simple orbite. Vers 1890, Henri Poincaré a "
                   "montré qu'on ne peut pas résoudre le problème des trois corps en général, et qu'il est chaotique : un "
                   "écart minuscule au départ devient énorme plus tard (effet papillon).\n\n"
                   "Les rares solutions exactes : le triangle de Lagrange (1772), où trois astres tournent en gardant la "
                   "forme d'un triangle équilatéral (mais pour trois masses égales il est instable : regardez le jumeau), "
                   "et le huit.";
        case Level::College:
            return "Gravitation universelle : chaque paire de corps s'attire avec F = G × m1 × m2 / d². Le corps 1 subit la "
                   "somme des attractions de tous les autres.\n"
                   "Ici les unités sont « normalisées » (G = 1, masses égales à 1) : seuls les rapports comptent. Le huit "
                   "se répète tous les T = 6,326.\n\n"
                   "La méthode simple (Euler, orange) crée de l'énergie : le huit se déforme vite. La méthode précise "
                   "(RK4, vert) le suit longtemps.";
        case Level::Lycee:
            return "Accélération du corps i : a_i = G Σ_j m_j (r_j − r_i) / |r_j − r_i|³ (somme vectorielle sur tous les autres).\n"
                   "Trois grandeurs se conservent : l'impulsion totale P = Σ m v (le centre de masse avance en ligne "
                   "droite, ici il reste immobile, 3e loi de Newton), le moment cinétique L = Σ m r × v, et l'énergie "
                   "E = Σ ½ m v² − Σ_{i<j} G m_i m_j / r_ij.\n"
                   "Le « jumeau » (cyan) démarre avec un écart de 1e-9 : s'il s'écarte exponentiellement de la référence, "
                   "le mouvement est chaotique ; si l'écart croît seulement proportionnellement au temps, il est stable.";
        case Level::Etudiant:
            return "H = Σ p_i²/(2 m_i) − Σ_{i<j} G m_i m_j / r_ij. Le théorème de Noether relie symétries et lois de "
                   "conservation : translation d'espace → P, rotation → L, translation du temps → E. 10 intégrales "
                   "premières classiques (P : 3, centre de masse : 3, L : 3, E : 1) ; Bruns (1887) et Poincaré (1890) ont montré "
                   "qu'il n'y en a pas d'autre « simple » : pas de solution générale à N ≥ 3.\n"
                   "Euler symplectique et Verlet conservent L exactement (à l'arrondi) et bornent E ; RK4 conserve P (invariant "
                   "linéaire) mais pas L (quadratique) ; Euler ne conserve rien.\n"
                   "Triangle de Lagrange : ω² = G M / s³ ; stable seulement si (Σm)² > 27 Σ_{i<j} m_i m_j (critère de Routh), "
                   "donc instable à masses égales.";
        case Level::Chercheur:
            return "Adoucissement de Plummer ε : a_i = G Σ_j m_j (r_j − r_i)/(r_ij² + ε²)^(3/2), dérivée de l'énergie "
                   "U = −G Σ m_i m_j/√(r_ij² + ε²) : E reste exactement conservée (une force tronquée ne la conserverait "
                   "pas). Boucle sur les paires i < j, la même force appliquée aux deux corps : P conservée à l'arrondi (1e-15).\n"
                   "Coût O(N²) par évaluation ; RK4 en demande 4, Verlet 1 : à précision d'énergie égale Verlet est bien moins "
                   "cher. Pas fixe et rencontres rapprochées : l'erreur d'énergie y saute (il faudrait un pas adaptatif, une "
                   "régularisation ou un pas partagé par particule).\n"
                   "Le jumeau mesure le plus grand exposant de Lyapunov λ ; le temps de prédictibilité vaut ~ ln(Δ/d0)/λ. En "
                   "simple précision (GPU, M7) d0 ≈ 1e-7 : l'horizon est plus court, à comparer avec ce calcul CPU en double.";
    }
    return "";
}

void NBodyModule::frameCamera(Camera& camera) const {
    camera.target[0] = camera.target[1] = camera.target[2] = 0.0f;
    camera.distance = 15.0f;
    camera.yaw = 0.0f;
    camera.pitch = 1.0f;
}

// ------------------------------ simulation -----------------------------

// Valeurs par défaut de chaque scénario (durée de 3 périodes pour les mouvements périodiques).
void NBodyModule::selectScenario(Scenario s) {
    scenario_ = s;
    mass3Factor_ = 1.0;
    switch (s) {
        case kEight:
            durationSec_ = 3.0 * nbody::kFigureEightPeriod;
            dt_ = 0.005;
            convergenceTime_ = 0.35 * nbody::kFigureEightPeriod;
            break;
        case kLagrange:
            durationSec_ = 3.0 * 2.0 * constants::pi / std::sqrt(3.0);
            dt_ = 0.005;
            convergenceTime_ = 0.35 * 2.0 * constants::pi / std::sqrt(3.0);
            break;
        case kCluster:
            durationSec_ = 20.0;
            dt_ = 0.002;
            convergenceTime_ = 2.0;
            break;
        default:
            break;
    }
    reset();
}

void NBodyModule::buildProblem() {
    switch (scenario_) {
        case kEight: problem_ = NBodyProblem::figureEight(); break;
        case kLagrange: problem_ = NBodyProblem::lagrangeTriangle(); break;
        default: problem_ = NBodyProblem::randomCluster(static_cast<int>(std::lround(clusterCount_)), seed_, 1.0, softening_); break;
    }
    if (scenario_ == kEight && mass3Factor_ != 1.0) {
        // Masses inégales : la chorégraphie n'existe plus. On recentre le système pour garder P = 0.
        problem_.mass[2] = mass3Factor_;
        Vec3 c, v;
        for (int i = 0; i < 3; ++i) { c += problem_.mass[i] * problem_.position[i]; v += problem_.mass[i] * problem_.velocity[i]; }
        c /= problem_.totalMass();
        v /= problem_.totalMass();
        for (int i = 0; i < 3; ++i) { problem_.position[i] -= c; problem_.velocity[i] -= v; }
    }
}

void NBodyModule::reset() {
    buildProblem();
    rhs_ = problem_.rhs();
    const int n = problem_.count();
    const State y0 = problem_.initialState();
    e0_ = problem_.energy(y0);
    l0_ = problem_.angularMomentum(y0);

    for (Run& r : runs_) {
        r.y = y0;
        r.energy.clear();
        r.error.clear();
        r.angular.clear();
        r.trails.assign(n, {});
        r.diverged = false;
        r.horizon = -1.0;
    }
    refY_ = y0;
    twinY_ = y0;
    twinY_[0] += perturbation_;
    refTrails_.assign(n, {});
    bodyX_.assign(std::min(n, 6), Series{});
    pairDistance_.assign(std::min(n * (n - 1) / 2, kMaxPairs), Series{});
    minDistance_.clear();
    twinT_.clear();
    twinD_.clear();
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void NBodyModule::sample() {
    const double t = clock_.time;
    const int n = problem_.count();

    for (int i = 0; i < n; ++i) {
        refTrails_[i].push_back(pt(refY_, i, kBodyColors[i % 10], 0.6f));
        if (refTrails_[i].size() > kTrailLength) refTrails_[i].erase(refTrails_[i].begin(), refTrails_[i].begin() + 100);
    }
    for (std::size_t i = 0; i < bodyX_.size(); ++i) bodyX_[i].add(t, refY_[3 * i]);

    int pair = 0;
    double minD = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d = std::sqrt(std::pow(refY_[3 * i] - refY_[3 * j], 2) + std::pow(refY_[3 * i + 1] - refY_[3 * j + 1], 2) +
                                       std::pow(refY_[3 * i + 2] - refY_[3 * j + 2], 2));
            minD = std::min(minD, d);
            if (pair < static_cast<int>(pairDistance_.size())) pairDistance_[pair++].add(t, d);
        }
    }
    if (std::isfinite(minD)) minDistance_.add(t, minD);
    twinT_.push_back(t);
    twinD_.push_back(std::max(problem_.distance(refY_, twinY_), 1e-300));

    for (int i = 0; i < SolverSet::kCount; ++i) {
        Run& r = runs_[i];
        if (r.diverged) continue;
        r.energy.add(t, (problem_.energy(r.y) - e0_) / std::abs(e0_));
        const double d = problem_.distance(r.y, refY_);
        r.error.add(t, std::max(d, 1e-12));
        r.angular.add(t, std::max((problem_.angularMomentum(r.y) - l0_).norm(), 1e-18));
        if (r.horizon < 0.0 && d > kHorizonThreshold) r.horizon = t;
        for (int b = 0; b < n; ++b) {
            r.trails[b].push_back(pt(r.y, b, solvers_.color(i), 0.7f));
            if (r.trails[b].size() > kTrailLength) r.trails[b].erase(r.trails[b].begin(), r.trails[b].begin() + 100);
        }
    }
}

void NBodyModule::update(double frameSeconds) {
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
                if (!isFinite(r.y)) r.diverged = true;
            }
            advance(refSolver_, rhs_, clock_.time, refY_, h);
            advance(twinSolver_, rhs_, clock_.time, twinY_, h);
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur de chaque solveur à pas fixe à `convergenceTime_` (avant que le chaos n'amplifie tout), pour différents pas.
void NBodyModule::computeConvergence() {
    static const int kStepCounts[] = {100, 200, 400, 800, 1600, 3200, 6400};
    const State reference = problem_.reference(convergenceTime_);
    const OdeFunction f = problem_.rhs();

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            State y = problem_.initialState();
            const double dt = convergenceTime_ / n;
            double t = 0.0;
            for (int k = 0; k < n; ++k) t += advance(*solver, f, t, y, dt);
            const double e = problem_.distance(y, reference);
            if (!std::isfinite(e) || e <= 1e-12) continue;  // sous ce seuil : précision de la référence
            c.x.push_back(dt);
            c.y.push_back(e);
            if (e < 0.1) {  // pente ajustée là où la méthode est dans son régime asymptotique
                const double lx = std::log10(dt), ly = std::log10(e);
                sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
                ++m;
            }
        }
        c.slope = std::numeric_limits<double>::quiet_NaN();  // hors du régime asymptotique
        if (m >= 2) {
            const double dm = static_cast<double>(m);
            c.slope = (dm * sxy - sx * sy) / (dm * sxx - sx * sx);
        }
    }
    convergenceDirty_ = false;
}

// --------------------------------- UI ----------------------------------

void NBodyModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    bool changed = false;

    ImGui::SeparatorText("Scénario");
    int s = scenario_;
    ImGui::RadioButton("Le « huit » (3 corps égaux)", &s, kEight);
    ImGui::RadioButton("Triangle de Lagrange", &s, kLagrange);
    ImGui::RadioButton("Amas au hasard", &s, kCluster);
    if (s != scenario_) {
        selectScenario(static_cast<Scenario>(s));
        return;  // l'état vient d'être réinitialisé : on redessine au prochain tour
    }

    if (scenario_ == kCluster) {
        if (interesse && ImGui::Button("Nouveau tirage")) {
            seed_ = seed_ * 1664525u + 1013904223u;  // graine suivante (congruentiel linéaire), reproductible
            changed = true;
        }
        pushSliderWidth();
        if (college) changed |= sliderD("Nombre de corps", &clusterCount_, 3.0, 10.0, "%.0f");
        popSliderWidth();
    }
    if (college && scenario_ == kEight) {
        pushSliderWidth();
        changed |= sliderD("Masse du corps 3 (x)", &mass3Factor_, 0.5, 2.0, "%.2f");
        popSliderWidth();
        if (mass3Factor_ != 1.0) ImGui::TextDisabled("Masses inégales : le huit n'existe plus.");
    }

    if (college) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt", &dt_, 1e-4, 0.05, "%.4f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Lycee)) changed |= sliderD("Durée", &durationSec_, 2.0, 60.0, "%.1f");
        if (atLeast(level, Level::Etudiant) && scenario_ == kCluster)
            changed |= sliderD("Adoucissement ε", &softening_, 0.005, 0.3, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Chercheur)) {
            changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
            changed |= sliderD("Écart du jumeau", &perturbation_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
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
    if (interesse) sliderD("Vitesse du temps", &timeScale_, 0.05, 3.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        ImGui::Text("t = %.2f / %.1f", clock_.time, durationSec_);
        if (scenario_ == kEight && mass3Factor_ == 1.0) ImGui::Text("(%.2f périodes du huit)", clock_.time / nbody::kFigureEightPeriod);
        ImGui::Text("E = %.4f", problem_.energy(refY_));
    }
    if (atLeast(level, Level::Lycee)) {
        // Verdict fondé sur l'amplification réelle de l'écart : le huit (stable) ne fait croître l'écart que linéairement en t
        // (x ~ 70 en 19 unités de temps), un amas chaotique le multiplie par des milliers et plus.
        const double amplification = twinD_.empty() ? 1.0 : twinD_.back() / std::max(twinD_.front(), 1e-300);
        ImGui::SeparatorText("Chaos");
        ImGui::Text("Écart du jumeau : %.1e (amplification x %.1e)", twinD_.empty() ? 0.0 : twinD_.back(), amplification);
        if (amplification >= kChaosAmplification) {
            int used = 0;
            const double lambda = lyapunovExponent(twinT_, twinD_, 20.0 * perturbation_, 0.1, &used);
            if (std::isfinite(lambda)) ImGui::Text("λ = %.3f   (%d points)", lambda, used);
            ImGui::TextDisabled("Chaotique : l'écart a été multiplié par plus de %.0f.", kChaosAmplification);
        } else {
            ImGui::TextDisabled("Pas d'explosion de l'écart (croissance lente : système stable, ou trop tôt).");
        }
    }
}

void NBodyModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant que les astres bougent, l'énergie totale, la « quantité de mouvement » et la « quantité de "
                           "rotation » de l'ensemble ne changent pas. C'est ce qui permet de vérifier que l'ordinateur ne "
                           "se trompe pas.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const State y0 = problem_.initialState();

    ImGui::SeparatorText("Système");
    wrapped(strf("%d corps, masse totale %.3g ; G = 1", problem_.count(), problem_.totalMass()));
    wrapped(strf("E = %.5f   (cinétique %.4f, potentielle %.4f)", e0_, problem_.kineticEnergy(y0), problem_.potentialEnergy(y0)));
    if (lycee) {
        wrapped(strf("|P| = %.1e   |L| = %.4f", problem_.momentum(y0).norm(), l0_.norm()));
        if (scenario_ == kEight && mass3Factor_ == 1.0) wrapped(strf("Période du huit T = %.8f", nbody::kFigureEightPeriod));
        if (scenario_ == kLagrange) wrapped(strf("Période de rotation T = 2π/ω = %.4f  (ω² = 3 G m / s³)", 2.0 * constants::pi / std::sqrt(3.0)));
        if (problem_.softening > 0.0) wrapped(strf("Adoucissement ε = %.3g", problem_.softening));
    }

    const bool split = etudiant;
    std::vector<std::string> headersA, headersB;
    headersA.push_back(etudiant ? "ΔE/E0" : "ΔE %");
    if (etudiant) { headersA.push_back("|P|"); headersA.push_back("|ΔL|"); }
    (split ? headersB : headersA).push_back(etudiant ? "dist" : "écart");
    if (lycee) (split ? headersB : headersA).push_back("t_H");

    std::vector<TableRow> rowsA, rowsB;
    rowsA.push_back({"Référence", kBlue, std::vector<std::string>(headersA.size(), "-")});
    if (split) rowsB.push_back({"Référence", kBlue, std::vector<std::string>(headersB.size(), "-")});

    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(level, i)) continue;
        const Run& r = runs_[i];
        TableRow a{solvers_.shortLabel(level, i), solvers_.color(i), {}};
        TableRow b = a;
        if (r.diverged) {
            a.cells.assign(headersA.size(), "-");
            a.cells[0] = "diverge";
            b.cells.assign(headersB.size(), "-");
        } else {
            const double dE = (problem_.energy(r.y) - e0_) / std::abs(e0_);
            a.cells.push_back(strf("%+.1e", etudiant ? dE : 100.0 * dE));
            if (etudiant) {
                a.cells.push_back(strf("%.1e", problem_.momentum(r.y).norm()));
                a.cells.push_back(strf("%.1e", (problem_.angularMomentum(r.y) - l0_).norm()));
            }
            std::vector<std::string>& prec = (split ? b : a).cells;
            prec.push_back(strf("%.1e", problem_.distance(r.y, refY_)));
            if (lycee) prec.push_back(r.horizon >= 0.0 ? strf("%.2f", r.horizon) : "-");
        }
        rowsA.push_back(a);
        if (split) rowsB.push_back(b);
    }

    ImGui::SeparatorText(split ? "Ce qui se conserve" : "Comparaison");
    drawResultTable("conservation", headersA, rowsA);
    if (split) {
        ImGui::SeparatorText("Précision de la trajectoire");
        drawResultTable("precision", headersB, rowsB);
    }
    wrapped(split ? "dist : écart à la référence dans l'espace des phases (sans unité). t_H : instant où il dépasse 0,1."
                  : "écart : distance à la référence (sans unité). ΔE : variation en %.",
            true);
    wrapped("« - » : la méthode suit encore la référence ; « diverge » : calcul perdu.", true);
}

void NBodyModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showEnergy = atLeast(level, Level::College);
    const bool showError = atLeast(level, Level::College);
    const bool showAngular = atLeast(level, Level::Lycee);
    const int cols = 1 + (showEnergy ? 1 : 0) + (showError ? 1 : 0) + (showAngular ? 1 : 0);

    auto plotRuns = [&](Series Run::*member) {
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i) || (runs_[i].*member).size() == 0) continue;
            const Series& s = runs_[i].*member;
            ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(solvers_.color(i), s.offset));
        }
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Position x de chaque corps")) {
            ImPlot::SetupAxes("t", "x", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            for (std::size_t i = 0; i < bodyX_.size(); ++i) {
                const Series& s = bodyX_[i];
                if (s.size() == 0) continue;
                ImPlot::PlotLine(strf("corps %d", static_cast<int>(i) + 1).c_str(), s.x.data(), s.y.data(), s.size(),
                                 lineSpec(kBodyColors[i % 10], s.offset));
            }
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie (E − E0) / |E0|")) {
            ImPlot::SetupAxes("t", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Écart à la référence")) {
            ImPlot::SetupAxes("t", "distance (espace des phases)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            const double threshold = kHorizonThreshold;
            ImPlot::PlotInfLines("seuil", &threshold, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kGrey), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            ImPlot::EndPlot();
        }
        if (showAngular && ImPlot::BeginPlot("Moment cinétique |L − L0|")) {
            ImPlot::SetupAxes("t", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::angular);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void NBodyModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool etudiant = atLeast(level, Level::Etudiant);
    if (etudiant && convergenceDirty_) computeConvergence();
    const int cols = etudiant ? 3 : 2;

    auto markerSpec = [](const float* c) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle,
                          ImPlotProp_MarkerSize, 4.0f);
    };

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        // Une collision rapprochée est la difficulté numérique : la distance minimale la montre.
        if (ImPlot::BeginPlot("Distances entre les corps (référence)")) {
            ImPlot::SetupAxes("t", "r", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            for (std::size_t k = 0; k < pairDistance_.size(); ++k) {
                const Series& s = pairDistance_[k];
                if (s.size() > 0) ImPlot::PlotLine(strf("paire %d", static_cast<int>(k) + 1).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(kBodyColors[k % 10], s.offset));
            }
            if (problem_.count() > 3 && minDistance_.size() > 0)
                ImPlot::PlotLine("minimale", minDistance_.x.data(), minDistance_.y.data(), minDistance_.size(), lineSpec(kGrey, minDistance_.offset));
            ImPlot::EndPlot();
        }

        // Sensibilité aux conditions initiales : exponentielle (chaos) ou linéaire (stable).
        if (ImPlot::BeginPlot("Sensibilité : écart référence - jumeau")) {
            ImPlot::SetupAxes("t", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (!twinT_.empty())
                ImPlot::PlotLine("Écart", twinT_.data(), twinD_.data(), static_cast<int>(twinT_.size()), lineSpec(kTwinColor, 0));
            int used = 0;
            double intercept = 0.0;
            const double lambda = lyapunovExponent(twinT_, twinD_, 20.0 * perturbation_, 0.1, &used, &intercept);
            const bool exponential = !twinD_.empty() && twinD_.back() / std::max(twinD_.front(), 1e-300) >= kChaosAmplification;
            if (exponential && std::isfinite(lambda)) {  // la droite n'a de sens que si la croissance est vraiment exponentielle
                std::vector<double> xs, ys;
                const double tEnd = std::max(clock_.time, 1.0);
                for (int i = 0; i <= 50; ++i) {
                    const double t = tEnd * i / 50.0;
                    xs.push_back(t);
                    ys.push_back(std::min(std::exp(intercept + lambda * t), 10.0));
                }
                const float orange[3] = {1.0f, 0.55f, 0.15f};
                const std::string name = strf("e^(λt), λ = %.2f", lambda);
                ImPlot::PlotLine(name.c_str(), xs.data(), ys.data(), static_cast<int>(xs.size()), lineSpec(orange, 0));
            }
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot(strf("Convergence : erreur à t = %.2f", convergenceTime_).c_str())) {
            ImPlot::SetupAxes("dt", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < SolverSet::kFixedStep; ++i) {
                const Curve& c = convergence_[i];
                if (!solvers_.show(i) || c.x.size() < 2) continue;
                const std::string name = std::isfinite(c.slope) ? strf("%s (pente %.2f)", solvers_.solver(i).name(), c.slope)
                                                                 : strf("%s (hors régime)", solvers_.solver(i).name());
                ImPlot::PlotLine(name.c_str(), c.x.data(), c.y.data(), static_cast<int>(c.x.size()), markerSpec(solvers_.color(i)));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

void NBodyModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    const int n = problem_.count();
    double massMax = 0.0;
    for (double m : problem_.mass) massMax = std::max(massMax, m);

    // Jumeau (niveau 4+), puis méthodes, puis la référence par-dessus (grosses boules colorées, une par corps).
    if (atLeast(ctx.level, Level::Lycee)) {
        std::vector<Vertex> twin;
        for (int i = 0; i < n; ++i) twin.push_back(pt(twinY_, i, kTwinColor));
        renderer.draw(Primitive::Points, twin, 8.0f * ctx.uiScale);
    }
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        const Run& r = runs_[i];
        if (r.diverged) continue;
        for (const auto& trail : r.trails)
            if (!trail.empty()) renderer.draw(Primitive::LineStrip, trail);
        std::vector<Vertex> bodies;
        for (int b = 0; b < n; ++b) bodies.push_back(pt(r.y, b, solvers_.color(i)));
        renderer.draw(Primitive::Points, bodies, 9.0f * ctx.uiScale);
    }
    for (int b = 0; b < n; ++b) {
        if (!refTrails_[b].empty()) renderer.draw(Primitive::LineStrip, refTrails_[b]);
        const float size = (11.0f + 11.0f * static_cast<float>(std::cbrt(problem_.mass[b] / massMax))) * ctx.uiScale;
        renderer.draw(Primitive::Points, {pt(refY_, b, kBodyColors[b % 10])}, size);
    }
}

}  // namespace pl
