#include "CradleModule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <implot.h>

#include "physicslab/core/Constants.hpp"

namespace pl {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kSlowMotion = 0.08;     // facteur de ralenti pendant le choc
constexpr double kLaneZ = 1.0;           // décalage des deux rangées quand elles sont affichées ensemble
const float kRail[3] = {0.45f, 0.45f, 0.5f};
const float kForceColor[3] = {1.0f, 0.82f, 0.30f};
const float kVelocityColor[3] = {0.35f, 1.0f, 0.45f};
const float kKineticColor[3] = {0.45f, 0.80f, 0.95f};
const float kTotalColor[3] = {0.90f, 0.90f, 0.90f};

// Une couleur par bille (de la gauche vers la droite).
const float kBallColors[7][3] = {{1.00f, 0.45f, 0.35f}, {1.00f, 0.75f, 0.25f}, {0.75f, 0.90f, 0.30f}, {0.35f, 0.85f, 0.55f},
                                 {0.30f, 0.80f, 0.95f}, {0.50f, 0.55f, 1.00f}, {0.85f, 0.50f, 0.95f}};

std::vector<Vertex> circleXZ(double cx, double cz, double r, const float* c, float dim = 1.0f) {
    std::vector<Vertex> v;
    for (int i = 0; i <= 40; ++i) {
        const double a = 2.0 * constants::pi * i / 40.0;
        v.push_back({static_cast<float>(cx + r * std::cos(a)), 0.0f, static_cast<float>(cz + r * std::sin(a)), c[0] * dim, c[1] * dim, c[2] * dim});
    }
    return v;
}

std::string velocityText(const std::vector<double>& v, double scale) {
    std::string s = "(";
    for (std::size_t i = 0; i < v.size(); ++i) s += strf(i ? " ; %+.3f" : "%+.3f", v[i] / scale);
    return s + ")";
}

}  // namespace

CradleModule::CradleModule() {
    ref_.relTol = 1e-12;
    ref_.absTol = 1e-14;
    reset();
}

const char* CradleModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Le berceau de Newton : on lâche une bille, et une seule bille repart de l'autre côté. On en lâche deux, deux repartent. "
                   "Étonnant : les lois de la physique (le mouvement total et l'énergie ne changent pas) autoriseraient d'autres "
                   "sorties, par exemple les deux dernières billes qui partent ensemble, plus lentement, pendant que la première recule.\n\n"
                   "Le secret : chaque bille s'écrase un tout petit peu au contact, comme un ressort très dur, et le choc traverse la rangée "
                   "comme une vague. Devant : les billes qui s'écrasent vraiment (la vague de force est en jaune). Derrière : un calcul "
                   "simplifié où les chocs sont instantanés, un par un. Regardez les premières billes de devant : elles reculent un peu, ce "
                   "que le calcul simplifié ne voit jamais. Le temps est ralenti pendant le choc.";
        case Level::Interesse:
            return "Dans un choc sans perte, deux choses ne changent jamais : la quantité de mouvement totale (masse × vitesse) et l'énergie de "
                   "mouvement. Avec 2 billes elles imposent le résultat : la bille lancée s'arrête, l'autre part. Avec 3 billes ou plus elles "
                   "ne suffisent plus : plusieurs sorties sont possibles. Seule la façon dont les billes se touchent départage.\n\n"
                   "Modèle réel (devant) : les billes se compriment comme un ressort de plus en plus dur, la force de Hertz. Modèle "
                   "simplifié (derrière) : des chocs instantanés, deux billes à la fois. Il donne un résultat différent : la première bille ne "
                   "recule jamais, et avec un rebond imparfait (curseur) de l'énergie se perd.";
        case Level::College:
            return "Quantité de mouvement p = m × v ; énergie de mouvement E = ½ m v². Une bille lancée à la vitesse v : p = m v et E = ½ m v². "
                   "Pourquoi pas deux billes qui repartent à v/2 ? p = 2 × m × v/2 = m v convient, mais E = 2 × ½ m (v/2)² = ¼ m v² : la moitié "
                   "de l'énergie aurait disparu.\n"
                   "Avec 3 billes, p et E ne donnent pourtant pas un seul résultat. Le tableau compare les vitesses (divisées par la vitesse "
                   "de lancement) du modèle réel et des chocs instantanés : leur différence vient de la façon dont le contact est décrit.";
        case Level::Lycee:
            return "Contact de Hertz : deux sphères qui s'écrasent d'une longueur δ se repoussent avec la force F = k δ^(3/2) (nulle si elles ne "
                   "se touchent pas). Chaque bille obéit à m a = F(gauche) − F(droite). Pour un choc de deux billes de vitesse relative v : "
                   "écrasement maximal δ_max = (5 μ v²/(4 k))^(2/5), avec μ = m/2, et durée du contact T = 2,94 δ_max / v.\n"
                   "L'énergie ½ m Σ v² + Σ (2/5) k δ^(5/2) reste constante (le contact est élastique). Les chocs instantanés conservent p mais "
                   "pas E si e < 1. Courbe « Force des contacts » : l'onde de compression qui traverse la rangée.";
        case Level::Etudiant:
            return "N billes : x_i'' = [F_(i−1) − F_i]/m, F_i = k δ_i^(3/2) pour δ_i = 2R − (x_(i+1) − x_i) > 0, énergie élastique U = (2/5) k δ^(5/2). "
                   "Avec l'échelle de longueur L = (m v²/k)^(2/5) et de temps L/v, les équations n'ont plus aucun paramètre : le résultat, en unités de "
                   "v, ne dépend ni de k ni de v (vérifié : essayez les curseurs).\n"
                   "Trois billes : Σ v_i = v et Σ v_i² = v² donnent v1 v2 + v1 v3 + v2 v3 = 0, un cercle de solutions (« Analyse »). Hertz "
                   "choisit x = −v1/v ≈ 0,071 ; les chocs instantanés (e = 1) donnent x = 0.\n"
                   "Schémas : la force n'est pas lisse au début du contact (δ^(3/2)), l'ordre de RK4 tombe vers 2,5 ; Euler symplectique et Verlet "
                   "donnent la même erreur ici (conjugués, hors contact).";
        case Level::Chercheur:
            return "Contact de Hertz (1881) F = k δ^(3/2), valable pour un écrasement δ très petit devant R ; amortissement optionnel de Hunt-Crossley F = k δ^(3/2) (1 + (3/2) α δ'), "
                   "restitution mesurée e ≈ 1/(1 + α v). F est C^1 mais pas C^2 en δ = 0 (début et fin de chaque contact) : le pas du schéma doit "
                   "résoudre le contact (T/100 ici) et l'ordre global de RK4 tombe à ≈ 2,5 (mesuré 2,3 à 2,75 selon la phase du contact sur la grille). "
                   "L'onde est de type solitaire dans une chaîne de Hertz non précomprimée (Nesterenko, 1983).\n"
                   "Impulsions séquentielles : chocs binaires de Newton. Un choc simultané à 3 corps n'existe pas dans ce modèle ; l'ordre de "
                   "résolution ne change rien dans le berceau (la causalité fixe la séquence : 245 cas vérifiés) mais change le résultat dès "
                   "qu'une bille est prise entre deux voisines qui s'approchent (démonstration dans « Invariants »). La dynamique du contact "
                   "lève l'indétermination sans paramètre ad hoc.";
    }
    return "";
}

void CradleModule::frameCamera(Camera& camera) const {
    camera.target[0] = static_cast<float>(0.5 * (viewMin_ + viewMax_));
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = static_cast<float>(0.8 * (viewMax_ - viewMin_) + 1.0);
    camera.yaw = 0.1f;
    camera.pitch = 0.65f;
}

// ------------------------------ simulation -----------------------------

void CradleModule::reset() {
    launched_ = std::clamp(launched_, 1, balls_ - 1);
    prob_ = CradleProblem{};
    prob_.balls = balls_;
    prob_.launched = launched_;
    prob_.stiffness = stiffness_;
    prob_.damping = damping_;
    prob_.speed = speed_;
    prob_.gap = gap_;
    rhs_ = prob_.rhs();
    y_ = yRef_ = prob_.initialState();

    const double mu = 0.5 * prob_.mass;
    contactTime_ = hertz::contactDuration(mu, speed_, stiffness_);
    dt_ = contactTime_ / stepsPerContact_;
    forceScale_ = hertz::force(stiffness_, hertz::maxCompression(mu, speed_, stiffness_));
    energy0_ = prob_.energy(y_);

    // chocs instantanés : vitesses en unités de la vitesse de lancement (le résultat est homogène)
    v0_.assign(balls_, 0.0);
    for (int i = 0; i < launched_; ++i) v0_[i] = 1.0;
    impulses_ = cradle::sequentialImpulses(v0_, restitution_, cradle::ResolveOrder::LeftToRight);

    // issue de référence de Hertz : durée de la collision, tableau et courbe des solutions
    RK45 tight;
    tight.relTol = 1e-12;
    tight.absTol = 1e-14;
    outcome_ = prob_.run(tight, prob_.suggestedStep(), gap_ / speed_ + 60.0);
    const double exitTime = 2.5 * 2.0 * prob_.radius / speed_;   // la bille qui part s'éloigne de 2,5 diamètres
    duration_ = (outcome_.finished ? outcome_.time : gap_ / speed_ + 3.0 * contactTime_) + exitTime;

    viewMin_ = -gap_ - 1.5;
    viewMax_ = 2.0 * prob_.radius * (balls_ - 1) + 1.0 + 2.5 * 2.0 * prob_.radius + 1.5;

    for (Series& s : vel_) s.clear();
    for (Series& s : impVel_) s.clear();
    for (Series& s : force_) s.clear();
    kinetic_.clear();
    potential_.clear();
    total_.clear();
    error_.clear();
    maxCompression_ = 0.0;
    contactStart_ = contactEnd_ = -1.0;
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    diverged_ = false;
    convergenceDirty_ = true;
    sample();
}

double CradleModule::contactForce(const State& y, int pair) const {
    const double delta = prob_.compression(y, pair);
    double f = hertz::force(prob_.stiffness, delta);
    if (prob_.damping > 0.0 && delta > 0.0) f *= std::max(0.0, 1.0 + 1.5 * prob_.damping * (y[balls_ + pair] - y[balls_ + pair + 1]));
    return f;
}

double CradleModule::phaseDistance(const State& a, const State& b) const {
    const double length = speed_ * contactTime_;
    double sum = 0.0;
    for (int i = 0; i < balls_; ++i) {
        const double dx = (a[i] - b[i]) / length, dv = (a[balls_ + i] - b[balls_ + i]) / speed_;
        sum += dx * dx + dv * dv;
    }
    return std::sqrt(sum);
}

double CradleModule::impulseVelocity(int ball, double t) const {
    if (t < gap_ / speed_) return ball < launched_ ? speed_ : 0.0;
    return impulses_.velocities[ball] * speed_;
}

double CradleModule::impulsePosition(int ball, double t) const {
    const double x0 = 2.0 * prob_.radius * ball, tc = gap_ / speed_;
    if (t < tc) return ball < launched_ ? x0 - gap_ + speed_ * t : x0;
    return x0 + impulseVelocity(ball, t) * (t - tc);
}

// Ralenti automatique : tant que la collision n'est pas finie et que la bille de tête est à moins de 0,05 du contact (ou déjà dessus).
bool CradleModule::slowMotion() const {
    if (prob_.collisionOver(y_)) return false;
    return -prob_.compression(y_, launched_ - 1) <= 0.05;
}

void CradleModule::step(double h) {
    const double t = clock_.time;
    advance(solvers_.solver(solverIndex_), rhs_, t, y_, h, 1000);
    advance(ref_, rhs_, t, yRef_, h, 1000);
    for (double v : y_)
        if (!std::isfinite(v) || std::abs(v) > 1e6) diverged_ = true;   // Euler explicite à pas trop grand : l'énergie explose
    if (diverged_) return;
    for (int i = 0; i + 1 < balls_; ++i) maxCompression_ = std::max(maxCompression_, prob_.compression(y_, i));
    const double front = prob_.compression(y_, launched_ - 1), t1 = t + h;
    if (contactStart_ < 0.0) {
        if (front > 0.0) contactStart_ = t1;
    } else if (contactEnd_ < 0.0 && front <= 0.0) {
        contactEnd_ = t1;
    }
}

void CradleModule::sample() {
    if (diverged_) return;
    const double t = clock_.time;
    for (int i = 0; i < balls_; ++i) {
        vel_[i].add(t, y_[balls_ + i] / speed_);
        impVel_[i].add(t, impulseVelocity(i, t) / speed_);
    }
    for (int i = 0; i + 1 < balls_; ++i) force_[i].add(t, contactForce(y_, i) / forceScale_);
    const double k = prob_.kineticEnergy(y_), u = prob_.potentialEnergy(y_);
    kinetic_.add(t, k / energy0_);
    potential_.add(t, u / energy0_);
    total_.add(t, (k + u) / energy0_);
    error_.add(t, std::max(phaseDistance(y_, yRef_), 1e-16));
}

void CradleModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;
    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    clock_.advance(
        frameSeconds, timeScale_ * (slowMotion() ? kSlowMotion : 1.0), dt_, duration_, [&](double h) { step(h); }, [&] { sample(); });
    if (clock_.finished || diverged_) running_ = false;
}

// Erreur de chaque solveur à pas fixe pour différents pas. Contact décalé de 0,0123 (jamais sur la grille des pas) comme dans les tests.
void CradleModule::computeConvergence() {
    static const int kStepCounts[] = {50, 100, 200, 400, 800, 1600, 3200};
    CradleProblem q = prob_;
    q.gap = 0.0123 * speed_;
    RK45 tight;
    tight.relTol = 1e-12;
    tight.absTol = 1e-14;
    const CradleOutcome ref = q.run(tight, q.suggestedStep(), q.gap / speed_ + 60.0);
    convergenceTime_ = (ref.finished ? ref.time : q.gap / speed_ + 3.0 * contactTime_) + 0.5 * contactTime_;

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            const double e = cradleError(q, *solver, n, convergenceTime_);
            if (!std::isfinite(e) || e <= 1e-12) continue;
            const double dt = convergenceTime_ / n;
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

// --------------------------------- UI ----------------------------------

void CradleModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (!etudiant && solverIndex_ != SolverSet::kRK4) {   // l'ordinateur « précis » par défaut tant que le choix du schéma n'est pas exposé
        solverIndex_ = SolverSet::kRK4;
        changed = true;
    }

    if (interesse) {
        ImGui::TextDisabled("Cas types");
        if (ImGui::SmallButton("1 sur 4")) { balls_ = 5; launched_ = 1; changed = true; }
        ImGui::SameLine();
        if (ImGui::SmallButton("2 sur 3")) { balls_ = 5; launched_ = 2; changed = true; }
        ImGui::SameLine();
        if (ImGui::SmallButton("1 sur 2 (le paradoxe)")) { balls_ = 3; launched_ = 1; changed = true; }
    }

    ImGui::SeparatorText("Berceau");
    pushSliderWidth();
    {
        int n = launched_;
        if (ImGui::SliderInt("Billes lancées", &n, 1, balls_ - 1)) { launched_ = n; changed = true; }
        if (interesse) {
            int b = balls_;
            if (ImGui::SliderInt("Nombre de billes", &b, 3, kMaxBalls)) { balls_ = b; changed = true; }
        }
        if (interesse) changed |= sliderD(college ? "Restitution e" : "Rebond (mou ↔ vif)", &restitution_, 0.0, 1.0, college ? "%.2f" : "");
        if (lycee) {
            changed |= sliderD("Dureté k des billes", &stiffness_, 1e2, 1e6, "%.0e", ImGuiSliderFlags_Logarithmic);
            changed |= sliderD("Distance avant le choc", &gap_, 0.0, 3.0, "%.2f");
        }
        if (etudiant) {
            changed |= sliderD("Vitesse de lancement", &speed_, 0.2, 3.0, "%.2f", ImGuiSliderFlags_Logarithmic);
            changed |= sliderD("Amortissement α (Hertz)", &damping_, 0.0, 0.3, "%.3f");
        }
    }
    popSliderWidth();

    if (college) {
        ImGui::TextDisabled("Afficher");
        int d = static_cast<int>(display_);
        ImGui::RadioButton("Les deux (comparaison)", &d, kBoth);
        ImGui::RadioButton("Billes réelles (contact de Hertz)", &d, kHertz);
        ImGui::RadioButton("Chocs instantanés", &d, kImpulses);
        display_ = static_cast<Display>(d);
    }

    if (etudiant) {
        ImGui::SeparatorText("Calcul");
        const char* names = "Euler explicite\0Euler symplectique\0Verlet (vitesses)\0RK4\0RK45 adaptatif\0";
        pushSliderWidth();
        if (ImGui::Combo("Schéma numérique", &solverIndex_, names)) changed = true;
        if (chercheur) {
            changed |= sliderD("Pas par contact", &stepsPerContact_, 10.0, 400.0, "%.0f", ImGuiSliderFlags_Logarithmic);
            if (solverIndex_ == SolverSet::kRK45) changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-10, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
        }
        popSliderWidth();
    }

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
    if (college) ImGui::Text("t = %.3f / %.2f s", clock_.time, duration_);
    if (diverged_) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.40f, 1.0f));
        ImGui::TextWrapped("Le calcul a divergé : le pas est trop grand pour ce schéma (l'énergie explose). Augmentez les pas par contact ou changez de schéma.");
        ImGui::PopStyleColor();
    }
    if (interesse) wrapped(slowMotion() ? "Ralenti pendant le choc." : "Vitesse normale.", true);
    if (interesse) wrapped("Devant : billes réelles (Hertz), barres jaunes = force de chaque contact. Derrière : chocs instantanés.", true);
}

void CradleModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Dans un choc, la « quantité de mouvement » totale ne change pas, et l'énergie de mouvement non plus quand les billes "
                           "sont élastiques. Mais ces deux règles ne disent pas toujours quelles billes repartent : c'est la manière dont les "
                           "billes se touchent, très brièvement et en se déformant, qui décide.");
        return;
    }
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    const double t = clock_.time;
    const int n = balls_;
    const double mu = 0.5 * prob_.mass;

    ImGui::SeparatorText("Théorie");
    wrapped(strf("%d bille(s) lancée(s) sur %d : impulsion p = %.2f, énergie E = %.3f (masse 1, vitesse %.2f)", launched_, n,
                 launched_ * speed_, 0.5 * launched_ * speed_ * speed_, speed_));
    if (lycee) {
        const double dMax = hertz::maxCompression(mu, speed_, stiffness_), tc = hertz::contactDuration(mu, speed_, stiffness_);
        wrapped(strf("Deux billes (exact) : δ_max = (5 μ v²/(4 k))^(2/5) = %.5f ; T = 2,943 δ_max / v = %.5f s", dMax, tc));
        wrapped(strf("Mesuré sur la chaîne : δ_max = %.5f", maxCompression_) +
                (contactStart_ >= 0.0 && contactEnd_ >= 0.0 ? strf(" ; 1er contact de %.5f s", contactEnd_ - contactStart_) : std::string(" ; contact en cours ou à venir")));
    }
    if (etudiant && n == 3 && launched_ == 1) {
        const double xH = -outcome_.velocities[0] / speed_, xI = std::abs(impulses_.velocities[0]) < 5e-5 ? 0.0 : -impulses_.velocities[0];
        wrapped(strf("3 billes, courbe des solutions (impulsion et énergie conservées) : 0 ≤ x ≤ 1/3 avec x = −v1/v. Hertz : x = %.4f ; "
                     "chocs instantanés : x = %.4f%s.",
                     xH, xI, restitution_ < 1.0 ? " (hors de la courbe : de l'énergie est perdue)" : ""));
    }

    // Vitesses de chaque bille à l'instant courant
    ImGui::SeparatorText("Vitesses des billes (v / v0)");
    std::vector<TableRow> rows;
    for (int i = 0; i < n; ++i) {
        const double h = y_[n + i] / speed_, im = impulseVelocity(i, t) / speed_;
        rows.push_back({strf("Bille %d", i + 1), kBallColors[i], {strf("%+.3f", h), strf("%+.3f", im), strf("%+.3f", h - im)}});
    }
    drawResultTable("vitesses", {"Hertz", "Impuls.", "Écart"}, rows);
    wrapped("Écart = Hertz − impulsions. Hertz : les premières billes reculent un peu et celles du milieu bougent ; chocs instantanés : une bille ne recule jamais.", true);

    if (lycee) {
        ImGui::SeparatorText("Bilan");
        std::vector<TableRow> bilan;
        double hp = 0.0, ip = 0.0, ie = 0.0, rp = 0.0;
        for (int i = 0; i < n; ++i) {
            const double iv = impulseVelocity(i, t);
            ip += iv;
            ie += 0.5 * iv * iv;
            hp += y_[n + i];
            rp += yRef_[n + i];
        }
        const double p0 = launched_ * speed_;
        const double hEnergy = prob_.energy(y_) / energy0_;
        bilan.push_back({solvers_.shortLabel(level, solverIndex_) + " (Hertz)", solvers_.color(solverIndex_),
                         {strf("%.6f", hp / p0), strf("%.6f", hEnergy), strf("%+.3f", y_[2 * n - 1] / speed_)}});
        bilan.push_back({"Impulsions", kForceColor, {strf("%.6f", ip / p0), strf("%.6f", ie / energy0_), strf("%+.3f", impulseVelocity(n - 1, t) / speed_)}});
        if (etudiant) {
            bilan.push_back({"Référence (RK45)", kBlue, {strf("%.6f", rp / p0), strf("%.6f", prob_.energy(yRef_) / energy0_), strf("%+.3f", yRef_[2 * n - 1] / speed_)}});
        }
        drawResultTable("bilan", {"P / P0", "E / E0", "dernière"}, bilan);
        wrapped("E comprend l'énergie élastique stockée dans les contacts (elle revient à zéro quand les billes se séparent). dernière : vitesse de la dernière bille.", true);
        if (damping_ > 0.0) wrapped(strf("Amortissement α = %.3f : l'énergie décroît, restitution d'un choc de deux billes e ≈ 1/(1 + α v) = %.3f.", damping_, 1.0 / (1.0 + damping_ * speed_)), true);
    }

    if (chercheur) {
        ImGui::SeparatorText("Ordre de résolution (chocs instantanés)");
        const std::vector<double> squeezed = {1.0, 0.5, 0.0, 0.0};
        const cradle::ImpulseResult a = cradle::sequentialImpulses(squeezed, restitution_, cradle::ResolveOrder::LeftToRight);
        const cradle::ImpulseResult b = cradle::sequentialImpulses(squeezed, restitution_, cradle::ResolveOrder::RightToLeft);
        wrapped(strf("Bille prise entre deux voisines qui s'approchent déjà, vitesses (1 ; 0,5 ; 0 ; 0), e = %.2f :", restitution_));
        wrapped("  gauche → droite : " + velocityText(a.velocities, 1.0));
        wrapped("  droite → gauche : " + velocityText(b.velocities, 1.0));
        wrapped("Avec e < 1 les deux ordres diffèrent. Dans le berceau (groupe lancé contre billes au repos) ils donnent exactement le même résultat : "
                "seule la paire de tête se rapproche à chaque instant.", true);
        if (solverIndex_ == SolverSet::kRK45)
            wrapped(strf("RK45 : %d pas acceptés, %d refusés, %d évaluations.", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps(), solvers_.rk45().evaluations()), true);
    }
}

void CradleModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showH = display_ != kImpulses, showI = display_ != kHertz;
    const bool showForce = atLeast(level, Level::Lycee) && showH;
    const bool showEnergy = atLeast(level, Level::Etudiant) && showH;
    const int cols = (showH ? 1 : 0) + (showI ? 1 : 0) + (showForce ? 1 : 0) + (showEnergy ? 1 : 0);
    const int n = balls_;

    if (!ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) return;
    if (showH && ImPlot::BeginPlot("Vitesses : contact de Hertz")) {
        ImPlot::SetupAxes("t (s)", "v / v0", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        for (int i = 0; i < n; ++i)
            if (vel_[i].size() > 0) ImPlot::PlotLine(strf("Bille %d", i + 1).c_str(), vel_[i].x.data(), vel_[i].y.data(), vel_[i].size(), lineSpec(kBallColors[i], vel_[i].offset));
        ImPlot::EndPlot();
    }
    if (showI && ImPlot::BeginPlot("Vitesses : chocs instantanés")) {
        ImPlot::SetupAxes("t (s)", "v / v0", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        for (int i = 0; i < n; ++i)
            if (impVel_[i].size() > 0) ImPlot::PlotLine(strf("Bille %d", i + 1).c_str(), impVel_[i].x.data(), impVel_[i].y.data(), impVel_[i].size(), lineSpec(kBallColors[i], impVel_[i].offset));
        ImPlot::EndPlot();
    }
    if (showForce && ImPlot::BeginPlot("Force des contacts (onde de compression)")) {
        ImPlot::SetupAxes("t (s)", "F / F_max(2 billes)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        for (int i = 0; i + 1 < n; ++i)
            if (force_[i].size() > 0) ImPlot::PlotLine(strf("Contact %d-%d", i + 1, i + 2).c_str(), force_[i].x.data(), force_[i].y.data(), force_[i].size(), lineSpec(kBallColors[i], force_[i].offset));
        ImPlot::EndPlot();
    }
    if (showEnergy && ImPlot::BeginPlot("Énergie (Hertz)")) {
        ImPlot::SetupAxes("t (s)", "E / E0", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        if (kinetic_.size() > 0) {
            ImPlot::PlotLine("Cinétique", kinetic_.x.data(), kinetic_.y.data(), kinetic_.size(), lineSpec(kKineticColor, kinetic_.offset));
            ImPlot::PlotLine("Élastique (contacts)", potential_.x.data(), potential_.y.data(), potential_.size(), lineSpec(kForceColor, potential_.offset));
            ImPlot::PlotLine("Totale", total_.x.data(), total_.y.data(), total_.size(), lineSpec(kTotalColor, total_.offset));
        }
        ImPlot::EndPlot();
    }
    ImPlot::EndSubplots();
}

void CradleModule::drawAnalysis(const UiContext& ctx) {
    const bool chercheur = atLeast(ctx.level, Level::Chercheur);
    if (convergenceDirty_) computeConvergence();
    const bool family = balls_ == 3 && launched_ == 1;

    auto markerSpec = [](const float* c, ImPlotMarker m = ImPlotMarker_Circle, float size = 4.0f) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, m, ImPlotProp_MarkerSize, size,
                          ImPlotProp_MarkerFillColor, toImVec4(c), ImPlotProp_MarkerLineColor, toImVec4(c));
    };
    const int cols = 1 + (family ? 1 : 0) + (chercheur ? 1 : 0);

    if (!ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) return;

    if (family) {
        // L'ensemble des sorties possibles de 3 billes : le résultat de chaque modèle y est placé.
        std::vector<double> xs, v2, v3;
        for (int k = 0; k <= 100; ++k) {
            const double x = k / 300.0;
            const std::array<double, 3> v = cradle::threeBallFamily(x, 1.0);
            xs.push_back(x);
            v2.push_back(v[1]);
            v3.push_back(v[2]);
        }
        const double xH = -outcome_.velocities[0] / speed_, xI = -impulses_.velocities[0];
        const double h2 = outcome_.velocities[1] / speed_, h3 = outcome_.velocities[2] / speed_;
        const double i2 = impulses_.velocities[1], i3 = impulses_.velocities[2];
        if (ImPlot::BeginPlot("Sorties possibles de 3 billes")) {
            ImPlot::SetupAxes("recul de la 1re bille x = −v1/v0", "v / v0");
            ImPlot::SetupAxesLimits(-0.06, 0.36, -0.1, 1.1);
            ImPlot::SetupLegend(ImPlotLocation_West);
            ImPlot::PlotLine("2e bille", xs.data(), v2.data(), static_cast<int>(xs.size()), lineSpec(kBallColors[1], 0));
            ImPlot::PlotLine("3e bille", xs.data(), v3.data(), static_cast<int>(xs.size()), lineSpec(kBallColors[2], 0));
            ImPlot::PlotLine("Hertz", &xH, &h2, 1, markerSpec(kTotalColor, ImPlotMarker_Circle, 8.0f));
            ImPlot::PlotLine("##hertz3", &xH, &h3, 1, markerSpec(kTotalColor, ImPlotMarker_Circle, 8.0f));
            ImPlot::PlotLine("Chocs instantanés", &xI, &i2, 1, markerSpec(kForceColor, ImPlotMarker_Square, 8.0f));
            ImPlot::PlotLine("##imp3", &xI, &i3, 1, markerSpec(kForceColor, ImPlotMarker_Square, 8.0f));
            ImPlot::EndPlot();
        }
    }

    if (ImPlot::BeginPlot(strf("Convergence : erreur à t = %.2f s", convergenceTime_).c_str())) {
        ImPlot::SetupAxes("dt (s)", "écart (phase)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        for (int i = 0; i < SolverSet::kFixedStep; ++i) {
            const Curve& c = convergence_[i];
            if (!solvers_.show(i) || c.x.size() < 2) continue;
            const std::string name = std::isfinite(c.slope) ? strf("%s (pente %.2f)", solvers_.solver(i).name(), c.slope) : strf("%s", solvers_.solver(i).name());
            ImPlot::PlotLine(name.c_str(), c.x.data(), c.y.data(), static_cast<int>(c.x.size()), markerSpec(solvers_.color(i)));
        }
        ImPlot::EndPlot();
    }

    if (chercheur && ImPlot::BeginPlot("Écart au calcul de référence (RK45 serré)")) {
        ImPlot::SetupAxes("t (s)", "écart (phase)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        if (error_.size() > 0) ImPlot::PlotLine(solvers_.solver(solverIndex_).name(), error_.x.data(), error_.y.data(), error_.size(), lineSpec(solvers_.color(solverIndex_), error_.offset));
        ImPlot::EndPlot();
    }
    ImPlot::EndSubplots();
}

// ------------------------------ rendu 3D --------------------------------

void CradleModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    const bool showH = display_ != kImpulses, showI = display_ != kHertz;
    const bool arrows = atLeast(ctx.level, Level::Lycee);
    const double R = prob_.radius, t = clock_.time;
    const int n = balls_;
    const float x0 = static_cast<float>(viewMin_), x1 = static_cast<float>(viewMax_);

    auto drawArrow = [&](double x, double z, double v) {
        const float fx = static_cast<float>(x), fz = static_cast<float>(z), len = static_cast<float>(v / speed_);
        const float y = 0.6f;
        renderer.draw(Primitive::Lines, {{fx, y, fz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]},
                                         {fx + len, y, fz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}});
        if (std::abs(len) > 0.08f) {   // pointe
            const float s = len > 0 ? -1.0f : 1.0f;
            renderer.draw(Primitive::Lines, {{fx + len, y, fz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]},
                                             {fx + len + 0.12f * s, y, fz + 0.08f, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]},
                                             {fx + len, y, fz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]},
                                             {fx + len + 0.12f * s, y, fz - 0.08f, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}});
        }
    };

    if (showH && !diverged_) {
        const double z = showI ? kLaneZ : 0.0;
        renderer.draw(Primitive::Lines, {{x0, 0.0f, static_cast<float>(z), kRail[0], kRail[1], kRail[2]}, {x1, 0.0f, static_cast<float>(z), kRail[0], kRail[1], kRail[2]}});
        std::vector<Vertex> centers;
        for (int i = 0; i < n; ++i) {
            const float* c = kBallColors[i];
            renderer.draw(Primitive::LineStrip, circleXZ(y_[i], z, R, c));
            renderer.draw(Primitive::LineStrip, circleXZ(y_[i], z, R * 0.93, c, 0.6f));
            centers.push_back({static_cast<float>(y_[i]), 0.0f, static_cast<float>(z), c[0], c[1], c[2]});
            if (arrows) drawArrow(y_[i], z, y_[n + i]);
        }
        renderer.draw(Primitive::Points, centers, 7.0f * ctx.uiScale);
        // force de chaque contact : barre verticale au point de contact (hauteur proportionnelle à F)
        for (int i = 0; i + 1 < n; ++i) {
            const double f = contactForce(y_, i) / forceScale_;
            if (f <= 1e-9) continue;
            const float px = static_cast<float>(0.5 * (y_[i] + y_[i + 1])), h = static_cast<float>(std::min(2.5 * f, 4.0));
            renderer.draw(Primitive::Lines, {{px, 0.0f, static_cast<float>(z), kForceColor[0], kForceColor[1], kForceColor[2]},
                                             {px, h, static_cast<float>(z), kForceColor[0], kForceColor[1], kForceColor[2]}});
            renderer.draw(Primitive::Points, {{px, h, static_cast<float>(z), kForceColor[0], kForceColor[1], kForceColor[2]}}, 6.0f * ctx.uiScale);
        }
    }

    if (showI) {
        const double z = showH ? -kLaneZ : 0.0;
        renderer.draw(Primitive::Lines, {{x0, 0.0f, static_cast<float>(z), kRail[0], kRail[1], kRail[2]}, {x1, 0.0f, static_cast<float>(z), kRail[0], kRail[1], kRail[2]}});
        std::vector<Vertex> centers;
        for (int i = 0; i < n; ++i) {
            const double x = impulsePosition(i, t);
            const float* c = kBallColors[i];
            renderer.draw(Primitive::LineStrip, circleXZ(x, z, R, c, 0.7f));
            centers.push_back({static_cast<float>(x), 0.0f, static_cast<float>(z), c[0] * 0.7f, c[1] * 0.7f, c[2] * 0.7f});
            if (arrows) drawArrow(x, z, impulseVelocity(i, t));
        }
        renderer.draw(Primitive::Points, centers, 6.0f * ctx.uiScale);
    }
}

}  // namespace pl
