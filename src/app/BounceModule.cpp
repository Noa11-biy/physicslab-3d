#include "BounceModule.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <implot.h>

namespace pl {
namespace {

constexpr double kPi = constants::pi;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kDeg = 180.0 / kPi;
constexpr std::size_t kTrailLength = 2500;
constexpr double kBallScaleBalls = 1.4;      // échelle de la vue de dessus des deux billes (unités de scène par mètre)
const float kGround[3] = {0.55f, 0.55f, 0.6f};
const float kBall2Color[3] = {0.60f, 0.80f, 1.00f};
const float kVelocityColor[3] = {0.35f, 1.0f, 0.45f};

std::vector<Vertex> circle(double cx, double cy, double r, double scale, const float* c, float dim = 1.0f) {
    std::vector<Vertex> v;
    for (int i = 0; i <= 40; ++i) {
        const double a = 2.0 * kPi * i / 40.0;
        v.push_back({static_cast<float>((cx + r * std::cos(a)) * scale), 0.0f, -static_cast<float>((cy + r * std::sin(a)) * scale),
                     c[0] * dim, c[1] * dim, c[2] * dim});
    }
    return v;
}

}  // namespace

BounceModule::BounceModule() {
    reset();
}

const char* BounceModule::explanation(Level level) const {
    const bool ball = scenario_ == kBall;
    switch (level) {
        case Level::Vulgarisation:
            return ball ? "Une balle lâchée rebondit de moins en moins haut : à chaque rebond elle perd un peu de son élan, qui se "
                          "transforme en chaleur et en bruit. Elle finit par s'arrêter, et pourtant elle a rebondi une infinité de "
                          "fois, de plus en plus vite : tout cela se passe en un temps fini !\n\n"
                          "La balle bleue suit les vraies lois, la verte est le calcul de l'ordinateur. Regardez-la au repos : "
                          "l'ordinateur la laisse trembler. Cochez « Calcul amélioré »."
                        : "Deux billes se heurtent. Si elles sont identiques et que l'une est immobile, la première s'arrête et la "
                          "seconde repart avec toute sa vitesse (choc de face) ; si le choc est de côté, elles repartent à angle droit.\n\n"
                          "Les billes bleues suivent les vraies lois, les vertes sont le calcul de l'ordinateur : il rate l'instant "
                          "du contact. Cochez « Calcul amélioré ».";
        case Level::Interesse:
            return ball ? "Le coefficient de restitution e dit quelle part de la vitesse la balle garde à chaque rebond : e = 1 pour un "
                          "rebond parfait, e = 0 pour une balle de pâte à modeler qui ne rebondit pas. La hauteur est multipliée "
                          "par e² à chaque rebond.\n\n"
                          "Tous les rebonds finissent par se serrer : la balle s'arrête après un temps fini, même s'il y en a une "
                          "infinité. C'est le paradoxe de Zénon, version balle."
                        : "Dans un choc, la « quantité de mouvement » totale (masse × vitesse) ne change jamais. L'énergie, elle, "
                          "se perd en partie (chaleur, bruit), sauf si le choc est parfaitement élastique (e = 1).\n\n"
                          "Essayez un choc de côté : avec e = 1 et deux billes identiques, les billes repartent exactement à angle "
                          "droit ; avec e plus petit, l'angle se referme.";
        case Level::College:
            return ball ? "Après chaque rebond la vitesse est multipliée par e : v' = e × v, donc la hauteur par e². Hauteur du n-ième "
                          "rebond : h_n = e^(2n) h0. Temps total avant l'arrêt : t_total = t0 (1 + e) / (1 − e), avec "
                          "t0 = √(2 h0 / g) le temps de la première chute (un nombre fini de secondes malgré une infinité de rebonds).\n"
                          "Le calcul naïf voit le sol trop tard (après le pas) et ne s'arrête jamais exactement."
                        : "Impulsion (quantité de mouvement) p = m v : conservée dans un choc. Restitution : la vitesse relative "
                          "d'éloignement vaut e fois celle d'approche. Énergie perdue : ½ (1 − e²) μ v², avec μ = m1 m2 / (m1 + m2).\n"
                          "Dans « Invariants » : le calcul naïf détecte le choc en retard et décale tout ce qui suit.";
        case Level::Lycee:
            return ball ? "Chute : y'' = −g, v = g t, h = ½ g t². Rebond : v' = −e v (vers le haut). Les vols successifs durent "
                          "2 e^n t0 (suite géométrique de raison e) : leur somme vaut t0 (1 + e)/(1 − e). Énergie : E = ½ v² + g y est "
                          "multipliée par e² à chaque rebond ; la différence est de la chaleur.\n"
                          "Dans « Analyse » : les hauteurs et les durées en échelle logarithmique sont des droites (suites géométriques)."
                        : "Au contact, la force est dirigée selon la normale n (ligne des centres), les composantes tangentielles ne "
                          "changent pas. Impulsion J = (1 + e) μ v_n le long de n ; v1' = v1 − (J/m1) n, v2' = v2 + (J/m2) n.\n"
                          "Billes égales, l'une au repos, e = 1 : angle de 90° ; si e < 1 : v1'·v2' = (1 − e²)/4 v_n², angle aigu "
                          "(voir « Analyse » : perte d'énergie et angle selon le décalage).";
        case Level::Etudiant:
            return ball ? "Système hybride : vol (EDO lisse) + saut à la surface y = 0 (reset vy ← −e vy). L'instant du saut se "
                          "détermine par une fonction d'événement g = y ; l'erreur d'un pas fixe est d'ordre 1 car le contact tombe "
                          "au milieu du pas. Avec résistance k > 0 les vols ne sont plus des polynômes : les ordres des solveurs "
                          "réapparaissent avec l'événement (graphe de convergence).\n"
                          "Accumulation de Zénon : les durées de vol tendent vers 0 ; on arrête la balle sous un seuil de vitesse."
                        : "Restitution de Newton : (v1' − v2')·n = −e (v1 − v2)·n. Hors choc le mouvement est rectiligne : TOUS les "
                          "solveurs sont exacts ; toute l'erreur vient de l'instant du contact, d'où le même écart pour les quatre schémas "
                          "(ordre 1, RK4 compris) avec le modèle naïf, nul avec l'événement.";
        case Level::Chercheur:
            return ball ? "Traitement événementiel (event-driven) : on localise l'impact par bissection sur le pas, on applique le "
                          "saut, puis on repart de la surface par un pas de décollage (g = 0 au départ : sans lui, un vol plus court "
                          "que le pas échappe à la détection et la balle traverse le sol). L'accumulation de Zénon n'est pas "
                          "dépassable en temps événementiel : seuil de vitesse, ou schéma time-stepping (Moreau-Jean) qui traite le "
                          "contact comme une complémentarité sur chaque pas. Restitution de Newton vs Poisson vs énergétique (Stronge) : "
                          "équivalentes sans frottement ; avec frottement elles divergent (paradoxe de Painlevé). Contact continu de "
                          "Hertz : M5c."
                        : "Même structure hybride. Le contact est détecté par la fonction d'événement max(distance − R, 0) (nulle après "
                          "le choc, donc pas de faux événement au redémarrage). Le choc est instantané : l'impulsion J est une masse "
                          "de Dirac de la force. Billes avec frottement tangentiel : cône de Coulomb dans le plan du contact, "
                          "restitution tangentielle (hors périmètre ici).";
    }
    return "";
}

void BounceModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 1.5f;
    camera.target[2] = -1.0f;
    camera.distance = 24.0f;
    camera.yaw = 0.25f;
    camera.pitch = 0.45f;
}

// ------------------------------ simulation -----------------------------

void BounceModule::reset() {
    ball_ = BounceProblem{};
    ball_.restitution = restitution_;
    ball_.y0 = height_;
    ball_.vx0 = horizontal_;
    ball_.drag = drag_;
    ball_.restSpeed = restSpeed_;

    balls_ = TwoBallProblem{};
    balls_.restitution = restitutionBalls_;
    balls_.p1 = {-4.0, impact_, 0.0};
    balls_.v1 = {approach_, 0.0, 0.0};
    balls_.m2 = massRatio_ * balls_.m1;

    impactX_.clear();
    impactT_.clear();
    if (scenario_ == kBall) {
        const double tRest = ball_.restTime();
        duration_ = std::isfinite(tRest) ? std::clamp(1.15 * tRest + 1.0, 4.0, 30.0) : 20.0;
        convergenceTime_ = std::clamp(0.37 * (std::isfinite(tRest) ? tRest : 10.0), 1.0, 12.0);
        impactT_ = ball_.impactTimes(40);
        for (double t : impactT_) impactX_.push_back(ball_.exact(t).x);
        const double xEnd = ball_.exact(duration_).x;
        const double range = std::max(xEnd - ball_.x0, 4.0);
        scale_ = std::clamp(18.0 / range, 0.35, 3.0);
        offsetX_ = -0.5 * (ball_.x0 + xEnd) * scale_;
    } else {
        const double tc = balls_.collisionTime();
        duration_ = std::isfinite(tc) ? tc + 2.5 : 5.0;
        convergenceTime_ = duration_;
        scale_ = kBallScaleBalls;
        offsetX_ = 0.0;
    }

    for (Run& r : runs_) {
        r.ball = scenario_ == kBall ? std::make_unique<BounceRun>(ball_, model_) : nullptr;
        r.balls = scenario_ == kBalls ? std::make_unique<TwoBallRun>(balls_, model_) : nullptr;
        r.primary.clear();
        r.secondary.clear();
        r.error.clear();
        r.energy.clear();
        r.trail.clear();
        r.trail2.clear();
    }
    exactPrimary_.clear();
    exactSecondary_.clear();
    exactEnergy_.clear();
    solvers_.rk45().resetStats();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    sample();
}

void BounceModule::sample() {
    const double t = clock_.time;
    if (scenario_ == kBall) {
        const BounceState ex = ball_.exact(t);
        exactPrimary_.add(t, ex.y);
        exactSecondary_.add(t, ex.vy);
        exactEnergy_.add(t, ball_.energy(ex.y, ex.vx, ex.vy));
        for (int i = 0; i < SolverSet::kCount; ++i) {
            Run& r = runs_[i];
            const State& y = r.ball->state();
            r.primary.add(t, y[1]);
            r.secondary.add(t, y[3]);
            r.error.add(t, std::max(std::hypot(y[0] - ex.x, y[1] - ex.y), 1e-12));
            r.energy.add(t, ball_.energy(y[1], y[2], y[3]));
            r.trail.push_back({static_cast<float>(y[0] * scale_ + offsetX_), static_cast<float>(y[1] * scale_), 0.0f, solvers_.color(i)[0] * 0.7f,
                               solvers_.color(i)[1] * 0.7f, solvers_.color(i)[2] * 0.7f});
            if (r.trail.size() > kTrailLength) r.trail.erase(r.trail.begin(), r.trail.begin() + 300);
        }
    } else {
        const State ex = balls_.exact(t);
        exactPrimary_.add(t, balls_.gap(ex));
        exactSecondary_.add(t, std::hypot(ex[4], ex[5]));
        exactEnergy_.add(t, balls_.kineticEnergy(ex));
        for (int i = 0; i < SolverSet::kCount; ++i) {
            Run& r = runs_[i];
            const State& y = r.balls->state();
            r.primary.add(t, balls_.gap(y));
            r.secondary.add(t, std::hypot(y[4], y[5]));
            double err = 0.0;
            for (int k = 0; k < 4; ++k) err = std::max(err, std::abs(y[k] - ex[k]));
            r.error.add(t, std::max(err, 1e-12));
            r.energy.add(t, balls_.kineticEnergy(y));
            const float* c = solvers_.color(i);
            r.trail.push_back({static_cast<float>(y[0] * scale_), 0.0f, -static_cast<float>(y[1] * scale_), c[0] * 0.7f, c[1] * 0.7f, c[2] * 0.7f});
            r.trail2.push_back({static_cast<float>(y[2] * scale_), 0.0f, -static_cast<float>(y[3] * scale_), c[0] * 0.7f, c[1] * 0.7f, c[2] * 0.7f});
            if (r.trail.size() > kTrailLength) { r.trail.erase(r.trail.begin(), r.trail.begin() + 300); r.trail2.erase(r.trail2.begin(), r.trail2.begin() + 300); }
        }
    }
}

void BounceModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;

    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    clock_.advance(
        frameSeconds, timeScale_, dt_, duration_,
        [&](double h) {
            for (int i = 0; i < SolverSet::kCount; ++i) {
                Run& r = runs_[i];
                if (r.ball) r.ball->advance(solvers_.solver(i), h);
                else r.balls->advance(solvers_.solver(i), h);
            }
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur de chaque solveur à pas fixe à `convergenceTime_`, pour différents pas, avec le modèle de contact courant.
void BounceModule::computeConvergence() {
    static const int kStepCounts[] = {50, 100, 200, 400, 800, 1600, 3200};
    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            const double e = scenario_ == kBall ? bounceError(ball_, model_, *solver, n, convergenceTime_)
                                                : twoBallError(balls_, model_, *solver, n, convergenceTime_);
            if (!std::isfinite(e) || e <= 1e-12) continue;  // sous ce seuil : arrondi
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

// Ce que la théorie prévoit, en mots.
std::string BounceModule::verdict() const {
    if (scenario_ == kBall) {
        const double tRest = ball_.restTime();
        if (!std::isfinite(tRest)) return "Rebond parfait : la balle rebondit indéfiniment à la même hauteur.";
        return strf("La balle s'arrête de rebondir après %d rebonds, à t = %.2f s : un temps FINI malgré ces rebonds de plus en plus rapprochés.",
                    ball_.exact(tRest + 1.0).bounces, tRest);
    }
    const double tc = balls_.collisionTime();
    if (!std::isfinite(tc)) return "Les billes se manquent : pas de choc.";
    const double angle = balls_.outgoingAngle();
    if (std::isfinite(angle)) return strf("Choc à t = %.3f s ; les billes repartent à %.1f° l'une de l'autre.", tc, angle * kDeg);
    return strf("Choc de face à t = %.3f s : une bille repart seule ou toutes deux dans la même direction.", tc);
}

// --------------------------------- UI ----------------------------------

void BounceModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool etudiant = atLeast(level, Level::Etudiant);
    bool changed = false;

    ImGui::SeparatorText("Scénario");
    int s = scenario_;
    ImGui::RadioButton("Balle qui rebondit", &s, kBall);
    ImGui::RadioButton("Choc de deux billes", &s, kBalls);
    if (s != scenario_) {
        scenario_ = static_cast<Scenario>(s);
        reset();
        return;
    }

    if (interesse) {
        ImGui::TextDisabled("Cas types");
        if (scenario_ == kBall) {
            if (ImGui::SmallButton("Caoutchouc")) { restitution_ = 0.8; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Balle molle")) { restitution_ = 0.3; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Parfaite")) { restitution_ = 1.0; changed = true; }
        } else {
            if (ImGui::SmallButton("De face")) { impact_ = 0.0; restitutionBalls_ = 1.0; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("De côté")) { impact_ = 0.6; restitutionBalls_ = 1.0; changed = true; }
            ImGui::SameLine();
            if (ImGui::SmallButton("Choc mou")) { impact_ = 0.0; restitutionBalls_ = 0.0; changed = true; }
        }
    }

    ImGui::SeparatorText("Paramètres");
    pushSliderWidth();
    if (scenario_ == kBall) {
        changed |= sliderD(interesse ? "Restitution e" : "Rebond (mou ↔ vif)", &restitution_, 0.0, 1.0, interesse ? "%.2f" : "");
        changed |= sliderD(interesse ? "Hauteur de chute (m)" : "Hauteur de chute", &height_, 0.5, 5.0, interesse ? "%.1f" : "");
        if (interesse) changed |= sliderD("Vitesse horizontale (m/s)", &horizontal_, 0.0, 4.0, "%.1f");
        if (etudiant) changed |= sliderD("Résistance k (1/s)", &drag_, 0.0, 2.0, "%.2f");
        if (atLeast(level, Level::Chercheur)) changed |= sliderD("Seuil d'arrêt (m/s)", &restSpeed_, 1e-6, 1e-2, "%.0e", ImGuiSliderFlags_Logarithmic);
    } else {
        changed |= sliderD(interesse ? "Restitution e" : "Choc (mou ↔ vif)", &restitutionBalls_, 0.0, 1.0, interesse ? "%.2f" : "");
        changed |= sliderD(interesse ? "Décalage (m)" : "Décalage", &impact_, -1.2, 1.2, interesse ? "%.2f" : "");
        if (college) {
            changed |= sliderD("Masse bille 2 (x)", &massRatio_, 0.25, 4.0, "%.2f", ImGuiSliderFlags_Logarithmic);
            changed |= sliderD("Vitesse bille 1 (m/s)", &approach_, 1.0, 6.0, "%.1f");
        }
    }
    popSliderWidth();

    ImGui::SeparatorText("Calcul");
    if (!college) {
        bool improved = model_ == ContactModel::EventDriven;
        if (ImGui::Checkbox("Calcul amélioré (détecte le contact)", &improved)) {
            model_ = improved ? ContactModel::EventDriven : ContactModel::Naive;
            changed = true;
        }
    } else {
        int m = static_cast<int>(model_);
        bool modelChanged = ImGui::RadioButton("Naïf : contact vu après le pas", &m, static_cast<int>(ContactModel::Naive));
        modelChanged |= ImGui::RadioButton("Événement : contact à l'instant exact", &m, static_cast<int>(ContactModel::EventDriven));
        if (modelChanged) {
            model_ = static_cast<ContactModel>(m);
            changed = true;
        }
        pushSliderWidth();
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 2e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Chercheur)) changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-10, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
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

    if (college) ImGui::Text("t = %.2f / %.1f s", clock_.time, duration_);
    if (interesse) wrapped("Théorie : " + verdict());
}

void BounceModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Dans un choc, la « quantité de mouvement » totale ne change pas, mais une partie de l'énergie du "
                           "mouvement se transforme en chaleur et en bruit : voilà pourquoi la balle rebondit de moins en moins haut.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const double t = clock_.time;
    const double g = constants::g0;

    ImGui::SeparatorText("Théorie");
    std::vector<std::string> headersA, headersB;
    std::vector<TableRow> rowsA, rowsB;

    if (scenario_ == kBall) {
        const double e = ball_.restitution, t0 = std::sqrt(2.0 * ball_.y0 / g), tRest = ball_.restTime();
        wrapped(strf("e = %.2f ; première chute t0 = %.3f s ; vitesse au 1er impact %.2f m/s", e, t0, g * t0));
        if (ball_.drag == 0.0 && e < 1.0)
            wrapped(strf("t_total = t0 (1 + e)/(1 − e) = %.3f s ; h_n = e^(2n) h0 : h_1 = %.3f m, h_2 = %.3f m", t0 * (1 + e) / (1 - e), e * e * ball_.y0, std::pow(e, 4) * ball_.y0));
        else if (std::isfinite(tRest))
            wrapped(strf("Repos à t = %.3f s (avec résistance : calculé numériquement, sans formule fermée)", tRest));
        wrapped("→ " + verdict());

        const BounceState ex = ball_.exact(t);
        const double exE = ball_.energy(ex.y, ex.vx, ex.vy);
        headersA = {"y (m)", "vy (m/s)", "écart (m)"};
        headersB = {"rebonds", "t_repos", "E−E_ex"};
        rowsA.push_back({"Exacte", kBlue, {strf("%.3f", ex.y), strf("%+.3f", ex.vy), "-"}});
        rowsB.push_back({"Exacte", kBlue, {strf("%d", ex.bounces), std::isfinite(tRest) ? strf("%.3f", tRest) : "-", "-"}});
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i)) continue;
            const BounceRun& run = *runs_[i].ball;
            const State& y = run.state();
            TableRow a{solvers_.shortLabel(level, i), solvers_.color(i), {strf("%.3f", y[1]), strf("%+.3f", y[3]), strf("%.1e", std::hypot(y[0] - ex.x, y[1] - ex.y))}};
            TableRow b{a.name, a.color, {strf("%d", run.bounces()), run.resting() ? strf("%.3f", run.restTime()) : "-", strf("%+.1e", ball_.energy(y[1], y[2], y[3]) - exE)}};
            rowsA.push_back(a);
            rowsB.push_back(b);
        }
    } else {
        const double mu = balls_.m1 * balls_.m2 / (balls_.m1 + balls_.m2), tc = balls_.collisionTime();
        wrapped(strf("m1 = %.2f, m2 = %.2f ; μ = m1 m2/(m1+m2) = %.3f ; e = %.2f", balls_.m1, balls_.m2, mu, balls_.restitution));
        if (std::isfinite(tc)) {
            const State at = balls_.exact(tc);
            const Vec3 n = Vec3{at[2] - at[0], at[3] - at[1], 0.0}.normalized();
            const double vn = dot(balls_.v1 - balls_.v2, n);
            wrapped(strf("Vitesse d'approche normale v_n = %.3f m/s ; impulsion J = (1+e) μ v_n = %.3f kg·m/s", vn, (1.0 + balls_.restitution) * mu * vn));
            wrapped(strf("Énergie perdue = ½(1−e²) μ v_n² = %.3f J", 0.5 * (1.0 - balls_.restitution * balls_.restitution) * mu * vn * vn));
        }
        wrapped("→ " + verdict());

        const State ex = balls_.exact(t);
        const double e0 = balls_.kineticEnergy(balls_.initialState());
        const Vec3 p0 = balls_.momentum(balls_.initialState());
        headersA = {"écart (m)", "t_choc", "angle (°)"};
        headersB = {"|ΔP|", "ΔE/E0", "chocs"};
        const double angle = balls_.outgoingAngle();
        rowsA.push_back({"Exacte", kBlue, {"-", std::isfinite(tc) ? strf("%.3f", tc) : "-", std::isfinite(angle) ? strf("%.1f", angle * kDeg) : "-"}});
        rowsB.push_back({"Exacte", kBlue, {strf("%.0e", (balls_.momentum(ex) - p0).norm()), strf("%+.2f", (balls_.kineticEnergy(ex) - e0) / e0), std::isfinite(tc) && t > tc ? "1" : "0"}});
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i)) continue;
            const TwoBallRun& run = *runs_[i].balls;
            const State& y = run.state();
            double err = 0.0;
            for (int k = 0; k < 4; ++k) err = std::max(err, std::abs(y[k] - ex[k]));
            const double n1 = std::hypot(y[4], y[5]), n2 = std::hypot(y[6], y[7]);
            const bool after = run.collisions() > 0 && n1 > 1e-9 && n2 > 1e-9;
            TableRow a{solvers_.shortLabel(level, i), solvers_.color(i),
                       {strf("%.1e", err), std::isfinite(run.collisionTime()) ? strf("%.3f", run.collisionTime()) : "-",
                        after ? strf("%.1f", std::acos(std::clamp((y[4] * y[6] + y[5] * y[7]) / (n1 * n2), -1.0, 1.0)) * kDeg) : "-"}};
            TableRow b{a.name, a.color,
                       {strf("%.0e", (balls_.momentum(y) - p0).norm()), strf("%+.2f", (balls_.kineticEnergy(y) - e0) / e0), strf("%d", run.collisions())}};
            rowsA.push_back(a);
            rowsB.push_back(b);
        }
    }

    ImGui::SeparatorText("Comparaison à l'instant courant");
    drawResultTable("comparaison", headersA, rowsA);
    if (lycee) {
        ImGui::SeparatorText(scenario_ == kBall ? "Rebonds et énergie" : "Impulsion et énergie");
        drawResultTable("bilan", headersB, rowsB);
    }
    if (scenario_ == kBall) {
        wrapped("écart : distance à la balle exacte. « - » : la balle naïve ne s'arrête jamais (t_repos non défini).", true);
        if (etudiant) wrapped("E−E_ex : énergie par unité de masse (½v² + g y) par rapport à l'exacte.", true);
    } else {
        wrapped("écart : plus grande erreur de position des deux billes. angle : entre les vitesses finales (à 90° si e = 1 et masses égales, bille au repos).", true);
        if (etudiant) wrapped("L'impulsion est conservée même par le calcul naïf (impulsions égales et opposées) ; l'erreur vient de l'instant du contact.", true);
    }
}

void BounceModule::drawGraphs(const UiContext& ctx) {
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
    const bool ball = scenario_ == kBall;

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot(ball ? "Hauteur y(t)" : "Distance entre les surfaces")) {
            ImPlot::SetupAxes("t (s)", ball ? "y (m)" : "écart (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactPrimary_);
            plotRuns(&Run::primary);
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot(ball ? "Vitesse verticale vy(t)" : "Vitesse de la bille 1")) {
            ImPlot::SetupAxes("t (s)", "m/s", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactSecondary_);
            plotRuns(&Run::secondary);
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Écart de position")) {
            ImPlot::SetupAxes("t (s)", "m", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plotRuns(&Run::error);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot(ball ? "Énergie ½v² + g y" : "Énergie cinétique")) {
            ImPlot::SetupAxes("t (s)", ball ? "J/kg" : "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotExact(exactEnergy_);
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void BounceModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool researcher = atLeast(level, Level::Chercheur);
    if (etudiant && convergenceDirty_) computeConvergence();
    const bool ball = scenario_ == kBall;
    const int cols = 2 + (etudiant ? 1 : 0) + (researcher && ball ? 1 : 0);

    auto markerSpec = [](const float* c) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle,
                          ImPlotProp_MarkerSize, 4.0f);
    };

    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ball) {
            // Hauteurs (sommets) et durées des vols : suites géométriques de raison e^2 et e : des droites en échelle log.
            std::vector<double> n, heights, durations;
            for (std::size_t k = 0; k + 1 < impactT_.size(); ++k) {
                const BounceState apex = ball_.exact(0.5 * (impactT_[k] + impactT_[k + 1]));
                n.push_back(static_cast<double>(k + 1));
                heights.push_back(std::max(apex.y, 1e-12));
                durations.push_back(impactT_[k + 1] - impactT_[k]);
            }
            if (ImPlot::BeginPlot("Hauteur des rebonds h_n")) {
                ImPlot::SetupAxes("rebond n", "h (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                if (!n.empty()) ImPlot::PlotLine("sommet", n.data(), heights.data(), static_cast<int>(n.size()), markerSpec(kBlue));
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Durée des vols (Zénon)")) {
                ImPlot::SetupAxes("rebond n", "durée (s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                if (!n.empty()) ImPlot::PlotLine("vol n", n.data(), durations.data(), static_cast<int>(n.size()), markerSpec(kBall2Color));
                ImPlot::EndPlot();
            }
        } else {
            // Perte d'énergie et angle de sortie selon le décalage (exact, mêmes masses, e courant).
            std::vector<double> bs, loss, bsAngle, angle;
            for (int k = 0; k <= 120; ++k) {
                TwoBallProblem q = balls_;
                const double R = q.r1 + q.r2, b = -0.98 * R + 1.96 * R * k / 120.0;
                q.p1.y = b;
                const double tc = q.collisionTime();
                if (!std::isfinite(tc)) continue;
                const State a = q.exact(tc);
                const Vec3 nn = Vec3{a[2] - a[0], a[3] - a[1], 0.0}.normalized();
                const double mu = q.m1 * q.m2 / (q.m1 + q.m2), vn = dot(q.v1 - q.v2, nn);
                bs.push_back(b);
                loss.push_back(0.5 * (1.0 - q.restitution * q.restitution) * mu * vn * vn);
                const double ang = q.outgoingAngle();
                if (std::isfinite(ang)) {  // angle indéfini quand une bille s'arrête (choc de face) : point omis
                    bsAngle.push_back(b);
                    angle.push_back(ang * kDeg);
                }
            }
            if (ImPlot::BeginPlot("Énergie perdue selon le décalage")) {
                ImPlot::SetupAxes("décalage b (m)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (!bs.empty()) ImPlot::PlotLine("½(1−e²) μ v_n²", bs.data(), loss.data(), static_cast<int>(bs.size()), lineSpec(kVelocityColor, 0));
                const double now = balls_.p1.y;
                ImPlot::PlotInfLines("décalage actuel", &now, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue)));
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Angle entre les vitesses finales")) {
                ImPlot::SetupAxes("décalage b (m)", "degrés", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (!bsAngle.empty()) ImPlot::PlotLine("angle", bsAngle.data(), angle.data(), static_cast<int>(bsAngle.size()), lineSpec(kBall2Color, 0));
                ImPlot::EndPlot();
            }
        }

        if (etudiant && ImPlot::BeginPlot(strf("Convergence : erreur à t = %.2f s", convergenceTime_).c_str())) {
            ImPlot::SetupAxes("dt (s)", "écart (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
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

        // Accumulation de Zénon : le temps restant avant le repos tombe comme e^n (droite en échelle log).
        if (researcher && ball && ImPlot::BeginPlot("Temps restant avant le repos")) {
            ImPlot::SetupAxes("rebond n", "t_repos − t_n (s)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            const double tRest = ball_.restTime();
            std::vector<double> n, remaining;
            if (std::isfinite(tRest))
                for (std::size_t k = 0; k < impactT_.size(); ++k)
                    if (tRest - impactT_[k] > 1e-9) { n.push_back(static_cast<double>(k + 1)); remaining.push_back(tRest - impactT_[k]); }
            if (!n.empty()) ImPlot::PlotLine("reste", n.data(), remaining.data(), static_cast<int>(n.size()), markerSpec(kGround));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

void BounceModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    const double t = clock_.time;
    int lane = 0;

    if (scenario_ == kBall) {
        // Sol (ligne épaisse) et impacts exacts ; trajectoire exacte lissée (vols échantillonnés) ; une balle par méthode, par couloir.
        const float x0 = static_cast<float>(offsetX_ + ball_.x0 * scale_ - 2.0), x1 = static_cast<float>(offsetX_ + ball_.exact(duration_).x * scale_ + 2.0);
        renderer.draw(Primitive::Lines, {{x0, 0.0f, 0.0f, kGround[0], kGround[1], kGround[2]}, {x1, 0.0f, 0.0f, kGround[0], kGround[1], kGround[2]}});
        std::vector<Vertex> marks;
        for (double x : impactX_) {
            marks.push_back({static_cast<float>(x * scale_ + offsetX_), 0.0f, 0.0f, 1.0f, 1.0f, 1.0f});
            marks.push_back({static_cast<float>(x * scale_ + offsetX_), 0.25f, 0.0f, 1.0f, 1.0f, 1.0f});
        }
        renderer.draw(Primitive::Lines, marks);

        std::vector<Vertex> path;
        double ta = 0.0;
        for (std::size_t k = 0; k <= impactT_.size(); ++k) {
            const double tb = k < impactT_.size() ? impactT_[k] : std::min(duration_, ta + 1.0);
            if (tb <= ta) continue;
            for (int j = 0; j <= 24; ++j) {
                const BounceState s = ball_.exact(ta + (tb - ta) * j / 24.0);
                path.push_back({static_cast<float>(s.x * scale_ + offsetX_), static_cast<float>(s.y * scale_), 0.0f, kBlue[0] * 0.5f, kBlue[1] * 0.5f, kBlue[2] * 0.5f});
            }
            ta = tb;
            if (ta >= duration_) break;
        }
        renderer.draw(Primitive::LineStrip, path);

        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(ctx.level, i)) continue;
            ++lane;
            const Run& r = runs_[i];
            const State& y = r.ball->state();
            const float z = -0.7f * static_cast<float>(lane);
            std::vector<Vertex> trail = r.trail;
            for (Vertex& v : trail) v.z = z;
            if (!trail.empty()) renderer.draw(Primitive::LineStrip, trail);
            renderer.draw(Primitive::Points, {{static_cast<float>(y[0] * scale_ + offsetX_), static_cast<float>(y[1] * scale_ + 0.0), z,
                                               solvers_.color(i)[0], solvers_.color(i)[1], solvers_.color(i)[2]}}, 13.0f * ctx.uiScale);
        }
        const BounceState ex = ball_.exact(t);
        renderer.draw(Primitive::Points, {{static_cast<float>(ex.x * scale_ + offsetX_), static_cast<float>(ex.y * scale_), 0.0f, kBlue[0], kBlue[1], kBlue[2]}}, 18.0f * ctx.uiScale);
        if (atLeast(ctx.level, Level::Lycee) && !ex.resting) {  // vitesse (flèche verte)
            const float px = static_cast<float>(ex.x * scale_ + offsetX_), py = static_cast<float>(ex.y * scale_);
            const float dx = static_cast<float>(ex.vx) * 0.5f, dy = static_cast<float>(ex.vy) * 0.5f;
            renderer.draw(Primitive::Lines, {{px, py, 0.0f, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}, {px + dx, py + dy, 0.0f, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}});
        }
        return;
    }

    // Deux billes : vue de dessus (x, y) -> (x, 0, -y).
    const double R1 = balls_.r1, R2 = balls_.r2;
    const State ex = balls_.exact(t);
    std::vector<Vertex> p1 = {{static_cast<float>(balls_.p1.x * scale_), 0.0f, -static_cast<float>(balls_.p1.y * scale_), kBlue[0] * 0.5f, kBlue[1] * 0.5f, kBlue[2] * 0.5f}};
    std::vector<Vertex> p2 = {{static_cast<float>(balls_.p2.x * scale_), 0.0f, -static_cast<float>(balls_.p2.y * scale_), kBall2Color[0] * 0.5f, kBall2Color[1] * 0.5f, kBall2Color[2] * 0.5f}};
    const double tc = balls_.collisionTime();
    const double tEndPath = std::max(duration_, 0.0);
    if (std::isfinite(tc)) {
        const State a = balls_.exact(tc);
        p1.push_back({static_cast<float>(a[0] * scale_), 0.0f, -static_cast<float>(a[1] * scale_), p1[0].r, p1[0].g, p1[0].b});
        p2.push_back({static_cast<float>(a[2] * scale_), 0.0f, -static_cast<float>(a[3] * scale_), p2[0].r, p2[0].g, p2[0].b});
    }
    const State last = balls_.exact(tEndPath);
    p1.push_back({static_cast<float>(last[0] * scale_), 0.0f, -static_cast<float>(last[1] * scale_), p1[0].r, p1[0].g, p1[0].b});
    p2.push_back({static_cast<float>(last[2] * scale_), 0.0f, -static_cast<float>(last[3] * scale_), p2[0].r, p2[0].g, p2[0].b});
    renderer.draw(Primitive::LineStrip, p1);
    renderer.draw(Primitive::LineStrip, p2);

    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        const Run& r = runs_[i];
        const State& y = r.balls->state();
        if (!r.trail.empty()) renderer.draw(Primitive::LineStrip, r.trail);
        if (!r.trail2.empty()) renderer.draw(Primitive::LineStrip, r.trail2);
        renderer.draw(Primitive::LineStrip, circle(y[0], y[1], R1, scale_, solvers_.color(i)));
        renderer.draw(Primitive::LineStrip, circle(y[2], y[3], R2, scale_, solvers_.color(i)));
    }
    renderer.draw(Primitive::LineStrip, circle(ex[0], ex[1], R1, scale_, kBlue));
    renderer.draw(Primitive::LineStrip, circle(ex[0], ex[1], R1 * 0.93, scale_, kBlue, 0.6f));
    renderer.draw(Primitive::LineStrip, circle(ex[2], ex[3], R2, scale_, kBall2Color));
    renderer.draw(Primitive::LineStrip, circle(ex[2], ex[3], R2 * 0.93, scale_, kBall2Color, 0.6f));
    renderer.draw(Primitive::Points, {{static_cast<float>(ex[0] * scale_), 0.0f, -static_cast<float>(ex[1] * scale_), kBlue[0], kBlue[1], kBlue[2]},
                                      {static_cast<float>(ex[2] * scale_), 0.0f, -static_cast<float>(ex[3] * scale_), kBall2Color[0], kBall2Color[1], kBall2Color[2]}},
                  8.0f * ctx.uiScale);
    if (atLeast(ctx.level, Level::Lycee)) {  // vitesses (flèches vertes)
        for (int b = 0; b < 2; ++b) {
            const float px = static_cast<float>(ex[2 * b] * scale_), pz = -static_cast<float>(ex[2 * b + 1] * scale_);
            const float dx = static_cast<float>(ex[4 + 2 * b]) * 0.8f, dz = -static_cast<float>(ex[5 + 2 * b]) * 0.8f;
            renderer.draw(Primitive::Lines, {{px, 0.0f, pz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}, {px + dx, 0.0f, pz + dz, kVelocityColor[0], kVelocityColor[1], kVelocityColor[2]}});
        }
    }
}

}  // namespace pl
