#include "KeplerModule.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>

#include <implot.h>

namespace pl {
namespace {

constexpr double kTwoPi = 2.0 * constants::pi;
constexpr double kDeg = 180.0 / constants::pi;
constexpr double kAU = 1.495978707e11;     // unité astronomique [m] (définition UAI 2012)
constexpr double kGMSun = 1.3271244e20;    // GM du Soleil [m^3/s^2] (valeur nominale UAI 2015)
constexpr double kYear = 365.25 * 86400.0; // année julienne [s]
constexpr float kScale = 6.0f;             // demi-grand axe = 6 unités de scène
constexpr std::size_t kTrailLength = 1500;
constexpr double kEscapeRadius = 40.0;     // au-delà (en demi-grands axes), une méthode est déclarée « échappée »

const float kSunColor[3] = {1.0f, 0.85f, 0.30f};
const float kVelocityColor[3] = {0.35f, 1.0f, 0.45f};
const float kForceColor[3] = {1.0f, 0.38f, 0.32f};
const float kGrey[3] = {0.55f, 0.55f, 0.6f};

struct Preset {
    const char* name;
    double semiMajorAU, eccentricity;
};
const Preset kPresets[] = {{"Mercure", 0.387, 0.2056}, {"Terre", 1.0, 0.0167}, {"Mars", 1.524, 0.0934},
                           {"Halley", 17.8, 0.967},    {"Pluton", 39.48, 0.2488}};

// Plan de l'orbite (x, y) -> plan horizontal de la scène (x, 0, -y) : vue du dessus, sens trigonométrique.
Vertex pt(double x, double y, const float* c, float dim = 1.0f) {
    return {static_cast<float>(x) * kScale, 0.0f, -static_cast<float>(y) * kScale, c[0] * dim, c[1] * dim, c[2] * dim};
}

double radiusOf(const State& y) { return std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]); }

std::string formatDuration(double seconds) {
    if (seconds >= 0.5 * kYear) return strf("%.3g an%s", seconds / kYear, seconds >= 1.5 * kYear ? "s" : "");
    return strf("%.1f jours", seconds / 86400.0);
}

// Texte qui passe à la ligne dans le panneau (grisé pour les notes).
void wrapped(const std::string& text, bool dim = false) {
    if (dim) ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s", text.c_str());
    if (dim) ImGui::PopStyleColor();
}

// Flèche d'origine (x, z), de composantes (dx, dz) limitées à `maxLength` unités de scène.
void drawArrow(Renderer& renderer, float x, float z, float dx, float dz, float maxLength, const float* c, float uiScale) {
    const float len = std::sqrt(dx * dx + dz * dz);
    if (len > maxLength) { dx *= maxLength / len; dz *= maxLength / len; }
    renderer.draw(Primitive::Lines, {{x, 0.0f, z, c[0], c[1], c[2]}, {x + dx, 0.0f, z + dz, c[0], c[1], c[2]}});
    renderer.draw(Primitive::Points, {{x + dx, 0.0f, z + dz, c[0], c[1], c[2]}}, 7.0f * uiScale);
}

// Précession moyenne par orbite [rad] mesurée sur `orbits` tours (NaN si l'orbite est perdue).
double measurePrecession(const KeplerProblem& p, Solver& solver, int perOrbit, int orbits) {
    State y = p.initialState();
    const OdeFunction f = p.rhs();
    PeriapsisTracker tracker;
    const double dt = p.period() / perOrbit;
    double t = 0.0;
    for (int i = 0; i < perOrbit * orbits; ++i) {
        t += advance(solver, f, t, y, dt);
        if (!std::isfinite(y[0]) || radiusOf(y) > kEscapeRadius) return std::numeric_limits<double>::quiet_NaN();
        tracker.update(t, y);
    }
    return tracker.count() >= 5 ? tracker.precessionPerOrbit() : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace

KeplerModule::KeplerModule() {
    reset();
}

double KeplerModule::lengthUnit() const { return semiMajorAU_ * kAU; }
double KeplerModule::timeUnit() const {
    const double l = lengthUnit();
    return std::sqrt(l * l * l / (massSun_ * kGMSun));
}

const char* KeplerModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Une planète tourne autour de son étoile parce que l'étoile l'attire. Elle ne tombe pas dessus : elle "
                   "avance aussi très vite sur le côté, alors elle « tombe » sans cesse en ratant l'étoile.\n\n"
                   "Son chemin est une ellipse, un cercle un peu aplati. Près de l'étoile la planète file à toute vitesse, "
                   "loin elle ralentit. Changez la forme de l'orbite !\n\n"
                   "Le point bleu est la vérité mathématique, le vert est le calcul de l'ordinateur.";
        case Level::Interesse:
            return "Les trois lois de Kepler (1609-1619) :\n"
                   "1. Chaque planète décrit une ellipse dont le Soleil occupe un foyer.\n"
                   "2. Elle balaie des aires égales en des temps égaux : les traits bleus sont tracés à intervalles de temps "
                   "réguliers, ils sont serrés loin du Soleil et écartés près de lui, car la planète y va plus vite.\n"
                   "3. Plus une planète est loin, plus son année est longue : l'année de Mars dure 687 jours, celle de la "
                   "Terre 365.\n\n"
                   "Newton a montré que la même force, qui fait tomber une pomme, explique tout cela. Essayez les planètes "
                   "du Système solaire, et la comète de Halley (orbite très aplatie, 76 ans).";
        case Level::College:
            return "Gravitation universelle : F = G × m1 × m2 / d². Si la distance double, la force est 4 fois plus petite.\n"
                   "Troisième loi de Kepler : T² / a³ = constante. Avec T en années, a en UA et une étoile d'une masse "
                   "solaire, T² = a³ : à 4 UA on trouve T = 8 ans. Vérifiez dans « Invariants ».\n"
                   "L'excentricité e mesure l'aplatissement (0 = cercle). Distance minimale (périhélie) : a (1 − e) ; "
                   "maximale (aphélie) : a (1 + e).\n\n"
                   "La méthode simple (Euler, orange) donne un peu d'énergie à la planète à chaque pas : elle s'éloigne en "
                   "spirale, ce qui n'arrive jamais dans la réalité. La méthode précise (RK4, vert) reste sur l'ellipse.";
        case Level::Lycee:
            return "Newton : l'accélération est a = −GM r / r³ (vecteur). La flèche rouge est la force, la verte la vitesse.\n"
                   "Énergie par unité de masse : E = v²/2 − GM/r < 0 pour une orbite liée ; E = −GM / (2a). Vitesse en "
                   "chaque point : v² = GM (2/r − 1/a). Vitesse de libération : v = √(2GM / r).\n"
                   "Moment cinétique : L = r × v est constant (force centrale) : c'est la loi des aires, dA/dt = L/2.\n"
                   "Période : T = 2π √(a³ / GM).\n\n"
                   "Euler fait croître E et L à chaque tour. RK4 les conserve presque.";
        case Level::Etudiant:
            return "Deux corps se ramènent à un seul, de masse réduite μ = m1 m2/(m1+m2), dans le champ de GM = G(m1+m2) : "
                   "r'' = −GM r/|r|³.\n"
                   "Conservés : E = v²/2 − GM/r ; L = r × v ; vecteur de Runge-Lenz A = v × L − GM r/|r| (flèche grise, "
                   "|A| = GM e, dirigé vers le périastre : l'orbite est fermée). Équation polaire (Binet) : "
                   "r = p / (1 + e cos θ), p = L²/GM. Le potentiel effectif U = L²/(2r²) − GM/r borne le mouvement entre "
                   "périastre et apoastre.\n"
                   "Solution exacte : M = E_a − e sin E_a (Kepler), résolue par Newton ; x = a (cos E_a − e), "
                   "y = a √(1−e²) sin E_a.\n"
                   "Euler symplectique et Verlet gardent une énergie bornée et L exact, mais FONT TOURNER l'ellipse (le "
                   "périastre précesse, flèche grise) ; RK4 dissipe lentement ; RK45 raccourcit son pas au périastre.";
        case Level::Chercheur:
            return "Verlet est le flot exact d'un hamiltonien modifié H~ = H + h² (1/12 pᵀ V'' p − 1/24 |grad V|²) + O(h⁴), "
                   "V'' étant la matrice des dérivées secondes de V : "
                   "d'où énergie bornée (sur des temps exponentiellement longs) mais une précession séculaire. La moyenne "
                   "de la perturbation sur l'orbite, injectée dans l'équation de Lagrange dω/dt = ∂<Δ>/∂G (variables de "
                   "Delaunay), donne par orbite :\n"
                   "Δω = −(π/8) (GM h²/a³) (4 + e²) / (1 − e²)³\n"
                   "(rétrograde, en h² ; comparée à la mesure dans « Invariants » et dans l'Analyse). Euler symplectique est "
                   "conjugué au schéma de Strang : même précession à l'ordre dominant, mesurée identique ici.\n"
                   "Forte excentricité : l'échelle de temps locale au périastre est ~ r^(3/2), ce qu'un pas fixe ne résout "
                   "pas ; RK45 adapte h (graphe « Pas RK45 »). Remèdes de la littérature : régularisation de "
                   "Kustaanheimo-Stiefel, changement de temps de Sundman (dt = r dE_a/(na)).";
    }
    return "";
}

void KeplerModule::frameCamera(Camera& camera) const {
    camera.target[0] = -3.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 24.0f;
    camera.yaw = 0.0f;
    camera.pitch = 1.0f;
}

// ------------------------------ simulation -----------------------------

void KeplerModule::buildGeometry() {
    ellipse_.clear();
    const double b = std::sqrt(1.0 - ecc_ * ecc_);
    for (int i = 0; i <= 240; ++i) {
        const double E = kTwoPi * i / 240.0;
        ellipse_.push_back(pt(std::cos(E) - ecc_, b * std::sin(E), kBlue, 0.55f));
    }
    // 12 rayons à intervalles de temps égaux (T/12) : la loi des aires se lit sur leur écartement.
    sectors_.clear();
    for (int k = 0; k < 12; ++k) {
        const State s = problem_.exact(problem_.period() * k / 12.0);
        sectors_.push_back(pt(0.0, 0.0, kGrey, 0.5f));
        sectors_.push_back(pt(s[0], s[1], kGrey, 0.5f));
    }
}

void KeplerModule::reset() {
    problem_.mu = 1.0;
    problem_.a = 1.0;
    problem_.e = ecc_;
    rhs_ = problem_.rhs();
    e0_ = problem_.exactEnergy();

    for (Run& r : runs_) {
        r.y = problem_.initialState();
        r.radius.clear();
        r.energy.clear();
        r.momentum.clear();
        r.stepSize.clear();
        r.precession.clear();
        r.trail.clear();
        r.tracker.reset();
        r.tracker.update(0.0, r.y);
        r.diverged = false;
        r.lastStep = 0.0;
    }
    exactRadius_.clear();
    exactKinetic_.clear();
    exactPotential_.clear();
    exactTotal_.clear();
    solvers_.rk45().resetStats();
    buildGeometry();

    clock_.reset();
    running_ = true;
    convergenceDirty_ = true;
    scanDirty_ = true;
    sample();
}

void KeplerModule::sample() {
    const double t = clock_.time, orbits = t / kTwoPi;
    const State ex = problem_.exact(t);
    const double rx = radiusOf(ex), scale = 1.0 / std::abs(e0_);
    exactRadius_.add(orbits, rx * semiMajorAU_);
    const double kinetic = 0.5 * (ex[3] * ex[3] + ex[4] * ex[4] + ex[5] * ex[5]);
    exactKinetic_.add(orbits, kinetic * scale);
    exactPotential_.add(orbits, -1.0 / rx * scale);
    exactTotal_.add(orbits, (kinetic - 1.0 / rx) * scale);

    for (int i = 0; i < SolverSet::kCount; ++i) {
        Run& r = runs_[i];
        if (r.diverged) continue;
        r.radius.add(orbits, radiusOf(r.y) * semiMajorAU_);
        r.energy.add(orbits, (problem_.energy(r.y) - e0_) * scale);
        r.momentum.add(orbits, problem_.angularMomentum(r.y).z / problem_.exactAngularMomentum() - 1.0);
        if (i == SolverSet::kRK45) r.stepSize.add(orbits, r.lastStep);
        r.trail.push_back(pt(r.y[0], r.y[1], solvers_.color(i), 0.7f));
        if (r.trail.size() > kTrailLength) r.trail.erase(r.trail.begin(), r.trail.begin() + 200);
    }
}

// Avance un solveur de h. RK45 est piloté ici pour garder la trace de son dernier pas (graphe « Pas RK45 »).
void KeplerModule::advanceRun(int i, Run& run, double t, double h) {
    if (i != SolverSet::kRK45) {
        advance(solvers_.solver(i), rhs_, t, run.y, h);
        return;
    }
    RK45& rk = solvers_.rk45();
    double done = 0.0;
    while (h - done > 1e-12 * h) {
        const double step = rk.step(rhs_, t + done, run.y, h - done);
        if (!(step > 0.0)) break;
        run.lastStep = step;
        done += step;
    }
}

void KeplerModule::update(double frameSeconds) {
    if (!running_ || clock_.finished) return;

    solvers_.rk45().relTol = relTol_;
    solvers_.rk45().absTol = relTol_ * 1e-2;
    const double dt = kTwoPi / stepsPerOrbit_;
    clock_.advance(
        frameSeconds, timeScale_, dt, orbits_ * kTwoPi,
        [&](double h) {
            for (int i = 0; i < SolverSet::kCount; ++i) {
                Run& r = runs_[i];
                if (r.diverged) continue;
                advanceRun(i, r, clock_.time, h);
                if (!std::isfinite(r.y[0]) || !std::isfinite(r.y[1]) || radiusOf(r.y) > kEscapeRadius) {
                    r.diverged = true;
                } else if (r.tracker.update(clock_.time + h, r.y)) {
                    r.precession.add((clock_.time + h) / kTwoPi, r.tracker.precession() * kDeg);
                }
            }
        },
        [&] { sample(); });
    if (clock_.finished) running_ = false;
}

// Erreur de chaque solveur à pas fixe après 0,35 période (les petites erreurs se mesurent plus vite et Euler y est
// asymptotique), pour différents pas. La pente est ajustée là où l'erreur reste sous 0,1.
void KeplerModule::computeConvergence() {
    static const int kStepCounts[] = {50, 100, 200, 400, 800, 1600, 3200, 6400};
    const double tEnd = 0.35 * problem_.period();

    for (int i = 0; i < SolverSet::kFixedStep; ++i) {
        Curve& c = convergence_[i];
        c.x.clear();
        c.y.clear();
        const std::unique_ptr<Solver> solver = SolverSet::makeFixedStep(i);
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int n : kStepCounts) {
            const double e = keplerError(problem_, *solver, n, tEnd);
            if (!std::isfinite(e) || e <= 1e-12) continue;  // sous ce seuil : arrondi
            c.x.push_back(tEnd / n);
            c.y.push_back(e);
            if (e < 0.1) {
                const double x = std::log10(tEnd / n), y = std::log10(e);
                sx += x; sy += y; sxx += x * x; sxy += x * y;
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

// Précession par orbite de Euler symplectique et Verlet selon le pas, face à la formule du hamiltonien modifié.
void KeplerModule::computePrecessionScan() {
    static const int kPerOrbit[] = {50, 100, 200, 400, 800};
    for (Curve* c : {&scanSymplectic_, &scanVerlet_, &scanTheory_}) { c->x.clear(); c->y.clear(); }
    SymplecticEuler symplectic;
    VelocityVerlet verlet;
    for (int n : kPerOrbit) {
        const double dt = problem_.period() / n;
        const double s = measurePrecession(problem_, symplectic, n, 30), v = measurePrecession(problem_, verlet, n, 30);
        if (std::isfinite(s) && s != 0.0) { scanSymplectic_.x.push_back(dt); scanSymplectic_.y.push_back(std::abs(s) * kDeg); }
        if (std::isfinite(v) && v != 0.0) { scanVerlet_.x.push_back(dt); scanVerlet_.y.push_back(std::abs(v) * kDeg); }
        scanTheory_.x.push_back(dt);
        scanTheory_.y.push_back(std::abs(problem_.verletPrecessionPerOrbit(dt)) * kDeg);
    }
    scanDirty_ = false;
}

// --------------------------------- UI ----------------------------------

void KeplerModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    bool changed = false;

    ImGui::SeparatorText("Orbite");
    if (interesse) {
        ImGui::TextDisabled("Astres du Système solaire");
        for (std::size_t i = 0; i < std::size(kPresets); ++i) {
            if (i > 0) ImGui::SameLine();
            if (ImGui::SmallButton(kPresets[i].name)) {
                semiMajorAU_ = kPresets[i].semiMajorAU;
                ecc_ = kPresets[i].eccentricity;
                changed = true;
            }
        }
    }
    pushSliderWidth();
    changed |= sliderD(interesse ? "Excentricité e (0 = cercle)" : "Forme de l'orbite", &ecc_, 0.0, 0.97,
                       interesse ? "%.3f" : "");
    if (interesse) changed |= sliderD("Demi-grand axe a (UA)", &semiMajorAU_, 0.3, 50.0, "%.2f", ImGuiSliderFlags_Logarithmic);
    if (college) changed |= sliderD("Masse étoile (masses sol.)", &massSun_, 0.2, 5.0, "%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= sliderD("Pas par orbite", &stepsPerOrbit_, 20.0, 3000.0, "%.0f", ImGuiSliderFlags_Logarithmic);
        if (atLeast(level, Level::Etudiant)) changed |= sliderD("Durée (orbites)", &orbits_, 1.0, 100.0, "%.0f");
        if (atLeast(level, Level::Chercheur))
            changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);
        popSliderWidth();
        ImGui::TextDisabled("dt = T / %.0f = %s", stepsPerOrbit_, formatDuration(kTwoPi * timeUnit() / stepsPerOrbit_).c_str());
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

    if (college) {
        const double orbitsDone = clock_.time / kTwoPi;
        ImGui::Text("t = %s   (%.2f orbite)", formatDuration(clock_.time * timeUnit()).c_str(), orbitsDone);
        if (atLeast(level, Level::Etudiant) && solvers_.show(SolverSet::kRK45))
            ImGui::Text("RK45 : %d pas acceptés, %d rejetés", solvers_.rk45().acceptedSteps(), solvers_.rk45().rejectedSteps());
    }
}

void KeplerModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant que la planète tourne, deux choses restent identiques : son énergie totale (mouvement + "
                           "attraction) et sa « quantité de rotation » autour de l'étoile. C'est pourquoi l'orbite se "
                           "referme toujours sur elle-même.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool researcher = atLeast(level, Level::Chercheur);

    // Orbite exacte en unités réelles.
    const double l = lengthUnit(), tau = timeUnit(), T = kTwoPi * tau, gm = massSun_ * kGMSun;
    const double rp = semiMajorAU_ * (1.0 - ecc_), ra = semiMajorAU_ * (1.0 + ecc_);
    const double vp = std::sqrt((1.0 + ecc_) / (1.0 - ecc_)) * l / tau / 1000.0;
    const double va = std::sqrt((1.0 - ecc_) / (1.0 + ecc_)) * l / tau / 1000.0;
    const double years = T / kYear;
    ImGui::SeparatorText("Orbite exacte");
    wrapped(strf("Période T = %s", formatDuration(T).c_str()));
    wrapped(strf("Périhélie %.3f UA (%.1f km/s) ; aphélie %.3f UA (%.1f km/s)", rp, vp, ra, va));
    wrapped(strf("T² × M / a³ = %.3f  (T en ans, a en UA, M en masses solaires)",
                 years * years * massSun_ / std::pow(semiMajorAU_, 3)));
    if (lycee) {
        wrapped(strf("E/m = −GM/2a = %.1f MJ/kg", -gm / (2.0 * l) / 1e6));
        wrapped(strf("L/m = √(GM a (1−e²)) = %.3e m²/s", std::sqrt(gm * l * (1.0 - ecc_ * ecc_))));
        wrapped(strf("Vitesse de libération au périhélie : %.1f km/s", std::sqrt(2.0 * gm / (rp * kAU)) / 1000.0));
    }

    // Tableaux étroits : à partir du niveau 5, ce qui se conserve d'un côté, la précision de la trajectoire de l'autre.
    const bool split = etudiant;
    std::vector<std::string> headersA, headersB;
    headersA.push_back(etudiant ? "ΔE/E0" : "ΔE %");
    if (lycee) headersA.push_back(etudiant ? "ΔL/L0" : "ΔL %");
    if (etudiant) headersA.push_back("ΔA/μ");
    (split ? headersB : headersA).push_back(etudiant ? "dist" : "écart");
    if (researcher) headersB.push_back("préc.°");

    std::vector<TableRow> rowsA, rowsB;
    rowsA.push_back({"Exacte", kBlue, std::vector<std::string>(headersA.size(), "-")});
    if (split) {
        rowsB.push_back({"Exacte", kBlue, std::vector<std::string>(headersB.size(), "-")});
        if (researcher) rowsB.back().cells.back() = "0";
    }

    const State exactNow = problem_.exact(clock_.time);
    const double lz0 = problem_.exactAngularMomentum();
    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(level, i)) continue;
        const Run& r = runs_[i];
        TableRow a{solvers_.shortLabel(level, i), solvers_.color(i), {}};
        TableRow b = a;
        if (r.diverged) {
            a.cells.assign(headersA.size(), "-");
            a.cells[0] = "échappe";
            b.cells.assign(headersB.size(), "-");
        } else {
            const double dE = (problem_.energy(r.y) - e0_) / std::abs(e0_);
            const double dL = problem_.angularMomentum(r.y).z / lz0 - 1.0;
            const double posErr = std::sqrt(std::pow(r.y[0] - exactNow[0], 2) + std::pow(r.y[1] - exactNow[1], 2));
            a.cells.push_back(strf("%+.1e", etudiant ? dE : 100.0 * dE));
            if (lycee) a.cells.push_back(strf("%+.1e", etudiant ? dL : 100.0 * dL));
            if (etudiant) a.cells.push_back(strf("%+.1e", problem_.rungeLenz(r.y).norm() - ecc_));
            const std::string err = strf("%.1e", etudiant ? problem_.distance(r.y, exactNow) : posErr * semiMajorAU_);
            (split ? b : a).cells.push_back(err);
            if (researcher) {
                const double p = r.tracker.precessionPerOrbit();
                b.cells.push_back(std::isfinite(p) ? strf("%+.4f", p * kDeg) : "-");
            }
        }
        rowsA.push_back(a);
        if (split) rowsB.push_back(b);
    }
    if (researcher) {
        TableRow th{"Théorie Verlet", nullptr, std::vector<std::string>(headersB.size(), "-")};
        th.cells.back() = strf("%+.4f", problem_.verletPrecessionPerOrbit(kTwoPi / stepsPerOrbit_) * kDeg);
        rowsB.push_back(th);
    }

    ImGui::SeparatorText(split ? "Ce qui se conserve" : "Comparaison");
    drawResultTable("conservation", headersA, rowsA);
    if (split) {
        ImGui::SeparatorText("Précision de l'orbite");
        drawResultTable("precision", headersB, rowsB);
    }
    wrapped(split ? "ΔA/μ = |A|/μ − e (Runge-Lenz). dist : écart à l'orbite exacte dans l'espace des phases."
                  : "écart : distance à la planète exacte, en UA. ΔE, ΔL : variation en %.",
            true);
    if (researcher) wrapped("préc.° : rotation du périastre par orbite (°/orbite), mesurée sur les passages au périastre.", true);
    wrapped("« - » : non défini ; « échappe » : la méthode a quitté l'orbite.", true);
}

void KeplerModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showEnergy = atLeast(level, Level::College);
    const bool showMomentum = atLeast(level, Level::Lycee);
    const int cols = 1 + (showEnergy ? 1 : 0) + (showMomentum ? 1 : 0);

    auto plotRuns = [&](Series Run::*member) {
        for (int i = 0; i < SolverSet::kCount; ++i) {
            if (!solvers_.isShown(level, i) || (runs_[i].*member).size() == 0) continue;
            const Series& s = runs_[i].*member;
            ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(solvers_.color(i), s.offset));
        }
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Distance à l'étoile r(t)")) {
            ImPlot::SetupAxes("t (orbites)", "r (UA)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (exactRadius_.size() > 0)
                ImPlot::PlotLine("Exacte", exactRadius_.x.data(), exactRadius_.y.data(), exactRadius_.size(), lineSpec(kBlue, exactRadius_.offset));
            plotRuns(&Run::radius);
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie (E − E0) / |E0|")) {
            ImPlot::SetupAxes("t (orbites)", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotRuns(&Run::energy);
            ImPlot::EndPlot();
        }
        if (showMomentum && ImPlot::BeginPlot("Moment cinétique L / L0 − 1")) {
            ImPlot::SetupAxes("t (orbites)", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            plotRuns(&Run::momentum);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void KeplerModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool researcher = atLeast(level, Level::Chercheur);
    if (etudiant && convergenceDirty_) computeConvergence();
    if (researcher && scanDirty_) computePrecessionScan();

    auto markerSpec = [](const float* c) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle,
                          ImPlotProp_MarkerSize, 4.0f);
    };

    // Page 1 : énergies (niveau 4), puis potentiel effectif et convergence (niveau 5).
    auto pageConservation = [&] {
        if (!ImPlot::BeginSubplots("##analyse", 1, etudiant ? 3 : 1, ImVec2(-1.0f, -1.0f))) return;
        // Ec + Ep = E reste constante, les deux autres s'échangent.
        if (ImPlot::BeginPlot("Énergies de la planète (exacte)")) {
            ImPlot::SetupAxes("t (orbites)", "énergie / |E|", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            const float red[3] = {1.0f, 0.4f, 0.4f};
            if (exactKinetic_.size() > 0) {
                ImPlot::PlotLine("Cinétique", exactKinetic_.x.data(), exactKinetic_.y.data(), exactKinetic_.size(), lineSpec(kVelocityColor, exactKinetic_.offset));
                ImPlot::PlotLine("Potentielle (attraction)", exactPotential_.x.data(), exactPotential_.y.data(), exactPotential_.size(), lineSpec(red, exactPotential_.offset));
                ImPlot::PlotLine("Totale", exactTotal_.x.data(), exactTotal_.y.data(), exactTotal_.size(), lineSpec(kBlue, exactTotal_.offset));
            }
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot("Potentiel effectif  U(r) = L²/2r² − GM/r")) {
            // en unités de |E0| ; r en UA ; le mouvement se limite à U <= E, entre périastre et apoastre
            const double L2 = std::pow(problem_.exactAngularMomentum(), 2);
            const double rp = 1.0 - ecc_, ra = 1.0 + ecc_, scale = 1.0 / std::abs(e0_);
            std::vector<double> xs, ys;
            for (int k = 0; k < 300; ++k) {
                const double r = 0.4 * rp + (1.25 * ra - 0.4 * rp) * k / 299.0;
                xs.push_back(r * semiMajorAU_);
                ys.push_back((0.5 * L2 / (r * r) - 1.0 / r) * scale);
            }
            const double uMin = -0.5 / L2 * scale;  // fond du puits : -GM²/(2 L²)
            ImPlot::SetupAxes("r (UA)", "U / |E|", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 1.35 * uMin, -0.55 * uMin, ImPlotCond_Always);
            ImPlot::PlotLine("U(r)", xs.data(), ys.data(), static_cast<int>(xs.size()), lineSpec(kVelocityColor, 0));
            const double eLevel = e0_ * scale;
            ImPlot::PlotInfLines("E", &eLevel, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            const double bounds[2] = {rp * semiMajorAU_, ra * semiMajorAU_};
            ImPlot::PlotInfLines("périastre / apoastre", bounds, 2, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kGrey)));
            const State ex = problem_.exact(clock_.time);
            const double cx = radiusOf(ex) * semiMajorAU_, cy = eLevel;
            ImPlot::PlotScatter("Planète", &cx, &cy, 1, markerSpec(kBlue));
            ImPlot::EndPlot();
        }

        if (etudiant && ImPlot::BeginPlot("Convergence : erreur à t = 0,35 T")) {
            ImPlot::SetupAxes("dt (normalisé)", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
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
    };

    // Page 2 (niveau 6) : pas adaptatif de RK45 et précession des schémas à pas fixe.
    auto pagePrecession = [&] {
        if (!ImPlot::BeginSubplots("##precession", 1, 3, ImVec2(-1.0f, -1.0f))) return;
        if (ImPlot::BeginPlot("Pas de RK45  h(t)")) {
            ImPlot::SetupAxes("t (orbites)", "h (normalisé)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            const Series& s = runs_[SolverSet::kRK45].stepSize;
            if (s.size() > 0) ImPlot::PlotLine("RK45", s.x.data(), s.y.data(), s.size(), lineSpec(solvers_.color(SolverSet::kRK45), s.offset));
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Précession : angle du périastre")) {
            ImPlot::SetupAxes("t (orbites)", "rotation (°)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            for (int i = 0; i < SolverSet::kCount; ++i) {
                const Series& s = runs_[i].precession;
                if (!solvers_.isShown(level, i) || s.size() == 0) continue;
                ImPlot::PlotLine(solvers_.label(level, i).c_str(), s.x.data(), s.y.data(), s.size(), lineSpec(solvers_.color(i), s.offset));
            }
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Précession par orbite selon le pas")) {
            ImPlot::SetupAxes("dt (normalisé)", "|rotation| (°/orbite)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (scanTheory_.x.size() >= 2)
                ImPlot::PlotLine("Théorie (en h²)", scanTheory_.x.data(), scanTheory_.y.data(), static_cast<int>(scanTheory_.x.size()), lineSpec(kGrey, 0));
            if (scanSymplectic_.x.size() >= 2)
                ImPlot::PlotLine("Euler symplectique", scanSymplectic_.x.data(), scanSymplectic_.y.data(), static_cast<int>(scanSymplectic_.x.size()), markerSpec(solvers_.color(SolverSet::kSymplectic)));
            if (scanVerlet_.x.size() >= 2)
                ImPlot::PlotLine("Verlet", scanVerlet_.x.data(), scanVerlet_.y.data(), static_cast<int>(scanVerlet_.x.size()), markerSpec(solvers_.color(SolverSet::kVerlet)));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    };

    if (!researcher) {
        pageConservation();
    } else if (ImGui::BeginTabBar("##analyseTabs")) {
        if (ImGui::BeginTabItem("Conservation et convergence")) { pageConservation(); ImGui::EndTabItem(); }
        if (ImGui::BeginTabItem("Pas adaptatif et précession")) { pagePrecession(); ImGui::EndTabItem(); }
        ImGui::EndTabBar();
    }
}

// ------------------------------ rendu 3D --------------------------------

void KeplerModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    renderer.draw(Primitive::LineStrip, ellipse_);
    if (atLeast(ctx.level, Level::Interesse)) renderer.draw(Primitive::Lines, sectors_);
    renderer.draw(Primitive::Points, {pt(0.0, 0.0, kSunColor)}, 28.0f * ctx.uiScale);

    for (int i = 0; i < SolverSet::kCount; ++i) {
        if (!solvers_.isShown(ctx.level, i)) continue;
        const Run& r = runs_[i];
        if (r.diverged) continue;
        if (!r.trail.empty()) renderer.draw(Primitive::LineStrip, r.trail);
        renderer.draw(Primitive::Points, {pt(r.y[0], r.y[1], solvers_.color(i))}, 11.0f * ctx.uiScale);

        // Vecteur de Runge-Lenz de la méthode (niveau 5+) : de l'étoile vers le périastre, de longueur e.
        // Pour l'orbite exacte il est fixe ; sa rotation est la précession.
        if (atLeast(ctx.level, Level::Etudiant)) {
            const Vec3 A = problem_.rungeLenz(r.y);
            renderer.draw(Primitive::Lines, {pt(0.0, 0.0, solvers_.color(i), 0.6f), pt(A.x, A.y, solvers_.color(i), 0.6f)});
        }
    }

    const State ex = problem_.exact(clock_.time);
    renderer.draw(Primitive::Points, {pt(ex[0], ex[1], kBlue)}, 15.0f * ctx.uiScale);

    // Niveau 4+ : vitesse (verte) et force (rouge) de la planète exacte. Longueurs proportionnelles, plafonnées.
    if (atLeast(ctx.level, Level::Lycee)) {
        const double r = radiusOf(ex), pull = 1.0 / (r * r * r);
        const float x = static_cast<float>(ex[0]) * kScale, z = -static_cast<float>(ex[1]) * kScale;
        drawArrow(renderer, x, z, static_cast<float>(ex[3]), -static_cast<float>(ex[4]), 7.0f, kVelocityColor, ctx.uiScale);
        drawArrow(renderer, x, z, static_cast<float>(-ex[0] * pull) * 0.5f, static_cast<float>(ex[1] * pull) * 0.5f, 6.0f, kForceColor, ctx.uiScale);
    }
}

}  // namespace pl
