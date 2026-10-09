#include "GpuNBodyModule.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include <implot.h>

namespace pl {
namespace {

using Clock = std::chrono::steady_clock;

constexpr float kScale = 4.0f;           // 1 unité de longueur = 4 unités de scène
constexpr double kBudgetMs = 14.0;       // temps de calcul visé par image : l'interface reste fluide
constexpr int kMaxStepsPerFrame = 2000;
constexpr int kMinBodies = 100;
constexpr int kMaxBodiesGpu = 20000;
constexpr int kMaxBodiesCpu = 4000;      // le CPU double (2 évaluations en moins d'une image) ne suit plus au-delà
constexpr int kExactEnergyMaxN = 4000;   // jusque-là E0 est calculée en double sur le CPU (O(N²) une fois)
constexpr int kAutoValidationMaxN = 1500;
const float kGpuColor[3] = {1.0f, 0.60f, 0.15f};   // orange : carte graphique
const float kCpuColor[3] = {0.35f, 0.65f, 1.0f};   // bleu : processeur (même bleu que la « solution exacte » des autres modules)
const float kGrey[3] = {0.55f, 0.55f, 0.6f};

double elapsedSeconds(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

double floorForLog(double v) { return std::max(v, 1e-18); }

// « 4950 paires », « 2 millions de paires », « 12 milliards de paires ».
std::string pairsText(double v) {
    if (v >= 1e9) return strf("%.2g milliards de paires", v * 1e-9);
    if (v >= 1e6) return strf("%.2g millions de paires", v * 1e-6);
    return strf("%.0f paires", v);
}

bool allFinite(const State& y) {
    for (double v : y) if (!std::isfinite(v)) return false;
    return true;
}

}  // namespace

void GpuNBodyModule::Track::clear() {
    y.clear();
    energy.clear();
    momentum.clear();
    angular.clear();
    stepMs.clear();
    stepMillis = kernelMillis = lastEnergy = lastMomentum = lastAngular = 0.0;
    diverged = false;
}

GpuNBodyModule::GpuNBodyModule(std::string shaderDir) : shaderDir_(std::move(shaderDir)) {
    gpuOk_ = gpu_.init(shaderDir_);
    if (!gpuOk_) backend_ = kCpu;
    lastCpuEnergy_ = lastValidation_ = Clock::now();
    reset();
}

GpuNBodyModule::~GpuNBodyModule() { gpu_.shutdown(); }

int GpuNBodyModule::maxBodies() const { return backend_ == kGpu && gpuOk_ ? kMaxBodiesGpu : kMaxBodiesCpu; }

const char* GpuNBodyModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Chaque étoile attire toutes les autres. Avec 2000 étoiles, cela fait deux millions de paires à calculer, "
                   "et il faut tout recommencer à chaque instant !\n\n"
                   "Le processeur de l'ordinateur calcule ces paires presque une par une. La carte graphique, celle des jeux "
                   "vidéo, a des milliers de petits calculateurs qui travaillent en même temps.\n\n"
                   "Essayez : augmentez le nombre d'étoiles avec le processeur, puis avec la carte graphique, et regardez "
                   "jusqu'où le mouvement reste fluide.";
        case Level::Interesse:
            return "Un amas d'étoiles est un « problème à N corps ». Pour N étoiles il y a N(N−1)/2 paires : doubler le nombre "
                   "d'étoiles multiplie le travail par quatre. Aucune formule ne donne la réponse, on calcule pas à pas.\n\n"
                   "La carte graphique (GPU) confie une étoile à chaque calculateur : tous avancent ensemble. Contrepartie : elle "
                   "calcule avec des nombres moins précis (7 chiffres au lieu de 16). Pour dessiner un amas c'est suffisant, mais "
                   "deux calculs de précisions différentes finissent par se séparer : c'est le chaos.\n\n"
                   "« Deux amas » : une collision, comme entre deux galaxies.";
        case Level::College:
            return "Chaque paire d'étoiles s'attire avec F = G × m1 × m2 / d². Une étoile subit la somme des attractions de toutes "
                   "les autres : pour N étoiles, N(N−1)/2 paires (2000 étoiles : environ 2 millions).\n"
                   "L'ordinateur avance le temps par petits pas dt : à chaque pas il recalcule toutes les forces, puis déplace les "
                   "étoiles. Le « temps par pas » (en millisecondes) dit qui est le plus rapide.\n"
                   "Le processeur calcule avec 16 chiffres (double), la carte graphique avec 7 (float). Regardez la colonne ΔE : "
                   "l'énergie totale doit rester constante, elle bouge à peine.";
        case Level::Lycee:
            return "Accélération de l'étoile i : a_i = G Σ_j m_j (r_j − r_i) / (|r_j − r_i|² + ε²)^(3/2). L'adoucissement ε évite "
                   "une force infinie quand deux étoiles passent très près l'une de l'autre.\n"
                   "Schéma « saute-mouton » : v ← v + a·dt/2 ; x ← x + v·dt ; on recalcule a ; v ← v + a·dt/2. Il conserve très "
                   "bien l'énergie (elle oscille sans dériver).\n"
                   "Trois grandeurs doivent rester constantes : l'énergie E = ½ Σ m v² − Σ_{i<j} G m_i m_j / r_ij, l'impulsion "
                   "P = Σ m v et le moment cinétique L = Σ m r × v. Leurs écarts mesurent l'erreur du calcul.";
        case Level::Etudiant:
            return "Compute shader : un thread par étoile cible i ; les sources j sont lues par tuiles de 256 corps chargées en "
                   "mémoire partagée, donc chaque étoile n'est lue qu'une fois par groupe et non une fois par thread. Coût O(N²), "
                   "sans arbre (Barnes-Hut) ni approximation.\n"
                   "Intégrateur kick-drift-kick (Verlet des vitesses), d'ordre 2 et symplectique, entièrement sur le GPU : "
                   "positions, vitesses et accélérations restent dans des tampons, seul l'affichage relit les positions.\n"
                   "Le mode « les deux » fait tourner le même schéma, avec le même pas, en double sur le CPU. L'écart entre les deux "
                   "(espace des phases) croît comme e^(λt) avec λ de l'ordre de 1 : l'arrondi float de 1e-7 est amplifié par le "
                   "chaos. On compare donc les accélérations à état fixé (panneau Analyse), pas les trajectoires longues.";
        case Level::Chercheur:
            return "Précision : float a une mantisse de 24 bits (6e-8). Une somme de N termes arrondis donne une erreur relative de "
                   "l'ordre de u√N : mesuré ici (Intel UHD) 1,6e-8 à N = 3 et 1,05e-6 à N = 5000 sur les accélérations. Les positions "
                   "sont centrées sur le barycentre en double avant la conversion : sans cela un amas décalé de 1000 perdrait un "
                   "facteur 345 de précision (mesuré).\n"
                   "Les forces par paire ne sont opposées qu'à l'arrondi : l'impulsion dérive (mesuré 5e-9 à 2e-6 selon N et le nombre "
                   "de pas) alors qu'elle est exacte à l'arrondi en double (1e-16). L'énergie, elle, dérive autant que l'erreur du "
                   "schéma d'ordre 2 : le float n'y ajoute presque rien. Le potentiel est sommé dans le même shader (terme m_j/r déjà "
                   "calculé), puis U = ½ Σ m_i φ_i en double.\n"
                   "Le double du shader existe mais il est émulé par le pilote de cette carte : de mille à plusieurs milliers de fois "
                   "plus lent que le float selon N (mesuré). Les envois sont découpés (~50 ms) pour ne pas déclencher la "
                   "réinitialisation du pilote (TDR) de Windows. Le banc d'essai (panneau Analyse) mesure débit et erreur en "
                   "fonction de N.";
    }
    return "";
}

void GpuNBodyModule::frameCamera(Camera& camera) const {
    camera.target[0] = camera.target[1] = camera.target[2] = 0.0f;
    camera.distance = 22.0f;
    camera.yaw = 0.0f;
    camera.pitch = 0.9f;
}

// ------------------------------ simulation -----------------------------

void GpuNBodyModule::buildProblem() {
    const int n = static_cast<int>(std::lround(bodies_));
    problem_ = scenario_ == kCluster ? NBodyProblem::plummer(n, seed_, 1.0, softening_)
                                     : NBodyProblem::plummerCollision(n, seed_, 6.0, 0.4, 0.6, softening_);
}

void GpuNBodyModule::reset() {
    pendingReset_ = false;
    bodies_ = std::clamp(bodies_, static_cast<double>(kMinBodies), static_cast<double>(maxBodies()));
    buildProblem();
    const int n = bodyCount();
    const State y0 = problem_.initialState();

    p0_ = problem_.momentum(y0);
    l0_ = problem_.angularMomentum(y0);
    pScale_ = lScale_ = 0.0;
    double v2 = 0.0;
    for (int i = 0; i < n; ++i) {
        const Vec3 r{y0[3 * i], y0[3 * i + 1], y0[3 * i + 2]}, v{y0[3 * n + 3 * i], y0[3 * n + 3 * i + 1], y0[3 * n + 3 * i + 2]};
        pScale_ += problem_.mass[i] * v.norm();
        lScale_ += problem_.mass[i] * cross(r, v).norm();
        v2 += problem_.mass[i] * v.norm2();
    }
    pScale_ = std::max(pScale_, 1e-300);
    lScale_ = std::max(lScale_, 1e-300);
    speedScale_ = static_cast<float>(2.0 * std::sqrt(v2 / std::max(problem_.totalMass(), 1e-300)));  // 2 x vitesse quadratique moyenne

    time_ = accumulator_ = 0.0;
    achieved_ = 1.0;
    gpuTrack_.clear();
    cpuTrack_.clear();
    gap_.clear();
    validation_ = Validation{};
    gpuTrack_.y = cpuTrack_.y = y0;
    cpuAccValid_ = false;

    double gpuPotential = 0.0;
    if (useGpu()) {
        gpu_.setState(y0.data(), y0.data() + 3 * n, problem_.mass.data(), n);
        gpu_.compute(problem_.G, problem_.softening);  // a(x0) : le premier pas ne le recalcule pas
        gpuPotential = gpu_.potentialEnergy();
    }
    exactEnergy_ = n <= kExactEnergyMaxN || !useGpu();
    const double u0 = exactEnergy_ ? problem_.potentialEnergy(y0) : gpuPotential;
    e0_ = problem_.kineticEnergy(y0) + u0;

    if (useGpu()) {
        recordConserved(gpuTrack_);
        recordEnergy(gpuTrack_, gpuPotential);
    }
    if (useCpu()) {
        recordConserved(cpuTrack_);
        recordEnergy(cpuTrack_, u0);
    }
    lastCpuEnergy_ = lastValidation_ = Clock::now();
}

void GpuNBodyModule::recordConserved(Track& t) {
    t.lastMomentum = floorForLog((problem_.momentum(t.y) - p0_).norm() / pScale_);
    t.lastAngular = floorForLog((problem_.angularMomentum(t.y) - l0_).norm() / lScale_);
    if (time_ > 0.0) {  // à t = 0 l'écart est nul : sur un axe logarithmique ce point écraserait toute la courbe
        t.momentum.add(time_, t.lastMomentum);
        t.angular.add(time_, t.lastAngular);
    }
}

void GpuNBodyModule::recordEnergy(Track& t, double potential) {
    t.lastEnergy = (problem_.kineticEnergy(t.y) + potential - e0_) / std::abs(e0_);
    t.energy.add(time_, t.lastEnergy);
}

void GpuNBodyModule::stepGpu(int steps) {
    const int n = bodyCount();
    const auto t0 = Clock::now();
    gpu_.step(problem_.G, problem_.softening, dt_, steps);
    const double ms = 1e3 * elapsedSeconds(t0, Clock::now()) / steps;
    gpuTrack_.stepMillis = gpuTrack_.stepMillis == 0.0 ? ms : 0.7 * gpuTrack_.stepMillis + 0.3 * ms;
    gpuTrack_.stepMs.add(time_, ms);
    const double kernel = 1e3 * gpu_.lastKernelSeconds() / steps;
    gpuTrack_.kernelMillis = gpuTrack_.kernelMillis == 0.0 ? kernel : 0.7 * gpuTrack_.kernelMillis + 0.3 * kernel;

    gpuTrack_.y.resize(6 * static_cast<std::size_t>(n));
    gpu_.readState(gpuTrack_.y.data(), gpuTrack_.y.data() + 3 * n);
    if (!allFinite(gpuTrack_.y)) {  // un float qui déborde : on s'arrête plutôt que de dessiner des valeurs infinies
        gpuTrack_.diverged = true;
        running_ = false;
        return;
    }
    recordConserved(gpuTrack_);
    recordEnergy(gpuTrack_, gpu_.potentialEnergy());
}

void GpuNBodyModule::stepCpu(int steps) {
    const int n = bodyCount();
    const std::size_t m = 3 * static_cast<std::size_t>(n);
    double* x = cpuTrack_.y.data();
    double* v = x + m;
    const double G = problem_.G, eps = problem_.softening;
    const auto t0 = Clock::now();
    if (!cpuAccValid_) {
        cpuAcc_.assign(m, 0.0);
        nbody::accelerations(x, problem_.mass.data(), n, G, eps, cpuAcc_.data());
        cpuAccValid_ = true;
    }
    const double half = 0.5 * dt_;
    for (int s = 0; s < steps; ++s) {  // kick-drift-kick avec l'accélération du dernier point en mémoire (1 évaluation par pas)
        for (std::size_t k = 0; k < m; ++k) v[k] += half * cpuAcc_[k];
        for (std::size_t k = 0; k < m; ++k) x[k] += dt_ * v[k];
        nbody::accelerations(x, problem_.mass.data(), n, G, eps, cpuAcc_.data());
        for (std::size_t k = 0; k < m; ++k) v[k] += half * cpuAcc_[k];
    }
    const double ms = 1e3 * elapsedSeconds(t0, Clock::now()) / steps;
    cpuTrack_.stepMillis = cpuTrack_.stepMillis == 0.0 ? ms : 0.7 * cpuTrack_.stepMillis + 0.3 * ms;
    cpuTrack_.stepMs.add(time_, ms);

    if (!allFinite(cpuTrack_.y)) {
        cpuTrack_.diverged = true;
        running_ = false;
        return;
    }
    recordConserved(cpuTrack_);
    // L'énergie potentielle du CPU est en O(N²) : au plus toutes les 0,3 s réelles quand N est grand.
    const auto now = Clock::now();
    if (n <= 300 || elapsedSeconds(lastCpuEnergy_, now) > 0.3) {
        recordEnergy(cpuTrack_, nbody::potentialEnergy(x, problem_.mass.data(), n, G, eps));
        lastCpuEnergy_ = now;
    }
}

void GpuNBodyModule::update(double frameSeconds) {
    if (!running_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    const int wanted = static_cast<int>(accumulator_ / dt_);
    if (wanted < 1) return;

    // Budget : on calcule autant de pas que le temps par pas mesuré le permet dans une image.
    const double perStep = (useGpu() ? gpuTrack_.stepMillis : 0.0) + (useCpu() ? cpuTrack_.stepMillis : 0.0);
    const int allowed = perStep > 0.0 ? std::max(1, static_cast<int>(kBudgetMs / perStep)) : 4;  // 1re image : sonde
    const int steps = std::min({wanted, allowed, kMaxStepsPerFrame});

    time_ += steps * dt_;
    if (useGpu()) stepGpu(steps);
    if (useCpu()) stepCpu(steps);
    if (useGpu() && useCpu() && !gpuTrack_.diverged && !cpuTrack_.diverged)
        gap_.add(time_, floorForLog(problem_.distance(gpuTrack_.y, cpuTrack_.y)));

    accumulator_ -= steps * dt_;
    if (steps < wanted) accumulator_ = 0.0;  // trop lent : on ralentit le temps plutôt que d'accumuler du retard
    achieved_ = 0.8 * achieved_ + 0.2 * static_cast<double>(steps) / wanted;
}

// --------------------------------- UI ----------------------------------

void GpuNBodyModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    level_ = level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    bool changed = false;
    if (backend_ == kBoth && !etudiant) {  // « les deux » n'existe qu'à partir du niveau 5
        backend_ = kGpu;
        changed = true;
    }

    ImGui::SeparatorText("Scénario");
    int s = scenario_;
    ImGui::RadioButton("Un amas d'étoiles", &s, kCluster);
    ImGui::RadioButton("Deux amas qui se rencontrent", &s, kCollision);
    if (s != scenario_) {
        scenario_ = static_cast<Scenario>(s);
        changed = true;
    }

    ImGui::SeparatorText("Qui calcule ?");
    int b = backend_;
    if (gpuOk_) {
        ImGui::RadioButton(college ? "Carte graphique (GPU, float)" : "La carte graphique", &b, kGpu);
        ImGui::RadioButton(college ? "Processeur (CPU, double)" : "Le processeur", &b, kCpu);
        if (etudiant) ImGui::RadioButton("Les deux, côte à côte", &b, kBoth);
    } else {
        wrapped("La carte graphique n'a pas pu lancer le calcul (compute shader refusé) : seul le processeur est disponible.", true);
        b = kCpu;
    }
    if (b != backend_) {
        backend_ = static_cast<Backend>(b);
        changed = true;
    }

    pushSliderWidth();
    bodies_ = std::min(bodies_, static_cast<double>(maxBodies()));
    changed |= sliderD("Nombre d'étoiles", &bodies_, kMinBodies, maxBodies(), "%.0f", ImGuiSliderFlags_Logarithmic);
    if (interesse) {
        if (ImGui::Button("Nouveau tirage")) {
            seed_ = seed_ * 1664525u + 1013904223u;  // graine suivante (congruentiel linéaire), reproductible
            changed = true;
        }
        sliderD("Vitesse du temps", &timeScale_, 0.1, 20.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    }
    if (college) changed |= sliderD("Pas de calcul dt", &dt_, 5e-4, 2e-2, "%.4f", ImGuiSliderFlags_Logarithmic);
    if (etudiant) changed |= sliderD("Adoucissement ε", &softening_, 0.01, 0.3, "%.3f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();
    if (backend_ != kGpu) ImGui::TextDisabled("Processeur : au plus %d étoiles.", kMaxBodiesCpu);

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ ? "Pause" : "Lecture")) {
        if (gpuTrack_.diverged || cpuTrack_.diverged) changed = true;
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) changed = true;

    if (changed) pendingReset_ = true;
    // Un curseur qu'on tire change à chaque image : on ne recrée l'amas qu'au relâchement.
    if (pendingReset_ && !ImGui::IsAnyItemActive()) {
        if (gpuTrack_.diverged || cpuTrack_.diverged) running_ = true;  // « Recommencer » après une divergence relance le calcul
        reset();
    }

    const int n = bodyCount();
    const double pairs = 0.5 * n * (n - 1.0);
    ImGui::SeparatorText("Calcul");
    if (gpuTrack_.diverged || cpuTrack_.diverged) wrapped("Le calcul a divergé (valeurs infinies) : « Recommencer ».", true);
    wrapped(strf("%d étoiles, %s à chaque pas", n, pairsText(pairs).c_str()));
    if (useGpu()) wrapped(strf("Carte graphique : %.2f ms par pas", gpuTrack_.stepMillis));
    if (useGpu() && college)
        wrapped(strf("   dont %.2f ms de calcul pur (le reste : attente, la carte dessine aussi)", gpuTrack_.kernelMillis), true);
    if (useCpu()) wrapped(strf("Processeur : %.2f ms par pas", cpuTrack_.stepMillis));
    if (useGpu() && useCpu() && gpuTrack_.stepMillis > 0.0)
        wrapped(strf("La carte graphique est %.1f fois plus rapide", cpuTrack_.stepMillis / gpuTrack_.stepMillis));
    if (running_ && achieved_ < 0.95) wrapped(strf("Trop lourd pour le temps réel : le temps est ralenti (x%.2f).", achieved_), true);
    if (college) {
        wrapped(strf("t = %.2f   (%.0f images/s)", time_, static_cast<double>(ImGui::GetIO().Framerate)), true);
        if (useGpu() && gpuTrack_.kernelMillis > 0.0)
            wrapped(strf("Débit du calcul GPU : %.1f milliards d'interactions/s", static_cast<double>(n) * n / (gpuTrack_.kernelMillis * 1e-3) * 1e-9), true);
    }
    if (lycee && !exactEnergy_) wrapped("E0 estimée par le shader (N trop grand pour le calcul exact sur le CPU).", true);
}

void GpuNBodyModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (!atLeast(level, Level::College)) {
        ImGui::TextWrapped("Pendant que les étoiles bougent, l'énergie totale de l'amas ne doit pas changer. Si l'ordinateur calcule "
                           "bien, elle reste constante : c'est ainsi qu'on vérifie qu'il ne se trompe pas.\n\n"
                           "Observez aussi le « temps par pas » dans la fenêtre Simulation : le plus petit gagne.");
        return;
    }

    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const int n = bodyCount();
    ImGui::SeparatorText("Système");
    wrapped(strf("%d étoiles, masse totale %.3g, G = 1, ε = %.3g", n, problem_.totalMass(), problem_.softening));
    wrapped(strf("E0 = %.5f   |P0| = %.1e   |L0| = %.4f", e0_, p0_.norm(), l0_.norm()));

    struct Entry { const char* name; const float* color; const Track* track; };
    std::vector<Entry> entries;
    if (useGpu()) entries.push_back({etudiant ? "GPU float" : "Carte graphique", kGpuColor, &gpuTrack_});
    if (useCpu()) entries.push_back({etudiant ? "CPU double" : "Processeur", kCpuColor, &cpuTrack_});

    std::vector<std::string> headers;
    if (lycee) headers = {"ΔE/E0", "|ΔP|", "|ΔL|"};
    else headers = {"ΔE en %"};
    std::vector<TableRow> rows;
    for (const Entry& e : entries) {
        TableRow r{e.name, e.color, {}};
        if (lycee) r.cells = {strf("%+.1e", e.track->lastEnergy), strf("%.1e", e.track->lastMomentum), strf("%.1e", e.track->lastAngular)};
        else r.cells = {strf("%+.1e", 100.0 * e.track->lastEnergy)};
        rows.push_back(r);
    }
    ImGui::SeparatorText("Ce qui se conserve");
    drawResultTable("conservation", headers, rows);

    std::vector<TableRow> timeRows;
    for (const Entry& e : entries)
        timeRows.push_back({e.name, e.color, {strf("%.3f", e.track->stepMillis), e.track->stepMillis > 0.0 ? strf("%.0f", 1e3 / e.track->stepMillis) : "-"}});
    ImGui::SeparatorText("Vitesse");
    drawResultTable("vitesse", {"ms/pas", "pas/s"}, timeRows);

    if (useGpu() && useCpu() && !gap_.x.empty()) {
        const int last = (gap_.offset + gap_.size() - 1) % std::max(gap_.size(), 1);
        wrapped(strf("Écart GPU − CPU (espace des phases) : %.1e", gap_.y[last]));
    }
    if (lycee)
        wrapped("ΔE/E0 : variation relative de l'énergie. |ΔP|, |ΔL| : variations de l'impulsion et du moment cinétique, relatives à "
                "Σ m|v| et Σ m|r × v|. Le processeur en double garde P et L à 1e-15 ; en float ils bougent un peu.", true);
    else
        wrapped("ΔE : variation de l'énergie totale depuis le départ, en % (plus c'est petit, plus le calcul est juste).", true);
}

void GpuNBodyModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool both = atLeast(level, Level::Etudiant) && useGpu() && useCpu();
    const int cols = 1 + (college ? 1 : 0) + (lycee ? 1 : 0) + (both ? 1 : 0);

    const auto line = [](const char* name, const Series& s, const float* color) {
        if (s.size() > 0) ImPlot::PlotLine(name, s.x.data(), s.y.data(), s.size(), lineSpec(color, s.offset));
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Temps de calcul par pas (ms)")) {
            ImPlot::SetupAxes("t", "ms", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (useGpu() && useCpu()) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);  // des ordres de grandeur d'écart
            if (useGpu()) line(college ? "GPU float" : "carte graphique", gpuTrack_.stepMs, kGpuColor);
            if (useCpu()) line(college ? "CPU double" : "processeur", cpuTrack_.stepMs, kCpuColor);
            ImPlot::EndPlot();
        }
        if (college && ImPlot::BeginPlot("Énergie (E − E0) / |E0|")) {
            ImPlot::SetupAxes("t", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (useGpu()) line("GPU float", gpuTrack_.energy, kGpuColor);
            if (useCpu()) line("CPU double", cpuTrack_.energy, kCpuColor);
            ImPlot::EndPlot();
        }
        if (lycee && ImPlot::BeginPlot("Impulsion et moment cinétique")) {
            ImPlot::SetupAxes("t", "écart relatif", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (useGpu()) {
                line("|ΔP| GPU", gpuTrack_.momentum, kGpuColor);
                const float darker[3] = {0.75f, 0.35f, 0.05f};
                line("|ΔL| GPU", gpuTrack_.angular, darker);
            }
            if (useCpu()) {
                line("|ΔP| CPU", cpuTrack_.momentum, kCpuColor);
                const float darker[3] = {0.20f, 0.40f, 0.75f};
                line("|ΔL| CPU", cpuTrack_.angular, darker);
            }
            ImPlot::EndPlot();
        }
        if (both && ImPlot::BeginPlot("Écart GPU − CPU (espace des phases)")) {
            ImPlot::SetupAxes("t", "distance", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            line("écart", gap_, kGrey);
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

// Compare les accélérations du shader à celles du CPU double sur l'état courant (le CPU vérifie le GPU).
void GpuNBodyModule::runValidation() {
    lastValidation_ = Clock::now();
    const int n = bodyCount();
    if (!useGpu() || n > kMaxBodiesCpu || gpuTrack_.diverged) return;

    std::vector<double> accGpu(3 * static_cast<std::size_t>(n)), accRef(accGpu.size());
    gpu_.readAccelerations(accGpu.data());
    const auto t0 = Clock::now();
    nbody::accelerations(gpuTrack_.y.data(), problem_.mass.data(), n, problem_.G, problem_.softening, accRef.data());
    Validation v;
    v.cpuMillis = 1e3 * elapsedSeconds(t0, Clock::now());
    v.potentialCpu = problem_.potentialEnergy(gpuTrack_.y);
    v.potentialGpu = gpu_.potentialEnergy();

    std::vector<double> rel;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        const Vec3 a{accGpu[3 * i], accGpu[3 * i + 1], accGpu[3 * i + 2]}, r{accRef[3 * i], accRef[3 * i + 1], accRef[3 * i + 2]};
        const double d = (a - r).norm(), rn = r.norm();
        num += d * d;
        den += rn * rn;
        if (rn > 1e-300) rel.push_back(d / rn);
    }
    if (rel.empty()) return;
    std::sort(rel.begin(), rel.end());
    v.valid = true;
    v.count = n;
    v.time = time_;
    v.globalRel = den > 0.0 ? std::sqrt(num / den) : 0.0;
    v.maxRel = rel.back();
    v.medianRel = rel[rel.size() / 2];
    v.p99Rel = rel[std::min(rel.size() - 1, static_cast<std::size_t>(0.99 * static_cast<double>(rel.size())))];

    constexpr int kBins = 17;  // log10(erreur relative) de -9 à -1 par pas de 0,5
    v.centers.resize(kBins);
    v.counts.assign(kBins, 0.0);
    for (int k = 0; k < kBins; ++k) v.centers[k] = -9.0 + 0.5 * k;
    for (double e : rel) {
        const int bin = std::clamp(static_cast<int>(std::floor((std::log10(std::max(e, 1e-12)) + 9.25) / 0.5)), 0, kBins - 1);
        v.counts[bin] += 1.0;
    }
    validation_ = std::move(v);
}

// Débit et erreur de l'accélération du shader float en fonction de N, sur des sphères de Plummer (quelques secondes).
void GpuNBodyModule::runBenchmark() {
    bench_.clear();
    GpuNBody g;
    if (!g.init(shaderDir_)) return;
    for (int n : {250, 500, 1000, 2000, 4000, 8000, 16000}) {
        const NBodyProblem p = NBodyProblem::plummer(n, 1, 1.0, 0.05);
        std::vector<double> pos(3 * static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            pos[3 * i] = p.position[i].x;
            pos[3 * i + 1] = p.position[i].y;
            pos[3 * i + 2] = p.position[i].z;
        }
        g.setBodies(pos.data(), p.mass.data(), n);
        g.compute(1.0, 0.05);  // échauffement
        double best = 1e30;
        for (int rep = 0; rep < 3; ++rep) {
            g.compute(1.0, 0.05);
            best = std::min(best, g.lastKernelSeconds());
        }
        BenchRow row;
        row.n = n;
        row.gpuGints = static_cast<double>(n) * n / best * 1e-9;
        if (n <= kMaxBodiesCpu) {
            std::vector<double> accCpu(pos.size()), accGpu(pos.size());
            const auto t0 = Clock::now();
            nbody::accelerations(pos.data(), p.mass.data(), n, 1.0, 0.05, accCpu.data());
            row.cpuGints = static_cast<double>(n) * n / elapsedSeconds(t0, Clock::now()) * 1e-9;
            g.readAccelerations(accGpu.data());
            double num = 0.0, den = 0.0;
            for (std::size_t k = 0; k < accCpu.size(); ++k) {
                num += (accGpu[k] - accCpu[k]) * (accGpu[k] - accCpu[k]);
                den += accCpu[k] * accCpu[k];
            }
            row.error = den > 0.0 ? std::sqrt(num / den) : 0.0;
        }
        bench_.push_back(row);
    }
    g.shutdown();
}

void GpuNBodyModule::drawAnalysis(const UiContext& ctx) {
    const bool chercheur = atLeast(ctx.level, Level::Chercheur);
    const int n = bodyCount();
    const bool canValidate = useGpu() && n <= kMaxBodiesCpu && !gpuTrack_.diverged;

    // Une première mesure après chaque recommencement, puis toutes les 2 s tant que le CPU la fait en quelques millisecondes ;
    // au-delà de 1500 étoiles, sur demande.
    if (canValidate && !validation_.valid) runValidation();
    else if (canValidate && n <= kAutoValidationMaxN && running_ && elapsedSeconds(lastValidation_, Clock::now()) > 2.0) runValidation();
    if (canValidate) {
        if (ImGui::Button("Comparer maintenant GPU et CPU")) runValidation();
    } else {
        wrapped(useGpu() ? "Comparaison au CPU : au plus 4000 étoiles." : "Comparaison disponible avec la carte graphique.", true);
    }
    if (chercheur) {
        ImGui::SameLine();
        if (ImGui::Button("Banc d'essai (quelques secondes)")) runBenchmark();
    }
    if (validation_.valid) {
        wrapped(strf("État t = %.2f, %d étoiles. Erreur relative de l'accélération : maximale %.1e, médiane %.1e, 99%% des étoiles sous %.1e ; "
                     "globale %.1e. Énergie potentielle GPU/CPU : écart %.1e. (CPU double : %.1f ms)",
                     validation_.time, validation_.count, validation_.maxRel, validation_.medianRel, validation_.p99Rel,
                     validation_.globalRel, std::abs(validation_.potentialGpu / validation_.potentialCpu - 1.0), validation_.cpuMillis));
    }

    const int cols = 1 + (chercheur && !bench_.empty() ? 2 : 0);
    if (ImPlot::BeginSubplots("##analyse", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Erreur de l'accélération GPU float, étoile par étoile")) {
            ImPlot::SetupAxes("log10(erreur relative)", "nombre d'étoiles", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (validation_.valid)
                ImPlot::PlotBars("étoiles", validation_.centers.data(), validation_.counts.data(), static_cast<int>(validation_.centers.size()), 0.45,
                                 ImPlotSpec(ImPlotProp_FillColor, toImVec4(kGpuColor, 0.8f)));
            ImPlot::EndPlot();
        }
        if (chercheur && !bench_.empty()) {
            const auto markers = [](const float* c) {
                return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle,
                                  ImPlotProp_MarkerSize, 4.0f);
            };
            std::vector<double> ns, gpu, nsCpu, cpu, nsErr, err, ref;
            for (const BenchRow& r : bench_) {
                ns.push_back(r.n);
                gpu.push_back(r.gpuGints);
                if (r.cpuGints > 0.0) { nsCpu.push_back(r.n); cpu.push_back(r.cpuGints); }
                if (r.error > 0.0) { nsErr.push_back(r.n); err.push_back(r.error); ref.push_back(6e-8 * std::sqrt(static_cast<double>(r.n))); }
            }
            if (ImPlot::BeginPlot("Débit (milliards d'interactions/s)")) {
                ImPlot::SetupAxes("N", "Ginteractions/s", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlot::PlotLine("GPU float", ns.data(), gpu.data(), static_cast<int>(ns.size()), markers(kGpuColor));
                if (!nsCpu.empty()) ImPlot::PlotLine("CPU double", nsCpu.data(), cpu.data(), static_cast<int>(nsCpu.size()), markers(kCpuColor));
                ImPlot::EndPlot();
            }
            if (ImPlot::BeginPlot("Erreur globale de l'accélération (GPU float / CPU double)")) {
                ImPlot::SetupAxes("N", "écart relatif", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                if (!nsErr.empty()) {
                    ImPlot::PlotLine("mesure", nsErr.data(), err.data(), static_cast<int>(nsErr.size()), markers(kGpuColor));
                    ImPlot::PlotLine("6e-8 √N", nsErr.data(), ref.data(), static_cast<int>(nsErr.size()), lineSpec(kGrey));
                }
                ImPlot::EndPlot();
            }
        }
        ImPlot::EndSubplots();
    }
}

// ------------------------------ rendu 3D --------------------------------

// Position de l'étoile i -> scène (le plan du mouvement devient l'horizontale), couleur selon la vitesse si `color` est nul.
Vertex GpuNBodyModule::starVertex(const State& y, int i, const float* color, float speedScale) const {
    const int n = bodyCount();
    float c[3];
    if (color) {
        c[0] = color[0]; c[1] = color[1]; c[2] = color[2];
    } else {
        const double vx = y[3 * n + 3 * i], vy = y[3 * n + 3 * i + 1], vz = y[3 * n + 3 * i + 2];
        const float t = std::clamp(static_cast<float>(std::sqrt(vx * vx + vy * vy + vz * vz)) / speedScale, 0.0f, 1.0f);
        c[0] = 0.35f + 0.65f * t;   // lent : bleu froid ; rapide : jaune chaud
        c[1] = 0.55f + 0.30f * t;
        c[2] = 1.00f - 0.70f * t;
    }
    return {static_cast<float>(y[3 * i]) * kScale, static_cast<float>(y[3 * i + 2]) * kScale,
            -static_cast<float>(y[3 * i + 1]) * kScale, c[0], c[1], c[2]};
}

void GpuNBodyModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer.draw(Primitive::Lines, grid);
    renderer.draw(Primitive::Lines, axes);

    const int n = bodyCount();
    const float size = (n <= 500 ? 5.0f : n <= 3000 ? 3.5f : 2.5f) * ctx.uiScale;
    std::vector<Vertex> stars;
    stars.reserve(static_cast<std::size_t>(n));
    const float speedScale = static_cast<float>(speedScale_);

    if (useGpu() && useCpu()) {  // côte à côte : le CPU en halo cyan dessous, le GPU en orange dessus
        const float halo[3] = {0.30f, 0.90f, 0.95f};
        for (int i = 0; i < n; ++i) stars.push_back(starVertex(cpuTrack_.y, i, halo, speedScale));
        renderer.draw(Primitive::Points, stars, size * 2.0f);
        stars.clear();
        for (int i = 0; i < n; ++i) stars.push_back(starVertex(gpuTrack_.y, i, kGpuColor, speedScale));
        renderer.draw(Primitive::Points, stars, size);
        return;
    }
    const Track& t = useGpu() ? gpuTrack_ : cpuTrack_;
    if (t.diverged) return;
    for (int i = 0; i < n; ++i) stars.push_back(starVertex(t.y, i, nullptr, speedScale));
    renderer.draw(Primitive::Points, stars, size);
}

}  // namespace pl
