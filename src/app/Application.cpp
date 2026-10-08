#include "Application.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <glad/gl.h>
// glad doit précéder GLFW
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_internal.h>  // DockBuilder (disposition par défaut)
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/core/Level.hpp"
#include "physicslab/core/Solver.hpp"
#include "physicslab/core/World.hpp"
#include "physicslab/mechanics/Projectile.hpp"
#include "physicslab/render/Camera.hpp"
#include "physicslab/render/Renderer.hpp"

#ifndef PHYSICSLAB_SHADER_DIR
#define PHYSICSLAB_SHADER_DIR "shaders"
#endif

namespace pl {
namespace {

// ---------------------------------------------------------------------------
// M1 : projectile avec frottement. Un même lancer est intégré en parallèle par
// 5 solveurs et comparé en direct à la solution analytique.
// ---------------------------------------------------------------------------

enum RunIndex { kEuler = 0, kSymplectic, kVerlet, kRK4, kRK45, kRunCount };
constexpr int kFixedStepSolvers = 4;  // Euler, symplectique, Verlet, RK4 (l'étude de convergence)

const float kBlue[3] = {0.35f, 0.65f, 1.0f};  // couleur de la solution exacte

// Courbe temps réel dans un tampon circulaire (compatible avec l'argument `offset` d'ImPlot).
struct Series {
    static constexpr int kCapacity = 4000;
    std::vector<double> x, y;
    int offset = 0;

    int size() const { return static_cast<int>(x.size()); }
    void clear() { x.clear(); y.clear(); offset = 0; }
    void add(double xx, double yy) {
        if (size() < kCapacity) {
            x.push_back(xx);
            y.push_back(yy);
        } else {
            x[offset] = xx;
            y[offset] = yy;
            offset = (offset + 1) % kCapacity;
        }
    }
};

// Un solveur, son monde et ses courbes.
struct Run {
    std::unique_ptr<Solver> solver;
    World world;
    Series y, posError, energyError;
    std::vector<Vertex> trail;
    float color[3] = {1, 1, 1};
    bool show = true;
};

// Erreur finale en fonction du pas, pour un solveur à pas fixe.
struct ConvergenceCurve {
    std::vector<double> dt, err;
    double slope = 0.0;  // ordre mesuré (pente en log-log)
};

const char* explanation(Level level) {
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
                   "Le panneau Convergence trace l'erreur finale en fonction de dt : en log-log, la pente est l'ordre.";
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

bool sliderD(const char* label, double* v, double lo, double hi, const char* fmt, ImGuiSliderFlags flags = 0) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &lo, &hi, fmt, flags);
}

ImVec4 toImVec4(const float c[3], float alpha = 1.0f) { return ImVec4(c[0], c[1], c[2], alpha); }

std::unique_ptr<Solver> makeFixedStepSolver(int index) {
    switch (index) {
        case kEuler: return std::make_unique<ExplicitEuler>();
        case kSymplectic: return std::make_unique<SymplecticEuler>();
        case kVerlet: return std::make_unique<VelocityVerlet>();
        default: return std::make_unique<RK4>();
    }
}

void applyTheme(float scale) {
    ImGui::StyleColorsDark();
    ImGuiStyle& s = ImGui::GetStyle();
    s.WindowRounding = 6.0f;
    s.ChildRounding = 4.0f;
    s.FrameRounding = 4.0f;
    s.PopupRounding = 4.0f;
    s.GrabRounding = 4.0f;
    s.TabRounding = 4.0f;
    s.WindowBorderSize = 0.0f;
    s.FrameBorderSize = 0.0f;
    s.WindowPadding = ImVec2(12.0f, 10.0f);
    s.FramePadding = ImVec2(8.0f, 5.0f);
    s.ItemSpacing = ImVec2(8.0f, 7.0f);

    const ImVec4 bg(0.105f, 0.115f, 0.135f, 0.96f);
    const ImVec4 panel(0.15f, 0.165f, 0.19f, 1.0f);
    const ImVec4 accent(0.20f, 0.60f, 0.78f, 1.0f);
    const ImVec4 accentHi(0.27f, 0.70f, 0.88f, 1.0f);
    ImVec4* c = s.Colors;
    c[ImGuiCol_WindowBg] = bg;
    c[ImGuiCol_PopupBg] = ImVec4(0.12f, 0.13f, 0.15f, 0.98f);
    c[ImGuiCol_MenuBarBg] = ImVec4(0.09f, 0.10f, 0.12f, 1.0f);
    c[ImGuiCol_TitleBg] = ImVec4(0.09f, 0.10f, 0.12f, 1.0f);
    c[ImGuiCol_TitleBgActive] = ImVec4(0.11f, 0.12f, 0.14f, 1.0f);
    c[ImGuiCol_FrameBg] = panel;
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.19f, 0.21f, 0.25f, 1.0f);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.22f, 0.24f, 0.29f, 1.0f);
    c[ImGuiCol_SliderGrab] = accent;
    c[ImGuiCol_SliderGrabActive] = accentHi;
    c[ImGuiCol_CheckMark] = accentHi;
    c[ImGuiCol_Button] = ImVec4(0.17f, 0.36f, 0.46f, 1.0f);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.20f, 0.46f, 0.58f, 1.0f);
    c[ImGuiCol_ButtonActive] = accent;
    c[ImGuiCol_Header] = ImVec4(0.17f, 0.36f, 0.46f, 0.55f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.20f, 0.46f, 0.58f, 0.80f);
    c[ImGuiCol_HeaderActive] = accent;
    c[ImGuiCol_Tab] = panel;
    c[ImGuiCol_TabHovered] = accent;
    c[ImGuiCol_TabSelected] = ImVec4(0.17f, 0.36f, 0.46f, 1.0f);
    c[ImGuiCol_DockingPreview] = ImVec4(0.20f, 0.60f, 0.78f, 0.45f);
    s.ScaleAllSizes(scale);
}

void loadFont(ImGuiIO& io, float size) {
    // Police système lisible (accents, symboles) ; à défaut, la police intégrée d'ImGui.
    const char* candidates[] = {
        "C:/Windows/Fonts/segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
    };
    for (const char* path : candidates) {
        if (std::filesystem::exists(path) && io.Fonts->AddFontFromFileTTF(path, size)) return;
    }
}

class App {
public:
    int run(const AppOptions& options);

private:
    bool initWindow(bool visible);
    void initImGui();
    void initRuns();
    void shutdown();

    // simulation
    void resetSimulation();
    void stepSimulation(double frameSeconds);
    void sample();
    void computeConvergence();

    // niveau pédagogique
    bool isShown(int run) const;
    std::string label(int run) const;

    // interface
    void buildUi();
    void buildDefaultLayout(ImGuiID dockId, const ImGuiViewport* viewport);
    void drawMenuBar();
    void drawSimulationPanel();
    void drawExplanationPanel();
    void drawInvariantsPanel();
    void drawGraphsPanel();
    void drawConvergencePanel();
    void drawScene(int fbW, int fbH);

    // fenêtre et rendu
    GLFWwindow* window_ = nullptr;
    Renderer renderer_;
    Camera camera_;
    float uiScale_ = 1.0f;
    ViewRect viewRect_{0, 0, 1, 1};
    bool resetLayout_ = false;

    // pédagogie
    Level level_ = Level::College;

    // paramètres du lancer (unités SI) et du calcul
    double speed_ = 14.0;               // [m/s]
    double angleDeg_ = 55.0;            // [deg]
    double height_ = 2.0;               // [m]
    double gravity_ = constants::g0;    // [m/s^2]
    double mass_ = 1.0;                 // [kg]
    double drag_ = 0.0;                 // b [kg/s]
    double dt_ = 1.0 / 30.0;            // pas de calcul [s] (grossier à dessein : les erreurs se voient)
    double relTol_ = 1e-8;              // tolérance relative de RK45
    double timeScale_ = 1.0;

    // simulation
    ProjectileProblem problem_;
    std::array<Run, kRunCount> runs_;
    RK45* rk45_ = nullptr;              // même objet que runs_[kRK45].solver
    Series refY_;                       // hauteur exacte
    std::vector<Vertex> refPath_;       // trajectoire exacte complète
    std::array<ConvergenceCurve, kFixedStepSolvers> convergence_;
    double convergenceEnd_ = 1.0;       // instant final de l'étude de convergence
    bool convergenceForcedDrag_ = false;  // l'étude impose un frottement (b nul => solution polynomiale)
    double landingTime_ = 0.0;
    double accumulator_ = 0.0;
    double lastSampleTime_ = 0.0;
    bool running_ = true;
    bool landed_ = false;
};

// --------------------------- fenêtre / ImGui ---------------------------

bool App::initWindow(bool visible) {
    if (!glfwInit()) {
        std::fprintf(stderr, "GLFW : initialisation impossible\n");
        return false;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, visible ? GLFW_TRUE : GLFW_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    // Fenêtre à 85 % de la zone de travail du moniteur principal, centrée (1600x900 si inconnue).
    int width = 1600, height = 900;
    int workX = 0, workY = 0, workW = 0, workH = 0;
    glfwGetMonitorWorkarea(glfwGetPrimaryMonitor(), &workX, &workY, &workW, &workH);
    if (workW > 0 && workH > 0) {
        width = workW * 85 / 100;
        height = workH * 85 / 100;
    }
    window_ = glfwCreateWindow(width, height, "PhysicsLab 3D", nullptr, nullptr);
    if (!window_) {
        std::fprintf(stderr, "GLFW : OpenGL 4.5 indisponible (pilote graphique trop ancien ?)\n");
        return false;
    }
    if (workW > 0 && workH > 0) glfwSetWindowPos(window_, workX + (workW - width) / 2, workY + (workH - height) / 2);
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1);

    if (!gladLoadGL(glfwGetProcAddress)) {
        std::fprintf(stderr, "GLAD : chargement d'OpenGL impossible\n");
        return false;
    }
    std::printf("OpenGL %s | %s\n", reinterpret_cast<const char*>(glGetString(GL_VERSION)),
                reinterpret_cast<const char*>(glGetString(GL_RENDERER)));
    glEnable(GL_MULTISAMPLE);

    return renderer_.init(PHYSICSLAB_SHADER_DIR);
}

void App::initImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable | ImGuiConfigFlags_NavEnableKeyboard;
    resetLayout_ = !std::filesystem::exists("imgui.ini");  // 1er lancement : disposition par défaut

    float xs = 1.0f, ys = 1.0f;
    glfwGetWindowContentScale(window_, &xs, &ys);
    uiScale_ = std::max(1.0f, xs);
    loadFont(io, 17.0f * uiScale_);
    applyTheme(uiScale_);

    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL3_Init("#version 450");
}

void App::shutdown() {
    if (ImGui::GetCurrentContext()) {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
    }
    renderer_.shutdown();
    if (window_) glfwDestroyWindow(window_);
    glfwTerminate();
}

void App::initRuns() {
    const float colors[kRunCount][3] = {{1.00f, 0.55f, 0.15f},   // Euler : orange
                                        {0.80f, 0.85f, 0.25f},   // symplectique : citron
                                        {0.85f, 0.45f, 0.95f},   // Verlet : violet
                                        {0.35f, 0.85f, 0.45f},   // RK4 : vert
                                        {1.00f, 0.40f, 0.50f}};  // RK45 : rose
    for (int i = 0; i < kRunCount; ++i) {
        if (i == kRK45) {
            auto rk = std::make_unique<RK45>();
            rk45_ = rk.get();
            runs_[i].solver = std::move(rk);
        } else {
            runs_[i].solver = makeFixedStepSolver(i);
        }
        std::copy(colors[i], colors[i] + 3, runs_[i].color);
    }
}

// ------------------------------ simulation -----------------------------

void App::resetSimulation() {
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
    rk45_->resetStats();

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

void App::sample() {
    const double t = runs_[0].world.time;
    lastSampleTime_ = t;
    const Vec3 ref = problem_.position(t);
    const double exactEnergy = problem_.energy(t);

    refY_.add(t, ref.y);
    for (Run& r : runs_) {
        const Particle& p = r.world.particles[0];
        r.y.add(t, p.position.y);
        r.posError.add(t, std::max((p.position - ref).norm(), 1e-12));  // plancher : l'axe est logarithmique
        r.energyError.add(t, r.world.invariants().total() - exactEnergy);
        r.trail.push_back({static_cast<float>(p.position.x), static_cast<float>(p.position.y),
                           static_cast<float>(p.position.z), r.color[0], r.color[1], r.color[2]});
    }
}

void App::stepSimulation(double frameSeconds) {
    if (!running_ || landed_) return;

    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;  // borne anti "spirale de la mort"
    rk45_->relTol = relTol_;
    rk45_->absTol = relTol_ * 1e-2;
    const double sampleInterval = std::max(dt_, 1.0 / 240.0);

    int guard = 0;
    while (!landed_ && guard++ < 5000) {
        const double remaining = landingTime_ - runs_[0].world.time;
        const double h = std::min(dt_, remaining);  // dernier pas raccourci : tous arrivent pile à l'impact exact
        if (accumulator_ < h) break;

        for (Run& r : runs_) r.world.step(*r.solver, h);
        accumulator_ -= h;
        landed_ = remaining <= dt_ * (1.0 + 1e-9);
        if (landed_ || runs_[0].world.time - lastSampleTime_ >= sampleInterval - 1e-12) sample();
    }
    if (landed_) running_ = false;
}

// Erreur finale de chaque solveur à pas fixe pour 7 pas différents, puis pente log-log (ordre mesuré).
void App::computeConvergence() {
    static const int kStepCounts[] = {5, 10, 20, 40, 80, 160, 320};
    convergenceEnd_ = std::max(0.2, landingTime_);

    // Sans frottement la solution est un polynôme de degré 2 que Verlet et RK4 intègrent exactement :
    // il n'y aurait rien à mesurer. On impose alors un frottement de référence (k = 0,8 /s).
    ProjectileProblem study = problem_;
    convergenceForcedDrag_ = study.linearDrag < 0.3 * study.mass;
    if (convergenceForcedDrag_) study.linearDrag = 0.8 * study.mass;

    for (int i = 0; i < kFixedStepSolvers; ++i) {
        ConvergenceCurve& c = convergence_[i];
        c.dt.clear();
        c.err.clear();
        const std::unique_ptr<Solver> solver = makeFixedStepSolver(i);
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

// ------------------------------ niveaux --------------------------------

bool App::isShown(int run) const {
    if (!atLeast(level_, Level::College)) return run == kRK4;
    if (!atLeast(level_, Level::Etudiant)) return run == kEuler || run == kRK4;
    return runs_[run].show;
}

std::string App::label(int run) const {
    if (!atLeast(level_, Level::College) && run == kRK4) return "Ordinateur";
    if (!atLeast(level_, Level::Etudiant)) {
        if (run == kEuler) return "Méthode simple (Euler)";
        if (run == kRK4) return "Méthode précise (RK4)";
    }
    return runs_[run].solver->name();
}

// --------------------------------- UI ----------------------------------

void App::drawMenuBar() {
    if (!ImGui::BeginMainMenuBar()) return;

    if (ImGui::BeginMenu("Fichier")) {
        if (ImGui::MenuItem("Quitter")) glfwSetWindowShouldClose(window_, GLFW_TRUE);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Fenêtre")) {
        if (ImGui::MenuItem("Réinitialiser la disposition")) resetLayout_ = true;
        ImGui::EndMenu();
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Niveau");
    ImGui::SetNextItemWidth(250.0f * uiScale_);
    const int idx = static_cast<int>(level_);
    if (ImGui::BeginCombo("##niveau", levelName(level_))) {
        for (int i = 0; i < kLevelCount; ++i) {
            const Level l = static_cast<Level>(i);
            if (ImGui::Selectable(levelName(l), i == idx)) level_ = l;
        }
        ImGui::EndCombo();
    }
    ImGui::EndMainMenuBar();
}

void App::buildDefaultLayout(ImGuiID dockId, const ImGuiViewport* viewport) {
    ImGui::DockBuilderRemoveNode(dockId);
    ImGui::DockBuilderAddNode(dockId, ImGuiDockNodeFlags_DockSpace | ImGuiDockNodeFlags_PassthruCentralNode);
    ImGui::DockBuilderSetNodeSize(dockId, viewport->Size);

    ImGuiID rest = dockId;
    const ImGuiID left = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Left, 0.27f, nullptr, &rest);
    const ImGuiID bottom = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Down, 0.34f, nullptr, &rest);
    const ImGuiID right = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Right, 0.34f, nullptr, &rest);

    ImGui::DockBuilderDockWindow("Simulation", left);
    ImGui::DockBuilderDockWindow("Explication", left);
    ImGui::DockBuilderDockWindow("Graphes", bottom);
    ImGui::DockBuilderDockWindow("Convergence", bottom);
    ImGui::DockBuilderDockWindow("Invariants", right);
    ImGui::DockBuilderFinish(dockId);
}

void App::drawSimulationPanel() {
    ImGui::Begin("Simulation");

    const bool college = atLeast(level_, Level::College);
    bool changed = false;

    // Colonne de libellés calée sur le plus long : curseurs alignés, textes jamais coupés.
    ImGui::PushItemWidth(-(ImGui::CalcTextSize("Tolérance relative RK45").x + ImGui::GetStyle().ItemInnerSpacing.x));

    ImGui::SeparatorText("Lancer");
    changed |= sliderD(college ? "Vitesse (m/s)" : "Force du lancer", &speed_, 1.0, 30.0, college ? "%.1f" : "");
    if (atLeast(level_, Level::Interesse)) changed |= sliderD("Angle (deg)", &angleDeg_, 5.0, 85.0, "%.0f");

    if (level_ == Level::Interesse) {
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
    if (atLeast(level_, Level::Lycee)) {
        changed |= sliderD("Hauteur initiale (m)", &height_, 0.0, 20.0, "%.1f");
        const double a = angleDeg_ * constants::pi / 180.0;
        ImGui::TextDisabled("v0 = (%.2f ; %.2f) m/s", speed_ * std::cos(a), speed_ * std::sin(a));
    }
    changed |= sliderD(college ? "Frottement b (kg/s)" : "Résistance de l'air", &drag_, 0.0, 3.0, college ? "%.2f" : "");

    if (college) {
        ImGui::SeparatorText("Calcul");
        changed |= sliderD("Pas de calcul dt (s)", &dt_, 1e-3, 0.1, "%.3f", ImGuiSliderFlags_Logarithmic);
    }
    if (atLeast(level_, Level::Chercheur))
        changed |= sliderD("Tolérance relative RK45", &relTol_, 1e-12, 1e-3, "%.0e", ImGuiSliderFlags_Logarithmic);

    if (atLeast(level_, Level::Etudiant)) {
        ImGui::TextDisabled("Méthodes affichées");
        for (int i = 0; i < kRunCount; ++i) {
            ImGui::PushStyleColor(ImGuiCol_CheckMark, toImVec4(runs_[i].color));
            ImGui::Checkbox(label(i).c_str(), &runs_[i].show);
            ImGui::PopStyleColor();
        }
    }
    ImGui::PopItemWidth();
    if (changed) resetSimulation();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !landed_ ? "Pause" : "Lecture")) {
        if (landed_) resetSimulation();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) resetSimulation();
    ImGui::PushItemWidth(-(ImGui::CalcTextSize("Tolérance relative RK45").x + ImGui::GetStyle().ItemInnerSpacing.x));
    if (atLeast(level_, Level::Interesse))
        sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    ImGui::PopItemWidth();

    if (college) {
        const Run& ref = runs_[kRK4];
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
    ImGui::End();
}

void App::drawExplanationPanel() {
    ImGui::Begin("Explication");
    ImGui::TextColored(ImVec4(0.45f, 0.80f, 0.95f, 1.0f), "%s", levelName(level_));
    ImGui::Separator();
    ImGui::TextWrapped("%s", explanation(level_));
    ImGui::End();
}

void App::drawInvariantsPanel() {
    ImGui::Begin("Invariants");

    if (!atLeast(level_, Level::College)) {
        ImGui::TextWrapped("Pendant le vol, l'énergie du mouvement et l'énergie de la hauteur s'échangent. "
                           "Sans air, leur total reste le même.");
        ImGui::End();
        return;
    }

    const double t = runs_[0].world.time;
    const Vec3 exactPos = problem_.position(t);
    const double exactEnergy = problem_.energy(t);
    const double e0 = problem_.energy(0.0);
    const bool lycee = atLeast(level_, Level::Lycee);
    const bool researcher = atLeast(level_, Level::Chercheur);

    // Chercheur : dE et dE/E0 remplacent Em (redondant) pour que le tableau tienne dans le panneau.
    const bool showEm = !researcher, showDE = lycee, showRel = researcher;
    const int columns = 2 + (showEm ? 1 : 0) + (showDE ? 1 : 0) + (showRel ? 1 : 0);
    if (ImGui::BeginTable("erreurs", columns, ImGuiTableFlags_RowBg)) {
        // Colonnes numériques à largeur fixe (le plus long nombre affiché) ; "Méthode" prend le reste.
        const float numW = ImGui::CalcTextSize("+0.0e+00").x + 2.0f * ImGui::GetStyle().CellPadding.x;
        ImGui::TableSetupColumn("Méthode", ImGuiTableColumnFlags_WidthStretch);
        if (showEm) ImGui::TableSetupColumn("Em (J)", ImGuiTableColumnFlags_WidthFixed, numW);
        if (showDE) ImGui::TableSetupColumn("dE (J)", ImGuiTableColumnFlags_WidthFixed, numW);
        if (showRel) ImGui::TableSetupColumn("dE/E0", ImGuiTableColumnFlags_WidthFixed, numW);
        ImGui::TableSetupColumn("dr (m)", ImGuiTableColumnFlags_WidthFixed, numW);
        ImGui::TableHeadersRow();

        auto dash = [] {
            ImGui::TableNextColumn();
            ImGui::TextDisabled("-");
        };

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::TextColored(toImVec4(kBlue), "Exacte");
        if (showEm) {
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", exactEnergy);
        }
        if (showDE) dash();
        if (showRel) dash();
        dash();

        for (int i = 0; i < kRunCount; ++i) {
            if (!isShown(i)) continue;
            const Run& r = runs_[i];
            const double energy = r.world.invariants().total();
            const double posErr = (r.world.particles[0].position - exactPos).norm();

            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextColored(toImVec4(r.color), "%s", label(i).c_str());
            if (showEm) {
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", energy);
            }
            if (showDE) {
                ImGui::TableNextColumn();
                ImGui::Text("%+.1e", energy - exactEnergy);
            }
            if (showRel) {
                ImGui::TableNextColumn();
                ImGui::Text("%+.1e", (energy - exactEnergy) / e0);
            }
            ImGui::TableNextColumn();
            ImGui::Text("%.1e", posErr);
        }
        ImGui::EndTable();
    }

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
        ImGui::Text("%d pas acceptés, %d rejetés", rk45_->acceptedSteps(), rk45_->rejectedSteps());
        ImGui::Text("%d évaluations de f", rk45_->evaluations());
        if (drag_ == 0.0) {
            const double predicted = 0.5 * mass_ * gravity_ * gravity_ * t * dt_;
            ImGui::TextDisabled("Euler : dE théorique = %+.3e J", predicted);
            ImGui::TextDisabled("(1/2 m g² t dt)");
        }
    }
    ImGui::End();
}

void App::drawGraphsPanel() {
    ImGui::Begin("Graphes");

    const bool showError = atLeast(level_, Level::College);
    const bool showEnergy = atLeast(level_, Level::Lycee);
    const int cols = 1 + (showError ? 1 : 0) + (showEnergy ? 1 : 0);

    auto spec = [](const float c[3], int offset) {
        return ImPlotSpec(ImPlotProp_LineColor, toImVec4(c), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Offset, offset);
    };

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Hauteur y(t)")) {
            ImPlot::SetupAxes("t (s)", "y (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            if (refY_.size() > 0) ImPlot::PlotLine("Exacte", refY_.x.data(), refY_.y.data(), refY_.size(), spec(kBlue, refY_.offset));
            for (int i = 0; i < kRunCount; ++i) {
                if (!isShown(i) || runs_[i].y.size() == 0) continue;
                const Series& s = runs_[i].y;
                ImPlot::PlotLine(label(i).c_str(), s.x.data(), s.y.data(), s.size(), spec(runs_[i].color, s.offset));
            }
            ImPlot::EndPlot();
        }
        if (showError && ImPlot::BeginPlot("Erreur de position |r - r_exact|")) {
            ImPlot::SetupAxes("t (s)", "m", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            for (int i = 0; i < kRunCount; ++i) {
                if (!isShown(i) || runs_[i].posError.size() == 0) continue;
                const Series& s = runs_[i].posError;
                ImPlot::PlotLine(label(i).c_str(), s.x.data(), s.y.data(), s.size(), spec(runs_[i].color, s.offset));
            }
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Erreur d'énergie Em - Em_exacte")) {
            ImPlot::SetupAxes("t (s)", "J", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            for (int i = 0; i < kRunCount; ++i) {
                if (!isShown(i) || runs_[i].energyError.size() == 0) continue;
                const Series& s = runs_[i].energyError;
                ImPlot::PlotLine(label(i).c_str(), s.x.data(), s.y.data(), s.size(), spec(runs_[i].color, s.offset));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
    ImGui::End();
}

void App::drawConvergencePanel() {
    if (!atLeast(level_, Level::Etudiant)) return;
    ImGui::Begin("Convergence");

    ImGui::TextWrapped("Erreur de position à t = %.2f s en fonction du pas dt. En échelle log-log, la pente est "
                       "l'ordre de la méthode (RK45, adaptatif, n'y figure pas).",
                       convergenceEnd_);
    if (convergenceForcedDrag_)
        ImGui::TextDisabled("Frottement k = 0,8 /s imposé : sans frottement, Verlet et RK4 sont exacts (polynôme).");
    if (ImPlot::BeginPlot("##convergence", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("dt (s)", "erreur (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        for (int i = 0; i < kFixedStepSolvers; ++i) {
            const ConvergenceCurve& c = convergence_[i];
            if (!runs_[i].show || c.dt.size() < 2) continue;
            char name[96];
            std::snprintf(name, sizeof(name), "%s (pente %.2f)", runs_[i].solver->name(), c.slope);
            const ImPlotSpec spec(ImPlotProp_LineColor, toImVec4(runs_[i].color), ImPlotProp_LineWeight, 2.0f,
                                  ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f);
            ImPlot::PlotLine(name, c.dt.data(), c.err.data(), static_cast<int>(c.dt.size()), spec);
        }
        ImPlot::EndPlot();
    }
    ImGui::End();
}

void App::buildUi() {
    drawMenuBar();

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const ImGuiID dockId = ImGui::DockSpaceOverViewport(0, vp, ImGuiDockNodeFlags_PassthruCentralNode);
    if (resetLayout_) {
        buildDefaultLayout(dockId, vp);
        resetLayout_ = false;
    }

    drawSimulationPanel();
    drawExplanationPanel();
    drawInvariantsPanel();
    drawGraphsPanel();
    drawConvergencePanel();

    // Zone laissée libre au centre par les panneaux : c'est là que se dessine la 3D.
    const ImGuiIO& io = ImGui::GetIO();
    ImVec2 pos = vp->Pos, size = vp->Size;
    if (const ImGuiDockNode* central = ImGui::DockBuilderGetCentralNode(dockId)) {
        if (central->Size.x > 1.0f && central->Size.y > 1.0f) {
            pos = central->Pos;
            size = central->Size;
        }
    }
    viewRect_ = {static_cast<int>(pos.x * io.DisplayFramebufferScale.x),
                 static_cast<int>(pos.y * io.DisplayFramebufferScale.y),
                 std::max(1, static_cast<int>(size.x * io.DisplayFramebufferScale.x)),
                 std::max(1, static_cast<int>(size.y * io.DisplayFramebufferScale.y))};
}

// ------------------------------ rendu 3D --------------------------------

void App::drawScene(int fbW, int fbH) {
    const float clear[3] = {0.055f, 0.06f, 0.075f};
    renderer_.beginFrame(fbW, fbH, viewRect_, camera_, clear);

    static const std::vector<Vertex> grid = makeGrid(25, 2.0f);
    static const std::vector<Vertex> axes = makeAxes(4.0f);
    renderer_.draw(Primitive::Lines, grid);
    renderer_.draw(Primitive::Lines, axes);

    renderer_.draw(Primitive::LineStrip, refPath_);  // solution exacte
    for (int i = 0; i < kRunCount; ++i)
        if (isShown(i)) renderer_.draw(Primitive::LineStrip, runs_[i].trail);

    // Billes : l'exacte en bleu (grosse), les numériques par-dessus (plus petites).
    const Vec3 exact = problem_.position(runs_[0].world.time);
    const std::vector<Vertex> exactBall = {{static_cast<float>(exact.x), static_cast<float>(exact.y),
                                            static_cast<float>(exact.z), kBlue[0], kBlue[1], kBlue[2]}};
    renderer_.draw(Primitive::Points, exactBall, 18.0f * uiScale_);
    for (int i = 0; i < kRunCount; ++i) {
        if (!isShown(i)) continue;
        const Vec3 p = runs_[i].world.particles[0].position;
        const std::vector<Vertex> ball = {{static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z),
                                           runs_[i].color[0], runs_[i].color[1], runs_[i].color[2]}};
        renderer_.draw(Primitive::Points, ball, 11.0f * uiScale_);
    }
    glDisable(GL_SCISSOR_TEST);
}

// ------------------------------ boucle ---------------------------------

int App::run(const AppOptions& options) {
    if (!initWindow(!options.smokeTest)) {
        shutdown();
        return 1;
    }
    level_ = static_cast<Level>(options.level - 1);
    initImGui();
    initRuns();
    resetSimulation();

    using Clock = std::chrono::steady_clock;
    auto last = Clock::now();
    int frames = 0;

    while (!glfwWindowShouldClose(window_)) {
        glfwPollEvents();

        const auto now = Clock::now();
        const double frameSeconds = std::chrono::duration<double>(now - last).count();
        last = now;

        int fbW = 0, fbH = 0;
        glfwGetFramebufferSize(window_, &fbW, &fbH);
        if (fbW == 0 || fbH == 0) {  // fenêtre réduite
            glfwWaitEventsTimeout(0.1);
            continue;
        }

        stepSimulation(frameSeconds);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        buildUi();

        // Caméra : la souris pilote la scène 3D tant qu'elle n'est pas sur un panneau.
        const ImGuiIO& io = ImGui::GetIO();
        if (!io.WantCaptureMouse) {
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) camera_.orbit(io.MouseDelta.x, io.MouseDelta.y);
            if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseDown(ImGuiMouseButton_Middle))
                camera_.pan(io.MouseDelta.x, io.MouseDelta.y);
            if (io.MouseWheel != 0.0f) camera_.zoom(io.MouseWheel);
        }

        ImGui::Render();
        drawScene(fbW, fbH);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window_);

        if (options.smokeTest && ++frames >= 30) break;
    }

    if (options.smokeTest) std::puts("smoke-test : OK");
    shutdown();
    return 0;
}

}  // namespace

int runApplication(const AppOptions& options) {
    App app;
    return app.run(options);
}

}  // namespace pl
