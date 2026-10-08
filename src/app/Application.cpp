#include "Application.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
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
#include "physicslab/render/Camera.hpp"
#include "physicslab/render/Renderer.hpp"

#ifndef PHYSICSLAB_SHADER_DIR
#define PHYSICSLAB_SHADER_DIR "shaders"
#endif

namespace pl {
namespace {

// ---------------------------------------------------------------------------
// Démo du socle (M0) : lancer d'une particule sous pesanteur uniforme, intégré par
// Euler explicite et comparé en direct à la solution analytique. M1 en fera un vrai module.
// ---------------------------------------------------------------------------

struct LaunchParams {
    double speed = 14.0;      // vitesse de lancer [m/s]
    double angleDeg = 55.0;   // angle au-dessus de l'horizontale [deg]
    double height = 2.0;      // hauteur de départ [m]
    double gravity = constants::g0;  // intensité de la pesanteur [m/s^2]
    double mass = 1.0;        // [kg]
};

// Courbes temps réel dans un tampon circulaire (compatible avec l'argument `offset` d'ImPlot).
struct Trace {
    static constexpr int kCapacity = 6000;
    std::vector<double> t, y, yRef, dE;
    int offset = 0;

    int size() const { return static_cast<int>(t.size()); }
    void clear() {
        t.clear(); y.clear(); yRef.clear(); dE.clear();
        offset = 0;
    }
    void add(double tt, double yy, double rr, double ee) {
        if (size() < kCapacity) {
            t.push_back(tt); y.push_back(yy); yRef.push_back(rr); dE.push_back(ee);
        } else {
            t[offset] = tt; y[offset] = yy; yRef[offset] = rr; dE[offset] = ee;
            offset = (offset + 1) % kCapacity;
        }
    }
};

const char* explanation(Level level) {
    switch (level) {
        case Level::Vulgarisation:
            return "Quand on lance une balle, elle monte, ralentit, puis retombe : la Terre l'attire vers le bas. "
                   "Lance plus fort et elle ira plus loin et plus haut.\n\n"
                   "La ligne bleue montre ce que fait vraiment la nature. Le point orange est ce que calcule "
                   "l'ordinateur, pas à pas.";
        case Level::Interesse:
            return "La gravité est une force qui tire tout vers le sol. Sur la Lune elle est six fois plus faible : "
                   "la même balle irait beaucoup plus loin.\n\n"
                   "Le mouvement vers l'avant et le mouvement vers le haut sont indépendants : la balle avance "
                   "à vitesse constante pendant qu'elle monte puis redescend.";
        case Level::College:
            return "Poids : P = m x g  (g = 9,81 N/kg sur Terre).\n"
                   "Vitesse : v = d / t.\n"
                   "Pendant la chute, la vitesse verticale change de g = 9,81 m/s chaque seconde.\n\n"
                   "Énergie de mouvement : Ec = 1/2 x m x v^2.\n"
                   "Énergie de hauteur : Ep = m x g x h.\n"
                   "Pendant le vol, Ec + Ep reste (presque) constante.";
        case Level::Lycee:
            return "2e loi de Newton : somme des forces = m a (vecteurs). Seule la pesanteur agit : a = (0 ; -g).\n\n"
                   "En intégrant :\n"
                   "  vx = v0x          x = x0 + v0x t\n"
                   "  vy = v0y - g t    y = y0 + v0y t - 1/2 g t^2\n\n"
                   "Énergie mécanique Em = Ec + Ep conservée ; quantité de mouvement horizontale conservée.";
        case Level::Etudiant:
            return "Lagrangien : L = T - V = 1/2 m (x'^2 + y'^2) - m g y.\n"
                   "Euler-Lagrange : d/dt (dL/dq') - dL/dq = 0  =>  x'' = 0, y'' = -g.\n"
                   "Hamiltonien : H = p^2 / 2m + m g y.\n\n"
                   "Le solveur est Euler explicite, qui n'est pas symplectique : H dérive au cours du temps. "
                   "Réduisez dt et observez la dérive diminuer proportionnellement.";
        case Level::Chercheur:
            return "Euler explicite : y_{n+1} = y_n + dt f(t_n, y_n). Erreur locale O(dt^2), globale O(dt).\n\n"
                   "Pour la chute libre, la dérive d'énergie se calcule à la main : "
                   "dE = 1/2 m g^2 t dt (vérifiée par le test unitaire à 1e-6 près).\n\n"
                   "Comparez cette valeur à la dérive mesurée dans le panneau Invariants. "
                   "Prochaine étape (M1) : Verlet, RK4 et RK45 adaptatif.";
    }
    return "";
}

bool sliderD(const char* label, double* v, double lo, double hi, const char* fmt, ImGuiSliderFlags flags = 0) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &lo, &hi, fmt, flags);
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
    void shutdown();

    void resetSimulation();
    void stepSimulation(double frameSeconds);
    Vec3 analyticPosition(double t) const;
    double analyticFlightTime() const;
    void recordSample();

    void buildUi();
    void buildDefaultLayout(ImGuiID dockId, const ImGuiViewport* viewport);
    void drawMenuBar();
    void drawSimulationPanel();
    void drawExplanationPanel();
    void drawInvariantsPanel();
    void drawGraphsPanel();
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

    // simulation
    LaunchParams params_;
    World world_;
    ExplicitEuler solver_;
    Invariants initial_;
    double dt_ = 1.0 / 240.0;
    double timeScale_ = 1.0;
    double accumulator_ = 0.0;
    bool running_ = true;
    bool landed_ = false;
    double landingTime_ = 0.0;
    Trace trace_;
    std::vector<Vertex> trail_;      // trajectoire numérique
    std::vector<Vertex> refPath_;    // trajectoire analytique
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

// ------------------------------ simulation -----------------------------

Vec3 App::analyticPosition(double t) const {
    const double a = params_.angleDeg * constants::pi / 180.0;
    const Vec3 r0{0.0, params_.height, 0.0};
    const Vec3 v0{params_.speed * std::cos(a), params_.speed * std::sin(a), 0.0};
    const Vec3 g{0.0, -params_.gravity, 0.0};
    return r0 + t * v0 + 0.5 * t * t * g;  // r(t) = r0 + v0 t + g t^2 / 2
}

double App::analyticFlightTime() const {
    const double vy = params_.speed * std::sin(params_.angleDeg * constants::pi / 180.0);
    return (vy + std::sqrt(vy * vy + 2.0 * params_.gravity * params_.height)) / params_.gravity;
}

void App::resetSimulation() {
    const double a = params_.angleDeg * constants::pi / 180.0;

    world_ = World{};
    world_.gravity = {0.0, -params_.gravity, 0.0};
    Particle p;
    p.mass = params_.mass;
    p.position = {0.0, params_.height, 0.0};
    p.velocity = {params_.speed * std::cos(a), params_.speed * std::sin(a), 0.0};
    world_.particles.push_back(p);

    initial_ = world_.invariants();
    accumulator_ = 0.0;
    landed_ = false;
    running_ = true;
    trace_.clear();
    trail_.clear();

    refPath_.clear();
    const double tEnd = analyticFlightTime();
    for (int i = 0; i <= 120; ++i) {
        const Vec3 r = analyticPosition(tEnd * i / 120.0);
        refPath_.push_back({static_cast<float>(r.x), static_cast<float>(r.y), static_cast<float>(r.z),
                            0.35f, 0.65f, 1.0f});
    }
    recordSample();
}

void App::recordSample() {
    const Vec3 p = world_.particles[0].position;
    const double yRef = analyticPosition(world_.time).y;
    trace_.add(world_.time, p.y, yRef, world_.invariants().total() - initial_.total());
    trail_.push_back({static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z),
                      1.0f, 0.55f, 0.15f});
}

void App::stepSimulation(double frameSeconds) {
    if (!running_ || landed_) return;

    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;  // borne anti "spirale de la mort"
    int guard = 0;
    while (accumulator_ >= dt_ && guard++ < 2000) {
        world_.step(solver_, dt_);
        accumulator_ -= dt_;
        recordSample();
        if (world_.particles[0].position.y < 0.0) {  // impact au sol
            landed_ = true;
            landingTime_ = world_.time;
            break;
        }
    }
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
    int idx = static_cast<int>(level_);
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
    const ImGuiID bottom = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Down, 0.32f, nullptr, &rest);
    const ImGuiID right = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Right, 0.24f, nullptr, &rest);

    ImGui::DockBuilderDockWindow("Simulation", left);
    ImGui::DockBuilderDockWindow("Explication", left);
    ImGui::DockBuilderDockWindow("Graphes", bottom);
    ImGui::DockBuilderDockWindow("Invariants", right);
    ImGui::DockBuilderFinish(dockId);
}

void App::drawSimulationPanel() {
    ImGui::Begin("Simulation");

    bool changed = false;
    ImGui::SeparatorText("Lancer");
    changed |= sliderD(atLeast(level_, Level::College) ? "Vitesse (m/s)" : "Force du lancer",
                       &params_.speed, 1.0, 30.0, atLeast(level_, Level::College) ? "%.1f" : "");
    if (atLeast(level_, Level::Interesse)) changed |= sliderD("Angle (deg)", &params_.angleDeg, 5.0, 85.0, "%.0f");

    if (level_ == Level::Interesse) {
        static const char* planets[] = {"Terre", "Lune", "Mars", "Jupiter"};
        static const double gs[] = {constants::g0, 1.62, 3.71, 24.79};
        static int planet = 0;
        if (ImGui::Combo("Planète", &planet, planets, 4)) {
            params_.gravity = gs[planet];
            changed = true;
        }
    }
    if (atLeast(level_, Level::College)) {
        changed |= sliderD("Pesanteur g (m/s²)", &params_.gravity, 0.5, 30.0, "%.2f");
        changed |= sliderD("Masse (kg)", &params_.mass, 0.1, 20.0, "%.1f", ImGuiSliderFlags_Logarithmic);
    }
    if (atLeast(level_, Level::Lycee)) {
        changed |= sliderD("Hauteur initiale (m)", &params_.height, 0.0, 20.0, "%.1f");
        const double a = params_.angleDeg * constants::pi / 180.0;
        ImGui::TextDisabled("v0 = (%.2f ; %.2f) m/s", params_.speed * std::cos(a), params_.speed * std::sin(a));
    }
    if (changed) resetSimulation();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !landed_ ? "Pause" : "Lecture")) running_ = !running_;
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) resetSimulation();

    if (atLeast(level_, Level::Interesse)) sliderD("Vitesse du temps", &timeScale_, 0.05, 2.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    if (atLeast(level_, Level::Etudiant)) {
        sliderD("Pas dt (s)", &dt_, 1e-4, 0.05, "%.4f", ImGuiSliderFlags_Logarithmic);
        ImGui::TextDisabled("Solveur : %s", solver_.name());
    }

    if (atLeast(level_, Level::College)) ImGui::Text("t = %.3f s", world_.time);
    const Particle& p = world_.particles[0];
    if (atLeast(level_, Level::College)) {
        ImGui::Text("x = %.2f m   y = %.2f m", p.position.x, p.position.y);
        ImGui::Text("vitesse = %.2f m/s", p.velocity.norm());
    }
    if (landed_) {
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Impact au sol !");
        if (atLeast(level_, Level::College)) {
            ImGui::Text("portée = %.2f m, durée = %.3f s", p.position.x, landingTime_);
            if (atLeast(level_, Level::Lycee)) ImGui::TextDisabled("durée analytique = %.3f s", analyticFlightTime());
        }
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
    const Invariants inv = world_.invariants();

    if (!atLeast(level_, Level::College)) {
        ImGui::TextWrapped("Pendant le vol, l'énergie du mouvement et l'énergie de la hauteur s'échangent, "
                           "mais leur total reste presque le même.");
        ImGui::End();
        return;
    }

    ImGui::Text("Ec  = %9.3f J", inv.kinetic);
    ImGui::Text("Ep  = %9.3f J", inv.potential);
    ImGui::Text("Em  = %9.3f J", inv.total());

    if (atLeast(level_, Level::Lycee)) {
        ImGui::SeparatorText("Quantité de mouvement");
        ImGui::Text("px = %8.3f kg.m/s", inv.momentum.x);
        ImGui::Text("py = %8.3f kg.m/s", inv.momentum.y);
        ImGui::SeparatorText("Moment cinétique / O");
        ImGui::Text("Lz = %8.3f kg.m²/s", inv.angularMomentum.z);
    }
    if (atLeast(level_, Level::Etudiant)) {
        const double drift = inv.total() - initial_.total();
        ImGui::SeparatorText("Dérive d'énergie");
        ImGui::Text("dE = %+.4e J", drift);
        ImGui::Text("dE/E0 = %+.3e", drift / initial_.total());
    }
    if (atLeast(level_, Level::Chercheur)) {
        const double g = params_.gravity;
        const double predicted = 0.5 * params_.mass * g * g * world_.time * dt_;  // dE = 1/2 m g^2 t dt
        ImGui::TextDisabled("dE théorique = %+.4e J", predicted);
        ImGui::TextDisabled("(1/2 m g² t dt)");
    }
    ImGui::End();
}

void App::drawGraphsPanel() {
    ImGui::Begin("Graphes");
    const bool showEnergy = atLeast(level_, Level::Lycee);
    const float h = ImGui::GetContentRegionAvail().y;
    const float plotH = showEnergy ? std::max(80.0f, h * 0.5f - 4.0f) : h;

    // Mêmes couleurs que dans la scène 3D : bleu = analytique, orange = numérique.
    const ImVec4 blue(0.35f, 0.65f, 1.0f, 1.0f), orange(1.0f, 0.55f, 0.15f, 1.0f);
    const ImPlotSpec analytic(ImPlotProp_LineColor, blue, ImPlotProp_LineWeight, 2.0f,
                              ImPlotProp_Offset, trace_.offset);
    const ImPlotSpec numeric(ImPlotProp_LineColor, orange, ImPlotProp_LineWeight, 2.0f,
                             ImPlotProp_Offset, trace_.offset);

    if (ImPlot::BeginPlot("Hauteur y(t)", ImVec2(-1.0f, plotH))) {
        ImPlot::SetupAxes("t (s)", "y (m)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        if (trace_.size() > 0) {
            ImPlot::PlotLine("Analytique", trace_.t.data(), trace_.yRef.data(), trace_.size(), analytic);
            ImPlot::PlotLine("Numérique", trace_.t.data(), trace_.y.data(), trace_.size(), numeric);
        }
        ImPlot::EndPlot();
    }
    if (showEnergy && ImPlot::BeginPlot("Dérive d'énergie Em(t) - Em(0)", ImVec2(-1.0f, plotH))) {
        ImPlot::SetupAxes("t (s)", "dE (J)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        if (trace_.size() > 0)
            ImPlot::PlotLine("dE", trace_.t.data(), trace_.dE.data(), trace_.size(), numeric);
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

    renderer_.draw(Primitive::LineStrip, refPath_);   // solution analytique
    renderer_.draw(Primitive::LineStrip, trail_);     // trajectoire numérique

    const Vec3 p = world_.particles[0].position;
    const std::vector<Vertex> ball = {{static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z),
                                       1.0f, 0.6f, 0.2f}};
    renderer_.draw(Primitive::Points, ball, 16.0f * uiScale_);
    glDisable(GL_SCISSOR_TEST);
}

// ------------------------------ boucle ---------------------------------

int App::run(const AppOptions& options) {
    if (!initWindow(!options.smokeTest)) {
        shutdown();
        return 1;
    }
    initImGui();
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
