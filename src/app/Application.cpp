#include "Application.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <vector>

#include <glad/gl.h>
// glad doit précéder GLFW
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_internal.h>  // DockBuilder (disposition par défaut)
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>

#include "DoublePendulumModule.hpp"
#include "KeplerModule.hpp"
#include "BounceModule.hpp"
#include "CradleModule.hpp"
#include "RigidBodyModule.hpp"
#include "FrictionModule.hpp"
#include "GpuNBodyModule.hpp"
#include "GpuTest.hpp"
#include "NBodyModule.hpp"
#include "OscillatorModule.hpp"
#include "PendulumModule.hpp"
#include "ProjectileModule.hpp"
#include "SimulationModule.hpp"
#include "physicslab/core/Level.hpp"
#include "physicslab/render/Camera.hpp"
#include "physicslab/render/Renderer.hpp"

#ifndef PHYSICSLAB_SHADER_DIR
#define PHYSICSLAB_SHADER_DIR "shaders"
#endif

namespace pl {
namespace {

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

    void buildUi();
    void buildDefaultLayout(ImGuiID dockId, const ImGuiViewport* viewport);
    void drawMenuBar();
    void activate(int index);

    // fenêtre et rendu
    GLFWwindow* window_ = nullptr;
    Renderer renderer_;
    Camera camera_;
    float uiScale_ = 1.0f;
    ViewRect viewRect_{0, 0, 1, 1};
    bool resetLayout_ = false;

    // pédagogie et simulations
    Level level_ = Level::College;
    std::vector<std::unique_ptr<SimulationModule>> modules_;
    int active_ = 0;
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

void App::activate(int index) {
    active_ = index;
    modules_[active_]->frameCamera(camera_);
}

// --------------------------------- UI ----------------------------------

void App::drawMenuBar() {
    if (!ImGui::BeginMainMenuBar()) return;

    if (ImGui::BeginMenu("Fichier")) {
        if (ImGui::MenuItem("Quitter")) glfwSetWindowShouldClose(window_, GLFW_TRUE);
        ImGui::EndMenu();
    }
    if (ImGui::BeginMenu("Simulation")) {
        ImGui::TextDisabled("Mécanique");
        for (int i = 0; i < static_cast<int>(modules_.size()); ++i)
            if (ImGui::MenuItem(modules_[i]->title(), nullptr, i == active_) && i != active_) activate(i);
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
    ImGui::DockBuilderDockWindow("Analyse", bottom);
    ImGui::DockBuilderDockWindow("Invariants", right);
    ImGui::DockBuilderFinish(dockId);
}

void App::buildUi() {
    drawMenuBar();

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const ImGuiID dockId = ImGui::DockSpaceOverViewport(0, vp, ImGuiDockNodeFlags_PassthruCentralNode);
    if (resetLayout_) {
        buildDefaultLayout(dockId, vp);
        resetLayout_ = false;
    }

    // L'application ouvre les fenêtres, le module actif les remplit.
    SimulationModule& sim = *modules_[active_];
    const UiContext ctx{level_, uiScale_};

    ImGui::Begin("Simulation");
    sim.drawControls(ctx);
    ImGui::End();

    ImGui::Begin("Explication");
    ImGui::TextColored(ImVec4(0.45f, 0.80f, 0.95f, 1.0f), "%s", sim.title());
    ImGui::TextDisabled("Niveau %s", levelName(level_));
    ImGui::Separator();
    ImGui::TextWrapped("%s", sim.explanation(level_));
    ImGui::End();

    ImGui::Begin("Invariants");
    sim.drawInvariants(ctx);
    ImGui::End();

    ImGui::Begin("Graphes");
    sim.drawGraphs(ctx);
    ImGui::End();

    if (sim.hasAnalysis(ctx)) {
        ImGui::Begin("Analyse");
        sim.drawAnalysis(ctx);
        ImGui::End();
    }

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

// ------------------------------ boucle ---------------------------------

int App::run(const AppOptions& options) {
    if (!initWindow(!options.smokeTest && !options.gpuTest)) {
        shutdown();
        return 1;
    }
    if (options.gpuTest) {  // pas d'interface : seulement le contexte OpenGL, pour le calcul GPU
        const int code = runGpuTest(PHYSICSLAB_SHADER_DIR, options.gpuMaxN);
        shutdown();
        return code;
    }
    level_ = static_cast<Level>(options.level - 1);
    initImGui();

    modules_.push_back(std::make_unique<ProjectileModule>());
    modules_.push_back(std::make_unique<OscillatorModule>());
    modules_.push_back(std::make_unique<PendulumModule>());
    modules_.push_back(std::make_unique<DoublePendulumModule>());
    modules_.push_back(std::make_unique<KeplerModule>());
    modules_.push_back(std::make_unique<NBodyModule>());
    modules_.push_back(std::make_unique<FrictionModule>());
    modules_.push_back(std::make_unique<BounceModule>());
    modules_.push_back(std::make_unique<CradleModule>());
    modules_.push_back(std::make_unique<RigidBodyModule>());
    modules_.push_back(std::make_unique<GpuNBodyModule>(PHYSICSLAB_SHADER_DIR));
    activate(std::clamp(options.simulation - 1, 0, static_cast<int>(modules_.size()) - 1));

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

        modules_[active_]->update(frameSeconds);

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

        const float clear[3] = {0.055f, 0.06f, 0.075f};
        renderer_.beginFrame(fbW, fbH, viewRect_, camera_, clear);
        modules_[active_]->drawScene(renderer_, UiContext{level_, uiScale_});
        glDisable(GL_SCISSOR_TEST);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window_);

        if (options.smokeTest && ++frames >= 30) break;
    }

    if (options.smokeTest) std::puts("smoke-test : OK");
    modules_.clear();
    shutdown();
    return 0;
}

}  // namespace

int runApplication(const AppOptions& options) {
    App app;
    return app.run(options);
}

}  // namespace pl
