#include "UiCommon.hpp"

#include <algorithm>

namespace pl {

bool sliderD(const char* label, double* v, double lo, double hi, const char* fmt, ImGuiSliderFlags flags) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &lo, &hi, fmt, flags);
}

void pushSliderWidth() {
    ImGui::PushItemWidth(-(ImGui::CalcTextSize("Tolérance relative RK45").x + ImGui::GetStyle().ItemInnerSpacing.x));
}
void popSliderWidth() { ImGui::PopItemWidth(); }

std::unique_ptr<Solver> SolverSet::makeFixedStep(int i) {
    switch (i) {
        case kEuler: return std::make_unique<ExplicitEuler>();
        case kSymplectic: return std::make_unique<SymplecticEuler>();
        case kVerlet: return std::make_unique<VelocityVerlet>();
        default: return std::make_unique<RK4>();
    }
}

SolverSet::SolverSet() {
    const float colors[kCount][3] = {{1.00f, 0.55f, 0.15f},   // Euler : orange
                                     {0.80f, 0.85f, 0.25f},   // symplectique : citron
                                     {0.85f, 0.45f, 0.95f},   // Verlet : violet
                                     {0.35f, 0.85f, 0.45f},   // RK4 : vert
                                     {1.00f, 0.40f, 0.50f}};  // RK45 : rose
    for (int i = 0; i < kCount; ++i) {
        if (i == kRK45) {
            auto rk = std::make_unique<RK45>();
            rk45_ = rk.get();
            entries_[i].solver = std::move(rk);
        } else {
            entries_[i].solver = makeFixedStep(i);
        }
        std::copy(colors[i], colors[i] + 3, entries_[i].color);
    }
}

bool SolverSet::isShown(Level level, int i) const {
    if (!atLeast(level, Level::College)) return i == kRK4;
    if (!atLeast(level, Level::Etudiant)) return i == kEuler || i == kRK4;
    return entries_[i].show;
}

std::string SolverSet::label(Level level, int i) const {
    if (!atLeast(level, Level::College) && i == kRK4) return "Ordinateur";
    if (!atLeast(level, Level::Etudiant)) {
        if (i == kEuler) return "Méthode simple (Euler)";
        if (i == kRK4) return "Méthode précise (RK4)";
    }
    return entries_[i].solver->name();
}

std::string SolverSet::shortLabel(Level level, int i) const {
    if (!atLeast(level, Level::Etudiant)) return label(level, i);
    static const char* names[kCount] = {"Euler", "Euler sympl.", "Verlet", "RK4", "RK45"};
    return names[i];
}

void SolverSet::drawToggles(Level level) {
    if (!atLeast(level, Level::Etudiant)) return;
    ImGui::TextDisabled("Méthodes affichées");
    for (int i = 0; i < kCount; ++i) {
        ImGui::PushStyleColor(ImGuiCol_CheckMark, toImVec4(entries_[i].color));
        ImGui::Checkbox(label(level, i).c_str(), &entries_[i].show);
        ImGui::PopStyleColor();
    }
}

}  // namespace pl
