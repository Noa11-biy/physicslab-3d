#include "StringModule.hpp"

#include <algorithm>
#include <cmath>

#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/waves/Fft.hpp"

namespace pl {
namespace {

constexpr double kTwoPi = 2.0 * constants::pi;
constexpr float kComputed[3] = {0.35f, 0.85f, 0.45f};   // vert : le calcul
constexpr float kInitial[3] = {0.45f, 0.48f, 0.55f};    // gris : la forme au départ
constexpr float kProbeColor[3] = {1.00f, 0.75f, 0.25f}; // orange : la masse observée
constexpr float kPostColor[3] = {0.92f, 0.92f, 0.95f};
constexpr float kRest[3] = {0.35f, 0.38f, 0.42f};
constexpr int kShownModes = 12;                         // modes tracés dans la courbe d'énergie par mode
constexpr double kDefaultDtFactor = 0.2;                // dt = 0.2 a / c : RK4 perd alors ~1e-4 d'énergie en 4 s (0,6 % à 0.4 : dissipation des modes aigus, ordre (w dt)^6)

const char* kSolverNames[] = {"Euler explicite", "Euler symplectique", "Verlet", "RK4"};

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

}  // namespace

StringModule::StringModule() { reset(); }

int StringModule::probeBead() const { return std::clamp(static_cast<int>(std::lround(0.29 * (beads_ + 1))), 1, beads_) - 1; }

double StringModule::maxAbsPositions() const {
    double m = 0.0;
    for (int j = 0; j < beads_; ++j) {
        const double a = std::abs(y_[j]);
        if (!(a <= m)) m = a;  // un NaN se propage : un calcul qui a explosé ne passe pas pour un calcul calme
    }
    return m;
}

// Limite de dt en unités de a / c = 2 / w_max : Verlet et Euler symplectique 1, RK4 sqrt(2) ; Euler explicite : aucune.
double StringModule::stabilityFactorLimit() const {
    switch (solverIndex_) {
        case SolverSet::kEuler: return 0.0;
        case SolverSet::kRK4: return std::sqrt(2.0);
        default: return 1.0;
    }
}

// ------------------------------ simulation -----------------------------

void StringModule::reset() {
    problem_.beads = beads_;
    problem_.length = length_;
    problem_.tension = tension_;
    problem_.density = density_;
    modes_ = std::make_unique<waves::StringModes>(problem_);
    u0_ = startMode_ > 0 ? problem_.modeShape(startMode_, height_) : problem_.pluck(pluckPos_ * length_, height_);
    v0_.assign(beads_, 0.0);
    y_ = waves::StringProblem::state(u0_, v0_);
    exact_ = y_;
    rhs_ = problem_.rhs();
    solver_ = SolverSet::makeFixedStep(solverIndex_);
    dt_ = dtFactor_ * problem_.verletDtLimit();
    energy0_ = std::max(problem_.energy(y_), 1e-300);
    endTime_ = 16.0 / problem_.fundamental();  // 16 périodes du mode 1 : une résolution de f1 / 16 sur le spectre

    time_ = accumulator_ = lastSample_ = 0.0;
    finished_ = diverged_ = false;
    probe_.clear();
    probeExact_.clear();
    drift_.clear();
    record_.clear();
    record_.push_back(y_[probeBead()]);
    spectrum_.clear();
    readFreq_.clear();
    spectrumAge_ = 1000;
    sample();
}

void StringModule::sample() {
    lastSample_ = time_;
    if (!(maxAbsPositions() < 1e3)) {  // la corde « explose » (pas trop grand) : on arrête sans polluer les courbes
        diverged_ = true;
        return;
    }
    exact_ = modes_->exact(u0_, v0_, time_);
    const int j = probeBead();
    probe_.add(time_, y_[j]);
    probeExact_.add(time_, exact_[j]);
    drift_.add(time_, std::max(std::abs(problem_.energy(y_) / energy0_ - 1.0), 1e-16));
}

void StringModule::update(double frameSeconds) {
    if (!running_ || finished_ || diverged_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    int guard = 0;
    while (accumulator_ >= dt_ && guard++ < 600 && !finished_ && !diverged_) {
        solver_->step(rhs_, time_, y_, dt_);
        time_ += dt_;
        accumulator_ -= dt_;
        if (record_.size() < kMaxRecord) record_.push_back(y_[probeBead()]);
        if (time_ - lastSample_ >= 1.0 / 120.0) sample();
        if (time_ >= endTime_) {
            finished_ = true;
            sample();
        }
    }
    // la solution exacte à l'instant EXACT de l'image (sample() n'est appelé qu'à 120 Hz : sans cela l'écart affiché mélangerait
    // l'erreur du schéma et un décalage de temps d'au plus 1/120 s)
    if (!diverged_ && time_ != lastSample_) exact_ = modes_->exact(u0_, v0_, time_);
}

void StringModule::refreshSpectrum() {
    if (++spectrumAge_ < 8 || record_.size() < 16) return;
    spectrumAge_ = 0;
    spectrum_ = waves::windowedSpectrum(record_, 1.0 / dt_, 4, &binWidth_);
    readFreq_.assign(6, 0.0);
    // modes 1 à 6 : lus seulement quand l'enregistrement couvre au moins 8 périodes du mode 1
    if (static_cast<double>(record_.size()) * dt_ >= 8.0 / problem_.fundamental())
        for (int n = 1; n <= std::min(beads_, 6); ++n)
            readFreq_[n - 1] = waves::peakFromSpectrum(spectrum_, binWidth_, problem_.chainOmega(n) / kTwoPi, 0.1);
}

// --------------------------------- UI ----------------------------------

void StringModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (!etudiant) {  // aux niveaux 1 à 4 : la corde continue (59 masses) et le calcul précis (RK4, pas moyen) ; les réglages de calcul viennent au niveau 5
        if (beads_ != 59 || solverIndex_ != SolverSet::kRK4 || dtFactor_ != kDefaultDtFactor) {
            beads_ = 59;
            solverIndex_ = SolverSet::kRK4;
            dtFactor_ = kDefaultDtFactor;
            changed = true;
        }
    }
    if (!interesse && startMode_ != 0) {
        startMode_ = 0;
        changed = true;
    }

    pushSliderWidth();
    ImGui::SeparatorText("Corde");
    if (startMode_ == 0)
        changed |= sliderD(college ? "Pincement x0 / L" : "Où pincer la corde", &pluckPos_, 0.05, 0.5, college ? "%.3f" : "");
    popSliderWidth();

    if (interesse) {
        ImGui::TextDisabled("Départ");
        int pure = startMode_ > 0 ? 1 : 0;
        changed |= ImGui::RadioButton("Pincer la corde", &pure, 0);
        ImGui::SameLine();
        changed |= ImGui::RadioButton("Un seul mode", &pure, 1);
        if (pure && startMode_ == 0) startMode_ = 1;
        if (!pure) startMode_ = 0;
        pushSliderWidth();
        if (startMode_ > 0) changed |= ImGui::SliderInt(college ? "Mode n (n arcs)" : "Nombre d'arcs", &startMode_, 1, kShownModes);
        popSliderWidth();

        if (startMode_ == 0) {
            ImGui::TextDisabled("Pincer à");
            const struct { const char* label; double value; } presets[] = {{"L/2", 0.5}, {"L/3", 1.0 / 3.0}, {"L/4", 0.25}, {"L/5", 0.2}};
            for (int i = 0; i < 4; ++i) {
                if (i > 0) ImGui::SameLine();
                if (ImGui::Button(presets[i].label)) {
                    pluckPos_ = presets[i].value;
                    changed = true;
                }
            }
        }
    }

    pushSliderWidth();
    if (interesse) {
        changed |= sliderD(college ? "Tension T (N)" : "Tension de la corde", &tension_, 0.5, 20.0, college ? "%.2f" : "", ImGuiSliderFlags_Logarithmic);
        changed |= sliderD(college ? "Masse linéique μ (kg/m)" : "Épaisseur (masse)", &density_, 0.2, 5.0, college ? "%.2f" : "", ImGuiSliderFlags_Logarithmic);
    }
    if (college) {
        changed |= sliderD("Longueur L (m)", &length_, 0.5, 2.0, "%.2f");
        changed |= sliderD(startMode_ > 0 ? "Amplitude (m)" : "Hauteur h (m)", &height_, 0.02, 0.25, "%.3f");
    }
    popSliderWidth();

    if (etudiant) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        changed |= ImGui::SliderInt("Nombre de masses N", &beads_, 8, 200);
        changed |= sliderD("Pas dt (× a / c)", &dtFactor_, 0.05, 1.8, "%.3f", ImGuiSliderFlags_Logarithmic);
        popSliderWidth();
        ImGui::TextDisabled("Méthode");
        for (int i = 0; i < 4; ++i) changed |= ImGui::RadioButton(kSolverNames[i], &solverIndex_, i);
        const double limit = stabilityFactorLimit();
        const bool unstable = limit == 0.0 || dtFactor_ > limit;
        if (unstable) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
            if (limit == 0.0) ImGui::TextWrapped("Euler explicite est instable pour tout pas : l'énergie grandit à chaque pas.");
            else ImGui::TextWrapped("dt = %.3f × a/c dépasse la limite %.3f de cette méthode : le calcul est instable.", dtFactor_, limit);
            ImGui::PopStyleColor();
        }
        if (beads_ % 60 != 59 && startMode_ == 0)
            wrapped("Astuce : avec N + 1 multiple de 60 (N = 59, 119, 179), un pincement à L/2, L/3, L/4, L/5 tombe exactement sur une masse : les harmoniques absents sont alors exactement nuls.", true);
    }
    if (changed) reset();

    ImGui::SeparatorText("Temps");
    if (ImGui::Button(running_ && !finished_ ? "Pause" : "Lecture")) {
        if (finished_) reset();
        else running_ = !running_;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recommencer")) reset();
    pushSliderWidth();
    if (interesse) sliderD("Vitesse du temps", &timeScale_, 0.05, 4.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (diverged_) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
        ImGui::TextWrapped("Le calcul a explosé : le pas de temps est trop grand pour cette méthode (voir « Calcul »).");
        ImGui::PopStyleColor();
    }

    // grandeurs
    if (interesse) {
        ImGui::SeparatorText("Grandeurs");
        const double f1 = problem_.fundamental();
        if (!college) {
            ImGui::Text("La corde vibre %.2f fois par seconde.", f1);
        } else {
            ImGui::Text("t = %.2f s / %.2f s", time_, endTime_);
            ImGui::Text("Vitesse de l'onde c = √(T/μ) = %.3f m/s", problem_.speed());
            ImGui::Text("Fréquence de base f1 = c/(2L) = %.3f Hz", f1);
            ImGui::Text("Harmoniques : f_n = n × f1");
            ImGui::TextDisabled("%.2f  %.2f  %.2f  %.2f  %.2f Hz", f1, 2 * f1, 3 * f1, 4 * f1, 5 * f1);
        }
        if (lycee) ImGui::Text("Longueur d'onde du mode n : λ_n = 2L/n (mode 1 : %.2f m)", 2.0 * length_);
        if (etudiant) {
            ImGui::Text("a = L/(N+1) = %.4f m   m = μ a = %.4f kg", problem_.spacing(), problem_.mass());
            ImGui::Text("k = T/a = %.2f N/m   dt = %.3f ms", problem_.stiffness(), 1e3 * dt_);
        }
        if (chercheur) {
            ImGui::Text("Coupure ω_max = 2c/a = %.1f rad/s (%.2f Hz)", problem_.maxOmega(), problem_.maxOmega() / kTwoPi);
            ImGui::Text("dt limite : Verlet a/c = %.3f ms, RK4 √2 a/c = %.3f ms", 1e3 * problem_.verletDtLimit(), 1e3 * problem_.rk4DtLimit());
        }
    }
}

void StringModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (diverged_) {
        wrapped("Calcul divergent : plus rien à mesurer. Baisser dt ou changer de méthode, puis « Recommencer ».", true);
        return;
    }
    if (!atLeast(level, Level::College)) {
        wrapped("Sans frottement, la corde ne perd aucune énergie : elle va et vient entre mouvement (rapide) et tension (étirée).", true);
        return;
    }

    const double energy = problem_.energy(y_);
    ImGui::Text("Énergie E = %.5f J", energy);
    ImGui::Text("Cinétique %.5f J   Élastique %.5f J", problem_.kineticEnergy(y_), problem_.potentialEnergy(y_));
    if (atLeast(level, Level::Lycee)) ImGui::Text("E / E0 - 1 = %+.1e", energy / energy0_ - 1.0);
    if (atLeast(level, Level::Etudiant)) {
        ImGui::SeparatorText("Contre la solution exacte");
        double err = 0.0;
        for (int j = 0; j < beads_; ++j) err = std::max(err, std::abs(y_[j] - exact_[j]));
        ImGui::Text("%s : max |u - exacte| = %.2e m", kSolverNames[solverIndex_], err);
        ImGui::Text("Énergie de la solution exacte : %.5f J", problem_.energy(exact_));
        wrapped("Euler explicite fait grandir l'énergie ; Verlet et Euler symplectique la gardent bornée ; RK4 la perd très lentement.", true);
    }
    if (atLeast(level, Level::Chercheur)) {
        ImGui::SeparatorText("Pas de temps");
        ImGui::Text("ω_max dt = %.3f  (Verlet : < 2, RK4 : < 2√2 = 2.828)", problem_.maxOmega() * dt_);
        const std::vector<double> e = modes_->modeEnergies(y_);
        double sum = 0.0;
        for (double v : e) sum += v;
        ImGui::Text("Σ E_n / E - 1 = %+.1e", sum / energy - 1.0);
        wrapped("L'énergie se range exactement par mode : la somme des E_n redonne E.", true);

        // fréquences des modes 1 à 6 lues dans le spectre (onglet Analyse) contre celles de la chaîne
        ImGui::SeparatorText("Fréquences lues dans le spectre");
        refreshSpectrum();
        std::vector<TableRow> rows;
        for (int n = 1; n <= std::min(beads_, 6); ++n) {
            const double chain = problem_.chainOmega(n) / kTwoPi;
            const double read = n - 1 < static_cast<int>(readFreq_.size()) ? readFreq_[n - 1] : 0.0;
            rows.push_back({strf("mode %d", n), nullptr,
                            {strf("%.5f", chain), read > 0.0 ? strf("%.5f", read) : "-", read > 0.0 ? strf("%+.0e", read / chain - 1.0) : "-"}});
        }
        drawResultTable("frequences", {"f chaîne", "f lue", "écart"}, rows);
        if (readFreq_.empty() || readFreq_[0] <= 0.0) wrapped("(lues après 8 périodes du mode 1)", true);
    }
}

void StringModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showModes = atLeast(level, Level::College);
    const bool showTheory = atLeast(level, Level::Lycee);
    const bool showDrift = atLeast(level, Level::Etudiant);
    const int cols = 1 + (showModes ? 1 : 0) + (showDrift ? 1 : 0);
    if (diverged_) return;

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Mouvement d'une masse u(t)")) {
            ImPlot::SetupAxes("t (s)", "u (m)", ImPlotAxisFlags_AutoFit, 0);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -1.3 * height_, 1.3 * height_, ImPlotCond_Always);
            if (showTheory && probeExact_.size() > 0)
                ImPlot::PlotLine("Exacte (modes)", probeExact_.x.data(), probeExact_.y.data(), probeExact_.size(), lineSpec(kBlue, probeExact_.offset));
            if (probe_.size() > 0)
                ImPlot::PlotLine(atLeast(level, Level::College) ? "Calcul" : "Corde", probe_.x.data(), probe_.y.data(), probe_.size(), lineSpec(kComputed, probe_.offset));
            ImPlot::EndPlot();
        }
        if (showModes && ImPlot::BeginPlot("Part de l'énergie dans chaque mode")) {
            ImPlot::SetupAxes("mode n", "E_n / E", ImPlotAxisFlags_AutoFit, 0);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.05, ImPlotCond_Always);
            const int shown = std::min(beads_, kShownModes);
            std::vector<double> xs(shown), ys(shown), theory(shown);
            const std::vector<double> e = modes_->modeEnergies(y_);
            double total = 0.0;
            for (double v : e) total += v;
            for (int n = 1; n <= shown; ++n) {
                xs[n - 1] = n;
                ys[n - 1] = total > 0.0 ? e[n - 1] / total : 0.0;
            }
            if (showTheory) {
                // corde continue : l'énergie du mode n vaut (mu L / 4) w_n² b_n², soit proportionnelle à n² b_n² pour un pincement
                double norm = 0.0;
                std::vector<double> share(2000);
                for (int n = 1; n <= 2000; ++n) {
                    const double b = startMode_ > 0 ? 0.0 : problem_.pluckCoefficient(n, pluckPos_ * length_, height_);
                    share[n - 1] = static_cast<double>(n) * n * b * b;
                    norm += share[n - 1];
                }
                for (int n = 1; n <= shown; ++n) theory[n - 1] = startMode_ > 0 ? (n == startMode_ ? 1.0 : 0.0) : share[n - 1] / norm;
            }
            ImPlot::PlotBars("Calcul", xs.data(), ys.data(), shown, 0.6, ImPlotSpec(ImPlotProp_FillColor, toImVec4(kComputed, 0.8f)));
            if (showTheory)
                ImPlot::PlotScatter("Théorie (corde continue)", xs.data(), theory.data(), shown,
                                    ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 6.0f, ImPlotProp_MarkerFillColor, toImVec4(kBlue),
                                               ImPlotProp_MarkerLineColor, toImVec4(kBlue)));
            ImPlot::EndPlot();
        }
        if (showDrift && ImPlot::BeginPlot("Dérive de l'énergie |E/E0 - 1|")) {
            ImPlot::SetupAxes("t (s)", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (drift_.size() > 0) ImPlot::PlotLine(kSolverNames[solverIndex_], drift_.x.data(), drift_.y.data(), drift_.size(), lineSpec(kComputed, drift_.offset));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void StringModule::drawAnalysis(const UiContext&) {
    if (diverged_) return;
    const double f1 = problem_.fundamental();

    refreshSpectrum();  // spectre du mouvement de la masse observée (la FFT de plusieurs dizaines de milliers de points coûte quelques ms)
    std::vector<double> fs(spectrum_.size());
    for (std::size_t k = 0; k < fs.size(); ++k) fs[k] = static_cast<double>(k) * binWidth_;
    double peak = 0.0;
    for (double v : spectrum_) peak = std::max(peak, v);
    std::vector<double> normalized(spectrum_.size());
    for (std::size_t k = 0; k < normalized.size(); ++k) normalized[k] = peak > 0.0 ? spectrum_[k] / peak : 0.0;

    const int shown = std::min(beads_, kShownModes);
    std::vector<double> chainF(shown), nn(shown);
    for (int n = 1; n <= shown; ++n) {
        chainF[n - 1] = problem_.chainOmega(n) / kTwoPi;
        nn[n - 1] = n;
    }
    std::vector<double> allN(beads_), allChain(beads_), allCont(beads_);
    for (int n = 1; n <= beads_; ++n) {
        allN[n - 1] = n;
        allChain[n - 1] = problem_.chainOmega(n) / kTwoPi;
        allCont[n - 1] = n * f1;
    }

    if (ImPlot::BeginSubplots("##analyse", 1, 2, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Spectre du mouvement (FFT)")) {
            ImPlot::SetupAxes("f (Hz)", "amplitude relative", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, (kShownModes + 0.5) * f1, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.05, ImPlotCond_Once);
            if (!fs.empty()) ImPlot::PlotLine("FFT de la masse observée", fs.data(), normalized.data(), static_cast<int>(fs.size()), lineSpec(kComputed));
            ImPlot::PlotInfLines("modes de la chaîne", chainF.data(), shown, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue, 0.7f)));
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Dispersion : fréquence du mode n")) {
            ImPlot::SetupAxes("n", "f_n (Hz)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("corde continue n c / 2L", allN.data(), allCont.data(), beads_, lineSpec(kBlue));
            ImPlot::PlotLine("chaîne de N masses", allN.data(), allChain.data(), beads_, lineSpec(kComputed));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void StringModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const float width = 6.0f;                          // largeur de la corde dessinée
    const float sx = width / static_cast<float>(length_), sy = 2.0f * sx;  // hauteurs exagérées d'un facteur 2
    const double a = problem_.spacing();
    auto X = [&](double x) { return static_cast<float>((x - 0.5 * length_) * sx); };

    std::vector<Vertex> lines;
    addSegment(lines, -0.5f * width, -0.5f, 0.0f, -0.5f * width, 0.5f, 0.0f, kPostColor);   // poteaux
    addSegment(lines, 0.5f * width, -0.5f, 0.0f, 0.5f * width, 0.5f, 0.0f, kPostColor);
    addSegment(lines, -0.5f * width, 0.0f, 0.0f, 0.5f * width, 0.0f, 0.0f, kRest);          // ligne de repos
    if (startMode_ > 0) {  // noeuds du mode pur : x = k L / n
        for (int k = 1; k < startMode_; ++k) {
            const float x = X(k * length_ / startMode_);
            addSegment(lines, x, -0.25f, 0.0f, x, 0.25f, 0.0f, kProbeColor);
        }
    }
    renderer.draw(Primitive::Lines, lines);

    auto polyline = [&](const std::vector<double>& u, float z, const float* color) {
        std::vector<Vertex> v;
        v.push_back({X(0.0), 0.0f, z, color[0], color[1], color[2]});
        for (int j = 0; j < beads_; ++j) v.push_back({X((j + 1) * a), sy * static_cast<float>(u[j]), z, color[0], color[1], color[2]});
        v.push_back({X(length_), 0.0f, z, color[0], color[1], color[2]});
        return v;
    };
    renderer.draw(Primitive::LineStrip, polyline(u0_, 0.0f, kInitial));  // forme de départ, en gris
    if (diverged_) return;

    if (atLeast(ctx.level, Level::Lycee))  // solution exacte (modes), bleue, un peu en retrait
        renderer.draw(Primitive::LineStrip, polyline(std::vector<double>(exact_.begin(), exact_.begin() + beads_), -0.04f, kBlue));
    renderer.draw(Primitive::LineStrip, polyline(std::vector<double>(y_.begin(), y_.begin() + beads_), 0.0f, kComputed));

    if (atLeast(ctx.level, Level::Etudiant)) {  // les masses de la chaîne
        std::vector<Vertex> points;
        for (int j = 0; j < beads_; ++j) points.push_back({X((j + 1) * a), sy * static_cast<float>(y_[j]), 0.0f, kComputed[0], kComputed[1], kComputed[2]});
        renderer.draw(Primitive::Points, points, 6.0f);
    }
    const int p = probeBead();  // la masse observée (courbe u(t)), en orange
    renderer.draw(Primitive::Points, {Vertex{X((p + 1) * a), sy * static_cast<float>(y_[p]), 0.0f, kProbeColor[0], kProbeColor[1], kProbeColor[2]}}, 9.0f);
}

void StringModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 7.5f;
    camera.yaw = 0.30f;
    camera.pitch = 0.22f;
}

// ------------------------------- textes --------------------------------

const char* StringModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Pincez une corde de guitare et lâchez-la : elle vibre et le son naît. Plus la corde est tendue, fine ou courte, plus la note est aiguë.\n\n"
                   "Déplacez le point de pincement : la même corde ne sonne pas pareil, car elle ne vibre pas d'une seule façon mais de plusieurs à la fois.\n\n"
                   "Le point orange est la petite masse dont on suit le mouvement dans la courbe du bas. La forme grise est la position de départ.";
        case Level::Interesse:
            return "Une corde peut vibrer de plusieurs façons simples, appelées modes : un arc (mode 1), deux arcs (mode 2), trois arcs... "
                   "Le mode 2 vibre deux fois plus vite que le mode 1, le mode 3 trois fois plus : ce sont les harmoniques.\n\n"
                   "Essayez « Un seul mode » : la corde dessine un motif régulier, avec des points qui ne bougent jamais (les noeuds, en orange).\n\n"
                   "Pincée à la main, la corde mélange tous les modes. Au milieu, il manque les modes pairs (2, 4, 6) ; au tiers, il manque les modes 3, 6, 9 : "
                   "essayez les boutons L/2, L/3... Une corde plus tendue ou plus légère vibre plus vite, une corde plus longue plus lentement.";
        case Level::College:
            return "Vitesse de l'onde sur la corde : c = √(T / μ) (T : tension en newtons, μ : masse par mètre en kg/m).\n"
                   "Le mode 1 est un arc de longueur d'onde λ1 = 2L ; sa fréquence est f1 = c / (2L) (en hertz, vibrations par seconde).\n"
                   "Le mode n a n arcs : λ_n = 2L / n et f_n = n × f1. Ce sont les harmoniques : 1, 2, 3, 4... fois la note de base.\n\n"
                   "La courbe « Part de l'énergie » montre combien chaque mode contient. Un pincement à L/3 n'envoie rien dans les modes 3, 6, 9...\n\n"
                   "L'ordinateur remplace la corde par 59 petites masses reliées par des ressorts et calcule le mouvement de chacune.";
        case Level::Lycee:
            return "Équation d'onde : ∂²u/∂t² = c² ∂²u/∂x², avec u(0, t) = u(L, t) = 0 (bouts fixes). Séparation des variables : u = sin(nπx / L) cos(ω_n t), "
                   "ω_n = nπc / L, soit f_n = n c / (2L). Toute vibration est une somme de modes (superposition).\n\n"
                   "Pincement en x0 (triangle de hauteur h) : u(x, 0) = Σ b_n sin(nπx / L) avec b_n = 2hL² sin(nπ x0 / L) / (n² π² x0 (L − x0)). "
                   "Si x0 = L / m, les modes n multiples de m disparaissent (sin(nπ / m) = 0).\n\n"
                   "Énergie : E = ½ Σ m v² + ½ Σ k (Δu)². Sans frottement elle est constante ; celle du mode n est proportionnelle à n² b_n² (points bleus). "
                   "La courbe bleue du mouvement est la solution exacte par modes.";
        case Level::Etudiant:
            return "Modèle : N masses m = μ a reliées par des ressorts k = T / a (a = L / (N+1)) : m u_j'' = k (u_{j+1} − 2 u_j + u_{j−1}). "
                   "C'est un ressort-masse de M2, N fois. Modes : u_j = sin(nπj / (N+1)) cos(ω_n t), ω_n = 2 √(k/m) sin(nπ / (2(N+1))).\n\n"
                   "La chaîne est DISPERSIVE : f_n < n f1 (graphique « Dispersion »), les harmoniques élevés sont un peu trop graves ; la corde continue est la limite N → ∞.\n\n"
                   "Les solveurs de la Mécanique intègrent la chaîne, comparés à la solution exacte par modes. Au-dessus d'un pas limite ils divergent : "
                   "Verlet et Euler symplectique pour dt < a / c, RK4 pour dt < √2 a / c, Euler explicite jamais stable. Le spectre (FFT) du mouvement d'une masse montre des raies aux f_n.";
        case Level::Chercheur:
            return "Dispersion de la chaîne : ω_n / (nπc / L) = sin(t) / t, t = nπ / (2(N+1)), soit 1 − (nπ)² / (24 (N+1)²) : le mode n est ralenti de ce rapport. Coupure ω_max = 2 √(k/m) = 2c / a.\n\n"
                   "Énergie par mode, exacte : E = Σ_n (m (N+1) / 4)(q_n'² + ω_n² q_n²), q_n = (2 / (N+1)) Σ_j u_j sin(nπj / (N+1)) (orthogonalité des sinus, transformée en sinus discrète).\n\n"
                   "Stabilité : le pire mode (ω_max) impose ω_max dt < 2 pour Verlet, ce qui donne C = c dt / a ≤ 1 : le Verlet de la chaîne EST le schéma saute-mouton de l'équation d'onde d'O0 (écart 4,6e-16). "
                   "RK4 : |ω dt| ≤ 2√2 sur l'axe imaginaire ; Euler explicite : |1 + iω dt| > 1 pour tout dt.\n\n"
                   "Fréquences lues par FFT (fenêtre de Hann, interpolation parabolique de ln|X|) : accord avec ω_n / 2π à 2e-6 (solution exacte) et 1e-5 (RK4). Les harmoniques absents d'un pincement à L/m sont exactement nuls si N+1 est multiple de m.";
    }
    return "";
}

}  // namespace pl
