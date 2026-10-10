#include "WaveTankModule.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/render/FieldPlot.hpp"
#include "physicslab/waves/Dispersion.hpp"

namespace pl {
namespace {

namespace dsp = waves::dispersion;

constexpr double kPi = constants::pi;
constexpr float kScene = 0.5f;                  // unités de scène par mètre
constexpr float kCpuColor[3] = {0.35f, 0.85f, 0.45f};
constexpr float kGuideGrey[3] = {0.55f, 0.57f, 0.62f};
constexpr float kGuideWhite[3] = {0.95f, 0.95f, 0.98f};
constexpr float kMedium[3] = {0.95f, 0.70f, 0.25f};   // contour des milieux
constexpr float kSponge[3] = {0.55f, 0.50f, 0.25f};
constexpr std::size_t kMaxCpuNodes = 110000;

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

// Fresnel d'une onde plane (u et c² du/dx continus) : fraction d'énergie réfléchie, et vrai s'il y a réflexion totale.
double fresnelReflected(double theta1, double c1, double c2, bool& total) {
    const double s2 = c2 / c1 * std::sin(theta1);
    total = s2 >= 1.0;
    if (total) return 1.0;
    const double th2 = std::asin(s2);
    const double r = (c1 * std::cos(theta1) - c2 * std::cos(th2)) / (c1 * std::cos(theta1) + c2 * std::cos(th2));
    return r * r;
}

}  // namespace

WaveTankModule::WaveTankModule(std::string shaderDir) : shaderDir_(std::move(shaderDir)) {
    gpuOk_ = gpu_.init(shaderDir_);
    reset();
}

WaveTankModule::~WaveTankModule() { gpu_.shutdown(); }

double WaveTankModule::cLeft() const { return medium_ == kInterface && swapMedia_ ? ratio_ : 1.0; }
double WaveTankModule::cRight() const { return medium_ == kInterface ? (swapMedia_ ? 1.0 : ratio_) : 1.0; }

double WaveTankModule::speedAt(double x, double y) const {
    switch (medium_) {
        case kInterface: return x < kInterfaceX ? cLeft() : cRight();
        case kLens: return (x - kLensX) * (x - kLensX) + (y - 5.0) * (y - 5.0) < kLensRadius * kLensRadius ? ratio_ : 1.0;
        default: return 1.0;
    }
}

// ------------------------------ simulation -----------------------------

void WaveTankModule::reset() {
    if (!gpuOk_ && backend_ != kCpu) backend_ = kCpu;
    const bool cpuNeeded = backend_ != kGpu;
    if (cpuNeeded) ppw_ = std::min(ppw_, 14);  // 365 x 281 noeuds au plus sur le CPU double

    const double dx = 0.5 / ppw_;
    waves::Wave2DParams p;
    p.cellsX = 26 * ppw_;
    p.cellsY = 20 * ppw_;
    p.dx = dx;
    p.speed = 1.0;
    p.cfl = 0.6;
    p.left = p.right = p.bottom = p.top = waves::Edge::Absorbing;
    p.spongeCells = 2 * ppw_;  // 1 m
    p.speedField = waves::Wave2D::makeSpeedField(p, [this](double x, double y) { return speedAt(x, y); });
    cpu_ = std::make_unique<waves::Wave2D>(p);

    if (source_ == kPoint) {
        waves::PointSource s;
        s.x = sourceX();
        s.y = sourceY();
        s.omega = 2.0 * kPi * freq_;
        s.amplitude = 2.0;
        cpu_->addSource(s);
    } else {
        const double th = angleDeg_ * kPi / 180.0, c1 = cLeft(), k = 2.0 * kPi * freq_ / c1;
        const double x0 = kInterfaceX - kPacketDistance, y0 = 5.0 - kPacketDistance * std::tan(th);
        const double nx = std::cos(th), ny = std::sin(th), dt = cpu_->dt(), sigma = kPacketSigma;
        const auto env = [sigma](double ex, double ey) { return std::exp(-0.5 * (ex * ex + ey * ey) / (sigma * sigma)); };
        // paquet qui avance vers l'interface à la vitesse c1 : u^n = F(r - c n t), donc u^{-1} = F(r + c n dt)
        cpu_->setStates(
            [=](double x, double y) { return env(x - x0 + c1 * nx * dt, y - y0 + c1 * ny * dt) * std::cos(k * ((x - x0) * nx + (y - y0) * ny) + c1 * k * dt); },
            [=](double x, double y) { return env(x - x0, y - y0) * std::cos(k * ((x - x0) * nx + (y - y0) * ny)); });
        hitTime_ = kPacketDistance / (c1 * nx);
    }

    if (gpuOk_ && backend_ != kCpu) gpu_.upload(*cpu_);
    accumulator_ = 0.0;
    finished_ = false;
    leftEnergy_.clear();
    rightEnergy_.clear();
    lastSeriesTime_ = -1.0;
    haveT_ = haveR_ = false;
    measuredT_ = measuredR_ = 0.0;
    gpuVsCpuMax_ = gpuVsCpuL2_ = 0.0;
    ampHistory_.clear();
    amplitude_.clear();
    ampRadii_.clear();
    for (int m = 0; m < 70; ++m) ampRadii_.push_back(0.4 + 3.6 * m / 69.0);

    refreshField();
    double e, cx, cy;
    regionStats(-1.0, 1e9, e, cx, cy);
    totalEnergy0_ = std::max(e, 1e-300);
    sampleSeries();
}

void WaveTankModule::refreshField() {
    if (backend_ == kCpu || !gpuOk_) {
        shown_ = &cpu_->u();
        return;
    }
    gpu_.readField(gpuField_);
    gpuFieldD_.assign(gpuField_.begin(), gpuField_.end());
    if (backend_ == kGpu) {
        shown_ = &gpuFieldD_;
        return;
    }
    shown_ = &cpu_->u();  // « les deux » : on dessine le CPU (référence) et on mesure l'écart du GPU
    double num = 0.0, den = 0.0, maxDiff = 0.0, maxCpu = 0.0;
    for (std::size_t i = 0; i < gpuFieldD_.size(); ++i) {
        const double diff = gpuFieldD_[i] - cpu_->u()[i];
        maxDiff = std::max(maxDiff, std::abs(diff));
        maxCpu = std::max(maxCpu, std::abs(cpu_->u()[i]));
        num += diff * diff;
        den += cpu_->u()[i] * cpu_->u()[i];
    }
    gpuVsCpuMax_ = maxCpu > 0.0 ? maxDiff / maxCpu : 0.0;
    gpuVsCpuL2_ = den > 0.0 ? std::sqrt(num / den) : 0.0;
}

// Σ u² dx² (énergie à un facteur ω² près), et centre de Σ x u², Σ y u², dans la zone xmin <= x < xmax.
void WaveTankModule::regionStats(double xmin, double xmax, double& energy, double& cx, double& cy) const {
    const int nx = nodesX(), ny = nodesY();
    const double dx = cpu_->dx();
    double s = 0.0, sx = 0.0, sy = 0.0;
    const std::vector<double>& u = field();
    for (int j = 0; j < ny; ++j) {
        const double y = j * dx;
        for (int i = 0; i < nx; ++i) {
            const double x = i * dx;
            if (x < xmin || x >= xmax) continue;
            const double u2 = u[static_cast<std::size_t>(j) * nx + i] * u[static_cast<std::size_t>(j) * nx + i];
            s += u2;
            sx += x * u2;
            sy += y * u2;
        }
    }
    energy = s * dx * dx;
    cx = s > 0.0 ? sx / s : 0.0;
    cy = s > 0.0 ? sy / s : 0.0;
}

void WaveTankModule::sampleSeries() {
    if (!packetCrossesInterface()) return;
    const double t = timeNow();
    if (lastSeriesTime_ >= 0.0 && t - lastSeriesTime_ < 0.1) return;
    lastSeriesTime_ = t;
    double el, er, cx, cy;
    regionStats(-1.0, kInterfaceX - 0.3, el, cx, cy);
    regionStats(kInterfaceX + 0.3, 1e9, er, cx, cy);
    leftEnergy_.add(t, el / totalEnergy0_);
    rightEnergy_.add(t, er / totalEnergy0_);
}

// Angles des faisceaux réfléchi et transmis : déplacement du centre de Σ u² entre un instant où la queue du paquet a fini de franchir l'interface
// (4 largeurs / c plus une marge) et l'instant courant (au moins 1 s plus tard). Voir tests/test_waves.cpp (testRefraction) pour la méthode et sa précision.
void WaveTankModule::measureAngles() {
    if (!packetCrossesInterface()) return;
    const double t = timeNow(), c1 = cLeft(), c2 = cRight();
    const double tT = hitTime_ + 4.0 * kPacketSigma / c2 + 0.5, tR = hitTime_ + 4.0 * kPacketSigma / c1 + 0.5;
    double e, cx, cy;
    if (!haveT_ && t >= tT) {
        regionStats(kInterfaceX + 0.05, 1e9, e, cx, cy);
        if (e > 0.02 * totalEnergy0_) { haveT_ = true; txA_ = cx; tyA_ = cy; tA_ = t; }
    } else if (haveT_ && t >= tA_ + 1.0) {
        regionStats(kInterfaceX + 0.05, 1e9, e, cx, cy);
        if (e > 0.02 * totalEnergy0_) { measuredT_ = std::atan2(cy - tyA_, cx - txA_) * 180.0 / kPi; measuredTEnergy_ = e / totalEnergy0_; }
    }
    if (!haveR_ && t >= tR) {
        regionStats(-1.0, kInterfaceX - 0.05, e, cx, cy);
        if (e > 0.02 * totalEnergy0_) { haveR_ = true; rxA_ = cx; ryA_ = cy; rA_ = t; }
    } else if (haveR_ && t >= rA_ + 1.0) {
        regionStats(-1.0, kInterfaceX - 0.05, e, cx, cy);
        if (e > 0.02 * totalEnergy0_) measuredR_ = std::atan2(cy - ryA_, -(cx - rxA_)) * 180.0 / kPi;
    }
}

// Amplitude de la source ponctuelle le long de +x : maximum de |u| sur les deux dernières périodes, en 70 points entre 0,4 et 4 m.
void WaveTankModule::updateAmplitudes() {
    if (source_ != kPoint) return;
    const double dx = cpu_->dx(), t = timeNow();
    const int nx = nodesX(), j = static_cast<int>(std::lround(sourceY() / dx)), i0 = static_cast<int>(std::lround(sourceX() / dx));
    std::vector<float> row(ampRadii_.size());
    for (std::size_t m = 0; m < ampRadii_.size(); ++m) {
        const int i = std::min(i0 + static_cast<int>(std::lround(ampRadii_[m] / dx)), nx - 1);
        row[m] = static_cast<float>(std::abs(field()[static_cast<std::size_t>(j) * nx + i]));
    }
    ampHistory_.emplace_back(t, std::move(row));
    while (!ampHistory_.empty() && ampHistory_.front().first < t - 2.0 / freq_) ampHistory_.erase(ampHistory_.begin());
    amplitude_.assign(ampRadii_.size(), 0.0);
    for (const auto& entry : ampHistory_)
        for (std::size_t m = 0; m < amplitude_.size(); ++m) amplitude_[m] = std::max(amplitude_[m], static_cast<double>(entry.second[m]));
}

void WaveTankModule::update(double frameSeconds) {
    if (!running_ || finished_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    const double dt = cpu_->dt();
    const int n = std::min(maxStepsPerFrame_, static_cast<int>(accumulator_ / dt));
    accumulator_ = std::min(accumulator_ - n * dt, 2.0 * dt * maxStepsPerFrame_);  // pas de « spirale de la mort » : on abandonne le retard
    if (n <= 0) return;

    double frameMs = 0.0;
    if (backend_ != kGpu) {
        const auto t0 = std::chrono::steady_clock::now();
        cpu_->advance(n);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        cpuMsPerStep_ = ms / n;
        frameMs += ms;
    }
    if (backend_ != kCpu) {
        const double seconds = gpu_.benchmark(n);
        gpuMsPerStep_ = 1e3 * seconds / n;
        frameMs += 1e3 * seconds;
    }
    // cadence : on vise 8 à 18 ms de calcul par image pour garder l'interface fluide (le CPU double est lent en build Debug)
    if (frameMs > 18.0) maxStepsPerFrame_ = std::max(1, static_cast<int>(maxStepsPerFrame_ * 0.8));
    else if (frameMs < 8.0) maxStepsPerFrame_ = std::min(400, maxStepsPerFrame_ + 2);

    refreshField();
    sampleSeries();
    measureAngles();
    updateAmplitudes();
    if (timeNow() >= endTime_) finished_ = true;
}

// --------------------------------- UI ----------------------------------

void WaveTankModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (!interesse && source_ != kPoint) {  // niveau 1 : des rides circulaires
        source_ = kPoint;
        changed = true;
    }

    ImGui::SeparatorText("Cuve");
    int m = medium_;
    changed |= ImGui::RadioButton("Eau partout pareille", &m, kUniform);
    changed |= ImGui::RadioButton(college ? "Interface : zone plus lente à droite" : "Une zone peu profonde", &m, kInterface);
    changed |= ImGui::RadioButton(college ? "Lentille circulaire plus lente" : "Une bosse au fond (lentille)", &m, kLens);
    medium_ = m;

    pushSliderWidth();
    changed |= sliderD(college ? "Fréquence f (Hz)" : "Rythme de la source", &freq_, 0.8, 4.0, college ? "%.2f" : "", ImGuiSliderFlags_Logarithmic);
    if (medium_ != kUniform)
        changed |= sliderD(college ? "Vitesse lente / rapide" : "Profondeur de la zone", &ratio_, 0.35, 1.0, college ? "%.2f" : "");
    popSliderWidth();

    if (interesse) {
        ImGui::SeparatorText("Source");
        int s = source_;
        changed |= ImGui::RadioButton("Une goutte (rides circulaires)", &s, kPoint);
        changed |= ImGui::RadioButton("Un train d'ondes plan", &s, kPacket);
        source_ = s;
        if (source_ == kPacket) {
            pushSliderWidth();
            changed |= sliderD(college ? "Angle d'incidence θ1 (°)" : "Inclinaison", &angleDeg_, 0.0, 55.0, college ? "%.0f" : "");
            popSliderWidth();
            if (medium_ != kInterface) wrapped("Le paquet plan se mesure surtout avec « zone plus lente » (interface plane).", true);
        }
        if (lycee && medium_ == kInterface) changed |= ImGui::Checkbox("Zone lente à GAUCHE (réflexion totale possible)", &swapMedia_);
    }

    if (etudiant) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        const int maxPpw = (backend_ == kGpu && gpuOk_) ? 24 : 14;
        changed |= ImGui::SliderInt("Finesse (cases / 0,5 m)", &ppw_, 8, maxPpw);
        popSliderWidth();
        if (gpuOk_) {
            ImGui::TextDisabled("Processeur de calcul");
            int b = backend_;
            changed |= ImGui::RadioButton("CPU (double)", &b, kCpu);
            changed |= ImGui::RadioButton("GPU (float)", &b, kGpu);
            changed |= ImGui::RadioButton("Les deux (le CPU vérifie le GPU)", &b, kBoth);
            backend_ = b;
        } else {
            wrapped("Pas de compute shader disponible : calcul sur le CPU.", true);
        }
        if (chercheur) {
            ImGui::Text("Grille %d x %d = %d noeuds", nodesX(), nodesY(), nodesX() * nodesY());
            ImGui::Text("dx = %.4f m   dt = %.5f s", cpu_->dx(), cpu_->dt());
        }
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

    if (college) {
        ImGui::SeparatorText("Grandeurs");
        ImGui::Text("t = %.2f s", timeNow());
        const double c1 = cLeft(), c2 = cRight();
        ImGui::Text("vitesse : %.2f m/s", c1);
        if (medium_ != kUniform) ImGui::Text("vitesse dans la zone lente : %.2f m/s", std::min(c1, c2) < 1.0 ? std::min(c1, c2) : ratio_);
        ImGui::Text("longueur d'onde λ = c / f = %.3f m", c1 / freq_);
        if (medium_ != kUniform) ImGui::Text("indice n = c_rapide / c_lent = %.2f", 1.0 / ratio_);
    }
}

void WaveTankModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);

    if (!college) {
        wrapped("L'eau ne se déplace pas : seule la forme de la vague avance. La vague accélère ou ralentit selon la profondeur, et tourne quand elle change de zone.", true);
        return;
    }

    double maxU = 0.0;
    for (double v : field()) maxU = std::max(maxU, std::abs(v));
    ImGui::Text("Hauteur maximale |u| = %.3f", maxU);

    if (packetCrossesInterface()) {
        const double th1 = angleDeg_ * kPi / 180.0, c1 = cLeft(), c2 = cRight();
        bool total = false;
        const double r2 = fresnelReflected(th1, c1, c2, total);
        const double s2 = c2 / c1 * std::sin(th1);
        ImGui::SeparatorText("Loi de Snell-Descartes");
        ImGui::Text("angle d'incidence θ1 = %.1f°", angleDeg_);
        const double lastLeft = leftEnergy_.size() > 0 ? leftEnergy_.y[(leftEnergy_.offset + leftEnergy_.size() - 1) % leftEnergy_.size()] : 0.0;
        if (total) ImGui::Text("sin θ2 = %.2f > 1 : réflexion TOTALE", s2);
        else ImGui::Text("θ2 prévu : %.1f°   mesuré : %s", std::asin(s2) * 180.0 / kPi, measuredT_ != 0.0 ? strf("%.1f°", measuredT_).c_str() : "(en attente)");
        // un faisceau réfléchi de quelques % d'énergie n'a pas de centre fiable (voir tests) : on ne l'affiche que s'il porte au moins 10 %
        ImGui::Text("angle réfléchi : prévu %.1f°   mesuré : %s", angleDeg_,
                    measuredR_ != 0.0 && lastLeft > 0.1 ? strf("%.1f°", measuredR_).c_str() : (haveR_ ? "(trop faible)" : "(en attente)"));
        if (!total && !etudiant) wrapped("Un écart de 1 à 2° vient de la grille de calcul, qui n'est pas parfaite (niveau 5).", true);
        if (c2 > c1) ImGui::Text("angle limite θc = asin(c_lent / c_rapide) = %.1f°", std::asin(c1 / c2) * 180.0 / kPi);
        if (lycee) {
            ImGui::SeparatorText("Énergie");
            ImGui::Text("réfléchie (onde plane) : %.1f %%   transmise : %.1f %%", 100.0 * r2, 100.0 * (1.0 - r2));
            const double l = leftEnergy_.size() > 0 ? leftEnergy_.y[(leftEnergy_.offset + leftEnergy_.size() - 1) % leftEnergy_.size()] : 0.0;
            const double r = rightEnergy_.size() > 0 ? rightEnergy_.y[(rightEnergy_.offset + rightEnergy_.size() - 1) % rightEnergy_.size()] : 0.0;
            ImGui::Text("mesurée : à gauche %.1f %%   à droite %.1f %%", 100.0 * l, 100.0 * r);
            if (etudiant) wrapped("Un paquet a une largeur angulaire d'environ ±10° : près de l'angle limite une partie des composantes se transmet, d'où quelques % de plus que l'onde plane.", true);
        }
        if (etudiant && !total) {
            // loi de Snell de la GRILLE : w et ky conservés, vitesse de groupe discrète
            const double dx = cpu_->dx(), dt = cpu_->dt(), k = 2.0 * kPi * freq_ / c1;
            const double a = k * std::cos(th1) * dx, b = k * std::sin(th1) * dx, cfl1 = c1 * dt / dx, cfl2 = c2 * dt / dx;
            const double sb = std::sin(0.5 * b), w = cfl1 * std::sqrt(std::sin(0.5 * a) * std::sin(0.5 * a) + sb * sb);
            const double sx2sq = (w / cfl2) * (w / cfl2) - sb * sb;
            if (sx2sq > 0.0 && sx2sq < 1.0) {
                double vx, vy;
                dsp::groupVelocity2D(2.0 * std::asin(std::sqrt(sx2sq)), b, cfl2, vx, vy);
                ImGui::Text("θ2 du faisceau sur la grille : %.2f° (%.1f cases / λ)", std::atan2(vy, vx) * 180.0 / kPi, c1 / freq_ / dx);
            }
        }
    } else if (source_ == kPoint) {
        ImGui::SeparatorText("Rides circulaires");
        ImGui::Text("λ = %.3f m   (%.1f cases par longueur d'onde)", cLeft() / freq_, cLeft() / freq_ / cpu_->dx());
        if (!amplitude_.empty() && medium_ == kUniform) {
            const double a1 = amplitude_[10], a2 = amplitude_[60];
            if (a1 > 0.0 && a2 > 0.0)
                ImGui::Text("pente de A(r) entre %.1f et %.1f m : %.3f (attendu -0.5)", ampRadii_[10], ampRadii_[60], std::log(a2 / a1) / std::log(ampRadii_[60] / ampRadii_[10]));
        }
    }

    if (etudiant && backend_ != kCpu && gpuOk_) {
        ImGui::SeparatorText("GPU");
        ImGui::Text("GPU : %.3f ms par pas", gpuMsPerStep_);
        if (backend_ == kBoth) {
            ImGui::Text("CPU : %.2f ms par pas (x%.0f plus lent)", cpuMsPerStep_, gpuMsPerStep_ > 0.0 ? cpuMsPerStep_ / gpuMsPerStep_ : 0.0);
            ImGui::Text("écart GPU / CPU : max %.1e, L2 %.1e", gpuVsCpuMax_, gpuVsCpuL2_);
            wrapped("Le GPU calcule en float, le CPU en double : un écart de l'ordre de 1e-6 est l'arrondi.", true);
        }
    } else if (etudiant) {
        ImGui::Text("CPU : %.2f ms par pas", cpuMsPerStep_);
    }
}

void WaveTankModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showTheory = atLeast(level, Level::Lycee);
    const double dx = cpu_->dx();
    const int nx = nodesX(), j = static_cast<int>(std::lround(5.0 / dx));

    if (ImPlot::BeginSubplots("##graphes", 1, 2, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Coupe de la cuve au milieu : u(x)")) {
            ImPlot::SetupAxes("x (m)", "u", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, kWidth, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, source_ == kPoint ? -0.4 : -1.1, source_ == kPoint ? 0.4 : 1.1, ImPlotCond_Once);
            std::vector<double> xs(nx), us(nx);
            for (int i = 0; i < nx; ++i) {
                xs[i] = i * dx;
                us[i] = field()[static_cast<std::size_t>(j) * nx + i];
            }
            ImPlot::PlotLine("u(x, y = 5 m)", xs.data(), us.data(), nx, lineSpec(kCpuColor));
            if (medium_ == kInterface) {
                const double xi = kInterfaceX;
                ImPlot::PlotInfLines("interface", &xi, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kMedium)));
            }
            ImPlot::EndPlot();
        }
        if (source_ == kPoint) {
            if (ImPlot::BeginPlot("Amplitude des rides A(r)")) {
                ImPlot::SetupAxes("r (m)", "A", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                if (!amplitude_.empty()) {
                    ImPlot::PlotLine("calcul", ampRadii_.data(), amplitude_.data(), static_cast<int>(ampRadii_.size()), lineSpec(kCpuColor));
                    if (showTheory && medium_ == kUniform && amplitude_[10] > 0.0) {
                        // 1 / sqrt(r) et |H0(k r)|, ancrés sur le calcul à r = ampRadii_[10]
                        const double k = 2.0 * kPi * freq_ / cLeft(), r0 = ampRadii_[10];
                        std::vector<double> sq(ampRadii_.size()), hk(ampRadii_.size());
                        const auto h0 = [&](double r) { return std::hypot(std::cyl_bessel_j(0.0, k * r), std::cyl_neumann(0.0, k * r)); };
                        for (std::size_t m = 0; m < sq.size(); ++m) {
                            sq[m] = amplitude_[10] * std::sqrt(r0 / ampRadii_[m]);
                            hk[m] = amplitude_[10] * h0(ampRadii_[m]) / h0(r0);
                        }
                        ImPlot::PlotLine("1 / √r", ampRadii_.data(), sq.data(), static_cast<int>(sq.size()), lineSpec(kBlue));
                        if (atLeast(level, Level::Etudiant)) ImPlot::PlotLine("|H0(k r)| (exacte)", ampRadii_.data(), hk.data(), static_cast<int>(hk.size()), lineSpec(kMedium));
                    }
                }
                ImPlot::EndPlot();
            }
        } else if (ImPlot::BeginPlot("Où va l'énergie du paquet")) {
            ImPlot::SetupAxes("t (s)", "part de l'énergie", ImPlotAxisFlags_AutoFit, 0);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.05, 1.1, ImPlotCond_Always);
            if (leftEnergy_.size() > 0) {
                ImPlot::PlotLine("à gauche de l'interface", leftEnergy_.x.data(), leftEnergy_.y.data(), leftEnergy_.size(), lineSpec(kCpuColor, leftEnergy_.offset));
                ImPlot::PlotLine("à droite de l'interface", rightEnergy_.x.data(), rightEnergy_.y.data(), rightEnergy_.size(), lineSpec(kBlue, rightEnergy_.offset));
            }
            if (showTheory && packetCrossesInterface()) {
                bool total = false;
                const double r2 = fresnelReflected(angleDeg_ * kPi / 180.0, cLeft(), cRight(), total);
                const double tr = 1.0 - r2;
                ImPlot::PlotInfLines("transmise (Fresnel, onde plane)", &tr, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue, 0.6f), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void WaveTankModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const float half = 0.5f * static_cast<float>(kWidth) * kScene;  // demi-largeur (x)
    const float halfZ = 0.5f * static_cast<float>(kHeight) * kScene;
    auto X = [](double x) { return static_cast<float>((x - 0.5 * kWidth) * kScene); };
    auto Z = [](double y) { return static_cast<float>((y - 0.5 * kHeight) * kScene); };

    std::vector<Vertex> lines;
    addSegment(lines, -half, 0, -halfZ, half, 0, -halfZ, kGuideWhite);  // cadre
    addSegment(lines, half, 0, -halfZ, half, 0, halfZ, kGuideWhite);
    addSegment(lines, half, 0, halfZ, -half, 0, halfZ, kGuideWhite);
    addSegment(lines, -half, 0, halfZ, -half, 0, -halfZ, kGuideWhite);
    const float layer = static_cast<float>(1.0) * kScene;  // éponge de 1 m
    addSegment(lines, -half + layer, 0, -halfZ + layer, half - layer, 0, -halfZ + layer, kSponge);
    addSegment(lines, half - layer, 0, -halfZ + layer, half - layer, 0, halfZ - layer, kSponge);
    addSegment(lines, half - layer, 0, halfZ - layer, -half + layer, 0, halfZ - layer, kSponge);
    addSegment(lines, -half + layer, 0, halfZ - layer, -half + layer, 0, -halfZ + layer, kSponge);

    if (medium_ == kInterface) addSegment(lines, X(kInterfaceX), 0.02f, -halfZ, X(kInterfaceX), 0.02f, halfZ, kMedium);
    if (medium_ == kLens) {
        const int seg = 72;
        for (int i = 0; i < seg; ++i) {
            const double a0 = 2.0 * kPi * i / seg, a1 = 2.0 * kPi * (i + 1) / seg;
            addSegment(lines, X(kLensX + kLensRadius * std::cos(a0)), 0.02f, Z(5.0 + kLensRadius * std::sin(a0)), X(kLensX + kLensRadius * std::cos(a1)), 0.02f,
                       Z(5.0 + kLensRadius * std::sin(a1)), kMedium);
        }
    }
    if (source_ == kPoint) {  // la source
        for (int i = 0; i < 12; ++i) {
            const double a0 = 2.0 * kPi * i / 12, a1 = 2.0 * kPi * (i + 1) / 12;
            addSegment(lines, X(sourceX() + 0.12 * std::cos(a0)), 0.03f, Z(sourceY() + 0.12 * std::sin(a0)), X(sourceX() + 0.12 * std::cos(a1)), 0.03f,
                       Z(sourceY() + 0.12 * std::sin(a1)), kGuideWhite);
        }
    }
    if (packetCrossesInterface() && atLeast(ctx.level, Level::Lycee)) {  // rayons : normale, incident, réfléchi et réfracté prévus (loi de Snell)
        const double th1 = angleDeg_ * kPi / 180.0, c1 = cLeft(), c2 = cRight(), len = 2.6;
        const double hx = kInterfaceX, hy = 5.0;
        addSegment(lines, X(hx - len), 0.04f, Z(hy), X(hx + len), 0.04f, Z(hy), kGuideGrey);
        addSegment(lines, X(hx), 0.04f, Z(hy), X(hx - len * std::cos(th1)), 0.04f, Z(hy - len * std::sin(th1)), kGuideWhite);
        addSegment(lines, X(hx), 0.04f, Z(hy), X(hx - len * std::cos(th1)), 0.04f, Z(hy + len * std::sin(th1)), kBlue);
        const double s2 = c2 / c1 * std::sin(th1);
        if (s2 < 1.0) {
            const double th2 = std::asin(s2);
            addSegment(lines, X(hx), 0.04f, Z(hy), X(hx + len * std::cos(th2)), 0.04f, Z(hy + len * std::sin(th2)), kBlue);
        }
    }
    renderer.draw(Primitive::Lines, lines);

    // Le champ est singulier au noeud de la source (une pointe 5 fois plus haute que les rides) : on écrête l'affichage, pas le calcul.
    const bool point = source_ == kPoint;
    const double clampValue = point ? 0.22 : 1.2;
    display_.resize(field().size());
    for (std::size_t i = 0; i < display_.size(); ++i) display_[i] = std::clamp(field()[i], -clampValue, clampValue);
    const int stride = std::max(1, nodesX() / 70);
    renderer.draw(Primitive::Lines, makeRelief(display_, nodesX(), nodesY(), cpu_->dx() * kScene, point ? 4.0f : 1.6f, point ? 0.15f : 0.5f, stride));
}

void WaveTankModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 10.0f;
    camera.yaw = 0.15f;
    camera.pitch = 0.95f;
}

// ------------------------------- textes --------------------------------

const char* WaveTankModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Dans une cuve à ondes, on voit comment une vague voyage. Laissez tomber une goutte : des rides circulaires s'étalent.\n\n"
                   "Rouge : l'eau monte. Bleu : elle descend.\n\n"
                   "Là où l'eau est peu profonde (zone orange), les vagues vont plus lentement. En changeant de zone, la vague change de direction, "
                   "comme une voiture dont une roue entre dans le sable. Une partie de la vague rebondit aussi sur la frontière.\n\n"
                   "Avec une « bosse au fond » en forme de disque, les vagues se rassemblent derrière, comme la lumière dans une loupe.";
        case Level::Interesse:
            return "La vitesse d'une vague dépend de la profondeur : lente en eau peu profonde, rapide en eau profonde. La fréquence (le rythme de la source) ne change pas, "
                   "donc dans la zone lente les vagues sont plus serrées : leur longueur d'onde diminue.\n\n"
                   "Choisissez « un train d'ondes plan » et inclinez-le : à la frontière, une partie rebondit (réflexion), une partie continue en changeant de direction (réfraction). "
                   "Quand la vague entre en zone lente, elle se redresse vers la perpendiculaire à la frontière. Les rayons blancs et bleus du niveau 4 montrent les directions prévues.\n\n"
                   "Autour d'une goutte, la hauteur des rides diminue avec la distance : l'énergie se répartit sur des cercles de plus en plus longs.";
        case Level::College:
            return "Fréquence f (en Hz), vitesse c, longueur d'onde λ : c = λ × f. Passer en zone lente (c plus petit) à f constante divise λ dans le même rapport.\n\n"
                   "Indice du milieu lent : n = c_rapide / c_lent. À la frontière, l'angle θ1 de l'onde incidente et l'angle θ2 de l'onde transmise (mesurés par rapport à la perpendiculaire à la frontière) "
                   "vérifient : sin θ1 / sin θ2 = c1 / c2 (loi de Snell-Descartes, la même que pour la lumière).\n\n"
                   "Réflexion : l'onde rebondit avec le même angle θ1. Le panneau Invariants compare l'angle mesuré dans la simulation à l'angle prévu.";
        case Level::Lycee:
            return "Sur la frontière x = 0 les crêtes doivent se raccorder : la composante k_y du vecteur d'onde est la même des deux côtés, et la pulsation ω aussi. Avec k = ω / c : "
                   "k1 sin θ1 = k2 sin θ2, soit sin θ1 / c1 = sin θ2 / c2 (Snell-Descartes).\n\n"
                   "Si on va de la zone lente vers la zone rapide (case « zone lente à gauche »), sin θ2 = (c2 / c1) sin θ1 peut dépasser 1 : plus d'onde transmise, c'est la réflexion totale, "
                   "au-delà de l'angle limite θc = asin(c_lent / c_rapide).\n\n"
                   "Les amplitudes se déduisent de la continuité de u et de c² ∂u/∂x : R = (c1 cos θ1 − c2 cos θ2) / (c1 cos θ1 + c2 cos θ2), et la part d'énergie réfléchie vaut R². "
                   "Autour d'une source ponctuelle l'amplitude décroît comme 1 / √r (courbe de droite).";
        case Level::Etudiant:
            return "Équation de la cuve à profondeur variable : ∂²u/∂t² = ∇·(c² ∇u) (c² = g h), discrétisée par le schéma saute-mouton d'O0 avec c² moyenné sur les faces entre cases ; "
                   "l'énergie discrète reste conservée à 4e-15 près (sans éponge). Le calcul se fait au choix en double sur le CPU ou en float sur le GPU (compute shader, un thread par case) : "
                   "le GPU est environ 100 fois plus rapide ici et reste à 1e-6 du CPU (mode « les deux »).\n\n"
                   "Une source ponctuelle sinusoïdale produit en 2D l'onde |H0(k r)|, de l'ordre de 1 / √r loin de la source : mesuré à 0,6 % près (pente −0,5017). "
                   "Pour un paquet, la direction du faisceau est celle de la vitesse de GROUPE : sur la grille, elle diffère de celle du vecteur d'onde de plus de 1° à 10 cases par longueur d'onde "
                   "(anisotropie, voir O2). La loi de Snell de la grille — conservation de ω et de k_y avec la relation de dispersion discrète — prévoit l'angle mesuré à 0,1° près, celle du continu se trompe de 1,5°.";
        case Level::Chercheur:
            return "Relation de dispersion de la grille : sin²(ω dt / 2) = C² [sin²(kx dx / 2) + sin²(ky dx / 2)], C = c dt / dx. Dans le milieu 2 : ω et ky conservés donnent sin²(kx2 dx / 2) = (sin(ω dt / 2) / C2)² − sin²(ky dx / 2) ; "
                   "pas de solution : onde évanescente (réflexion totale). Le faisceau suit grad_k ω = c (sx cos(kx dx / 2), sy cos(ky dx / 2)) / (S √(1 − C² S²)).\n\n"
                   "Un paquet de largeur σ a un spectre angulaire d'écart-type 1 / (k σ) (≈ 10° ici) : l'énergie réfléchie est la moyenne de |R(θ)|² sur ce spectre, pondérée par ω² exp(−σ² |k − k0|²). "
                   "Mesuré à θ1 = 30°, c2 / c1 = 1 / 0,7 : 0,0227 contre 0,0247 prévu pour le paquet (0,0051 pour une onde plane). Près de l'angle limite la prédiction perd en précision (1,3 % d'énergie et 0,9° d'angle en réflexion totale).\n\n"
                   "Interface d'une case : c²_face = moyenne arithmétique des c² voisins ; une lentille circulaire est un escalier de cases (rugosité de l'ordre de dx).";
    }
    return "";
}

}  // namespace pl
