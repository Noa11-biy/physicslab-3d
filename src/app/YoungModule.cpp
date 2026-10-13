#include "YoungModule.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/render/FieldPlot.hpp"
#include "physicslab/waves/Fft.hpp"

namespace pl {
namespace {

namespace dif = waves::diffraction;

constexpr double kPi = constants::pi;
constexpr float kScene = 0.5f;                  // unités de scène par mètre
constexpr float kCpuColor[3] = {0.35f, 0.85f, 0.45f};
constexpr float kGuideGrey[3] = {0.55f, 0.57f, 0.62f};
constexpr float kGuideWhite[3] = {0.95f, 0.95f, 0.98f};
constexpr float kScreenColor[3] = {0.95f, 0.70f, 0.25f};   // écran percé de fentes
constexpr float kSpongeColor[3] = {0.55f, 0.50f, 0.25f};
constexpr float kObsColor[3] = {0.45f, 0.90f, 0.85f};      // écran d'observation

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

}  // namespace

YoungModule::YoungModule(std::string shaderDir) : shaderDir_(std::move(shaderDir)) {
    gpuOk_ = gpu_.init(shaderDir_);
    if (!gpuOk_) backend_ = kCpu;
    reset();
}

YoungModule::~YoungModule() { gpu_.shutdown(); }

// ------------------------------ simulation -----------------------------

void YoungModule::reset() {
    if (!gpuOk_ && backend_ != kCpu) backend_ = kCpu;
    // un groupe de fentes doit tenir entre les éponges (19 m utiles)
    width_ = std::clamp(width_, 0.5, 4.0);
    if (slits_ > 1) pitch_ = std::clamp(pitch_, width_ + 0.5, std::max(width_ + 0.5, (kHeight - 2.0 * kSponge - 0.5 - width_) / (slits_ - 1)));

    spec_ = dif::SlitScreenSpec();
    spec_.dx = kDx;
    spec_.cfl = kCfl;
    spec_.lengthX = kWidth;
    spec_.lengthY = kHeight;
    spec_.sourceX = kSourceX;
    spec_.screenX = kScreenX;
    spec_.sponge = kSponge;
    spec_.slits = slits_;
    spec_.slitWidth = width_;
    spec_.pitch = pitch_;
    spec_.lambda = lambda();
    spec_.amplitude = 6.0;   // onde plane d'amplitude ~ 0,5 (A / (2 w) fois 1,06 sur la grille)
    screen_ = dif::makeSlitScreen(spec_);

    cpu_ = std::make_unique<waves::Wave2D>(screen_.params);
    cpu_->addLineSource(screen_.source);
    if (gpuOk_ && backend_ != kCpu) gpu_.upload(*cpu_);

    obsColumn_ = static_cast<int>(std::lround(obsX() / kDx));
    accumulator_ = 0.0;
    finished_ = false;
    snapshots_.clear();
    profileY_.clear();
    profileI_.clear();
    minimaY_.clear();
    measuredFringe_ = 0.0;
    gridValid_ = false;
    ray_.reset();
    gpuVsCpuMax_ = gpuVsCpuL2_ = 0.0;

    // courbes théoriques sur l'écran d'observation : 361 points de -9 à 9 m
    theoryY_.resize(361);
    fraunhoferI_.resize(361);
    for (int k = 0; k < 361; ++k) {
        const double y = -kViewHalf + kViewHalf * 2.0 * k / 360.0;
        theoryY_[k] = y;
        fraunhoferI_[k] = dif::slitIntensity(screen_.slitWidth, screen_.pitch, slits_, lambda(), y / std::hypot(kObsDistance, y));
    }
    gridI_.assign(361, 0.0);

    refreshField(true);
}

void YoungModule::refreshField(bool compare) {
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
    shown_ = &cpu_->u();  // « les deux » : on dessine et on mesure le CPU (référence), on compare le GPU
    if (!compare) return;
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

// Colonne d'observation : un instantané tous les quarts de période (les quatre derniers servent à lire l'amplitude).
void YoungModule::takeSnapshot() {
    const int nx = nodesX(), ny = nodesY();
    std::vector<double> column(static_cast<std::size_t>(ny));
    for (int j = 0; j < ny; ++j) column[static_cast<std::size_t>(j)] = field()[static_cast<std::size_t>(j) * nx + obsColumn_];
    snapshots_.push_back(std::move(column));
    if (snapshots_.size() > 4) snapshots_.pop_front();
}

// Pour u = A cos(phase + k pi/2) : u0 - u2 = 2 A cos(phase), u1 - u3 = -2 A sin(phase), donc A = hypot(u0 - u2, u1 - u3) / 2 : exact pour une sinusoïde, sans réglage.
void YoungModule::computeProfile() {
    profileY_.clear();
    profileI_.clear();
    minimaY_.clear();
    measuredFringe_ = 0.0;
    if (snapshots_.size() < 4 || !steady()) return;
    const double yc = screen_.centerY;
    double top = 0.0;
    for (int j = 0; j < nodesY(); ++j) {
        const double y = j * kDx - yc;
        if (std::abs(y) > kViewHalf) continue;
        const double a = 0.5 * std::hypot(snapshots_[0][j] - snapshots_[2][j], snapshots_[1][j] - snapshots_[3][j]);
        profileY_.push_back(y);
        profileI_.push_back(a * a);
        top = std::max(top, a * a);
    }
    if (top <= 0.0) return;
    for (double& v : profileI_) v /= top;

    // minima profonds (< 30 % du maximum) : ceux qui entourent le maximum central donnent l'interfrange (deux fentes) ou la demi-largeur de la tache (une fente)
    for (double m : dif::findMinima(profileY_, profileI_)) {
        const std::size_t k = static_cast<std::size_t>(std::lround((m - profileY_.front()) / kDx));
        if (k < profileI_.size() && profileI_[k] < 0.3) minimaY_.push_back(m);
    }
    double below = -1e9, above = 1e9;
    for (double m : minimaY_) {
        if (m <= 0.0) below = std::max(below, m);
        else above = std::min(above, m);
    }
    if (below > -1e8 && above < 1e8) measuredFringe_ = above - below;
}

void YoungModule::ensureGridTheory() {
    if (gridValid_) return;
    ray_ = std::make_unique<waves::dispersion::RayPhase2D>(2.0 * kPi / (4.0 * quarter_), kCfl, kDx, 91);
    double top = 0.0;
    for (int k = 0; k < 361; ++k) {
        const double a = dif::apertureAmplitudeGrid(screen_.aperture, kDx, *ray_, screen_.centerX, obsX(), screen_.centerY + theoryY_[k]);
        gridI_[k] = a * a;
        top = std::max(top, gridI_[k]);
    }
    if (top > 0.0)
        for (double& v : gridI_) v /= top;
    gridValid_ = true;
}

void YoungModule::update(double frameSeconds) {
    if (!running_ || finished_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    const double dt = cpu_->dt();
    const int cap = backend_ == kCpu ? maxStepsPerFrame_ : std::min(maxStepsPerFrame_, 160);   // une relecture du champ par quart de période
    int n = std::min(cap, static_cast<int>(accumulator_ / dt));
    accumulator_ = std::min(accumulator_ - n * dt, 2.0 * dt * cap);
    if (n <= 0) return;

    const int q = quarter_;
    double frameMs = 0.0;
    while (n > 0) {
        const long long now = backend_ == kGpu ? gpu_.steps() : cpu_->steps();
        const int chunk = static_cast<int>(std::min<long long>(n, q - now % q));  // jusqu'au prochain quart de période
        if (backend_ != kGpu) {
            const auto t0 = std::chrono::steady_clock::now();
            cpu_->advance(chunk);
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            cpuMsPerStep_ = ms / chunk;
            frameMs += ms;
        }
        if (backend_ != kCpu) {
            const double seconds = gpu_.benchmark(chunk);
            gpuMsPerStep_ = 1e3 * seconds / chunk;
            frameMs += 1e3 * seconds;
        }
        n -= chunk;
        if ((now + chunk) % q == 0) {
            refreshField(false);
            takeSnapshot();
        }
    }
    if (frameMs > 18.0) maxStepsPerFrame_ = std::max(2, static_cast<int>(maxStepsPerFrame_ * 0.8));
    else if (frameMs < 8.0) maxStepsPerFrame_ = std::min(400, maxStepsPerFrame_ + 2);

    refreshField(true);
    computeProfile();
    if (timeNow() >= endTime_) finished_ = true;
}

// --------------------------------- UI ----------------------------------

void YoungModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (!etudiant && slits_ > 2) {  // jusqu'au niveau 4 : une ou deux fentes
        slits_ = 2;
        changed = true;
    }

    ImGui::SeparatorText("Fentes");
    if (etudiant) {
        pushSliderWidth();
        changed |= ImGui::SliderInt("Nombre de fentes", &slits_, 1, 5);
        popSliderWidth();
    } else {
        int s = slits_;
        changed |= ImGui::RadioButton("Une fente", &s, 1);
        changed |= ImGui::RadioButton(college ? "Deux fentes (expérience de Young)" : "Deux fentes", &s, 2);
        slits_ = s;
    }
    pushSliderWidth();
    if (interesse) changed |= sliderD(college ? "Largeur b (m)" : "Largeur des fentes", &width_, 0.5, 4.0, college ? "%.1f" : "");
    if (slits_ > 1) changed |= sliderD(college ? "Écartement a (m)" : "Écartement", &pitch_, 1.5, 8.0, college ? "%.1f" : "");
    if (interesse) {
        double lam = lambda();  // la longueur d'onde est un multiple de 0,25 m (la période dure un nombre entier de quarts de 4 pas)
        if (sliderD(college ? "Longueur d'onde λ (m)" : "Taille des vagues", &lam, 0.75, 2.0, college ? "%.2f" : "")) {
            const int q = std::clamp(static_cast<int>(std::lround(lam / 0.25)), 3, 8);
            if (q != quarter_) {
                quarter_ = q;
                changed = true;
            }
        }
    }
    popSliderWidth();
    if (college && interesse && quarter_ == 3) wrapped("λ = 0,75 m : 7,5 cases par longueur d'onde, la grille est grossière (niveau 5).", true);

    if (etudiant) {
        ImGui::SeparatorText("Calcul");
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
            ImGui::Text("dx = %.3f m   dt = %.5f s   C = %.3f", cpu_->dx(), cpu_->dt(), kCfl);
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
    if (interesse) sliderD("Vitesse du temps", &timeScale_, 0.1, 8.0, "x%.2f", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (college) {
        ImGui::SeparatorText("Grandeurs");
        ImGui::Text("t = %.1f s   (mesures prêtes à t = %.0f s)", timeNow(), readyTime());
        ImGui::Text("longueur d'onde λ = %.2f m", lambda());
        ImGui::Text("fentes : b = %.1f m%s", screen_.slitWidth, slits_ > 1 ? strf(", a = %.1f m", screen_.pitch).c_str() : "");
        ImGui::Text("distance à l'écran D = %.0f m", kObsDistance);
        if (slits_ > 1) ImGui::Text("interfrange prévu i = λ D / a = %.2f m", lambda() * kObsDistance / screen_.pitch);
    }
}

void YoungModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);

    if (!college) {
        wrapped("Les vagues passent par les fentes et s'étalent de l'autre côté. Avec deux fentes, les vagues de chaque fente se croisent : là où deux crêtes se rencontrent la mer monte beaucoup, "
                "là où une crête rencontre un creux elle reste calme. Cela dessine des bandes (les franges) sur l'écran d'observation.",
                true);
        if (!steady()) ImGui::TextDisabled("Attendez quelques secondes : les vagues n'ont pas encore atteint l'écran d'observation.");
        return;
    }

    const double lam = lambda(), D = kObsDistance, a = screen_.pitch, b = screen_.slitWidth;
    if (!steady()) {
        ImGui::TextDisabled("Mesure en attente : t = %.0f s sur %.0f s.", timeNow(), readyTime());
    }

    if (slits_ == 1) {
        ImGui::SeparatorText("Diffraction par une fente");
        const double s1 = lam / b;
        if (s1 < 1.0) {
            const double y1 = D * std::tan(std::asin(s1));
            wrapped(strf("premier minimum : sin θ = λ / b = %.3f  (θ = %.1f°)", s1, std::asin(s1) * 180.0 / kPi));
            wrapped(strf("sur l'écran : y = ± %.2f m (tache centrale large de %.2f m)", y1, 2.0 * y1));
            if (measuredFringe_ > 0.0) wrapped(strf("mesuré : tache centrale large de %.2f m", measuredFringe_));
            else if (steady()) ImGui::Text("mesuré : (la tache déborde de l'écran tracé)");
        } else {
            ImGui::Text("λ > b : pas de minimum, les vagues s'étalent dans toutes les directions");
        }
        if (lycee) wrapped("Plus la fente est étroite, plus la tache est large : la diffraction est d'autant plus forte que l'ouverture est petite devant λ.", true);
    } else if (slits_ == 2) {
        ImGui::SeparatorText("Interférences de Young");
        const double i = lam * D / a;
        ImGui::Text("interfrange prévu : i = λ D / a = %.2f m", i);
        if (lycee) {
            const double s = lam / (2.0 * a);
            wrapped(strf("exact (premier minimum sin θ = λ / (2a)) : %.2f m", 2.0 * D * std::tan(std::asin(s))));
        }
        if (measuredFringe_ > 0.0) wrapped(strf("mesuré (entre les deux minima autour du centre) : %.2f m, soit %+.1f %% par rapport à λ D / a", measuredFringe_, 100.0 * (measuredFringe_ / i - 1.0)));
        else if (steady()) ImGui::Text("mesuré : (minima hors de l'écran tracé)");
        if (lycee && lam / b < 1.0) wrapped(strf("la tache d'une fente (enveloppe) s'annule à sin θ = λ / b = %.3f", lam / b));
        if (college && !etudiant) wrapped("Un écart de quelques % avec la formule vient de la grille de calcul, qui n'est pas parfaite : 10 cases seulement par longueur d'onde (niveau 5).", true);
    } else {
        ImGui::SeparatorText("Réseau de fentes");
        wrapped(strf("champ lointain : pics principaux à sin θ = m λ / a (m = 0, ±1, ...), le premier à sin θ = %.3f, soit y = %.2f m", lam / a, D * std::tan(std::asin(std::min(1.0, lam / a)))));
        const double fresnel = std::pow((slits_ - 1) * a + b, 2.0) / (lam * D);
        if (fresnel > 1.0)
            wrapped(strf("Ici le nombre de Fresnel du réseau vaut %.0f : on est en champ PROCHE, les pics ne sont pas encore séparés (on voit une image des fentes qui se brouille). "
                         "La courbe de Fraunhofer ne s'applique pas à cette distance.", fresnel), true);
        wrapped(strf("entre deux pics principaux : %d minima et %d pics secondaires", slits_ - 1, slits_ - 2));
        wrapped(strf("largeur d'un pic : environ λ / (N a) = %.3f en sin θ", lam / (slits_ * a)));
        wrapped("Plus il y a de fentes, plus les pics principaux sont fins et brillants (ils sont N² fois plus intenses en crête).", true);
    }

    if (etudiant) {
        ImGui::SeparatorText("Grille et théorie");
        const double lg = dif::gridWavelength(2.0 * kPi / lam, cpu_->dt(), kDx, kCfl);
        wrapped(strf("λ de la grille %.4f m (%.2f cases) pour λ = %.3f m", lg, lg / kDx, lam));
        if (b / lam < 1.5) wrapped("La fente est étroite (b < 1,5 λ) : l'approximation de Fraunhofer (ouverture éclairée uniformément) est moins bonne.", true);
        if (chercheur) {
            wrapped("La grille est anisotrope : elle décale légèrement les franges loin des fentes (courbe orange = Rayleigh-Sommerfeld avec la phase de la grille).", true);
            wrapped(strf("nombre de Fresnel b² / (λ D) = %.2f   (a² / (λ D) = %.2f)", b * b / (lam * D), a * a / (lam * D)));
        }
    }

    if (etudiant && backend_ != kCpu && gpuOk_) {
        ImGui::SeparatorText("GPU");
        ImGui::Text("GPU : %.3f ms par pas", gpuMsPerStep_);
        if (backend_ == kBoth) {
            ImGui::Text("CPU : %.2f ms par pas (x%.0f plus lent)", cpuMsPerStep_, gpuMsPerStep_ > 0.0 ? cpuMsPerStep_ / gpuMsPerStep_ : 0.0);
            ImGui::Text("écart GPU / CPU : max %.1e, L2 %.1e", gpuVsCpuMax_, gpuVsCpuL2_);
            wrapped("Le GPU calcule en float, le CPU en double : un écart de l'ordre de 1e-6 est l'arrondi. L'écran est un simple coefficient nul sur les faces : le shader n'a rien de nouveau.", true);
        }
    } else if (etudiant) {
        ImGui::Text("CPU : %.2f ms par pas", cpuMsPerStep_);
    }
}

void YoungModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool theory = atLeast(level, Level::College);
    const bool factors = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);

    const int plots = factors ? 2 : 1;
    if (ImPlot::BeginSubplots("##graphes", 1, plots, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Intensité sur l'écran d'observation")) {
            ImPlot::SetupAxes("y - yc (m)", "intensité / max", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, -kViewHalf, kViewHalf, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.05, 1.1, ImPlotCond_Always);
            if (!profileY_.empty())
                ImPlot::PlotLine("calcul (mesure)", profileY_.data(), profileI_.data(), static_cast<int>(profileY_.size()), lineSpec(kCpuColor));
            if (theory) ImPlot::PlotLine("théorie de Fraunhofer", theoryY_.data(), fraunhoferI_.data(), 361, lineSpec(kBlue));
            if (etudiant) {
                ensureGridTheory();
                ImPlot::PlotLine("Rayleigh-Sommerfeld, phase de la grille", theoryY_.data(), gridI_.data(), 361, lineSpec(kScreenColor));
            }
            ImPlot::EndPlot();
        }
        if (factors && ImPlot::BeginPlot("Théorie : fente × interférences, I(sin θ)")) {
            ImPlot::SetupAxes("sin θ", "I / I(0)", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, -1.0, 1.0, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.05, 1.1, ImPlotCond_Always);
            std::vector<double> s(401), single(401), grating(401), product(401);
            const double lam = lambda(), b = screen_.slitWidth, a = screen_.pitch;
            for (int k = 0; k < 401; ++k) {
                s[k] = -1.0 + 2.0 * k / 400.0;
                single[k] = dif::slitIntensity(b, a, 1, lam, s[k]);
                product[k] = dif::slitIntensity(b, a, slits_, lam, s[k]);
                grating[k] = single[k] > 1e-12 ? product[k] / single[k] : dif::slitIntensity(1e-9, a, slits_, lam, s[k]);
            }
            ImPlot::PlotLine("une fente : sinc²(π b sin θ / λ)", s.data(), single.data(), 401, lineSpec(kBlue));
            if (slits_ > 1) {
                ImPlot::PlotLine("interférences des fentes", s.data(), grating.data(), 401, lineSpec(kScreenColor));
                ImPlot::PlotLine("produit : ce qu'on voit", s.data(), product.data(), 401, lineSpec(kCpuColor));
            }
            if (chercheur) {
                // champ lointain par FFT de l'ouverture (courbe discrète, points)
                const dif::FarField ff = dif::fraunhofer(screen_.aperture, kDx, lam, 16);
                ImPlot::PlotLine("FFT de l'ouverture", ff.sinTheta.data(), ff.intensity.data(), static_cast<int>(ff.sinTheta.size()), lineSpec(kGuideWhite));
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void YoungModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const float half = 0.5f * static_cast<float>(kWidth) * kScene;
    const float halfZ = 0.5f * static_cast<float>(kHeight) * kScene;
    auto X = [](double x) { return static_cast<float>((x - 0.5 * kWidth) * kScene); };
    auto Z = [](double y) { return static_cast<float>((y - 0.5 * kHeight) * kScene); };

    std::vector<Vertex> lines;
    addSegment(lines, -half, 0, -halfZ, half, 0, -halfZ, kGuideWhite);  // cadre
    addSegment(lines, half, 0, -halfZ, half, 0, halfZ, kGuideWhite);
    addSegment(lines, half, 0, halfZ, -half, 0, halfZ, kGuideWhite);
    addSegment(lines, -half, 0, halfZ, -half, 0, -halfZ, kGuideWhite);
    const float layer = static_cast<float>(kSponge) * kScene;
    addSegment(lines, -half + layer, 0, -halfZ + layer, half - layer, 0, -halfZ + layer, kSpongeColor);
    addSegment(lines, half - layer, 0, -halfZ + layer, half - layer, 0, halfZ - layer, kSpongeColor);
    addSegment(lines, half - layer, 0, halfZ - layer, -half + layer, 0, halfZ - layer, kSpongeColor);
    addSegment(lines, -half + layer, 0, halfZ - layer, -half + layer, 0, -halfZ + layer, kSpongeColor);

    // la source en ligne
    addSegment(lines, X(kSourceX), 0.03f, -halfZ, X(kSourceX), 0.03f, halfZ, kGuideGrey);

    // l'écran : les tronçons pleins (contour des cases fermées)
    const double x0 = (screen_.screenColumn - 0.5) * kDx, x1 = (screen_.screenColumn + screen_.screenColumns - 0.5) * kDx;
    const int ny = static_cast<int>(screen_.aperture.size());
    int start = -1;
    for (int j = 0; j <= ny; ++j) {
        const bool solid = j < ny && screen_.aperture[static_cast<std::size_t>(j)] < 0.5;
        if (solid && start < 0) start = j;
        if (!solid && start >= 0) {
            const double ya = (start - 0.5) * kDx, yb = (j - 0.5) * kDx;
            addSegment(lines, X(x0), 0.04f, Z(ya), X(x1), 0.04f, Z(ya), kScreenColor);
            addSegment(lines, X(x1), 0.04f, Z(ya), X(x1), 0.04f, Z(yb), kScreenColor);
            addSegment(lines, X(x1), 0.04f, Z(yb), X(x0), 0.04f, Z(yb), kScreenColor);
            addSegment(lines, X(x0), 0.04f, Z(yb), X(x0), 0.04f, Z(ya), kScreenColor);
            start = -1;
        }
    }
    // l'écran d'observation (|y - yc| <= 9 m)
    const double yc = screen_.centerY;
    addSegment(lines, X(obsX()), 0.04f, Z(yc - kViewHalf), X(obsX()), 0.04f, Z(yc + kViewHalf), kObsColor);

    if (atLeast(ctx.level, Level::Lycee)) {  // directions prévues : premier minimum de la fente (bleu) et premier minimum des franges (blanc)
        const double lam = lambda(), len = kObsDistance * 1.05, cx = screen_.centerX;
        auto ray = [&](double s, const float* c) {
            if (s <= 0.0 || s >= 1.0) return;
            const double th = std::asin(s);
            for (double sign : {1.0, -1.0})
                addSegment(lines, X(cx), 0.05f, Z(yc), X(cx + len * std::cos(th)), 0.05f, Z(yc + sign * len * std::sin(th)), c);
        };
        ray(lam / screen_.slitWidth, kBlue);
        if (slits_ > 1) ray(lam / (2.0 * screen_.pitch), kGuideWhite);
    }
    renderer.draw(Primitive::Lines, lines);

    display_.resize(field().size());
    for (std::size_t i = 0; i < display_.size(); ++i) display_[i] = std::clamp(field()[i], -1.2, 1.2);
    const int stride = std::max(1, nodesX() / 90);
    renderer.draw(Primitive::Lines, makeRelief(display_, nodesX(), nodesY(), kDx * kScene, 1.0f, 0.5f, stride));
}

void YoungModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 11.0f;
    camera.yaw = 0.0f;
    camera.pitch = 1.1f;
}

// ------------------------------- textes --------------------------------

const char* YoungModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Des vagues régulières arrivent de la gauche sur une barrière (orange) percée d'une ou deux fentes. De l'autre côté, elles ne continuent pas tout droit : elles s'étalent en demi-cercles. C'est la diffraction.\n\n"
                   "Avec deux fentes, les vagues de chaque fente se croisent. Là où deux crêtes se rencontrent, la mer monte beaucoup. Là où une crête rencontre un creux, elle reste calme. "
                   "Sur le mur d'observation (trait turquoise), cela dessine des bandes régulières, les franges (courbe verte du graphe).\n\n"
                   "Rouge : l'eau monte. Bleu : elle descend. Essayez d'écarter les fentes : les bandes se resserrent.";
        case Level::Interesse:
            return "Une fente étroite étale les vagues dans toutes les directions : plus la fente est petite devant la longueur d'onde, plus l'étalement est fort. "
                   "Une fente large laisse passer un faisceau qui s'élargit à peine.\n\n"
                   "Avec deux fentes, chaque fente est une nouvelle source de vagues. En un point de l'écran d'observation, les deux vagues arrivent avec un retard différent : "
                   "si le retard est un nombre entier de périodes, elles s'ajoutent (frange brillante) ; s'il est d'une demi-période en plus, elles s'annulent (frange sombre).\n\n"
                   "C'est l'expérience que Thomas Young a faite avec de la lumière en 1801 : elle a montré que la lumière est une onde. "
                   "Essayez de changer la taille des vagues (la longueur d'onde) : des vagues plus longues donnent des franges plus larges.";
        case Level::College:
            return "Longueur d'onde λ (en m), écartement des fentes a, largeur b, distance D à l'écran d'observation. "
                   "Deux fentes donnent des franges régulières espacées de i = λ × D / a (interfrange) : doubler a divise i par deux, doubler λ le double. "
                   "Exemple : λ = 1 m, D = 12 m, a = 4 m donne i = 3 m.\n\n"
                   "Une fente de largeur b donne une tache centrale dont le bord (premier minimum) se trouve à l'angle θ tel que sin θ = λ / b. Si b est plus petit que λ, il n'y a plus de minimum : l'onde part dans toutes les directions.\n\n"
                   "Le graphe de gauche compare la mesure faite sur la simulation (vert) à la théorie (bleu). Le panneau Invariants donne l'interfrange prévu et mesuré.";
        case Level::Lycee:
            return "Deux fentes séparées de a : pour la direction θ, la différence de marche entre les deux ondes vaut δ = a sin θ. Interférence constructive si δ = m λ, destructive si δ = (m + 1/2) λ. "
                   "Premier minimum : sin θ = λ / (2a). Sur un écran à la distance D (petits angles, tan θ ≈ sin θ ≈ y / D) : i = λ D / a.\n\n"
                   "Chaque fente de largeur b diffracte : son intensité vaut sinc²(π b sin θ / λ), avec un premier zéro en sin θ = λ / b. L'intensité des deux fentes est le PRODUIT : "
                   "sinc²(π b sin θ / λ) × cos²(π a sin θ / λ). Le graphe de droite montre les deux facteurs et leur produit : la tache de diffraction (bleu) sert d'enveloppe aux franges (orange).\n\n"
                   "Les traits bleus et blancs sur la scène sont les directions des premiers minima.";
        case Level::Etudiant:
            return "Champ lointain (Fraunhofer) = transformée de Fourier de l'ouverture : A(θ) ~ somme des e^{i k y sin θ} sur les fentes. Pour N fentes de largeur b et de pas a : "
                   "I(θ) / I(0) = sinc²(π b sin θ / λ) × [sin(N φ) / (N sin φ)]², φ = π a sin θ / λ. Pics principaux en sin θ = m λ / a, largeur ~ λ / (N a), N − 2 pics secondaires entre eux.\n\n"
                   "Le calcul est la même équation d'onde qu'en O3 : l'écran est un masque de noeuds dont toutes les faces ont un coefficient nul (mur rigide, ∂u/∂n = 0), l'énergie discrète reste conservée à 1e-14 près. "
                   "La source est une LIGNE de sources ponctuelles en phase : onde plane d'amplitude A dx / (2 c² sin k dx) sur la grille (mesuré à 0,06 % près). "
                   "L'intensité se lit par détection synchrone sur la grille, ou ici sur quatre instantanés espacés d'un quart de période : A = √((u0 − u2)² + (u1 − u3)²) / 2.\n\n"
                   "Choisissez le GPU, le CPU, ou les deux : le GPU a besoin du même shader qu'en O3 plus une ligne de source, et reste à 1e-6 du CPU. "
                   "Limite de Fraunhofer : b²/(λ D) petit ; ici D est de l'ordre de 12 λ, donc le champ n'est que « presque » lointain.";
        case Level::Chercheur:
            return "Au-delà de Fraunhofer : intégrale de Rayleigh-Sommerfeld (écran rigide, 2e espèce) u(P) = −(k/2) × somme sur l'ouverture de H0(k r) dy', valable à toute distance (proche ou lointaine). "
                   "Elle redonne Fraunhofer quand r tend vers l'infini et l'amplitude 1 sans écran (la somme de H0 sur un plan vaut 2 e^{i k x} / k).\n\n"
                   "La courbe orange remplace la phase continue k r par celle de la grille : la grille de calcul est anisotrope (le nombre d'onde varie d'environ 1 % avec la direction à 10 cases par λ), "
                   "et sur 12 à 22 longueurs d'onde de trajet cela déphase les deux fentes de quelques % de λ. Mesuré ici (D = 12 m, a = 4 m, λ = 1 m) : interfrange 3,21 m contre 3,02 m pour le continu (+7 %) et 3,20 m avec la phase de la grille. "
                   "Sur un arc de 22 λ (tests/test_waves.cpp) le premier minimum des franges est à sin θ = 0,1338 contre 0,1238 prévu (+8,1 %) à 10 cases par λ, 0,1282 contre 0,1248 (+2,7 %) à 20 cases par λ ; "
                   "avec la phase de la grille l'écart tombe à 0,006 en sin θ au plus. Fente simple b = 2 m : premier minimum à 0,5135 mesuré, 0,4948 continu, 0,5111 grille.\n\n"
                   "Nombre de Fresnel F = b² / (λ D) : F << 1 Fraunhofer, F ~ 1 Fresnel (la tache change de forme avec D). Pour N fentes, la FFT de l'ouverture (courbe blanche) retrouve la formule sinc² × réseau à 8e-4 près.";
    }
    return "";
}

}  // namespace pl
