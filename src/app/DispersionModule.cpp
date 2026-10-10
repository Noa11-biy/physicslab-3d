#include "DispersionModule.hpp"

#include <algorithm>
#include <cmath>

#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/render/FieldPlot.hpp"
#include "physicslab/waves/Dispersion.hpp"

namespace pl {
namespace {

namespace dsp = waves::dispersion;

constexpr double kWavelength = 0.8;      // longueur d'onde de la porteuse du paquet, et échelle des « cases par longueur d'onde » [m]
constexpr double kSpeed = 1.0;           // c [m/s]
constexpr double kLength1D = 16.0;       // 20 longueurs d'onde
constexpr double kSide2D = 6.4;          // 8 longueurs d'onde
constexpr double kPulseSigma = 0.4;      // largeur de la bosse 1D [m]
constexpr double kPulseCenter = 8.0;
constexpr double kPacketStart = 6.0;     // centre initial du paquet [m]
constexpr double kPacketWidth = 1.0;     // largeur de son enveloppe [m]
constexpr double kSigma2D = 0.25;        // largeur de la bosse 2D [m]
constexpr double kEndTime1D = 6.5;
constexpr double kEndTime2D = 4.0;
constexpr float kHeight1D = 1.8f;
constexpr float kHeight2D = 1.4f;
constexpr float kScale1D = 6.0f / 16.0f; // unités de scène par mètre (1D)
constexpr float kScale2D = 0.9f;         // (2D)
constexpr double kColorScale = 0.4;
constexpr float kComputed[3] = {0.35f, 0.85f, 0.45f};   // vert : le calcul
constexpr float kOrange[3] = {1.00f, 0.70f, 0.25f};
constexpr float kSpongeColor[3] = {0.95f, 0.78f, 0.30f};
constexpr float kWall[3] = {0.92f, 0.92f, 0.95f};

double gaussian(double x, double center, double sigma) {
    const double z = (x - center) / sigma;
    return std::exp(-0.5 * z * z);
}

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

}  // namespace

DispersionModule::DispersionModule() { reset(); }

double DispersionModule::cflLimit() const { return dsp::cflLimit(dimension_); }

double DispersionModule::carrierPhase() const { return 2.0 * constants::pi / ppw_; }

// Onde exacte du continu. Bosse lâchée au repos : deux moitiés. Paquet : une porteuse qui va vers la droite sous une enveloppe gaussienne.
double DispersionModule::exact1D(double x, double t) const {
    if (shape_ == kPulse) return 0.5 * (gaussian(x - kSpeed * t, kPulseCenter, kPulseSigma) + gaussian(x + kSpeed * t, kPulseCenter, kPulseSigma));
    const double k0 = 2.0 * constants::pi / kWavelength, y = x - kSpeed * t;
    return gaussian(y, kPacketStart, kPacketWidth) * std::cos(k0 * (y - kPacketStart));
}

// ------------------------------ simulation -----------------------------

void DispersionModule::reset() {
    const double dx = kWavelength / ppw_;
    if (dimension_ == 1) {
        waves::Wave1DParams p;
        p.cells = 20 * ppw_;
        p.length = kLength1D;
        p.speed = kSpeed;
        p.cfl = cfl1D_;
        p.left = p.right = waves::Edge::Absorbing;
        p.spongeCells = 3 * ppw_;  // 2,4 m
        wave1_ = std::make_unique<waves::Wave1D>(p);
        wave2_.reset();
        if (shape_ == kPulse) {
            wave1_->setInitial([](double x) { return gaussian(x, kPulseCenter, kPulseSigma); });
        } else {
            // paquet purement progressif, cohérent avec la relation de dispersion de la grille (sinon : onde de retour parasite, voir Wave1D::setStates)
            const double k0 = 2.0 * constants::pi / kWavelength, x = k0 * dx;
            const double vg = dsp::groupRatio1D(x, cfl1D_) * kSpeed;
            const double arg = cfl1D_ * std::sin(0.5 * x);
            const double omegaDt = arg <= 1.0 ? 2.0 * std::asin(arg) : x * cfl1D_;  // au-delà de la limite l'onde n'existe pas : on garde le continu
            const double dt = wave1_->dt();
            const double groupShift = std::isnan(vg) ? kSpeed * dt : vg * dt;
            wave1_->setStates(
                [=](double xx) { return gaussian(xx + groupShift, kPacketStart, kPacketWidth) * std::cos(k0 * (xx - kPacketStart) + omegaDt); },
                [=](double xx) { return gaussian(xx, kPacketStart, kPacketWidth) * std::cos(k0 * (xx - kPacketStart)); });
        }
        endTime_ = kEndTime1D;
    } else {
        waves::Wave2DParams p;
        p.cellsX = p.cellsY = 8 * ppw_;
        p.dx = dx;
        p.speed = kSpeed;
        p.cfl = cfl2D_;
        p.left = p.right = p.bottom = p.top = waves::Edge::Absorbing;
        p.spongeCells = 2 * ppw_;  // 1,6 m
        wave2_ = std::make_unique<waves::Wave2D>(p);
        wave1_.reset();
        const double half = 0.5 * kSide2D;
        wave2_->setInitial([half](double x, double y) {
            const double r2 = (x - half) * (x - half) + (y - half) * (y - half);
            return std::exp(-0.5 * r2 / (kSigma2D * kSigma2D));
        });
        endTime_ = kEndTime2D;
    }
    accumulator_ = 0.0;
    finished_ = false;
    diverged_ = false;
    peaks_.clear();
}

void DispersionModule::update(double frameSeconds) {
    if (!running_ || finished_ || diverged_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    const double dt = dtNow();
    int guard = 0;
    while (accumulator_ >= dt && guard++ < 400 && !diverged_) {
        if (dimension_ == 1) wave1_->step(); else wave2_->step();
        accumulator_ -= dt;
        const double m = maxAbsNow();
        if (!(m < 1e6)) {  // la grille explose (C trop grand) : on arrête, les valeurs ne se dessinent plus
            diverged_ = true;
            break;
        }
        if (peaks_.size() < 20000) peaks_.push_back(m);
        // un calcul instable continue jusqu'à l'explosion (ou 4000 pas) : c'est la croissance qu'on veut mesurer
        if (timeNow() >= endTime_ && (cflNow() <= cflLimit() + 1e-12 || peaks_.size() >= 4000)) {
            finished_ = true;
            break;
        }
    }
}

// Centre d'énergie (somme x u² / somme u²) du calcul et de la solution exacte, dans la même fenêtre : leur différence est le retard de la grille.
double DispersionModule::lag1D() const {
    const double lo = (shape_ == kPulse ? kPulseCenter : kPacketStart) - 0.5, hi = kLength1D - 2.4;  // hors de l'éponge
    double sn = 0.0, sxn = 0.0, se = 0.0, sxe = 0.0;
    for (int i = 0; i < wave1_->points(); ++i) {
        const double x = wave1_->x(i);
        if (x < lo || x > hi) continue;
        const double un = wave1_->u()[i], ue = exact1D(x, wave1_->time());
        sn += un * un;
        sxn += x * un * un;
        se += ue * ue;
        sxe += x * ue * ue;
    }
    if (sn <= 1e-12 || se <= 1e-12) return std::nan("");
    return sxe / se - sxn / sn;
}

// Front de l'anneau (maximum de u) le long de l'axe x et de la diagonale, avec interpolation parabolique : distance au centre [m].
void DispersionModule::frontRadii2D(double& axis, double& diagonal) const {
    axis = diagonal = 0.0;
    const int n = wave2_->pointsX(), c = n / 2;
    auto crest = [&](int stepI, int stepJ, double unit) {
        double best = 0.0;
        int bestM = 0;
        for (int m = 2; c + stepI * m + 1 < n && c + stepJ * m + 1 < n; ++m) {
            const double v = wave2_->u()[wave2_->index(c + stepI * m, c + stepJ * m)];
            if (v > best) { best = v; bestM = m; }
        }
        if (best < 0.01 || bestM < 2) return 0.0;
        const double a = wave2_->u()[wave2_->index(c + stepI * (bestM - 1), c + stepJ * (bestM - 1))], b = best,
                     d = wave2_->u()[wave2_->index(c + stepI * (bestM + 1), c + stepJ * (bestM + 1))];
        const double denom = a - 2.0 * b + d;
        const double delta = denom != 0.0 ? 0.5 * (a - d) / denom : 0.0;
        return (bestM + delta) * unit;
    };
    axis = crest(1, 0, wave2_->dx());
    diagonal = crest(1, 1, wave2_->dx() * std::sqrt(2.0));
}

// Erreur max d'une bosse gaussienne (largeur 0,4 m) contre d'Alembert, domaine fermé de 8 m, au temps `time`.
double DispersionModule::pulseError(int cellsPerWavelength, double cfl, double time) const {
    waves::Wave1DParams p;
    p.cells = 10 * cellsPerWavelength;
    p.length = 8.0;
    p.speed = kSpeed;
    p.cfl = cfl;
    waves::Wave1D w(p);
    const auto g = [](double x) { return gaussian(x, 4.0, kPulseSigma); };
    w.setInitial(g);
    w.advance(std::lround(time / w.dt()));
    double err = 0.0;
    for (int i = 0; i < w.points(); ++i) {
        const double ct = kSpeed * w.time();
        err = std::max(err, std::abs(w.u()[i] - 0.5 * (g(w.x(i) - ct) + g(w.x(i) + ct))));
    }
    return err;
}

void DispersionModule::computeConvergence() {
    conv_ = Convergence();
    const double cfl = cfl1D_, time = 1.2;
    convCfl_ = cfl;
    convDirty_ = false;
    if (cfl > 1.0 + 1e-12) return;  // instable : pas d'étude
    conv_.exactScheme = cfl >= 1.0 - 1e-9;
    for (int ppw : {6, 12, 24, 48}) {
        const double e = pulseError(ppw, cfl, time);
        conv_.dx.push_back(kWavelength / ppw);
        conv_.err.push_back(std::max(e, 1e-17));
    }
    // pente log-log (moindres carrés) là où l'erreur dépasse l'arrondi
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int m = 0;
    for (std::size_t i = 0; i < conv_.dx.size(); ++i) {
        if (conv_.err[i] < 1e-13) continue;
        const double x = std::log(conv_.dx[i]), y = std::log(conv_.err[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
        ++m;
    }
    if (m >= 2) conv_.slope = (m * sxy - sx * sy) / (m * sxx - sx * sx);
    for (int i = 1; i <= 20; ++i) {
        const double c = 0.05 * i;
        conv_.cfl.push_back(c);
        conv_.errC.push_back(std::max(pulseError(24, c, time), 1e-17));
    }
}

// --------------------------------- UI ----------------------------------

void DispersionModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (!interesse && (dimension_ != 1 || shape_ != kPulse)) {  // niveau 1 : une bosse sur une corde
        dimension_ = 1;
        shape_ = kPulse;
        changed = true;
    }

    if (interesse) {
        ImGui::SeparatorText("Forme");
        int d = dimension_;
        changed |= ImGui::RadioButton("Corde (1D)", &d, 1);
        ImGui::SameLine();
        changed |= ImGui::RadioButton("Surface (2D)", &d, 2);
        dimension_ = d;
        if (dimension_ == 2 && ppw_ > 32) {  // en 2D le nombre de cases croît comme ppw² : bornons à 256 × 256
            ppw_ = 32;
            changed = true;
        }
        if (dimension_ == 1) {
            int s = shape_;
            changed |= ImGui::RadioButton(college ? "Bosse (toutes les longueurs d'onde)" : "Une bosse", &s, kPulse);
            changed |= ImGui::RadioButton(college ? "Paquet d'ondes (une longueur d'onde)" : "Un train d'ondes", &s, kPacket);
            shape_ = s;
        }
    }

    pushSliderWidth();
    ImGui::SeparatorText("Calcul");
    changed |= ImGui::SliderInt(college ? "Cases / longueur d'onde" : "Finesse des cases", &ppw_, 3, dimension_ == 1 ? 60 : 32);
    double& cfl = dimension_ == 1 ? cfl1D_ : cfl2D_;
    if (interesse) changed |= sliderD(college ? "Nombre de Courant C" : "Pas de temps", &cfl, 0.05, 1.2, college ? "%.3f" : "");
    popSliderWidth();
    if (lycee) {
        if (ImGui::Button(dimension_ == 1 ? "C = 1" : "C = 1/√2")) {
            cfl = cflLimit();
            changed = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("C = 0,5")) {
            cfl = 0.5;
            changed = true;
        }
    }
    if (interesse && cfl > cflLimit() + 1e-12) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
        if (college) ImGui::TextWrapped("C = %.3f dépasse la limite %.4f : l'onde saute plus d'une case par pas, le calcul explose.", cfl, cflLimit());
        else ImGui::TextWrapped("Pas de temps trop grand : le calcul va exploser.");
        ImGui::PopStyleColor();
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
        ImGui::TextWrapped("Le calcul a explosé : le pas de temps est trop grand pour la taille des cases (voir « Stabilité » dans Invariants au niveau 6).");
        ImGui::PopStyleColor();
    }

    if (college) {
        ImGui::SeparatorText("Grandeurs");
        ImGui::Text("t = %.2f s / %.2f s", timeNow(), endTime_);
        ImGui::Text("Case dx = %.4f m   pas dt = %.5f s", dxNow(), dtNow());
        ImGui::Text("Cases par longueur d'onde : %d", ppw_);
        if (lycee) ImGui::Text("Phase par case k dx = 2π / %d = %.3f rad", ppw_, carrierPhase());
        if (chercheur) ImGui::Text("Limite CFL : %.4f", cflLimit());
    }
}

void DispersionModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool college = atLeast(level, Level::College);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);

    if (diverged_) {
        wrapped("Calcul divergent : les valeurs sont hors d'échelle. Baisser C, puis « Recommencer ».", true);
    } else if (!college) {
        wrapped(shape_ == kPulse && dimension_ == 1
                    ? "Comparez la courbe verte (l'ordinateur) à la bleue (la nature) : avec des cases grossières, la bosse calculée est un peu déformée et laisse de petits creux derrière elle."
                    : "Comparez ce que fait l'ordinateur à ce que ferait la vraie vague.",
                true);
    } else if (dimension_ == 1) {
        double err = 0.0;
        for (int i = 0; i < wave1_->points(); ++i) err = std::max(err, std::abs(wave1_->u()[i] - exact1D(wave1_->x(i), wave1_->time())));
        ImGui::Text("Écart max |u - exacte| = %.3e", err);
        const double lag = lag1D();
        if (!std::isnan(lag)) ImGui::Text("Retard de la vague calculée : %.2f cm", 100.0 * lag);
        if (etudiant && shape_ == kPacket) {
            const double x = carrierPhase();
            const double vg = dsp::groupRatio1D(x, cfl1D_), vp = dsp::phaseRatio1D(x, cfl1D_);
            ImGui::SeparatorText("Paquet d'ondes");
            ImGui::Text("v_phase / c = %.5f   v_groupe / c = %.5f", vp, vg);
            if (!std::isnan(lag) && wave1_->time() > 0.5) {
                ImGui::Text("retard prévu c t (1 − v_g/c) = %.2f cm", 100.0 * kSpeed * wave1_->time() * (1.0 - vg));
            }
            wrapped("L'énergie voyage à la vitesse de groupe : c'est elle qui fixe le retard.", true);
        }
    } else {
        double axis, diag;
        frontRadii2D(axis, diag);
        // La crête d'une bosse 2D n'est pas exactement en c t même dans le continu (queue de l'onde) : on ne compare donc pas à c t,
        // seulement les deux directions entre elles. Dans le continu elles seraient égales.
        ImGui::Text("Cercle bleu : c t = %.3f m", kSpeed * wave2_->time());
        if (axis > 0.0 && diag > 0.0) {
            ImGui::Text("crête le long d'un axe : %.3f m", axis);
            ImGui::Text("crête en diagonale : %.3f m", diag);
            ImGui::Text("écart diagonale / axe : %+.2f %%", 100.0 * (diag / axis - 1.0));
            if (etudiant) wrapped("Dans le continu ces deux rayons seraient égaux : l'écart mesure l'anisotropie de la grille.", true);
        } else {
            wrapped("(front hors du domaine ou trop faible)", true);
        }
    }
    if (chercheur) {
        const double pi = constants::pi, c = cflNow();
        const double predicted = dimension_ == 1 ? dsp::growthPerStep1D(pi, c) : dsp::growthPerStep2D(pi, pi, c);
        ImGui::SeparatorText("Stabilité");
        ImGui::Text("C = %.3f, limite %.4f (marge %.1f %%)", c, cflLimit(), 100.0 * (1.0 - c / cflLimit()));
        ImGui::Text("croissance par pas du pire mode (prévue) : %.4f", predicted);
        // mesure : moyenne géométrique de max|u| sur les derniers pas, une fois l'instabilité sortie du bruit d'arrondi (valeurs > 1e-3)
        if (predicted > 1.0 && peaks_.size() > 20 && peaks_.back() > 1e-3) {
            const std::size_t n = peaks_.size(), last = std::min<std::size_t>(10, n - 1);
            ImGui::Text("croissance mesurée sur les %d derniers pas : %.4f", static_cast<int>(last), std::pow(peaks_[n - 1] / peaks_[n - 1 - last], 1.0 / static_cast<double>(last)));
            wrapped(dimension_ == 1 ? "Prévu : (C + √(C² − 1))² pour le mode de Nyquist. La mesure est un peu plus faible ici : l'éponge amortit les ondes courtes et le "
                                      "mode dominant n'est pas exactement celui de Nyquist (avec des bords fixes et un mode exact, l'accord est de 6 chiffres : voir test_waves)."
                                    : "Prévu : plus grande racine de λ + 1/λ = 2 − 4 C² (sin² + sin²), phases égales à π. La mesure est plus faible : éponge et modes voisins.", true);
        }
    }
}

void DispersionModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showTheory = atLeast(level, Level::College);
    const int cols = showTheory ? 2 : 1;
    const double pi = constants::pi;
    if (diverged_) return;

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (dimension_ == 1) {
            if (ImPlot::BeginPlot("Profil u(x) : la nature et le calcul")) {
                ImPlot::SetupAxes("x (m)", "u", 0, 0);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, kLength1D, ImPlotCond_Once);
                ImPlot::SetupAxisLimits(ImAxis_Y1, -0.5, 1.1, ImPlotCond_Once);
                const int n = wave1_->points();
                std::vector<double> xs(n), un(n), ue(n);
                for (int i = 0; i < n; ++i) {
                    xs[i] = wave1_->x(i);
                    un[i] = wave1_->u()[i];
                    ue[i] = exact1D(xs[i], wave1_->time());
                }
                ImPlot::PlotLine("Exacte (nature)", xs.data(), ue.data(), n, lineSpec(kBlue));
                ImPlot::PlotLine("Calcul (grille)", xs.data(), un.data(), n, lineSpec(kComputed));
                ImPlot::EndPlot();
            }
        } else if (ImPlot::BeginPlot("Coupes de l'anneau : axe et diagonale")) {
            ImPlot::SetupAxes("r (m)", "u", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, 0.5 * kSide2D, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.3, 0.7, ImPlotCond_Once);
            const int n = wave2_->pointsX(), c = n / 2;
            std::vector<double> ra, ua, rd, ud;
            for (int m = 0; c + m < n; ++m) {
                ra.push_back(m * wave2_->dx());
                ua.push_back(wave2_->u()[wave2_->index(c + m, c)]);
                rd.push_back(m * wave2_->dx() * std::sqrt(2.0));
                ud.push_back(wave2_->u()[wave2_->index(c + m, c + m)]);
            }
            ImPlot::PlotLine("le long d'un axe", ra.data(), ua.data(), static_cast<int>(ra.size()), lineSpec(kComputed));
            ImPlot::PlotLine("en diagonale", rd.data(), ud.data(), static_cast<int>(rd.size()), lineSpec(kOrange));
            const double ct = kSpeed * wave2_->time();
            ImPlot::PlotInfLines("front attendu c t", &ct, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue)));
            ImPlot::EndPlot();
        }

        if (showTheory) {
            const double c = cflNow();
            if (dimension_ == 1 && ImPlot::BeginPlot("Vitesses de la grille / c")) {
                ImPlot::SetupAxes("k dx (rad par case)", "v / c", 0, 0);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, pi, ImPlotCond_Once);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.1, ImPlotCond_Once);
                const int n = 160;
                std::vector<double> xs(n), phase(n), group(n), chain(n);
                for (int i = 0; i < n; ++i) {
                    const double x = pi * (i + 1) / n;
                    xs[i] = x;
                    phase[i] = dsp::phaseRatio1D(x, c);
                    group[i] = dsp::groupRatio1D(x, c);
                    chain[i] = dsp::phaseRatio1D(x, 0.0);
                }
                ImPlot::PlotLine("phase (C actuel)", xs.data(), phase.data(), n, lineSpec(kComputed));
                ImPlot::PlotLine("groupe (C actuel)", xs.data(), group.data(), n, lineSpec(kOrange));
                if (atLeast(level, Level::Lycee)) ImPlot::PlotLine("phase, espace seul (C → 0)", xs.data(), chain.data(), n, lineSpec(kBlue));
                const double k0 = carrierPhase();
                ImPlot::PlotInfLines(shape_ == kPacket ? "porteuse du paquet" : "porteuse de référence", &k0, 1,
                                     ImPlotSpec(ImPlotProp_LineColor, ImVec4(0.8f, 0.8f, 0.8f, 0.8f)));
                ImPlot::EndPlot();
            } else if (dimension_ == 2 && ImPlot::BeginPlot("Vitesse de phase selon la direction")) {
                ImPlot::SetupAxes("angle θ du vecteur d'onde (°)", "v_phase / c", 0, 0);
                ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, 90.0, ImPlotCond_Once);
                ImPlot::SetupAxisLimits(ImAxis_Y1, 0.8, 1.05, ImPlotCond_Once);
                const int n = 91;
                const double x = carrierPhase();
                std::vector<double> th(n), v(n);
                for (int i = 0; i < n; ++i) {
                    th[i] = i;
                    v[i] = dsp::phaseRatio2D(x, i * pi / 180.0, c);
                }
                ImPlot::PlotLine("longueur d'onde 0,8 m", th.data(), v.data(), n, lineSpec(kComputed));
                const double one = 1.0;
                ImPlot::PlotInfLines("vraie vitesse c", &one, 1, ImPlotSpec(ImPlotProp_LineColor, toImVec4(kBlue), ImPlotProp_Flags, ImPlotInfLinesFlags_Horizontal));
                ImPlot::EndPlot();
            }
        }
        ImPlot::EndSubplots();
    }
}

void DispersionModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    if (convDirty_ || std::abs(convCfl_ - cfl1D_) > 1e-12) computeConvergence();

    wrapped("Étude en 1D : erreur d'une bosse (largeur 0,4 m) contre la solution exacte, à t = 1,2 s.", true);
    if (conv_.dx.empty()) {
        wrapped("C > 1 en 1D : le calcul est instable, il n'y a pas d'erreur de dispersion à mesurer (voir Invariants, niveau 6).");
        return;
    }
    if (ImPlot::BeginSubplots("##analyse", 1, 2, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot("Erreur selon la taille des cases")) {
            ImPlot::SetupAxes("dx (m)", "erreur max", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (!conv_.dx.empty()) {
                ImPlot::PlotLine(strf("calcul, C = %.2f", cfl1D_).c_str(), conv_.dx.data(), conv_.err.data(), static_cast<int>(conv_.dx.size()),
                                 ImPlotSpec(ImPlotProp_LineColor, toImVec4(kComputed), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 5.0f));
                if (!conv_.exactScheme) {  // droite de pente 2 ancrée sur le premier point
                    std::vector<double> ref(conv_.dx.size());
                    for (std::size_t i = 0; i < ref.size(); ++i) ref[i] = conv_.err[0] * std::pow(conv_.dx[i] / conv_.dx[0], 2.0);
                    ImPlot::PlotLine("pente 2", conv_.dx.data(), ref.data(), static_cast<int>(ref.size()), lineSpec(kBlue));
                }
            }
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Erreur selon C (24 cases par longueur d'onde)")) {
            ImPlot::SetupAxes("C", "erreur max", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            if (!conv_.cfl.empty()) {
                ImPlot::PlotLine("calcul", conv_.cfl.data(), conv_.errC.data(), static_cast<int>(conv_.cfl.size()),
                                 ImPlotSpec(ImPlotProp_LineColor, toImVec4(kComputed), ImPlotProp_LineWeight, 2.0f, ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 4.0f));
                if (atLeast(level, Level::Chercheur)) {  // prévu : proportionnelle à (1 - C²), ancrée à C = 0,5
                    std::vector<double> ref(conv_.cfl.size());
                    const double anchor = conv_.errC[9] / 0.75;  // C = 0,5 est le 10e point
                    for (std::size_t i = 0; i < ref.size(); ++i) ref[i] = std::max(anchor * (1.0 - conv_.cfl[i] * conv_.cfl[i]), 1e-17);
                    ImPlot::PlotLine("∝ 1 − C²", conv_.cfl.data(), ref.data(), static_cast<int>(ref.size()), lineSpec(kBlue));
                }
            }
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void DispersionModule::drawScene(Renderer& renderer, const UiContext&) {
    std::vector<Vertex> lines;
    if (dimension_ == 1) {
        const float half = 0.5f * static_cast<float>(kLength1D) * kScale1D;
        const float rest[3] = {0.35f, 0.38f, 0.42f};
        addSegment(lines, -half, 0.0f, 0.0f, half, 0.0f, 0.0f, rest);
        addSegment(lines, -half, -0.6f, 0.0f, -half, 1.0f, 0.0f, kSpongeColor);
        addSegment(lines, half, -0.6f, 0.0f, half, 1.0f, 0.0f, kSpongeColor);
        const float layer = static_cast<float>(3 * ppw_ * wave1_->dx()) * kScale1D;  // éponge
        addSegment(lines, -half + layer, -0.25f, 0.0f, -half + layer, 0.25f, 0.0f, kSpongeColor);
        addSegment(lines, half - layer, -0.25f, 0.0f, half - layer, 0.25f, 0.0f, kSpongeColor);
        renderer.draw(Primitive::Lines, lines);
        if (diverged_) return;

        std::vector<Vertex> exact, numeric, curtain;
        const int n = wave1_->points();
        for (int i = 0; i < n; ++i) {
            const double x = wave1_->x(i);
            const float X = static_cast<float>(x - 0.5 * kLength1D) * kScale1D;
            exact.push_back({X, kHeight1D * static_cast<float>(exact1D(x, wave1_->time())), -0.04f, kBlue[0], kBlue[1], kBlue[2]});
            float rgb[3];
            divergingColor(wave1_->u()[i] / kColorScale, rgb);
            numeric.push_back({X, kHeight1D * static_cast<float>(wave1_->u()[i]), 0.0f, kComputed[0], kComputed[1], kComputed[2]});
            if (i % std::max(1, n / 200) == 0) {
                for (float& comp : rgb) comp *= 0.5f;
                addSegment(curtain, X, 0.0f, 0.0f, X, kHeight1D * static_cast<float>(wave1_->u()[i]), 0.0f, rgb);
            }
        }
        renderer.draw(Primitive::Lines, curtain);
        renderer.draw(Primitive::LineStrip, exact);
        renderer.draw(Primitive::LineStrip, numeric);
        return;
    }

    const float half = 0.5f * static_cast<float>(kSide2D) * kScale2D;
    addSegment(lines, -half, 0, -half, half, 0, -half, kWall);
    addSegment(lines, half, 0, -half, half, 0, half, kWall);
    addSegment(lines, half, 0, half, -half, 0, half, kWall);
    addSegment(lines, -half, 0, half, -half, 0, -half, kWall);
    const float layer = static_cast<float>(2 * ppw_ * wave2_->dx()) * kScale2D;
    const float inner = half - layer;
    addSegment(lines, -inner, 0, -inner, inner, 0, -inner, kSpongeColor);
    addSegment(lines, inner, 0, -inner, inner, 0, inner, kSpongeColor);
    addSegment(lines, inner, 0, inner, -inner, 0, inner, kSpongeColor);
    addSegment(lines, -inner, 0, inner, -inner, 0, -inner, kSpongeColor);
    // le front exact, un cercle de rayon c t, en bleu : l'anneau calculé s'en écarte, plus en diagonale que sur un axe
    const float r = static_cast<float>(kSpeed * wave2_->time()) * kScale2D;
    if (r > 0.01f && r < half) {
        const int seg = 96;
        for (int i = 0; i < seg; ++i) {
            const float a0 = 2.0f * static_cast<float>(constants::pi) * i / seg, a1 = 2.0f * static_cast<float>(constants::pi) * (i + 1) / seg;
            addSegment(lines, r * std::cos(a0), 0.0f, r * std::sin(a0), r * std::cos(a1), 0.0f, r * std::sin(a1), kBlue);
        }
    }
    renderer.draw(Primitive::Lines, lines);
    if (diverged_) return;
    const int stride = std::max(1, wave2_->pointsX() / 64);
    renderer.draw(Primitive::Lines, makeRelief(wave2_->u(), wave2_->pointsX(), wave2_->pointsY(), wave2_->dx() * kScale2D, kHeight2D,
                                               static_cast<float>(kColorScale), stride));
}

void DispersionModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 8.5f;
    camera.yaw = 0.35f;
    camera.pitch = 0.45f;
}

// ------------------------------- textes --------------------------------

const char* DispersionModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Un ordinateur ne connaît pas l'eau « continue » : il la découpe en petites cases et avance par petits sauts de temps. "
                   "La vague qu'il calcule est donc un peu fausse.\n\n"
                   "Courbe bleue : ce que fait la vraie nature. Courbe verte : ce que calcule l'ordinateur.\n\n"
                   "Avec des cases grossières (curseur tout à gauche), la bosse calculée se déforme un peu, arrive en retard et laisse de petits creux derrière elle : "
                   "l'effet est discret, regardez bien le graphique. Plus les cases sont fines, plus les deux courbes se confondent.";
        case Level::Interesse:
            return "Ce sont surtout les petites vagues (les « ondes courtes ») que la grille déforme : elles sont ralenties, alors que les grandes vagues passent presque sans erreur. "
                   "Une bosse étroite contient beaucoup d'ondes courtes ; un « train d'ondes » n'en contient qu'une sorte.\n\n"
                   "Le curseur « pas de temps » est étonnant : avec un pas de temps bien choisi, le calcul devient exact sur la corde. Mais un pas trop grand fait exploser la simulation.\n\n"
                   "Sur la surface (2D), la grille a des directions privilégiées : l'anneau calculé n'est pas un cercle parfait (le cercle bleu).";
        case Level::College:
            return "On appelle nombre de Courant C = c × dt / dx : c'est la distance parcourue par l'onde en un pas de temps, mesurée en cases. "
                   "Si C dépasse 1 (limite CFL), l'onde devrait « sauter » plus d'une case par pas : le calcul explose.\n\n"
                   "La longueur d'onde est λ ; il y en a λ / dx cases (curseur « cases par longueur d'onde »). Moins il y en a, plus l'onde est ralentie : "
                   "avec 6 cases par longueur d'onde, la vitesse calculée est de l'ordre de 3 % inférieure à c.\n\n"
                   "Deux vitesses différentes : celle des crêtes (vitesse de phase) et celle du paquet d'ondes tout entier (vitesse de groupe). Le graphique de droite les donne.";
        case Level::Lycee:
            return "Une onde plane exp(i (k x − ω t)) est solution exacte de la grille si sin(ω dt / 2) = C sin(k dx / 2).\n"
                   "Vitesse de phase : v = ω / k. Vitesse de groupe : dω / dk = c cos(k dx / 2) / √(1 − C² sin²(k dx / 2)).\n\n"
                   "Quand les cases tendent vers 0 : ω = c k (pas de dispersion). Sinon v < c et v dépend de k : les ondes courtes sont plus lentes, "
                   "la bosse (somme d'ondes de tous les k) se déforme : c'est la dispersion numérique.\n\n"
                   "À C = 1 en 1D, v_phase = v_groupe = c : le calcul est exact (l'erreur d'espace et l'erreur de temps se compensent). "
                   "L'erreur d'une bosse est proportionnelle à (1 − C²) dx² : divisée par 4 quand les cases sont deux fois plus petites (ordre 2).";
        case Level::Etudiant:
            return "Développement : v_phase / c = 1 − (1 − C²)(k dx)² / 24 + O(k⁴ dx⁴) ; sans le pas de temps (C → 0, la chaîne de masses d'O1) : sin(k dx / 2) / (k dx / 2). "
                   "L'erreur est du second ordre, de coefficient ∝ (1 − C²) : mesuré, les erreurs à C = 0,2 / 0,9 / 0,99 valent 1,28 / 0,25 / 0,03 fois celle à C = 0,5, comme (1 − C²) / 0,75.\n\n"
                   "Un paquet d'ondes avance à la vitesse de groupe : retard prévu c t (1 − v_g / c), mesuré dans l'onglet Invariants.\n\n"
                   "2D : sin²(ω dt / 2) = C² [sin²(kx dx / 2) + sin²(ky dx / 2)] : v dépend de la direction. L'écart au continu est (k dx)² (cos⁴θ + sin⁴θ − C²) / 24 : "
                   "plus faible en diagonale que sur un axe, nul en diagonale à C = 1/√2.\n\n"
                   "Instabilité : au-delà de C = 1 (1D) ou 1/√2 (2D) les modes de la grille proches de Nyquist grandissent à chaque pas.";
        case Level::Chercheur:
            return "Stabilité de von Neumann : pour chaque mode, T^{m+1} = (2 − 4 C² s²) T^m − T^{m−1}, s² = sin²(k dx / 2) (1D) ou la somme des deux sin² (2D). "
                   "Les racines λ + 1/λ = 2 − 4 C² s² sont de module 1 si C s ≤ 1. Sinon la racine de plus grand module |b| + √(b² − 1), b = 1 − 2 C² s², "
                   "amplifie le mode à chaque pas : au pire mode (s² = 1 en 1D, 2 en 2D), (C + √(C² − 1))² en 1D. Mesuré sur la grille : 1,442964 à C = 1,02 (1D, mode n = 190), "
                   "3,433275 à C = 1,2, 2,450844 en 2D à C = 0,78 : égal à la formule à 6 chiffres.\n\n"
                   "Condition CFL : le domaine de dépendance numérique doit contenir le domaine physique. En 1D : c dt ≤ dx. En 2D, un pas ne porte l'information que d'une case en x ou en y : "
                   "le domaine numérique est un losange |x| + |y| ≤ n dx ; le cercle physique de rayon c n dt y tient si c dt ≤ dx / √2 (cercle inscrit) : c'est exactement la limite de von Neumann.\n\n"
                   "Vitesse de groupe mesurée sur un paquet (6,3 cases par longueur d'onde, 6 s) : égale à la formule à 1e-15 près à C = 1 et à 2,5e-5 pour C ≤ 0,9 "
                   "(à condition de démarrer par un état cohérent avec la grille : la vitesse initiale du continu ajoute une onde de retour de 8 %).";
    }
    return "";
}

}  // namespace pl
