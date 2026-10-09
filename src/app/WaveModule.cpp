#include "WaveModule.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <implot.h>

#include "physicslab/core/Constants.hpp"
#include "physicslab/render/FieldPlot.hpp"
#include "physicslab/waves/Fft.hpp"

namespace pl {
namespace {

constexpr double kLength1D = 8.0;   // longueur de la corde [m]
constexpr double kSide2D = 4.0;     // côté de la surface [m]
constexpr float kHeight1D = 1.6f;   // hauteur dessinée pour u = 1 [m de scène]
constexpr float kHeight2D = 1.4f;
constexpr double kColorScale = 0.6; // |u| qui sature la palette
constexpr float kComputed[3] = {0.35f, 0.85f, 0.45f};  // vert : le calcul (comme RK4 dans la Mécanique)
constexpr float kWallColor[3] = {0.92f, 0.92f, 0.95f};
constexpr float kFreeColor[3] = {0.45f, 0.90f, 0.55f};
constexpr float kSpongeColor[3] = {0.95f, 0.78f, 0.30f};

waves::Edge toEdge(int kind) {
    return kind == 0 ? waves::Edge::Absorbing : (kind == 1 ? waves::Edge::Fixed : waves::Edge::Free);
}

const float* edgeColor(int kind) { return kind == 0 ? kSpongeColor : (kind == 1 ? kWallColor : kFreeColor); }

void addSegment(std::vector<Vertex>& v, float x0, float y0, float z0, float x1, float y1, float z1, const float* c) {
    v.push_back({x0, y0, z0, c[0], c[1], c[2]});
    v.push_back({x1, y1, z1, c[0], c[1], c[2]});
}

// Rectangle dans le plan y = 0, demi-côtés hx et hz.
void addFrame(std::vector<Vertex>& v, float hx, float hz, const float* c) {
    addSegment(v, -hx, 0, -hz, hx, 0, -hz, c);
    addSegment(v, hx, 0, -hz, hx, 0, hz, c);
    addSegment(v, hx, 0, hz, -hx, 0, hz, c);
    addSegment(v, -hx, 0, hz, -hx, 0, -hz, c);
}

}  // namespace

WaveModule::WaveModule() { reset(); }

double WaveModule::length() const { return dimension_ == 1 ? kLength1D : kSide2D; }
double WaveModule::cflLimit() const { return dimension_ == 1 ? 1.0 : waves::Wave2D::kCflLimit; }
double WaveModule::dtNow() const { return dimension_ == 1 ? wave1_->dt() : wave2_->dt(); }
double WaveModule::dxNow() const { return dimension_ == 1 ? wave1_->dx() : wave2_->dx(); }
double WaveModule::timeNow() const { return dimension_ == 1 ? wave1_->time() : wave2_->time(); }
double WaveModule::energyNow() const { return dimension_ == 1 ? wave1_->energy() : wave2_->energy(); }
double WaveModule::maxAbsNow() const { return dimension_ == 1 ? wave1_->maxAbs() : wave2_->maxAbs(); }

double WaveModule::bump1D(double x) const {
    const double z = (x - 0.5 * kLength1D) / width_;
    return std::exp(-0.5 * z * z);
}

// Méthode des images : un mur fixe renvoie l'impulsion inversée, un bout libre la renvoie telle quelle, un bord absorbant la laisse
// partir (milieu infini, l'éponge réelle réfléchit un peu : l'écart avec cette courbe mesure ce reste). Le profil prolongé F vérifie
// F(-y) = s F(y) en 0 et F(L + z) = s F(L - z) en L ; u(x, t) = ½ [F(x - ct) + F(x + ct)] (au repos initial).
double WaveModule::exact1D(double x, double t) const {
    const double c = speed_, big = kLength1D;
    auto folded = [&](double y) {
        if (edge_ == kAbsorbing) return bump1D(y);
        const double sign = edge_ == kWall ? -1.0 : 1.0;
        double s = 1.0;
        for (int k = 0; k < 200; ++k) {  // chaque réflexion ramène y vers [0, L]
            if (y < 0.0) { y = -y; s *= sign; }
            else if (y > big) { y = 2.0 * big - y; s *= sign; }
            else break;
        }
        return s * bump1D(y);
    };
    return 0.5 * (folded(x - c * t) + folded(x + c * t));
}

// ------------------------------ simulation -----------------------------

void WaveModule::reset() {
    const waves::Edge edge = toEdge(edge_);
    if (dimension_ == 1) {
        waves::Wave1DParams p;
        p.cells = cells1D_;
        p.length = kLength1D;
        p.speed = speed_;
        p.cfl = cfl1D_;
        p.left = p.right = edge;
        p.spongeCells = std::max(4, static_cast<int>(std::lround(spongeFraction_ * cells1D_)));
        p.spongeReflection = spongeStrength_;
        p.spongeOrder = spongeOrder_;
        wave1_ = std::make_unique<waves::Wave1D>(p);
        wave2_.reset();
        wave1_->setInitial([this](double x) { return bump1D(x); });
    } else {
        waves::Wave2DParams p;
        p.cellsX = p.cellsY = cells2D_;
        p.dx = kSide2D / cells2D_;
        p.speed = speed_;
        p.cfl = cfl2D_;
        p.left = p.right = p.bottom = p.top = edge;
        p.spongeCells = std::max(4, static_cast<int>(std::lround(spongeFraction_ * cells2D_)));
        p.spongeReflection = spongeStrength_;
        p.spongeOrder = spongeOrder_;
        wave2_ = std::make_unique<waves::Wave2D>(p);
        wave1_.reset();
        const double w = width_, half = 0.5 * kSide2D;
        wave2_->setInitial([w, half](double x, double y) {
            const double r2 = (x - half) * (x - half) + (y - half) * (y - half);
            return std::exp(-0.5 * r2 / (w * w));
        });
    }
    initialEnergy_ = std::max(energyNow(), 1e-300);
    // trois traversées du domaine : les murs ramènent l'impulsion au point de départ au bout de 2 L / c
    endTime_ = 3.0 * length() / speed_;
    accumulator_ = 0.0;
    finished_ = false;
    diverged_ = false;
    energy_.clear();
    lastSample_ = 0.0;
    sample();
}

void WaveModule::sample() {
    lastSample_ = timeNow();
    if (!(maxAbsNow() < 1e6)) {  // C trop grand : la grille explose (NaN compris) ; on arrête sans polluer les courbes
        diverged_ = true;
        return;
    }
    energy_.add(lastSample_, energyNow() / initialEnergy_);
}

void WaveModule::update(double frameSeconds) {
    if (!running_ || finished_ || diverged_) return;
    accumulator_ += std::min(frameSeconds, 0.1) * timeScale_;
    const double dt = dtNow();
    int guard = 0;
    while (accumulator_ >= dt && guard++ < 400 && !diverged_) {
        if (dimension_ == 1) wave1_->step(); else wave2_->step();
        accumulator_ -= dt;
        if (timeNow() - lastSample_ >= 1.0 / 60.0) sample();
        if (timeNow() >= endTime_) {
            finished_ = true;
            sample();
            break;
        }
    }
}

std::vector<double> WaveModule::cut() const {
    if (dimension_ == 1) return wave1_->u();
    const int n = wave2_->pointsX(), j = wave2_->pointsY() / 2;
    std::vector<double> row(n);
    for (int i = 0; i < n; ++i) row[i] = wave2_->u()[wave2_->index(i, j)];
    return row;
}

// --------------------------------- UI ----------------------------------

void WaveModule::drawControls(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool interesse = atLeast(level, Level::Interesse);
    const bool college = atLeast(level, Level::College);
    const bool lycee = atLeast(level, Level::Lycee);
    const bool etudiant = atLeast(level, Level::Etudiant);
    const bool chercheur = atLeast(level, Level::Chercheur);
    bool changed = false;

    if (level == Level::Vulgarisation && dimension_ != 2) {  // au niveau 1 on ne montre que la surface de l'eau
        dimension_ = 2;
        changed = true;
    }

    pushSliderWidth();
    ImGui::SeparatorText(dimension_ == 1 ? "Corde" : "Surface de l'eau");
    changed |= sliderD(college ? "Largeur de la bosse (m)" : "Taille de la bosse", &width_, 0.1, 0.6, college ? "%.2f" : "");
    if (interesse) changed |= sliderD(college ? "Vitesse c (m/s)" : "Vitesse de la vague", &speed_, 0.25, 4.0, college ? "%.2f" : "", ImGuiSliderFlags_Logarithmic);
    popSliderWidth();

    if (interesse) {
        ImGui::TextDisabled("Forme");
        int d = dimension_;
        changed |= ImGui::RadioButton("Corde (1D)", &d, 1);
        ImGui::SameLine();
        changed |= ImGui::RadioButton("Surface (2D)", &d, 2);
        dimension_ = d;
    }

    ImGui::SeparatorText("Bords");
    int e = edge_;
    changed |= ImGui::RadioButton(college ? "Absorbants (éponge)" : "Absorbants : la vague s'en va", &e, kAbsorbing);
    changed |= ImGui::RadioButton(college ? "Mur fixe : l'onde revient inversée" : "Mur : la vague rebondit", &e, kWall);
    if (lycee) changed |= ImGui::RadioButton("Bout libre : l'onde revient de même signe", &e, kFree);
    edge_ = e;

    if (etudiant) {
        ImGui::SeparatorText("Calcul");
        pushSliderWidth();
        if (dimension_ == 1) changed |= ImGui::SliderInt("Nombre de cases", &cells1D_, 100, 1600);
        else changed |= ImGui::SliderInt("Cases par côté", &cells2D_, 40, 200);
        double& cfl = dimension_ == 1 ? cfl1D_ : cfl2D_;
        changed |= sliderD("Nombre de Courant C", &cfl, 0.1, 1.2, "%.3f");
        if (chercheur && edge_ == kAbsorbing) {
            ImGui::SeparatorText("Éponge");
            changed |= sliderD("Épaisseur (part du côté)", &spongeFraction_, 0.05, 0.4, "%.2f");
            changed |= sliderD("Force R (ondes courtes)", &spongeStrength_, 1e-4, 0.5, "%.4f", ImGuiSliderFlags_Logarithmic);
            changed |= ImGui::SliderInt("Ordre du profil", &spongeOrder_, 1, 4);
        }
        popSliderWidth();
        if (cfl > cflLimit()) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
            ImGui::TextWrapped("C = %.3f dépasse la limite %.4f de cette grille : le calcul est instable.", cfl, cflLimit());
            ImGui::PopStyleColor();
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

    if (diverged_) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
        ImGui::TextWrapped("Le calcul a explosé : le pas de temps est trop grand pour la taille des cases (C trop grand).");
        ImGui::PopStyleColor();
    }

    if (college) {
        ImGui::SeparatorText("Grandeurs");
        ImGui::Text("t = %.2f s / %.2f s", timeNow(), endTime_);
        ImGui::Text("Traverser le domaine : L / c = %.2f s", length() / speed_);
        if (lycee) ImGui::Text("Largeur de la bosse : %.2f m, durée de passage ≈ %.2f s", width_, width_ / speed_);
    }
}

void WaveModule::drawInvariants(const UiContext& ctx) {
    const Level level = ctx.level;
    if (diverged_) {
        wrapped("Calcul divergent : plus rien à mesurer. Baisser le nombre de Courant C, puis « Recommencer ».", true);
        return;
    }
    const double fraction = energyNow() / initialEnergy_;
    ImGui::TextDisabled("Énergie restante dans le domaine");
    ImGui::ProgressBar(static_cast<float>(std::clamp(fraction, 0.0, 1.0)), ImVec2(-1.0f, 0.0f), strf("%.2f %%", 100.0 * fraction).c_str());

    if (!atLeast(level, Level::College)) {
        wrapped(edge_ == kAbsorbing
                    ? "Sur les bords absorbants la vague s'éteint : l'énergie sort du domaine."
                    : "Contre un mur, la vague rebondit sans rien perdre : l'énergie reste la même.",
                true);
        return;
    }

    ImGui::Text("Hauteur maximale |u| = %.3f", maxAbsNow());
    ImGui::Text("Énergie E / E0 = %.6f", fraction);
    if (atLeast(level, Level::Etudiant)) {
        ImGui::SeparatorText("Conservation");
        if (edge_ == kAbsorbing) {
            wrapped("Avec l'éponge E décroît (c'est voulu).", true);
        } else {
            ImGui::Text("E / E0 - 1 = %+.1e", fraction - 1.0);
            wrapped("Le schéma saute-mouton conserve exactement une énergie discrète : l'écart est de l'arrondi (1e-15).", true);
        }
    }
    if (atLeast(level, Level::Chercheur)) {
        ImGui::SeparatorText("Grille");
        ImGui::Text("dx = %.4f m   dt = %.5f s", dxNow(), dtNow());
        ImGui::Text("C = c dt / dx = %.3f  (limite %.4f)", cflNow(), cflLimit());
        const double k0 = 1.0 / width_, c = speed_;
        const double s = cflNow() * std::sin(0.5 * k0 * dxNow());
        if (std::abs(s) < 1.0) {
            const double omega = 2.0 / dtNow() * std::asin(s);
            ImGui::Text("v_phase / c à k = 1/largeur : %.5f", omega / (c * k0));
        }
        if (dimension_ == 1) {
            double err = 0.0;
            for (int i = 0; i < wave1_->points(); ++i) err = std::max(err, std::abs(wave1_->u()[i] - exact1D(wave1_->x(i), wave1_->time())));
            ImGui::Text("max |u - exacte| = %.2e", err);
        }
    }
}

void WaveModule::drawGraphs(const UiContext& ctx) {
    const Level level = ctx.level;
    const bool showEnergy = atLeast(level, Level::College);
    const bool showExact = dimension_ == 1 && atLeast(level, Level::Lycee);
    const int cols = showEnergy ? 2 : 1;
    if (diverged_) return;  // valeurs hors d'échelle : rien de lisible à tracer

    const std::vector<double> u = cut();
    const double dx = dxNow();
    std::vector<double> xs(u.size()), exact;
    for (std::size_t i = 0; i < u.size(); ++i) xs[i] = static_cast<double>(i) * dx;
    if (showExact) {
        exact.resize(u.size());
        for (std::size_t i = 0; i < u.size(); ++i) exact[i] = exact1D(xs[i], wave1_->time());
    }

    if (ImPlot::BeginSubplots("##graphes", 1, cols, ImVec2(-1.0f, -1.0f))) {
        if (ImPlot::BeginPlot(dimension_ == 1 ? "Profil de la corde u(x)" : "Coupe au milieu de la surface u(x)")) {
            ImPlot::SetupAxes("x (m)", "u", 0, 0);
            ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, length(), ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.6, 1.1, ImPlotCond_Once);
            if (showExact) ImPlot::PlotLine("Exacte (d'Alembert)", xs.data(), exact.data(), static_cast<int>(xs.size()), lineSpec(kBlue));
            ImPlot::PlotLine(atLeast(level, Level::College) ? "Calcul" : "Vague", xs.data(), u.data(), static_cast<int>(xs.size()), lineSpec(kComputed));
            ImPlot::EndPlot();
        }
        if (showEnergy && ImPlot::BeginPlot("Énergie restante E(t) / E0")) {
            ImPlot::SetupAxes("t (s)", "E / E0", ImPlotAxisFlags_AutoFit, 0);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.05, 1.1, ImPlotCond_Once);
            if (energy_.size() > 0)
                ImPlot::PlotLine("E / E0", energy_.x.data(), energy_.y.data(), energy_.size(), lineSpec(kComputed, energy_.offset));
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
}

void WaveModule::drawAnalysis(const UiContext& ctx) {
    const Level level = ctx.level;
    if (diverged_) return;
    const std::vector<double> u = cut();
    const std::vector<double> amp = waves::amplitudeSpectrum(u);
    const std::size_t n = (amp.size() - 1) * 2;  // longueur du signal complété par des zéros
    const double dk = 2.0 * constants::pi / (static_cast<double>(n) * dxNow());
    std::vector<double> ks(amp.size()), theory;
    for (std::size_t m = 0; m < amp.size(); ++m) ks[m] = static_cast<double>(m) * dk;

    // 1D, domaine ouvert, avant l'arrivée aux bords : deux moitiés d'impulsion qui s'éloignent. Spectre = spectre de la bosse
    // × |cos(k c t)|, où la bosse (de hauteur 1) a pour amplitude 2 w sqrt(2 pi) exp(-k² w² / 2) / (n dx) à l'échantillonnage près.
    const bool showTheory = dimension_ == 1 && edge_ == kAbsorbing && atLeast(level, Level::Chercheur);
    if (showTheory) {
        theory.resize(amp.size());
        const double t = wave1_->time();
        for (std::size_t m = 0; m < amp.size(); ++m) {
            const double k = ks[m];
            theory[m] = 2.0 * width_ * std::sqrt(2.0 * constants::pi) / (static_cast<double>(n) * dxNow()) *
                        std::exp(-0.5 * k * k * width_ * width_) * std::abs(std::cos(k * speed_ * t));
        }
    }

    wrapped("Spectre de Fourier de la coupe (FFT) : combien de chaque longueur d'onde la forme contient.", true);
    if (ImPlot::BeginPlot("Spectre d'amplitude", ImVec2(-1.0f, -1.0f))) {
        ImPlot::SetupAxes("k = 2π / λ (rad/m)", "amplitude", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, std::min(ks.back(), 40.0), ImPlotCond_Once);
        if (showTheory) ImPlot::PlotLine("Théorie (bosse × |cos kct|)", ks.data(), theory.data(), static_cast<int>(ks.size()), lineSpec(kBlue));
        ImPlot::PlotLine("FFT du calcul", ks.data(), amp.data(), static_cast<int>(ks.size()), lineSpec(kComputed));
        ImPlot::EndPlot();
    }
}

void WaveModule::drawScene(Renderer& renderer, const UiContext& ctx) {
    const float half = static_cast<float>(0.5 * length());
    const float* edgeCol = edgeColor(edge_);
    std::vector<Vertex> lines;

    if (dimension_ == 1) {
        // poteaux aux extrémités, ligne de repos, début de l'éponge
        addSegment(lines, -half, -1.0f, 0, -half, 1.4f, 0, edgeCol);
        addSegment(lines, half, -1.0f, 0, half, 1.4f, 0, edgeCol);
        const float rest[3] = {0.35f, 0.38f, 0.42f};
        addSegment(lines, -half, 0, 0, half, 0, 0, rest);
        if (edge_ == kAbsorbing) {
            const float layer = static_cast<float>(std::max(4.0, std::round(spongeFraction_ * cells1D_)) * wave1_->dx());
            addSegment(lines, -half + layer, -0.35f, 0, -half + layer, 0.35f, 0, kSpongeColor);
            addSegment(lines, half - layer, -0.35f, 0, half - layer, 0.35f, 0, kSpongeColor);
        }
        renderer.draw(Primitive::Lines, lines);
        if (diverged_) return;

        if (atLeast(ctx.level, Level::Lycee)) {  // solution exacte en bleu, un peu en retrait
            std::vector<Vertex> exact;
            const int n = wave1_->points();
            for (int i = 0; i < n; ++i)
                exact.push_back({static_cast<float>(wave1_->x(i)) - half, kHeight1D * static_cast<float>(exact1D(wave1_->x(i), wave1_->time())), -0.04f,
                                 kBlue[0], kBlue[1], kBlue[2]});
            renderer.draw(Primitive::LineStrip, exact);
        }
        // « rideau » : un trait vertical tous les 4 points, du repos à la courbe, pour lire la hauteur quel que soit l'angle de vue
        std::vector<Vertex> curtain;
        for (int i = 0; i < wave1_->points(); i += 4) {
            float rgb[3];
            divergingColor(wave1_->u()[i] / kColorScale, rgb);
            for (float& comp : rgb) comp *= 0.55f;
            const float x = static_cast<float>(wave1_->x(i)) - half;
            addSegment(curtain, x, 0.0f, 0.0f, x, kHeight1D * static_cast<float>(wave1_->u()[i]), 0.0f, rgb);
        }
        renderer.draw(Primitive::Lines, curtain);
        renderer.draw(Primitive::LineStrip, makeProfile(wave1_->u(), wave1_->dx(), kHeight1D, static_cast<float>(kColorScale)));
        return;
    }

    addFrame(lines, half, half, edgeCol);
    if (edge_ == kAbsorbing) {
        const float layer = static_cast<float>(std::max(4.0, std::round(spongeFraction_ * cells2D_)) * wave2_->dx());
        addFrame(lines, half - layer, half - layer, kSpongeColor);
    }
    renderer.draw(Primitive::Lines, lines);
    if (diverged_) return;
    const int stride = std::max(1, cells2D_ / 60);  // ~60 lignes de maillage dessinées : lisible, léger
    renderer.draw(Primitive::Lines, makeRelief(wave2_->u(), wave2_->pointsX(), wave2_->pointsY(), wave2_->dx(), kHeight2D,
                                               static_cast<float>(kColorScale), stride));
}

void WaveModule::frameCamera(Camera& camera) const {
    camera.target[0] = 0.0f;
    camera.target[1] = 0.0f;
    camera.target[2] = 0.0f;
    camera.distance = 8.5f;
    camera.yaw = 0.55f;
    camera.pitch = 0.62f;
}

// ------------------------------- textes --------------------------------

const char* WaveModule::explanation(Level level) const {
    switch (level) {
        case Level::Vulgarisation:
            return "Jetez un caillou dans l'eau : une bosse se forme, puis un cercle s'étale tout autour. L'eau ne voyage pas "
                   "avec la vague : chaque goutte monte et redescend sur place, c'est la forme qui avance.\n\n"
                   "Rouge : l'eau monte. Bleu : elle descend.\n\n"
                   "Bords absorbants : la vague s'éteint doucement, comme sur une plage en pente douce. Avec un mur : elle rebondit.\n\n"
                   "Un détail : le cercle perd en hauteur en grandissant, parce que la même énergie se répartit sur un cercle de "
                   "plus en plus long.";
        case Level::Interesse:
            return "Une onde est une perturbation qui se propage sans transporter de matière. Sa vitesse dépend du milieu, pas de la "
                   "forme de la bosse : essayez une bosse étroite et une large, elles avancent aussi vite l'une que l'autre.\n\n"
                   "Corde : la bosse se coupe en deux moitiés qui partent chacune de leur côté, deux fois moins hautes. "
                   "Surface : l'onde s'étale en cercle et s'affaiblit.\n\n"
                   "Aux bords : un mur renvoie l'onde retournée (une bosse devient un creux), un bout libre la renvoie dans le "
                   "même sens, un bord absorbant la laisse partir. L'énergie ne se perd pas, elle sort ou elle rebondit.";
        case Level::College:
            return "Vitesse d'une onde : v = distance / temps. Ici la vague traverse le domaine en L / c secondes.\n"
                   "Longueur d'onde λ (distance entre deux bosses), période T (temps entre deux passages) et fréquence f = 1 / T : "
                   "v = λ / T = λ × f.\n\n"
                   "La bosse lâchée au repos se coupe en deux ondes de hauteur ½ qui vont en sens contraires. Contre un mur fixe, "
                   "l'onde revient avec une hauteur opposée ; contre un bout libre, elle revient identique.\n\n"
                   "L'ordinateur découpe la corde en petites cases et le temps en petits pas, puis calcule la hauteur de chaque case "
                   "à partir de ses voisines. La courbe « Énergie restante » montre ce qui reste dans le domaine.";
        case Level::Lycee:
            return "Équation d'onde : ∂²u/∂t² = c² ∂²u/∂x² (en 2D : c² Δu). u est la hauteur, c la vitesse.\n\n"
                   "Solution de d'Alembert pour une bosse f lâchée au repos : u(x, t) = ½ [f(x − ct) + f(x + ct)]. Deux demi-bosses "
                   "qui se séparent à la vitesse c. La courbe bleue de la corde est cette formule (avec la méthode des images pour les murs : "
                   "un mur fixe inverse l'onde, un bout libre non) ; la verte est le calcul.\n\n"
                   "Énergie (par unité de masse) : E = ½ Σ (u_t² + c² u_x²) Δx, moitié cinétique, moitié élastique. "
                   "Elle est constante avec des murs. En 2D, une onde circulaire s'affaiblit comme 1 / √r (r : distance à la source) : "
                   "la même énergie se répartit sur un cercle de longueur 2πr.";
        case Level::Etudiant:
            return "Schéma saute-mouton (leapfrog, ordre 2) : u(i, n+1) = 2 u(i, n) − u(i, n−1) + C² [u(i+1, n) − 2 u(i, n) + u(i−1, n)], "
                   "avec C = c Δt / Δx. C'est le Verlet de la Mécanique appliqué à chaque case. Stable si C ≤ 1 en 1D et C ≤ 1/√2 en 2D "
                   "(essayez C = 1,2 : la grille explose).\n\n"
                   "Il conserve exactement une énergie discrète ½ Σ (u_t² + c² u_x u_x') : l'écart reste à l'arrondi (1e-15). "
                   "À C = 1 en 1D le schéma est exact (l'onde avance d'une case par pas) ; pour mesurer son ordre 2 il faut C < 1.\n\n"
                   "Le panneau Analyse montre la transformée de Fourier rapide (FFT) de la forme. Bords absorbants : "
                   "u_tt + 2σ(x) u_t = c² Δu devant un mur, avec σ qui croît dans la couche (éponge).";
        case Level::Chercheur:
            return "Relation de dispersion de la grille (1D) : sin(ω Δt / 2) = C sin(k Δx / 2) ; en 2D : sin²(ω Δt / 2) = C² [sin²(kx Δx / 2) + "
                   "sin²(ky Δx / 2)]. Un mode propre de la grille est une solution EXACTE du schéma (testé à 1e-14). Au développement : "
                   "ω / (c k) ≈ 1 − (1 − C²)(k Δx)² / 24 : la grille ralentit les courtes longueurs d'onde (dispersion numérique), "
                   "sauf à C = 1 en 1D où elle disparaît.\n\n"
                   "Stabilité de von Neumann : |sin(ω Δt / 2)| ≤ 1 pour tout k donne C ≤ 1 (1D), C ≤ 1/√2 (2D).\n\n"
                   "Éponge : u_tt + 2σ u_t = c² Δu, σ = σmax s^p. Le réglage « force R » fixe σmax pour qu'une onde COURTE soit réduite à R "
                   "après un aller-retour. Mesuré : pour une grande bosse (longueur d'onde du même ordre que la couche) la zone est suramortie "
                   "et réfléchit ; une éponge plus dure ou plus lisse fait pire, seule une couche plus épaisse aide (résidu divisé par 4 quand "
                   "l'épaisseur double). C'est la raison d'être des couches parfaitement adaptées (PML).";
    }
    return "";
}

}  // namespace pl
