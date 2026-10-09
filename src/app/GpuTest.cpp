#include "GpuTest.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "physicslab/core/Solver.hpp"
#include "physicslab/core/Vec3.hpp"
#include "physicslab/mechanics/NBody.hpp"
#include "physicslab/render/GpuNBody.hpp"

namespace pl {
namespace {

using Clock = std::chrono::steady_clock;

int g_failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        ++g_failures;
        std::printf("    ÉCHEC : %s\n", what);
    }
}

// Corps de masses égales dans une boule unité (positions seulement : pour les mesures de temps à grand N).
struct Cloud {
    std::vector<double> pos, mass;
    int n = 0;
};

Cloud makeBall(int n, unsigned seed) {
    std::mt19937 rng(seed);
    auto uniform = [&] { return (rng() + 0.5) / 4294967296.0; };
    Cloud c;
    c.n = n;
    c.mass.assign(n, 1.0 / n);
    c.pos.reserve(3 * static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        double x, y, z;
        do {
            x = 2.0 * uniform() - 1.0;
            y = 2.0 * uniform() - 1.0;
            z = 2.0 * uniform() - 1.0;
        } while (x * x + y * y + z * z > 1.0);
        c.pos.push_back(x);
        c.pos.push_back(y);
        c.pos.push_back(z);
    }
    return c;
}

Cloud fromProblem(const NBodyProblem& p) {
    Cloud c;
    c.n = p.count();
    c.mass = p.mass;
    for (const Vec3& r : p.position) {
        c.pos.push_back(r.x);
        c.pos.push_back(r.y);
        c.pos.push_back(r.z);
    }
    return c;
}

double norm3(const double* v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

// Écart entre deux tableaux d'accélérations.
struct ErrorStats {
    double maxRel = 0.0, medianRel = 0.0, p99Rel = 0.0;  // |Δa_i| / |a_i|, corps par corps
    double globalRel = 0.0;                              // ||Δa|| / ||a|| sur tout le système
    double maxVsScale = 0.0;                             // max_i |Δa_i| / S_i, S_i = somme des modules des termes de la somme
};

// S_i = G sum_j m_j |r_j - r_i| / (|r|^2 + eps^2)^(3/2) : l'échelle naturelle de l'erreur d'une somme de N termes
// (un arrondi relatif u sur chaque terme donne une erreur de l'ordre de u S_i, même si les termes se compensent dans a_i).
void termScales(const double* pos, const double* mass, int n, double G, double eps, std::vector<double>& s) {
    s.assign(n, 0.0);
    const double eps2 = eps * eps;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double d[3] = {pos[3 * j] - pos[3 * i], pos[3 * j + 1] - pos[3 * i + 1], pos[3 * j + 2] - pos[3 * i + 2]};
            const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + eps2;
            if (r2 > 0.0) s[i] += G * mass[j] * norm3(d) / (r2 * std::sqrt(r2));
        }
}

ErrorStats compare(const std::vector<double>& ref, const std::vector<double>& test, const std::vector<double>& scale) {
    const int n = static_cast<int>(ref.size() / 3);
    ErrorStats e;
    std::vector<double> rel;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        const double d[3] = {test[3 * i] - ref[3 * i], test[3 * i + 1] - ref[3 * i + 1], test[3 * i + 2] - ref[3 * i + 2]};
        const double dn = norm3(d), rn = norm3(&ref[3 * i]);
        num += dn * dn;
        den += rn * rn;
        if (rn > 1e-300) rel.push_back(dn / rn);
        if (!scale.empty() && scale[i] > 0.0) e.maxVsScale = std::max(e.maxVsScale, dn / scale[i]);
    }
    e.globalRel = den > 0.0 ? std::sqrt(num / den) : 0.0;
    if (!rel.empty()) {
        std::sort(rel.begin(), rel.end());
        e.maxRel = rel.back();
        e.medianRel = rel[rel.size() / 2];
        e.p99Rel = rel[std::min(rel.size() - 1, static_cast<std::size_t>(0.99 * static_cast<double>(rel.size())))];
    }
    return e;
}

// Force totale relative |sum m_i a_i| / sum m_i |a_i| : nulle si les forces par paire sont exactement opposées.
double netForce(const std::vector<double>& acc, const double* mass) {
    const int n = static_cast<int>(acc.size() / 3);
    double fx = 0.0, fy = 0.0, fz = 0.0, scale = 0.0;
    for (int i = 0; i < n; ++i) {
        fx += mass[i] * acc[3 * i];
        fy += mass[i] * acc[3 * i + 1];
        fz += mass[i] * acc[3 * i + 2];
        scale += mass[i] * norm3(&acc[3 * i]);
    }
    return scale > 0.0 ? std::sqrt(fx * fx + fy * fy + fz * fz) / scale : 0.0;
}

// L'entrée telle que le GPU float la voit : positions centrées, positions et masses arrondies au float, eps^2 arrondi.
struct RoundedInput {
    std::vector<double> pos, mass;
    double eps = 0.0;
};

RoundedInput roundToFloat(const Cloud& c, double eps, bool centered) {
    RoundedInput r;
    r.pos.resize(3 * static_cast<std::size_t>(c.n));
    if (centered) {
        centerOnBarycenter(c.pos.data(), c.mass.data(), c.n, r.pos.data());
    } else {
        r.pos = c.pos;
    }
    for (double& x : r.pos) x = static_cast<double>(static_cast<float>(x));
    for (double m : c.mass) r.mass.push_back(static_cast<double>(static_cast<float>(m)));
    r.eps = std::sqrt(static_cast<double>(static_cast<float>(eps * eps)));
    return r;
}

void printStats(const char* label, const ErrorStats& e) {
    std::printf("    %-34s max %.2e  médiane %.2e  p99 %.2e  globale %.2e  max/S_i %.2e\n", label, e.maxRel, e.medianRel,
                e.p99Rel, e.globalRel, e.maxVsScale);
}

// Seuil de l'écart global en float : un arrondi relatif u = 6e-8 par terme, N termes sommés au hasard -> u sqrt(N) (mesuré : de
// 1,6e-8 pour N = 3 à 1,05e-6 pour N = 5000, soit environ le tiers du seuil ; l'arrondi des entrées ajoute la constante).
double toleranceFloat(int n) { return 5e-8 * (3.0 + std::sqrt(static_cast<double>(n))); }
constexpr double kToleranceNetForce = 1e-6;  // |sum m a| / sum m |a| en float (mesuré : 1e-8 à 2e-8)
constexpr double kToleranceDouble = 1e-12;   // GPU double contre CPU double (mesuré : 3e-16)
// Le double du shader est émulé par le pilote de cette carte (Intel UHD) : environ 3e6 interactions/s mesurées, soit des milliers de
// fois moins que le float. Au-delà de ce nombre de corps on ne le lance plus (5000 corps : 9 s d'un seul envoi, risque de redémarrage du pilote).
constexpr int kMaxBodiesDouble = 1000;
// Intégration (mesuré en float sur Intel UHD, voir les tables du test) :
//  * huit, 5 périodes : écart au CPU 2,0e-4 (croissance linéaire : le bruit d'arrondi des mises à jour de v et de x est amplifié par
//    le cisaillement de l'orbite) ;
//  * amas : l'écart croît comme exp(lambda t), lambda ~ 1 : 3,9e-5 à t = 2 (N = 100), 2,0e-5 à t = 1 (N = 500) ;
//  * dérives GPU float au plus 2,2e-6 (E, P, L) ; le CPU double ne dérive que de l'erreur du schéma (E 1e-10 à 8e-7 selon dt).
constexpr double kMaxTrajectoryGapEight = 2e-3;
constexpr double kMaxEarlyGapCluster = 1e-3;   // amas, avant que le chaos n'amplifie l'arrondi
constexpr double kMaxDriftFloat = 3e-5;        // E (au-delà du CPU), P et L en float, sur 2e4 pas

// Un état fixé : compare GPU float (et GPU double si disponible) au CPU double. Renvoie l'écart GPU float / CPU double.
// `strict` faux : on ne vérifie que l'écart à l'entrée arrondie (cas sans centrage, où l'erreur vient de l'arrondi des positions).
ErrorStats accuracyCase(const char* name, const Cloud& c, double G, double eps, GpuNBody& gpuF, GpuNBody* gpuD, bool centering,
                        bool strict = true) {
    std::printf("  %s : N = %d, ε = %g%s\n", name, c.n, eps, centering ? "" : " (SANS centrage)");
    std::vector<double> ref(3 * static_cast<std::size_t>(c.n)), refRounded(ref.size()), gpu(ref.size());
    nbody::accelerations(c.pos.data(), c.mass.data(), c.n, G, eps, ref.data());

    std::vector<double> scale;
    if (c.n <= 6000) termScales(c.pos.data(), c.mass.data(), c.n, G, eps, scale);

    gpuF.setCentering(centering);
    gpuF.accelerations(c.pos.data(), c.mass.data(), c.n, G, eps, gpu.data());

    bool finite = true;
    for (double v : gpu) finite = finite && std::isfinite(v);
    check(finite, "valeurs non finies dans les accélérations GPU");

    const RoundedInput in = roundToFloat(c, eps, centering);
    nbody::accelerations(in.pos.data(), in.mass.data(), c.n, G, in.eps, refRounded.data());

    const ErrorStats vsRaw = compare(ref, gpu, scale);
    const ErrorStats vsRounded = compare(refRounded, gpu, scale);
    const ErrorStats roundingOnly = compare(ref, refRounded, scale);
    std::printf("    temps du noyau GPU float : %.3f ms\n", 1e3 * gpuF.lastKernelSeconds());
    printStats("GPU float / CPU double", vsRaw);
    printStats("GPU float / CPU (entrée arrondie)", vsRounded);
    printStats("arrondi de l'entrée seul (CPU/CPU)", roundingOnly);
    std::printf("    force totale relative : CPU %.2e  GPU float %.2e\n", netForce(ref, c.mass.data()), netForce(gpu, c.mass.data()));
    if (strict) check(vsRaw.globalRel <= toleranceFloat(c.n), "écart global GPU float / CPU double trop grand");
    check(vsRounded.globalRel <= toleranceFloat(c.n), "écart global GPU float / CPU (entrée arrondie) trop grand");
    check(netForce(gpu, c.mass.data()) <= kToleranceNetForce, "force totale GPU float non nulle (3e loi de Newton)");

    if (gpuD && c.n <= kMaxBodiesDouble) {
        gpuD->setCentering(centering);
        std::vector<double> gd(ref.size());
        gpuD->accelerations(c.pos.data(), c.mass.data(), c.n, G, eps, gd.data());
        const ErrorStats e = compare(ref, gd, scale);
        std::printf("    temps du noyau GPU double : %.3f ms\n", 1e3 * gpuD->lastKernelSeconds());
        printStats("GPU double / CPU double", e);
        check(e.globalRel <= kToleranceDouble, "écart global GPU double / CPU double trop grand");
    }
    return vsRaw;
}

double seconds(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

// ------------------------- intégration : GPU kick-drift-kick contre CPU Verlet -------------------------

// Dérive maximale des grandeurs conservées sur les échantillons d'une trajectoire, en relatif :
//   énergie |E - E0| / |E0| ; impulsion |P - P0| / sum m|v| ; moment cinétique |L - L0| / sum m|r x v| (états initiaux).
struct Drift {
    double energy = 0.0, momentum = 0.0, angular = 0.0;
};

struct Trace {
    std::vector<State> states;   // état (positions | vitesses) après chaque lot de pas
    std::vector<Drift> drifts;   // dérive cumulée maximale jusqu'à ce lot
    double wallSeconds = 0.0;    // temps passé dans les pas (hors relectures et énergie)
    double kernelSeconds = 0.0;  // temps GPU des accélérations
};

class Conservation {
public:
    explicit Conservation(const NBodyProblem& p) : p_(p) {
        const State y0 = p.initialState();
        e0_ = p.energy(y0);
        p0_ = p.momentum(y0);
        l0_ = p.angularMomentum(y0);
        const int n = p.count();
        for (int i = 0; i < n; ++i) {
            const Vec3 r{y0[3 * i], y0[3 * i + 1], y0[3 * i + 2]}, v{y0[3 * n + 3 * i], y0[3 * n + 3 * i + 1], y0[3 * n + 3 * i + 2]};
            pScale_ += p.mass[i] * v.norm();
            lScale_ += p.mass[i] * cross(r, v).norm();
        }
    }
    void update(const State& y) {
        drift_.energy = std::max(drift_.energy, std::abs(p_.energy(y) - e0_) / std::abs(e0_));
        if (pScale_ > 0.0) drift_.momentum = std::max(drift_.momentum, (p_.momentum(y) - p0_).norm() / pScale_);
        if (lScale_ > 0.0) drift_.angular = std::max(drift_.angular, (p_.angularMomentum(y) - l0_).norm() / lScale_);
    }
    const Drift& drift() const { return drift_; }

private:
    const NBodyProblem& p_;
    double e0_ = 0.0, pScale_ = 0.0, lScale_ = 0.0;
    Vec3 p0_, l0_;
    Drift drift_;
};

// `samples` lots de `stepsPerSample` pas de `dt`. CPU : Verlet des vitesses de Solver.hpp en double.
Trace traceCpu(const NBodyProblem& p, double dt, int stepsPerSample, int samples) {
    Trace tr;
    VelocityVerlet solver;
    const OdeFunction f = p.rhs();
    State y = p.initialState();
    Conservation cons(p);
    double t = 0.0;
    for (int s = 0; s < samples; ++s) {
        const auto t0 = Clock::now();
        for (int k = 0; k < stepsPerSample; ++k) t += solver.step(f, t, y, dt);
        tr.wallSeconds += seconds(t0, Clock::now());
        cons.update(y);
        tr.states.push_back(y);
        tr.drifts.push_back(cons.drift());
    }
    return tr;
}

// GPU : l'état reste dans les SSBO pendant un lot ; il est relu entre deux lots seulement.
Trace traceGpu(GpuNBody& gpu, const NBodyProblem& p, double dt, int stepsPerSample, int samples) {
    Trace tr;
    const int n = p.count();
    const State y0 = p.initialState();
    gpu.setState(y0.data(), y0.data() + 3 * n, p.mass.data(), n);
    Conservation cons(p);
    State y(6 * static_cast<std::size_t>(n));
    for (int s = 0; s < samples; ++s) {
        const auto t0 = Clock::now();
        gpu.step(p.G, p.softening, dt, stepsPerSample);
        tr.kernelSeconds += gpu.lastKernelSeconds();
        gpu.readState(y.data(), y.data() + 3 * n);  // la relecture attend la fin des envois : le temps est celui du lot
        tr.wallSeconds += seconds(t0, Clock::now());
        cons.update(y);
        tr.states.push_back(y);
        tr.drifts.push_back(cons.drift());
    }
    return tr;
}

// Un cas : table GPU float contre CPU double (même schéma, même pas) à chaque échantillon. `withCpu` faux : GPU seul.
struct IntegrationResult {
    std::vector<double> distances;  // GPU float / CPU double à chaque échantillon (vide sans CPU)
    double finalDistance = 0.0;     // au dernier échantillon
    Drift gpu, cpu;
    double microsecondsPerStep = 0.0;
};

// Les dérives du GPU float doivent rester petites : E au-delà de l'erreur propre du schéma (celle du CPU), P et L en valeur absolue.
void checkDrifts(const IntegrationResult& r) {
    check(r.gpu.energy <= r.cpu.energy + kMaxDriftFloat, "dérive d'énergie du GPU float trop grande");
    check(r.gpu.momentum <= kMaxDriftFloat, "dérive d'impulsion du GPU float trop grande");
    check(r.gpu.angular <= kMaxDriftFloat, "dérive de moment cinétique du GPU float trop grande");
}

IntegrationResult integrationCase(const char* name, const NBodyProblem& p, GpuNBody& gpu, double dt, int stepsPerSample, int samples,
                                  bool withCpu) {
    std::printf("  %s : N = %d, ε = %g, dt = %g, %d pas (t = %g)\n", name, p.count(), p.softening, dt, stepsPerSample * samples,
                dt * stepsPerSample * samples);
    const Trace g = traceGpu(gpu, p, dt, stepsPerSample, samples);
    Trace c;
    if (withCpu) c = traceCpu(p, dt, stepsPerSample, samples);

    IntegrationResult r;
    r.gpu = g.drifts.back();
    r.microsecondsPerStep = 1e6 * g.wallSeconds / (static_cast<double>(stepsPerSample) * samples);
    std::printf("    %8s %14s %12s %12s %12s %12s\n", "t", withCpu ? "|GPU - CPU|" : "", "dE/E GPU", "dE/E CPU", "dP GPU", "dL GPU");
    // On affiche une dizaine de lignes au plus.
    const int stride = std::max(1, samples / 10);
    for (int s = 0; s < samples; ++s) {
        if ((s + 1) % stride != 0 && s != samples - 1) continue;
        const double t = dt * stepsPerSample * (s + 1);
        char dist[32] = "", ec[32] = "-";
        if (withCpu) {
            std::snprintf(dist, sizeof(dist), "%.3e", p.distance(g.states[s], c.states[s]));
            std::snprintf(ec, sizeof(ec), "%.3e", c.drifts[s].energy);
        }
        std::printf("    %8.3f %14s %12.3e %12s %12.3e %12.3e\n", t, dist, g.drifts[s].energy, ec, g.drifts[s].momentum,
                    g.drifts[s].angular);
    }
    if (withCpu) {
        for (int s = 0; s < samples; ++s) r.distances.push_back(p.distance(g.states[s], c.states[s]));
        r.finalDistance = r.distances.back();
        r.cpu = c.drifts.back();
        std::printf("    dérives maximales (GPU float | CPU double) : E %.2e | %.2e   P %.2e | %.2e   L %.2e | %.2e\n", r.gpu.energy,
                    r.cpu.energy, r.gpu.momentum, r.cpu.momentum, r.gpu.angular, r.cpu.angular);
        std::printf("    temps par pas : GPU %.0f µs (dont noyau %.0f µs), CPU %.0f µs\n", r.microsecondsPerStep,
                    1e6 * g.kernelSeconds / (static_cast<double>(stepsPerSample) * samples),
                    1e6 * c.wallSeconds / (static_cast<double>(stepsPerSample) * samples));
    } else {
        std::printf("    dérives maximales GPU float : E %.2e   P %.2e   L %.2e\n", r.gpu.energy, r.gpu.momentum, r.gpu.angular);
        std::printf("    temps par pas : GPU %.0f µs (dont noyau %.0f µs)\n", r.microsecondsPerStep,
                    1e6 * g.kernelSeconds / (static_cast<double>(stepsPerSample) * samples));
    }
    return r;
}

}  // namespace

int runGpuTest(const std::string& shaderDir, int maxN) {
    std::setvbuf(stdout, nullptr, _IONBF, 0);  // sortie immédiate, même redirigée (utile si le pilote bloque)
    const GpuLimits lim = queryGpuLimits();
    std::printf("=== Test GPU (M7) ===\n");
    std::printf("Carte : %s | %s\n", lim.renderer.c_str(), lim.version.c_str());
    std::printf("Groupe de travail : %d threads max (%d x %d x %d), mémoire partagée %d octets, SSBO %lld octets, fp64 annoncé : %s\n",
                lim.maxGroupInvocations, lim.maxGroupSize[0], lim.maxGroupSize[1], lim.maxGroupSize[2], lim.maxSharedBytes,
                lim.maxStorageBlockBytes, lim.doublePrecision ? "oui" : "non");

    GpuNBody gpuF;
    if (!gpuF.init(shaderDir, GpuPrecision::Float)) {
        std::printf("ÉCHEC : le compute shader float ne se compile pas.\n");
        return 1;
    }
    GpuNBody gpuD;
    const bool haveDouble = gpuD.init(shaderDir, GpuPrecision::Double);
    std::printf("Calcul en double dans le shader : %s\n", haveDouble ? "disponible" : "refusé par le pilote");

    // ---- 1. Précision à état fixé ----
    std::printf("\n[1] Accélérations à état fixé, GPU contre CPU double (écarts relatifs)\n");
    const double G = 1.0;
    const NBodyProblem eight = NBodyProblem::figureEight();
    GpuNBody* dbl = haveDouble ? &gpuD : nullptr;
    accuracyCase("huit de Chenciner-Montgomery", fromProblem(eight), G, 0.0, gpuF, dbl, true);
    const NBodyProblem lagrange = NBodyProblem::lagrangeTriangle(1.0, 1.0, G);
    accuracyCase("triangle de Lagrange", fromProblem(lagrange), G, 0.0, gpuF, dbl, true);
    for (int n : {100, 1000, 5000}) {
        const NBodyProblem cluster = NBodyProblem::randomCluster(n, 42, 1.0, 0.05);
        accuracyCase("amas aléatoire (graine 42)", fromProblem(cluster), G, 0.05, gpuF, dbl, true);
    }

    // Tailles qui ne sont pas un multiple du groupe de travail (256) : dernière tuile incomplète ; un seul corps : accélération nulle.
    for (int n : {1, 2, 127, 129, 255, 257}) accuracyCase("cas limite", makeBall(n, 7), G, 0.05, gpuF, dbl, true);

    // Amas éloigné de l'origine : le centrage en double avant la conversion en float doit absorber le décalage.
    {
        Cloud c = fromProblem(NBodyProblem::randomCluster(1000, 42, 1.0, 0.05));
        for (int i = 0; i < c.n; ++i) {
            c.pos[3 * i] += 1000.0;
            c.pos[3 * i + 1] -= 2000.0;
            c.pos[3 * i + 2] += 500.0;
        }
        const ErrorStats centered = accuracyCase("amas décalé de (1000, -2000, 500)", c, G, 0.05, gpuF, dbl, true);
        const ErrorStats raw = accuracyCase("amas décalé de (1000, -2000, 500)", c, G, 0.05, gpuF, nullptr, false, false);
        std::printf("    gain du centrage : écart global divisé par %.0f\n", raw.globalRel / centered.globalRel);
        check(raw.globalRel > 20.0 * centered.globalRel, "le centrage devrait réduire l'écart d'un facteur d'au moins 20");
    }

    // ---- 2. Temps ----
    std::printf("\n[2] Tailles du groupe de travail (N = 16384, GPU float, meilleur de 8 envois)\n");
    {
        const Cloud c = makeBall(16384, 3);
        for (int wg : {32, 64, 128, 256, 512, 1024}) {
            if (wg > lim.maxGroupInvocations || wg * 16 > lim.maxSharedBytes) continue;
            GpuNBody g;
            if (!g.init(shaderDir, GpuPrecision::Float, wg)) continue;
            g.setBodies(c.pos.data(), c.mass.data(), c.n);
            g.compute(G, 0.05);  // échauffement
            double best = 1e30;
            for (int rep = 0; rep < 8; ++rep) {
                g.compute(G, 0.05);
                best = std::min(best, g.lastKernelSeconds());
            }
            std::printf("    groupe %4d : %.3f ms, %.2f Ginteractions/s\n", wg, 1e3 * best,
                        static_cast<double>(c.n) * c.n / best * 1e-9);
            g.shutdown();
        }
    }

    std::printf("\n[3] Temps et débit (interactions ordonnées : N^2 ; une paire = 2 interactions)\n");
    std::printf("    %8s %12s %12s %14s %12s %12s %10s\n", "N", "GPU noyau", "GPU double", "GPU total", "CPU double", "GPU Gint/s",
                "CPU Gint/s");
    const auto cell = [](double s) {  // durée en ms, ou « - » si non mesurée
        char buf[32];
        if (s < 0.0) std::snprintf(buf, sizeof(buf), "-");
        else std::snprintf(buf, sizeof(buf), "%.2f ms", 1e3 * s);
        return std::string(buf);
    };
    for (int n : {1000, 2000, 4000, 8000, 16000, 32000, 64000, 100000, 200000}) {
        if (n > maxN) break;
        const Cloud c = makeBall(n, 5);
        std::vector<double> acc(3 * static_cast<std::size_t>(n)), accCpu(acc.size());
        gpuF.setBodies(c.pos.data(), c.mass.data(), n);
        gpuF.compute(G, 0.05);  // échauffement (la fréquence de la carte monte et le débit est mesuré)
        double total = 1e30, kernel = 1e30;  // meilleur de 3 : la fréquence d'une carte intégrée fluctue
        for (int rep = 0; rep < 3; ++rep) {
            const auto t0 = Clock::now();
            gpuF.setBodies(c.pos.data(), c.mass.data(), n);
            gpuF.compute(G, 0.05);
            gpuF.readAccelerations(acc.data());
            total = std::min(total, seconds(t0, Clock::now()));
            kernel = std::min(kernel, gpuF.lastKernelSeconds());
        }
        const double interactions = static_cast<double>(n) * static_cast<double>(n);

        double kernelDouble = -1.0;
        if (haveDouble && n <= 2 * kMaxBodiesDouble) {
            gpuD.setBodies(c.pos.data(), c.mass.data(), n);
            gpuD.compute(G, 0.05);
            kernelDouble = gpuD.lastKernelSeconds();
        }
        double cpu = -1.0;
        if (n <= 16000) {
            const auto c0 = Clock::now();
            nbody::accelerations(c.pos.data(), c.mass.data(), n, G, 0.05, accCpu.data());
            cpu = seconds(c0, Clock::now());
        }
        char cpuRate[32] = "-";
        if (cpu >= 0.0) std::snprintf(cpuRate, sizeof(cpuRate), "%.3f", interactions / cpu * 1e-9);
        std::printf("    %8d %12s %12s %14s %12s %12.2f %10s\n", n, cell(kernel).c_str(), cell(kernelDouble).c_str(),
                    cell(total).c_str(), cell(cpu).c_str(), interactions / kernel * 1e-9, cpuRate);
        std::fflush(stdout);
    }

    // ---- 4. Intégration ----
    std::printf("\n[4] Intégration kick-drift-kick sur le GPU contre Verlet des vitesses du CPU (même schéma, même pas)\n");
    {
        const NBodyProblem huit = NBodyProblem::figureEight();
        const double dtEight = nbody::kFigureEightPeriod / 2000.0;

        // 4a. Le même schéma en double doit coïncider avec le CPU à l'arrondi près : valide la logique du shader sans l'arrondi float.
        if (haveDouble) {
            const double dt = nbody::kFigureEightPeriod / 1000.0;
            const Trace d = traceGpu(gpuD, huit, dt, 10, 2);
            const Trace f = traceGpu(gpuF, huit, dt, 10, 2);
            const Trace c = traceCpu(huit, dt, 10, 2);
            const double dd = huit.distance(d.states.back(), c.states.back());
            const double df = huit.distance(f.states.back(), c.states.back());
            std::printf("  huit, 20 pas : écart au CPU double : GPU double %.2e, GPU float %.2e\n", dd, df);
            check(dd <= 1e-12, "le schéma GPU en double devrait coïncider avec le CPU à l'arrondi près");
            check(df <= 1e-5, "le schéma GPU en float s'écarte trop du CPU après 20 pas");

            // 4b. Le même nombre de pas en lots de 10 ou de 5 : l'accélération mémorisée et la coupure des demi-coups ne changent rien.
            const Trace split = traceGpu(gpuD, huit, dt, 5, 4);
            const double ds = huit.distance(d.states.back(), split.states.back());
            std::printf("  huit, 20 pas en 2 lots de 10 ou 4 lots de 5 (GPU double) : écart %.2e\n", ds);
            check(ds <= 1e-13, "découper le calcul en lots ne devrait rien changer (GPU double)");
        }

        const IntegrationResult eight = integrationCase("huit de Chenciner-Montgomery (5 périodes)", huit, gpuF, dtEight, 2000, 5, true);
        check(eight.finalDistance <= kMaxTrajectoryGapEight, "huit : écart GPU float / CPU double trop grand après 5 périodes");
        checkDrifts(eight);

        // Amas chaotique (lambda ~ 1) : l'arrondi float de 1e-7 est amplifié ; seul le début de la trajectoire est comparable.
        const NBodyProblem c100 = NBodyProblem::randomCluster(100, 42, 1.0, 0.05);
        const IntegrationResult r100 = integrationCase("amas aléatoire (graine 42)", c100, gpuF, 1e-3, 1000, 20, true);
        check(r100.distances[1] <= kMaxEarlyGapCluster, "amas N = 100 : écart GPU float / CPU double trop grand à t = 2");
        checkDrifts(r100);

        const NBodyProblem c500 = NBodyProblem::randomCluster(500, 42, 1.0, 0.05);
        const IntegrationResult r500 = integrationCase("amas aléatoire (graine 42)", c500, gpuF, 2e-3, 100, 5, true);
        check(r500.finalDistance <= kMaxEarlyGapCluster, "amas N = 500 : écart GPU float / CPU double trop grand à t = 1");
        checkDrifts(r500);

        // 5000 corps : le CPU (2 évaluations de force par pas) serait trop long ; on ne regarde que les grandeurs conservées.
        const NBodyProblem c5000 = NBodyProblem::randomCluster(5000, 42, 1.0, 0.05);
        const IntegrationResult r5000 = integrationCase("amas aléatoire (graine 42)", c5000, gpuF, 2e-3, 100, 5, false);
        check(r5000.gpu.energy <= kMaxDriftFloat, "amas N = 5000 : dérive d'énergie du GPU float trop grande");
        check(r5000.gpu.momentum <= kMaxDriftFloat, "amas N = 5000 : dérive d'impulsion du GPU float trop grande");
        check(r5000.gpu.angular <= kMaxDriftFloat, "amas N = 5000 : dérive de moment cinétique du GPU float trop grande");
    }

    gpuF.shutdown();
    gpuD.shutdown();
    std::printf("\n%s\n", g_failures == 0 ? "gpu-test : OK" : "gpu-test : ÉCHEC");
    return g_failures == 0 ? 0 : 1;
}

}  // namespace pl
