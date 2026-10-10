// O1 : corde vibrante, vue comme une chaîne de N masses reliées par des ressorts (le ressort-masse de M2, N fois), fixée aux deux bouts.
//
// Théorie. Corde de longueur L, tension T, masse linéique mu, vitesse c = sqrt(T / mu). On la découpe en N masses m = mu a, a = L / (N + 1),
// aux abscisses x_j = j a (j = 1 .. N), reliées par des ressorts de raideur k = T / a (petits déplacements transversaux u_j, u_0 = u_{N+1} = 0) :
//     m u_j'' = k (u_{j+1} - 2 u_j + u_{j-1}).
// Modes propres : u_j = sin(n pi j / (N + 1)) cos(w_n t) donne  -m w² = k (2 cos(theta) - 2),  theta = n pi / (N + 1),  d'où
//     w_n = 2 sqrt(k / m) sin(n pi / (2 (N + 1)))     (n = 1 .. N),  pulsation de coupure 2 sqrt(k / m).
// Limite continue (N -> infini, a -> 0, sqrt(k / m) = c / a) : w_n -> n pi c / L, soit f_n = n c / (2 L) : la corde de la guitare, des harmoniques
// entiers. À N fini la chaîne est DISPERSIVE : w_n / (n pi c / L) = sin(theta / 2) / (theta / 2) < 1, les harmoniques élevés sont trop graves.
// Énergie : E = sum ½ m v_j² + sum_{j=0..N} ½ k (u_{j+1} - u_j)². Avec q_n = (2 / (N + 1)) sum_j u_j sin(n pi j / (N + 1)) (amplitude du mode n),
// E = sum_n (m (N + 1) / 4) (q_n'² + w_n² q_n²) : l'énergie se range exactement par mode (orthogonalité des sinus).
// Pincement en x0 (triangle de hauteur h), limite continue : u(x, 0) = sum_n b_n sin(n pi x / L),
//     b_n = 2 h L² / (n² pi² x0 (L - x0)) sin(n pi x0 / L)    (les harmoniques multiples de L / x0 sont absents : x0 = L / 2 supprime les pairs).
// Solution exacte de la chaîne pour une position et une vitesse initiales quelconques : transformée en sinus discrète (StringModes).
#pragma once

#include <vector>

#include "physicslab/core/Solver.hpp"

namespace pl::waves {

struct StringProblem {
    int beads = 50;           // N, nombre de masses
    double length = 1.0;      // L [m]
    double tension = 4.0;     // T [N]
    double density = 1.0;     // mu [kg/m]

    // --- grandeurs caractéristiques ---
    double speed() const;               // c = sqrt(T / mu) [m/s]
    double spacing() const;             // a = L / (N + 1) [m]
    double mass() const;                // m = mu a [kg]
    double stiffness() const;           // k = T / a [N/m]
    double fundamental() const;         // f_1 = c / (2 L) de la corde continue [Hz]
    double continuumOmega(int n) const; // n pi c / L [rad/s]
    double chainOmega(int n) const;     // 2 sqrt(k / m) sin(n pi / (2 (N + 1))) [rad/s]
    double maxOmega() const;            // 2 sqrt(k / m) = 2 c / a, pulsation de coupure de la chaîne [rad/s]
    double verletDtLimit() const;       // 2 / maxOmega() = a / c : au-delà Verlet (et le saute-mouton d'O0 avec C = c dt / a) diverge [s]
    double rk4DtLimit() const;          // 2 sqrt(2) / maxOmega() : RK4 est stable sur l'axe imaginaire jusqu'à |w dt| = 2 sqrt(2) [s]

    // --- conditions initiales ---
    std::vector<double> pluck(double x0, double height) const;                 // triangle échantillonné sur les masses
    std::vector<double> modeShape(int n, double amplitude) const;              // amplitude sin(n pi x_j / L)

    // --- intégration numérique : état y = [u_1..u_N | v_1..v_N] (convention des solveurs, n = N) ---
    static State state(const std::vector<double>& u, const std::vector<double>& v);
    OdeFunction rhs() const;
    double kineticEnergy(const State& y) const;
    double potentialEnergy(const State& y) const;
    double energy(const State& y) const { return kineticEnergy(y) + potentialEnergy(y); }

    // --- corde continue (série de Fourier du pincement) ---
    double pluckCoefficient(int n, double x0, double height) const;            // b_n
    double continuumPluck(double x, double t, double x0, double height, int terms) const;
};

// Modes propres de la chaîne : solution exacte (aucune erreur de schéma) et énergie par mode.
class StringModes {
public:
    explicit StringModes(const StringProblem& problem);

    int count() const { return n_; }
    double omega(int n) const { return omega_[n - 1]; }                        // n = 1 .. N [rad/s]

    // Amplitudes q_n = (2 / (N + 1)) sum_j f_j sin(n pi j / (N + 1)), n = 1 .. N (rangées dans un tableau d'indice n - 1).
    std::vector<double> amplitudes(const std::vector<double>& f) const;
    // État exact [u | v] à l'instant t pour la position u0 et la vitesse v0 initiales.
    State exact(const std::vector<double>& u0, const std::vector<double>& v0, double t) const;
    // Énergie de chaque mode pour l'état y ; la somme vaut StringProblem::energy(y) (à l'arrondi).
    std::vector<double> modeEnergies(const State& y) const;

private:
    StringProblem p_;
    int n_;
    std::vector<double> omega_;
    std::vector<double> sine_;  // sine_[(n - 1) * N + (j - 1)] = sin(n pi j / (N + 1))
};

// Erreur finale (norme max sur les positions) d'un solveur sur la chaîne à tEnd, contre la solution exacte, en `steps` pas égaux.
// Mesurer un ordre : tEnd qui n'est pas un multiple de la période du mode principal (voir oscillatorError).
double stringError(const StringProblem& problem, const std::vector<double>& u0, Solver& solver, int steps, double tEnd);

}  // namespace pl::waves
