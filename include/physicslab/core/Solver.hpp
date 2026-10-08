// Interface commune des intégrateurs d'EDO du premier ordre : dy/dt = f(t, y).
//
// Convention pour les systèmes mécaniques : y = [positions | vitesses], chaque moitié de
// taille n = y.size() / 2. La seconde moitié de f(t, y) est l'accélération.
// Euler symplectique et Verlet s'appuient sur cette convention ; Euler, RK4 et RK45 sont
// des intégrateurs génériques qui n'en dépendent pas.
#pragma once

#include <array>
#include <functional>
#include <vector>

namespace pl {

using State = std::vector<double>;

// Écrit dy/dt dans `dydt` (déjà dimensionné comme `y`).
using OdeFunction = std::function<void(double t, const State& y, State& dydt)>;

class Solver {
public:
    virtual ~Solver() = default;

    virtual const char* name() const = 0;

    // Avance `y` de t à t + h et renvoie h, avec 0 < h <= dt.
    // Un pas fixe renvoie toujours dt ; RK45 adaptatif peut renvoyer moins : l'appelant
    // (World::step) recommence jusqu'à avoir couvert tout l'intervalle voulu.
    virtual double step(const OdeFunction& f, double t, State& y, double dt) = 0;
};

// Euler explicite : y_{n+1} = y_n + dt f(t_n, y_n). Ordre 1, non symplectique.
// Référence "naïve" pour montrer l'accumulation d'erreur et la dérive d'énergie.
class ExplicitEuler final : public Solver {
public:
    const char* name() const override { return "Euler explicite"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

private:
    State k_;  // tampon réutilisé pour éviter une allocation à chaque pas
};

// Euler symplectique (semi-implicite) : v_{n+1} = v_n + dt a(x_n, v_n) puis x_{n+1} = x_n + dt v_{n+1}.
// Ordre 1, mais conserve la structure hamiltonienne (énergie bornée sur les systèmes conservatifs).
class SymplecticEuler final : public Solver {
public:
    const char* name() const override { return "Euler symplectique"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

private:
    State k_;
};

// Verlet des vitesses : x_{n+1} = x_n + dt v_n + dt^2/2 a_n ; v_{n+1} = v_n + dt/2 (a_n + a_{n+1}).
// Ordre 2, symplectique et réversible quand l'accélération ne dépend que de la position.
// Si elle dépend de la vitesse (frottement), v* = v_n + dt a_n sert de prédicteur pour a_{n+1}.
class VelocityVerlet final : public Solver {
public:
    const char* name() const override { return "Verlet (vitesses)"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

private:
    State k1_, k2_, tmp_;
};

// Runge-Kutta classique d'ordre 4 (4 évaluations de f par pas).
class RK4 final : public Solver {
public:
    const char* name() const override { return "RK4"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

private:
    State k1_, k2_, k3_, k4_, tmp_;
};

// Dormand-Prince 5(4) à pas adaptatif. Le pas est choisi pour que l'erreur locale estimée
// reste sous  absTol + relTol * |y|  (norme RMS) ; `dt` n'est qu'un plafond.
class RK45 final : public Solver {
public:
    double relTol = 1e-8;
    double absTol = 1e-10;

    const char* name() const override { return "RK45 adaptatif"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

    // Statistiques cumulées depuis le dernier resetStats().
    int acceptedSteps() const { return accepted_; }
    int rejectedSteps() const { return rejected_; }
    int evaluations() const { return evaluations_; }
    void resetStats() { accepted_ = rejected_ = evaluations_ = 0; }

private:
    std::array<State, 7> k_;
    State tmp_;
    double hNext_ = 0.0;  // pas proposé pour la prochaine tentative (0 = pas encore connu)
    int accepted_ = 0, rejected_ = 0, evaluations_ = 0;
};

}  // namespace pl
