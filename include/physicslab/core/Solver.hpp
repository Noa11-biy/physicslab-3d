// Interface commune des intégrateurs d'EDO du premier ordre : dy/dt = f(t, y).
//
// Convention pour les systèmes mécaniques : y = [positions | vitesses].
// Verlet/symplectique (M1) s'appuieront sur cette convention : la seconde moitié
// de f(t, y) donne l'accélération, qui ne doit pas dépendre des vitesses.
#pragma once

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

    // Avance `y` de t à t + dt et renvoie la durée réellement avancée
    // (égale à dt pour un pas fixe ; peut être plus courte pour RK45 adaptatif).
    virtual double step(const OdeFunction& f, double t, State& y, double dt) = 0;
};

// Euler explicite : y_{n+1} = y_n + dt f(t_n, y_n). Ordre 1, non symplectique.
// Sert de référence "naïve" pour montrer la dérive d'énergie.
class ExplicitEuler final : public Solver {
public:
    const char* name() const override { return "Euler explicite"; }
    double step(const OdeFunction& f, double t, State& y, double dt) override;

private:
    State k_;  // tampon réutilisé pour éviter une allocation à chaque pas
};

}  // namespace pl
