# M2 : ressort-masse (libre, amorti, forcé) et application modulaire

Commit : `4cc1c7b`. Date : 2026-10-08. `--sim 2`.

## Objectif

Le système le plus étudié de la physique, et le premier où l'on voit **ce qui distingue les schémas sur la durée** : l'énergie. Plus un refactoring : avec deux simulations à héberger, l'application est devenue modulaire.

## Application modulaire

Interface `SimulationModule` : `title`, `explanation(niveau)`, `reset`, `update`, `drawControls`, `drawInvariants`, `drawGraphs`, `hasAnalysis` / `drawAnalysis`, `drawScene`, `frameCamera`. L'application ouvre les fenêtres ImGui (Simulation, Explication, Invariants, Graphes, Analyse) et le module les remplit. Chaque nouveau module n'a plus qu'à réaliser cette interface et à s'inscrire dans `Application.cpp`. C'est le commit qui a le plus de suppressions (609) : le code du projectile a été déplacé dans son module.

## Physique

`m x'' = -k x - c x' + F cos(Ωt)`. Pulsation propre `ω0 = √(k/m)`, taux d'amortissement `ζ = c / (2√(km))` : trois régimes (sous-critique `ζ < 1`, critique `ζ = 1`, sur-critique `ζ > 1`) avec leurs solutions exactes ; cas forcé : amplitude et déphasage en régime établi, résonance près de `Ω = ω0 √(1 - 2ζ²)`.

## Résultats mesurés

- **Euler explicite multiplie l'énergie** par `(1 + ω0² dt²)` à chaque pas : mesuré 7,5543, prédit 7,5543.
- **Verlet et Euler symplectique** : énergie bornée, sans dérive, sur 100 périodes (l'énergie oscille autour de la valeur exacte).
- **RK4** : dissipation lente, d'ordre `(ω0 dt)⁶` par pas.

## Pièges rencontrés

- Euler symplectique et Verlet ne sont symplectiques que pour un hamiltonien **séparable** (`H = T(p) + V(q)`) : vrai ici, faux pour le pendule double (M3).
- Pour juger des ordres : même précaution qu'en M1 (instant hors multiple de période, espace des phases).

## Vérification

`testOscillatorExact` (trois régimes, résidu de l'équation), `testOscillatorConvergence`, `testOscillatorEnergy`.

## Taille

18 fichiers, 1823 lignes ajoutées ; `Oscillator.cpp` 114 lignes, `OscillatorModule.cpp` 579 lignes.
