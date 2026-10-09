# M3 : pendule simple (solution elliptique exacte) puis pendule double (chaos)

Commits : `62a34a5` (physique), `b393ecd` (interface). Date : 2026-10-08. `--sim 3` (simple) et `--sim 4` (double).

## Objectif

Sortir du linéaire. Le pendule simple a une solution exacte **non élémentaire** (fonctions elliptiques) : excellente occasion de montrer que « exacte » ne veut pas dire « formule à une ligne ». Le pendule double n'a aucune solution fermée et sert à introduire le **chaos** et son exposant de Lyapunov.

## Pendule simple

`θ'' = -(g/L) sin θ`. Période exacte : `T = 4 K(k)/ω0` avec `k = sin(θ0/2)` et `K` l'intégrale elliptique complète de première espèce, calculée par moyenne arithmético-géométrique (`ellipticK`). La solution `θ(t)` s'écrit avec les fonctions de Jacobi `sn`, `cn`, `dn`.

- `T/T0 = 1,0732` à 60° (exact).
- Solution elliptique validée contre un RK45 serré (`relTol` 1e-13) : écart de l'ordre de 10⁻¹².

## Pendule double

Équations de Lagrange à deux degrés de liberté (angles `θ1`, `θ2`). Cas de référence : `m1 = m2 = l1 = l2 = 1`, `θ1 = 120°`, `θ2 = -10°`.

- **Exposant de Lyapunov** : `λ ≈ 1,1 /s`. Un écart initial de 10⁻⁹ rad est amplifié de **×4400 en 8 s** (la fonction générique `lyapunovExponent` ajuste une exponentielle sur l'écart de deux trajectoires voisines).
- **Horizon de prédictibilité** (premier instant où l'écart dépasse 0,1) à `dt = 1/60 s` : Euler 0,43 s, Euler symplectique 0,52 s, Verlet 3 s, RK4 8,6 s. C'est le résultat pédagogique central : le choix du schéma décide de combien de temps on peut croire la simulation.

## Pièges rencontrés

- Euler symplectique et Verlet ne sont symplectiques que pour un hamiltonien **séparable** ; le pendule double ne l'est pas (la masse effective dépend des angles). Ils gardent leur ordre mais pas leur garantie sur l'énergie.
- Ajuster une exponentielle (Lyapunov) sur une croissance seulement linéaire donne un faux `λ` : cette erreur est revenue en M4b (voir `M4-gravitation.md`) et c'est la raison de la règle « décider du chaos sur l'amplification réelle de l'écart ».

## Vérification

`testElliptic`, `testPendulumExact`, `testPendulumConvergence`, `testPendulumEnergy`, `testDoublePendulumEquations`, `testDoublePendulumNumerics`, `testChaosVersusRegular` (le pendule double amplifie un écart de 10⁻⁹ d'au moins cinq ordres de grandeur avec `λ` entre 0,5 et 5, le pendule simple régulier non), `testLyapunovFit`.

## Taille

619 + 1313 lignes ajoutées ; `Pendulum.cpp` 118 lignes, `DoublePendulum.cpp` 86 lignes, `PendulumModule.cpp` 508 lignes, `DoublePendulumModule.cpp` 497 lignes.

## Suite immédiate

Le commit `efd75ef` (« fichier de passation et script de capture d'écran ») date de cette période : `docs/PASSATION.md` pour reprendre le projet dans une nouvelle conversation, et `tools/screenshot.ps1` pour vérifier l'interface sur captures.
