# M4 : gravitation (M4a : Kepler à deux corps, M4b : problème à N corps)

Commits : `7d26dd3` (M4a), `dc9a04b` (M4b et passation M4). Date : 2026-10-08. `--sim 5` (Kepler) et `--sim 6` (N corps).

## M4a : orbite de Kepler

### Objectif

Le premier problème où l'on peut **prédire** l'erreur d'un schéma, pas seulement la mesurer : un schéma symplectique d'ordre 2 fait précesser l'orbite de Kepler d'une quantité qu'on sait calculer.

### Physique

`r'' = -GM r/|r|³`. Solution exacte par l'équation de Kepler `M = E - e sin E` (résolue par Newton, avec le cas particulier `E = 0` pour `M = 0`), puis `x = a (cos E - e)`, `y = b sin E`. Unités normalisées : `GM = a = 1`, `T = 2π`.

### Résultats mesurés (e = 0,5)

- **Équation de Kepler** : résidu ≤ 4,4·10⁻¹⁶ ; 4 à 9 itérations de Newton pour `e ≤ 0,99`, 20 à `e = 0,9999` près du périastre.
- **Rapports d'erreur quand on double le nombre de pas** (mesurés à 0,35 T et dans l'espace des phases) : Euler 1,98 ; Euler symplectique 1,98 ; Verlet 4,00 ; RK4 16,95. (Euler à 2,7 T n'est encore qu'à 1,76 avec 25 600 pas : il est pré-asymptotique.)
- **RK45** (tolérance 10⁻¹⁰) : erreur 8,5·10⁻⁹ en 654 pas sur 2,7 T.
- **100 orbites à 200 pas par orbite** : Euler explicite s'échappe (`E > 0`) ; Euler symplectique : `dE/E` borné à 5,2·10⁻² ; Verlet : 2,7·10⁻³ ; RK4 : 9,9·10⁻⁶ puis 9,3·10⁻⁵ (dissipation lente). Le moment cinétique `L` est conservé à l'arrondi par les schémas symplectiques (test < 10⁻⁹), pas par RK4 (1,4·10⁻⁵).
- **Précession numérique** (rétrograde) : Verlet −0,2232 °/orbite, Euler symplectique −0,2234 (identiques à l'ordre dominant : Lie-Trotter est conjugué à Strang), **RK4 +0,00029, soit 770 fois moins**. En `h²` : −0,886 → −0,0559 de 100 à 400 pas par orbite (rapport 15,85).
- **Formule du hamiltonien modifié** `Δω = −(π/8)(GM h²/a³)(4+e²)/(1−e²)³` : écart à la mesure de 0,02 à 0,2 % pour `e` de 0,1 à 0,9.
- **`e = 0,9`** : à 200 pas par orbite, RK4 perd 753 fois l'énergie ; RK45 prend 698 pas pour 2,3·10⁻⁶ sur 2 orbites, avec un pas 300 fois plus court au périastre qu'à l'apoastre (RK4 à pas fixe : 32 000 pas pour 7,8·10⁻⁶).

### Pièges

- Newton seul échoue en `M = 0` (la racine est au bord du crochet) : cas particulier `E = 0`.
- Mesurer la précession par l'angle du vecteur de Runge-Lenz à des instants fixes est **faux** (il oscille au sein d'une orbite et la période numérique n'est pas `T`) : on détecte les passages au périastre (`PeriapsisTracker`).
- À `e = 0,9` et 200 pas par orbite tous les schémas à pas fixe sont détruits (le périastre n'est pas résolu) : c'est précisément le sujet de RK45.

Tests : `testKeplerEquation`, `testKeplerExact`, `testKeplerConvergence`, `testKeplerEnergyAndPrecession`, `testKeplerPrecessionTheory`, `testKeplerAdaptiveStep`, `testPeriapsisTracker`.

## M4b : problème à N corps (CPU, double)

### Objectif

La référence de précision du futur calcul GPU (M7) et le premier système chaotique **conservatif** à plusieurs corps. Interface unique « force sur chaque particule » (`nbody::accelerations`), découplée de l'intégrateur : c'est elle que le GPU reprend en M7.

### Physique

`a_i = G Σ_j m_j (r_j − r_i) / (|r_j − r_i|² + ε²)^(3/2)` avec adoucissement de Plummer `ε` (dérive de l'énergie `U = −G Σ m_i m_j / √(r_ij² + ε²)` : avec `ε > 0` l'énergie est exactement conservée). La force sur `i` par `j` est l'opposée de celle sur `j` par `i` : l'impulsion est conservée à l'arrondi. Cas de référence : 2 corps (retour à Kepler), triangle de Lagrange, « huit » de Chenciner-Montgomery (données de Simó), amas aléatoire à l'équilibre du viriel.

### Résultats mesurés

- **Huit** : `E = −1,28714199`, `P = L = 0`, fermeture à 8·10⁻⁸ à `T = 6,32591398` (meilleure période retrouvée à 10⁻⁷ près). Une perturbation de 10⁻³ y croît **linéairement** (0,058 à 2 T → 0,275 à 10 T) : le huit est stable.
- **Amplification d'un écart de 10⁻⁹ entre t = 2 et t = 16** : huit ×7,3 ; amas de 6 corps (graine 42) ×6,6·10⁴ (`λ ≈ 0,8`) ; graine 7 : ×5·10⁷.
- **Rapports d'erreur au doublement des pas sur le huit** (0,35 T) : Euler 1,97 ; Euler symplectique 2,00 ; Verlet 4,00 ; RK4 16,1 (ordres 1, 1, 2 et 4). Triangle de Lagrange par RK4 : 8,0·10⁻⁶ (100 pas par période) → 8,3·10⁻¹¹ (1600), ordre 4.
- **Amas de 6 corps, dt = 10⁻³, t = 20** : `|P| ≤ 3·10⁻¹⁵` pour tous les solveurs (forces opposées) ; `|ΔL|` : Euler 3,4·10⁻³, symplectiques ≤ 4·10⁻¹⁵ (exact à l'arrondi), RK4 9,8·10⁻¹¹ (invariant quadratique non conservé) ; `max |dE/E|` : Euler 0,475, symplectiques 2,4·10⁻², Verlet 8,0·10⁻⁴, RK4 1,8·10⁻⁷. Deux corps : accord avec Kepler à 1,5·10⁻¹¹.

### Pièges

- Le triangle de Lagrange à masses égales est **instable** (critère de Routh) : le jumeau perturbé s'en écarte exponentiellement, c'est normal.
- Un ajustement exponentiel (Lyapunov) sur une croissance seulement linéaire donne un **faux `λ`** (0,12 pour le huit, qui est stable) : on décide du chaos sur l'amplification réelle de l'écart (> 1000), pas sur `λ`.
- L'énergie des schémas d'ordre élevé ne reste pas « bornée » de façon lisse sur un amas (les rencontres rapprochées la font sauter) : on compare des maxima, on n'exige pas un plateau.

Tests : `testNBodyAccelerations`, `testNBodyDiagnostics`, `testLagrangeTriangle`, `testTwoBodyMatchesKepler`, `testFigureEight`, `testNBodyConvergence`, `testNBodyConservation`, `testNBodyChaos`.

## Taille

M4a : 10 fichiers, 1307 lignes ; M4b : 12 fichiers, 1351 lignes. `Kepler.cpp` 159 lignes, `NBody.cpp` 281 lignes (avec les sphères de Plummer de M7), `KeplerModule.cpp` 654 lignes, `NBodyModule.cpp` 585 lignes.
