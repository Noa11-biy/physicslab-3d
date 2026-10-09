# M1 : projectile avec frottement, et les cinq solveurs

Commits : `529f722` (solveurs et projectile avec solution analytique), `6cd12b9` (interface de comparaison). Date : 2026-10-08. `--sim 1`.

## Objectif

Le premier phénomène sert surtout à installer **la méthode** : dérivation, solution exacte, plusieurs schémas numériques, mesure de leur erreur. Le cas : un projectile avec frottement linéaire `F = -k v`, dont la solution exacte existe.

## Physique

`m dv/dt = m g - k v` (g vers le bas). Avec `τ = m/k` : la vitesse tend vers la vitesse limite `m g / k`, et
`v(t) = v_lim + (v0 - v_lim) e^(-t/τ)`, `x(t) = x0 + v_lim t + τ (v0 - v_lim)(1 - e^(-t/τ))`.
Sans frottement la chute libre est un polynôme de degré 2 : Verlet et RK4 l'intègrent **exactement**, il n'y a alors rien à mesurer.

## Les cinq solveurs

Euler explicite, Euler symplectique, Verlet des vitesses, RK4, et RK45 de Dormand-Prince adaptatif (`relTol`, `absTol`, statistiques de pas). `step()` renvoie le pas réellement avancé ; `advance(solver, f, t, y, dt)` enchaîne les pas pour couvrir tout `dt` (nécessaire avec un pas adaptatif). Couleurs fixes pour toute l'application (Euler orange, symplectique citron, Verlet violet, RK4 vert, RK45 rose, exacte bleu) et visibilité selon le niveau : niveaux 1-2 un seul « ordinateur », 3-4 Euler contre RK4, 5-6 cases à cocher.

## Résultats mesurés

- **Ordres de convergence** (pente log-log de l'erreur finale en fonction de dt) : Euler 1,0 ; Euler symplectique 1,0 ; Verlet des vitesses 2,0 ; RK4 4,0. RK45 : erreur 9·10⁻¹¹ pour une tolérance de 10⁻¹⁰.
- **Dérive d'énergie d'Euler sans frottement** (chute libre) : `dE = ½ m g² t dt`, vérifiée à 10⁻⁶.

## Pièges rencontrés (valables pour tout le projet)

- Mesurer un ordre **à un instant qui n'est pas un multiple de la période**, et **dans l'espace des phases** : à `t = n T` l'erreur de phase d'ordre 1 s'annule (sur-convergence apparente : Verlet « d'ordre 4 »).
- Un schéma peut être **pré-asymptotique** : RK4 sur un grand angle n'atteint un rapport d'erreur de 16 qu'à 1280 pas ou plus. On ajuste le point de mesure, jamais le seuil du test, sans avoir regardé les valeurs.
- Un test qui passe du premier coup n'est pas forcément bon : toujours afficher les valeurs réelles.

## Vérification

`testFreeFall`, `testProjectileAnalytic`, `testConvergenceOrders`, `testPolynomialExactness`, `testRK45`. Chaque problème expose `rhs()`, l'énergie, la solution exacte et une fonction `xxxError(problem, solver, steps, tEnd)` : le même schéma sert à tous les modules suivants.

## Taille

446 + 469 lignes ajoutées sur les deux commits ; `Projectile.cpp` 80 lignes. À M1 l'interface vivait encore dans `Application.cpp` ; elle a été déplacée dans `ProjectileModule.cpp` (420 lignes aujourd'hui) au moment de la modularisation de M2.

## Limites

Frottement linéaire seulement (le frottement quadratique, plus réaliste à haute vitesse, n'a pas de solution simple : noté en dette).
