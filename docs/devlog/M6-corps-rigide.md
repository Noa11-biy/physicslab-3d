# M6 : corps rigide (solides libres, toupie de Lagrange, intégrateurs d'orientation)

Commit : `5935a4f` (suivi de `e0380f7` pour la passation). Date : 2026-10-09. `--sim 10`.

## Objectif

Passer de la translation à la **rotation**, où l'état vit sur une variété (le groupe des rotations) et non dans un espace plat : un quaternion qui perd sa norme, ou un schéma qui perd le moment cinétique, trahit le calcul. Trois situations : solide libre symétrique, solide libre asymétrique (instabilité de l'axe intermédiaire), toupie de Lagrange (corps pesant fixé en un point).

## Physique

Équations d'Euler dans le repère du corps : `I ω' + ω × (I ω) = τ`, couplées à `q' = ½ q ⊗ (0, ω)` (convention : `q` envoie le repère du corps vers le repère fixe, `ω` exprimé dans le corps). Invariants du solide libre : énergie `E` et vecteur moment cinétique `L` (fixe dans l'espace). Tenseurs d'inertie en formules fermées (sphère, boîte, cylindre).

- **Solide symétrique libre** (`I1 = I2`) : `q(t) = q_L(L t / I1) q0 q_3(−Ω t)` avec `Ω = (I3 − I1) ω3 / I1`.
- **Solide asymétrique** : solution par les fonctions elliptiques de Jacobi (séries de Fourier en nome `q`, DLMF 22.11), formules de Landau-Lifchitz § 37, pour un départ `(a, 0, c)`.
- **Toupie de Lagrange** (`m = 1`, `I = (1,2 ; 1,2 ; 0,4)` autour du pivot, `l = 0,5`) : précession régulière exacte, nutation par la cubique `u'² = f(u)` (racines `u1`, `u2`, `u3`), spin critique de la toupie « endormie » `2 √(I1 m g l) / I3`.

## Intégrateurs d'orientation comparés

Euler explicite, RK4 + renormalisation du quaternion, **Lie-Heun** (groupe de Lie, norme exacte), **découpage symplectique** (ordre 2, conserve `L` par construction).

## Résultats mesurés

- Tenseurs d'inertie : formules fermées contre une intégration directe sur grille, à 10⁻³ (boîte) et 5·10⁻³ (sphère, cylindre : erreur de bord en 1/n).
- Solide symétrique : accord avec RK45 serré à 10⁻¹⁴ … 2·10⁻¹² sur `ω` et 7·10⁻¹⁴ … 4·10⁻¹² sur l'orientation (t = 3 et 25) ; avec le signe de la rotation propre faux le même test donne 0,56 … 2,8 : le test est sensible. Invariants le long de la trajectoire (`I = (1, 2, 3)`) : `E` à 2·10⁻¹³ et `L` à 10⁻¹¹ jusqu'à t = 100.
- **Fonctions de Jacobi** (séries en nome) contre RK4 : 10⁻¹⁰ pour `m` de 10⁻¹⁴ à 0,999 ; identités `sn² + cn² = 1` à 4·10⁻¹⁵, jusqu'à `m = 1 − 10⁻¹⁰`.
- Solide asymétrique : périodes `4K/λ` = 8,62606 / 5,72746 / 11,08236 s ; `|ω_exact − ω_RK45| ≤ 5,7·10⁻¹³` pour t ≤ 40. **Axe intermédiaire** : `λ` prévu 0,577350, mesuré 0,577411 (écart relatif 10⁻⁴) ; axes extrêmes stables (perturbation 1,4·10⁻³ → au plus 2,0·10⁻³ sur 100 s).
- **Ordres** (`I = (1,2,3)`, `ω0 = (0,9 ; 0,5 ; 1,1)`, t = 6,3) : Euler ordre 1 (×2,01), RK4 + renormalisation ordre 4 (×16,06 ; erreur 3,5·10⁻¹¹ à 1600 pas), Lie-Heun ordre 2 (×3,99), découpage ordre 2 (×4,00) et environ 3 fois plus précis que Heun à pas égal. Euler explicite sur une rotation pure : `|q|ᴺ = (1 + h² ω²/4)^(N/2)` exactement.
- **Dérives sur 1000 s** à `h = 0,05` (20 000 pas, proche de l'axe instable) : découpage `max |ΔL|/L = 6,2·10⁻¹⁴`, `|q| − 1 = 3·10⁻¹⁶`, `max |ΔE/E| = 6,9·10⁻⁵` (en `h²`, **aucune dérive**) ; RK4 + renormalisation : `ΔE/E` final −3,1·10⁻⁶, `ΔL/L` 1,8·10⁻⁶ ; RK4 sans renormalisation : `|q| − 1 = 4·10⁻⁶` ; **Lie-Heun conserve `|q| = 1` mais dérive** (`ΔE/E` +1,4·10⁻², `ΔL/L` 7·10⁻³) et diverge à `h = 0,2` ; Euler diverge au pas 2251.
- **Toupie** : `E`, `L_z`, `L_3` conservés à 10⁻¹⁰ ; précession régulière contre RK45 à 10⁻⁸ ; exemple `θ = 0,6`, `ω3 = 25` : `φ' = 0,517 rad/s` (un tour en 12,2 s). Nutation (`θ0 = 0,9`) : `u1 = 0,621609968` (= cos 0,9), `u2 = 0,926152996`, période 1,069185 s, `u(T) − u0 = 1,6·10⁻¹⁴`. Spin critique 12,12847 rad/s ; au-dessus (1,15 / 1,5 / 3 fois) une inclinaison de 0,01 reste ≤ 0,0203 / 0,0134 / 0,0106 ; en dessous (0,7 fois) la croissance est `e^(γt)`, `γ` prévu 1,44358, mesuré 1,44357.
- **Intégrateurs avec couple** (nutation, 100 s à `h = 0,005`) : découpage `ΔL_z/L_z = 2,8·10⁻¹³` (les « coups » du couple n'ont pas de composante verticale), `ΔE/E` 3,9·10⁻⁵ borné ; RK4 `ΔE/E` 7,5·10⁻⁷ (à ce pas il conserve mieux `E`) ; Heun 2,3·10⁻³ ; Euler diverge.

## Correction apportée par la mesure

Hypothèse du brouillon : « une toupie sous le spin critique tombe ». **Faux** : de 10⁻⁶ elle monte jusqu'à 1,2457 rad (71°) puis **remonte** (`θ(20 s) ≈ 10⁻⁵`) : grande nutation d'un système sans perte, fixé au pivot. Pour qu'elle tombe vraiment il faudrait une table (contact, frottement), non écrite.

Autre idée reçue testée : « garder `|q| = 1` suffit ». **Faux** : Lie-Heun l'assure et dérive en énergie et en `L` ; c'est la structure symplectique qui conserve `L`.

## Pièges

- `Quaternion::rotate` n'est valable que pour un quaternion **unitaire** : pour dessiner un état dont la norme a dérivé, utiliser le produit complet `q v q*` (le solide gonfle de `|q|²`, c'est l'erreur à montrer) ; pour les diagnostics (`E`, `L`) normaliser d'abord.
- `ellipticK(√(1 − m))` perd la précision pour `m` petit (nome faux de 0,5 % à `m = 10⁻¹⁴`) : le nome se forme avec `K = π / (2 AGM(1, k'))` et `K' = π / (2 AGM(1, k))` directement.
- Appeler `reference(t)` dans une boucle sur `t` ré-intègre depuis 0 à chaque fois : coût **quadratique** (un test a pris 76 s en Debug). Intégrer une fois et échantillonner.
- La solution de Jacobi exige `I1 < I2 < I3` et un départ `(a, 0, c)` avec `a, c > 0` ; sans `ω2(0) = 0` il faudrait l'intégrale elliptique incomplète de première espèce (non écrite).
- Un test de l'« ordre » de Euler n'est asymptotique qu'à pas fin (avec le couple de la toupie, 50 pas donnent NaN).
- Interface : le symbole ⊥ ne s'affiche pas dans la police Segoe UI (Ω, ⁻, ² et ³ passent).

## Vérification

`testInertiaTensors`, `testEulerEquations`, `testFreeBodyInvariants`, `testSymmetricTopExact`, `testJacobiFunctions`, `testAsymmetricExact`, `testIntermediateAxisInstability`, `testRotationIntegratorBasics`, `testSplittingConservation`, `testRotationOrders`, `testRotationDrift`, `testHeavyTopInvariants`, `testHeavyTopSteadyPrecession`, `testNutationPeriod`, `testSleepingTop`, `testRotationIntegratorsWithTorque`.

## Taille et limites

10 fichiers, 2235 lignes ; `RigidBody.cpp` 361 lignes, `RigidBodyModule.cpp` 808 lignes (le plus gros module). Limites : boîte de dimensions fixes (3 × 2 × 1) ; pas de solution exacte d'orientation du solide asymétrique (fonctions thêta) ; toupie à pivot fixe seulement ; pas de choc de corps rigides (cône de Coulomb, paradoxe de Painlevé) ; découpage d'ordre 2 seulement ; corps dessinés en fil de fer.
