# M5 : collisions et frottements (M5a frottement sec, M5b chocs et rebonds, M5c berceau de Newton)

Commits : `07795c2` (M5a), `9171015` (M5b), `e32c7ea` (passation et feuille de route de M5c), `48a1c12` (M5c), `dc7dbde` (passation de fin de M5c). Date : 2026-10-08. `--sim 7`, `--sim 8`, `--sim 9`.

Fil conducteur : **la physique devient discontinue**, et les schémas d'ordre élevé s'effondrent s'ils ignorent la discontinuité. Trois façons de la traiter : détecter l'événement (M5a, M5b), le régulariser (M5a), ou ne pas en avoir parce que la force est continue (M5c, contact de Hertz).

## Brique commune : détection d'événement (`core/Events`)

`advanceToEvent` : instant exact où `g(t, y)` change de signe, par bissection sur le pas. Elle sert à l'arrêt d'un bloc, à l'impact d'une balle, au contact de deux billes.

## M5a : frottement sec de Coulomb sur un plan incliné

Trois modèles numériques : **événement + adhérence** (on arrête exactement quand la vitesse s'annule, puis on décide s'il y a adhérence : `tan θ ≤ μs`), **naïf** (`sgn(v)` dans l'équation différentielle), **régularisé** (`tanh(v/ε)`). Une solution exacte par morceaux sert de référence ; la résistance `k = 0,7` évite que la solution soit un polynôme.

Résultats mesurés :
- Événement + adhérence : on retrouve les ordres des solveurs, rapports d'erreur au doublement des pas Euler 2,0, Euler symplectique 2,0, Verlet 4,0, RK4 16,2.
- **Naïf** : Verlet tombe à 2,0 et RK4 à 2,05 puis 1,94 (16,1 sur le même problème sans frottement sec) : l'ordre 1 pour tous. À 400 pas : erreur 6,7·10⁻² (naïf) contre 3,3·10⁻¹⁰ (événement), soit **2·10⁸ fois mieux**.
- Erreurs finales de l'événement à `dt = 0,05` : RK4 4·10⁻⁸ à 9·10⁻⁸, Verlet 1,4·10⁻³, RK45 10⁻¹⁰ ; sans résistance Verlet et RK4 sont exacts à 2·10⁻¹⁴. Instant d'arrêt trouvé par bissection : 10⁻⁶ (RK4), 10⁻⁸ (RK45).
- **Piège de l'adhérence** (`μd < tan θ ≤ μs`, bloc au repos) : exact = reste collé ; le naïf RK4 est parti à −9,2 m en 5 s (le modèle ignore `μs` : erreur de **modèle**, qui ne tend pas vers 0 avec `dt`). Un bloc qui doit rester collé après un arrêt garde une vitesse résiduelle `1,18 |A| dt`, proportionnelle au pas, jamais nulle.
- Régularisation `tanh(v/ε)` : erreur de position à 5 s 0,017 / 0,085 / 0,35 m pour `ε` = 0,01 / 0,05 / 0,2, **indépendante de `dt`** (vérifié à 0,01 et 0,001).
- **RK45 naïf sur un bloc qui doit rester collé ne termine pas** : l'erreur locale reste d'ordre `h` sur la surface `v = 0`, le pas descend jusqu'à `hMin = 10⁻¹⁴`. D'où `advance(..., maxSteps)`, qui borne le travail (400 pas par pas de calcul dans l'interface, état « bloqué »).

Pièges : un plateau d'erreur indépendant du pas n'est pas forcément un défaut du schéma (l'arrêt tombait au même décalage de 0,18 ms de la frontière de pas pour tous les `dt` multiples de 0,0025 s) : choisir des cas dont l'arrêt n'est pas aligné sur la grille. L'événement doit ignorer un départ exactement sur la surface (`g = 0`) : sinon faux arrêt à chaque redémarrage.

Tests : `testEventDetection`, `testAdvanceBudget`, `testInclineExact`, `testInclineEventDriven`, `testInclineNaive`, `testInclineRegularized`. Taille : 19 fichiers, 1464 lignes ; `Friction.cpp` 160 lignes, `FrictionModule.cpp` 614 lignes.

## M5b : chocs et rebonds

Restitution de Newton (`collide1D`, `collideSpheres` : seule la composante normale est inversée et amortie). Cas : balle rebondissante (accumulation de **Zénon** : une infinité de rebonds en temps fini) et choc de deux billes (disques de rayon 0,5, `e = 1`, masses égales, décalage 0,6 m). Deux modèles : **événement** (l'instant du choc est résolu exactement) et **naïf** (le contact n'est vu qu'après le pas).

Résultats mesurés :
- Balle `e = 0,8` lâchée de 2 m : instants d'impact de la solution exacte égaux à la formule à 1,3·10⁻¹⁵, sommets `e^(2n) h0` à 7·10⁻¹⁶ ; repos à 5,747939132 s contre `t0 (1+e)/(1−e) = 5,747939222 s` (écart −9·10⁻⁸ dû au seuil de 10⁻⁷) ; 81 rebonds (50 avec le seuil par défaut 10⁻⁴), calcul exact en 0,2 ms.
- Modèle événement : sans résistance RK4 et Verlet sont exacts à 4·10⁻¹³ jusqu'après le repos ; avec `k = 0,7` rapports Euler 2,0, Euler symplectique 2,0, Verlet 4,0, RK4 16,1 ; instant de repos retrouvé à 2·10⁻⁶ (`dt = 0,05`), 8·10⁻⁸ (0,02), 2·10⁻⁹ (0,007), 1,5·10⁻¹² (0,001) : ordre 4, 46 rebonds.
- Modèle **naïf** : pente log-log 0,94 pour RK4 contre 4,0 pour l'événement ; à 400 pas 1,5·10⁻² contre 7,3·10⁻¹¹ ; premier impact vu en retard de moins d'un pas (0,6390 s pour 0,63868 s à `dt` = 1 ms) ; **la balle ne s'arrête jamais** : 0,5 s après le repos exact, `|vy|` max = 0,78 / 0,28 / 0,11 m/s pour `dt` = 0,02 / 0,01 / 0,005.
- Deux billes : contact à `t = 1,0667 s`, angle de sortie 90° ; événement : erreur ≤ 8·10⁻¹⁴ pour les quatre solveurs, instant du choc exact à 10⁻¹² ; naïf : **même erreur pour les quatre solveurs** (RK4 compris), proportionnelle au pas (0,43 puis 0,10 puis 0,025 quand le pas est divisé par 4), choc vu 3,3·10⁻³ s trop tard à `dt = 0,01` ; l'impulsion reste conservée même en naïf.

Pièges : une fonction d'événement qui vaut **exactement 0** au départ (balle posée sur le sol juste après un rebond) ne voit pas un vol entier contenu dans un pas : la balle traverse le sol (observé : `y = −26 m` vers `t = 3,7 s`). Remèdes : `g = max(y, 0)` puis un petit pas de décollage ; même idée pour les billes (`g = max(distance, 0)`). Le seuil d'arrêt doit rester ≥ 10⁻⁸ m/s (en dessous, `landingTime` renvoie 0). Ne pas asserter le nombre de rebonds du modèle naïf (24, 52, 99 selon `dt`) ni les ordres mesurés **après** le repos.

Tests : `testCollide1D`, `testCollideSpheres`, `testTwoBallExact`, `testBounceExact`, `testBounceEventDriven`, `testBounceNaive`, `testTwoBallContact`. Taille : 13 fichiers, 1777 lignes ; `Collision.cpp` 142 lignes, `Bounce.cpp` 186 lignes, `BounceModule.cpp` 705 lignes.

## M5c : berceau de Newton, contact de Hertz contre impulsions séquentielles

Deux philosophies. **Hertz** : la force de contact `F = k δ^(3/2)` (avec amortissement de Hunt-Crossley en option) est continue, donc intégrable sans événement ; le résultat d'une chaîne de N billes émerge de la dynamique. **Impulsions séquentielles** : chaque choc binaire est instantané, résolu dans un ordre ; avec `e = 1` et masses égales chaque choc échange les vitesses.

Résultats mesurés, unités normalisées (`m = 1`, `R = 1/2`, `k = 10⁴`, `v = 1`) :
- Constante du temps de contact `(4/5) Γ(2/5) Γ(1/2) / Γ(9/10) = 2,943275` ; deux billes (RK45 serré) : `δ_max = 0,0208138302` et `T = 0,0612608299` retrouvés à 6·10⁻¹³ et 1,6·10⁻¹³ (relatif), sortie = échange des vitesses à 5·10⁻¹³.
- **Le résultat d'une chaîne ne dépend ni de `k` ni de `v`** (une seule échelle de longueur `(m v²/k)^(2/5)`) : 3 billes avec `(k, v)` = (10⁴, 1), (10⁶, 1), (10⁴, 3), (2·10³, 0,2) donnent les mêmes six chiffres ; seule la durée change (rapport 6,3 = 100^(2/5) pour `k` × 100).
- **Hertz, 3 billes, 1 lancée** : (−0,070952 ; +0,076403 ; +0,994549) : les billes de tête **reculent** d'environ 0,07 v, contrairement au modèle d'impulsions (« une entre, une sort » n'est pas exact). La dernière bille part à 0,9945 / 0,9910 / 0,9888 / 0,9875 / 0,9867 pour N = 3 à 7. Deux lancées sur cinq : (−0,112615 ; −0,041960 ; +0,214486 ; +0,800367 ; +1,139722). Symétrie exacte (miroir + changement de repère galiléen) entre `n` lancées et `N−n` lancées : écart ≤ 10⁻⁸.
- **Impulsions séquentielles** : `n` entrent, `n` sortent. **L'ordre de résolution ne change rien dans le berceau, même avec `e < 1`** : écart exactement 0 sur 245 cas (N = 3 à 9). Il compte seulement pour une bille prise entre deux voisines qui s'approchent déjà (écart 0,070 pour `e = 0,5`). Avec `e < 1` aucune bille ne recule jamais (différence qualitative avec Hertz).
- **Ordres des schémas sur le contact** (contact décalé de 0,0123 pour ne pas tomber sur la grille) : Euler 1 (rapport 4,07), Verlet 2 (15,4), Euler symplectique = Verlet à 5·10⁻¹³ (conjugués ; la mesure est faite hors contact), **RK4 n'atteint pas l'ordre 4** : ordre apparent 2,29 à 2,75 selon la phase du contact sur la grille, cause probable (non testée avec une force lisse) : `k δ^(3/2)` n'a pas de dérivée seconde en `δ = 0`. RK45 par défaut (tol 10⁻⁸) : erreur 2,4·10⁻⁷ à 3,4·10⁻⁷, 88 pas acceptés, 26 refusés.
- **Stabilité** (3 lancées sur 7) : Euler explicite **diverge** à 10 pas par contact ; Euler symplectique et Verlet : `E/E0 = 1,0004` à 10 pas, 1,0000 à 30 pas. **Hunt-Crossley** `F = k δ^(3/2)(1 + (3/2) α δ')` : `e ≈ 1/(1 + α v)` (0,990099 / 0,952370 / 0,909016 / 0,768 pour `α` = 0,01 / 0,05 / 0,1 / 0,3).

Corrections apportées par la mesure au brouillon de la feuille de route : (1) l'ordre de résolution des impulsions « compte » dans le berceau : **faux**, il est sans effet ; (2) RK4 serait d'ordre 4 sur Hertz : **faux**, 2,3 à 2,75.

Pièges : ne jamais asserter « ×16 » ni « ×5,6 » sur Hertz (utiliser une pente sur 200 → 3200 pas et des bornes larges) ; mesurer `δ_max` et `T` par événement (vitesse relative nulle, puis compression nulle), jamais en échantillonnant le maximum aux pas de RK45 ; la fin de la collision est `collisionOver` : plus aucune compression **et** vitesses rangées. Un schéma qui diverge met des valeurs infinies dans l'état : détecter, arrêter, ne pas dessiner.

Tests : `testHertzReference`, `testHertzTwoBalls`, `testCradleInvariants`, `testThreeBallFamily`, `testHertzThreeBalls`, `testHertzChain`, `testSequentialImpulses`, `testSequentialOrder`, `testCradleGap`, `testCradleConvergence`, `testHertzDamping`. Taille : 10 fichiers, 1446 lignes ; `Cradle.cpp` 183 lignes, `CradleModule.cpp` 610 lignes.

## Limites de M5

Chaîne 1D seulement (pas de choc oblique) ; billes dessinées en fil de fer ; billes non suspendues (pas de force de rappel `g/L`) ; frottement de roulement et frottement visqueux quadratique absents ; un seul modèle d'amortissement (Hunt-Crossley).
