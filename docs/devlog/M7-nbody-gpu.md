# M7 : N corps sur carte graphique (compute shader), et écart CPU / GPU

Commits : `c8c4e3b` (étapes 1 et 2 : accélérations puis kick-drift-kick sur GPU, `--gpu-test`), `a7ee166` (étape 3 : module `--sim 11`, sphères de Plummer, passation). Date : 2026-10-09. Matériel de toutes les mesures : **Intel UHD Graphics** (carte intégrée), OpenGL 4.5.0, pilote 31.0.101.3962, 1024 threads par groupe au plus, 32 Kio de mémoire partagée par groupe.

## Objectif

Reprendre le calcul O(N²) de `nbody::accelerations` (M4b, CPU `double`) dans un compute shader en `float` et **mesurer l'écart** : « le CPU vérifie le GPU ». Contrat : la même interface (positions, masses, n, G, adoucissement) que le CPU.

## Conception

- **Un thread par étoile cible.** Les sources sont lues par tuiles de 256 corps chargées en mémoire partagée (chaque étoile n'est lue qu'une fois par groupe au lieu d'une fois par thread), accumulation de `a_i = G Σ_j m_j (r_j − r_i)/(|r|² + ε²)^(3/2)`. Cases de remplissage de masse nulle pour la dernière tuile.
- **Centrage sur le barycentre en `double` avant la conversion en `float`** : les accélérations ne dépendent pas de l'origine, mais l'arrondi des positions, lui, en dépend.
- **Intégration kick-drift-kick sur GPU** (Verlet des vitesses, ordre 2, symplectique) : positions, vitesses et accélérations restent dans des SSBO. Deux demi-coups consécutifs se fondent en un coup entier : par pas, un envoi d'accélérations et une mise à jour élément par élément.
- **Énergie potentielle dans le même shader** : `φ_i` rangé dans la quatrième composante des accélérations (le terme `m_j/r` est déjà calculé), `U = ½ Σ m_i φ_i` sommée en double. Le suivi de l'énergie ne coûte donc plus de O(N²) sur le CPU.
- **Précision** : `float` ou `double` choisis par un `#define` injecté à la compilation (même source de shader). Taille du groupe de travail réglable.
- **Envois courts** : le calcul est découpé en envois d'environ 50 ms, dimensionnés sur le débit mesuré (voir « TDR » plus bas).

## Résultats mesurés

Précision des accélérations à état fixé (GPU `float` contre CPU `double`, écart global `‖Δa‖/‖a‖`) : 1,6·10⁻⁸ (huit, N = 3), 1,3·10⁻⁸ (Lagrange), 2,7·10⁻⁷ (N = 100), 4,7·10⁻⁷ (1000), 1,05·10⁻⁶ (5000), pire étoile 3,4·10⁻⁶ à N = 5000. **Environ `5·10⁻⁸ √N`** : une somme de N termes arrondis donne `u √N` avec `u = 6·10⁻⁸`. Le GPU `double` coïncide à 3·10⁻¹⁶ (N ≤ 1000). La force totale relative `|Σ m a| / Σ m|a|` vaut 1 à 2·10⁻⁸ en float contre 10⁻¹⁷ au CPU : la 3ᵉ loi n'est vraie qu'à l'arrondi près.

**Centrage** : amas de 1000 corps décalé de (1000, −2000, 500) : écart 4,7·10⁻⁷ avec centrage contre 1,6·10⁻⁴ sans, soit un facteur **345** ; sans centrage l'erreur vient de l'arrondi des positions et non du calcul.

**Taille du groupe** (N = 16384, meilleur de 8 envois) : 32 → 17,7 ms ; 64 → 10,0 ; 128 → 9,2 ; 256 → 8,35 ; 512 → 8,33 ; 1024 → 9,0. On retient 256.

**Débit** (interactions ordonnées `N²` par seconde) : GPU float 5 Ginteractions/s à N = 1000, 17,7 à 4000, 28 à 8000, 33 à 32 000, 30 à 100 000 (0,34 s), 21 à 200 000 (1,9 s). CPU `double` en Release : 0,33 Ginteractions/s (0,11 à 0,14 les jours où la machine est moins performante). À N = 16 000 : GPU 8,1 à 13 ms contre CPU 0,77 à 1,98 s, soit **80 à 150 fois plus rapide**.

**Le `double` du shader est émulé par le pilote** : il est accepté (et exact : 3·10⁻¹⁶) mais tourne à 2,9·10⁶ interactions/s (0,35 s à N = 1000, 1,4 s à N = 2000, 22 s à N = 8000), de mille à plusieurs milliers de fois moins que le float et 100 fois moins que le CPU. Il ne sert qu'à valider la logique du shader.

**Intégration** (GPU contre Verlet CPU, même pas) : en `double` le GPU coïncide avec le CPU à 7,2·10⁻¹⁶ après 20 pas du huit (8,3·10⁻¹⁷ entre un calcul en un lot et un calcul en lots plus petits) ; en `float` 2,2·10⁻⁷. Dérives maximales GPU float (CPU double entre parenthèses) :

| Cas | Énergie | Impulsion | Moment cinétique |
|---|---|---|---|
| Huit, 5 périodes (10 000 pas) | 2,0·10⁻⁶ (1,4·10⁻¹⁰) | 2,2·10⁻⁶ (2·10⁻¹⁵) | 1,3·10⁻⁶ (10⁻¹⁴) |
| Amas N = 100, t = 20 (20 000 pas) | 1,0·10⁻⁶ (7,9·10⁻⁷) | 4,9·10⁻⁷ (6·10⁻¹⁶) | 2,9·10⁻⁶ (3·10⁻¹⁵) |
| Amas N = 500, t = 1 | 4,1·10⁻⁷ (3,8·10⁻⁷) | 1,8·10⁻⁸ | 3,3·10⁻⁸ |
| Amas N = 5000, t = 1 | 2,5·10⁻⁷ | 5,5·10⁻⁹ | 4,9·10⁻⁹ |

L'énergie dérive autant que l'**erreur propre du schéma** d'ordre 2 (le float n'y ajoute presque rien). L'impulsion et le moment cinétique, exacts à l'arrondi en double, dérivent de 5·10⁻⁹ à 3·10⁻⁶ en float.

**Écart de trajectoire** : sur le huit (régulier) l'écart GPU/CPU croît **linéairement** (4,3·10⁻⁵ après 1 période, 2,0·10⁻⁴ après 5) : le bruit d'arrondi des mises à jour de `v` et `x` est amplifié par le cisaillement de l'orbite. Sur l'amas N = 100 (chaotique, `λ ≈ 1`) il croît exponentiellement : 3,9·10⁻⁵ (t = 2), 2,7·10⁻⁴ (4), 2,9·10⁻³ (6), 2,5·10⁻² (8), 0,16 (10), 1,1 (12), 4,0 (14), 9,0 (16), 13,8 (20, saturé). C'est pourquoi on compare les accélérations **à état fixé** puis des horizons courts, jamais des trajectoires longues.

**Temps par pas** (une évaluation de force par pas côté GPU) : 53 à 89 µs à N = 3 et 100 (CPU : 60 à 152 µs), 181 à 283 µs à N = 500 (CPU 1,5 à 3,8 ms), 1,5 ms à N = 5000. Avant de regrouper les lectures de requêtes de temps : 308, 717 et 1920 µs.

**Sphère de Plummer** (`NBodyProblem::plummer`, méthode d'Aarseth, Hénon et Wielen, rayons limités à 10 a) : N = 4000, `ε = 0` : `E = −0,1505` (7 graines : −0,147 à −0,159, moyenne −0,152 contre `−3π/64 = −0,1473` : couper les rayons au-delà de 10 a creuse un peu le puits), `2T/|U| = 0,995`, rayon de demi-masse 1,307 (théorie 1,3048). Changement d'échelle exact (`a = 2` : `E/2`, rayons ×2), équilibre tenu en dynamique (N = 500, t = 5 : rayon de demi-masse +3,4 %, `dE/E` 4·10⁻⁷), collision de deux sphères lié (`E = −0,148`, `P = 10⁻¹⁶`).

**Module `--sim 11`** (N = 2000, build Debug) : GPU 1,5 ms par pas dont 1,2 ms de calcul pur (le reste : attente, la carte dessine aussi), 655 pas/s ; mode « les deux » : CPU 63,6 ms par pas (Debug) contre 2,5 ms, soit un rapport 25 ; écart GPU − CPU 2,9·10⁻⁵ après quelques unités de temps ; histogramme d'erreur centré sur 5·10⁻⁷ (maximum 2,1·10⁻⁶).

## Ce que la mesure a corrigé

- Le brouillon prévoyait une dérive d'impulsion de l'ordre de 10⁻⁷ **par pas** : **faux**, bien moins (5·10⁻⁹ à 3·10⁻⁶ cumulés).
- L'idée d'un GPU `double` pour séparer l'arrondi du parallélisme : réalisable mais **inutilisable** sur cette carte (émulé).
- La dérive d'énergie « plus forte qu'en double » : vraie seulement pour P et L ; l'énergie suit le schéma.

## Incidents, dans l'ordre où ils sont arrivés

1. **TDR : un seul envoi GPU de plus de ~2 s peut geler le pilote.** Le `double` à N = 5000, d'un seul envoi, a bloqué le test plus de deux minutes. Remède : envois découpés (~50 ms), taille calculée sur le débit mesuré (un seul groupe au premier envoi). Le test limite maintenant le `double` à N ≤ 1000.
2. **Une lecture bloquante de la requête de temps à chaque envoi coûtait 130 à 350 µs par pas** : à N = 100 le GPU perdait contre le CPU (308 µs contre 157 µs). Remède : empiler les envois, relire les temps par lots (au plus 0,25 s de travail empilé, 256 requêtes) : 89 µs.
3. **NaN en ε = 0 (régression que j'avais introduite).** En ajoutant le potentiel j'ai écrit le terme d'une étoile sur elle-même hors de la garde : `0 · ∞ = NaN`. `--gpu-test` l'a signalé aussitôt (accélérations non finies sur le huit et le triangle de Lagrange), ce qui justifie ses cas à ε = 0.
4. **Les étoiles disparaissaient** dès que l'on changeait N ou lançait le banc d'essai (grille intacte, calcul et données sains, aucune erreur GL). Deux hypothèses fausses (mon calcul ; la création de tampons pendant la trame ImGui). La dichotomie (même `setBodies` seul le déclenchait) a désigné le **`Renderer` partagé**, qui ré-allouait son VBO à chaque dessin (`glNamedBufferData`) : le pilote Intel le supporte mal dès que d'autres tampons existent. Remède : tampon de flux persistant à décalage (4 Mio au moins, décalage remis à 0 par image). M4b et M6 vérifiés sur captures.

## Validation de la validation

`--gpu-test` doit échouer quand le shader est faux : injection volontaire de fautes (facteur 1,0001 sur la force, facteur 1,002 sur le coup de pied, remplissage de tuile avec une masse non nulle) : détectées à chaque fois (code de sortie 1 ; 35 échecs pour le facteur 1,0001 sur la force, sept contrôles différents en échec pour le facteur 1,002 sur le coup de pied), puis fautes retirées. Les shaders étant lus à l'exécution, l'essai ne demande pas de recompilation.

## Vérification

- `testPlummer` (dans `tests/test_core.cpp`).
- `physicslab --gpu-test [--gpu-max-n N]` : 1. précision à état fixé, 2. tailles de groupe, 3. temps et débit, 4. intégration. Code de sortie non nul au-delà des seuils (`5·10⁻⁸ (3 + √N)` sur l'écart global float, 10⁻¹² en double, 3·10⁻⁵ sur les dérives, etc.). Option CMake `PHYSICSLAB_GPU_TESTS` pour l'ajouter à `ctest` (désactivée : exige un pilote OpenGL 4.5).
- `physicslab --smoke-test --sim 11 --level N` pour les six niveaux ; captures `tools/screenshot.ps1` aux niveaux 1, 3, 4, 5 et 6.

## Taille

11 fichiers (1154 lignes) puis 17 fichiers (1046 lignes). `GpuNBody.cpp` 404 lignes, `GpuTest.cpp` 539 lignes, `GpuNBodyModule.cpp` 672 lignes, shaders 59 et 37 lignes.

## Limites

Mises à jour de `v` et `x` en float sans sommation compensée ; pas de relecture asynchrone (les positions repassent par le CPU pour le dessin) ; pas d'arbre de Barnes-Hut (O(N²)) ; pas de disque galactique dans le module ; N ≤ 20 000 sur GPU et 4000 sur CPU ; le banc d'essai bloque l'interface quelques secondes ; énergie estimée en float par le shader au-delà de 4000 étoiles ; textes des niveaux 2 et 3 non relus à l'écran.
