# Devlog du domaine Mécanique (M0 à M7)

Journal de développement de PhysicsLab 3D pour tout le domaine **Mécanique**, écrit en fin de domaine (2026-10-09). Un fichier par module ; ce fichier donne la vue d'ensemble, la méthode, ce que les mesures ont **corrigé**, et comment tout revérifier.

Provenance des chiffres : ils viennent de l'exécution réelle, relevés module par module (valeurs de `docs/PASSATION.md`, section 2) et, pour M7, de la session de développement du 2026-10-09 (journaux de `--gpu-test`). Les bornes correspondantes sont **vérifiées en continu** par `tests/test_core.cpp` (76 fonctions de test, 3311 lignes, tout passe). Au moment d'écrire ce devlog quatre valeurs clés ont été recalculées sur le code actuel et coïncident : `T/T0 = 1,0732` du pendule à 60°, l'énergie du huit `−1,28714199` (période 6,32591398), la constante du temps de contact de Hertz `2,943275`, le spin critique de la toupie `12,12847 rad/s`.

## Les modules

| Module | Contenu | `--sim` | Commits | Fichier |
|---|---|---|---|---|
| M0 | socle : CMake, fenêtre OpenGL 4.5, ImGui docking, maths, interface `Solver`, 6 niveaux | - | `8ce9ea3` | [M0-socle.md](M0-socle.md) |
| M1 | projectile avec frottement, cinq solveurs, solution exacte | 1 | `529f722`, `6cd12b9` | [M1-projectile-et-solveurs.md](M1-projectile-et-solveurs.md) |
| M2 | ressort-masse libre, amorti, forcé ; application modulaire | 2 | `4cc1c7b` | [M2-ressort-masse.md](M2-ressort-masse.md) |
| M3 | pendule simple (solution elliptique), pendule double (chaos) | 3, 4 | `62a34a5`, `b393ecd` | [M3-pendules.md](M3-pendules.md) |
| M4 | gravitation : Kepler (précession numérique prédite), N corps (huit, Lagrange, amas) | 5, 6 | `7d26dd3`, `dc9a04b` | [M4-gravitation.md](M4-gravitation.md) |
| M5 | frottement sec, chocs et rebonds, berceau de Newton | 7, 8, 9 | `07795c2`, `9171015`, `48a1c12` | [M5-collisions-et-frottements.md](M5-collisions-et-frottements.md) |
| M6 | corps rigide : solides libres, toupie de Lagrange, intégrateurs d'orientation | 10 | `5935a4f` | [M6-corps-rigide.md](M6-corps-rigide.md) |
| M7 | N corps sur GPU (compute shader), écart CPU / GPU | 11 | `c8c4e3b`, `a7ee166` | [M7-nbody-gpu.md](M7-nbody-gpu.md) |

## Chronologie

- **2026-10-08** : initialisation, licence MIT, dépendances, M0 à M5c (dont le fichier de passation `efd75ef` après M3). Dix-sept commits de code et de documentation.
- **2026-10-09** : M6 (corps rigide) puis M7 (N corps sur GPU), chacun suivi d'une mise à jour de la passation. Le devlog et la suite des livrables de fin de domaine viennent ensuite.

Le dépôt compte 21 commits, 13 493 lignes d'en-têtes, de sources et de shaders (dont 11 modules d'interface), 3311 lignes de tests, 11 simulations, 6 niveaux pédagogiques, 5 solveurs à pas fixe ou adaptatif.

## La méthode, appliquée à chaque phénomène

1. **Dérivation théorique** à la main.
2. **Solution analytique** d'un cas de référence (ou référence RK45 à `relTol` 10⁻¹³ quand il n'y en a pas).
3. **Discrétisation** numérique : plusieurs schémas, comparés.
4. **Implémentation** CPU `double`, puis GPU `float` quand le cas est lourd (seulement M7).
5. **Validation** : invariants, comparaison à l'analytique, ordres de convergence, écart CPU / GPU. Les **tests précèdent l'interface**.
6. **Déclinaison pédagogique** sur les 6 niveaux : l'explication, les paramètres exposés, les équations et le détail des courbes changent avec le niveau.

Règles de travail qui ont payé : ne jamais afficher une mesure sans l'avoir vérifiée ; afficher les valeurs réelles avant d'être satisfait d'un test qui passe ; ajuster le point de mesure, jamais le seuil d'un test ; une feuille de route est un **brouillon à valider par le calcul**.

## Ordres de convergence mesurés (pente log-log de l'erreur en fonction de dt)

| Schéma | Ordre | Remarques |
|---|---|---|
| Euler explicite | 1,0 | crée de l'énergie : multiplie l'énergie de l'oscillateur par `1 + ω0² dt²` à chaque pas (M2) |
| Euler symplectique | 1,0 | énergie bornée sur un hamiltonien séparable ; identique à Verlet à un décalage de demi-pas près |
| Verlet des vitesses | 2,0 | symplectique ; fait précesser l'orbite de Kepler de −0,2232 °/orbite à 200 pas/orbite (M4a) |
| RK4 | 4,0 (4,03 ressort, 3,98 pendule) | précis mais non symplectique : dissipation lente, perd `L` ; **ordre 2,3 à 2,75 seulement sur un contact de Hertz** (force non lisse) |
| RK45 adaptatif | erreur 9·10⁻¹¹ pour une tolérance 10⁻¹⁰ | indispensable près d'un périastre (pas 300 fois plus court qu'à l'apoastre) ; **ne termine pas** sur un modèle de frottement sec naïf sans budget |

Dès que la physique est **discontinue** (frottement sec, rebonds), les ordres s'effondrent à 1 pour tous les schémas si l'on ignore l'événement : 2·10⁸ fois moins précis à 400 pas pour le frottement sec (M5a), et la balle « ne s'arrête jamais » (M5b).

## Ce que la mesure a corrigé

Presque chaque module a démenti une affirmation du brouillon ou une intuition. C'est le principal acquis méthodologique du domaine.

| Module | Affirmation | Réalité mesurée |
|---|---|---|
| M3 / M4b | un ajustement exponentiel suffit à dire si c'est chaotique | sur une croissance linéaire (le huit, stable) il donne un faux `λ` = 0,12 : décider sur l'amplification réelle (> 1000) |
| M4a | la précession se lit sur le vecteur de Runge-Lenz à des instants fixes | faux (il oscille dans l'orbite) : détecter les passages au périastre |
| M1 à M4 | mesurer l'ordre à `t = nT` | sur-convergence apparente (Verlet « d'ordre 4 ») : mesurer à 2,7 T et dans l'espace des phases |
| M5a | le frottement sec se met dans l'équation avec `sgn(v)` | ordre 1 pour tous les schémas, et un bloc qui devrait rester collé s'en va de 9 m en 5 s |
| M5b | détecter le contact seulement après le pas (modèle naïf) | la balle ne s'arrête jamais (`|vy|` résiduel proportionnel à `dt`) |
| M5c | l'ordre de résolution des impulsions compte dans le berceau | **faux** : écart exactement 0 sur 245 cas ; seule la causalité compte |
| M5c | RK4 est d'ordre 4 sur un contact de Hertz | **faux** : 2,3 à 2,75 (la force n'est pas lisse en `δ = 0`) |
| M6 | une toupie sous le spin critique tombe | **faux** : grande nutation, elle remonte (système sans perte fixé à un pivot) |
| M6 | garder `|q| = 1` suffit à bien intégrer une rotation | **faux** : Lie-Heun conserve la norme et dérive en énergie et en `L` |
| M7 | l'impulsion dérive de 10⁻⁷ par pas en float | **faux** : 5·10⁻⁹ à 3·10⁻⁶ cumulés |
| M7 | un GPU `double` séparerait l'arrondi du parallélisme | émulé par le pilote, des milliers de fois plus lent : inutilisable |
| M7 | le calcul GPU est correct puisque l'on voit les étoiles | `Renderer` : les points disparaissaient dès que d'autres tampons GL existaient ; calcul sain, rendu cassé |

## Un fil rouge : les quatre pièges numériques qui reviennent

1. **Mesurer au bon endroit** (instant, espace des phases, régime asymptotique).
2. **Les discontinuités** (événements, adhérence, `hMin` atteint par RK45) : toujours un budget de pas (`advance(..., maxSteps)`).
3. **Les invariants ne se valent pas** : l'énergie (Verlet borné), l'impulsion (forces opposées), le moment cinétique (structure symplectique), la norme d'un quaternion (renormalisation) se conservent pour des raisons différentes ; chaque schéma en conserve certains et pas d'autres.
4. **La précision finie** (M7) : une somme de N termes arrondis donne `u √N` ; centrer en `double` avant de convertir ; le chaos interdit de comparer des trajectoires longues, seulement des accélérations à état fixé et des grandeurs d'ensemble.

## Comment tout revérifier

```bash
cmake -S . -B build -G Ninja && cmake --build build           # aucun avertissement attendu (-Wall -Wextra -Wpedantic)
ctest --test-dir build --output-on-failure                     # « test_core : OK »
for s in 1 2 3 4 5 6 7 8 9 10 11; do ./build/physicslab --smoke-test --sim $s --level 6; done
./build/physicslab --gpu-test                                  # M7 : « gpu-test : OK » (pilote OpenGL 4.5 avec compute shaders)
```

Les temps du M7 fluctuent d'une exécution à l'autre sur une carte intégrée ; les précisions, non. Pour mesurer un temps CPU, construire en Release dans un dossier jetable (le dossier `build` de cette machine est en Debug, 5 à 10 fois plus lent).

## Ce qui reste (détail dans `docs/PASSATION.md`, section 7)

- `computeConvergence()` est dupliqué dans une dizaine de modules : à factoriser. `tests/test_core.cpp` (3311 lignes) est à scinder par domaine.
- Orbites non liées de Kepler ; triangle de Lagrange à masses inégales ; frottement quadratique ; chocs obliques et corps rigides avec contact ; toupie sur une table ; solution thêta du solide asymétrique.
- M7 : sommation compensée, relecture asynchrone, arbre de Barnes-Hut, disque galactique.
- Une fenêtre console s'ouvre à côté de l'application sous Windows ; identité Git et copyright de la licence à confirmer avec l'utilisateur.

## Suite

Fin de domaine : cours compilé pour non-initiés (`docs/cours/`), roue des domaines pour choisir le module suivant, prompt de reprise court, puis merge dans `main` et tag `mecanique-1` (à confirmer avec l'utilisateur).
