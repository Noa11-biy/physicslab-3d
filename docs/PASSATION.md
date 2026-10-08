# Passation : PhysicsLab 3D

Ce fichier permet de reprendre le projet dans une nouvelle conversation sans rien perdre. À lire en entier avant de coder.
Dernière mise à jour : fin de M4 (Mécanique : M0 à M4 terminés), branche `module/mecanique`.

## 1. Rôle et règles de travail

**Rôle** : professeur de physique et mentor C++/GPU. On construit ensemble, module par module, un moteur de simulation de
toute la physique, du collège au doctorat. L'utilisateur est **professeur** : le projet est aussi un outil pédagogique pour
ses élèves, et il faut pouvoir l'expliquer à des non-initiés. Parler **français**, ton simple et direct.

**Règles**
- Pendant un module : réponses courtes, pas de longues explications, pas d'images ni de SVG dans le chat.
- Code propre et commenté, sans sur-architecture du module en cours. Physique proche du réel (pas d'approximation « arcade »).
- Avant chaque nouveau domaine : une feuille de route courte (concepts, équations, difficultés, simulations).
- **Méthode obligatoire pour chaque phénomène** : (1) dérivation théorique à la main, (2) solution analytique d'un cas de référence,
  (3) discrétisation numérique, (4) implémentation CPU `double` puis GPU `float` si le cas est lourd, (5) validation (invariants,
  comparaison à l'analytique, écart CPU/GPU), (6) déclinaison pédagogique sur les 6 niveaux.
- **Fin de domaine** (ex. : quand toute la Mécanique est finie), livrer : (1) un devlog complet, (2) un prompt court de reprise
  prêt à copier, (3) le cours compilé en PDF ou Word pour expliquer à des « moldus » (l'utilisateur est prof et s'en sert avec
  ses élèves), (4) la « roue des domaines » pour choisir le module suivant.
- Ne jamais afficher une mesure sans l'avoir vérifiée : les chiffres du devlog doivent venir de l'exécution réelle.

**Les 6 niveaux** (sélecteur global, `enum class Level`) : 1 Vulgarisation (images, analogies, zéro équation, 1-2 curseurs),
2 Intéressé (concepts, formules en mots), 3 Collège (formules simples, unités), 4 Lycée (vecteurs, dérivées, énergie, Newton),
5 Étudiant L1-L3 (EDO, Lagrange/Hamilton), 6 Doctorant/Chercheur (formulations complètes, schémas avancés, analyse d'erreur).
Le niveau change : l'explication, les paramètres exposés, les équations, le détail des courbes.

**Stack** : C++17, CMake + Ninja, Windows (MSYS2) et Linux, OpenGL 4.5 + compute shaders (GLFW + GLAD), Dear ImGui (docking)
+ ImPlot. Deux backends visés pour les simulations lourdes : CPU `double` de référence, GPU `float` (le CPU vérifie le GPU).

**Ordre des domaines** : Mécanique → Ondes (acoustique, optique) → Thermodynamique → Électrodynamique → Fluides → Plasma →
Atomique/Quantique/Nucléaire → Relativité/Astro/Cosmologie → Domaines appliqués. Liste complète des domaines dans le brief
initial (Mathematical Physics, Chaos, Milieux continus, Acoustique, Optique, Magnétisme, Électronique, Solide, Cryogénie,
Moléculaire, Information quantique, Particules, Géophysique, Météo, Biophysique, Médical, etc.).

## 2. État d'avancement

| Module | Contenu | État |
|---|---|---|
| M0 | Socle : CMake, fenêtre, ImGui docking, maths, `Solver`, sélecteur de niveau | fait |
| M1 | Projectile avec frottement linéaire, comparaison de 5 solveurs, solution exacte | fait |
| M2 | Ressort-masse libre / amorti / forcé, résonance, 3 régimes d'amortissement | fait |
| M3 | Pendule simple (solution elliptique exacte) puis pendule double (chaos, Lyapunov) | fait |
| M4 | Gravitation : M4a Kepler à 2 corps (précession numérique prédite), M4b N-corps CPU (huit, Lagrange, amas) | fait |
| **M5** | **Collisions et frottements** | **à faire (suite logique)** |
| M6 | Corps rigide (quaternions, toupie, solide libre) | à faire |
| M7 | N-corps GPU en compute shader (écart CPU/GPU) | à faire |

Dépôt : https://github.com/Noa11-biy/physicslab-3d (public, licence MIT). Branche de travail : `module/mecanique`
(fin de domaine : merge dans `main` puis tag `mecanique-1`). Conventions de commit : `feat(mecanique): ...`, `fix:`, `docs:`,
`build:`, `test:` en français, avec le trailer `Co-Authored-By` demandé par l'environnement. Pousser sur la branche après chaque
module quand l'utilisateur le demande ou continue la chaîne (c'était son souhait jusqu'ici).

### Résultats mesurés (réels, à réutiliser dans le cours et le devlog)

- **Ordres de convergence** (pente log-log de l'erreur finale en fonction de dt) : Euler 1,0 ; Euler symplectique 1,0 ;
  Verlet des vitesses 2,0 ; RK4 4,0 (4,03 sur le ressort, 3,98 sur le pendule). RK45 adaptatif : erreur 9·10⁻¹¹ pour une tolérance 10⁻¹⁰.
- **Dérive d'énergie d'Euler sans frottement** (chute libre) : dE = ½ m g² t dt, vérifié à 10⁻⁶.
- **Oscillateur** : Euler multiplie l'énergie par (1 + ω0² dt²) à chaque pas (mesuré 7,5543, prédit 7,5543).
  Verlet et Euler symplectique : énergie bornée sans dérive sur 100 périodes ; RK4 : dissipation lente, ordre (ω0 dt)⁶.
- **Pendule simple** : T/T0 = 1,0732 à 60° (exact, par K(k)) ; solution elliptique validée contre un RK45 à 10⁻¹³ (écart ~10⁻¹²).
- **Pendule double** (m1=m2=l1=l2=1, θ1=120°, θ2=−10°) : λ ≈ 1,1 /s, amplification ×4400 en 8 s pour un écart initial 10⁻⁹ rad.
  Horizon de prédictibilité (écart > 0,1) à dt = 1/60 s : Euler 0,43 s, Euler symplectique 0,52 s, Verlet 3 s, RK4 8,6 s.
- **Kepler (M4a)**, e = 0,5, unités normalisées (GM = a = 1, T = 2π) : rapports d'erreur lors du doublement des pas : Euler 1,98
  (mesuré à 0,35 T : à 2,7 T il est encore à 1,76 avec 25600 pas), Euler symplectique 1,98 (25600 → 51200 pas), Verlet 4,00, RK4 16,95
  (3200 → 6400 pas). RK45 (tol 1e-10) : erreur 8,5e-9 en 654 pas sur 2,7 T.
  Équation de Kepler : résidu ≤ 4,4e-16, 4 à 9 itérations de Newton pour e ≤ 0,99 (20 à e = 0,9999 près du périastre).
  100 orbites à 200 pas/orbite : Euler explicite s'échappe (E > 0) ; Euler symplectique : dE/E borné à 5,2e-2 ; Verlet : borné à 2,7e-3 ;
  RK4 : 9,9e-6 puis 9,3e-5 (dissipation lente) ; L conservé à l'arrondi par les symplectiques (testé < 1e-9) mais pas par RK4 (1,4e-5).
  **Précession numérique** (rétrograde) : Verlet −0,2232 °/orbite, Euler symplectique −0,2234 (identiques à l'ordre dominant : Lie-Trotter est
  conjugué à Strang), RK4 +0,00029 (770 fois moins). En h² : −0,886 → −0,0559 de 100 à 400 pas/orbite (rapport 15,85).
  Formule du hamiltonien modifié, Δω = −(π/8)(GM h²/a³)(4+e²)/(1−e²)³ : écart à la mesure 0,02 à 0,2 % pour e de 0,1 à 0,9.
  e = 0,9 : à 200 pas/orbite RK4 perd 753 fois l'énergie ; RK45 prend 698 pas pour 2,3e-6 sur 2 orbites, avec un pas 300 fois plus court
  au périastre qu'à l'apoastre (RK4 fixe : 32000 pas pour 7,8e-6).
- **N corps (M4b)** : « huit » de Chenciner-Montgomery (données de Simó) : E = −1,28714199, P = L = 0, fermeture à 8e-8 à T = 6,32591398
  (meilleure période retrouvée à 1e-7 près). Une perturbation de 1e-3 y croît **linéairement** (0,058 à 2 T → 0,275 à 10 T) : stable.
  Amplification d'un écart de 1e-9 entre t = 2 et t = 16 : huit 7,3 ; amas de 6 corps (graine 42) 6,6e4 (λ ≈ 0,8) ; graine 7 : 5e7.
  Ordres sur le huit (0,35 T) : Euler 1,97, Euler symplectique 2,00, Verlet 4,00, RK4 16,1. Triangle de Lagrange RK4 : 8,0e-6 (100 pas/période)
  → 8,3e-11 (1600), ordre 4. Amas de 6 corps, dt = 1e-3, t = 20 : |P| ≤ 3e-15 pour TOUS les solveurs (forces opposées) ; |ΔL| : Euler 3,4e-3,
  symplectiques ≤ 4e-15 (exact à l'arrondi), RK4 9,8e-11 (invariant quadratique non conservé) ; max |dE/E| : Euler 0,475, sympl. 2,4e-2,
  Verlet 8,0e-4, RK4 1,8e-7. Deux corps : accord avec Kepler à 1,5e-11.

## 3. Architecture du code

```
include/physicslab/
  core/        Vec3, Mat3, Quaternion, Constants (SI, CODATA 2022), Level, Solver, World
  mechanics/   problèmes de référence avec solution exacte : Projectile, Oscillator, Pendulum, DoublePendulum, Kepler (+ PeriapsisTracker),
               NBody (nbody::accelerations : interface « force sur chaque particule » découplée de l'intégrateur, reprise par le GPU en M7)
  render/      Camera (orbitale, float), Renderer (OpenGL 4.5 DSA : lignes et points colorés)
src/core, src/mechanics, src/render   implémentations
src/app/       Application (fenêtre, thème, disposition, menu, boucle), SimulationModule (interface),
               UiCommon (SolverSet, StepClock, drawResultTable, Series, sliders), un module par simulation
shaders/       line.vert, line.frag
tests/         test_core.cpp (un seul exécutable, CHECK/CHECK_NEAR maison)
tools/         screenshot.ps1 (capture automatique pour vérifier l'interface)
third_party/   glfw 3.4, imgui v1.92.9b-docking, implot v1.0 (sous-modules) ; glad généré (GL 4.5 core) versionné
```

**Briques à réutiliser**
- `Solver` (premier ordre `dy/dt = f(t,y)`, convention `y = [positions | vitesses]`) : `ExplicitEuler`, `SymplecticEuler`,
  `VelocityVerlet`, `RK4`, `RK45` (Dormand-Prince adaptatif, `relTol`, `absTol`, statistiques). `step()` renvoie le pas réellement
  avancé ; `advance(solver, f, t, y, dt)` enchaîne les pas pour couvrir tout `dt` (nécessaire avec RK45).
- Chaque fichier de `mechanics/` expose un `XxxProblem` (paramètres, `rhs()`, énergie, **solution exacte ou référence RK45
  1e-13**) et une fonction `xxxError(problem, solver, steps, tEnd)` pour mesurer l'ordre de convergence.
- `SimulationModule` (src/app) : interface commune (`title`, `explanation(level)`, `reset`, `update`, `drawControls`,
  `drawInvariants`, `drawGraphs`, `hasAnalysis`/`drawAnalysis`, `drawScene`, `frameCamera`). L'application ouvre les fenêtres
  ImGui « Simulation », « Explication », « Invariants », « Graphes », « Analyse » et le module les remplit.
- `SolverSet` : les 5 solveurs, couleurs fixes (Euler orange, symplectique citron, Verlet violet, RK4 vert, RK45 rose, exacte bleu),
  visibilité par niveau (1-2 : seul RK4 « Ordinateur » ; 3-4 : Euler contre RK4 ; 5+ : cases à cocher), libellés par niveau.
- `StepClock` (pas fixe, accumulateur, dernier pas raccourci pour finir pile, échantillonnage ≤ 240 Hz), `drawResultTable`,
  `strf`, `lineSpec`, `Series` (tampon circulaire pour ImPlot).
- Renderer : `draw(Primitive::{Lines, LineStrip, Points}, vertices, pointSize)` ; la grille et les axes sont des helpers.

**Ajouter un module (liste de contrôle)**
1. `include/physicslab/<domaine>/Xxx.hpp` + `src/<domaine>/Xxx.cpp` : problème, solution exacte/référence, erreur ; les ajouter à
   `physicslab_core` dans `CMakeLists.txt`.
2. Tests dans `tests/test_core.cpp` **avant** l'interface : solution exacte (résidu de l'EDO par différences finies, conditions
   initiales), invariants, ordres de convergence, comparaison à une référence indépendante.
3. `src/app/XxxModule.{hpp,cpp}` réalisant `SimulationModule`, ajouté aux sources de `physicslab` et enregistré dans
   `Application.cpp` (`modules_.push_back`). Mettre à jour `--sim N` (`main.cpp`, `Application.hpp`) et le README.
4. Contenu des 6 niveaux : textes d'explication, paramètres exposés, colonnes du tableau d'invariants, graphes.
5. Vérifier (section 5), commiter, pousser.

## 4. Pièges rencontrés (à ne pas refaire)

**Numérique**
- Mesurer un ordre de convergence **à un instant qui n'est pas un multiple de la période** (on utilise 2,7 T) et **dans l'espace des
  phases** : à t = n T l'erreur de phase d'ordre 1 s'annule (sur-convergence apparente : Verlet « ordre 4 », symplectique « ordre 2 »).
- Un schéma peut être **pré-asymptotique** : RK4 sur un grand angle n'atteint un rapport d'erreur de 16 qu'à ≥ 1280 pas ; Euler
  n'atteint 2 qu'à ≥ 1600 pas. Ajuster le point de mesure, jamais le seuil de tolérance du test, sans avoir regardé les valeurs réelles.
- Sans frottement, la chute libre est un polynôme de degré 2 : Verlet et RK4 l'intègrent exactement, il n'y a rien à mesurer (l'étude
  de convergence impose alors un frottement).
- Euler symplectique et Verlet ne sont symplectiques que pour un hamiltonien **séparable** ; ce n'est pas le cas du pendule double.
- Un test qui passe du premier coup n'est pas forcément bon : afficher les valeurs réelles avant d'être satisfait.
- **Kepler** : Newton seul échoue en M = 0 (la racine est au bord du crochet) : cas particulier E = 0. Mesurer la précession par l'angle du vecteur de
  Runge-Lenz à des instants fixes est faux (il oscille au sein d'une orbite et la période numérique n'est pas T) : détecter les passages au périastre
  (`PeriapsisTracker`). À e = 0,9 et 200 pas/orbite tous les schémas à pas fixe sont détruits (le périastre n'est pas résolu) : c'est le sujet de RK45.
- **N corps** : le triangle de Lagrange à masses égales est **instable** (critère de Routh) ; un ajustement exponentiel (Lyapunov) sur une croissance
  seulement linéaire donne un faux λ (0,12 pour le huit, stable) : décider du chaos sur l'amplification réelle de l'écart (> 1000), pas sur λ.
  L'énergie des schémas d'ordre élevé ne reste pas « bornée » de façon lisse sur un amas (les rencontres rapprochées la font sauter) : comparer des maxima,
  ne pas exiger un plateau.
- **Interface** : la police Segoe UI n'a pas ∇ ni ∝ (affichés « � ») ; ∂, ᵀ, ω, √, ≈, Δ passent. Plus de 3 colonnes numériques dans « Invariants » (~400 px)
  écrasent la colonne « Méthode » : scinder en deux tableaux. `TextDisabled` ne passe pas à la ligne : utiliser `TextWrapped` colorée pour les notes.
- **PowerShell 5.1** : `Get-Content -Raw | Set-Content -Encoding utf8` ré-encode les accents (mojibake) et ajoute un BOM. Éditer avec l'outil Edit, ou avec
  `python -I` et `encoding="utf-8", newline=""`. En cas d'accident : `git checkout -- <fichier>`.

**Outils**
- Toolchain réelle : MSYS2 **mingw64** (GCC 16.1), pas `ucrt64` (non installé). Configure : `cmake -S . -B build -G Ninja`.
- ImPlot v1.0 : `PlotLine(label, xs, ys, n, ImPlotSpec(...))`, pas d'argument `offset` séparé (`ImPlotProp_Offset`). GLFW 3.4 est un tag, pas
  une branche : `git clone --depth 1 --branch 3.4 ...` puis `git submodule add`.
- Dans l'outil Bash, un **gros heredoc** (scripts Python, C++) peut échouer sur le quoting : écrire le script dans un fichier (outil
  Write, dans le dossier temporaire de la session) puis l'exécuter avec `python -I fichier.py`. Ne pas utiliser `/tmp`.
- Captures d'écran Windows : le processus PowerShell doit être déclaré DPI-aware (`SetProcessDPIAware`), sinon la capture est tronquée
  sur un écran mis à l'échelle. `tools/screenshot.ps1` le fait. Écran de l'utilisateur : 1920×1080 physique.
- ImGui crée `imgui.ini` dans le dossier courant (ignoré par git) ; s'il n'existe pas, la disposition par défaut est reconstruite.
- Avertissements « CRLF will be replaced by LF » : normaux, `.gitattributes` force LF dans le dépôt.
- Les chaînes de l'interface sont en UTF-8 avec lettres grecques (ω, ζ, θ, λ…) : elles s'affichent grâce à Segoe UI (Windows) ou DejaVu
  (Linux). Sans ces polices, la police intégrée d'ImGui n'affiche pas le grec.

## 5. Vérifier avant de commiter

```bash
cmake -S . -B build -G Ninja && cmake --build build           # aucun avertissement attendu (-Wall -Wextra -Wpedantic)
ctest --test-dir build --output-on-failure                     # "test_core : OK"
for s in 1 2 3 4 5 6; do ./build/physicslab --smoke-test --sim $s --level 6; done   # démarrage hors écran, GL 4.5
```
Build Debug propre depuis zéro de temps en temps (`-DCMAKE_BUILD_TYPE=Debug` dans un dossier jetable : vérifie les `assert`).
Vérification visuelle (PowerShell) : `.\tools\screenshot.ps1 -Level 5 -Sim 4 -Wait 8` puis ouvrir le PNG indiqué. Regarder au moins les
niveaux 1, 3, 5 et 6 d'un nouveau module, et un ancien module pour la non-régression. `-Clicks "x,y;x,y"` simule des clics.

Options de l'application : `--level 1..6`, `--sim 1..6` (1 = M1, 2 = M2, 3 = M3 pendule simple, 4 = M3b pendule double, 5 = M4a Kepler,
6 = M4b N corps), `--smoke-test`. Navigation 3D : clic gauche tourner, clic droit/milieu déplacer, molette zoomer.

## 6. Suite proposée : M5 collisions et frottements

Feuille de route à présenter d'abord (règle du domaine), puis (idées, à affiner) :
1. **Frottement sec de Coulomb** sur plan incliné : seuil d'adhérence tan θ = μs, glissement uniformément accéléré a = g (sin θ − μd cos θ),
   solution exacte. Difficulté numérique : le modèle est discontinu en v = 0 (broutement « stick-slip » des schémas naïfs) : prévoir un état
   « adhérence » explicite. Comparer aux solveurs de M1 (RK45 doit détecter l'événement).
2. **Collisions** : coefficient de restitution e (impulsion conservée, énergie perdue de 1 − e²), choc 1D élastique/inélastique avec solution
   exacte, choc oblique 2D (billard). Balle rebondissante : hauteurs en e^(2n) h0, **durée totale finie** (suite géométrique, paradoxe de Zénon),
   rebonds infinis en temps fini à traiter à part. Détection d'événement (instant exact du contact par bissection / RK45 « dense output »)
   plutôt que test à chaque pas fixe (sinon : effet tunnel, énergie fausse).
3. **Résolution par impulsions** de plusieurs corps (cône de friction de Coulomb), puis le berceau de Newton (cas piégeux : les chocs
   successifs et simultanés). Prépare M6 (corps rigide) et M7 (GPU).
4. Réutiliser : `Solver`/`RK45`, `SolverSet`, `StepClock`, `drawResultTable` (≤ 3 colonnes numériques par tableau), `PeriapsisTracker` n'est pas
   générique. `World` n'a toujours que pesanteur uniforme + frottement linéaire (la gravitation a son propre `NBodyProblem`) : y ajouter les
   contacts ou créer un `ContactWorld` séparé, à décider en feuille de route.

## 7. Dette technique et idées

- `computeConvergence()` est dupliqué dans 6 modules (deux variantes : erreur à `tEnd` fixe, pente ajustée là où l'erreur < 0,1 pour M4) :
  à factoriser (fonction générique prenant une fonction d'erreur). Idem `wrapped()` (texte enveloppé) dans KeplerModule et NBodyModule : à
  mettre dans UiCommon.
- Euler explicite et symplectique restent mal comparés à pas égal (ordre 1 chacun mais constantes très différentes) ; le pas par particule ou adaptatif
  n'existe pas pour les N corps (rencontres rapprochées : l'énergie saute). À reprendre avec M7 (GPU float).
- Orbites non liées (e ≥ 1, hyperboles) absentes de Kepler ; triangle de Lagrange à masses inégales (stable pour les Troyens) absent de NBody.
- `ProjectileModule` n'utilise pas encore `StepClock` ni `World` de la même façon que les autres modules.
- Absents : type `Tensor` (relativité), backend GPU / compute shaders, export de données (niveau 6), tableau périodique, constantes
  fondamentales et convertisseur d'unités en interface.
- Une fenêtre console s'ouvre à côté de l'application (à masquer en Release sous Windows).
- Dans le tableau d'invariants du niveau 6, certains libellés sont abrégés (« Euler sympl. »).
- Identité Git : `user.name` vaut `Noa-biy11` alors que le login GitHub est `Noa11-biy` (à corriger si l'utilisateur le souhaite).
- Copyright du fichier LICENSE : « Undermania » (nom du profil GitHub) ; à confirmer.

## 8. Prompt de reprise (à coller dans la nouvelle conversation)

> Reprends le projet PhysicsLab 3D dans ce dossier. Lis d'abord `docs/PASSATION.md` en entier (rôle, règles, état, architecture, pièges),
> puis `README.md`. Vérifie que ça compile et que les tests passent (section 5), puis présente-moi la feuille de route courte du
> module M5 (collisions et frottements) avant de coder. Réponses courtes pendant le module, en français.
