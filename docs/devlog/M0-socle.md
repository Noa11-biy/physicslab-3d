# M0 : le socle

Commits : `9c46f44` (initialisation), `9f6fa9b` (licence MIT), `16f667f` (dépendances), `8ce9ea3` (socle). Date : 2026-10-08.

## Objectif

Avant la moindre physique : un projet qui compile partout, une fenêtre OpenGL 4.5, une interface à panneaux, des briques de calcul réutilisables par tous les modules, et la notion de **niveau pédagogique** (6 niveaux, de la vulgarisation au doctorat) présente dès le départ plutôt que collée après coup.

## Ce qui a été construit

| Brique | Rôle |
|---|---|
| CMake (C++17, `-Wall -Wextra -Wpedantic` sur notre code seulement) | trois bibliothèques (`physicslab_core` sans aucune dépendance graphique, `physicslab_render`, l'application) et un exécutable de tests |
| GLFW 3.4, GLAD (GL 4.5 core, généré et versionné), Dear ImGui v1.92 (docking), ImPlot v1.0 | fenêtre, chargeur OpenGL, interface, courbes temps réel (sous-modules git, sauf GLAD) |
| `Vec3`, `Mat3`, `Quaternion`, `Constants` (SI, CODATA 2022) | maths de base ; le cœur reste en `double` |
| `Solver` | interface « premier ordre » `dy/dt = f(t, y)`, avec la convention `y = [positions \| vitesses]` ; toute la suite (Euler, Verlet, RK4, RK45) s'y branche |
| `World`, `Level` | horloge de simulation ; énumération des 6 niveaux et `atLeast(niveau, minimum)` |
| `Camera` (orbitale), `Renderer` (OpenGL 4.5 DSA) | rendu de lignes et de points colorés ; le rendu ne connaît rien de la physique |

Règle d'architecture posée ici et tenue jusqu'à M7 : **le cœur physique ne dépend pas d'OpenGL**. Les tests (`tests/test_core.cpp`, un seul exécutable avec des macros `CHECK` / `CHECK_NEAR` maison) ne lient que `physicslab_core` et tournent donc sans écran ni pilote graphique.

## Taille

25 fichiers, 1678 lignes ajoutées dans le commit du socle (les dépendances tierces sont dans un commit à part).

## Vérification

`testVec3`, `testMat3`, `testQuaternion`, `testEulerOrder` (tests unitaires des maths). Démarrage hors écran : `physicslab --smoke-test`.

## Pièges rencontrés

- GLFW 3.4 est un **tag**, pas une branche : `git clone --depth 1 --branch 3.4`, puis `git submodule add`.
- ImPlot v1.0 a changé d'API : `PlotLine(label, xs, ys, n, ImPlotSpec(...))`, le décalage d'un tampon circulaire passe par `ImPlotProp_Offset`.
- Chaîne de compilation réelle : MSYS2 **mingw64** (GCC 16.1), pas `ucrt64`.
- Fins de ligne : `.gitattributes` force LF dans le dépôt ; les avertissements « CRLF will be replaced by LF » sont normaux.

## Limites

Rien de physique à ce stade. L'interface du module (`SimulationModule`) n'existe pas encore : elle apparaît à M2, quand il y a deux simulations à héberger.
