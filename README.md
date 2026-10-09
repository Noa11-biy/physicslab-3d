# PhysicsLab 3D

Moteur de simulation physique 3D temps réel, du collège au doctorat.
Chaque phénomène s'explique à **6 niveaux** : vulgarisation, intéressé, collège, lycée, étudiant (L1–L3), doctorant/chercheur.

## Objectifs

- Un « Unity/Unreal de la physique » : une catégorie par domaine, chacune avec ses simulations.
- Visualisation 3D temps réel (champs, particules, ondes, fluides, courbure de l'espace-temps).
- Équations, constantes et conditions initiales modifiables en direct.
- Physique proche du réel, sans approximation « arcade ».
- Outils annexes : tableau périodique, constantes fondamentales, convertisseur d'unités, graphes temps réel.

## Stack

| Élément    | Choix                                                   |
|------------|---------------------------------------------------------|
| Langage    | C++17, CMake                                            |
| Plateformes| Windows (MSYS2 mingw64) et Linux                        |
| Rendu / GPU| OpenGL 4.5 + compute shaders (GLFW + GLAD)              |
| Interface  | Dear ImGui (docking) + ImPlot                           |

Chaque simulation lourde a deux backends : un **CPU de référence en `double`** (validation) et un **GPU en `float`** (gros volumes). Le CPU vérifie que le GPU reste dans une tolérance définie.

## Compilation

Prérequis : compilateur C++17 (GCC via MSYS2 sous Windows), CMake ≥ 3.20, Ninja, pilote graphique OpenGL 4.5.

```bash
git clone --recurse-submodules https://github.com/Noa11-biy/physicslab-3d.git
cd physicslab-3d
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure   # tests de validation
./build/physicslab                           # lance l'application
./build/physicslab --level 5 --sim 4        # niveau pédagogique 5 (1 à 6), simulation 4 (1 à 11 : M1, M2, M3, M3b, M4a, M4b, M5a, M5b, M5c, M6, M7)
./build/physicslab --gpu-test               # M7 : compare le calcul N-corps du GPU (compute shader) au CPU double, mesure les temps
```

Si le dépôt a été cloné sans `--recurse-submodules` : `git submodule update --init --recursive`.

## Structure du dépôt

```
include/physicslab/   en-têtes publics : core/ (maths, solveurs, monde), mechanics/ (problèmes de référence
                      avec solution exacte), render/ (OpenGL)
src/                  sources : core/, mechanics/, render/, app/ (fenêtre, interface, un module par simulation)
shaders/              *.comp, *.vert, *.frag
tests/                validations (invariants, solutions analytiques, écart CPU/GPU)
tools/                scripts utilitaires (capture d'écran automatique pour vérifier l'interface)
third_party/          GLFW 3.4, ImGui 1.92 (docking), ImPlot 1.0 en sous-modules ; GLAD généré (GL 4.5 core)
docs/
  PASSATION.md        état du projet, règles de travail, pièges et prochaine étape (à lire pour reprendre)
  cours/              cours compilés par module (PDF / Word)
  devlog/             journal de développement par module
```

## Feuille de route

Ordre de progression (cocher au fil des modules) :

- [ ] Mécanique
  - [x] M0 : socle (CMake, fenêtre, ImGui docking, maths, interface Solver, sélecteur de niveau)
  - [x] M1 : projectile avec frottement, comparaison Euler / Euler symplectique / Verlet / RK4 / RK45
  - [x] M2 : ressort-masse (libre, amorti, forcé, résonance)
  - [x] M3 : pendule simple (solution exacte elliptique) puis pendule double (chaos, exposant de Lyapunov)
  - [x] M4 : gravitation (Kepler à 2 corps avec précession numérique prédite ; N-corps : huit, triangle de Lagrange, amas chaotique)
  - [x] M5 : collisions et frottements (frottement sec de Coulomb ; chocs et rebonds ; berceau de Newton : contact de Hertz contre impulsions séquentielles)
  - [x] M6 : corps rigide (solides libres symétrique et asymétrique, toupie de Lagrange, intégrateurs d'orientation : Euler, RK4, groupe de Lie, découpage symplectique)
  - [x] M7 : N-corps sur GPU (compute shader, float ou double, kick-drift-kick sur GPU, comparaison au CPU double ; `--sim 11`, `--gpu-test`)
- [ ] Ondes (acoustique, optique)
- [ ] Thermodynamique
- [ ] Électrodynamique
- [ ] Fluides
- [ ] Plasma
- [ ] Atomique / Quantique / Nucléaire
- [ ] Relativité / Astro / Cosmologie
- [ ] Domaines appliqués (vivant, santé, géophysique, climat)

## Méthode par phénomène

1. Dérivation théorique à la main
2. Solution analytique d'un cas de référence
3. Discrétisation numérique (Euler / Verlet / RK4 ou schéma adapté)
4. Implémentation CPU `double`, puis GPU si le cas est lourd
5. Validation : invariants, comparaison à l'analytique, écart CPU/GPU
6. Déclinaison pédagogique sur les 6 niveaux

## Conventions Git

- Branche principale : `main`, toujours dans un état qui compile.
- Une branche par module : `module/<nom>` (ex. `module/mecanique`).
- Commits courts à l'impératif, préfixés : `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `build:`.
- Fin de module : merge dans `main` puis tag `v<module>-<n>` (ex. `mecanique-1`).

## Licence

Projet open source sous licence [MIT](LICENSE).
