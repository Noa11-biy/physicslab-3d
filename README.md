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
| Plateformes| Windows (MSYS2 UCRT64) et Linux                         |
| Rendu / GPU| OpenGL 4.5 + compute shaders (GLFW + GLAD)              |
| Interface  | Dear ImGui (docking) + ImPlot                           |

Chaque simulation lourde a deux backends : un **CPU de référence en `double`** (validation) et un **GPU en `float`** (gros volumes). Le CPU vérifie que le GPU reste dans une tolérance définie.

## Structure du dépôt

```
include/      en-têtes (.hpp), un sous-dossier par module
src/          sources (.cpp), un sous-dossier par module
shaders/      *.comp, *.vert, *.frag
tests/        validations (invariants, solutions analytiques, écart CPU/GPU)
third_party/  dépendances externes (GLFW, GLAD, ImGui, ImPlot)
docs/
  cours/      cours compilés par module (PDF / Word)
  devlog/     journal de développement par module
```

## Feuille de route

Ordre de progression (cocher au fil des modules) :

- [ ] Mécanique
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
