// Point d'entrée de l'application graphique (fenêtre, interface, boucle principale).
#pragma once

namespace pl {

struct AppOptions {
    // Ouvre une fenêtre cachée, rend quelques images et quitte (0 = succès). Sert aux vérifications automatiques.
    bool smokeTest = false;
    // Ouvre une fenêtre cachée, compare le calcul GPU des N corps (M7) au CPU et mesure les temps (0 = succès).
    bool gpuTest = false;
    // Plus grand nombre de corps des mesures de temps de --gpu-test.
    int gpuMaxN = 16000;
    // Niveau pédagogique au démarrage (1 à 6).
    int level = 3;
    // Simulation affichée au démarrage (1 = M1 projectile, 2 = M2 ressort-masse, 3 = M3 pendule simple, 4 = M3b pendule double,
    // 5 = M4a orbite de Kepler, 6 = M4b problème à N corps, 7 = M5a frottement sec, 8 = M5b chocs et rebonds, 9 = M5c berceau de Newton, 10 = M6 corps rigide, 11 = M7 N corps sur GPU,
    // 12 = O0 une impulsion sur une grille).
    int simulation = 1;
};

int runApplication(const AppOptions& options);

}  // namespace pl
