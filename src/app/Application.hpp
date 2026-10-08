// Point d'entrée de l'application graphique (fenêtre, interface, boucle principale).
#pragma once

namespace pl {

struct AppOptions {
    // Ouvre une fenêtre cachée, rend quelques images et quitte (0 = succès). Sert aux vérifications automatiques.
    bool smokeTest = false;
    // Niveau pédagogique au démarrage (1 à 6).
    int level = 3;
    // Simulation affichée au démarrage (1 = M1 projectile, 2 = M2 ressort-masse, 3 = M3 pendule simple, 4 = M3b pendule double,
    // 5 = M4a orbite de Kepler, 6 = M4b problème à N corps,
    // 7 = M5a frottement sec).
    int simulation = 1;
};

int runApplication(const AppOptions& options);

}  // namespace pl
