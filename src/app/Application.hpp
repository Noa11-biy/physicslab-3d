// Point d'entrée de l'application graphique (fenêtre, interface, boucle principale).
#pragma once

namespace pl {

struct AppOptions {
    // Ouvre une fenêtre cachée, rend quelques images et quitte (0 = succès). Sert aux vérifications automatiques.
    bool smokeTest = false;
};

int runApplication(const AppOptions& options);

}  // namespace pl
