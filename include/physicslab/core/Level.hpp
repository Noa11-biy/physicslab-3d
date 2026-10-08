// Les 6 niveaux pédagogiques, sélecteur global de l'interface.
#pragma once

namespace pl {

enum class Level {
    Vulgarisation = 0,  // images, analogies, aucune équation
    Interesse,          // concepts et formules en mots
    College,            // formules simples, unités, calculs numériques
    Lycee,              // vecteurs, dérivées, énergie, lois de Newton
    Etudiant,           // EDO, Lagrange/Hamilton, Maxwell, Schrödinger (L1-L3)
    Chercheur           // formulations complètes, schémas avancés, analyse d'erreur
};

inline constexpr int kLevelCount = 6;

constexpr const char* levelName(Level l) {
    switch (l) {
        case Level::Vulgarisation: return "1 - Vulgarisation";
        case Level::Interesse:     return "2 - Intéressé";
        case Level::College:       return "3 - Collège";
        case Level::Lycee:         return "4 - Lycée";
        case Level::Etudiant:      return "5 - Étudiant (L1-L3)";
        case Level::Chercheur:     return "6 - Doctorant / Chercheur";
    }
    return "?";
}

// Vrai si le niveau courant est au moins `minimum` (ex. atLeast(level, Level::Lycee)).
constexpr bool atLeast(Level current, Level minimum) {
    return static_cast<int>(current) >= static_cast<int>(minimum);
}

}  // namespace pl
