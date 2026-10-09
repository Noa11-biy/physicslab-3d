// Constantes physiques en unités SI.
// Valeurs exactes depuis la redéfinition du SI de 2019 ; les autres sont issues de CODATA 2022.
#pragma once

namespace pl::constants {

inline constexpr double pi = 3.14159265358979323846;

// --- Constantes exactes (définition du SI) ---
inline constexpr double c       = 299792458.0;        // vitesse de la lumière [m/s]
inline constexpr double h       = 6.62607015e-34;     // constante de Planck [J.s]
inline constexpr double hbar    = h / (2.0 * pi);     // constante de Planck réduite [J.s]
inline constexpr double e       = 1.602176634e-19;    // charge élémentaire [C]
inline constexpr double k_B     = 1.380649e-23;       // constante de Boltzmann [J/K]
inline constexpr double N_A     = 6.02214076e23;      // nombre d'Avogadro [1/mol]
inline constexpr double g0      = 9.80665;            // pesanteur normale [m/s^2]

// --- Constantes mesurées (CODATA 2022) ---
inline constexpr double G       = 6.67430e-11;        // constante gravitationnelle [m^3/(kg.s^2)]
inline constexpr double epsilon0 = 8.8541878188e-12;  // permittivité du vide [F/m]
inline constexpr double mu0     = 1.25663706127e-6;   // perméabilité du vide [N/A^2]
inline constexpr double m_e     = 9.1093837139e-31;   // masse de l'électron [kg]
inline constexpr double m_p     = 1.67262192595e-27;  // masse du proton [kg]

// --- Dérivées ---
inline constexpr double R_gas   = k_B * N_A;          // constante des gaz parfaits [J/(mol.K)]

}  // namespace pl::constants
