// O2 : dispersion numérique et condition CFL du schéma saute-mouton (Wave1D, Wave2D) : formules exactes, confrontées aux simulations par les tests.
//
// Théorie. Une onde plane u = exp(i (k x - w t)) est solution EXACTE du schéma (la grille est linéaire et invariante par translation) à condition que
//     1D :  sin(w dt / 2) = C sin(k dx / 2)                                   C = c dt / dx  (nombre de Courant)
//     2D :  sin²(w dt / 2) = C² [sin²(kx dx / 2) + sin²(ky dx / 2)]
// Notations ci-dessous : x = k dx (« phase par case » : 2 pi / x cases par longueur d'onde).
//   vitesse de phase     v_phase / c = w / (c k) = 2 asin(C sin(x/2)) / (C x)        (1 - (1 - C²) x² / 24 + O(x^4) : la grille ralentit les ondes courtes)
//   vitesse de groupe    v_group / c = dw/dk / c = cos(x/2) / sqrt(1 - C² sin²(x/2)) (vitesse d'un paquet d'ondes : celle de l'énergie)
//   C -> 0 (espace seul, c'est la chaîne de masses d'O1) : v_phase / c = sin(x/2) / (x/2), v_group / c = cos(x/2)
//   C = 1 : v_phase = v_group = c exactement : l'erreur d'espace et l'erreur de temps (de signes opposés) se compensent (1D seulement).
// Stabilité (von Neumann) : pour chaque k, la récurrence en temps T^{m+1} = (2 - 4 C² s²) T^m - T^{m-1} (s² = sin²(k dx / 2) en 1D, somme en 2D) a pour
// racines lambda + 1/lambda = 2 - 4 C² s². Si |2 - 4 C² s²| <= 2 elles sont de module 1 (onde qui se propage) ; sinon la racine de plus grand module,
//     |lambda| = |b| + sqrt(b² - 1),  b = 1 - 2 C² s²,
// grandit à chaque pas (instabilité). Le pire mode est celui de plus grand s² : s² = 1 en 1D (C > 1 instable), s² = 2 en 2D (C > 1/sqrt(2)).
// La condition CFL (Courant-Friedrichs-Lewy) est la version géométrique : le domaine de dépendance numérique doit contenir le domaine physique.
#pragma once

namespace pl::waves::dispersion {

constexpr double kPi = 3.14159265358979323846;

// Cases par longueur d'onde pour une phase x = k dx par case.
inline double pointsPerWavelength(double x) { return 2.0 * kPi / x; }

// Vitesses rapportées à c, pour la phase x = k dx et le nombre de Courant cfl (0 : espace seul, la chaîne). NaN si C sin(x/2) > 1 (onde non propagée).
double phaseRatio1D(double x, double cfl);
double groupRatio1D(double x, double cfl);

// 2D : onde de nombre d'onde k = x / dx dont la direction fait l'angle theta avec l'axe x (kx = k cos theta, ky = k sin theta). v_phase / c.
// Elle dépend de theta : la grille est ANISOTROPE (en x ou y, comme 1D ; sur la diagonale, plus lente aux courtes longueurs d'onde).
double phaseRatio2D(double x, double theta, double cfl);

// Facteur de croissance par pas |lambda| >= 1 du mode (kx dx, ky dx) ; 1 si le mode est stable. 1D : ky dx = 0 et s² = sin²(kx dx / 2).
double growthPerStep1D(double x, double cfl);
double growthPerStep2D(double kxdx, double kydx, double cfl);

// Limites de stabilité : C <= 1 (1D), C <= 1/sqrt(2) (2D).
double cflLimit(int dimension);

}  // namespace pl::waves::dispersion
