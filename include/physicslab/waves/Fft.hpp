// Transformée de Fourier rapide (Cooley-Tukey, base 2) : spectres d'un champ, modes d'une corde.
// Convention : X_k = somme_j x_j exp(-2 i pi j k / n) (sans facteur devant la directe, 1/n devant l'inverse).
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace pl::waves {

using Complex = std::complex<double>;

bool isPowerOfTwo(std::size_t n);

// Plus petite puissance de 2 supérieure ou égale à n (1 pour n = 0).
std::size_t nextPowerOfTwo(std::size_t n);

// Transformée directe en place ; a.size() doit être une puissance de 2 (assert en Debug).
void fft(std::vector<Complex>& a);

// Transformée inverse en place (divisée par n) : ifft(fft(x)) = x.
void ifft(std::vector<Complex>& a);

// Transformée directe naïve en O(n²) : référence indépendante pour les tests, tout n est permis.
std::vector<Complex> dft(const std::vector<Complex>& x);

// Spectre d'amplitude d'un signal réel : A_k = |X_k| 2 / n pour k = 1 .. n/2 - 1 (une sinusoïde d'amplitude A tombant
// exactement sur la fréquence k donne A_k = A), A_0 = |X_0| / n (la moyenne) et A_{n/2} = |X_{n/2}| / n.
// Le signal est complété par des zéros jusqu'à la prochaine puissance de 2. Renvoie n/2 + 1 valeurs, k = 0 .. n/2.
std::vector<double> amplitudeSpectrum(const std::vector<double>& x);

// Spectre d'amplitude avec fenêtre de Hann (réduit la fuite spectrale d'un enregistrement qui ne contient pas un nombre entier de
// périodes), complété par des zéros jusqu'à `padding` fois la longueur (puissance de 2). Renvoie n/2 + 1 valeurs ; `binWidth` reçoit
// la largeur d'une case en Hz (sampleRate / n_complété). L'amplitude d'une sinusoïde est atténuée d'environ 2 par la fenêtre :
// utile pour la FORME et les positions des raies, pas pour leur hauteur absolue.
std::vector<double> windowedSpectrum(const std::vector<double>& samples, double sampleRate, int padding, double* binWidth);

// Fréquence [Hz] du pic le plus haut d'un spectre déjà calculé (cases de largeur binWidth) dans [guess (1 - searchFraction),
// guess (1 + searchFraction)], par interpolation parabolique du logarithme de l'amplitude autour du maximum. 0 s'il n'y a pas de pic.
double peakFromSpectrum(const std::vector<double>& spectrum, double binWidth, double guess, double searchFraction = 0.25);

// Idem depuis les échantillons : fenêtre de Hann, remplissage x4, puis peakFromSpectrum.
double peakFrequency(const std::vector<double>& samples, double sampleRate, double guess, double searchFraction = 0.25);

}  // namespace pl::waves
