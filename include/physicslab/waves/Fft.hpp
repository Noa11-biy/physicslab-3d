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

}  // namespace pl::waves
