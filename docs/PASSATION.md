# Passation : PhysicsLab 3D

Ce fichier permet de reprendre le projet dans une nouvelle conversation sans rien perdre. À lire en entier avant de coder.
Dernière mise à jour : fin de M6 (Mécanique : M0 à M6 terminés : frottement sec, chocs et rebonds, berceau de Newton, corps rigide), branche `module/mecanique`. Prochaine étape : M7 N-corps GPU, dernier module de la Mécanique (brouillon de feuille de route en section 6, l'utilisateur ne l'a pas encore validé) ; ensuite les livrables de fin de domaine.

## 1. Rôle et règles de travail

**Rôle** : professeur de physique et mentor C++/GPU. On construit ensemble, module par module, un moteur de simulation de
toute la physique, du collège au doctorat. L'utilisateur est **professeur** : le projet est aussi un outil pédagogique pour
ses élèves, et il faut pouvoir l'expliquer à des non-initiés. Parler **français**, ton simple et direct.

**Règles**
- Pendant un module : réponses courtes, pas de longues explications, pas d'images ni de SVG dans le chat.
- Code propre et commenté, sans sur-architecture du module en cours. Physique proche du réel (pas d'approximation « arcade »).
- Avant chaque nouveau domaine : une feuille de route courte (concepts, équations, difficultés, simulations).
- **Méthode obligatoire pour chaque phénomène** : (1) dérivation théorique à la main, (2) solution analytique d'un cas de référence,
  (3) discrétisation numérique, (4) implémentation CPU `double` puis GPU `float` si le cas est lourd, (5) validation (invariants,
  comparaison à l'analytique, écart CPU/GPU), (6) déclinaison pédagogique sur les 6 niveaux.
- **Fin de domaine** (ex. : quand toute la Mécanique est finie), livrer : (1) un devlog complet, (2) un prompt court de reprise
  prêt à copier, (3) le cours compilé en PDF ou Word pour expliquer à des « moldus » (l'utilisateur est prof et s'en sert avec
  ses élèves), (4) la « roue des domaines » pour choisir le module suivant.
- Ne jamais afficher une mesure sans l'avoir vérifiée : les chiffres du devlog doivent venir de l'exécution réelle.

**Les 6 niveaux** (sélecteur global, `enum class Level`) : 1 Vulgarisation (images, analogies, zéro équation, 1-2 curseurs),
2 Intéressé (concepts, formules en mots), 3 Collège (formules simples, unités), 4 Lycée (vecteurs, dérivées, énergie, Newton),
5 Étudiant L1-L3 (EDO, Lagrange/Hamilton), 6 Doctorant/Chercheur (formulations complètes, schémas avancés, analyse d'erreur).
Le niveau change : l'explication, les paramètres exposés, les équations, le détail des courbes.

**Stack** : C++17, CMake + Ninja, Windows (MSYS2) et Linux, OpenGL 4.5 + compute shaders (GLFW + GLAD), Dear ImGui (docking)
+ ImPlot. Deux backends visés pour les simulations lourdes : CPU `double` de référence, GPU `float` (le CPU vérifie le GPU).

**Ordre des domaines** : Mécanique → Ondes (acoustique, optique) → Thermodynamique → Électrodynamique → Fluides → Plasma →
Atomique/Quantique/Nucléaire → Relativité/Astro/Cosmologie → Domaines appliqués. Liste complète des domaines dans le brief
initial (Mathematical Physics, Chaos, Milieux continus, Acoustique, Optique, Magnétisme, Électronique, Solide, Cryogénie,
Moléculaire, Information quantique, Particules, Géophysique, Météo, Biophysique, Médical, etc.).

## 2. État d'avancement

| Module | Contenu | État |
|---|---|---|
| M0 | Socle : CMake, fenêtre, ImGui docking, maths, `Solver`, sélecteur de niveau | fait |
| M1 | Projectile avec frottement linéaire, comparaison de 5 solveurs, solution exacte | fait |
| M2 | Ressort-masse libre / amorti / forcé, résonance, 3 régimes d'amortissement | fait |
| M3 | Pendule simple (solution elliptique exacte) puis pendule double (chaos, Lyapunov) | fait |
| M4 | Gravitation : M4a Kepler à 2 corps (précession numérique prédite), M4b N-corps CPU (huit, Lagrange, amas) | fait |
| M5a | Frottement sec de Coulomb sur plan incliné (événements, adhérence, 3 modèles numériques) | fait |
| M5b | Chocs et rebonds (restitution, balle rebondissante et accumulation de Zénon, choc de deux billes) | fait |
| M5c | Berceau de Newton (contact de Hertz contre impulsions séquentielles, `--sim 9`) | fait |
| M6 | Corps rigide (solide libre symétrique et asymétrique, toupie de Lagrange, 4 intégrateurs d'orientation, `--sim 10`) | fait |
| **M7** | **N-corps GPU en compute shader (écart CPU/GPU)** | **à faire (suite logique, dernier module de la Mécanique)** |

Dépôt : https://github.com/Noa11-biy/physicslab-3d (public, licence MIT). Branche de travail : `module/mecanique`
(fin de domaine : merge dans `main` puis tag `mecanique-1`). Conventions de commit : `feat(mecanique): ...`, `fix:`, `docs:`,
`build:`, `test:` en français, avec le trailer `Co-Authored-By` demandé par l'environnement. Pousser sur la branche après chaque
module quand l'utilisateur le demande ou continue la chaîne (c'était son souhait jusqu'ici).

### Résultats mesurés (réels, à réutiliser dans le cours et le devlog)

- **Ordres de convergence** (pente log-log de l'erreur finale en fonction de dt) : Euler 1,0 ; Euler symplectique 1,0 ;
  Verlet des vitesses 2,0 ; RK4 4,0 (4,03 sur le ressort, 3,98 sur le pendule). RK45 adaptatif : erreur 9·10⁻¹¹ pour une tolérance 10⁻¹⁰.
- **Dérive d'énergie d'Euler sans frottement** (chute libre) : dE = ½ m g² t dt, vérifié à 10⁻⁶.
- **Oscillateur** : Euler multiplie l'énergie par (1 + ω0² dt²) à chaque pas (mesuré 7,5543, prédit 7,5543).
  Verlet et Euler symplectique : énergie bornée sans dérive sur 100 périodes ; RK4 : dissipation lente, ordre (ω0 dt)⁶.
- **Pendule simple** : T/T0 = 1,0732 à 60° (exact, par K(k)) ; solution elliptique validée contre un RK45 à 10⁻¹³ (écart ~10⁻¹²).
- **Pendule double** (m1=m2=l1=l2=1, θ1=120°, θ2=−10°) : λ ≈ 1,1 /s, amplification ×4400 en 8 s pour un écart initial 10⁻⁹ rad.
  Horizon de prédictibilité (écart > 0,1) à dt = 1/60 s : Euler 0,43 s, Euler symplectique 0,52 s, Verlet 3 s, RK4 8,6 s.
- **Kepler (M4a)**, e = 0,5, unités normalisées (GM = a = 1, T = 2π) : rapports d'erreur lors du doublement des pas : Euler 1,98
  (mesuré à 0,35 T : à 2,7 T il est encore à 1,76 avec 25600 pas), Euler symplectique 1,98 (25600 → 51200 pas), Verlet 4,00, RK4 16,95
  (3200 → 6400 pas). RK45 (tol 1e-10) : erreur 8,5e-9 en 654 pas sur 2,7 T.
  Équation de Kepler : résidu ≤ 4,4e-16, 4 à 9 itérations de Newton pour e ≤ 0,99 (20 à e = 0,9999 près du périastre).
  100 orbites à 200 pas/orbite : Euler explicite s'échappe (E > 0) ; Euler symplectique : dE/E borné à 5,2e-2 ; Verlet : borné à 2,7e-3 ;
  RK4 : 9,9e-6 puis 9,3e-5 (dissipation lente) ; L conservé à l'arrondi par les symplectiques (testé < 1e-9) mais pas par RK4 (1,4e-5).
  **Précession numérique** (rétrograde) : Verlet −0,2232 °/orbite, Euler symplectique −0,2234 (identiques à l'ordre dominant : Lie-Trotter est
  conjugué à Strang), RK4 +0,00029 (770 fois moins). En h² : −0,886 → −0,0559 de 100 à 400 pas/orbite (rapport 15,85).
  Formule du hamiltonien modifié, Δω = −(π/8)(GM h²/a³)(4+e²)/(1−e²)³ : écart à la mesure 0,02 à 0,2 % pour e de 0,1 à 0,9.
  e = 0,9 : à 200 pas/orbite RK4 perd 753 fois l'énergie ; RK45 prend 698 pas pour 2,3e-6 sur 2 orbites, avec un pas 300 fois plus court
  au périastre qu'à l'apoastre (RK4 fixe : 32000 pas pour 7,8e-6).
- **N corps (M4b)** : « huit » de Chenciner-Montgomery (données de Simó) : E = −1,28714199, P = L = 0, fermeture à 8e-8 à T = 6,32591398
  (meilleure période retrouvée à 1e-7 près). Une perturbation de 1e-3 y croît **linéairement** (0,058 à 2 T → 0,275 à 10 T) : stable.
  Amplification d'un écart de 1e-9 entre t = 2 et t = 16 : huit 7,3 ; amas de 6 corps (graine 42) 6,6e4 (λ ≈ 0,8) ; graine 7 : 5e7.
  Ordres sur le huit (0,35 T) : Euler 1,97, Euler symplectique 2,00, Verlet 4,00, RK4 16,1. Triangle de Lagrange RK4 : 8,0e-6 (100 pas/période)
  → 8,3e-11 (1600), ordre 4. Amas de 6 corps, dt = 1e-3, t = 20 : |P| ≤ 3e-15 pour TOUS les solveurs (forces opposées) ; |ΔL| : Euler 3,4e-3,
  symplectiques ≤ 4e-15 (exact à l'arrondi), RK4 9,8e-11 (invariant quadratique non conservé) ; max |dE/E| : Euler 0,475, sympl. 2,4e-2,
  Verlet 8,0e-4, RK4 1,8e-7. Deux corps : accord avec Kepler à 1,5e-11.
- **Frottement sec (M5a)**, résistance k = 0,7 pour que la solution ne soit pas un polynôme : avec « événement + adhérence » on retrouve les ordres des
  solveurs (rapports d'erreur au doublement des pas : Euler 2,0, Euler symplectique 2,0, Verlet 4,0, RK4 16,2) ; le modèle **naïf** (sgn dans l'EDO) fait
  tomber Verlet à 2,0 et RK4 à 2,05 puis 1,94 (16,1 sur le même problème sans frottement sec) : l'ordre 1 pour tous. À n = 400 pas : erreur 6,7e-2
  (naïf) contre 3,3e-10 (événement), soit 2e8 fois mieux. Erreurs finales de l'événement à dt = 0,05 : RK4 4e-8 à 9e-8, Verlet 1,4e-3, RK45 1e-10 ; sans
  résistance Verlet et RK4 sont exacts à 2e-14 (phases polynomiales). Instant d'arrêt trouvé par bissection : 1e-6 (RK4), 1e-8 (RK45).
  **Piège de l'adhérence** (μd < tan θ ≤ μs, bloc au repos) : exact = reste collé ; naïf RK4 = parti à −9,2 m en 5 s (le modèle ignore μs : erreur de modèle,
  l'erreur ne tend pas vers 0 avec dt). Bloc qui doit rester collé après un arrêt (k = 0) : le naïf garde une vitesse résiduelle 1,18 |A| dt (∝ dt, jamais nulle).
  Régularisation tanh(v/ε) : erreur de position à 5 s 0,017 / 0,085 / 0,35 m pour ε = 0,01 / 0,05 / 0,2, **indépendante de dt** (mesuré dt = 0,01 et 0,001) ;
  le bloc qui doit rester collé rampe de −9,4 m en 5 s. **RK45 naïf sur un bloc qui doit rester collé : ne termine pas** (l'erreur locale reste d'ordre h sur
  la surface v = 0, le pas descend jusqu'à hMin = 1e-14) ; `advance(..., maxSteps)` le rend borné (400 pas par pas de calcul dans l'interface, « bloqué »).
- **Chocs et rebonds (M5b)** : balle e = 0,8 lâchée de 2 m : instants d'impact de la solution exacte égaux à la formule à 1,3e-15, sommets e^(2n) h0 à 7e-16 ;
  repos à 5,747939132 s contre t0 (1+e)/(1−e) = 5,747939222 s (écart −9e-8 dû au seuil de 1e-7) ; 81 rebonds (50 avec le seuil par défaut 1e-4), calcul exact en
  0,2 ms. Modèle **événement** : sans résistance RK4 et Verlet sont exacts à 4e-13 jusqu'après le repos ; avec k = 0,7 (erreur de position à t = 2,37 s) rapports
  Euler 2,0, Euler symplectique 2,0, Verlet 4,0, RK4 16,1 ; instant de repos retrouvé à 2e-6 (dt = 0,05), 8e-8 (0,02), 2e-9 (0,007), 1,5e-12 (0,001) : ordre 4,
  46 rebonds avec résistance. Modèle **naïf** (contact vu après le pas) : pente log-log 0,94 pour RK4 contre 4,0 pour l'événement ; à 400 pas 1,5e-2 contre 7,3e-11 ;
  premier impact vu en retard de moins d'un pas (0,6390 s pour 0,63868 s à dt = 1 ms) ; la balle **ne s'arrête jamais** : 0,5 s après le repos exact, |vy| max =
  0,78 / 0,28 / 0,11 m/s pour dt = 0,02 / 0,01 / 0,005 (2,3 à 4 g dt). Sans résistance Verlet et RK4 naïfs donnent exactement la même erreur (vols exacts) : toute
  l'erreur vient de la détection. Deux billes (disques de rayon 0,5, e = 1, masses égales, décalage 0,6 m) : contact à t = 1,0667 s, angle de sortie 90° ; événement :
  erreur ≤ 8e-14 pour les quatre solveurs, instant du choc exact à 1e-12 ; naïf : **même erreur pour les quatre solveurs** (RK4 compris), proportionnelle au pas
  (0,43 puis 0,10 puis 0,025 quand le pas est divisé par 4 à chaque fois), choc vu 3,3e-3 s trop tard à dt = 0,01 ; l'impulsion reste conservée même en naïf.
- **Berceau de Newton (M5c)**, unités normalisées (m = 1, R = 1/2, k = 1e4, v = 1) : constante du temps de contact (4/5) Γ(2/5) Γ(1/2) / Γ(9/10) = 2,943275 ;
  deux billes (RK45 serré) : δ_max = 0,0208138302 et T = 0,0612608299 retrouvés à 6e-13 et 1,6e-13 (relatif), sortie = échange des vitesses à 5e-13, P conservée à 1e-16.
  **Le résultat d'une chaîne ne dépend ni de k ni de v** (une seule échelle de longueur (m v²/k)^(2/5)) : 3 billes (k, v) = (1e4, 1), (1e6, 1), (1e4, 3), (2e3, 0,2) donnent les
  mêmes six chiffres ; seule la durée change (0,0943 s contre 0,01495 s pour k × 100, rapport 6,3 = 100^(2/5)).
  **Hertz, 3 billes, 1 lancée** : (−0,070952 ; +0,076403 ; +0,994549), soit x = −v1/v = 0,0710 sur la courbe des solutions [0 ; 1/3] (impulsions e = 1 : x = 0).
  Chaîne de N billes, 1 lancée : N = 4 (−0,0711 ; −0,0296 ; +0,1097 ; +0,9910) ; N = 5 (−0,071085 ; −0,030274 ; −0,014464 ; +0,127045 ; +0,988777) ; la dernière bille part à
  0,9945 / 0,9910 / 0,9888 / 0,9875 / 0,9867 pour N = 3 à 7 ; les billes de tête reculent d'environ 0,07 v (pas « une entre, une sort » exact). 2 lancées sur 5 :
  (−0,112615 ; −0,041960 ; +0,214486 ; +0,800367 ; +1,139722). Symétrie exacte (miroir + changement de repère galiléen) entre n lancées et N−n lancées : écart ≤ 1e-8 (N = 3 à 6).
  RK4 au pas T/100 : accord avec RK45 serré à 1e-5, énergie finale à 4e-7. Durée de la collision : 0,094 s (N = 3), 0,156 s (N = 5), 0,216 s (N = 7).
  **Impulsions séquentielles** (e = 1, masses égales, chaque choc échange les vitesses : résultat = tri des vitesses) : n entrent, n sortent, quel que soit l'ordre. **L'ordre de
  résolution ne change rien dans le berceau, même avec e < 1** : écart exactement 0 sur 245 cas (N = 3 à 9, n, e de 0,99 à 0) ; la causalité impose la séquence des chocs. Il
  compte pour une bille prise entre deux voisines qui s'approchent déjà : vitesses (1 ; 0,5 ; 0 ; 0), e = 0,5 : écart 0,070 entre gauche → droite et droite → gauche
  (nul pour e = 1). Avec e < 1, 3 billes, 1 lancée : e = 0,9 donne (0,0476 ; 0,0499 ; 0,9025) et perd 0,090 ; e = 0,5 donne (0,203 ; 0,234 ; 0,563) ; aucune bille ne recule
  jamais (différence qualitative avec Hertz).
  **Ordres des schémas sur le contact** (3 billes, contact décalé de 0,0123 pour ne pas tomber sur la grille, erreur dans l'espace des phases à t = 0,2) : Euler 1 (rapport 4,07 de
  800 à 3200 pas), Verlet 2 (15,4 de 400 à 1600 pas) ; **Euler symplectique = Verlet à 5e-13** (conjugués, et la mesure est faite hors contact : pas de force au début ni à la fin) ;
  **RK4 n'atteint pas l'ordre 4** : ordre apparent 2,29 / 2,38 / 2,51 / 2,75 selon la phase du contact sur la grille (décalages 0,0123 / 0 / 0,0231 / 0,0377 ; pas de 200 à 3200
  pas), rapports par doublement erratiques (73, 1, 8, 6, 1,4, 8) ; cause probable (non testée avec une force lisse) : la force k δ^(3/2) n'a pas de dérivée seconde en δ = 0.
  RK45 par défaut (tol 1e-8) : erreur 2,4e-7 (3 billes) à 3,4e-7 (5 billes), 88 pas acceptés, 26 refusés, 798 évaluations sur [0 ; 0,2]. Stabilité (3 lancées sur 7) : Euler
  explicite **diverge** à 10 pas par contact (valeurs infinies pour k = 1e4 et 1e6 ; E/E0 = 1672 pour k = 1e2), E/E0 = 1,881 à 30 pas et 1,1135 à 100 pas, identique pour k de 1e2 à 1e6 ; Euler
  symplectique et Verlet : 1,0004 à 10 pas, 1,0000 à 30 pas. **Amortissement de Hunt-Crossley** F = k δ^(3/2) (1 + (3/2) α δ'), v = 1 : e = 0,990099 / 0,952370 / 0,909016 / 0,768 pour α = 0,01 / 0,05 / 0,1 / 0,3, soit
  e ≈ 1/(1 + α v) (écart à 1 − α v : 0,99 (α v)² pour α = 0,01, 0,90 (α v)² pour α = 0,1), perte relative d'énergie 2 α v au premier ordre (0,0197 pour α = 0,01) ; chaîne de 5, α = 0,05 :
  (−0,0553 ; −0,0183 ; −0,0026 ; +0,1488 ; +0,9273), moins de recul qu'en contact conservatif.
- **Corps rigide (M6)**, convention q : repère du corps → repère fixe, ω dans le repère du corps. **Tenseurs d'inertie** : formules fermées contre une intégration directe sur grille, à 1e-3 (boîte) et
  5e-3 (sphère, cylindre : erreur de bord en 1/n). **Solide symétrique libre** (I1 = I2) : q(t) = q_L(L t / I1) q0 q_3(−Ω t), Ω = (I3 − I1) ω3 / I1 ; accord avec RK45 serré (relTol 1e-13) à 1e-14 ... 2e-12 sur ω et
  7e-14 ... 4e-12 sur l'orientation (t = 3 et 25) ; la même formule avec le signe de la rotation propre faux donne 0,56 ... 2,8 (O(1)) : le test est sensible (sauf pour la sphère, Ω = 0). Invariants le long de la
  trajectoire de référence (I = (1, 2, 3)) : E à 2e-13 et L (vecteur fixe) à 1e-11 jusqu'à t = 100, |q| − 1 à 9e-13. **Fonctions de Jacobi** (séries de Fourier en nome q, DLMF 22.11) contre une intégration RK4 de
  sn' = cn dn, cn' = −sn dn, dn' = −m sn cn : 1e-10 pour m de 1e-14 à 0,999 ; identités sn² + cn² = 1 à 4e-15 (1,4e-14 pour m = 1e-14), jusqu'à m = 1 − 1e-10 à 4e-15. **Solide asymétrique** (I = (1, 2, 3), départ (a, 0, c),
  formules de Landau-Lifchitz § 37) : m = 0,75 / 0,1875 / 0,853333 pour (a, c) = (1,5 ; 1) / (2 ; 0,5) / (1,5 ; 0,8), périodes 4K/λ = 8,62606 / 5,72746 / 11,08236 s ; |ω_exact − ω_RK45| ≤ 5,7e-13 pour t ≤ 40,
  |ω(T) − ω(0)| ≤ 1e-13. **Axe intermédiaire** : λ prévu 0,577350 (ω2 = 1), mesuré 0,577411 (pente de ln|ω1| de t = 8 à 14, écart relatif 1e-4) ; axes extrêmes stables (perturbation de 1,4e-3 → au plus 2,0e-3 sur 100 s).
  **Intégrateurs d'orientation** (corps I = (1, 2, 3), ω0 = (0,9 ; 0,5 ; 1,1), t = 6,3) : Euler ordre 1 (×2,01 par doublement à 1600 pas), RK4 + renormalisation ordre 4 (×16,06 ; ×263,5 de 100 à 400 pas ; erreur 3,5e-11 à
  1600 pas), Lie-Heun ordre 2 (×3,99), découpage symplectique ordre 2 (×4,00) et environ 3 fois plus précis que Heun à pas égal (3,7e-2 contre 1,07e-1 à 25 pas) ; avec le couple de la toupie (t = 2) : RK4 ×16,04, Heun ×4,00,
  découpage ×4,00 (Euler : non asymptotique, NaN à 50 pas). **Dérives sur 1000 s à h = 0,05** (20000 pas, ω0 = (0,1 ; 2 ; 0,1), proche de l'axe instable) : découpage max |ΔL|/L = 6,2e-14 (3,0e-14 à h = 0,2 : indépendant du pas),
  |q| − 1 = 3e-16, max |ΔE/E| = 6,9e-5 (h = 0,05), 2,74e-4 (0,1), 1,07e-3 (0,2) : en h² (×3,98 puis ×3,9) et le MÊME maximum dans chaque tiers du calcul (aucune dérive) ; RK4 + renormalisation : ΔE/E final −3,1e-6, ΔL/L 1,8e-6,
  |q| − 1 = 3e-16 ; RK4 sans renormalisation : |q| − 1 = 4e-6, ΔL/L 1,7e-5 ; **Lie-Heun** : |q| = 1 à l'arrondi MAIS ΔE/E +1,4e-2 et ΔL/L 7e-3 (dérive séculaire), **diverge à h = 0,2** (pas 2691) ; Euler diverge au pas 2251 (h = 0,05) et 148
  (h = 0,2). Euler explicite sur une rotation pure : |q|N = (1 + h² ω²/4)^(N/2) exactement (testé à 1e-12). **Toupie de Lagrange** (m = 1, I = (1,2 ; 1,2 ; 0,4) autour du pivot, l = 0,5) : E, L_z, L_3 conservés à 1e-10 ; précession
  régulière exacte q(t) = q_z(φ' t) q_x(θ) q_z(ψ' t) contre RK45 : accord à 1e-8 (borne du test) ; exemple θ = 0,6, ω3 = 25 : φ' = 0,517 rad/s (un tour en 12,2 s). **Nutation** par la cubique u'² = f(u) : θ0 = 0,9, ω0 = (0, 2, 20) :
  u1 = 0,621609968 (= cos 0,9), u2 = 0,926152996, u3 = 5,001797, T = 4K/√(β (u3 − u1)) = 1,069185 s ; u(T/2) − u2 = −2,3e-14, u(T) − u0 = 1,6e-14, u(3T) − u0 = 5e-14 (cos θ de l'interface touche u1 et u2) ; θ0 = 0,6, ω3 = 25, ω0 sans
  précession : u1 = 0,7791, u2 = 0,8253, T = 0,8357 s. **Toupie endormie** : spin critique 2 √(I1 m g l) / I3 = 12,12847 rad/s ; au-dessus (1,15 / 1,5 / 3 fois) une inclinaison de 0,01 reste ≤ 0,0203 / 0,0134 / 0,0106 ; en dessous (0,7 fois) la
  croissance est e^(γ t), γ prévu 1,44358, mesuré 1,44357 ; **la toupie ne tombe pas** : de 1e-6 elle monte jusqu'à 1,2457 rad (71°, échantillonné à t = 11 s) puis REMONTE (θ(20 s) ≈ 1e-5) : grande nutation d'un système sans perte
  fixé au pivot. **Intégrateurs avec couple** (nutation, 100 s à h = 0,005) : découpage ΔL_z/L_z = 2,8e-13 (le couple n'a pas de composante verticale : les « coups » ne la changent pas), ΔE/E 3,9e-5 borné, ΔL_3/L_3 2,8e-5 (L_3 n'est PAS exact) ;
  RK4 ΔL_z/L_z 5,1e-6, ΔE/E 7,5e-7 (à ce pas RK4 conserve mieux E sur 100 s) ; Heun ΔL_z/L_z 7,5e-3, ΔE/E 2,3e-3 ; Euler diverge. Précession régulière (θ = 0,8, ω3 = 25) sur 10 s en 2000 pas : distance à l'exact RK4 5,0e-5, découpage 2,0e-4,
  Heun 8,2e-3 ; cos θ : écart 1,3e-6 / 8e-6 / 8,8e-4.

## 3. Architecture du code

```
include/physicslab/
  core/        Vec3, Mat3, Quaternion, Constants (SI, CODATA 2022), Level, Solver, World
  core/        ... + Events (advanceToEvent : instant exact où g(t, y) change de signe, bissection sur le pas ; brique pour tous les chocs)
  mechanics/   problèmes de référence avec solution exacte : Projectile, Oscillator, Pendulum, DoublePendulum, Kepler (+ PeriapsisTracker),
               Friction (InclineProblem : solution exacte par morceaux ; InclineRun : événement + adhérence ; modèles Naive / Regularized),
               Collision (collide1D, collideSpheres : restitution de Newton ; TwoBallProblem / TwoBallRun), Bounce (BounceProblem : solution exacte par vols
               successifs réutilisant ProjectileProblem ; BounceRun : modèles Naive / EventDriven, ContactModel),
               Cradle (hertz:: force, énergie, références exactes de deux billes ; CradleProblem : chaîne de N billes de Hertz avec amortissement de Hunt-Crossley et
               écart initial, run() jusqu'à la fin de la collision, reference() et cradleError() ; cradle:: sequentialImpulses, threeBallFamily),
               RigidBody (inertia:: sphere/box/cylinder ; RotationState = (q, ω corps) et état plat [q | ω] de 7 nombres ; rotationRhs avec couple optionnel ; FreeBodyProblem : exactSymmetric,
               exactAsymmetricOmega par jacobiSnCnDn, asymmetricPeriod, intermediateAxisGrowthRate, reference ; HeavyTopProblem : torque, steadyPrecession / exactSteady, nutation() (racines de la cubique) et
               nutationPeriod, sleepingCriticalSpin ; RotationIntegrator : EulerRotation, RK4Rotation(renormalize), LieHeunRotation, SplittingRotation ; integrateRotation, rotationDistance),
               NBody (nbody::accelerations : interface « force sur chaque particule » découplée de l'intégrateur, reprise par le GPU en M7)
  render/      Camera (orbitale, float), Renderer (OpenGL 4.5 DSA : lignes et points colorés)
src/core, src/mechanics, src/render   implémentations
src/app/       Application (fenêtre, thème, disposition, menu, boucle), SimulationModule (interface),
               UiCommon (SolverSet, StepClock, drawResultTable, Series, sliders), un module par simulation
shaders/       line.vert, line.frag
tests/         test_core.cpp (un seul exécutable, CHECK/CHECK_NEAR maison)
tools/         screenshot.ps1 (capture automatique pour vérifier l'interface)
third_party/   glfw 3.4, imgui v1.92.9b-docking, implot v1.0 (sous-modules) ; glad généré (GL 4.5 core) versionné
```

**Briques à réutiliser**
- `Solver` (premier ordre `dy/dt = f(t,y)`, convention `y = [positions | vitesses]`) : `ExplicitEuler`, `SymplecticEuler`,
  `VelocityVerlet`, `RK4`, `RK45` (Dormand-Prince adaptatif, `relTol`, `absTol`, statistiques). `step()` renvoie le pas réellement
  avancé ; `advance(solver, f, t, y, dt)` enchaîne les pas pour couvrir tout `dt` (nécessaire avec RK45).
- Chaque fichier de `mechanics/` expose un `XxxProblem` (paramètres, `rhs()`, énergie, **solution exacte ou référence RK45
  1e-13**) et une fonction `xxxError(problem, solver, steps, tEnd)` pour mesurer l'ordre de convergence.
- `SimulationModule` (src/app) : interface commune (`title`, `explanation(level)`, `reset`, `update`, `drawControls`,
  `drawInvariants`, `drawGraphs`, `hasAnalysis`/`drawAnalysis`, `drawScene`, `frameCamera`). L'application ouvre les fenêtres
  ImGui « Simulation », « Explication », « Invariants », « Graphes », « Analyse » et le module les remplit.
- `SolverSet` : les 5 solveurs, couleurs fixes (Euler orange, symplectique citron, Verlet violet, RK4 vert, RK45 rose, exacte bleu),
  visibilité par niveau (1-2 : seul RK4 « Ordinateur » ; 3-4 : Euler contre RK4 ; 5+ : cases à cocher), libellés par niveau.
- `StepClock` (pas fixe, accumulateur, dernier pas raccourci pour finir pile, échantillonnage ≤ 240 Hz), `drawResultTable`,
  `strf`, `lineSpec`, `Series` (tampon circulaire pour ImPlot).
- Renderer : `draw(Primitive::{Lines, LineStrip, Points}, vertices, pointSize)` ; la grille et les axes sont des helpers.

**Ajouter un module (liste de contrôle)**
1. `include/physicslab/<domaine>/Xxx.hpp` + `src/<domaine>/Xxx.cpp` : problème, solution exacte/référence, erreur ; les ajouter à
   `physicslab_core` dans `CMakeLists.txt`.
2. Tests dans `tests/test_core.cpp` **avant** l'interface : solution exacte (résidu de l'EDO par différences finies, conditions
   initiales), invariants, ordres de convergence, comparaison à une référence indépendante.
3. `src/app/XxxModule.{hpp,cpp}` réalisant `SimulationModule`, ajouté aux sources de `physicslab` et enregistré dans
   `Application.cpp` (`modules_.push_back`). Mettre à jour `--sim N` (`main.cpp`, `Application.hpp`) et le README.
4. Contenu des 6 niveaux : textes d'explication, paramètres exposés, colonnes du tableau d'invariants, graphes.
5. Vérifier (section 5), commiter, pousser.

## 4. Pièges rencontrés (à ne pas refaire)

**Numérique**
- Mesurer un ordre de convergence **à un instant qui n'est pas un multiple de la période** (on utilise 2,7 T) et **dans l'espace des
  phases** : à t = n T l'erreur de phase d'ordre 1 s'annule (sur-convergence apparente : Verlet « ordre 4 », symplectique « ordre 2 »).
- Un schéma peut être **pré-asymptotique** : RK4 sur un grand angle n'atteint un rapport d'erreur de 16 qu'à ≥ 1280 pas ; Euler
  n'atteint 2 qu'à ≥ 1600 pas. Ajuster le point de mesure, jamais le seuil de tolérance du test, sans avoir regardé les valeurs réelles.
- Sans frottement, la chute libre est un polynôme de degré 2 : Verlet et RK4 l'intègrent exactement, il n'y a rien à mesurer (l'étude
  de convergence impose alors un frottement).
- Euler symplectique et Verlet ne sont symplectiques que pour un hamiltonien **séparable** ; ce n'est pas le cas du pendule double.
- Un test qui passe du premier coup n'est pas forcément bon : afficher les valeurs réelles avant d'être satisfait.
- **Kepler** : Newton seul échoue en M = 0 (la racine est au bord du crochet) : cas particulier E = 0. Mesurer la précession par l'angle du vecteur de
  Runge-Lenz à des instants fixes est faux (il oscille au sein d'une orbite et la période numérique n'est pas T) : détecter les passages au périastre
  (`PeriapsisTracker`). À e = 0,9 et 200 pas/orbite tous les schémas à pas fixe sont détruits (le périastre n'est pas résolu) : c'est le sujet de RK45.
- **N corps** : le triangle de Lagrange à masses égales est **instable** (critère de Routh) ; un ajustement exponentiel (Lyapunov) sur une croissance
  seulement linéaire donne un faux λ (0,12 pour le huit, stable) : décider du chaos sur l'amplification réelle de l'écart (> 1000), pas sur λ.
  L'énergie des schémas d'ordre élevé ne reste pas « bornée » de façon lisse sur un amas (les rencontres rapprochées la font sauter) : comparer des maxima,
  ne pas exiger un plateau.
- **Frottement sec / discontinuités** : un solveur adaptatif (RK45) peut ne jamais terminer : il accepte des pas jusqu'à 1e-14 quand l'erreur reste trop grande ;
  toujours appeler `advance(..., maxSteps)` quand la dynamique peut être discontinue, et ne pas tester un RK45 sur le modèle naïf sans budget (un test peut se
  figer). Un plateau d'erreur indépendant du pas n'est pas forcément un défaut du schéma : ici l'arrêt tombait au même décalage (0,18 ms) de la frontière
  de pas pour tous les dt multiples de 0,0025 s (erreur de vitesse = ΔA × décalage). Choisir des cas dont l'arrêt n'est pas aligné sur la grille pour mesurer un
  ordre. Sans résistance, les phases sont des polynômes de degré 2 : Verlet et RK4 exacts, rien à mesurer (ajouter k > 0). L'événement doit ignorer un
  départ exactement sur la surface (g = 0) : sinon faux arrêt à chaque redémarrage.
- **Chocs / événements** : une fonction d'événement qui vaut EXACTEMENT 0 au départ (balle posée sur le sol juste après un rebond) ne voit pas un vol entier contenu dans
  un pas : la balle traverse le sol (observé : y = −26 m dès que les rebonds deviennent plus courts que le pas, vers t = 3,7 s). Remèdes en place : g = max(y, 0) (nulle
  après le choc), puis un pas de décollage d'1/1000 du temps de montée (`launch_`). Même idée pour les billes : g = max(distance, 0), sinon faux événement au redémarrage.
  Le seuil d'arrêt doit rester ≥ 1e-8 m/s : en dessous `ProjectileProblem::landingTime` (bracket depuis 1e-9 s) renvoie 0. Ne pas asserter le nombre de rebonds du
  modèle naïf (24, 52, 99 selon dt) ni les erreurs d'ordre mesurées APRÈS le repos (les rapports « x4 » y sont trompeurs).
- **Hertz / berceau** : une force continue n'a pas besoin d'événement, mais elle n'est pas lisse en δ = 0 : l'ordre de RK4 tombe vers 2,5 et ses rapports d'erreur par doublement
  sont erratiques (ne jamais asserter « x16 » ni « x5,6 » : utiliser une pente sur 200 → 3200 pas et des bornes larges). Le début du contact doit être décalé (`gap` non
  aligné sur la grille) pour mesurer un ordre ; sinon les erreurs sont celles d'une phase alignée (0 : ordre apparent 2,38). **Euler symplectique et Verlet donnent la même erreur ici parce que l'état
  initial et l'état final n'ont pas de force** (conjugaison par un demi-pas de vitesse) : ne pas en conclure que Euler symplectique est d'ordre 2 en général. Mesurer δ_max et T
  par événement (vitesse relative nulle, puis compression nulle), jamais en échantillonnant le maximum aux pas de RK45 (les pas adaptatifs sont longs : erreur d'échantillonnage non négligeable). La fin de la collision est
  `collisionOver` : plus aucune compression ET vitesses rangées (gauche pas plus rapide que droite) ; la compression seule ne suffit pas. L'ordre de résolution des impulsions
  ne change rien dans le berceau (le brouillon de M5c affirmait le contraire : faux, mesuré). `cradleError` recalcule la référence RK45 à chaque appel (28 appels pour l'étude de convergence : sans à-coup visible, durée non mesurée).
- **Corps rigide (M6)** : `Quaternion::rotate` n'est valable que pour un quaternion UNITAIRE (sa forme développée diffère de q v q* dès que |q| ≠ 1) : pour dessiner un état dont la norme a dérivé, utiliser le produit complet
  q v q* (le solide gonfle de |q|², c'est l'erreur à montrer) ; pour les diagnostics (E, L) normaliser d'abord. **Garder |q| = 1 ne suffit pas** : Lie-Heun l'assure et dérive en énergie et en L (mesuré) ; la structure symplectique
  est ce qui conserve L. `ellipticK(sqrt(1 − m))` perd la précision pour m petit (nome faux de 0,5 % à m = 1e-14) : le nome se forme avec K = π / (2 AGM(1, k')) et K' = π / (2 AGM(1, k)) directement. Appeler `reference(t)` dans une
  boucle sur t ré-intègre depuis 0 à chaque fois (coût QUADRATIQUE : un seul test a pris 76 s en Debug) : intégrer une fois et échantillonner (`walkReference` dans les tests). L'ODE de la rotation est lisse : les ordres sont propres
  (×16,06 pour RK4, ×4,00 pour les schémas d'ordre 2), contrairement au contact de Hertz de M5c. Un test de l'« ordre » de Euler n'est asymptotique qu'à pas fin (avec le couple de la toupie, 50 pas donnent NaN). La solution de Jacobi exige
  I1 < I2 < I3 et un départ (a, 0, c) avec a, c > 0 ; sans ω2(0) = 0 il faudrait l'intégrale elliptique incomplète de première espèce pour la phase (non écrite). Hypothèse fausse corrigée : une toupie sous le spin critique ne tombe
  pas (grande nutation, elle remonte).
- **Interface** : la police Segoe UI n'a pas ∇ ni ∝ (affichés « � ») ; ∂, ᵀ, ω, √, ≈, Δ passent. Plus de 3 colonnes numériques dans « Invariants » (~400 px)
  écrasent la colonne « Méthode » : scinder en deux tableaux. `TextDisabled` ne passe pas à la ligne : utiliser `TextWrapped` colorée pour les notes.
- **PowerShell 5.1** : `Get-Content -Raw | Set-Content -Encoding utf8` ré-encode les accents (mojibake) et ajoute un BOM. Éditer avec l'outil Edit, ou avec
  `python -I` et `encoding="utf-8", newline=""`. En cas d'accident : `git checkout -- <fichier>`.

**Outils**
- Toolchain réelle : MSYS2 **mingw64** (GCC 16.1), pas `ucrt64` (non installé). Configure : `cmake -S . -B build -G Ninja`.
- ImPlot v1.0 : `PlotLine(label, xs, ys, n, ImPlotSpec(...))`, pas d'argument `offset` séparé (`ImPlotProp_Offset`). GLFW 3.4 est un tag, pas
  une branche : `git clone --depth 1 --branch 3.4 ...` puis `git submodule add`.
- Dans l'outil Bash, un **gros heredoc** (scripts Python, C++) peut échouer sur le quoting : écrire le script dans un fichier (outil
  Write, dans le dossier temporaire de la session) puis l'exécuter avec `python -I fichier.py`. Ne pas utiliser `/tmp`.
- Captures d'écran Windows : le processus PowerShell doit être déclaré DPI-aware (`SetProcessDPIAware`), sinon la capture est tronquée
  sur un écran mis à l'échelle. `tools/screenshot.ps1` le fait. Écran de l'utilisateur : 1920×1080 physique.
- ImGui crée `imgui.ini` dans le dossier courant (ignoré par git) ; s'il n'existe pas, la disposition par défaut est reconstruite.
- Avertissements « CRLF will be replaced by LF » : normaux, `.gitattributes` force LF dans le dépôt.
- Les chaînes de l'interface sont en UTF-8 avec lettres grecques (ω, ζ, θ, λ…) : elles s'affichent grâce à Segoe UI (Windows) ou DejaVu
  (Linux). Sans ces polices, la police intégrée d'ImGui n'affiche pas le grec.
- M5c : ℓ (U+2113) et ≪ (U+226A) remplacés par « L » et un texte par précaution (non testés) ; ¼, α, ≤, ≈, ², ½, −, →, ↔ passent. Un schéma qui diverge (Euler explicite à 10 pas par contact) met des
  valeurs infinies dans l'état : détecter (`diverged_`), arrêter, ne pas dessiner, ne pas échantillonner. Le ralenti du choc se décide par image (`slowMotion()`), pas par pas de temps.
  `tests/test_core.cpp` est en CRLF dans la copie de travail (un `replace` Python doit utiliser \r\n) ; `Application.cpp` a des fins de ligne mixtes : utiliser l'outil Edit ; les nouveaux fichiers sont en LF.
- M6 : ⊥ (U+22A5) NE s'affiche PAS (« � ») ; Ω, ⁻, ² et ³ passent (vérifié sur captures) ; β, τ, ¹, ⁵, ⁶ utilisés dans les textes des niveaux 4 à 6 mais non vérifiés à l'écran. Les gros heredocs `python - <<'EOF'` avec des apostrophes échouent dans l'outil Bash (« unexpected EOF ») : écrire le script avec l'outil Write dans
  le dossier temporaire puis `python -I script.py`. Ne PAS repasser `-G Ninja` sur un dossier `build` déjà configuré avec un autre générateur (erreur « Does not match the generator used previously ») : le dossier `build` de cette machine a été
  reconfiguré en MinGW Makefiles / Debug par autre chose que moi (les tests y durent 4,2 s contre 1 s en Release) ; construire avec `cmake --build build` sans `-G`. Pour une vérification Release : dossier jetable dans le dossier temporaire de la session.

## 5. Vérifier avant de commiter

```bash
cmake -S . -B build -G Ninja && cmake --build build           # aucun avertissement attendu (-Wall -Wextra -Wpedantic)
ctest --test-dir build --output-on-failure                     # "test_core : OK"
for s in 1 2 3 4 5 6 7 8 9 10; do ./build/physicslab --smoke-test --sim $s --level 6; done   # démarrage hors écran, GL 4.5
```
Build Debug propre depuis zéro de temps en temps (`-DCMAKE_BUILD_TYPE=Debug` dans un dossier jetable : vérifie les `assert`).
Vérification visuelle (PowerShell) : `.\tools\screenshot.ps1 -Level 5 -Sim 4 -Wait 8` puis ouvrir le PNG indiqué. Regarder au moins les
niveaux 1, 3, 5 et 6 d'un nouveau module, et un ancien module pour la non-régression. `-Clicks "x,y;x,y"` simule des clics.

Options de l'application : `--level 1..6`, `--sim 1..10` (1 = M1, 2 = M2, 3 = M3 pendule simple, 4 = M3b pendule double, 5 = M4a Kepler,
6 = M4b N corps, 7 = M5a frottement sec, 8 = M5b chocs et rebonds, 9 = M5c berceau de Newton, 10 = M6 corps rigide), `--smoke-test`. Navigation 3D : clic gauche tourner, clic droit/milieu déplacer, molette zoomer.

## 6. Suite : M7 N-corps GPU (brouillon, à re-présenter brièvement au début de la prochaine conversation, puis commencer par l'étape 1 si l'utilisateur valide)

Brouillon rédigé de mémoire, à valider par le calcul et par des mesures (règle : ne rien affirmer sans l'avoir mesuré). Dernier module de la Mécanique.

**Le problème.** Reprendre le calcul O(N²) de `nbody::accelerations` (M4b, CPU double) dans un compute shader OpenGL 4.5, en `float`, et mesurer l'écart au CPU : « le CPU vérifie le GPU ». Le contrat est l'interface existante
(positions, masses, n, G, adoucissement, accélérations). L'adoucissement de Plummer est obligatoire (en `float` une rencontre rapprochée sature).

**Algorithme.** Un thread par corps cible, boucle sur les sources par tuiles chargées en mémoire partagée (taille du groupe de travail), accumulation de a_i = G sum_j m_j (r_j − r_i) / (|r|² + eps²)^(3/2) ; intégration kick-drift-kick (Verlet)
dans le shader sur des SSBO (double tampon des positions). Centrer les positions sur le barycentre pour limiter l'annulation dans r_j − r_i.

**Difficultés (à mesurer, pas à affirmer).**
- Précision `float` (mantisse 24 bits, 6e-8 relatif) : erreur d'accélération en fonction de N ; dérive d'énergie plus forte qu'en `double` ; ordre de sommation non déterministe d'un GPU à l'autre.
- Le chaos limite la comparaison : après un temps de l'ordre de 1/λ (M4b : λ ≈ 0,8 pour l'amas de 6 corps, amplification ×6,6e4 sur 14 unités de temps) une trajectoire `float` et une trajectoire `double` n'ont plus rien de commun. Comparer
  les ACCÉLÉRATIONS à état fixé (erreur relative corps par corps), puis des horizons courts et des grandeurs d'ensemble (E, P, L), pas des trajectoires longues.
- Les forces par paire sont antisymétriques à l'arrondi près seulement : l'impulsion totale dérive (à mesurer, ~1e-7 relatif attendu par pas).
- Les tests de `test_core` n'ont pas de contexte GL : prévoir une option `--gpu-test` de l'application (fenêtre cachée GL 4.5 comme `--smoke-test`) qui compare GPU et CPU et renvoie un code de sortie. Vérifier que le glad généré expose
  `glDispatchCompute`, `glMemoryBarrier`, les SSBO et les requêtes de temps (GL 4.5 core devrait suffire : à contrôler).
- Mesurer le temps avec `GL_TIME_ELAPSED` : N = 1e3 ... 1e5 (1e5² = 1e10 interactions par pas), débit en interactions par seconde contre le CPU `double`. Barnes-Hut hors périmètre.

**Simulation (`--sim 11`).** Amas / galaxie de N réglable (de 100 à plusieurs dizaines de milliers : à mesurer), bascule CPU / GPU, panneau d'erreur (distribution des écarts d'accélération, dérive de E et P), temps par pas ; niveaux 1-6 (au niveau 1 :
« l'ordinateur calcule toutes les paires de corps, la carte graphique les fait en parallèle »). Rendu : `Renderer::draw(Points)` avec renvoi des positions au CPU à chaque image au début ; un rendu direct depuis le SSBO demanderait de modifier le renderer.

**Plan.** (1) Cœur et tests d'abord : enveloppe de programme de calcul et SSBO dans `render/` (a besoin d'un contexte GL : hors `physicslab_core`), shader `shaders/nbody.comp`, comparaison GPU/CPU à état fixé via `--gpu-test`. (2) Interface. (3) Commit, push, passation.
**Fin de domaine Mécanique (après M7)** : devlog complet (`docs/devlog/`), prompt de reprise court, cours compilé en PDF ou Word pour non-initiés (`docs/cours/`), « roue des domaines » pour choisir le module suivant, merge de
`module/mecanique` dans `main` et tag `mecanique-1`.

**À réutiliser :** `nbody::accelerations` / `NBodyProblem` (CPU de référence), `SolverSet`, `StepClock`, `drawResultTable`, `wrapped`, `CradleModule` et `RigidBodyModule` comme exemples récents de modules (modèles côte à côte, détection de divergence,
fenêtre « Analyse »), `tools/screenshot.ps1`.

## 7. Dette technique et idées

- `computeConvergence()` est dupliqué dans 10 modules (variantes : erreur à `tEnd` fixe ; pente ajustée là où l'erreur < 0,1 pour M4 ; tous les points pour M5a) :
  à factoriser (fonction générique prenant une fonction d'erreur).
- M5a : tableau des modèles d'orbite/chocs sans sélecteur de solveur par modèle (les 5 solveurs partagent le modèle choisi) ; le plan est dessiné en fil de fer,
  une rampe pleine (triangles) donnerait une meilleure lecture. Frottement de roulement, frottement visqueux non linéaire (quadratique) : absents.
- Euler explicite et symplectique restent mal comparés à pas égal (ordre 1 chacun mais constantes très différentes) ; le pas par particule ou adaptatif
  n'existe pas pour les N corps (rencontres rapprochées : l'énergie saute). À reprendre avec M7 (GPU float).
- Orbites non liées (e ≥ 1, hyperboles) absentes de Kepler ; triangle de Lagrange à masses inégales (stable pour les Troyens) absent de NBody.
- `ProjectileModule` n'utilise pas encore `StepClock` ni `World` de la même façon que les autres modules.
- Absents : type `Tensor` (relativité), backend GPU / compute shaders, export de données (niveau 6), tableau périodique, constantes
  fondamentales et convertisseur d'unités en interface.
- M5c : billes suspendues (force de rappel g/L, vrai berceau) et contact linéaire « ressort » pour comparaison : non faits (bonus du plan) ; chaîne 1D seulement (pas de choc oblique) ;
  billes dessinées en fil de fer (cercles) ; un seul modèle d'amortissement (Hunt-Crossley) ; l'issue de référence de Hertz est calculée à chaque `reset()` (RK45 serré, quelques ms) ;
  l'étude de convergence du module recalcule 28 références à chaque changement de paramètre quand « Analyse » est visible ; le texte du niveau 6 cite Nesterenko (1983, onde solitaire de
  la chaîne de Hertz) de mémoire : non vérifié ici ; les textes des niveaux 2 à 4 n'ont pas été relus à l'écran. `tests/test_core.cpp` dépasse 2500 lignes : à scinder par domaine.
- M6 : boîte de dimensions fixes (3 × 2 × 1) ; solution de Jacobi seulement pour un départ (a, 0, c) (pas de phase quelconque : intégrale elliptique incomplète à écrire) ; pas d'orientation exacte du solide asymétrique (fonctions thêta) ;
  toupie à pivot fixe seulement (pas de toupie sur une table : contact, frottement, pointe qui glisse) ; pas de choc de corps rigides (cône de Coulomb, paradoxe de Painlevé) ; découpage d'ordre 2 seulement (pas de composition d'ordre 4 ni
  de RKMK d'ordre 4, alors que RK4 + renormalisation est d'ordre 4 mais perd L) ; textes des niveaux 2 et 3 non relus à l'écran ; valeurs extrêmes des curseurs et chemin « précession régulière impossible » non déclenchés dans l'application ;
  corps dessinés en fil de fer. `tests/test_core.cpp` : plus de 3100 lignes (à scinder par domaine).
- Une fenêtre console s'ouvre à côté de l'application (à masquer en Release sous Windows).
- Dans le tableau d'invariants du niveau 6, certains libellés sont abrégés (« Euler sympl. »).
- Identité Git : `user.name` vaut `Noa-biy11` alors que le login GitHub est `Noa11-biy` (à corriger si l'utilisateur le souhaite).
- Copyright du fichier LICENSE : « Undermania » (nom du profil GitHub) ; à confirmer.

## 8. Prompt de reprise (à coller dans la nouvelle conversation)

> Reprends le projet PhysicsLab 3D dans ce dossier. Lis d'abord `docs/PASSATION.md` en entier (rôle, règles, état, architecture, pièges),
> puis `README.md`. Vérifie que ça compile et que les tests passent (section 5), puis présente-moi la feuille de route courte du
> module M7 (N-corps GPU en compute shader, écart CPU/GPU ; M6 corps rigide est terminé ; la section 6 en donne le brouillon) avant de coder.
> Réponses courtes pendant le module, en français. Après M7 : livrables de fin de domaine Mécanique (devlog, cours, roue des domaines, merge, tag).
