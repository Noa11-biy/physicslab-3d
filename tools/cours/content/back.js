const D = require("docx");
const B = require("../lib/blocks");

// Glossaire (ordre alphabétique)
const glossary = [
  ["Accélération", "Variation de la vitesse par unité de temps (m/s²). La deuxième loi de Newton la relie à la force : $a = F/m$."],
  ["Adoucissement (de Plummer)", "Petite longueur ε ajoutée dans la force gravitationnelle, $1/(r² + ε²)^{3/2}$, pour éviter une force infinie quand deux étoiles passent très près l’une de l’autre."],
  ["Chaos", "Sensibilité extrême aux conditions initiales : un écart minuscule au départ est amplifié exponentiellement, ce qui limite la prévision."],
  ["Coefficient de restitution", "Rapport $e$ entre les vitesses relatives après et avant un choc ($e = 1$ : choc parfaitement élastique, $e = 0$ : choc mou)."],
  ["Compute shader", "Petit programme exécuté en parallèle par des milliers de calculateurs de la carte graphique, utilisé ici pour calculer les forces entre étoiles."],
  ["Conservation", "Une grandeur (énergie, impulsion, moment cinétique) qui ne change pas au cours du mouvement ; ses variations numériques mesurent l’erreur d’un calcul."],
  ["Convergence (ordre de)", "Exposant $p$ tel que l’erreur globale se comporte comme $Δt^p$. Ordre 1 : Euler ; 2 : Verlet ; 4 : RK4."],
  ["CPU, processeur", "Cerveau généraliste de l’ordinateur ; calcule en double précision dans ce cours (la référence)."],
  ["Double précision (double)", "Nombres à 64 bits, environ 16 chiffres significatifs (erreur d’arrondi relative de l’ordre de $10^{−16}$)."],
  ["Énergie cinétique, potentielle", "Énergie de mouvement ($½ m v²$) et énergie stockée dans une position (poids, ressort, gravitation). Leur somme se conserve sans frottement."],
  ["Équation différentielle (EDO)", "Équation qui relie une grandeur à ses dérivées par rapport au temps, par exemple $m x'' = −k x$."],
  ["Erreur locale, globale", "Erreur commise en un pas, et erreur accumulée sur toute la simulation."],
  ["Espace des phases", "Espace dont les coordonnées sont les positions et les vitesses (ou impulsions) ; l’état d’un système y est un point."],
  ["Euler (méthode d’)", "La plus simple des méthodes pas à pas : on avance avec la pente au début du pas. Ordre 1."],
  ["Événement", "Instant où une condition change (rebond, arrêt d’un bloc). On le localise par bissection plutôt que de le subir."],
  ["Exposant de Lyapunov", "Taux $λ$ de croissance exponentielle d’un petit écart, $δ(t) ≈ δ_0 e^{λ t}$ ; positif pour un système chaotique."],
  ["Float, simple précision", "Nombres à 32 bits, environ 7 chiffres significatifs (erreur relative de l’ordre de $6 × 10^{−8}$) ; format des cartes graphiques."],
  ["Frottement sec, visqueux", "Force de Coulomb indépendante de la vitesse ($μ N$), ou proportionnelle à la vitesse ($k v$)."],
  ["GPU, carte graphique", "Processeur massivement parallèle ; ici il calcule les forces entre étoiles en float."],
  ["Hertz (contact de)", "Modèle de contact élastique entre deux sphères : force $k δ^{3/2}$ en fonction de l’écrasement $δ$."],
  ["Impulsion", "Quantité de mouvement $P = Σ m v$ ; se conserve si aucune force extérieure n’agit."],
  ["Kick-drift-kick", "Forme du schéma de Verlet des vitesses : demi-coup sur la vitesse, dérive de la position, demi-coup."],
  ["Moment cinétique", "Grandeur de rotation $L = Σ m r × v$ ; se conserve pour une force centrale."],
  ["Moment d’inertie", "Résistance d’un solide à la mise en rotation autour d’un axe, $I = Σ m r²$ ; dépend de l’axe."],
  ["Nutation", "Oscillation de l’inclinaison de l’axe d’une toupie, superposée à sa précession."],
  ["Pas de temps", "Durée $Δt$ d’un petit calcul élémentaire d’une simulation."],
  ["Précession", "Rotation lente de l’axe d’un objet (toupie) ou de l’orbite d’une planète."],
  ["Quaternion", "Quadruplet de nombres qui représente une rotation sans ambiguïté ; il doit rester de norme 1."],
  ["Résonance", "Réponse très forte d’un système excité à sa fréquence propre."],
  ["RK4, RK45", "Méthodes de Runge-Kutta : d’ordre 4 à pas fixe, ou paire d’ordres 5 et 4 à pas adaptatif (Dormand et Prince)."],
  ["Solution exacte", "Formule qui donne le mouvement sans approximation ; rare, elle sert de juge aux méthodes numériques."],
  ["Stabilité", "Une méthode est stable si les erreurs ne s’amplifient pas au fil des pas ; Euler sur $y' = −y/τ$ l’est si $Δt < 2τ$."],
  ["Sphère de Plummer", "Modèle d’amas d’étoiles à l’équilibre, de densité $(1 + r²/a²)^{−5/2}$, utilisé pour les conditions initiales du chapitre 7."],
  ["Symplectique", "Se dit d’un schéma qui conserve la structure géométrique du mouvement hamiltonien : l’énergie oscille sans dériver."],
  ["TDR", "Réinitialisation du pilote graphique de Windows quand une commande dure plus de 2 s ; impose de découper les calculs."],
  ["Tuile (tile)", "Paquet de 256 étoiles chargé une fois en mémoire rapide et partagé par tous les calculateurs d’un groupe."],
  ["Verlet (des vitesses)", "Schéma d’ordre 2, symplectique, très employé en dynamique moléculaire et en astronomie."],
  ["Viriel (équilibre du)", "Relation $2T = −U$ entre énergie cinétique $T$ et potentielle $U$ d’un amas à l’équilibre."],
  ["Zénon (accumulation de)", "Une balle qui rebondit une infinité de fois en un temps fini : $t_0 (1+e)/(1−e)$."],
];

const formulas = [
  ["1", "Euler, énergie fabriquée en chute libre", "$ΔE = ½ m g² t Δt$", "100 J pour $m = 2$ kg, $t = 2$ s, $Δt = 0,5$ s"],
  ["1", "Vitesse limite, temps caractéristique", "$v_{∞} = m g / k$ ; $τ = m / k$", "40 m/s et 4 s (parachutiste de 80 kg, $k = 20$ kg/s)"],
  ["2", "Période et énergie du ressort", "$T = 2π sqrt{m/k}$ ; $E = ½ k A²$", "0,314 s et 0,16 J ($m = 0,5$ kg, $k = 200$ N/m, $A = 4$ cm)"],
  ["2", "Énergie par pas d’Euler (oscillateur)", "×$(1 + ω_0² Δt²)$", "1,48 par période ($T = 1$ s, $Δt = 0,01$ s)"],
  ["2", "Énergie par pas de RK4", "×$(1 − y⁶/72 + y⁸/576)$, $y = ω_0 Δt$", "0,9921 après 30 périodes à 20 pas par période (mesuré : 0,9921)"],
  ["3", "Période du pendule", "$T_0 = 2π sqrt{L/g}$ ; $T = 4 K(k)/ω_0$, $k = sin(θ_0 / 2)$", "$T/T_0 = 1,0732$ à 60°"],
  ["3", "Vitesse en bas du pendule", "$v = sqrt{2 g L (1 − cos θ_0)}$", "3,84 m/s ($L = 1,5$ m, 60°)"],
  ["3", "Amplification d’un écart (chaos)", "$δ(t) = δ_0 e^{λ t}$", "$λ ≈ 1$ /s pour le pendule double ; ×4400 en 8 s"],
  ["4", "Troisième loi de Kepler", "$T² = (4π² / G M) a³$", "Mars : 1,88 an ; Jupiter : 11,87 ans"],
  ["4", "Vitesse circulaire, libération", "$v = sqrt{G M / r}$ ; $v_{lib} = sqrt{2} v$", "7,67 km/s à 400 km d’altitude ; période 92,6 min"],
  ["4", "Triangle de Lagrange", "$ω² = 3 G m / s³$", "$T = 2π/sqrt{3} = 3,628$ ($G = m = s = 1$)"],
  ["4", "Précession numérique de Verlet", "$Δω = −frac{π}{8} frac{G M h²}{a³}$\n$× frac{4 + e²}{(1 − e²)³}$", "−0,2237° par orbite à 200 pas/orbite, $e = 0,5$ (mesuré : −0,2232°)"],
  ["5", "Plan incliné avec frottement", "$a = g (sin θ − μ cos θ)$ ; glisse si $tan θ > μ$", "0,59 m/s² ($θ = 20°$, $μ = 0,3$)"],
  ["5", "Rebonds : durée totale", "$t_0 (1 + e)/(1 − e)$, $t_0 = sqrt{2 h_0 / g}$", "5,748 s ($h_0 = 2$ m, $e = 0,8$)"],
  ["5", "Contact de Hertz (deux billes égales)", "$δ_{max} = (5 m_r v² / 4k)^{2/5}$ ; $T = 2,943275 δ_{max} / v$", "0,0208138 et 0,0612608 ($k = 10⁴$, $m = 1$, $v = 1$)"],
  ["6", "Moments d’inertie", "disque $½ M R²$ ; sphère pleine $2/5 M R²$", "0,010 et 0,008 kg·m² ($M = 2$ kg, $R = 0,1$ m)"],
  ["6", "Précession rapide d’une toupie", "$Ω ≈ m g l / (I_3 ω_3)$", "0,490 rad/s (exacte : 0,517 rad/s à $θ = 0,6$)"],
  ["6", "Spin critique de la toupie endormie", "$ω_c = 2 sqrt{I_1 m g l} / I_3$", "12,128 rad/s"],
  ["6", "Axe intermédiaire", "$λ = ω_2 sqrt{frac{(I_2 − I_1)(I_3 − I_2)}{I_1 I_3}}$", "0,5774 /s pour $I = (1, 2, 3)$, $ω_2 = 1$"],
  ["7", "Nombre de paires", "$N(N − 1)/2$", "1 999 000 pour $N = 2000$"],
  ["7", "Erreur d’une somme en float", " ≈ $u sqrt{N}$, $u = 2⁻²⁴ = 6 × 10⁻⁸$", "mesuré : 1,6·10⁻⁸ ($N = 3$) à 1,05·10⁻⁶ ($N = 5000$)"],
];

function sections({ chapters, heading, makeRenderer, T }) {
  const out = [];
  const R = (n) => makeRenderer({ fig: { chapter: n, n: 0 } });

  // ---------------------------------------------------------------- conclusion
  {
    const r = R("C");
    const kids = [heading(1, "Conclusion · dix idées à retenir")];
    kids.push(...r.render(B.lead("Sept chapitres, onze simulations et une seule méthode : dériver, résoudre exactement, discrétiser, programmer, valider, expliquer. Voici ce qui reste quand on a tout oublié.")));
    kids.push(
      ...r.render(
        B.ol([
          "**Un calcul numérique avance pas à pas.** Plus le pas est fin, plus on est précis et plus c’est cher (chapitre 1).",
          "**L’ordre d’une méthode** dit comment l’erreur baisse avec le pas : ÷2 en ordre 1, ÷4 en ordre 2, ÷16 en ordre 4.",
          "**Une solution exacte est un trésor** : même rare, elle sert de juge à tous les calculs (chapitres 1, 2, 3, 4, 6).",
          "**Les grandeurs conservées sont des détecteurs d’erreurs** : énergie, impulsion, moment cinétique (chapitres 2 et 4).",
          "**Précis n’est pas conservatif** : RK4 est plus précis que Verlet à pas égal mais perd de l’énergie lentement ; Verlet, symplectique, oscille sans dériver (chapitres 2 et 4).",
          "**Le chaos limite la prévision** : un écart $δ_0$ devient $δ$ en un temps $ln(δ/δ_0)/λ$ (chapitre 3).",
          "**Quand la physique est discontinue, il faut détecter l’événement** ; sinon tous les schémas retombent à l’ordre 1 (chapitre 5).",
          "**En rotation, conserver la norme ne suffit pas** : c’est la structure du schéma qui conserve le moment cinétique (chapitre 6).",
          "**La précision finie compte** : une somme de $N$ termes en float se trompe de $u sqrt{N}$ ; on centre les données en double avant de les convertir (chapitre 7).",
          "**Mesurer avant d’affirmer** : une toupie sous son spin critique ne tombe pas, RK4 n’est pas d’ordre 4 sur un contact de Hertz, l’ordre de résolution des chocs du berceau ne change rien.",
        ])
      )
    );
    kids.push(
      ...r.render(
        B.box("plus", "Et ensuite ?", [
          B.p("Le domaine **Mécanique** est terminé ; il a posé les outils (solveurs, invariants, chaos, événements, GPU) dont les autres domaines auront besoin : Ondes (acoustique, optique), Thermodynamique, Électrodynamique, Fluides, Plasma, Atomique, quantique et nucléaire, Relativité, astrophysique et cosmologie, puis les domaines appliqués. Chacun reprendra la même méthode en six étapes et les mêmes six niveaux."),
        ])
      )
    );
    out.push({ header: "Conclusion", children: kids });
  }

  // ---------------------------------------------------------------- annexe A : corrigés
  {
    const r = R("A");
    const kids = [heading(1, "Annexe A · Corrigés des exercices")];
    kids.push(...r.render({ k: "p", x: "Chaque corrigé donne le raisonnement et le résultat. Les valeurs numériques ont été calculées par programme et, quand c’est possible, recoupées avec des mesures du logiciel.", opt: { italic: true, color: T.colors.grey } }));
    for (const ch of chapters) {
      kids.push(heading(2, `Chapitre ${ch.num} · ${ch.title}`, { color: T.colors.navy }));
      ch.exercises.forEach((ex, i) => {
        const lv = T.levels[ex.level - 1];
        kids.push(
          new D.Paragraph({
            keepNext: true,
            spacing: { before: 200, after: 60 },
            children: [
              new D.TextRun({ text: `Exercice ${ch.num}.${i + 1}`, bold: true, color: T.colors.navy, size: 22 }),
              new D.TextRun({ text: `   N${ex.level} · ${lv.name}`, bold: true, color: lv.color, size: 18 }),
              ...(ex.title ? [new D.TextRun({ text: `   ${ex.title}`, italics: true, color: T.colors.grey, size: 20 })] : []),
            ],
          })
        );
        for (const b of ex.solution) kids.push(...r.render(b));
      });
    }
    out.push({ header: "Annexe A", children: kids });
  }

  // ---------------------------------------------------------------- annexe B : glossaire
  {
    const r = R("B");
    const kids = [heading(1, "Annexe B · Glossaire")];
    kids.push(...r.render(B.tbl(["Terme", "Définition"], glossary, [2700, 6598], { boldFirst: true, size: 19 })));
    out.push({ header: "Annexe B", children: kids });
  }

  // ---------------------------------------------------------------- annexe C : formulaire
  {
    const r = R("C");
    const kids = [heading(1, "Annexe C · Formulaire et valeurs de référence")];
    kids.push(...r.render({ k: "p", x: "Les formules principales de chaque chapitre, avec une valeur de référence vérifiable dans le logiciel ou par un calcul à la main.", opt: { italic: true, color: T.colors.grey } }));
    kids.push(...r.render(B.tbl(["Chap.", "Grandeur", "Formule", "Valeur de référence"], formulas, [760, 2000, 3340, 3198], { align: ["c", "l", "l", "l"], size: 18 })));
    kids.push(heading(2, "Constantes utilisées", { color: T.colors.navy }));
    kids.push(
      ...r.render(
        B.tbl(
          ["Constante", "Valeur"],
          [
            ["Pesanteur terrestre $g$", "9,80665 m/s² (arrondie à 9,81, ou à 10 dans les exercices à la main)"],
            ["Constante de gravitation $G$", "6,674 × 10⁻¹¹ m³/(kg·s²) (valeur normalisée $G = 1$ dans les simulations de N corps)"],
            ["$G M$ de la Terre", "3,986004 × 10¹⁴ m³/s²"],
            ["$G M$ du Soleil", "1,32712 × 10²⁰ m³/s²"],
            ["Unité astronomique", "1,495979 × 10¹¹ m"],
            ["Précision d’un float / d’un double", "6 × 10⁻⁸ (2⁻²⁴) / 1,1 × 10⁻¹⁶ (2⁻⁵³)"],
          ],
          [3600, 5698],
          { boldFirst: true, size: 19 }
        )
      )
    );
    out.push({ header: "Annexe C", children: kids });
  }

  // ---------------------------------------------------------------- annexe D : le logiciel
  {
    const r = R("D");
    const kids = [heading(1, "Annexe D · Utiliser le logiciel")];
    kids.push(...r.render(B.p("PhysicsLab 3D est un programme en ligne de commande qui ouvre une fenêtre. Voici les options utiles :")));
    kids.push(
      ...r.render(
        B.tbl(
          ["Option", "Effet"],
          [
            ["{{--level 1..6}}", "Choisit le niveau pédagogique au démarrage (1 vulgarisation, 2 intéressé, 3 collège, 4 lycée, 5 étudiant, 6 doctorant). Il se change aussi dans le menu."],
            ["{{--sim 1..11}}", "Choisit la simulation au démarrage (tableau ci-dessous). Elle se change aussi dans le menu Simulation."],
            ["{{--smoke-test}}", "Ouvre une fenêtre cachée, dessine quelques images et quitte : vérifie que tout démarre."],
            ["{{--gpu-test}}", "Compare le calcul de la carte graphique au processeur (précision, temps, intégration) et donne un code de sortie ; option {{--gpu-max-n N}} pour la taille maximale."],
          ],
          [2400, 6898],
          { size: 19 }
        )
      )
    );
    kids.push(...r.render(B.p("Navigation dans la scène 3D : clic gauche pour tourner, clic droit ou milieu pour déplacer, molette pour zoomer.", { before: 120 })));
    kids.push(heading(2, "Les onze simulations", { color: T.colors.navy }));
    kids.push(
      ...r.render(
        B.tbl(
          ["--sim", "Module", "Contenu", "Chapitre"],
          [
            ["1", "M1", "Projectile avec frottement, cinq méthodes contre la solution exacte", "1"],
            ["2", "M2", "Ressort-masse libre, amorti, forcé ; résonance", "2"],
            ["3", "M3", "Pendule simple, solution exacte elliptique", "3"],
            ["4", "M3b", "Pendule double : chaos, exposant de Lyapunov", "3"],
            ["5", "M4a", "Orbite de Kepler, précession numérique", "4"],
            ["6", "M4b", "Problème à N corps : le huit, le triangle de Lagrange, un amas", "4"],
            ["7", "M5a", "Frottement sec de Coulomb sur un plan incliné", "5"],
            ["8", "M5b", "Chocs et rebonds, accumulation de Zénon", "5"],
            ["9", "M5c", "Berceau de Newton : contact de Hertz contre impulsions", "5"],
            ["10", "M6", "Corps rigide : solides libres, toupie de Lagrange", "6"],
            ["11", "M7", "N corps sur carte graphique, écart processeur / GPU", "7"],
          ],
          [900, 1000, 6200, 1198],
          { align: ["c", "c", "l", "c"], size: 19 }
        )
      )
    );
    kids.push(heading(2, "Le sélecteur de niveau", { color: T.colors.navy }));
    kids.push(
      ...r.render(
        B.p("Le menu **Niveau**, en haut de la fenêtre, change à la fois le texte de la fenêtre **Explication**, les curseurs offerts dans **Simulation**, les colonnes du tableau d’**Invariants** et les courbes de **Graphes** et d’**Analyse**. Les niveaux de ce cours correspondent un à un à ceux du logiciel.")
      )
    );
    kids.push(
      ...r.render(
        B.box("logiciel", "Compiler le logiciel", [
          B.p("Il faut un compilateur C++17, CMake 3.20 ou plus, et un pilote graphique OpenGL 4.5 :"),
          B.p("{{git clone --recurse-submodules <dépôt PhysicsLab 3D>}}", { size: 20 }),
          B.p("{{cmake -S . -B build -G Ninja}}  puis  {{cmake --build build}}", { size: 20 }),
          B.p("{{ctest --test-dir build --output-on-failure}}  vérifie les 76 tests automatiques.", { size: 20 }),
        ])
      )
    );
    out.push({ header: "Annexe D", children: kids });
  }

  // ---------------------------------------------------------------- annexe E : références
  {
    const r = R("E");
    const kids = [heading(1, "Annexe E · Pour aller plus loin")];
    kids.push(...r.render(B.p("Quelques ouvrages et articles qui ont servi de base ou qui prolongent le cours.", { italic: true, color: T.colors.grey })));
    kids.push(
      ...r.render(
        B.ul([
          "L. Landau et E. Lifchitz, //Physique théorique, tome 1 : Mécanique//, Mir (solide libre, § 37 ; toupie ; petites oscillations).",
          "H. Goldstein, C. Poole et J. Safko, //Classical Mechanics//, Addison-Wesley (formulations de Lagrange et Hamilton, solide rigide).",
          "V. I. Arnold, //Méthodes mathématiques de la mécanique classique//, Mir (mécanique hamiltonienne, géométrie des orbites).",
          "E. Hairer, C. Lubich et G. Wanner, //Geometric Numerical Integration//, Springer (schémas symplectiques, analyse par erreur inverse).",
          "J. R. Dormand et P. J. Prince, « A family of embedded Runge-Kutta formulae », //Journal of Computational and Applied Mathematics// 6 (1980), 19-26 (méthode RK45).",
          "H. Poincaré, « Sur le problème des trois corps et les équations de la dynamique », //Acta Mathematica// 13 (1890) (premier théorème d’impossibilité et origine du chaos).",
          "C. Moore, « Braids in classical dynamics », //Physical Review Letters// 70 (1993), 3675 ; A. Chenciner et R. Montgomery, « A remarkable periodic solution of the three-body problem in the case of equal masses », //Annals of Mathematics// 152 (2000), 881-901 (le « huit »).",
          "S. J. Aarseth, M. Hénon et R. Wielen, //Astronomy and Astrophysics// 37 (1974), 183 (conditions initiales de Plummer) ; H. C. Plummer, //Monthly Notices of the RAS// 71 (1911), 460.",
          "H. Hertz, « Über die Berührung fester elastischer Körper », //Journal für die reine und angewandte Mathematik// 92 (1882), 156-171 (contact de Hertz).",
          "Documentation du dépôt PhysicsLab 3D : {{docs/PASSATION.md}} (état et pièges) et {{docs/devlog/}} (journal de développement, avec toutes les mesures).",
        ])
      )
    );
    out.push({ header: "Annexe E", children: kids });
  }
  return out;
}

module.exports = { sections };
