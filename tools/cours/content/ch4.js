const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 4,
  title: "Gravitation : de Kepler aux N corps",
  lead: "Deux astres qui s’attirent dansent une valse parfaitement prévisible depuis Newton. Avec trois astres il n’y a plus de formule, et c’est l’ordinateur qui doit danser à leur place.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** pourquoi les planètes tournent, comment l’on calcule une orbite, pourquoi un ordinateur peut faire « tourner » une ellipse qui devrait rester fixe, et pourquoi trois corps suffisent à rendre l’avenir imprévisible.",
      "**Les résultats marquants :** la précession numérique de l’orbite se **prédit par une formule** ; le « huit » de trois étoiles se poursuit pour toujours ; un amas de six étoiles est chaotique.",
      "**Dans le logiciel :** module M4, simulations 5 (orbite de Kepler, {{--sim 5}}) et 6 (problème à N corps, {{--sim 6}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`La Lune ne reste pas là-haut par magie : la gravité de la Terre la tire sans arrêt vers nous, elle « tombe » en permanence. Mais elle avance aussi très vite sur le côté, et la Terre étant ronde, elle manque sans cesse son but. Elle tombe en tournant. C’est l’idée du canon de Newton : tirez un boulet assez fort horizontalement du haut d’une montagne, et il fera le tour de la Terre.`),
      p(`Pour **deux** astres, Newton a su tout prédire : chaque planète décrit une ellipse, pour toujours, comme sur des rails invisibles. Mais ajoutez un troisième astre et le problème devient si difficile que **personne ne sait écrire la réponse**. Vers 1890 Henri Poincaré a même montré qu’il est souvent impossible de prévoir ce qui se passera très longtemps après.`),
      p(`Il existe pourtant de rares danses parfaites. La plus célèbre : le « **huit** ». Trois étoiles de même masse se poursuivent sans fin le long d’une unique courbe en forme de 8, sans jamais se heurter. Elle a été trouvée en 1993 par ordinateur, et c’est le logiciel qui nous la fait voir.`),
      fig("fig_kepler.png", "Une même orbite suivie par trois méthodes (calcul du logiciel, dix tours, 200 pas par tour). Euler : l’orbite gonfle et s’éloigne ; Verlet : la forme est préservée ; RK4 : ellipse presque parfaite. Énergie finale (exacte −0,5) : Euler −0,03 ; Verlet −0,499 ; RK4 −0,500005.", 1.0),
    ],
    2: [
      p(`Johannes Kepler a énoncé trois lois pour les planètes : elles décrivent des **ellipses** dont le Soleil occupe un foyer ; elles balaient des **aires égales en des temps égaux** (donc vont plus vite près du Soleil) ; et le carré de leur période est proportionnel au cube de la taille de leur orbite. Isaac Newton les a toutes expliquées d’un coup avec une seule idée : l’attraction entre deux masses diminue comme l’inverse du carré de la distance.`),
      p(`Le problème à **deux corps** se résout exactement. Avec un ordinateur on le simule pour apprendre à juger les méthodes, car on connaît la bonne réponse. On y découvre un phénomène inattendu : la **précession numérique**. Une ellipse devrait rester fixe, mais la méthode de Verlet la fait tourner lentement sur elle-même. On sait même **prédire** de combien (niveau 5).`),
      p(`À **trois corps ou plus**, plus de formule : le calcul pas à pas est le seul moyen. Le mouvement est en général **chaotique**. Le logiciel distingue un système chaotique d’un système simplement compliqué en lançant deux copies presque identiques : si l’écart entre elles est multiplié par plus de mille, le système est chaotique. Le « huit » ne l’est pas : c’est un mouvement stable.`),
    ],
    3: [
      p(`Deux masses $m_1$ et $m_2$ distantes de $d$ s’attirent avec une force :`),
      eq(`F = G frac{m_1 m_2}{d²}`),
      p(`où $G = 6,674 × 10^{−11}$ N·m²/kg². Pour une masse de 1 kg à la surface de la Terre (rayon 6371 km, masse 5,972 × 10²⁴ kg) on retrouve le poids : F = 9,82 N.`),
      p(`La troisième loi de Kepler s’écrit très simplement quand on mesure le temps en années et les distances en unités astronomiques (1 UA = distance Terre-Soleil) : **T² = a³**.`),
      tbl(
        ["Astre", "Demi-grand axe a", "Période (années)"],
        [
          ["Terre", "1 UA", "1 an"],
          ["Mars", "1,524 UA", "1,88 an"],
          ["Jupiter", "5,203 UA", "11,87 ans"],
          ["Comète imaginaire", "4 UA", "8 ans"],
        ],
        [3100, 3100, 3098],
        { align: ["l", "c", "c"] }
      ),
      box("exemple", "Exemple : l’apesanteur n’est pas l’absence de gravité", [
        p(`À 400 km d’altitude (la Station spatiale), la distance au centre de la Terre est de 6771 km au lieu de 6371 km, donc la force de gravité est $(6371/6771)² = 88,5 %$ de celle du sol : 8,69 N par kilogramme. Les astronautes ne flottent pas parce que la gravité a disparu, mais parce qu’ils **tombent** en permanence avec la station.`),
      ]),
    ],
    4: [
      p(`Sous forme vectorielle, l’accélération d’une planète autour d’un astre de masse $M$ est $a = − G M r / |r|³$. Par unité de masse, deux grandeurs se conservent : l’**énergie** et le **moment cinétique** $L = r × v$ (dont la conservation est la loi des aires) :`),
      eq(`E = frac{1}{2} v² − frac{G M}{r}`),
      p(`Pour une orbite circulaire de rayon $r$, la force gravitationnelle joue le rôle de force centripète, d’où :`),
      eq(`v = sqrt{frac{G M}{r}}, T = frac{2π r}{v}, v_{lib} = sqrt{2} v`),
      p(`$v_{lib}$ est la vitesse de libération. Pour une ellipse de demi-grand axe $a$, l’énergie vaut $E = −G M/(2a)$ : elle ne dépend pas de l’excentricité.`),
      box("exemple", "Exemples", [
        ul([
          "**Station spatiale** (r = 6378 + 400 = 6778 km, $G M_{Terre} = 3,986 × 10^{14}$ m³/s²) : v = 7,67 km/s et T = 5553 s = 92,6 min ; vitesse de libération à cette altitude : 10,85 km/s.",
          "**Satellite géostationnaire** (T = 86 164 s, un jour sidéral) : $r = (G M T²/4π²)^{1/3} = 42\u00a0164$ km, soit une altitude de 35 786 km.",
        ]),
      ]),
      p(`Dans le logiciel (simulation 5), la fenêtre **Invariants** donne l’énergie, le moment cinétique et la période pour chaque méthode : ce sont les deux quantités qu’un calcul exact conserverait.`),
    ],
    5: [
      p(`**Orbite exacte.** Avec l’anomalie moyenne $M = n t$ (où $n = 2π/T$) et l’anomalie excentrique $E$, l’**équation de Kepler** $M = E − e sin E$ se résout par la méthode de Newton (4 à 9 itérations pour $e ≤ 0,99$, résidu $4 × 10^{−16}$) ; puis $x = a (cos E − e)$ et $y = a sqrt{1 − e²} sin E$. C’est la référence du logiciel (unités normalisées $G M = a = 1$, $T = 2π$).`),
      p(`**Précession numérique.** Un schéma symplectique d’ordre 2 de pas $h$ résout exactement un hamiltonien légèrement perturbé ; cette perturbation fait précesser l’ellipse, à chaque orbite, de :`),
      eq(`Δω = − frac{π}{8} frac{G M h²}{a³} frac{4 + e²}{(1 − e²)³}`),
      p(`Pour $e = 0,5$ et 200 pas par orbite ($h = T/200$), la formule donne $Δω = −0,2237°$ par orbite ; le logiciel mesure −0,2232° avec Verlet et −0,2234° avec Euler symplectique (écart 0,02 à 0,2 % sur l’ensemble des $e$ de 0,1 à 0,9). Avec 400 pas par orbite, elle donne −0,0559° (mesuré −0,0559°) : en $h²$, comme prévu. RK4, **non symplectique**, ne fait presque pas précesser (+0,00029°, soit 770 fois moins).`),
      tbl(
        ["Schéma (200 pas/orbite, 100 orbites)", "Énergie", "Moment cinétique"],
        [
          ["Euler explicite", "diverge (E > 0, l’astre s’échappe)", "dérive"],
          ["Euler symplectique", "écart relatif borné à 5,2 × 10⁻²", "exact à l’arrondi"],
          ["Verlet", "écart relatif borné à 2,7 × 10⁻³", "exact à l’arrondi"],
          ["RK4", "9,9 × 10⁻⁶ puis 9,3 × 10⁻⁵ (dissipation lente)", "1,4 × 10⁻⁵ (non conservé)"],
        ],
        [3000, 3500, 2798],
        { align: ["l", "l", "l"] }
      ),
      p(`**Problème à N corps.** Avec un adoucissement $ε$ (qui évite une force infinie lors d’une rencontre), l’accélération de la masse $i$ est :`),
      eq(`a_i = G sum_{j≠i} frac{m_j (r_j − r_i)}{(|r_j − r_i|² + ε²)^{3/2}}`),
      p(`La force sur $i$ due à $j$ est l’opposée de celle sur $j$ due à $i$ (troisième loi de Newton) : l’impulsion totale $P = Σ m v$ se conserve (à l’arrondi). Le moment cinétique $L = Σ m r × v$ se conserve exactement avec les schémas symplectiques. Triangle de Lagrange : trois masses égales aux sommets d’un triangle équilatéral de côté $s$ tournent rigidement avec $ω² = 3 G m / s³$, **mais** pour trois masses égales cet équilibre est instable (critère de Routh : stabilité si $(Σ m)² > 27 Σ_{i<j} m_i m_j$, ce qui n’est pas le cas ici).`),
    ],
    6: [
      p(`**Intégrales premières.** Le problème à $N$ corps possède dix intégrales premières classiques (énergie, impulsion à 3 composantes, mouvement du centre de masse à 3 composantes, moment cinétique à 3 composantes) liées par le théorème de Noether aux symétries : translation du temps, translation d’espace, galiléennes, rotation. Bruns (1887) puis Poincaré (1890) ont montré qu’il n’y en a pas d’autres algébriques simples pour $N ≥ 3$ : pas de solution générale.`),
      p(`**Le « huit ».** Trois masses égales parcourent la même courbe, décalées d’un tiers de période. Avec les conditions initiales numériques de Simó ($G = m = 1$) : $E = −1,28714199$, $P = L = 0$ et période $T = 6,32591398$ ; le logiciel retrouve la fermeture à $8 × 10^{−8}$. Une perturbation de $10^{−3}$ y croît **linéairement** (0,058 à $2T$, 0,275 à $10T$) : l’orbite est stable.`),
      p(`**Décider du chaos.** Entre $t = 2$ et $t = 16$ un écart de $10^{−9}$ est multiplié par 7,3 pour le huit (la croissance linéaire donnerait $16/2 = 8$), par $6,6 × 10^{4}$ pour un amas de six corps (graine 42, $λ ≈ 0,8$) et par $5 × 10^{7}$ pour la graine 7. Un ajustement exponentiel sur la croissance linéaire du huit donnerait un faux $λ ≈ 0,12 > 0$ : on décide du chaos sur l’**amplification réelle** (seuil 1000), non sur $λ$.`),
      p(`**Pas adaptatif.** À $e = 0,9$, 200 pas par orbite détruisent tous les schémas à pas fixe (RK4 perd 753 fois l’énergie du cas précédent) car le périastre n’est pas résolu. RK45 résout ce cas en 698 pas pour $2,3 × 10^{−6}$ sur deux orbites, avec un pas 300 fois plus court au périastre qu’à l’apoastre ; RK4 à pas fixe en demande 32 000 pour $7,8 × 10^{−6}$. Pour un amas, l’énergie des schémas d’ordre élevé n’est pas bornée de façon lisse (les rencontres rapprochées la font sauter) : on compare des maxima, on n’exige pas un plateau.`),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Orbite de Kepler : {{physicslab --level 5 --sim 5}}. N corps : {{physicslab --level 5 --sim 6}}.`),
      ul([
        "Sur l’orbite, choisissez une planète du Système solaire avec les boutons (Mercure, Terre, Mars, Halley, Pluton), puis augmentez l’**excentricité** : les pas deviennent insuffisants au périastre, RK45 s’adapte.",
        "Cochez Verlet et RK4 et regardez la colonne **Précession** du tableau : elle correspond à la formule du niveau 5.",
        "Sur les N corps, comparez **Le huit**, **Le triangle de Lagrange** et **Amas au hasard**, puis lisez le verdict « chaotique » ou « stable » sous **Chaos**, fondé sur l’amplification de l’écart du jumeau.",
      ]),
    ]),
    fig("shot_s06.png", "Le « huit » dans la simulation 6, au niveau 5 : les trois étoiles et les traces de chaque méthode, le tableau de conservation à droite, les distances entre les corps en bas.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "La Lune tombe-t-elle ?",
      statement: [
        p(`La Lune tombe-t-elle sur la Terre ? Expliquez pourquoi elle ne s’écrase pas, en vous aidant du canon de Newton.`),
      ],
      solution: [
        p(`Oui, elle tombe **en permanence** : la gravité la tire sans arrêt vers la Terre. Mais elle a aussi une grande vitesse perpendiculaire à cette direction. Pendant qu’elle tombe, la surface de la Terre s’« éloigne » en s’incurvant, et la Lune manque le sol : elle tombe indéfiniment en tournant. Un boulet tiré assez fort fait de même.`),
      ],
    },
    {
      level: 3,
      title: "La troisième loi de Kepler",
      statement: [
        p(`Avec T² = a³ (T en années, a en UA), calculez la période de Mars (a = 1,524 UA), de Jupiter (a = 5,203 UA), et le demi-grand axe d’une comète de période 8 ans.`),
      ],
      solution: [
        p(`Mars : $T = 1,524^{3/2}$ = **1,88 an**. Jupiter : $T = 5,203^{3/2}$ = **11,87 ans**. Comète : $a = T^{2/3} = 8^{2/3} = 4$, donc **4 UA**.`),
      ],
    },
    {
      level: 3,
      title: "Peser à 400 km d’altitude",
      statement: [
        p(`Au sol (distance 6371 km du centre de la Terre) un kilogramme pèse 9,82 N. Combien pèse-t-il à 400 km d’altitude ? Exprimez le résultat en pourcentage du poids au sol et expliquez pourquoi les astronautes flottent.`),
      ],
      solution: [
        p(`La force varie en $1/d²$ : $9,82 × (6371/6771)² = 9,82 × 0,885$ = **8,69 N**, soit **88,5 %** du poids au sol.`),
        p(`Les astronautes flottent parce que la station et eux tombent ensemble à la même accélération : en chute libre il n’y a pas d’appui, donc pas de sensation de poids.`),
      ],
    },
    {
      level: 4,
      title: "La Station spatiale",
      statement: [
        p(`La Station spatiale circule à 400 km d’altitude (r = 6778 km du centre de la Terre ; G M = 3,986 × 10¹⁴ m³/s²). Calculez sa vitesse, sa période et la vitesse de libération à cette altitude.`),
      ],
      solution: [
        p(`$v = sqrt{G M / r} = sqrt{3,986 × 10^{14} / 6,778 × 10^{6}}$ = **7,67 km/s**.`),
        p(`$T = 2π r / v = 2π × 6,778 × 10^{6} / 7669$ = **5553 s**, soit **92,6 minutes** : la station fait environ 15,5 tours par jour.`),
        p(`$v_{lib} = sqrt{2} × 7,67$ = **10,85 km/s**.`),
      ],
    },
    {
      level: 4,
      title: "Le satellite géostationnaire",
      statement: [
        p(`Un satellite géostationnaire reste au-dessus du même point de l’équateur : sa période est un jour sidéral (86 164 s). Calculez le rayon de son orbite et son altitude (rayon terrestre 6378 km).`),
      ],
      solution: [
        p(`Pour une orbite circulaire $v² = G M/r$ et $v = 2π r / T$, d’où $r³ = G M T² / (4π²) = 3,986 × 10^{14} × 86164² / 39,478 = 7,496 × 10^{22}$ m³ et $r$ = **42 164 km**.`),
        p(`Altitude : 42 164 − 6378 = **35 786 km**.`),
      ],
    },
    {
      level: 5,
      title: "Le triangle de Lagrange",
      statement: [
        p(`Trois masses égales $m$ sont aux sommets d’un triangle équilatéral de côté $s$. (a) Trouvez la pulsation de rotation $ω$ autour du centre de masse. (b) Calculez $ω$ et la période pour $G = m = s = 1$. (c) Le triangle est-il stable ?`),
      ],
      solution: [
        p(`(a) Chaque masse est à la distance $s/sqrt{3}$ du centre. La force des deux autres, projetée vers le centre, vaut $2 × (G m² / s²) cos 30° = sqrt{3} G m² / s²$. Elle est égale à $m ω² (s/sqrt{3})$, donc $ω² = 3 G m / s³$.`),
        p(`(b) $ω = sqrt{3}$ = **1,732 rad/s** et $T = 2π/ω$ = **3,628**.`),
        p(`(c) Critère de Routh : stable si $(Σ m)² > 27 Σ_{i<j} m_i m_j$. Ici $(3m)² = 9 m²$ et $27 × 3 m² = 81 m²$ : **instable** (à masses égales ; il devient stable si l’une des masses est très grande devant les autres, comme dans le cas des Troyens de Jupiter).`),
      ],
    },
    {
      level: 5,
      title: "Prédire la précession numérique",
      statement: [
        p(`Pour une orbite d’excentricité e = 0,5 (G M = a = 1, T = 2π), calculez par la formule $Δω = −(π/8)(G M h²/a³)(4 + e²)/(1 − e²)³$ la précession par orbite de Verlet avec 200 puis 400 pas par orbite, et comparez aux mesures −0,2232° et −0,0559°. Combien de pas par orbite faut-il pour que la précession soit inférieure à 0,01° par orbite ?`),
      ],
      solution: [
        p(`$(4 + e²)/(1 − e²)³ = 4,25/0,421875 = 10,074$. À 200 pas : $h = 2π/200 = 0,031416$, $h² = 9,8696 × 10^{−4}$ et $Δω = −(π/8) × 9,8696 × 10^{−4} × 10,074 = −3,904 × 10^{−3}$ rad = **−0,2237°** (mesuré −0,2232°, écart 0,2 %).`),
        p(`À 400 pas, $h²$ est divisé par 4 : **−0,0559°** (mesuré −0,0559°).`),
        p(`Pour 0,01° il faut $Δω ∝ 1/N²$ : $N = 200 × sqrt{0,2237/0,01} = 200 × 4,73$ ≈ **946 pas par orbite**.`),
      ],
    },
    {
      level: 6,
      title: "Énergie du « huit »",
      statement: [
        p(`Le « huit » (G = m = 1) a pour positions initiales $(0,97000436 ; −0,24308753)$, $(−0,97000436 ; 0,24308753)$ et $(0 ; 0)$, et pour vitesses $(0,466203685 ; 0,43236573)$, la même, et $(−0,93240737 ; −0,86473146)$. Calculez l’énergie cinétique, l’énergie potentielle, l’énergie totale et l’impulsion totale.`),
      ],
      solution: [
        p(`Énergie cinétique : $½ (|v_1|² + |v_2|² + |v_3|²) = ½ (0,4043 + 0,4043 + 1,6171)$ = **1,21286**.`),
        p(`Distances : $r_{12} = 2,0000$, $r_{13} = r_{23} = 1,0000$. Énergie potentielle : $−(1/2 + 1 + 1)$ = **−2,5000**. Énergie totale : $1,21286 − 2,5$ = **−1,28714**, la valeur du logiciel.`),
        p(`Impulsion : $0,466203685 + 0,466203685 − 0,93240737 = 0$ ; de même pour la seconde composante : $P = 0$ (le centre de masse reste au repos à l’origine).`),
      ],
    },
    {
      level: 6,
      title: "Chaotique ou non ?",
      statement: [
        p(`Entre t = 2 et t = 16, le logiciel mesure une amplification d’un écart initial de 10⁻⁹ de 7,3 pour le « huit » et de 6,6 × 10⁴ pour un amas de six corps. (a) Que prédit une croissance linéaire du type δ ∝ t ? (b) Calculez l’exposant λ = ln(amplification)/14 dans chaque cas. (c) Quel critère permet de conclure sans ambiguïté ?`),
      ],
      solution: [
        p(`(a) Une croissance linéaire donne $δ(16)/δ(2) = 16/2$ = **8**, très proche des 7,3 mesurés : le huit est stable.`),
        p(`(b) Huit : $ln(7,3)/14 = 0,14$ ; amas : $ln(6,6 × 10^{4})/14$ = **0,79**. Le « 0,14 » n’a aucun sens physique ici : il provient d’un ajustement exponentiel appliqué à une croissance linéaire.`),
        p(`(c) On compare l’**amplification** à un seuil : au-dessus de 1000 le système est déclaré chaotique (amas : $6,6 × 10^{4} ≫ 1000$), en dessous non (huit : 7,3).`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Deux corps : ellipse exacte, **lois de Kepler** ($T² = a³$ en années et UA). Trois corps ou plus : plus de solution générale, calcul pas à pas et chaos.",
      "Le calcul numérique d’une orbite se juge à ce qu’il conserve : **énergie** et **moment cinétique** ; les schémas symplectiques font précesser l’ellipse d’une quantité **prédite** ($−(π/8)(G M h²/a³)(4+e²)/(1−e²)³$ par orbite).",
      "Un système est déclaré chaotique quand un petit écart est amplifié d’un facteur supérieur à 1000 ; le « huit » ($E = −1,28714$) ne l’est pas.",
      "À forte excentricité, un **pas adaptatif** (RK45) est indispensable près du périastre.",
    ]),
  ]),
};
