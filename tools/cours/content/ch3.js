const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 3,
  title: "Le pendule, et le chaos",
  lead: "Un pendule règle les horloges depuis Galilée. Accrochez-en un second au bout du premier et tout change : le mouvement devient imprévisible. Entre les deux se joue l’une des grandes idées du XXᵉ siècle, le chaos déterministe.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** la période d’un pendule, une solution exacte qui n’est pas une formule simple, et ce que veut dire « chaotique » : un mouvement parfaitement déterminé mais impossible à prévoir longtemps.",
      "**Le résultat marquant :** un écart de départ d’un milliardième de radian est multiplié par 4400 en 8 secondes.",
      "**Dans le logiciel :** module M3, simulations 3 (pendule simple, {{--sim 3}}) et 4 (pendule double, {{--sim 4}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`Le balancier d’une horloge va et vient avec une régularité parfaite : si vous savez où il est maintenant, vous savez où il sera dans une heure. Le mouvement est **régulier**.`),
      p(`Accrochez maintenant une seconde tige au bout de la première : vous obtenez un **double pendule**, un jouet fascinant dont le mouvement semble agité sans règle. Pourtant il obéit à des lois parfaitement précises, aucun hasard n’y intervient.`),
      p(`Lâchez deux double-pendules presque identiques, qui ne diffèrent que d’un milliardième de radian au départ. Pendant quelques secondes on ne voit aucune différence. Puis les deux mouvements se séparent, et très vite n’ont plus rien de commun. C’est l’**effet papillon**, popularisé par le météorologue Edward Lorenz : une toute petite cause peut avoir de grands effets, à condition d’attendre assez longtemps.`),
      p(`La leçon est étonnante : prévoir loin est impossible, **non pas parce que l’ordinateur calcule mal, mais parce que le système est ainsi fait**. Mesurer le départ plus précisément ne repousse l’horizon que très peu.`),
    ],
    2: [
      p(`Un pendule simple oscille avec une période qui, pour de **petits angles**, ne dépend pas de l’amplitude (c’est l’isochronisme découvert par Galilée). Pour de plus grands angles, la période s’allonge : elle dépasse de 7 % celle des petits angles à 60°, de 18 % à 90° et de 37 % à 120°. Il existe une solution exacte, mais elle utilise des fonctions mathématiques plus savantes (fonctions elliptiques) : « exacte » ne veut pas dire « formule simple ».`),
      p(`Le **double pendule** n’a pas de solution : seul le calcul pas à pas permet de le suivre. Il est **chaotique**, ce qui a trois conséquences : (1) il est déterministe, le futur est entièrement fixé par le présent ; (2) il est extrêmement sensible aux conditions initiales ; (3) son mouvement n’est jamais périodique.`),
      p(`Pour l’ordinateur, cela fixe un **horizon de prédiction** : au-delà d’une certaine durée, ce qu’il affiche n’a plus de rapport avec la réalité. Cet horizon dépend de la précision du départ, mais aussi de la qualité de la méthode de calcul. Dans le logiciel, avec un pas d’un soixantième de seconde, on suit la vraie trajectoire jusqu’à 0,43 s avec la méthode d’Euler, 0,52 s avec Euler symplectique, 3 s avec Verlet et 8,6 s avec RK4 (premier instant où l’écart à la référence dépasse 0,1).`),
      fig("fig_chaos.png", "Écart entre deux pendules doubles qui diffèrent de 10⁻⁹ rad au départ (calcul très précis, RK45). L’échelle verticale est logarithmique : la droite est une croissance exponentielle.", 0.82),
    ],
    3: [
      p(`La période d’un pendule de longueur $L$ (en m) pour de petits angles, avec $g = 9,81$ m/s², vaut :`),
      eq(`T = 2π sqrt{frac{L}{g}}`),
      p(`Elle ne dépend **ni de la masse, ni de l’amplitude** (si elle est petite). Voici quelques valeurs :`),
      tbl(
        ["Longueur L", "Période T", "Remarque"],
        [
          ["0,25 m", "1,00 s", "petit pendule de bureau"],
          ["0,994 m", "2,00 s", "« pendule qui bat la seconde » (une seconde par demi-oscillation)"],
          ["4 m", "4,01 s", "quatre fois plus long, deux fois plus lent"],
        ],
        [1900, 1900, 5498],
        { align: ["c", "c", "l"] }
      ),
      box("exemple", "Exemple : sur la Lune", [
        p(`La pesanteur lunaire est $g = 1,62$ m/s². Un pendule de 1 m y bat avec une période $T = 2π × √(1 / 1,62)$ = **4,94 s**, contre 2,01 s sur Terre : à longueur égale, il va 2,5 fois plus lentement.`),
      ]),
      p(`Quand l’amplitude est grande, la période s’allonge. Le pendule de 1 m lâché de 60° met 2,153 s au lieu de 2,006 s : une horloge réglée pour de petites oscillations mais qui balance à 60° prendrait **7 % de retard**, soit plus d’une heure et demie par jour (98 minutes).`),
    ],
    4: [
      p(`Le poids et la tension du fil donnent, pour l’angle $θ$ par rapport à la verticale, l’équation :`),
      eq(`θ'' = − frac{g}{L} sin θ`),
      p(`Pour de petits angles, $sin θ ≈ θ$ et l’on retrouve l’oscillateur du chapitre 2, avec $ω_0 = sqrt{g/L}$. L’énergie vaut $E = ½ m L² θ'² + m g L (1 − cos θ)$ et se conserve. Lâché sans vitesse d’un angle $θ_0$, le pendule atteint en bas la vitesse :`),
      eq(`v = sqrt{2 g L (1 − cos θ_0)}`),
      p(`et la tension du fil y vaut $m g + m v²/L = m g (3 − 2 cos θ_0)$.`),
      box("exemple", "Exemple", [
        p(`Pendule de longueur L = 1,5 m, de masse 0,2 kg, lâché de 60° : $v = √(2 × 9,81 × 1,5 × 0,5)$ = **3,84 m/s** en bas, et la tension vaut $0,2 × 9,81 × (3 − 2 × 0,5)$ = **3,92 N**, soit deux fois le poids.`),
      ]),
      p(`Dans le logiciel, le niveau 4 affiche pour le pendule simple la solution exacte (période réelle) à côté de la méthode simple (Euler) et de la méthode précise (RK4) : regardez Euler s’écarter de plus en plus.`),
    ],
    5: [
      p(`**Solution exacte.** L’énergie donne $θ'^2 = (2g/L)(cos θ − cos θ_0)$ ; en séparant les variables on obtient une intégrale elliptique, et la période exacte vaut :`),
      eq(`T = frac{4}{ω_0} K(k), k = sin frac{θ_0}{2}, K(k) = int_{0}^{π/2} frac{dφ}{sqrt{1 − k² sin² φ}}`),
      p(`$K$ se calcule très vite par la moyenne arithmético-géométrique : $K(k) = π / (2 AGM(1, sqrt{1 − k²}))$. Le rapport $T/T_0$ vaut 1,0019 à 10°, 1,0174 à 30°, 1,0732 à 60°, 1,1803 à 90°, 1,3729 à 120° et 2,4394 à 170°. Le développement $T/T_0 = 1 + θ_0²/16 + 11 θ_0⁴/3072 + …$ donne 1,07285 à 60° (erreur de 3 × 10⁻⁴) mais seulement 1,343 à 120° (erreur de 2 %).`),
      p(`**Pendule double.** Avec $m_1, m_2$, $l_1, l_2$ et les angles $θ_1, θ_2$, le lagrangien $ℒ = T − V$ est :`),
      eq(`T = frac{1}{2}(m_1 + m_2) l_1² θ_1'² + frac{1}{2} m_2 l_2² θ_2'² + m_2 l_1 l_2 θ_1' θ_2' cos(θ_1 − θ_2)`),
      eq(`V = −(m_1 + m_2) g l_1 cos θ_1 − m_2 g l_2 cos θ_2`),
      p(`Les équations d’Euler-Lagrange forment un système de deux équations du second ordre couplées, que le logiciel met sous forme de quatre équations du premier ordre $(θ_1, θ_2, ω_1, ω_2)$.`),
      p(`**Chaos.** Deux trajectoires voisines s’écartent comme :`),
      eq(`δ(t) = δ_0 e^{λ t}`),
      p(`où $λ$ est l’**exposant de Lyapunov**. Dans le logiciel on suit deux copies du système décalées de $10^{−9}$ rad et on ajuste $ln δ(t)$ par une droite sur la partie exponentielle (avant la saturation). On mesure $λ ≈ 1$ par seconde (0,98 dans la figure du niveau 2, 1,1 dans le logiciel avec d’autres seuils d’ajustement). L’écart étant multiplié par 4400 en 8 s, $λ = ln(4400)/8 = 1,05$ s⁻¹. **Horizon** :`),
      eq(`t_H = frac{1}{λ} ln frac{δ_{max}}{δ_0}`),
    ],
    6: [
      p(`**Pourquoi le pendule simple ne peut pas être chaotique.** Son espace des phases a deux dimensions $(θ, θ')$. Pour un système autonome continu, le théorème de Poincaré-Bendixson interdit le chaos en dimension 2 : les trajectoires ne peuvent pas se croiser, donc ne peuvent pas se replier sur elles-mêmes. Le pendule double vit dans un espace à quatre dimensions ; l’énergie conservée le réduit à une hypersurface de dimension 3, ce qui suffit au chaos (rupture des tores invariants du théorème KAM à grande énergie).`),
      p(`**Définition.** L’exposant de Lyapunov maximal est la limite de $(1/t) ln(|δ(t)|/|δ_0|)$ quand le temps $t$ devient grand et l’écart initial $δ_0$ infiniment petit. Numériquement on ne dispose que d’une fenêtre : il faut $δ$ au-dessus de l’arrondi ($δ ≥ 20 δ_0$ dans le logiciel) et en dessous de la saturation ($δ ≤ 0,1$) ; hors de la fenêtre l’ajustement donne un faux $λ$. Un faux $λ$ peut aussi apparaître sur une croissance **linéaire** (système régulier mais instable), d’où la règle du chapitre 4 : on décide du chaos sur l’amplification réelle de l’écart.`),
      p(`**Ce que vaut une trajectoire numérique chaotique.** Une trajectoire calculée avec une erreur $ε$ par unité de temps ne reste proche de la vraie trajectoire de même départ que pendant $t_H ≈ ln(δ_{max}/ε)/λ$ : au-delà, elle est sans rapport avec elle. Elle peut rester proche d’une autre vraie trajectoire (lemme de pistage, rigoureux pour les systèmes hyperboliques), ce qui justifie les statistiques mais pas les prévisions à long terme. Pour RK4 à $Δt = 1/60$ s, l’horizon mesuré (8,6 s) correspond à une erreur effective $δ_{max} e^{−λ t_H} ≈ 10^{−5}$, à comparer aux $10^{−9}$ de l’écart physique de la figure du niveau 2.`),
      p(`**Conséquence pour le calcul en précision réduite.** L’arrondi d’un nombre en simple précision (6 × 10⁻⁸) fixe $δ_0$ : avec $λ = 1$ s⁻¹ et $δ_{max} = 0,1$ l’horizon est de 14 s en simple précision contre 34 s en double précision (1,1 × 10⁻¹⁶). Gagner 8 chiffres significatifs ne rallonge l’horizon que de $ln(10⁸)/λ ≈ 18$ s : cette idée est au cœur du chapitre 7.`),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Pendule simple : {{physicslab --level 4 --sim 3}}. Pendule double : {{physicslab --level 5 --sim 4}}.`),
      ul([
        "Au pendule simple, tirez le curseur de l’**angle initial** de 10° à 170° et lisez la courbe « période selon l’amplitude » : elle monte de 1 à 2,4.",
        "Au pendule double, observez le bouton **Départ au hasard** et le graphique **Sensibilité** : l’écart entre la référence et son jumeau (décalé de 10⁻⁹ rad) croît en ligne droite sur une échelle logarithmique. Le niveau 5 affiche l’exposant λ et l’amplification.",
        "Cochez les cinq méthodes et notez, dans le tableau, l’instant où chacune « décroche » de la référence (colonne t_H) : Euler 0,43 s, RK4 8,6 s au pas de 1/60 s.",
      ]),
    ]),
    fig("shot_s04.png", "Le pendule double au niveau 5 : les trajectoires en couleur (une par méthode), le tableau des horizons de prédiction à droite, la sensibilité aux conditions initiales en bas.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "Qui est responsable ?",
      statement: [
        p(`Un ordinateur prédit très bien la position d’un balancier d’horloge dans une heure, mais pas celle d’un double pendule dans une minute. Est-ce parce que l’ordinateur calcule moins bien dans le second cas ? Expliquez avec vos mots.`),
      ],
      solution: [
        p(`Non. L’ordinateur calcule aussi bien dans les deux cas. La différence vient du **système** : le balancier est régulier (une petite erreur reste petite) alors que le double pendule est chaotique (une petite erreur grossit exponentiellement). Même avec un ordinateur parfait, un départ connu avec une précision finie donne une prévision qui cesse d’être fiable après un certain temps.`),
      ],
    },
    {
      level: 2,
      title: "Vrai ou faux ?",
      statement: [
        p(`(a) Un pendule d’un mètre met environ deux secondes pour un aller-retour sur Terre. (b) Si l’on double la masse du pendule, sa période double. (c) Un pendule lâché de très haut a une période plus longue que lâché de petit angle. (d) Le mouvement d’un double pendule est aléatoire.`),
      ],
      solution: [
        p(`(a) **Vrai** : T = 2π √(1/9,81) = 2,006 s. (b) **Faux** : la période ne dépend pas de la masse. (c) **Vrai** : à 60° elle est plus longue de 7 %, à 120° de 37 %. (d) **Faux** : il est déterministe ; il est imprévisible à long terme, ce qui est tout autre chose que le hasard.`),
      ],
    },
    {
      level: 3,
      title: "Périodes",
      statement: [
        p(`(a) Quelle est la période d’un pendule de longueur 0,5 m ? (b) Quelle longueur donne une période de 2 s ? (c) Quelle est la période d’un pendule de 1 m sur la Lune (g = 1,62 m/s²) ?`),
      ],
      solution: [
        p(`(a) T = 2π √(0,5/9,81) = **1,42 s**. (b) $L = g T²/(4π²) = 9,81 × 4 / 39,48$ = **0,994 m**. (c) T = 2π √(1/1,62) = **4,94 s**.`),
      ],
    },
    {
      level: 3,
      title: "Une horloge en retard",
      statement: [
        p(`Une horloge à pendule de 1 m est réglée pour de petites oscillations (T₀ = 2,006 s). Par erreur on la lance avec une amplitude de 60°, pour laquelle la période est 1,0732 fois plus longue. Quelle est la nouvelle période ? De combien l’horloge retarde-t-elle en 24 heures ?`),
      ],
      solution: [
        p(`Nouvelle période : 2,006 × 1,0732 = **2,153 s**. L’horloge fait moins d’oscillations que prévu, dans le rapport 1/1,0732 = 0,9318 : elle retarde de $1 − 0,9318 = 6,8 %$ du temps, soit $24 × 60 × 0,068$ = **98 minutes** par jour.`),
      ],
    },
    {
      level: 4,
      title: "Vitesse et tension en bas",
      statement: [
        p(`Un pendule de longueur L = 1,5 m et de masse 0,2 kg est lâché sans vitesse de 60° (g = 9,81 m/s²). Calculez sa vitesse et la tension du fil au passage par la verticale.`),
      ],
      solution: [
        p(`Conservation de l’énergie : $m g L (1 − cos 60°) = ½ m v²$, donc $v = sqrt{2 g L (1 − cos 60°)} = sqrt{2 × 9,81 × 1,5 × 0,5}$ = **3,84 m/s**.`),
        p(`En bas, la seconde loi de Newton donne $T − m g = m v²/L$, d’où $T = m g (3 − 2 cos 60°) = 0,2 × 9,81 × 2$ = **3,92 N**, soit deux fois le poids.`),
      ],
    },
    {
      level: 5,
      title: "Série et solution exacte",
      statement: [
        p(`On approche $T/T_0$ par $1 + θ_0²/16 + 11 θ_0⁴/3072$. Calculez cette approximation à 60° et à 120° et comparez aux valeurs exactes 1,07318 et 1,37288. Commentez.`),
      ],
      solution: [
        p(`À 60° ($θ_0 = 1,0472$ rad) : $1 + 0,06854 + 0,00431$ = **1,07285**, erreur relative 3 × 10⁻⁴. À 120° ($θ_0 = 2,0944$ rad) : $1 + 0,27416 + 0,06889$ = **1,34305**, erreur relative **2,2 %**.`),
        p(`La série est excellente pour des angles modérés mais converge mal près de 180° (sa variable naturelle est $k = sin(θ_0/2)$, qui tend vers 1 où $K$ diverge) : il faut alors la moyenne arithmético-géométrique, comme le fait le logiciel.`),
      ],
    },
    {
      level: 5,
      title: "Temps de prévision",
      statement: [
        p(`Dans le logiciel un écart de 10⁻⁹ rad est amplifié d’un facteur 4400 en 8 s. (a) Calculez l’exposant λ. (b) Quelle durée de prévision obtient-on pour un écart initial de 10⁻⁹, de 10⁻¹², de 10⁻¹⁵, en acceptant une erreur finale de 0,1 ? (c) De combien gagne-t-on en améliorant la mesure initiale d’un facteur 1000 ?`),
      ],
      solution: [
        p(`(a) $λ = ln(4400)/8 = 8,39/8$ = **1,05 s⁻¹**.`),
        p(`(b) $t_H = ln(0,1/δ_0)/λ$ : **17,6 s** pour 10⁻⁹ ; **24,2 s** pour 10⁻¹² ; **30,7 s** pour 10⁻¹⁵.`),
        p(`(c) Gagner un facteur 1000 ajoute seulement $ln(1000)/λ$ = **6,6 s** : la prévision n’augmente que de façon **logarithmique** avec la précision. Passer de 10⁻⁹ à 10⁻¹⁵ (un million de fois plus précis) ne rallonge l’horizon que de 13 s.`),
      ],
    },
    {
      level: 6,
      title: "Dimensions et chaos",
      statement: [
        p(`(a) Quelle est la dimension de l’espace des phases du pendule simple ? du pendule double ? (b) Pourquoi le pendule simple ne peut-il pas être chaotique ? (c) Quelle est la dimension de la surface d’énergie constante du pendule double ?`),
      ],
      solution: [
        p(`(a) 2 pour le pendule simple $(θ, θ')$ ; 4 pour le double $(θ_1, θ_2, θ_1', θ_2')$.`),
        p(`(b) Pour un système autonome continu en dimension 2, les trajectoires d’un plan ne peuvent pas se croiser : le théorème de Poincaré-Bendixson montre alors que le comportement asymptotique est simple (point fixe, cycle limite), jamais chaotique. Il faut au moins trois dimensions.`),
        p(`(c) L’énergie est une intégrale première : le mouvement a lieu sur une hypersurface de dimension $4 − 1$ = **3**, ce qui suffit à la dynamique chaotique.`),
      ],
    },
    {
      level: 6,
      title: "Simple ou double précision ?",
      statement: [
        p(`On suppose λ = 1 s⁻¹ et on accepte une erreur finale de 0,1. L’arrondi initial est de 6 × 10⁻⁸ en simple précision et de 1,1 × 10⁻¹⁶ en double précision. Calculez l’horizon dans chaque cas. De combien le double rallonge-t-il l’horizon, et que dire du coût de cette précision supplémentaire ?`),
      ],
      solution: [
        p(`$t_H = ln(0,1/ε)/λ$ : simple précision **14,3 s** ; double précision **34,4 s**.`),
        p(`Gain : 20 s seulement, alors que l’on a gagné **huit chiffres** significatifs. Le chaos transforme un énorme gain de précision en petit gain de temps : doubler la précision n’est pas un investissement rentable pour prévoir plus longtemps, d’où l’intérêt de comparer des **accélérations à état fixé** plutôt que des trajectoires longues (chapitre 7).`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Le pendule simple a une période $T = 2π sqrt{L/g}$ quasi indépendante de l’amplitude ; sa solution exacte existe mais fait intervenir une fonction elliptique (+7 % à 60°, +37 % à 120°).",
      "Un système **chaotique** est déterministe mais amplifie les écarts : $δ(t) = δ_0 e^{λ t}$ ; il lui faut au moins trois dimensions (le double pendule en a quatre).",
      "L’**horizon de prédiction** $t_H = ln(δ_{max}/δ_0)/λ$ ne croît que comme le **logarithme** de la précision de départ.",
      "À pas égal, la qualité de la méthode décide de la durée pendant laquelle on peut croire la simulation : 0,43 s (Euler) à 8,6 s (RK4) pour le double pendule.",
    ]),
  ]),
};
