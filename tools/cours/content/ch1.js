const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl, sp } = B;

module.exports = {
  num: 1,
  title: "Calculer un mouvement pas à pas",
  lead: "Un ordinateur ne « sait » pas où sera une balle dans une seconde. Il le découvre petit à petit, en avançant par tout petits pas de temps. Tout le reste du cours repose sur cette idée.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** pourquoi une simulation avance « pas à pas », pourquoi toutes les méthodes ne se valent pas, et comment mesurer la qualité d’une méthode (erreur, ordre).",
      "**Le cas d’étude :** une balle lancée, freinée par l’air, dont la solution exacte est connue.",
      "**Dans le logiciel :** modules M0 et M1, simulation 1 ({{--sim 1}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`Imaginez que vous marchez dans le brouillard, avec seulement une boussole et un compteur de pas. Toutes les dix secondes, vous regardez la boussole et vous notez la direction. Puis vous marchez droit devant vous pendant dix secondes, sans rien regarder. Si le sentier tourne, vous vous écartez un peu. Vous recommencez : regarder, avancer, regarder, avancer.`),
      p(`Un ordinateur calcule un mouvement exactement comme cela. « Regarder », c’est calculer la force qui s’exerce sur l’objet à cet instant. « Avancer », c’est déplacer l’objet un tout petit peu, pendant un court instant que l’on appelle le **pas de temps**. Puis on recommence, des milliers de fois par seconde.`),
      p(`Si l’on regarde plus souvent (des pas plus courts), on suit mieux le sentier, mais il faut faire beaucoup plus de calculs. Si l’on regarde moins souvent, c’est plus rapide, mais on s’écarte du vrai chemin. **Il n’existe pas de calcul parfait : il existe des calculs dont on sait de combien ils se trompent.**`),
      p(`Il y a aussi des méthodes plus malignes que d’autres. La plus simple, due à Léonhard Euler (XVIIIᵉ siècle), ne regarde qu’au début de chaque pas. Une méthode plus maligne regarde au début, au milieu et à la fin du pas pour deviner où va le sentier : avec le même nombre de pas, elle se trompe beaucoup moins. Le logiciel en met cinq en concurrence sur le même problème, et c’est ce qui rend le résultat si instructif.`),
    ],
    2: [
      p(`Pour simuler une balle lancée, on construit un **modèle** : la balle est réduite à un point, la pesanteur la tire vers le bas, l’air la freine un peu. Les lois de Newton fournissent alors un mode d’emploi : connaissant la position et la vitesse à un instant, on sait comment elles vont changer l’instant d’après.`),
      p(`Dans de rares cas on peut écrire la réponse sous forme de formule : c’est la **solution exacte**. Notre balle freinée par un frottement proportionnel à sa vitesse en a une : elle joue le rôle de juge de paix, car on peut y comparer chaque méthode numérique. Dans la plupart des situations réelles il n’y a pas de formule, et le calcul pas à pas est le seul moyen.`),
      p(`Ce qui distingue les méthodes, c’est la façon dont l’erreur diminue quand on raccourcit le pas. On parle de l’**ordre** d’une méthode : avec une méthode d’ordre 1 (Euler), diviser le pas par 2 divise l’erreur par 2 ; avec l’ordre 2 (Verlet), par 4 ; avec l’ordre 4 (RK4), par 16. En passant de 100 à 200 pas, RK4 gagne un facteur 16 sur l’erreur, alors qu’Euler, pour de petits pas, gagne seulement un facteur 2.`),
      p(`Il y a enfin le **coût**. RK4 évalue la force quatre fois par pas, Euler une seule fois. Mais RK4 atteint la même précision avec des pas bien plus longs : il est largement gagnant. La méthode RK45 va plus loin encore : elle choisit elle-même son pas, long quand tout est calme, court quand le mouvement devient délicat.`),
      fig("fig_convergence.png", "Erreur finale de chaque méthode en fonction du nombre de pas, relevée avec le logiciel sur un ressort (chapitre 2). Plus la droite est inclinée, plus la méthode est « maligne » ; les pentes mesurées sont proches de −1, −1, −2 et −4.", 0.78),
    ],
    3: [
      p(`On lâche une balle sans vitesse. Sa vitesse augmente d’environ g = 10 m/s chaque seconde (on prend 10 au lieu de 9,81 pour simplifier). La méthode d’Euler consiste à faire deux petits calculs à chaque pas :`),
      ul([
        "nouvelle vitesse = ancienne vitesse + g × pas ;",
        "nouvelle distance de chute = ancienne distance + ancienne vitesse × pas.",
      ]),
      p(`Prenons un pas de 0,5 s. Le tableau compare ce que calcule l’ordinateur à la vraie réponse (vitesse = g × t, distance = ½ × g × t²) :`),
      tbl(
        ["Temps (s)", "Vitesse Euler (m/s)", "Vitesse exacte (m/s)", "Distance Euler (m)", "Distance exacte (m)"],
        [
          ["0,5", "5", "5", "0", "1,25"],
          ["1,0", "10", "10", "2,5", "5"],
          ["1,5", "15", "15", "7,5", "11,25"],
          ["2,0", "20", "20", "15", "20"],
        ],
        [1500, 1900, 1900, 1900, 2098],
        { align: ["c", "c", "c", "c", "c"] }
      ),
      p(`La vitesse est juste, mais la distance est toujours trop petite : l’ordinateur utilise la vitesse du **début** de chaque pas, alors que la balle accélère pendant le pas. Après 2 s, il manque 5 m sur 20 m. Avec un pas de 0,1 s il ne manquerait plus que 1 m : **l’erreur est proportionnelle au pas**.`),
      box("exemple", "Exemple : le parachutiste", [
        p(`L’air exerce une force de frottement F = k × v (k en kg/s). Pour un parachutiste de masse m = 80 kg et k = 20 kg/s, la vitesse augmente jusqu’à ce que le frottement équilibre le poids : k × v = m × g, soit une **vitesse limite** $v_{∞}$ = m × g / k = 80 × 10 / 20 = **40 m/s** (144 km/h).`),
        p(`Le temps caractéristique est τ = m / k = 4 s. Au bout de τ, il atteint environ 63 % de sa vitesse limite (25 m/s) ; au bout de 3τ = 12 s, environ 95 % (38 m/s).`),
      ]),
    ],
    4: [
      p(`Pour une balle de masse $m$ soumise à son poids et à un frottement linéaire, la deuxième loi de Newton, projetée sur la verticale orientée vers le bas, s’écrit :`),
      eq(`m frac{dv}{dt} = m g − k v`),
      p(`En posant $v_{∞} = m g / k$ et $τ = m / k$, l’équation devient $v' = (v_{∞} − v) / τ$. Avec $v(0) = 0$, la solution est :`),
      eq(`v(t) = v_{∞} paren{1 − e^{−t/τ}}`),
      eq(`x(t) = v_{∞} paren{t − τ paren{1 − e^{−t/τ}}}`),
      p(`C’est la **solution exacte** qui sert de référence dans le logiciel (simulation 1). Avec un frottement d’une autre forme elle n’existe plus et seul le calcul numérique reste.`),
      p(`**Pourquoi Euler « triche » en chute libre.** Sans frottement, l’énergie mécanique $E = ½ m v² − m g x$ ($x$ est la distance de chute) doit rester nulle. Avec Euler, après $n$ pas de durée $Δt$, on a $v = g n Δt$ et $x = g Δt² n(n−1)/2$. On en déduit :`),
      eq(`ΔE = frac{1}{2} m g² Δt² n = frac{1}{2} m g² t Δt`),
      p(`L’énergie fabriquée par la méthode croît avec le temps et reste proportionnelle au pas : c’est la signature d’une méthode d’ordre 1 (exercice 1.5).`),
      box("attention", "Attention : sans frottement, rien à mesurer", [
        p(`Sans frottement, la chute libre est un polynôme du second degré en $t$ : les méthodes de Verlet et RK4 l’intègrent **exactement** et ne se distinguent plus. C’est pour cela que le logiciel ajoute toujours un frottement quand il veut comparer les méthodes.`),
      ]),
      p(`Distinguez enfin l’**erreur de méthode** (l’ordinateur calcule mal le bon modèle) de l’**erreur de modèle** (le modèle lui-même oublie une force). Le calcul numérique ne corrige que la première.`),
    ],
    5: [
      p(`Le problème général est une équation différentielle ordinaire (EDO) $y' = f(t, y)$ avec $y(t_0) = y_0$, où $y$ regroupe positions et vitesses. La méthode d’Euler explicite s’écrit :`),
      eq(`y_{n+1} = y_n + Δt f(t_n, y_n)`),
      p(`Comparons-la au développement de Taylor $y(t+Δt) = y + Δt y' + ½ Δt² y'' + …$ : Euler en garde les deux premiers termes, donc l’**erreur locale** est $O(Δt²)$ à chaque pas. Après $T/Δt$ pas elle s’accumule en une **erreur globale** en $O(Δt)$ : la méthode est d’**ordre 1**. Une méthode est d’ordre $p$ si son erreur globale est en $O(Δt^p)$.`),
      p(`La méthode de Runge-Kutta classique (RK4) combine quatre évaluations de la pente par pas :`),
      eq(`k_1 = f(t_n, y_n)`),
      eq(`k_2 = f(t_n + frac{Δt}{2}, y_n + frac{Δt}{2} k_1)`),
      eq(`k_3 = f(t_n + frac{Δt}{2}, y_n + frac{Δt}{2} k_2)`),
      eq(`k_4 = f(t_n + Δt, y_n + Δt k_3)`),
      eq(`y_{n+1} = y_n + frac{Δt}{6} (k_1 + 2 k_2 + 2 k_3 + k_4)`),
      p(`Elle coïncide avec le développement de Taylor jusqu’au terme en $Δt⁴$ : ordre 4. La méthode de Verlet des vitesses et celle d’Euler symplectique, d’ordres 2 et 1, ont une propriété que RK4 n’a pas : elles conservent une structure géométrique du mouvement (chapitres 2 et 4).`),
      p(`**Stabilité.** Sur l’équation test $y' = −y/τ$, Euler donne $y_{n+1} = (1 − Δt/τ) y_n$. Le facteur est de module inférieur à 1 si et seulement si $Δt < 2τ$ : au-delà, la solution numérique oscille avec une amplitude croissante alors que la vraie solution décroît. Pour notre balle freinée la solution numérique d’Euler est exactement $v_n = v_{∞}(1 − (1 − Δt/τ)^n)$, qui tend bien vers $v_{∞}(1 − e^{−t/τ})$ quand $Δt → 0$.`),
      tbl(
        ["Méthode", "Évaluations de la force par pas", "Ordre théorique", "Pente mesurée"],
        [
          ["Euler explicite", "1", "1", "1,04"],
          ["Euler symplectique", "1", "1", "1,01"],
          ["Verlet des vitesses", "1 (2 sans mémoire)", "2", "2,00"],
          ["Runge-Kutta 4 (RK4)", "4", "4", "4,06"],
        ],
        [3000, 2700, 1800, 1798],
        { align: ["l", "c", "c", "c"], size: 20 }
      ),
      p(`Les pentes mesurées sont celles de la figure du niveau 2, relevées sur un oscillateur avec le logiciel (pentes sur les quatre derniers points).`),
      box("attention", "Attention : où mesurer l’ordre ?", [
        p(`Mesurer l’erreur à un instant multiple de la période ($t = nT$) donne des ordres **faussement élevés** : l’erreur de phase d’ordre 1 s’annule (Verlet paraît d’ordre 4). Le logiciel mesure à $t = 2,7 T$ et dans l’espace des phases (position et vitesse ensemble).`),
      ]),
    ],
    6: [
      p(`**Contrôle de pas.** RK45 est la paire emboîtée de Dormand et Prince, d’ordres 5 et 4, à sept évaluations dont la dernière est réutilisée au pas suivant (propriété FSAL) : six évaluations nouvelles par pas accepté. La différence entre les deux approximations estime l’erreur locale $ε$, normalisée par la tolérance combinée (absolue + relative × |y|). Le pas suivant est proposé par :`),
      eq(`Δt_{nouveau} = 0,9 Δt paren{frac{1}{ε}}^{1/5}`),
      p(`L’exposant $1/5$ vient de l’ordre 4 de l’estimateur d’erreur, le facteur 0,9 est une marge de sécurité, et le nouveau pas est limité entre 0,2 et 5 fois l’ancien pour éviter les changements brutaux. Un pas est refusé et recommencé si $ε > 1$. Avec une tolérance de $10^{−10}$ le logiciel mesure une erreur de $9 × 10^{−11}$ sur la balle freinée.`),
      p(`**Mesurer proprement.** L’erreur est une distance dans l’espace des phases, $sqrt{Δx² + (Δv/ω_0)²}$ pour un oscillateur, afin que position et vitesse soient commensurables, évaluée à $t = 2,7 T$ (pas un multiple de $T$). Un schéma peut être **pré-asymptotique** : RK4 sur un grand angle n’atteint un rapport d’erreur de 16 qu’à partir de 1280 pas, Euler un rapport 2 qu’à partir de 1600 pas ; on ajuste le point de mesure, jamais le seuil du test.`),
      p(`**Ordre contre coût.** Ce qui compte est l’erreur atteinte pour un nombre donné d’évaluations de la force. Mesures sur l’oscillateur de la figure du niveau 2 :`),
      tbl(
        ["Méthode", "Pas", "Évaluations", "Erreur finale"],
        [
          ["Euler", "6400", "6400", "2,3 × 10⁻²"],
          ["Verlet (avec mémoire)", "6400", "6400", "4,8 × 10⁻⁶"],
          ["RK4", "400", "1600", "4,6 × 10⁻⁷"],
        ],
        [3000, 1700, 2000, 2598],
        { align: ["l", "c", "c", "c"] }
      ),
      p(`RK4, avec quatre fois moins d’évaluations que Verlet, est dix fois plus précis ; face à Euler le rapport dépasse $10^{4}$. En revanche RK4 n’est pas symplectique : l’énergie d’un système hamiltonien y dérive lentement (en $(ω_0 Δt)^6$ pour l’oscillateur), alors qu’un schéma symplectique conserve exactement un hamiltonien voisin $text{H̃} = H + O(Δt^p)$ et reste proche de l’énergie vraie sur de très longues durées (analyse par erreur inverse, voir Hairer, Lubich et Wanner en annexe E).`),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Lancez {{physicslab --level 3 --sim 1}} : un projectile avec frottement, comparé à sa solution exacte.`),
      ul([
        "Changez le curseur **Pas de calcul dt** : l’erreur de la « méthode simple » (Euler) diminue proportionnellement au pas, celle de la « méthode précise » (RK4) bien plus vite.",
        "Au niveau 5, cochez les cinq méthodes : chaque couleur est fixe dans tout le logiciel (Euler orange, Euler symplectique citron, Verlet violet, RK4 vert, RK45 rose, solution exacte bleue).",
        "Augmentez le frottement jusqu’à 0 : la méthode précise devient exacte, car la chute libre est un polynôme du second degré (voir l’encadré du niveau 4).",
      ]),
    ]),
    fig("shot_s01.png", "La simulation 1 au niveau 3 : à gauche les paramètres, au centre la trajectoire, à droite le tableau des erreurs, en bas les courbes de hauteur et d’erreur.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "Plus petit ou plus grand ?",
      statement: [p(`Une balle est suivie par ordinateur une première fois avec un pas de 1 seconde, une seconde fois avec un pas de 0,1 seconde. Laquelle des deux simulations est la plus proche de la réalité ? Quel prix paie-t-on ?`)],
      solution: [p(`La seconde : en raccourcissant le pas, on « regarde » plus souvent et l’on suit mieux le mouvement. Le prix est le nombre de calculs : pour simuler la même durée il en faut **dix fois plus** (0,1 s au lieu de 1 s par pas).`)],
    },
    {
      level: 2,
      title: "Les ordres d’une méthode",
      statement: [
        p(`Une méthode d’ordre 2 commet une erreur de 8 cm avec un pas de 0,2 s. Quelle erreur attendre avec un pas de 0,1 s ? avec 0,05 s ? Même question pour une méthode d’ordre 4 qui commet 16 cm avec le pas de 0,2 s.`),
      ],
      solution: [
        p(`Ordre 2 : diviser le pas par 2 divise l’erreur par $2² = 4$ : 8 cm → **2 cm** (pas 0,1 s) → **0,5 cm** (pas 0,05 s).`),
        p(`Ordre 4 : division par $2⁴ = 16$ à chaque fois : 16 cm → **1 cm** → **0,0625 cm**, soit 0,6 mm.`),
      ],
    },
    {
      level: 3,
      title: "Euler à la main",
      statement: [
        p(`Une balle tombe sans frottement, avec g = 10 m/s². Calculez avec la méthode d’Euler, en deux pas de 1 s, la vitesse et la distance de chute à t = 1 s et à t = 2 s. Comparez à la distance exacte à t = 2 s. Que remarquez-vous par rapport au tableau du cours (pas de 0,5 s) ?`),
      ],
      solution: [
        p(`Pas 1 (t = 1 s) : v = 0 + 10 × 1 = **10 m/s** ; x = 0 + 0 × 1 = **0 m**.`),
        p(`Pas 2 (t = 2 s) : v = 10 + 10 × 1 = **20 m/s** ; x = 0 + 10 × 1 = **10 m**.`),
        p(`Distance exacte : ½ × 10 × 2² = 20 m. L’erreur est de **10 m**, le double de celle du tableau (5 m avec un pas deux fois plus court) : l’erreur d’Euler est proportionnelle au pas.`),
      ],
    },
    {
      level: 3,
      title: "Le parachutiste (bis)",
      statement: [
        p(`Un parachutiste de masse m = 70 kg est freiné par une force F = k × v avec k = 14 kg/s (g = 10 m/s²). Calculez sa vitesse limite $v_{∞}$ = m g / k, le temps caractéristique τ = m / k, puis sa vitesse après 5 s et après 10 s sachant que v(t) = $v_{∞}$ × (1 − e^(−t/τ)).`),
      ],
      solution: [
        p(`$v_{∞}$ = 70 × 10 / 14 = **50 m/s** ; τ = 70 / 14 = **5 s**.`),
        p(`À t = 5 s = τ : v = 50 × (1 − e⁻¹) = 50 × 0,632 = **31,6 m/s**.`),
        p(`À t = 10 s = 2τ : v = 50 × (1 − e⁻²) = 50 × 0,865 = **43,2 m/s**.`),
      ],
    },
    {
      level: 4,
      title: "L’énergie fabriquée par Euler",
      statement: [
        p(`Une balle de masse m = 2 kg tombe sans frottement (g = 10 m/s²). On l’intègre par Euler avec un pas Δt = 0,5 s. (a) Calculez l’énergie mécanique E = ½ m v² − m g x à t = 2 s à partir des valeurs v = 20 m/s et x = 15 m du tableau du cours. (b) Vérifiez-la avec la formule ΔE = ½ m g² t Δt. (c) Que devient ΔE avec Δt = 0,1 s ?`),
      ],
      solution: [
        p(`(a) E = ½ × 2 × 20² − 2 × 10 × 15 = 400 − 300 = **100 J** (alors que l’énergie exacte est nulle).`),
        p(`(b) ΔE = ½ × 2 × 10² × 2 × 0,5 = **100 J** : la formule redonne bien le résultat.`),
        p(`(c) ΔE est proportionnelle au pas : avec Δt = 0,1 s, ΔE = 100 × 0,1 / 0,5 = **20 J**.`),
      ],
    },
    {
      level: 4,
      title: "Rejoindre 90 % de la vitesse limite",
      statement: [
        p(`Pour le parachutiste du cours (m = 80 kg, k = 20 kg/s, donc $v_{∞}$ = 40 m/s et τ = 4 s), au bout de combien de temps atteint-il 90 % de sa vitesse limite ? Quelle distance a-t-il alors parcourue ? On rappelle x(t) = $v_{∞}$ (t − τ (1 − e^(−t/τ))).`),
      ],
      solution: [
        p(`On résout $1 − e^{−t/τ} = 0,9$, soit $e^{−t/τ} = 0,1$ et $t = τ ln 10 = 4 × 2,303$ = **9,21 s**.`),
        p(`Distance : x = 40 × (9,21 − 4 × 0,9) = 40 × 5,61 = **224 m** environ.`),
      ],
    },
    {
      level: 5,
      title: "Stabilité d’Euler",
      statement: [
        p(`On étudie $v' = g − v/τ$ avec g = 10 m/s² et τ = 0,5 s (donc $v_{∞} = 5$ m/s), $v(0) = 0$. (1) Écrivez la récurrence d’Euler pour Δt = 0,1 s et calculez v à t = 0,5 s. (2) Comparez à la solution exacte. (3) À partir de quel pas la méthode devient-elle instable ? Que se passe-t-il pour Δt = 1,2 s ?`),
      ],
      solution: [
        p(`(1) $v_{n+1} = v_n + 0,1 (10 − 2 v_n) = 0,8 v_n + 1$, donc $v_n = 5 (1 − 0,8^n)$ et $v_5 = 5 (1 − 0,32768)$ = **3,3616 m/s**.`),
        p(`(2) Exacte : $5 (1 − e^{−1}) = 3,1606$ m/s. L’erreur est de +0,2010 m/s, soit **+6,4 %**.`),
        p(`(3) Le facteur d’amplification est $1 − Δt/τ$ ; il est de module inférieur à 1 si $Δt < 2τ = 1$ s. Pour Δt = 1,2 s le facteur vaut −1,4 : les valeurs successives sont 12 ; −4,8 ; 18,7 ; −14,2 ; 31,9 ; … La solution numérique oscille et diverge alors que la vraie solution tend vers 5 m/s.`),
      ],
    },
    {
      level: 5,
      title: "Lire un ordre sur des mesures",
      statement: [
        p(`Le logiciel a relevé les erreurs suivantes sur un oscillateur : Verlet 1,965 × 10⁻² (100 pas) et 4,896 × 10⁻³ (200 pas) ; RK4 1,171 × 10⁻⁴ (100 pas) et 7,318 × 10⁻⁶ (200 pas) ; Euler 4,599 × 10⁻² (3200 pas) et 2,274 × 10⁻² (6400 pas). Calculez l’ordre de chaque méthode par $p = log_{2}(e_1 / e_2)$ quand on double le nombre de pas. Combien de pas Euler faudrait-il pour égaler RK4 à 200 pas ?`),
      ],
      solution: [
        p(`Verlet : $e_1/e_2 = 4,015$, $p = 2,005 ≈ 2$. RK4 : rapport 15,996, $p = 4,00$. Euler : rapport 2,023, $p = 1,016 ≈ 1$.`),
        p(`L’erreur d’Euler varie comme 1/N. Pour atteindre 7,318 × 10⁻⁶ en partant de 2,274 × 10⁻² à 6400 pas, il faut multiplier le nombre de pas par 2,274 × 10⁻² / 7,318 × 10⁻⁶ = 3107, soit environ **2 × 10⁷ pas** (20 millions), contre 200 pas (800 évaluations) pour RK4 : un rapport voisin de 25 000.`),
      ],
    },
    {
      level: 6,
      title: "Contrôle de pas de RK45",
      statement: [
        p(`Un pas Δt = 0,2 s de RK45 produit une erreur locale estimée de 3,2 × 10⁻⁵ pour une tolérance de 10⁻⁶. (a) Le pas est-il accepté ? Quel pas suivant est proposé par la formule $Δt' = 0,9 Δt (tol/ε)^{1/5}$ ? (b) Même question si l’erreur estimée vaut 10⁻⁸. (c) Pourquoi le facteur 0,9 ?`),
      ],
      solution: [
        p(`(a) ε > tol : le pas est **refusé** et recommencé avec $Δt' = 0,9 × 0,2 × (10⁻⁶ / 3,2 × 10⁻⁵)^{1/5} = 0,9 × 0,2 × (1/32)^{1/5} = 0,9 × 0,2 × 0,5$ = **0,09 s**.`),
        p(`(b) ε < tol : le pas est accepté et le suivant peut être allongé : $Δt' = 0,9 × 0,2 × 100^{1/5} = 0,9 × 0,2 × 2,512$ = **0,452 s** (inférieur au plafond de 5 fois le pas).`),
        p(`(c) Le facteur 0,9 est une marge de sécurité : sans lui, le pas proposé serait juste à la limite de la tolérance et serait souvent refusé au coup suivant, ce qui coûte plus cher qu’un pas un peu prudent.`),
      ],
    },
    {
      level: 6,
      title: "Ordre contre coût",
      statement: [
        p(`À budget égal de 1600 évaluations de la force, comparez : Verlet avec mémoire (1 évaluation par pas) et RK4 (4 évaluations par pas), sachant que les erreurs mesurées sont 7,64 × 10⁻⁵ pour Verlet à 1600 pas et 4,57 × 10⁻⁷ pour RK4 à 400 pas. Quel est le rapport ? Pourquoi l’avantage de RK4 s’accroît-il quand on exige une précision plus grande ?`),
      ],
      solution: [
        p(`Rapport : 7,64 × 10⁻⁵ / 4,57 × 10⁻⁷ = **167** en faveur de RK4.`),
        p(`Pour une erreur visée ε, un schéma d’ordre $p$ demande $N ∝ ε^{−1/p}$ pas. Quand ε diminue, $N$ croît beaucoup plus vite pour $p = 2$ que pour $p = 4$ : par exemple pour gagner un facteur 10⁴ sur l’erreur il faut 100 fois plus de pas avec Verlet mais seulement 10 fois plus avec RK4. Le surcoût constant de 4 évaluations par pas est vite amorti.`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Un ordinateur calcule un mouvement **pas à pas** : il évalue la force, avance d’un petit pas de temps, et recommence.",
      "L’erreur dépend du pas et de la méthode : l’**ordre** dit de combien elle diminue quand on raccourcit le pas (÷2, ÷4, ÷16 pour les ordres 1, 2, 4).",
      "À précision égale une méthode d’ordre élevé coûte bien moins cher, mais il faut vérifier ce qu’elle conserve (énergie, structure) : chapitres suivants.",
      "Euler **fabrique de l’énergie** en chute libre (ΔE = ½ m g² t Δt) ; une méthode instable (Δt > 2τ) peut diverger.",
      "On mesure l’ordre **hors période et dans l’espace des phases**, jamais à un instant « trop gentil ».",
    ]),
  ]),
};
