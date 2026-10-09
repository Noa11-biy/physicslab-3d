const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 2,
  title: "Le ressort : oscillations, énergie, résonance",
  lead: "Presque tout ce qui vibre, d’une corde de guitare à un pont, se comporte près de son équilibre comme une masse au bout d’un ressort. C’est le système le plus étudié de la physique, et le meilleur banc d’essai pour juger ce que chaque méthode de calcul fait de l’énergie.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** ce qu’est une oscillation (période, amplitude), les trois façons de s’amortir, la résonance, et pourquoi l’énergie permet de juger un calcul.",
      "**Le résultat marquant :** Euler fait exploser l’énergie, Verlet la fait osciller sans dériver, RK4 la perd lentement.",
      "**Dans le logiciel :** module M2, simulation 2 ({{--sim 2}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`Accrochez une masse au bout d’un ressort, tirez-la vers le bas, lâchez : elle monte, redescend, remonte, encore et encore. À chaque instant de l’énergie change de forme. Quand la masse file au milieu, c’est de l’énergie de mouvement. Quand elle s’arrête tout en haut ou tout en bas, le ressort est étiré ou comprimé et l’énergie est « en réserve ». C’est exactement une balançoire.`),
      p(`Les frottements finissent par éteindre le mouvement. Mais si l’on pousse au bon rythme, comme un parent sur une balançoire, de petites poussées s’additionnent et le mouvement devient énorme : c’est la **résonance**. Les ingénieurs la redoutent pour les ponts et les bâtiments, et la recherchent pour les instruments de musique.`),
      p(`Pour l’ordinateur, le ressort est un excellent test, parce que l’on sait ce que l’énergie doit faire : **rester constante** (quand il n’y a pas de frottement). Si le calcul est mal fait, la masse monte un peu plus haut à chaque va-et-vient, comme si quelqu’un la poussait en cachette. Les méthodes les plus simples font exactement cela.`),
    ],
    2: [
      p(`Trois mots décrivent une oscillation : la **période** (durée d’un aller-retour), la **fréquence** (nombre d’allers-retours par seconde) et l’**amplitude** (de combien on s’écarte). Une masse plus lourde oscille plus lentement ; un ressort plus raide, plus vite. Fait remarquable : la période **ne dépend pas de l’amplitude** (les horloges à pendule reposent sur la même propriété).`),
      p(`Quand il y a du frottement, trois comportements sont possibles. **Peu de frottement** : le ressort oscille de moins en moins fort avant de s’arrêter (corde de guitare). **Frottement fort** : il revient lentement vers le repos sans jamais dépasser (porte lourde avec un ferme-porte trop dur). Entre les deux, le **régime critique** : retour au repos le plus vite possible, sans oscillation ; l’amortisseur d’une voiture est réglé tout près de ce régime.`),
      p(`Quand une force extérieure répétitive agit sur le ressort, la réponse est maximale si le rythme est celui du ressort : c’est la **résonance**. Sans frottement, l’amplitude grandirait indéfiniment.`),
      p(`Côté calcul, le ressort révèle la personnalité de chaque méthode. **Euler** ajoute un peu d’énergie à chaque pas : le mouvement enfle. **Verlet** fait osciller l’énergie autour de sa vraie valeur sans jamais dériver. **RK4** est très précis mais perd très lentement de l’énergie.`),
      fig("fig_energie.png", "Énergie d’un ressort suivie pendant 30 périodes avec 20 pas par période (mesures du logiciel). Euler multiplie l’énergie par plus de 10²⁴ ; Euler symplectique et Verlet la font osciller sans dériver ; RK4 en perd 0,8 % en 30 périodes.", 1.0),
    ],
    3: [
      p(`Pour une masse $m$ (en kg) accrochée à un ressort de raideur $k$ (en N/m), la période d’oscillation vaut :`),
      eq(`T = 2π sqrt{frac{m}{k}}`),
      p(`La fréquence est $f = 1/T$ (en hertz). L’énergie emmagasinée pour une amplitude $A$ est $E = ½ k A²$. Voici comment les choses varient :`),
      tbl(
        ["Si on…", "La période…", "L’énergie…"],
        [
          ["multiplie la masse par 4", "est multipliée par 2", "ne change pas"],
          ["multiplie la raideur par 4", "est divisée par 2", "est multipliée par 4"],
          ["double l’amplitude", "ne change pas", "est multipliée par 4"],
        ],
        [3600, 2900, 2798],
        { align: ["l", "c", "c"] }
      ),
      box("exemple", "Exemple", [
        p(`Masse m = 0,5 kg, ressort k = 200 N/m, amplitude A = 4 cm = 0,04 m.`),
        ul([
          "Période : T = 2π × √(0,5 / 200) = 2π × 0,05 = **0,314 s**, soit une fréquence de 3,18 Hz.",
          "Énergie : E = ½ × 200 × 0,04² = **0,16 J**.",
          "Au passage au milieu, toute l’énergie est cinétique : ½ × 0,5 × v² = 0,16, donc v = **0,8 m/s**.",
        ]),
      ]),
    ],
    4: [
      p(`La force de rappel est $F = −k x$ et la deuxième loi de Newton donne $m x'' = −k x$. Avec la pulsation propre $ω_0 = sqrt{k/m}$ (en rad/s), la solution est :`),
      eq(`x(t) = A cos(ω_0 t + φ)`),
      p(`avec $T = 2π/ω_0$. La vitesse maximale est $v_{max} = A ω_0$ et l’**énergie mécanique** se conserve :`),
      eq(`E = frac{1}{2} m v² + frac{1}{2} k x² = frac{1}{2} k A²`),
      p(`Avec un frottement visqueux $F = −c v$ l’équation devient $m x'' = −k x − c x'$ ; avec une force excitatrice $F_0 cos(Ω t)$, elle devient $m x'' = −k x − c x' + F_0 cos(Ω t)$. Sans frottement et à $Ω = ω_0$, la solution particulière $x = (F_0 / 2 m ω_0) t sin(ω_0 t)$ grandit **proportionnellement au temps** : c’est la résonance parfaite.`),
      p(`Lire les graphes : la courbe de $x(t)$ est une sinusoïde ; le **portrait de phase** (la vitesse en fonction de la position) est une ellipse parcourue en une période, qui rétrécit en spirale quand il y a du frottement.`),
      box("attention", "Attention", [
        p(`Dans un ressort vertical il faut distinguer la longueur à vide et la position d’équilibre : la pesanteur décale simplement l’équilibre de $m g / k$ et ne change ni la période ni la forme des oscillations autour de ce nouvel équilibre.`),
      ]),
    ],
    5: [
      p(`On pose $γ = c/(2m)$ et le **taux d’amortissement réduit** $ζ = c / (2 sqrt{k m}) = γ/ω_0$. Le régime est sous-critique si $ζ < 1$, critique si $ζ = 1$, sur-critique si $ζ > 1$. Dans le premier cas :`),
      eq(`x(t) = e^{−γ t} (P cos(ω_d t) + Q sin(ω_d t)), ω_d = sqrt{ω_0² − γ²}`),
      p(`En régime forcé permanent, l’amplitude est :`),
      eq(`X(Ω) = frac{F_0 / m}{sqrt{(ω_0² − Ω²)² + (2γΩ)²}}`),
      p(`Elle passe par un maximum en $Ω = ω_0 sqrt{1 − 2ζ²}$ lorsque $ζ < 1/sqrt{2}$.`),
      p(`**Ce que fait chaque schéma.** Écrivons la méthode d’Euler sur le ressort ($y = (x, v)$) :`),
      eq(`x_{n+1} = x_n + v_n Δt, v_{n+1} = v_n − ω_0² x_n Δt`),
      p(`On vérifie par un calcul direct que $ω_0² x² + v²$ est multipliée à chaque pas par $1 + ω_0² Δt²$ : l’énergie croît **toujours**, quel que soit le pas. Pour 20 pas par période ($ω_0 Δt = 2π/20$) et 30 périodes, la formule $(1 + ω_0² Δt²)^{600}$ prédit $3,3615 × 10²⁴$ ; le logiciel mesure $3,3615 × 10²⁴$ (figure du niveau 2).`),
      p(`Pour Euler symplectique ($v$ d’abord, puis $x$ avec la nouvelle vitesse), la quantité $ω_0² x² + v² − Δt ω_0² x v$ est conservée **exactement** (vérifié à $6 × 10^{−14}$ près sur 2000 pas) : l’énergie vraie $½ (v² + ω_0² x²)$ n’en diffère que par le terme croisé, d’où une oscillation de 0,864 à 1,186 fois sa valeur initiale, sans dérive. Verlet des vitesses se comporte de même, avec une oscillation bien plus petite (0,975 à 1) car il est d’ordre 2.`),
      tbl(
        ["Schéma", "Énergie sur 30 périodes (20 pas/période)", "Comportement"],
        [
          ["Euler explicite", "×3,4 × 10²⁴", "croissance exponentielle"],
          ["Euler symplectique", "de 0,864 à 1,186", "oscille, sans dérive"],
          ["Verlet des vitesses", "de 0,975 à 1,000", "oscille, sans dérive"],
          ["RK4", "0,9921 à la fin", "décroissance lente et régulière"],
        ],
        [2600, 3500, 3198],
        { align: ["l", "c", "l"] }
      ),
    ],
    6: [
      p(`**RK4 sur l’oscillateur.** Le facteur d’amplification est le polynôme de stabilité $R(z) = 1 + z + z²/2 + z³/6 + z⁴/24$ évalué en $z = i y$, avec $y = ω_0 Δt$. Un calcul direct donne :`),
      eq(`|R(i y)|² = 1 − frac{y⁶}{72} + frac{y⁸}{576}`),
      p(`L’énergie est donc multipliée à chaque pas par ce facteur, un peu inférieur à 1 : la dissipation numérique est d’ordre $(ω_0 Δt)^6$. Pour 20 pas par période ($y = 2π/20$), $|R|² = 1 − 1,3188 × 10^{−5}$ par pas, soit $0,99212$ après 600 pas : exactement ce que mesure le logiciel (0,992118). Le schéma est stable tant que $|R| ≤ 1$, c’est-à-dire $ω_0 Δt ≤ 2 sqrt{2} ≈ 2,83$ ; Euler explicite, lui, n’est stable pour aucun pas sur l’axe imaginaire ($|1 + i y|² = 1 + y² > 1$).`),
      p(`**Verlet et fréquence numérique.** Le schéma de Verlet des vitesses est équivalent au « saute-mouton » : il tourne dans l’espace des phases avec une pulsation numérique $ω_{num}$ telle que $cos(ω_{num} Δt) = 1 − (ω_0 Δt)²/2$, soit $ω_{num} = (2/Δt) arcsin(ω_0 Δt / 2)$, définie tant que $ω_0 Δt < 2$. Le rapport $ω_{num}/ω_0 ≈ 1 + (ω_0 Δt)²/24$ vaut 1,0042 à 20 pas par période : l’orbite numérique tourne **un peu trop vite** (erreur de phase de 0,78 rad après 30 périodes), mais l’énergie ne dérive pas.`),
      p(`**Analyse par erreur inverse.** Un schéma symplectique d’ordre $p$ est, pour un problème analytique, la solution exacte (à des termes exponentiellement petits près) d’un hamiltonien modifié $text{H̃} = H + Δt^p H_p + …$. L’énergie vraie $H$ ne s’écarte donc de la constante $H̃$ que de $O(Δt^p)$ pendant des durées très longues. C’est ce qui sépare Verlet (ordre 2, borné) de RK4 (ordre 4 mais dissipatif) pour les simulations longues d’un système hamiltonien (chapitre 4).`),
      box("plus", "Pour aller plus loin", [
        p(`Le critère de choix d’un schéma n’est donc pas « le plus grand ordre » mais « ce qu’il faut conserver sur la durée voulue » : une trajectoire précise sur quelques périodes (RK4, RK45) ou une énergie qui ne dérive pas sur des millions de périodes (schémas symplectiques).`),
      ]),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Lancez {{physicslab --level 4 --sim 2}} : un ressort-masse avec ses méthodes numériques.`),
      ul([
        "Choisissez **Libre**, **Sous-amorti**, **Critique** ou **Sur-amorti** pour comparer les trois régimes ; la fenêtre **Invariants** donne l’énergie de chaque méthode.",
        "Activez la **force périodique** et faites varier sa fréquence autour de la pulsation propre ω₀ : la courbe de résonance se dessine.",
        "Avec le niveau 5 ou 6, augmentez le pas de calcul : regardez l’énergie d’Euler s’emballer, celle de Verlet osciller et celle de RK4 s’éroder.",
      ]),
    ]),
    fig("shot_s02.png", "La simulation 2 au niveau 4 : le ressort en trois dimensions, le tableau des énergies et, en bas, le portrait de phase des quatre méthodes.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "La balançoire",
      statement: [
        p(`Un enfant est assis sur une balançoire. (a) Quelle est la bonne façon de le pousser pour qu’il aille de plus en plus haut ? (b) Que se passe-t-il si on le pousse au hasard, sans tenir compte de son mouvement ? (c) Comment s’appelle le phénomène de la question (a) ?`),
      ],
      solution: [
        p(`(a) Pousser **à chaque aller-retour, au bon moment**, au rythme naturel de la balançoire. (b) Les poussées tombent tantôt dans le sens du mouvement, tantôt contre lui : elles se compensent à peu près et l’amplitude reste petite. (c) La **résonance**.`),
      ],
    },
    {
      level: 2,
      title: "Trois régimes",
      statement: [
        p(`Associez chaque situation à un régime d’amortissement (sous-amorti, critique, sur-amorti) : (a) une corde de guitare pincée ; (b) un ferme-porte réglé trop fort, qui ramène la porte très lentement sans la faire claquer ; (c) un amortisseur de voiture neuf qui ramène la caisse à sa position sans oscillation notable.`),
      ],
      solution: [
        p(`(a) **Sous-amorti** : la corde oscille longtemps en s’éteignant. (b) **Sur-amorti** : retour très lent, sans dépassement. (c) **Proche du régime critique** : retour au repos le plus rapide sans oscillation.`),
      ],
    },
    {
      level: 3,
      title: "Période et énergie",
      statement: [
        p(`Un ressort de raideur k = 200 N/m porte une masse m = 0,5 kg et on l’écarte de A = 4 cm. Calculez la période, la fréquence, l’énergie et la vitesse maximale. Que deviennent la période et l’énergie si l’on écarte de 8 cm ?`),
      ],
      solution: [
        p(`T = 2π √(m/k) = 2π × √(0,5/200) = **0,314 s** ; f = 1/T = **3,18 Hz** ; E = ½ k A² = ½ × 200 × 0,04² = **0,16 J** ; $ω_0 = sqrt{k/m} = 20$ rad/s et $v_{max} = A ω_0 = 0,04 × 20$ = **0,8 m/s**.`),
        p(`Avec A = 8 cm la période est **inchangée** (0,314 s) et l’énergie est multipliée par 4 : **0,64 J** ; la vitesse maximale double (1,6 m/s).`),
      ],
    },
    {
      level: 3,
      title: "Trouver la raideur",
      statement: [
        p(`Une masse de 0,2 kg allonge un ressort vertical de 5 cm à l’équilibre (g = 9,81 m/s²). Calculez k, puis la période d’oscillation de cette masse.`),
      ],
      solution: [
        p(`À l’équilibre le poids compense la force de rappel : $k x = m g$, donc $k = 0,2 × 9,81 / 0,05$ = **39,2 N/m**.`),
        p(`Période : $T = 2π sqrt{m/k} = 2π sqrt{0,2 / 39,24}$ = **0,449 s**.`),
      ],
    },
    {
      level: 4,
      title: "Énergie à mi-chemin",
      statement: [
        p(`Pour le ressort de l’exercice 2.3 (A = 4 cm, ω₀ = 20 rad/s), calculez la vitesse quand la masse passe à x = A/2, puis la position où l’énergie cinétique est égale à l’énergie potentielle.`),
      ],
      solution: [
        p(`Conservation de l’énergie : $½ m v² + ½ k x² = ½ k A²$ donne $v = ω_0 sqrt{A² − x²}$. Pour $x = A/2$ : $v = 20 × 0,04 × sqrt{1 − 1/4} = 0,8 × 0,866$ = **0,693 m/s**.`),
        p(`Énergies égales : $½ k x² = ¼ k A²$, donc $x = A/sqrt{2} = 0,04/1,414$ = **2,83 cm**.`),
      ],
    },
    {
      level: 5,
      title: "Un oscillateur amorti",
      statement: [
        p(`Un oscillateur a m = 1 kg, k = 10 N/m et c = 2 kg/s. Calculez le taux d’amortissement réduit ζ, le régime, la pseudo-pulsation ω_d, la pseudo-période, et le facteur par lequel l’amplitude est multipliée en une pseudo-période.`),
      ],
      solution: [
        p(`$ω_0 = sqrt{10} = 3,162$ rad/s, $γ = c/(2m) = 1$ s⁻¹, $ζ = γ/ω_0$ = **0,316** : régime **sous-critique** (ζ < 1).`),
        p(`$ω_d = sqrt{ω_0² − γ²} = sqrt{10 − 1}$ = **3 rad/s** ; pseudo-période $T_d = 2π/3$ = **2,094 s**.`),
        p(`L’amplitude décroît en $e^{−γ t}$ : en une pseudo-période le facteur vaut $e^{−2,094}$ = **0,123**, soit une perte de près de 88 %.`),
      ],
    },
    {
      level: 5,
      title: "Euler et l’énergie",
      statement: [
        p(`On simule un oscillateur de période T = 1 s (ω₀ = 2π rad/s) avec Euler explicite et un pas Δt = 0,01 s. (a) De quel facteur l’énergie est-elle multipliée à chaque pas, en une période, en dix périodes ? (b) Quel pas faudrait-il pour que l’énergie ne gagne pas plus de 1 % en dix périodes ? Combien de pas cela représente-t-il ?`),
      ],
      solution: [
        p(`(a) Facteur par pas : $1 + ω_0² Δt² = 1 + (2π × 0,01)²$ = **1,003948**. En une période (100 pas) : 1,003948¹⁰⁰ = **1,483** ; en dix périodes (1000 pas) : **51,4**. L’énergie est multipliée par plus de 50 !`),
        p(`(b) Pour un gain de 1 % sur 10 s : $(1 + ω_0² Δt²)^{10/Δt} ≈ e^{10 ω_0² Δt} ≤ 1,01$, donc $Δt ≤ ln(1,01)/(10 ω_0²)$ = **2,5 × 10⁻⁵ s**, soit environ **400 000 pas** pour dix périodes (au lieu de 1000) : Euler est un très mauvais choix pour un oscillateur.`),
      ],
    },
    {
      level: 6,
      title: "Prédire l’énergie de RK4",
      statement: [
        p(`(a) Pour un oscillateur simulé par RK4 avec 20 pas par période pendant 30 périodes, prédisez le rapport E/E₀ final à partir de $|R(i y)|² = 1 − y⁶/72 + y⁸/576$ et comparez à la valeur mesurée 0,992118. (b) Combien de pas par période faut-il pour que la perte d’énergie ne dépasse pas 10⁻⁶ par période ?`),
      ],
      solution: [
        p(`(a) $y = 2π/20 = 0,31416$ : $y⁶ = 9,614 × 10^{−4}$, $y⁸ = 9,49 × 10^{−5}$, donc $|R|² = 1 − 1,3353 × 10^{−5} + 1,65 × 10^{−7} = 1 − 1,3188 × 10^{−5}$. Après 600 pas : $(1 − 1,3188 × 10^{−5})^{600} = e^{−7,913 × 10^{−3}}$ = **0,99212**, en accord avec la valeur mesurée (0,992118).`),
        p(`(b) Perte par période ≈ $(2π/y) × y⁶/72 = 2π y⁵/72 ≤ 10^{−6}$, soit $y⁵ ≤ 1,146 × 10^{−5}$ et $y ≤ 0,1028$ : il faut $2π/y$ ≈ **61 pas par période** au minimum. Avec Euler, en suivant le même raisonnement, il en faudrait des centaines de milliers.`),
      ],
    },
    {
      level: 6,
      title: "La fréquence numérique de Verlet",
      statement: [
        p(`Le schéma de Verlet des vitesses fait tourner l’oscillateur à la pulsation numérique $ω_{num} = (2/Δt) arcsin(ω_0 Δt / 2)$. Calculez le rapport $ω_{num}/ω_0$ pour 20 pas par période, l’erreur de phase accumulée en 30 périodes, puis la condition de stabilité.`),
      ],
      solution: [
        p(`Avec $ω_0 Δt = 2π/20 = 0,31416$ : $arcsin(0,15708) = 0,157726$ et $ω_{num}/ω_0 = 2 × 0,157726 / 0,31416$ = **1,00416**.`),
        p(`En 30 périodes, la phase exacte vaut $30 × 2π = 188,5$ rad et la phase numérique la dépasse de $188,5 × 0,00416$ = **0,78 rad** (environ 45°) : l’orbite numérique est en avance, mais son énergie reste bornée.`),
        p(`La pulsation numérique n’existe (le sinus inverse est défini) que si $ω_0 Δt/2 ≤ 1$, soit **$ω_0 Δt ≤ 2$** : c’est la condition de stabilité de Verlet (soit au moins π ≈ 3,14 pas par période).`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Un ressort oscille avec la période $T = 2π sqrt{m/k}$, **indépendante de l’amplitude** ; son énergie $½ k A²$ se conserve s’il n’y a pas de frottement.",
      "Avec frottement : trois régimes (sous-critique, critique, sur-critique) selon $ζ = c / (2 sqrt{k m})$ ; en régime forcé, l’amplitude est maximale près de la pulsation propre (**résonance**).",
      "L’énergie est un **détecteur d’erreur** : Euler la multiplie par $1 + ω_0² Δt²$ à chaque pas, Verlet et Euler symplectique la font osciller sans dérive, RK4 la perd lentement (en $(ω_0 Δt)^6$ par pas).",
      "Une bonne méthode n’est pas seulement précise : elle doit **respecter ce qui compte** pour la durée visée.",
    ]),
  ]),
};
