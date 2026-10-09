const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 5,
  title: "Frottements, chocs et rebonds",
  lead: "Quand une balle touche le sol ou qu’un bloc se met à glisser, la physique change brutalement d’un instant à l’autre. Ces ruptures sont le cauchemar des méthodes numériques les plus élégantes.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** le frottement sec, les rebonds avec perte d’énergie, le paradoxe de Zénon (une infinité de rebonds en un temps fini) et le berceau de Newton.",
      "**Le résultat marquant :** dès que la physique est discontinue, toutes les méthodes retombent à la précision d’Euler, sauf si l’on **détecte l’événement** : à 400 pas, l’erreur passe de $6,7 × 10^{−2}$ à $3,3 × 10^{−10}$, soit 200 millions de fois mieux.",
      "**Dans le logiciel :** module M5, simulations 7 (frottement sec, {{--sim 7}}), 8 (chocs et rebonds, {{--sim 8}}) et 9 (berceau de Newton, {{--sim 9}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`**Le bloc sur la planche.** Posez un bloc sur une planche et soulevez-la doucement : rien ne bouge, le bloc « colle ». À un certain angle, d’un coup, il se met à glisser. Il y a donc deux mondes, l’adhérence et le glissement, et un instant précis où l’on passe de l’un à l’autre. L’ordinateur doit trouver cet instant.`),
      p(`**La balle qui rebondit.** Lâchée de haut, elle rebondit de moins en moins haut, de plus en plus vite, une infinité de fois, et pourtant elle s’arrête au bout d’un temps **fini**. C’est le paradoxe de Zénon, version balle de tennis. Un calcul trop naïf, lui, n’arrive jamais à l’arrêter : la balle continue à trembloter.`),
      p(`**Le berceau de Newton.** Vous connaissez ces cinq billes suspendues : on lâche une bille, une bille repart de l’autre côté. Presque : si l’on regarde de très près, les premières billes **reculent** un peu et la dernière repart un peu moins vite. L’idée simple « les billes s’échangent leur mouvement » est une bonne approximation, pas la vérité exacte.`),
      fig("fig_rebonds.png", "Une balle lâchée de 2 m avec un coefficient de rebond de 0,8 : une infinité de rebonds de plus en plus petits, mais tout est terminé au bout de 5,75 s (formule exacte).", 0.82),
    ],
    2: [
      p(`**Le frottement sec** (celui d’un solide qui glisse) obéit à une règle simple : tant que la force qui cherche à faire bouger le bloc reste inférieure à une limite, il **adhère** ; au-delà il glisse, freiné par une force constante. Pour calculer proprement, il faut repérer l’instant exact où le bloc s’arrête ou repart : on appelle cela un **événement**.`),
      p(`**Le rebond** est décrit par un seul nombre, le **coefficient de restitution** $e$ : après le choc, la vitesse relative est $e$ fois celle d’avant. $e = 1$ : choc parfaitement élastique ; $e = 0$ : choc mou, sans rebond. Une balle lâchée de 2 m avec $e = 0,8$ remonte à 64 % de sa hauteur au premier rebond ($e²$), puis à 41 %, et ainsi de suite.`),
      p(`**Deux philosophies pour traiter un choc.** On peut le considérer comme **instantané** : à l’instant du contact on change brutalement les vitesses (la méthode des « impulsions »). Ou bien on suppose que les billes s’écrasent un peu comme des ressorts très raides (le **contact de Hertz**) et l’on laisse le calcul suivre cette courte phase. La seconde est plus réaliste et plus coûteuse ; la première est plus simple et suffit souvent. Dans le berceau de Newton elles donnent des résultats **différents** : avec les impulsions, « autant de billes qui entrent que de billes qui sortent » ; avec Hertz, les billes de tête reculent d’environ 7 % de la vitesse de lancement.`),
    ],
    3: [
      p(`**Frottement.** La force de frottement qui s’oppose au glissement vaut $F = μ N$, où $N$ est la force d’appui du plan sur le bloc et $μ$ le coefficient de frottement (sans unité). Sur un plan incliné d’un angle $θ$, un bloc commence à glisser si $tan θ > μ$ ; sinon il reste immobile. Quand il glisse, son accélération vaut $a = g (sin θ − μ cos θ)$.`),
      box("exemple", "Exemple", [
        p(`Plan incliné de 20° (tan 20° = 0,364). Avec μ = 0,3, le bloc glisse (0,364 > 0,3) avec a = 9,81 × (sin 20° − 0,3 × cos 20°) = **0,59 m/s²**. Avec μ = 0,4 il ne bouge pas (0,364 < 0,4).`),
      ]),
      p(`**Rebonds.** Si une balle arrive au sol à la vitesse $v$ elle repart à $e × v$, donc remonte à la hauteur $e² × h$. Après $n$ rebonds la hauteur maximale vaut $h_n = e^{2n} × h_0$.`),
      box("exemple", "Exemple : une balle lâchée de 1,6 m, e = 0,5", [
        ul([
          "Vitesse à l’arrivée : v = √(2 × 9,81 × 1,6) = **5,60 m/s**.",
          "Vitesse après le premier rebond : 0,5 × 5,60 = **2,80 m/s**.",
          "Hauteur après le premier rebond : 0,5² × 1,6 = **0,4 m** ; après le deuxième, 0,1 m.",
        ]),
      ]),
    ],
    4: [
      p(`**Plan incliné.** Pour un bloc lancé vers le haut du plan à la vitesse $v_0$, la pesanteur et le frottement freinent ensemble : $a = g (sin θ + μ cos θ)$. Il s’arrête après la distance $d = v_0² / (2a)$ et la durée $v_0 / a$ ; il redescend ensuite seulement si $tan θ > μ_s$ (frottement d’adhérence).`),
      box("exemple", "Exemple", [
        p(`Plan de 30°, μ = 0,2, v₀ = 4 m/s : a = 9,81 × (0,5 + 0,2 × 0,866) = **6,60 m/s²** ; le bloc monte de d = 16 / (2 × 6,60) = **1,21 m** pendant 0,61 s.`),
      ]),
      p(`**Durée totale des rebonds.** Si la balle tombe de la hauteur $h_0$ sa première chute dure $t_0 = sqrt{2 h_0 / g}$. Un vol complet après le $n$-ième choc dure $2 t_0 e^n$. La somme est une série géométrique :`),
      eq(`t_{total} = t_0 + 2 t_0 sum_{n=1}^{∞} e^n = t_0 frac{1 + e}{1 − e}`),
      p(`Pour $h_0 = 2$ m et $e = 0,8$ : $t_0 = 0,6386$ s et $t_{total} = 0,6386 × 9 = 5,748$ s, exactement la durée du graphique du niveau 1. Pour $h_0 = 1,6$ m et $e = 0,5$ : $t_0 = 0,571$ s et $t_{total} = 1,713$ s. La hauteur est inférieure à 1 mm après 6 rebonds.`),
      p(`**Conservation de l’impulsion.** Dans un choc entre deux billes, l’impulsion totale $m_1 v_1 + m_2 v_2$ se conserve toujours, même lorsque l’énergie cinétique ne se conserve pas ($e < 1$).`),
    ],
    5: [
      p(`**Choc à une dimension avec restitution $e$.** L’impulsion se conserve et $v_2' − v_1' = −e (v_2 − v_1)$, d’où :`),
      eq(`v_1' = frac{(m_1 − e m_2) v_1 + (1 + e) m_2 v_2}{m_1 + m_2}, v_2' = frac{(m_2 − e m_1) v_2 + (1 + e) m_1 v_1}{m_1 + m_2}`),
      p(`L’énergie perdue vaut $ΔE = ½ μ_r (1 − e²) (v_1 − v_2)²$ avec $μ_r = m_1 m_2 / (m_1 + m_2)$. Pour $m_1 = 2$ kg, $v_1 = 3$ m/s, $m_2 = 1$ kg au repos et $e = 0,5$ : $v_1' = 1,5$ m/s, $v_2' = 3$ m/s, énergie perdue 2,25 J sur 9 J (25 %).`),
      p(`**Contact de Hertz.** Deux billes élastiques se repoussent avec la force $F = k δ^{3/2}$ où $δ$ est l’écrasement. Pour deux billes égales de masse $m$ de vitesse relative $v$ et de masse réduite $m_r = m/2$, la conservation de l’énergie donne l’écrasement maximal, puis la durée du contact :`),
      eq(`δ_{max} = paren{frac{5 m_r v²}{4 k}}^{2/5}, T = frac{4}{5} frac{Γ(2/5) Γ(1/2)}{Γ(9/10)} frac{δ_{max}}{v} = 2,943275 frac{δ_{max}}{v}`),
      p(`Pour $m = 1$, $k = 10^4$, $v = 1$ (unités normalisées du logiciel) on obtient $δ_{max} = 0,0208138$ et $T = 0,0612608$, retrouvés par le calcul numérique à $6 × 10^{−13}$ et $1,6 × 10^{−13}$ près. La durée du contact varie comme $v^{−1/5}$ : presque indépendante de la vitesse. Les résultats d’une chaîne de billes ne dépendent **ni de $k$ ni de $v$** : une seule échelle de longueur $(m v²/k)^{2/5}$ intervient.`),
      p(`**Trois façons de traiter le frottement sec** (erreur à 400 pas, RK4) :`),
      tbl(
        ["Modèle", "Principe", "Erreur à 400 pas", "Ordre observé"],
        [
          ["Naïf", "{{sgn(v)}} dans l’équation", "6,7 × 10⁻²", "1 pour tous les schémas"],
          ["Régularisé", "$tanh(v/ε)$", "0,017 à 0,35 m selon ε", "plafond indépendant de $Δt$"],
          ["Événement + adhérence", "on s’arrête à l’instant exact", "3,3 × 10⁻¹⁰", "ordre normal du schéma"],
        ],
        [1900, 2900, 2300, 2198],
        { align: ["l", "l", "c", "l"], size: 19 }
      ),
      box("attention", "Attention : le piège de l’adhérence", [
        p(`Un bloc qui devrait **rester collé** après un arrêt (μ_d < tan θ ≤ μ_s) est parti à −9,2 m en 5 s avec le modèle naïf : il ignore tout bonnement l’adhérence. Ce n’est pas une erreur de calcul qui disparaît avec le pas : c’est une erreur de **modèle**. Même quand la vitesse s’annule exactement, le modèle naïf garde une vitesse résiduelle proportionnelle au pas, jamais nulle.`),
      ]),
    ],
    6: [
      p(`**Détection d’événement.** On cherche l’instant où une fonction $g(t, y)$ change de signe. Dès qu’un pas voit le changement de signe, on bissecte le pas ; partant d’un pas de 0,05 s, il faut $log_{2}(0,05/10^{−10}) = 28,9$, soit 29 itérations (on arrondit à l’entier supérieur), pour localiser l’instant à $10^{−10}$ près. L’instant d’arrêt d’un bloc est ainsi trouvé à $10^{−6}$ (RK4) et $10^{−8}$ (RK45). La fonction doit ignorer un départ exactement sur la surface ($g = 0$) : on prend $g = y$ au-dessus du sol et $g = 0$ une fois posée, sinon un vol entier contenu dans un pas est manqué et la balle traverse le sol (observé jusqu’à −26 m).`),
      p(`**Pourquoi la discontinuité détruit l’ordre.** Un schéma d’ordre $p$ suppose une solution régulière sur le pas ; un saut de la force à l’intérieur du pas donne une erreur locale en $O(h)$ ; accumulée sur les événements, l’erreur globale est d’ordre 1 pour tous les schémas : pente log-log de 0,94 pour RK4 sur les rebonds naïfs, contre 4,0 avec l’événement, et un instant de repos retrouvé à 2·10⁻⁶ ($Δt = 0,05$) puis 1,5·10⁻¹² ($Δt = 0,001$). RK45 sur un bloc qui doit rester collé **ne termine pas** : l’erreur locale reste d’ordre $h$ sur la surface $v = 0$, le pas tombe jusqu’à $h_{min} = 10^{−14}$ ; d’où la règle de toujours borner le nombre de pas.`),
      p(`**Zénon numérique.** Avec le modèle naïf la balle ne s’arrête jamais : 0,5 s après le repos exact, $|v_y|$ atteint encore 0,78 / 0,28 / 0,11 m/s pour $Δt$ = 0,02 / 0,01 / 0,005 s, soit 2,3 à 4 fois $g Δt$. Avec événement, le seuil d’arrêt (vitesse de repos) doit rester ≥ $10^{−8}$ m/s.`),
      p(`**Contact de Hertz et ordre.** La force $k δ^{3/2}$ n’a pas de dérivée seconde en $δ = 0$ : RK4 n’atteint pas l’ordre 4 mais 2,3 à 2,75 selon la phase du contact sur la grille de pas ; Euler (ordre 1) et Verlet (ordre 2, avec Euler symplectique à 5 × 10⁻¹³ près) sont retrouvés. Euler explicite **diverge** à 10 pas par contact.`),
      p(`**Amortissement de Hunt-Crossley** $F = k δ^{3/2} (1 + (3/2) α δ')$ : le coefficient de restitution vaut $e ≈ 1/(1 + α v)$ (0,990099 / 0,952370 / 0,909016 / 0,768 pour $α$ = 0,01 / 0,05 / 0,1 / 0,3 à $v = 1$), et la perte relative d’énergie ≈ $2 α v$ au premier ordre.`),
      box("attention", "Attention : une idée fausse corrigée par la mesure", [
        p(`On aurait pu croire que l’**ordre** dans lequel on résout les chocs successifs d’un berceau change le résultat. Sur 245 cas (de 3 à 9 billes, avec $e$ de 0,99 à 0), l’écart est **exactement zéro** : la causalité impose la séquence. L’ordre ne compte que pour une bille prise entre deux voisines qui s’approchent déjà (écart 0,070 pour $e = 0,5$).`),
      ]),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Frottement sec : {{physicslab --level 4 --sim 7}}. Rebonds : {{physicslab --level 4 --sim 8}}. Berceau de Newton : {{physicslab --level 5 --sim 9}}.`),
      ul([
        "Sur le plan incliné, choisissez les cas types (**Pente douce**, **Pente forte**, **Descente freinée**, **Piège : adhérence**) et comparez les trois modèles : naïf, régularisé et événement.",
        "Pour les rebonds, comparez le modèle **Naïf** (contact vu après le pas) à l’**Événement** : la balle du modèle naïf ne s’arrête jamais. Le graphe « Durée des vols (Zénon) » montre la suite des vols.",
        "Au berceau, sélectionnez **1 sur 4**, **2 sur 3** ou **2 sur 2 (paradoxe)**, puis **Les deux (comparaison)** pour voir les impulsions contre Hertz : repérez le recul des billes de tête.",
      ]),
    ]),
    fig("shot_s09.png", "Le berceau de Newton au niveau 5 : les billes, leurs vitesses finales (Hertz contre impulsions) à droite, la convergence des schémas en bas.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "Zénon et la balle",
      statement: [
        p(`Une balle rebondit une infinité de fois. (a) Rebondit-elle pour toujours ? (b) Un ordinateur naïf, qui ne regarde le sol qu’à la fin de chaque pas, la voit-il s’arrêter ?`),
      ],
      solution: [
        p(`(a) **Non** : bien qu’il y ait une infinité de rebonds, ils sont de plus en plus courts, et leur durée totale est finie : la balle s’arrête en un temps fini. (b) **Non** : à la fin de chaque pas il manque la subtilité des tout petits rebonds, la balle garde une petite vitesse parasite et continue à trembloter indéfiniment.`),
      ],
    },
    {
      level: 3,
      title: "Le plan incliné",
      statement: [
        p(`Un bloc est posé sur un plan incliné de 20°. Avec un coefficient de frottement μ = 0,3, glisse-t-il ? Avec quelle accélération ? Que se passe-t-il avec μ = 0,4 ?`),
      ],
      solution: [
        p(`tan 20° = 0,364 > 0,3 : le bloc **glisse**, avec $a = g (sin θ − μ cos θ) = 9,81 × (0,342 − 0,3 × 0,940)$ = **0,59 m/s²**.`),
        p(`Avec μ = 0,4, tan 20° = 0,364 < 0,4 : il **reste immobile** (le frottement suffit à le retenir).`),
      ],
    },
    {
      level: 3,
      title: "Rebond à 50 %",
      statement: [
        p(`Une balle de coefficient de restitution e = 0,5 est lâchée de 1,6 m. Calculez sa vitesse à l’arrivée, sa vitesse après le premier rebond, la hauteur atteinte après le premier puis le deuxième rebond.`),
      ],
      solution: [
        p(`$v = sqrt{2 g h} = sqrt{2 × 9,81 × 1,6}$ = **5,60 m/s** ; après le rebond : $e v$ = **2,80 m/s**.`),
        p(`Hauteurs : $e² h = 0,25 × 1,6$ = **0,4 m** puis $e⁴ h$ = **0,1 m**.`),
      ],
    },
    {
      level: 4,
      title: "La durée du paradoxe",
      statement: [
        p(`La balle de l’exercice précédent (h₀ = 1,6 m, e = 0,5) rebondit jusqu’à s’arrêter. (a) Calculez la durée totale par la formule $t_0 (1 + e)/(1 − e)$. (b) Après combien de rebonds la hauteur maximale est-elle inférieure à 1 mm ?`),
      ],
      solution: [
        p(`(a) $t_0 = sqrt{2 h_0/g} = sqrt{2 × 1,6/9,81} = 0,571$ s ; $t_{total} = 0,571 × (1,5/0,5)$ = **1,71 s**.`),
        p(`(b) On veut $1,6 × 0,25^n < 10^{−3}$, donc $n > ln(1600)/ln 4 = 5,32$ : **6 rebonds** suffisent.`),
      ],
    },
    {
      level: 4,
      title: "Un bloc lancé vers le haut",
      statement: [
        p(`Un bloc est lancé à 4 m/s vers le haut d’un plan incliné de 30° (μ = 0,2, g = 9,81 m/s²). Calculez la décélération, la distance parcourue avant l’arrêt et la durée de la montée.`),
      ],
      solution: [
        p(`$a = g (sin 30° + 0,2 × cos 30°) = 9,81 × (0,5 + 0,1732)$ = **6,60 m/s²**.`),
        p(`Distance : $d = v_0²/(2a) = 16/13,21$ = **1,21 m**. Durée : $v_0/a = 4/6,60$ = **0,61 s**.`),
      ],
    },
    {
      level: 5,
      title: "Un choc avec perte d’énergie",
      statement: [
        p(`Une bille de 2 kg à 3 m/s heurte une bille de 1 kg au repos (e = 0,5). Calculez les vitesses finales, vérifiez la conservation de l’impulsion, et calculez l’énergie perdue.`),
      ],
      solution: [
        p(`$v_1' = ((2 − 0,5 × 1) × 3)/3$ = **1,5 m/s** ; $v_2' = ((1 + e) m_1 v_1)/(m_1 + m_2) = 1,5 × 2 × 3 / 3$ = **3,0 m/s**.`),
        p(`Impulsion : avant $2 × 3 = 6$ ; après $2 × 1,5 + 1 × 3 = 6$ : conservée. Énergie : avant $½ × 2 × 9 = 9$ J, après $½ × 2 × 2,25 + ½ × 1 × 9 = 6,75$ J : **2,25 J perdus** (25 %). Vérification par la formule : $½ × (2/3) × (1 − 0,25) × 3² = 2,25$ J.`),
      ],
    },
    {
      level: 5,
      title: "Le contact de Hertz",
      statement: [
        p(`Deux billes égales (m = 1, k = 10⁴ en unités normalisées) se heurtent avec la vitesse relative v. Calculez l’écrasement maximal et la durée du contact pour v = 1, puis pour v = 2. Que remarquez-vous sur la durée ? On donne $T = 2,943275 δ_{max}/v$ et $m_r = 1/2$.`),
      ],
      solution: [
        p(`$v = 1$ : $δ_{max} = (5 × 0,5 × 1/(4 × 10^4))^{2/5} = (6,25 × 10^{−5})^{0,4}$ = **0,020814** ; $T = 2,943275 × 0,020814$ = **0,061261**.`),
        p(`$v = 2$ : $δ_{max} = 2^{4/5} × 0,020814$ = **0,036239** ; $T = 2,943275 × 0,036239/2$ = **0,053331**.`),
        p(`En doublant la vitesse la durée ne baisse que de 13 % (facteur $2^{−1/5} = 0,87$) : le contact est **presque indépendant de la vitesse**, comme un ressort. Pour un vrai ressort linéaire il serait exactement indépendant.`),
      ],
    },
    {
      level: 6,
      title: "Événement et amortissement",
      statement: [
        p(`(a) Combien d’itérations de bissection faut-il pour localiser un événement à 10⁻¹⁰ s près à partir d’un pas de 0,05 s ? (b) Avec l’amortissement de Hunt-Crossley, e ≈ 1/(1 + α v) : calculez e et la perte relative d’énergie pour α = 0,1 et v = 1. Comparez avec l’approximation 2αv.`),
      ],
      solution: [
        p(`(a) Chaque itération divise l’intervalle par 2 : $log_{2}(0,05/10^{−10}) = 28,9$, arrondi à l’entier supérieur : **29 itérations**.`),
        p(`(b) $e = 1/(1 + 0,1)$ = **0,909** (mesuré 0,909016). Énergie après / avant = $e² = 0,826$ : perte **17,4 %**. L’approximation $2αv = 20 %$ est du bon ordre (elle sert pour $α v ≪ 1$ : à $α = 0,01$ elle donne 2 % pour une perte réelle de 1,97 %).`),
      ],
    },
    {
      level: 6,
      title: "La balle qui ne s’arrête jamais",
      statement: [
        p(`Avec le modèle naïf, la vitesse verticale résiduelle observée 0,5 s après le repos exact vaut 0,78 / 0,28 / 0,11 m/s pour Δt = 0,02 / 0,01 / 0,005 s. (a) Calculez le rapport de chaque valeur à g Δt. (b) Que dire du comportement quand Δt → 0 ? (c) Pourquoi l’événement règle-t-il le problème ?`),
      ],
      solution: [
        p(`(a) $g Δt$ = 0,196 / 0,0981 / 0,0491 m/s ; rapports **3,98 / 2,85 / 2,24**.`),
        p(`(b) La vitesse résiduelle est **proportionnelle au pas** (de l’ordre de 2 à 4 fois $g Δt$) : elle tend vers 0 quand $Δt → 0$, mais pour tout pas fini la balle tremble indéfiniment.`),
        p(`(c) Avec l’événement, on détecte exactement l’instant où la hauteur s’annule, on applique le coefficient de restitution à cet instant, et on arrête le calcul quand la vitesse de rebond passe sous un seuil : il n’y a plus de rebonds parasites.`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Frottement sec : le bloc glisse si $tan θ > μ$ ; sinon il adhère. Choc : $v' = e v$ sur la composante normale, $h_n = e^{2n} h_0$, durée totale des rebonds $t_0 (1+e)/(1−e)$.",
      "Une discontinuité à l’intérieur d’un pas ramène **tous** les schémas à l’ordre 1 : il faut **détecter l’événement** (bissection) ; un solveur adaptatif doit être **borné** en nombre de pas.",
      "Le contact de **Hertz** ($F = k δ^{3/2}$) évite l’événement mais donne un autre ordre de convergence (RK4 : 2,3 à 2,75), et dans le berceau de Newton les billes de tête **reculent**, contrairement au modèle d’impulsions.",
      "Une erreur de **modèle** (ignorer l’adhérence) ne disparaît pas quand on raccourcit le pas.",
    ]),
  ]),
};
