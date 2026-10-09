const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 6,
  title: "Les solides qui tournent",
  lead: "Lancez un livre en l’air en le faisant tourner : selon l’axe, il tourne bien droit ou fait un demi-tour surprenant. Une toupie, elle, semble défier la pesanteur. La rotation d’un solide réserve des surprises qu’aucun point matériel ne connaît.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** le moment d’inertie, les trois axes d’un solide dont un est instable, la précession d’une toupie, et pourquoi représenter une orientation demande un objet particulier, le quaternion.",
      "**Les résultats marquants :** l’axe « du milieu » croît bien comme la théorie le prédit (0,5774 mesuré à 0,5774) ; une toupie sous sa vitesse critique **ne tombe pas** ; un schéma qui garde la norme du quaternion peut quand même perdre le moment cinétique.",
      "**Dans le logiciel :** module M6, simulation 10 ({{--sim 10}}).",
    ]),
  ]),
  levels: {
    1: [
      p(`Prenez un livre fermé par un élastique et lancez-le en l’air en le faisant tourner. Il a trois axes de rotation possibles : autour de sa hauteur, de sa largeur, de son épaisseur. Pour deux d’entre eux, le livre tourne bien droit. Pour le troisième, il fait un demi-tour surprenant au milieu de la rotation. Un cosmonaute l’a observé en 1985 en apesanteur avec un écrou à ailettes (effet Djanibekov) ; une raquette de tennis lancée en tournant fait la même chose.`),
      p(`Une **toupie** lancée vite reste debout, alors qu’à l’arrêt elle tomberait. Elle ne tombe pas : elle se met à tourner lentement autour de la verticale, comme une aiguille qui décrirait un cône. C’est la **précession**. Le patineur qui rentre les bras tourne plus vite : la rotation « se souvient » de son élan.`),
      p(`Pour l’ordinateur, il faut en plus retenir l’orientation du solide à chaque instant, avec quatre nombres appelés un **quaternion**. Si la méthode de calcul est mal choisie, ces nombres « gonflent » et le solide lui-même grossit à l’écran : c’est visible dans le logiciel, on y voit le cylindre d’Euler enfler.`),
      fig("fig_axes.png", "Un solide libre dont les trois moments d’inertie sont différents (1, 2 et 3), mis en rotation autour chacun de ses trois axes avec une très petite perturbation. À gauche et à droite, la rotation reste régulière ; au milieu, la perturbation grossit et le solide « bascule ».", 1.0),
    ],
    2: [
      p(`Le **moment d’inertie** mesure la résistance d’un solide à la mise en rotation, comme la masse mesure la résistance à la mise en mouvement. Il dépend de la répartition de la masse **et de l’axe** : un disque est plus facile à faire tourner autour de son axe qu’autour d’un diamètre. Tout solide a trois axes particuliers (les **axes principaux**) avec trois moments d’inertie : un petit, un moyen, un grand.`),
      p(`La rotation autour de l’axe du **plus petit** ou du **plus grand** moment d’inertie est stable. Celle autour de l’axe **intermédiaire** est instable : une perturbation infime grossit exponentiellement et le solide bascule. C’est ce que montre la figure du niveau 1.`),
      p(`En l’absence de force extérieure, la rotation conserve le **moment cinétique** : un patineur qui rentre les bras (moment d’inertie plus petit) tourne plus vite. Pour une **toupie**, la pesanteur exerce un couple qui ne la fait pas tomber mais précesser. Une toupie verticale « endormie » n’est stable que si elle tourne assez vite : sous une **vitesse critique**, la moindre inclinaison grandit ; mais, contre l’intuition, l’inclinaison finit par se refermer au lieu de se terminer par une chute (système sans frottement fixé au pivot).`),
    ],
    3: [
      p(`Le moment d’inertie d’un disque plein de masse $M$ et de rayon $R$ qui tourne autour de son axe est $I = ½ M R²$ ; celui d’une sphère pleine, $I = (2/5) M R²$ (en kg·m²). Pour une vitesse de rotation $ω$ en rad/s (1 tour = 2π rad), le **moment cinétique** et l’**énergie** valent :`),
      eq(`L = I ω, E = frac{1}{2} I ω²`),
      box("exemple", "Exemple : un disque et une sphère", [
        p(`Disque plein de 2 kg et de rayon 0,1 m : I = ½ × 2 × 0,1² = **0,010 kg·m²** ; sphère pleine de mêmes masse et rayon : I = 0,4 × 2 × 0,1² = **0,008 kg·m²**. À ω = 20 rad/s : pour le disque L = 0,2 kg·m²/s et E = **2 J**.`),
      ]),
      box("exemple", "Exemple : le patineur", [
        p(`Un patineur de moment d’inertie I₁ = 3 kg·m² tourne à 2 tours par seconde (ω₁ = 4π = 12,57 rad/s). Il rentre les bras : I₂ = 1,5 kg·m². Le moment cinétique se conserve : I₁ ω₁ = I₂ ω₂, donc ω₂ = 2 ω₁ : il tourne à **4 tours par seconde**. Son énergie passe de ½ × 3 × 12,57² = **237 J** à **474 J** : elle double, grâce au travail de ses muscles (la force centripète des bras).`),
      ]),
    ],
    4: [
      p(`La rotation obéit à une loi analogue à la deuxième loi de Newton : la variation du moment cinétique est égale au **couple** (moment des forces) $τ = r × F$ :`),
      eq(`frac{dL}{dt} = τ`),
      p(`Pour une **toupie** de masse $m$ dont le centre de masse est à la distance $l$ du pivot, le poids exerce un couple horizontal de module $m g l sin θ$. Quand la toupie tourne vite autour de son axe (moment d’inertie $I_3$, vitesse $ω_3$), l’axe tourne autour de la verticale à la vitesse de précession :`),
      eq(`Ω ≈ frac{m g l}{I_3 ω_3}`),
      box("exemple", "Exemples", [
        ul([
          "Toupie du logiciel : m = 1 kg, l = 0,5 m, I₃ = 0,4 kg·m², ω₃ = 25 rad/s : Ω ≈ 9,81 × 0,5 / (0,4 × 25) = **0,49 rad/s**, soit un tour en 12,8 s. La valeur exacte pour θ = 0,6 rad est 0,517 rad/s (un tour en 12,2 s).",
          "Toupie jouet : m = 0,5 kg, l = 4 cm, I₃ = 2 × 10⁻⁴ kg·m², ω₃ = 300 rad/s : Ω = 0,5 × 9,81 × 0,04 / (2 × 10⁻⁴ × 300) = **3,27 rad/s**, soit un tour en 1,9 s.",
        ]),
      ]),
      p(`Deux grandeurs se conservent pour une toupie à pivot fixe : l’énergie, et les composantes du moment cinétique selon la verticale et selon l’axe de la toupie. Dans le logiciel (niveau 4), les invariants sont tracés pour les quatre méthodes de calcul.`),
    ],
    5: [
      p(`**Équations d’Euler.** Dans le repère lié au solide, avec les moments d’inertie principaux $I_1, I_2, I_3$ et la rotation $ω = (ω_1, ω_2, ω_3)$, sans couple :`),
      eq(`I_1 ω_1' = (I_2 − I_3) ω_2 ω_3`),
      eq(`I_2 ω_2' = (I_3 − I_1) ω_3 ω_1`),
      eq(`I_3 ω_3' = (I_1 − I_2) ω_1 ω_2`),
      p(`L’énergie $E = ½ Σ I_i ω_i²$ et la norme du moment cinétique $|L|² = Σ (I_i ω_i)²$ sont conservées. **Solide symétrique** ($I_1 = I_2$) : $ω_3$ est constante et ($ω_1, ω_2$) tourne dans le repère du corps à la vitesse $Ω = (I_3 − I_1) ω_3 / I_1$. Pour $I_1 = I_2 = 1$, $I_3 = 0,5$ et $ω_3 = 3$ : $Ω = −1,5$ rad/s.`),
      p(`**Axe intermédiaire.** Linéarisons autour de $ω = (0, ω_2, 0)$ avec $I_1 < I_2 < I_3$ : $ω_1'' = λ² ω_1$ avec :`),
      eq(`λ = ω_2 sqrt{frac{(I_2 − I_1)(I_3 − I_2)}{I_1 I_3}}`),
      p(`La perturbation croît comme $e^{λ t}$. Pour $I = (1, 2, 3)$ et $ω_2 = 1$ : $λ = 1/sqrt{3} = 0,5774$ ; le logiciel mesure 0,577411 (écart relatif $10^{−4}$) et la figure du niveau 1, recalculée avec un RK4 indépendant, 0,584. Autour des axes extrêmes, la linéarisation donne des oscillations stables (perturbation de $1,4 × 10^{−3}$ qui reste inférieure à $2 × 10^{−3}$ en 100 s).`),
      p(`**Toupie de Lagrange.** Avec $u = cos θ$, la conservation de l’énergie et des deux moments cinétiques ramène la nutation à une cubique $u'² = f(u)$. La toupie verticale est stable si $ω_3$ dépasse la vitesse critique :`),
      eq(`ω_c = frac{2 sqrt{I_1 m g l}}{I_3}`),
      p(`Pour $m = 1$, $I = (1,2 ; 1,2 ; 0,4)$ et $l = 0,5$ : $ω_c = 12,128$ rad/s. Au-dessus (1,15 / 1,5 / 3 fois) une inclinaison initiale de 0,01 reste inférieure à 0,0203 / 0,0134 / 0,0106. À 0,7 $ω_c$ elle croît comme $e^{γ t}$ avec $γ$ prévu 1,44358 et mesuré 1,44357.`),
    ],
    6: [
      p(`**Orientation.** L’état est $(q, ω)$ où le quaternion unitaire $q$ envoie le repère du corps sur le repère fixe et $ω$ est exprimée dans le corps ; $q' = ½ q ⊗ (0, ω)$. La solution exacte du solide symétrique est $q(t) = q_L(L t/I_1) q_0 q_3(−Ω t)$ ; elle coïncide avec un RK45 serré à $10^{−14}$ … $2 × 10^{−12}$ sur $ω$ et $7 × 10^{−14}$ … $4 × 10^{−12}$ sur l’orientation ; avec le signe de la rotation propre faux le même test donne 0,56 à 2,8 (le test est sensible). Le solide asymétrique se résout par les fonctions elliptiques de Jacobi (Landau et Lifchitz, § 37) : pour $(a, c)$ = (1,5 ; 1), (2 ; 0,5), (1,5 ; 0,8) le paramètre vaut $m$ = 0,75 / 0,1875 / 0,8533 et la période $4K/λ$ = 8,626 / 5,727 / 11,082 s, avec un écart à RK45 inférieur à $5,7 × 10^{−13}$ pour $t ≤ 40$ s.`),
      p(`**Quatre intégrateurs d’orientation** (corps $I = (1, 2, 3)$, $ω_0 = (0,9 ; 0,5 ; 1,1)$, $t = 6,3$) : Euler ordre 1 (rapport d’erreur 2,01 par doublement des pas) ; RK4 avec renormalisation du quaternion, ordre 4 (16,06) ; Lie-Heun (méthode de groupe de Lie, norme exacte), ordre 2 (3,99) ; découpage symplectique, ordre 2 (4,00) et environ trois fois plus précis que Heun à pas égal.`),
      tbl(
        ["Schéma (1000 s, pas 0,05)", "|q| − 1", "Énergie", "Moment cinétique"],
        [
          ["Découpage symplectique", "3 × 10⁻¹⁶", "max 6,9 × 10⁻⁵, sans dérive", "6,2 × 10⁻¹⁴"],
          ["RK4 + renormalisation", "3 × 10⁻¹⁶", "−3,1 × 10⁻⁶ à la fin", "1,8 × 10⁻⁶"],
          ["RK4 sans renormalisation", "4 × 10⁻⁶", "—", "1,7 × 10⁻⁵"],
          ["Lie-Heun", "à l’arrondi", "+1,4 × 10⁻² (dérive)", "7 × 10⁻³ (dérive)"],
          ["Euler explicite", "diverge au pas 2251", "—", "—"],
        ],
        [2900, 1700, 2500, 2198],
        { align: ["l", "c", "c", "c"], size: 19 }
      ),
      p(`L’énergie du découpage symplectique évolue en $h²$ ($6,9 × 10^{−5}$ à $h = 0,05$, 2,74 × 10⁻⁴ à 0,1, 1,07 × 10⁻³ à 0,2, soit des rapports 3,98 puis 3,9) et son maximum est le même dans chaque tiers du calcul : **aucune dérive**. Pour Euler explicite sur une rotation pure, la norme croît exactement comme $|q|_N = (1 + h² ω²/4)^{N/2}$.`),
      box("attention", "Attention : garder la norme ne suffit pas", [
        p(`Lie-Heun conserve la norme du quaternion à l’arrondi près, et pourtant il dérive en énergie et en moment cinétique : c’est la **structure symplectique** du découpage qui conserve le moment cinétique, non le simple respect de la contrainte $|q| = 1$.`),
      ]),
      p(`**Toupie, nutation.** Pour $θ_0 = 0,9$ et $ω_0 = (0, 2, 20)$, les racines de la cubique valent $u_1 = 0,621609968 = cos 0,9$ et $u_2 = 0,926152996$ ; la période de nutation $T = 4K/sqrt{β (u_3 − u_1)}$ vaut 1,069185 s et $u(T) − u_0 = 1,6 × 10^{−14}$. Avec un couple, le découpage conserve $L_z$ à $2,8 × 10^{−13}$ (le couple n’a pas de composante verticale) et $E$ à $3,9 × 10^{−5}$ sur 100 s.`),
      box("attention", "Attention : une toupie sous la vitesse critique ne tombe pas", [
        p(`Hypothèse de départ : « sous $ω_c$ la toupie tombe ». La mesure la contredit : de $10^{−6}$ rad elle monte jusqu’à 1,2457 rad (71°), puis **remonte** (θ ≈ $10^{−5}$ à 20 s) : grande nutation d’un système sans perte fixé à son pivot. Pour qu’elle tombe réellement il faudrait modéliser une table (contact, frottement), ce que le logiciel ne fait pas.`),
      ]),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Lancez {{physicslab --level 5 --sim 10}} : trois scénarios (**Solide symétrique libre**, **Solide asymétrique libre**, **Toupie pesante**) et des cas types (**Cigare**, **Sphère**, **Disque**).`),
      ul([
        "Solide asymétrique : lancez la rotation autour de l’axe du milieu avec une très petite perturbation (**Écart de l’axe**) et regardez la **polhodie** (le graphe de ω₂ en fonction de ω₁) ; essayez les deux autres axes.",
        "Toupie : réglez le **spin ω₃** ; la fenêtre **Théorie** indique la vitesse de précession exacte et la vitesse critique ; passez sous la vitesse critique et observez la nutation.",
        "Comparez les cinq schémas d’orientation : le cylindre d’**Euler** enfle, **Lie-Heun** garde sa taille mais dérive en énergie, le **découpage symplectique** conserve le moment cinétique.",
      ]),
    ]),
    fig("shot_s10.png", "Le solide libre symétrique au niveau 5 : les enveloppes dessinées par chaque méthode, la théorie à droite, le moment cinétique et la norme du quaternion en bas.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "Le livre qui bascule",
      statement: [
        p(`Pour un livre lancé en tournant, lequel des trois axes (hauteur, largeur, épaisseur) est instable ? Quelle est la caractéristique de cet axe ? Pourquoi le patineur tourne-t-il plus vite en rentrant les bras ?`),
      ],
      solution: [
        p(`L’axe instable est celui du **moment d’inertie intermédiaire** : pour un livre c’est l’axe de la largeur (ni le plus long, ni le plus court). Les deux autres axes sont stables.`),
        p(`Le patineur conserve son moment cinétique $I ω$. En rentrant les bras, il diminue $I$ ; $ω$ doit donc augmenter.`),
      ],
    },
    {
      level: 3,
      title: "Moments d’inertie",
      statement: [
        p(`Un disque plein et une sphère pleine ont la même masse 2 kg et le même rayon 0,1 m. Calculez leurs moments d’inertie, puis le moment cinétique et l’énergie du disque tournant à 20 rad/s.`),
      ],
      solution: [
        p(`Disque : $I = ½ M R² = ½ × 2 × 0,01$ = **0,010 kg·m²**. Sphère : $I = (2/5) M R²$ = **0,008 kg·m²**.`),
        p(`Disque à 20 rad/s : $L = I ω$ = **0,2 kg·m²/s** et $E = ½ I ω² = ½ × 0,01 × 400$ = **2 J**.`),
      ],
    },
    {
      level: 3,
      title: "Le patineur",
      statement: [
        p(`Un patineur tourne à 2 tours par seconde avec I₁ = 3 kg·m². Il rentre les bras et son moment d’inertie devient I₂ = 1,5 kg·m². Quelle est sa nouvelle vitesse en tours par seconde ? Comparez les énergies avant et après. D’où vient la différence ?`),
      ],
      solution: [
        p(`$I_1 ω_1 = I_2 ω_2$ donne $ω_2 = 3/1,5 × ω_1 = 2 ω_1$ : **4 tours par seconde**.`),
        p(`Énergies : avant $½ × 3 × (4π)²$ = **237 J** ; après $½ × 1,5 × (8π)²$ = **474 J**. L’énergie a **doublé** : la différence vient du travail fourni par le patineur pour ramener ses bras vers l’axe, contre la force centrifuge.`),
      ],
    },
    {
      level: 4,
      title: "Précession d’une toupie",
      statement: [
        p(`La toupie du logiciel a m = 1 kg, l = 0,5 m, I₃ = 0,4 kg·m² et ω₃ = 25 rad/s. Calculez la vitesse de précession approchée Ω = m g l/(I₃ ω₃) et la durée d’un tour de précession. La valeur exacte pour θ = 0,6 rad est 0,517 rad/s : quel est l’écart relatif ?`),
      ],
      solution: [
        p(`$Ω = 1 × 9,81 × 0,5/(0,4 × 25) = 4,905/10$ = **0,49 rad/s** ; un tour dure $2π/Ω$ = **12,8 s**.`),
        p(`Écart avec 0,517 : $(0,517 − 0,490)/0,517 = 5,2 %$. L’approximation « toupie rapide » néglige un terme en $I_1 Ω cos θ$ devant $I_3 ω_3$.`),
      ],
    },
    {
      level: 4,
      title: "Une toupie jouet",
      statement: [
        p(`Une toupie de m = 0,5 kg a son centre de masse à 4 cm du pivot, un moment d’inertie I₃ = 2 × 10⁻⁴ kg·m² et tourne à ω₃ = 300 rad/s. Calculez sa vitesse de précession et le nombre de tours de précession par minute.`),
      ],
      solution: [
        p(`$Ω = m g l/(I_3 ω_3) = 0,5 × 9,81 × 0,04/(2 × 10^{−4} × 300) = 0,1962/0,06$ = **3,27 rad/s**.`),
        p(`Nombre de tours par minute : $3,27/(2π) × 60$ ≈ **31 tours/min** (un tour en 1,9 s).`),
      ],
    },
    {
      level: 5,
      title: "L’axe intermédiaire",
      statement: [
        p(`Pour un solide de moments d’inertie (1, 2, 3) qui tourne autour de l’axe 2 à ω₂, calculez le taux de croissance λ d’une perturbation pour ω₂ = 1 puis ω₂ = 2, le facteur d’amplification en 3 s, et le temps de doublement.`),
      ],
      solution: [
        p(`$λ = ω_2 sqrt{(I_2 − I_1)(I_3 − I_2)/(I_1 I_3)} = ω_2 sqrt{1 × 1/(1 × 3)} = ω_2 × 0,5774$ : **0,577 s⁻¹** pour $ω_2 = 1$ et **1,155 s⁻¹** pour $ω_2 = 2$.`),
        p(`Amplification en 3 s : $e^{1,732}$ = **5,65** et $e^{3,464}$ = **31,9**. Temps de doublement $ln 2/λ$ : **1,20 s** et **0,60 s**.`),
        p(`Plus le solide tourne vite, plus le basculement est rapide : à vitesse double, il est deux fois plus rapide.`),
      ],
    },
    {
      level: 5,
      title: "Corps libre symétrique et vitesse critique",
      statement: [
        p(`(a) Un solide libre a I₁ = I₂ = 1, I₃ = 0,5 et ω₃ = 3 rad/s. Calculez la vitesse Ω de rotation de (ω₁, ω₂) dans le repère du corps et sa période. (b) La toupie du logiciel (I₁ = 1,2, I₃ = 0,4, m = 1, l = 0,5) a ω_c = 2√(I₁ m g l)/I₃. Calculez ω_c ; la toupie lancée à 25 rad/s est-elle stable ? À 0,7 ω_c, que vaut la croissance, avec γ = 1,4436 s⁻¹ : temps de doublement ?`),
      ],
      solution: [
        p(`(a) $Ω = (I_3 − I_1) ω_3/I_1 = (0,5 − 1) × 3/1$ = **−1,5 rad/s** (le signe indique le sens de rotation) ; période $2π/1,5$ = **4,19 s**.`),
        p(`(b) $ω_c = 2 sqrt{1,2 × 9,81 × 0,5}/0,4 = 2 × 2,426/0,4$ = **12,13 rad/s**. À 25 rad/s, soit 2,06 $ω_c$, la toupie est **stable**. À 0,7 $ω_c = 8,49$ rad/s, une inclinaison double en $ln 2/γ$ = **0,48 s**.`),
      ],
    },
    {
      level: 6,
      title: "La norme du quaternion avec Euler",
      statement: [
        p(`Euler explicite fait croître la norme du quaternion comme $|q|_N = (1 + h² ω²/4)^{N/2}$. (a) Pour ω = 2 rad/s et h = 0,05 s, calculez |q| après N = 100 puis N = 2000 pas. (b) Quel pas faut-il pour que |q| − 1 reste inférieur à 10⁻⁶ sur 100 pas ?`),
      ],
      solution: [
        p(`(a) $1 + h² ω²/4 = 1 + 0,0025 = 1,0025$. $N = 100$ : $1,0025^{50}$ = **1,133** (le solide a grossi de 13 %). $N = 2000$ : $1,0025^{1000}$ = **12,1** : le solide est devenu douze fois plus grand.`),
        p(`(b) On veut $(1 + x)^{50} − 1 < 10^{−6}$, donc $x ≈ 2 × 10^{−8}$ avec $x = h² ω²/4 = h²$, soit $h$ ≈ **1,4 × 10⁻⁴ s**. Il faudrait 350 fois plus de pas : une renormalisation ou un meilleur schéma est bien plus rentable.`),
      ],
    },
    {
      level: 6,
      title: "Quel schéma pour mille secondes ?",
      statement: [
        p(`Sur 1000 s à h = 0,05 s, le découpage symplectique garde ΔL/L = 6,2 × 10⁻¹⁴, RK4 avec renormalisation 1,8 × 10⁻⁶ et Lie-Heun 7 × 10⁻³. Calculez les rapports entre eux. Quel schéma choisir pour (a) une animation de quelques secondes, (b) une simulation de plusieurs heures ? Justifiez.`),
      ],
      solution: [
        p(`Rapports : RK4 / découpage = $1,8 × 10^{−6}/6,2 × 10^{−14} ≈ 3 × 10^{7}$ ; Lie-Heun / RK4 = $7 × 10^{−3}/1,8 × 10^{−6} ≈ 4 × 10^{3}$.`),
        p(`(a) Pour quelques secondes, **RK4 avec renormalisation** : ordre 4 donc très précis à court terme, la dérive est invisible. (b) Pour de longues durées, le **découpage symplectique** : sa précision à pas égal est moindre (ordre 2) mais le moment cinétique reste exact à l’arrondi et l’énergie ne dérive pas ; les erreurs des autres s’accumulent. **Lie-Heun** est à éviter malgré sa norme exacte.`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Le moment d’inertie dépend de l’axe ; $L = I ω$ et $E = ½ I ω²$. Un solide libre est stable autour de ses axes extrêmes, **instable autour de l’axe intermédiaire** (croissance en $e^{λ t}$).",
      "Une toupie rapide précesse à $Ω ≈ m g l/(I_3 ω_3)$ ; la verticale est stable si $ω_3 > ω_c = 2 sqrt{I_1 m g l}/I_3$ ; sous $ω_c$ elle **ne tombe pas** (grande nutation sans perte).",
      "L’orientation se représente par un **quaternion unitaire** ; garder $|q| = 1$ ne suffit pas à conserver le moment cinétique : un schéma **symplectique** le fait.",
      "Euler explicite fait enfler le solide : $|q|_N = (1 + h² ω²/4)^{N/2}$.",
    ]),
  ]),
};
