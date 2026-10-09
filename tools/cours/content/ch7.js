const B = require("../lib/blocks");
const { p, eq, ul, ol, box, fig, tbl } = B;

module.exports = {
  num: 7,
  title: "Mille étoiles sur une carte graphique",
  lead: "Chaque étoile attire toutes les autres : doubler le nombre d’étoiles quadruple le travail. Pour suivre des milliers d’étoiles en temps réel, on confie ce travail à des milliers de petits calculateurs qui avancent ensemble : ceux de la carte graphique.",
  enBref: box("enbref", "En bref", [
    ul([
      "**Ce que vous allez comprendre :** pourquoi un amas d’étoiles coûte cher à calculer, comment une carte graphique le fait en parallèle, et ce que l’on perd en précision (7 chiffres au lieu de 16).",
      "**Les résultats marquants :** la carte graphique est **80 à 150 fois plus rapide** que le processeur sur 16 000 étoiles ; son erreur croît comme $6 × 10^{−8} sqrt{N}$ ; l’énergie dérive autant que l’erreur du schéma lui-même.",
      "**Dans le logiciel :** module M7, simulation 11 ({{--sim 11}}) et test {{--gpu-test}}.",
    ]),
  ]),
  levels: {
    1: [
      p(`Imaginez une cuisine qui doit préparer deux millions de plats. Avec un seul grand chef, qui fait tout très bien et très vite, cela prend un certain temps. Avec mille petits commis, chacun moins brillant mais tous travaillant en même temps, c’est fini bien plus tôt. Un **processeur** d’ordinateur est le grand chef. Une **carte graphique**, celle qui dessine les jeux vidéo, est la brigade de mille commis.`),
      p(`Notre problème tombe à pic : avec 2000 étoiles, chacune attire toutes les autres, soit près de **deux millions de paires** à calculer, et il faut tout recommencer à chaque instant. Si l’on double le nombre d’étoiles, il y a quatre fois plus de paires. Pour ce genre de travail, la brigade gagne largement.`),
      p(`Il y a un prix : la carte graphique calcule avec des nombres moins précis (7 chiffres au lieu de 16). Pour dessiner un amas c’est largement suffisant, mais deux calculs de précisions différentes finissent par se séparer, comme les deux double-pendules du chapitre 3. Alors pour **vérifier** la carte graphique, on la fait calculer et l’on compare avec le processeur : « le processeur vérifie la carte graphique ».`),
      fig("shot_s11_l5.png", "Deux mille étoiles calculées en même temps par le processeur (en bleu clair) et par la carte graphique (en orange) : au début on ne voit qu’un seul nuage. À droite, les écarts mesurés.", 1.0),
    ],
    2: [
      p(`Le calcul des forces entre toutes les étoiles s’appelle un calcul « **N corps** » : pour $N$ étoiles il y a $N(N−1)/2$ paires, donc un travail qui croît comme le **carré** de $N$. C’est ce qui fait qu’on ne peut pas simplement acheter un ordinateur dix fois plus rapide pour traiter dix fois plus d’étoiles : il en faudrait cent fois plus.`),
      p(`La carte graphique contient des milliers de petits calculateurs qui exécutent **le même programme** en même temps, chacun sur une étoile différente. C’est le **parallélisme**. Le logiciel mesure, avec 16 000 étoiles, une accélération de 80 à 150 fois par rapport au processeur de la même machine, et c’est ce qui permet de voir bouger des milliers d’étoiles à l’écran.`),
      p(`Le prix à payer est la précision : en **double précision** (processeur) un nombre a environ 16 chiffres significatifs, en **simple précision** (carte graphique) environ 7. On vérifie donc deux choses : (1) à un instant donné, les forces calculées par la carte sont-elles les mêmes que celles du processeur ? (2) sur quelques unités de temps, l’énergie, l’impulsion et le moment cinétique restent-ils constants ? Les réponses : oui à environ un millionième près, et oui. Au-delà de quelques unités de temps, les trajectoires se séparent, comme dans tout système chaotique.`),
      p(`Le logiciel simule un **amas de Plummer** (un modèle réaliste d’amas d’étoiles à l’équilibre) ou la **collision de deux amas**. Si le calcul est trop lourd pour fonctionner en temps réel, il ralentit le temps plutôt que de bloquer l’écran.`),
    ],
    3: [
      p(`Le nombre de paires d’étoiles est $N(N−1)/2$. Voici quelques valeurs :`),
      tbl(
        ["Nombre d’étoiles N", "Paires N(N−1)/2", "Interactions par pas (N²)"],
        [
          ["500", "124 750", "250 000"],
          ["2 000", "1 999 000", "4 millions"],
          ["10 000", "49 995 000", "100 millions"],
          ["100 000", "5 milliards", "10 milliards"],
        ],
        [2800, 3000, 3498],
        { align: ["c", "c", "c"] }
      ),
      p(`Le logiciel mesure sur la carte graphique de l’ordinateur de développement (une carte intégrée, modeste) environ **30 milliards d’interactions par seconde**, et sur le processeur environ **0,33 milliard** par seconde. Avec 100 000 étoiles il y a 10 milliards d’interactions par pas :`),
      ul([
        "carte graphique : 10 / 30 = **0,33 s** par pas (mesuré : 0,34 s) ;",
        "processeur : 10 / 0,33 = **30 s** par pas ;",
        "rapport : environ **90 fois**.",
      ]),
      box("exemple", "Combien de chiffres ?", [
        p(`Un nombre à 7 chiffres significatifs comme 0,1234567 n’a pas la précision de 0,12345678901234567 (16 chiffres). Quand on additionne deux millions de forces, chaque addition ajoute une toute petite erreur d’arrondi ; elles ne s’annulent pas toutes. Le niveau 5 explique à combien elles se montent.`),
      ]),
    ],
    4: [
      p(`**La force.** L’accélération de l’étoile $i$ est la somme des attractions de toutes les autres, avec un petit adoucissement $ε$ pour qu’elle ne soit jamais infinie :`),
      eq(`a_i = G sum_{j≠i} frac{m_j (r_j − r_i)}{(|r_j − r_i|² + ε²)^{3/2}}`),
      p(`**Le schéma « saute-mouton »** (Verlet, vu aux chapitres 2 et 4) avance en quatre petites opérations par pas $Δt$ :`),
      ol([
        "demi-coup sur la vitesse : $v ← v + a Δt/2$ ;",
        "dérive de la position : $x ← x + v Δt$ ;",
        "recalcul de toutes les accélérations (la partie lourde) ;",
        "second demi-coup : $v ← v + a Δt/2$.",
      ]),
      p(`Sur le GPU les positions, les vitesses et les accélérations restent en mémoire de la carte ; seul l’affichage récupère les positions à chaque image.`),
      p(`**Trois grandeurs doivent rester constantes** : l’énergie $E = ½ Σ m v² − G Σ_{i<j} m_i m_j / r_{ij}$, l’impulsion $P = Σ m v$ et le moment cinétique $L = Σ m r × v$. Dans le logiciel, avec 2000 étoiles et la carte graphique, elles varient d’environ $10^{−7}$ ou $10^{−8}$ en valeur relative après quelques secondes de simulation.`),
      p(`**L’amas de Plummer** a une densité $(1 + r²/a²)^{−5/2}$ ; pour une masse totale $M$ et $G = 1$ son énergie vaut $E = −3π/64 × M²/a = −0,147$ pour $a = 1$. À l’équilibre l’énergie cinétique $T$ et l’énergie potentielle $U$ vérifient $2T = −U$ (théorème du viriel), ce que le logiciel mesure à 0,5 % près.`),
    ],
    5: [
      p(`**Précision.** Un nombre en simple précision est arrondi à $u = 2^{−24} = 6 × 10^{−8}$ en valeur relative. La somme de $N$ termes arrondis accumule une erreur qui, si les arrondis sont indépendants, croît comme une marche au hasard :`),
      eq(`Δ_{rel} ≈ u sqrt{N}`),
      p(`(et non comme $u N$, le pire cas où tous les arrondis iraient dans le même sens). Le logiciel mesure, pour l’écart global entre accélérations de la carte et du processeur à un instant fixé, 1,6 × 10⁻⁸ ($N$ = 3), 2,7 × 10⁻⁷ (100), 4,7 × 10⁻⁷ (1000) et 1,05 × 10⁻⁶ (5000) : environ le quart de la borne $u sqrt{N}$.`),
      fig("fig_gpu.png", "À gauche : milliards d’interactions calculées par seconde, selon le nombre d’étoiles, pour la carte graphique et pour le processeur (mesures du logiciel sur une carte intégrée Intel UHD). À droite : erreur de la carte graphique en simple précision.", 1.0),
      p(`**Comment la carte fait ses calculs.** Un groupe de 256 calculateurs prend en charge 256 étoiles cibles. Les étoiles sources sont lues **par tuiles de 256** : chaque tuile est chargée une seule fois dans une mémoire rapide partagée par le groupe, puis tous les calculateurs l’utilisent. Chaque étoile n’est donc lue qu’une fois par groupe au lieu de 256 fois ; pour 2000 étoiles il y a 8 tuiles. Le choix de la taille du groupe compte : pour 16 384 étoiles, 32 calculateurs par groupe donnent 17,7 ms, 256 en donnent 8,35 ms.`),
      p(`**Le piège du décalage.** Un amas centré sur l’origine garde de bonnes précisions ; un amas de mêmes étoiles situé vers (1000, −2000, 500) en perd, car en simple précision un nombre proche de 1000 ne sait plus distinguer des positions séparées de moins de $6 × 10^{−5}$. On **centre** donc les positions sur le barycentre, en double précision, avant de les convertir : l’erreur passe de $1,6 × 10^{−4}$ à $4,7 × 10^{−7}$, soit un gain d’un facteur 345.`),
      tbl(
        ["Cas (GPU float)", "Énergie", "Impulsion", "Moment cinétique"],
        [
          ["Amas 100 étoiles, t = 20 (20 000 pas)", "1,0 × 10⁻⁶", "4,9 × 10⁻⁷", "2,9 × 10⁻⁶"],
          ["Amas 500 étoiles, t = 1", "4,1 × 10⁻⁷", "1,8 × 10⁻⁸", "3,3 × 10⁻⁸"],
          ["Amas 5000 étoiles, t = 1", "2,5 × 10⁻⁷", "5,5 × 10⁻⁹", "4,9 × 10⁻⁹"],
        ],
        [3600, 1900, 1900, 1898],
        { align: ["l", "c", "c", "c"], size: 19 }
      ),
      p(`À titre de comparaison, le processeur en double précision garde l’impulsion et le moment cinétique à $10^{−15}$ ; son énergie dérive autant que celle de la carte (7,9 × 10⁻⁷ contre 1,0 × 10⁻⁶ pour 100 étoiles) : l’énergie est limitée par l’**erreur du schéma** d’ordre 2, non par la simple précision.`),
    ],
    6: [
      p(`**Erreur d’arrondi des forces.** Chaque terme $m_j d / (r² + ε²)^{3/2}$ est évalué avec quelques unités de dernier chiffre ; en sommant $N$ termes, l’erreur absolue est de l’ordre de $u sqrt{N} S_i$ où $S_i = Σ |m_j d/(r²+ε²)^{3/2}|$ est la somme des modules, même quand les termes se compensent dans la somme vectorielle (étoile au centre d’un amas). Les forces par paire ne sont opposées qu’à l’arrondi : $|Σ m_i a_i| / Σ m_i |a_i|$ vaut 1 à 2 × 10⁻⁸ en simple précision, contre 10⁻¹⁷ au processeur. L’impulsion dérive donc de 5 × 10⁻⁹ à 2 × 10⁻⁶ selon $N$ et le nombre de pas (et non $10^{−7}$ par pas comme on l’avait prévu).`),
      p(`**Chaos et comparaison.** Les trajectoires GPU et CPU de l’amas de 100 étoiles ($λ ≈ 1$) s’écartent de $3,9 × 10^{−5}$ (t = 2), $2,7 × 10^{−4}$ (4), $2,9 × 10^{−3}$ (6), $2,5 × 10^{−2}$ (8), 0,16 (10), 1,1 (12), 4,0 (14), puis saturent vers 14 : il n’y a **aucun sens** à comparer des trajectoires longues. Sur le « huit », système régulier, l’écart croît linéairement (4,3 × 10⁻⁵ après une période, 2,0 × 10⁻⁴ après cinq) : le bruit d’arrondi des mises à jour de la vitesse et de la position est amplifié par le cisaillement de l’orbite. On compare donc les **accélérations à état fixé**, puis des durées courtes et des grandeurs d’ensemble.`),
      p(`**Contraintes de la carte.** (1) Une commande qui dure plus de deux secondes déclenche la réinitialisation du pilote de Windows (TDR) : le logiciel découpe le calcul en envois d’environ 50 ms, dimensionnés d’après le débit mesuré. (2) Le calcul en double précision du shader est accepté mais **émulé** par le pilote de cette carte : 3 millions d’interactions par seconde, mille à plusieurs milliers de fois moins que le float, et cent fois moins que le processeur ; il ne sert qu’à valider la logique du programme (accord à $3 × 10^{−16}$). (3) Lire le chronomètre de la carte à chaque envoi coûte 130 à 350 µs ; en les regroupant par lots le temps par pas est tombé de 308 à 89 µs pour 100 étoiles.`),
      p(`**Énergie potentielle.** Le programme de la carte calcule en passant le potentiel $φ_i = −G Σ_{j≠i} m_j / sqrt{r² + ε²}$ de chaque étoile (le terme $m_j/r$ est déjà disponible, un seul ajout par interaction) ; $U = ½ Σ m_i φ_i$ est sommée ensuite en double précision. L’écart au calcul du processeur est de 5 × 10⁻⁹ à 4,3 × 10⁻⁷ en valeur relative.`),
      box("attention", "Attention : un défaut qui n’était pas celui qu’on croyait", [
        p(`Quand on changeait le nombre d’étoiles, les étoiles disparaissaient presque toutes de l’écran alors que le calcul était juste : les mesures d’énergie étaient normales. Après avoir soupçonné le calcul sur carte, puis l’ordre des opérations de l’interface, la dichotomie a désigné le module de dessin : il réallouait sa mémoire de sommets à chaque dessin, ce que le pilote de cette carte supporte mal dès que d’autres tampons existent. Un tampon de flux persistant a réglé le problème : **un calcul correct peut cacher un affichage faux**, et inversement.`),
      ]),
      box("plus", "Pour aller plus loin", [
        p(`Le calcul par paires est en $O(N²)$ ; pour des millions d’étoiles on utilise des méthodes d’arbre (Barnes-Hut, en $O(N log N)$) ou multipolaires rapides, qui sacrifient un peu de précision. Les sommations compensées (Kahan) ou en « double-float » récupèrent des chiffres sans passer en double précision. Ces deux pistes ne sont pas implémentées dans le logiciel.`),
      ]),
    ],
  },
  logiciel: [
    box("logiciel", "", [
      p(`Lancez {{physicslab --level 5 --sim 11}} : un amas de Plummer de 2000 étoiles calculé par la carte graphique.`),
      ul([
        "Choisissez **Les deux, côte à côte** : le processeur est en halo bleu clair, la carte graphique en orange ; le tableau **Invariants** montre ΔE/E₀, |ΔP| et |ΔL| de chacun et le **temps par pas** (la carte graphique est environ 25 fois plus rapide avec un build non optimisé).",
        "Dans l’onglet **Analyse**, le bouton **Comparer maintenant GPU et CPU** trace l’histogramme de l’erreur des accélérations, étoile par étoile ; au niveau 6, le **Banc d’essai** trace le débit et l’erreur en fonction de N.",
        "Changez le **nombre d’étoiles** (jusqu’à 20 000 sur la carte, 4000 sur le processeur) : si le calcul est trop lourd, le logiciel affiche « le temps est ralenti » plutôt que de se bloquer.",
        "Hors interface, la commande {{physicslab --gpu-test}} refait les comparaisons du chapitre et signale toute dérive anormale par un code de sortie non nul.",
      ]),
    ]),
    fig("shot_s11_l1.png", "La simulation 11 au niveau 1 : seulement le choix du calculateur (carte graphique ou processeur), le nombre d’étoiles et le temps de calcul par pas.", 1.0),
  ],
  exercises: [
    {
      level: 1,
      title: "Pourquoi la carte graphique ?",
      statement: [
        p(`(a) Pourquoi une carte graphique est-elle plus rapide qu’un processeur pour calculer l’attraction entre des milliers d’étoiles ? (b) Pourquoi ne gagne-t-on pas de temps avec 3 étoiles seulement ?`),
      ],
      solution: [
        p(`(a) Parce qu’elle contient des milliers de petits calculateurs qui exécutent le même programme en parallèle, chacun pour une étoile différente ; le travail se partage.`),
        p(`(b) Avec 3 étoiles le travail est minuscule : le temps est dominé par le coût de lancer le calcul sur la carte (environ 50 à 90 µs par pas), pas par le calcul lui-même. Le logiciel mesure d’ailleurs que, pour 3 à 100 étoiles, processeur et carte sont du même ordre (60 à 150 µs contre 53 à 89 µs).`),
      ],
    },
    {
      level: 3,
      title: "Compter les paires",
      statement: [
        p(`Calculez le nombre de paires N(N−1)/2 pour N = 500, 2000 et 10 000. Par quel facteur le travail est-il multiplié quand on passe de 2000 à 10 000 étoiles ?`),
      ],
      solution: [
        p(`N = 500 : 500 × 499 / 2 = **124 750**. N = 2000 : 2000 × 1999 / 2 = **1 999 000**. N = 10 000 : 10 000 × 9999 / 2 = **49 995 000**.`),
        p(`De 2000 à 10 000 étoiles (facteur 5), le travail est multiplié par $(5)²$ = **25** (49 995 000 / 1 999 000 = 25,0).`),
      ],
    },
    {
      level: 3,
      title: "Temps de calcul",
      statement: [
        p(`Un pas de simulation de 100 000 étoiles demande 10 milliards d’interactions. La carte graphique en calcule 30 milliards par seconde, le processeur 0,33 milliard. Combien de temps dure un pas dans chaque cas ? Combien de fois la carte est-elle plus rapide ?`),
      ],
      solution: [
        p(`Carte : 10 / 30 = **0,33 s**. Processeur : 10 / 0,33 = **30 s**. Rapport : 30 / 0,33 ≈ **91**.`),
        p(`À 0,33 s par pas la simulation avance de 3 pas par seconde : trop lent pour de l’animation fluide ; à 30 s elle est inutilisable. C’est pourquoi le logiciel limite le nombre d’étoiles du processeur à 4000 et celui de la carte à 20 000.`),
      ],
    },
    {
      level: 4,
      title: "Combien de pas par image ?",
      statement: [
        p(`Pour garder l’écran fluide, le logiciel consacre au plus 14 ms de calcul par image. Un pas de 5000 étoiles prend 1,5 ms sur la carte. Combien de pas peut-il faire par image ? Avec un pas de temps Δt = 0,004 et 60 images par seconde, quel temps simulé s’écoule-t-il par seconde réelle ?`),
      ],
      solution: [
        p(`14 / 1,5 = 9,3 : **9 pas par image**.`),
        p(`Par seconde : 9 × 60 = 540 pas, soit $540 × 0,004$ = **2,16 unités de temps simulé par seconde**. (La traversée d’un amas dure environ 1 à 2 unités de temps : on voit donc bouger l’amas à un rythme confortable.)`),
      ],
    },
    {
      level: 5,
      title: "L’erreur d’une somme",
      statement: [
        p(`Avec u = 2⁻²⁴ = 6 × 10⁻⁸, calculez la borne u√N pour N = 100, 1000, 5000 et comparez aux erreurs mesurées 2,7 × 10⁻⁷, 4,7 × 10⁻⁷, 1,05 × 10⁻⁶. Pour quel N l’erreur atteindrait-elle 10⁻⁴ si la loi se maintenait ?`),
      ],
      solution: [
        p(`$u sqrt{N}$ : N = 100 → **6,0 × 10⁻⁷** ; N = 1000 → **1,9 × 10⁻⁶** ; N = 5000 → **4,2 × 10⁻⁶**. Les erreurs mesurées sont de 2 à 4 fois plus petites : la loi en $sqrt{N}$ est une **borne pessimiste** de la bonne forme.`),
        p(`Pour $u sqrt{N} = 10^{−4}$ : $sqrt{N} = 10^{−4}/(6 × 10^{−8}) = 1678$ et $N$ = **2,8 millions** d’étoiles (en dessous de cette taille, la précision de la carte reste meilleure que 10⁻⁴).`),
      ],
    },
    {
      level: 5,
      title: "Les tuiles",
      statement: [
        p(`Un groupe de 256 calculateurs lit les étoiles sources par tuiles de 256. (a) Combien de tuiles pour N = 2000 et pour N = 100 000 ? (b) Par quel facteur diminue le nombre de lectures de la mémoire principale par rapport à une lecture individuelle de chaque calculateur ?`),
      ],
      solution: [
        p(`(a) 2000/256 = 7,8 soit **8** tuiles ; 100 000/256 = 390,6 soit **391** tuiles (on arrondit à l’entier supérieur).`),
        p(`(b) Sans tuile, chacun des 256 calculateurs lit les $N$ sources : $256 N$ lectures par groupe. Avec les tuiles, le groupe lit chaque source une fois : $N$ lectures. Le facteur est **256**, ce qui explique en grande partie la vitesse du calcul (le résultat reste identique : seule la façon de lire change).`),
      ],
    },
    {
      level: 6,
      title: "L’horizon en simple précision",
      statement: [
        p(`Le système a λ ≈ 1 s⁻¹. En simple précision l’arrondi initial est 6 × 10⁻⁸ ; en double précision 1,1 × 10⁻¹⁶. En acceptant un écart de 0,1, calculez l’horizon de prédiction dans chaque cas. L’écart mesuré entre GPU et CPU atteint 0,1 vers t = 9,5 pour 100 étoiles : pourquoi plus tôt que ln(0,1/6 × 10⁻⁸) ?`),
      ],
      solution: [
        p(`$t_H = ln(0,1/ε)/λ$ : simple précision $ln(1,67 × 10^{6})$ = **14,3 s** ; double précision $ln(9,1 × 10^{14})$ = **34,4 s**.`),
        p(`L’écart mesuré atteint 0,1 vers 9,5 s, **avant** les 14 s prévues, parce que l’arrondi n’est pas commis une seule fois au départ : il est recommencé à **chacun des 20 000 pas** et l’écart est réalimenté en continu. L’erreur effective de départ est donc plus grande que $u$ : d’après le tableau (3,9 × 10⁻⁵ à t = 2 pour une croissance de rapport $e^{λ t}$), de l’ordre de quelques $10^{−6}$.`),
      ],
    },
    {
      level: 6,
      title: "Le piège du décalage",
      statement: [
        p(`Un float a 24 bits de mantisse : le plus petit écart entre deux nombres voisins vaut $2^{e−23}$ si le nombre est entre $2^e$ et $2^{e+1}$. (a) Calculez cet écart pour un nombre voisin de 1000 (entre 512 et 1024) et pour un nombre voisin de 1. (b) Comparez à l’adoucissement ε = 0,05 des étoiles. (c) Que prévoyez-vous pour le gain du centrage ? Le logiciel mesure 345.`),
      ],
      solution: [
        p(`(a) Pour 1000 : $e = 9$, écart $2^{9−23} = 2^{−14}$ = **6,1 × 10⁻⁵**. Pour 1 : $e = 0$, écart $2^{−23}$ = **1,2 × 10⁻⁷**.`),
        p(`(b) Rapportés à 0,05 : **1,2 × 10⁻³** (0,12 %) pour le nombre voisin de 1000, contre **2,4 × 10⁻⁶** pour le nombre voisin de 1 : les positions décalées sont arrondies avec une erreur relative 512 fois plus grande.`),
        p(`(c) Le gain attendu est d’un ordre de **500**, voisin du **345** mesuré (les erreurs ne s’ajoutent pas toutes de la même façon).`),
      ],
    },
    {
      level: 6,
      title: "Énergie d’une collision d’amas",
      statement: [
        p(`Deux amas de Plummer de masse 0,5 chacun et de rayon d’échelle a = 0,6 sont à 6 unités l’un de l’autre et se rapprochent à la vitesse relative 0,4 (G = 1). Estimez l’énergie totale avec $E_{propre} = −(3π/64) M²/a$ pour chacun, une énergie d’interaction $−M_1 M_2/d$ et l’énergie cinétique du mouvement relatif. Le système est-il lié ?`),
      ],
      solution: [
        p(`Énergie propre de chaque amas : $−(3π/64) × 0,25/0,6 = −0,0614$, soit −0,1227 pour les deux.`),
        p(`Interaction : $−0,5 × 0,5/6 = −0,0417$. Cinétique du mouvement relatif : deux masses de 0,5 à 0,2 chacune : $2 × ½ × 0,5 × 0,2² = 0,020$.`),
        p(`Total : $−0,1227 − 0,0417 + 0,020$ = **−0,144**, valeur négative : le système est **lié** et les deux amas fusionneront. Le logiciel mesure −0,148 (écart de 3 % dû à la troncature des rayons à 10 a).`),
      ],
    },
  ],
  retenir: box("retenir", "", [
    ul([
      "Un calcul **N corps** coûte $N(N−1)/2$ paires : le travail croît comme $N²$. Une carte graphique le partage entre des milliers de calculateurs (ici 80 à 150 fois plus vite que le processeur sur 16 000 étoiles).",
      "En **simple précision** l’erreur d’une somme de $N$ termes est de l’ordre de $u sqrt{N}$, avec $u = 6 × 10^{−8}$ : on **centre** les positions en double avant de convertir (gain 345 mesuré).",
      "L’énergie dérive autant que l’erreur du schéma ; l’impulsion et le moment cinétique, exacts à l’arrondi en double, dérivent de $10^{−9}$ à $10^{−6}$ en simple précision.",
      "Le **chaos** interdit de comparer des trajectoires longues : on compare des accélérations à état fixé et des grandeurs d’ensemble ; « le processeur vérifie la carte graphique ».",
    ]),
  ]),
};
