const fs = require("fs");
const path = require("path");
const D = require("docx");
const T = require("../lib/theme");
const B = require("../lib/blocks");

const none = { style: D.BorderStyle.NONE, size: 0, color: "FFFFFF" };

function cover(R) {
  const w = T.textWidth;
  const band = (children) =>
    new D.Table({
      width: { size: w, type: D.WidthType.DXA },
      columnWidths: [w],
      borders: { top: none, bottom: none, left: none, right: none, insideHorizontal: none, insideVertical: none },
      rows: [
        new D.TableRow({
          children: [
            new D.TableCell({
              width: { size: w, type: D.WidthType.DXA },
              shading: { type: D.ShadingType.CLEAR, fill: T.colors.navy, color: "auto" },
              margins: { top: 420, bottom: 420, left: 520, right: 520 },
              borders: { top: none, bottom: { style: D.BorderStyle.SINGLE, size: 36, color: T.colors.teal }, left: none, right: none },
              children,
            }),
          ],
        }),
      ],
    });
  const img = fs.readFileSync(path.join(__dirname, "..", "img", "cover_scene.png"));
  const kids = [
    new D.Paragraph({ spacing: { before: 200, after: 0 }, children: [] }),
    band([
      new D.Paragraph({ spacing: { after: 160 }, children: [new D.TextRun({ text: "PHYSICSLAB 3D  ·  DOMAINE MÉCANIQUE", size: 22, color: "9FD6D7", bold: true, characterSpacing: 40 })] }),
      new D.Paragraph({ spacing: { after: 200, line: 240 }, children: [new D.TextRun({ text: "La mécanique par la simulation", size: 74, bold: true, color: "FFFFFF" })] }),
      new D.Paragraph({ spacing: { after: 0 }, children: [new D.TextRun({ text: "Un cours à six niveaux, de la vulgarisation au doctorat, avec exercices corrigés", size: 30, color: "DCE6F2" })] }),
    ]),
    new D.Paragraph({ spacing: { after: 280 }, children: [] }),
    new D.Paragraph({
      alignment: D.AlignmentType.CENTER,
      spacing: { after: 100 },
      children: [new D.ImageRun({ type: "png", data: img, transformation: { width: 450, height: 326 }, altText: { title: "Double pendule", description: "Trajectoires d’un pendule double chaotique, dans le logiciel PhysicsLab 3D", name: "couverture" } })],
    }),
    new D.Paragraph({ alignment: D.AlignmentType.CENTER, spacing: { after: 360 }, children: [new D.TextRun({ text: "Un pendule double, dessiné par le logiciel : deux tiges, un mouvement imprévisible.", italics: true, size: 19, color: T.colors.grey })] }),
  ];
  const stat = (n, label) => [
    new D.Paragraph({ alignment: D.AlignmentType.CENTER, spacing: { after: 0 }, children: [new D.TextRun({ text: n, bold: true, size: 44, color: T.colors.teal })] }),
    new D.Paragraph({ alignment: D.AlignmentType.CENTER, spacing: { after: 0 }, children: [new D.TextRun({ text: label, size: 19, color: T.colors.grey })] }),
  ];
  const cw = [Math.floor(w / 4), Math.floor(w / 4), Math.floor(w / 4), w - 3 * Math.floor(w / 4)];
  const cell = (c, i) => new D.TableCell({ width: { size: cw[i], type: D.WidthType.DXA }, borders: { top: none, bottom: none, left: none, right: none }, margins: { top: 60, bottom: 60, left: 60, right: 60 }, children: c });
  kids.push(
    new D.Table({
      width: { size: w, type: D.WidthType.DXA },
      columnWidths: cw,
      borders: { top: none, bottom: none, left: none, right: none, insideHorizontal: none, insideVertical: none },
      rows: [new D.TableRow({ children: [cell(stat("7", "chapitres"), 0), cell(stat("6", "niveaux de lecture"), 1), cell(stat("11", "simulations"), 2), cell(stat("64", "exercices corrigés"), 3)] })],
    }),
    new D.Paragraph({ spacing: { before: 360, after: 40 }, alignment: D.AlignmentType.CENTER, children: [new D.TextRun({ text: "Édition 1.0  ·  octobre 2026", size: 22, bold: true, color: T.colors.navy })] }),
    new D.Paragraph({ alignment: D.AlignmentType.CENTER, spacing: { after: 0 }, children: [new D.TextRun({ text: "Document compagnon du logiciel libre PhysicsLab 3D (licence MIT)", size: 19, color: T.colors.grey })] })
  );
  return kids;
}

const avantPropos = [
  B.lead("Ce cours raconte comment un ordinateur calcule le mouvement des choses, depuis la balle qu’on lance jusqu’à mille étoiles qui s’attirent. Il s’appuie sur un logiciel que l’on peut manipuler, et chaque idée y est expliquée six fois, du plus simple au plus technique."),
  B.h3("À qui s’adresse ce cours ?"),
  B.p("À tout le monde, à condition de choisir le bon niveau de lecture. Un professeur peut s’en servir pour préparer un cours, un élève pour découvrir, un étudiant pour s’entraîner, un curieux pour comprendre ce que fait vraiment un « moteur physique » de jeu vidéo ou de film d’animation. **Aucun prérequis** n’est nécessaire pour les deux premiers niveaux."),
  B.tbl(
    ["Niveau", "Public visé", "Ce que l’on y trouve", "Outils mathématiques"],
    [
      ["**N1** Vulgarisation", "Tout public, sans bagage scientifique", "Images, analogies, histoires ; aucune équation", "Aucun"],
      ["**N2** Intéressé", "Curieux, lecteurs de science", "Les concepts, les formules dites en mots, les ordres de grandeur", "Proportionnalité"],
      ["**N3** Collège", "Élèves de collège", "Formules simples, unités, calculs numériques à la main", "Quatre opérations, puissances"],
      ["**N4** Lycée", "Élèves de lycée", "Vecteurs, lois de Newton, énergie, équations différentielles simples", "Dérivées, exponentielle, trigonométrie"],
      ["**N5** Étudiant (L1 à L3)", "Licence de physique, maths, ingénierie", "Équations différentielles, Lagrange et Hamilton, analyse de stabilité", "Algèbre linéaire, développements limités"],
      ["**N6** Doctorant / Chercheur", "Master, doctorat, recherche", "Formulations complètes, schémas avancés, analyse d’erreur, limites connues", "Analyse numérique, mécanique analytique"],
    ],
    [2000, 2000, 3098, 2200],
    { size: 19 }
  ),
  B.h3("Comment lire ce cours"),
  B.p("Chaque chapitre traite un phénomène, du plus simple (une balle) au plus riche (mille étoiles sur une carte graphique). Dans chaque chapitre, les **six niveaux** se suivent : chacun reprend la même idée avec plus de précision que le précédent. On peut donc :"),
  B.ul([
    "**lire en largeur** : les niveaux N1 et N2 de tous les chapitres forment un petit livre de culture générale (une vingtaine de pages) ;",
    "**lire en profondeur** : choisir un chapitre et monter de N1 à N6 jusqu’à ce que cela devienne trop technique, puis s’arrêter, sans rien perdre ;",
    "**s’entraîner** : les exercices sont rangés par niveau et leurs corrigés détaillés sont à l’annexe A.",
  ]),
  B.h3("Les encadrés"),
  B.box("plus", "", [
    B.p("[[1F3A5F|À retenir]] : le résumé d’une idée essentielle, à connaître.", { after: 60 }),
    B.p("[[0F8B8D|Dans le logiciel]] : quoi cliquer dans PhysicsLab 3D, et ce qu’il faut observer.", { after: 60 }),
    B.p("[[3C8D5A|Exemple]] : un calcul complet, chiffres à l’appui.", { after: 60 }),
    B.p("[[B03A2E|Attention]] : une erreur classique, ou un piège que l’ordinateur nous a vraiment tendu pendant la construction du logiciel.", { after: 60 }),
    B.p("[[7B4B94|Pour aller plus loin]] : une piste pour les plus curieux (références en annexe E).", { after: 0 }),
  ]),
  B.h3("La méthode du cours"),
  B.p("Chaque phénomène a été étudié avec la même méthode en six étapes, que l’on retrouve d’un chapitre à l’autre :"),
  B.ol([
    "**dériver** les équations à la main ;",
    "**résoudre exactement** un cas de référence, quand c’est possible ;",
    "**discrétiser** : transformer l’équation en une suite de petits calculs ;",
    "**programmer** (en double précision d’abord, sur carte graphique ensuite si le cas est lourd) ;",
    "**valider** : conservation de l’énergie, comparaison à la solution exacte, ordre de convergence, comparaison processeur contre carte graphique ;",
    "**décliner** l’explication sur les six niveaux.",
  ]),
  B.h3("Le logiciel"),
  B.p("Tout ce que le cours décrit peut être observé dans **PhysicsLab 3D**, un logiciel libre (licence MIT) écrit en C++ avec OpenGL. On le lance avec le niveau et la simulation voulus :"),
  B.p("{{physicslab --level 4 --sim 2}}", { align: D.AlignmentType.CENTER }),
  B.p("Le sélecteur de niveau du logiciel correspond exactement aux six niveaux de ce cours : changer de niveau change les explications, les curseurs offerts, les équations affichées et le détail des courbes. L’annexe D récapitule les options et la liste des onze simulations."),
  B.box("attention", "Une promesse : ce qui est écrit a été mesuré", [
    B.p("Les chiffres de ce cours (erreurs, ordres de convergence, temps de calcul, écarts entre processeur et carte graphique) ont été relevés en exécutant le logiciel, et les propriétés les plus importantes sont vérifiées automatiquement par 76 tests. Les temps de calcul dépendent de la machine : ils ont été mesurés ici sur un ordinateur portable à carte graphique intégrée (Intel UHD), et varient d’un jour à l’autre ; les précisions, elles, se reproduisent à l’identique."),
    B.p("Plusieurs affirmations qui « allaient de soi » se sont révélées fausses à la mesure (une toupie qui ne tombe pas, une méthode d’ordre 4 qui n’en est pas une, un ordre de calcul qui ne change rien). Elles sont signalées par un encadré **Attention** : l’erreur est souvent plus instructive que la règle."),
  ]),
];

module.exports = { cover, avantPropos };
