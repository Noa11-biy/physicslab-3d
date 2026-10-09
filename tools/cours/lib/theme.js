module.exports = {
  colors: {
    navy: "1F3A5F",
    teal: "0F8B8D",
    orange: "E07A1F",
    red: "B03A2E",
    purple: "7B4B94",
    green: "3C8D5A",
    grey: "6B7280",
    light: "F3F6F9",
    line: "D5DBE3",
    text: "1F2933",
    white: "FFFFFF",
  },
  // couleur et libellé de chaque niveau pédagogique (mêmes que dans le logiciel : 1 à 6)
  levels: [
    { n: 1, name: "Vulgarisation", color: "3C8D5A" },
    { n: 2, name: "Intéressé", color: "0F8B8D" },
    { n: 3, name: "Collège", color: "2E6DA4" },
    { n: 4, name: "Lycée", color: "C77D0A" },
    { n: 5, name: "Étudiant (L1 à L3)", color: "7B4B94" },
    { n: 6, name: "Doctorant / Chercheur", color: "A23B3B" },
  ],
  // types d'encadrés : couleur de la barre, fond clair
  boxes: {
    retenir: { color: "1F3A5F", fill: "EAF0F7", label: "À retenir" },
    logiciel: { color: "0F8B8D", fill: "E6F4F4", label: "Dans le logiciel" },
    attention: { color: "B03A2E", fill: "FBECEA", label: "Attention" },
    exemple: { color: "3C8D5A", fill: "EAF5EE", label: "Exemple" },
    plus: { color: "7B4B94", fill: "F3EDF7", label: "Pour aller plus loin" },
    enbref: { color: "E07A1F", fill: "FDF2E7", label: "En bref" },
  },
  fonts: { body: "Calibri", head: "Calibri", mono: "Consolas" },
  // A4, marges 2,3 cm à gauche et à droite : largeur de texte 9298 DXA (16,4 cm)
  page: { width: 11906, height: 16838, top: 1418, bottom: 1304, left: 1304, right: 1304 },
  textWidth: 9298,
};
