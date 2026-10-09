# Générateur du cours « La mécanique par la simulation »

Ces fichiers fabriquent le cours Word (`docs/cours/La-mecanique-par-la-simulation.docx`) et sa version PDF. Le Word livré peut être modifié directement ; ce dossier sert à **régénérer** le cours (corriger un chiffre, ajouter un chapitre ou un exercice) en gardant la même mise en page.

## Contenu

| Fichier | Rôle |
|---|---|
| `content/front.js`, `ch1.js` … `ch7.js`, `back.js` | le texte du cours : couverture et avant-propos, un fichier par chapitre (six niveaux, logiciel, exercices et corrigés, « À retenir »), conclusion et annexes (corrigés, glossaire, formulaire, logiciel, références) |
| `lib/blocks.js`, `lib/render.js`, `lib/tex.js`, `lib/theme.js` | mise en page (encadrés, tableaux, figures, listes), mini-langage de formules converti en vraies équations Word (OMML), couleurs des niveaux |
| `build.js` | assemble le `.docx`, y compris la table des matières (champ Word avec numéros de page en cache) |
| `finalize.py`, `convert.py` | convertit en PDF avec LibreOffice, relève la page de chaque titre et recommence jusqu'à ce que la table des matières soit stable |
| `figdata.cpp`, `figures.py` | figures tirées de **vraies exécutions** des solveurs du projet (convergence, énergie de l'oscillateur, orbites de Kepler, chaos du pendule double) et des mesures de `--gpu-test` |
| `exos_check.py` | recalcule toutes les réponses numériques des exercices (aucune valeur n'est écrite à la main) |
| `crop_shots.py` | découpe les captures de `tools/screenshot.ps1` pour les figures du logiciel |

## Fabriquer le cours

Prérequis : Node.js, Python 3 avec `pymupdf`, `matplotlib`, `pillow`, `numpy`, LibreOffice (pour le PDF), et le dépôt compilé (`build/libphysicslab_core.a`).

```bash
cd tools/cours
npm install                                   # docx, jszip
# 1. données et figures
g++ -std=c++17 -O1 -I ../../include figdata.cpp ../../build/libphysicslab_core.a -o figdata && mkdir -p data && ./figdata data
python figures.py                             # écrit img/fig_*.png
python exos_check.py                          # réponses numériques des exercices
# 2. captures du logiciel (PowerShell, depuis la racine du dépôt) :
#    .\tools\screenshot.ps1 -Level 3 -Sim 1 -Wait 9 -Out s01.png       (de même s02 : niveau 4 sim 2 ; s03 : 4 / 3 ; s04 : 5 / 4 ;
#    s05 : 5 / 5 ; s06 : 5 / 6 ; s07 : 4 / 7 ; s08 : 4 / 8 ; s09 : 5 / 9 ; s10 : 5 / 10 ;
#    s11_l1 : niveau 1 sim 11 ; s11_l5 : niveau 5 sim 11 -Clicks "175,411" -Wait 12)
python crop_shots.py
# 3. le cours
python finalize.py                            # cours.docx + cours.pdf, table des matières calculée en deux passes
```

Le chemin de LibreOffice est écrit dans `convert.py` (`C:\Program Files\LibreOffice\program\soffice.exe`). Pour valider le `.docx` : `python <compétence docx>/scripts/office/validate.py cours.docx`.

## Syntaxe du contenu

Marqueurs dans le texte : `**gras**`, `//italique//`, `{{code}}`, `[[RRGGBB|texte coloré]]`, `$formule$`. Les formules s'écrivent **sans antislash** : `frac{a}{b}`, `sqrt{x}`, `paren{…}`, `sum_{i}^{n}`, `x_{i}`, `x^{2}`, `text{…}`, les fonctions `sin cos tan ln log exp` suivies de leur argument ; lettres grecques et symboles en Unicode.

## Pièges rencontrés (à ne pas refaire)

- **LibreOffice sert à fabriquer le PDF** et interprète certains mots des formules comme des mots-clés de son propre langage : un indice `lim`, un `max`, un `min` s'affichent mal (un « ¿ » apparaît). Éviter ces mots dans les formules (on note la vitesse limite `v_{∞}`). Une formule qui **commence ou finit par « = »** donne aussi un « ¿ » : mettre le signe égal dans le texte, hors de `$…$`. La flèche de vecteur (caractère combinant) se détache de sa lettre : ne pas l'utiliser.
- Le repère de la table des matières doit être suivi d'un paragraphe vide : la rupture de section est portée par le dernier paragraphe de la section, et le remplacer fusionnerait les sections.
- Vérifier le PDF à l'œil (PyMuPDF : `python render.py cours.pdf pg`) et par une recherche de glyphes suspects (« ¿ », U+FFFD) dans le texte extrait.
