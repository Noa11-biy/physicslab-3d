"""Construit le .docx, le convertit en PDF avec LibreOffice, relève la page de chaque titre et recommence
jusqu'à ce que les numéros de la table des matières soient stables."""
import json, os, re, subprocess, sys
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(HERE, "cours.docx")
PDF = os.path.join(HERE, "cours.pdf")
PAGES = os.path.join(HERE, "pages.json")


def run(cmd):
    r = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True, shell=isinstance(cmd, str))
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        raise SystemExit("échec : " + str(cmd))
    return r.stdout


def roman(n):
    out = ""
    for r, v in (("x", 10), ("ix", 9), ("v", 5), ("iv", 4), ("i", 1)):
        while n >= v:
            out += r
            n -= v
    return out


def norm(t):
    return re.sub(r"\s+", " ", t.replace(" ", " ")).strip()


def locate():
    headings = json.load(open(os.path.join(HERE, "headings.json"), encoding="utf-8"))
    doc = pymupdf.open(PDF)
    big = []     # titres (lignes en grande police) de chaque page
    kind = []    # type de page selon l'en-tête : toc, front, body
    for page in doc:
        lines = []
        for b in page.get_text("dict")["blocks"]:
            for l in b.get("lines", []):
                size = max(s["size"] for s in l["spans"])
                if size > 13:
                    lines.append(norm("".join(s["text"] for s in l["spans"])))
        big.append(lines)
        txt = page.get_text()
        head = norm(txt[:200])
        if "Table des matières" in head and "Avant-propos" not in head[:60]:
            kind.append("toc")
        else:
            kind.append("other")
    # première page du premier chapitre
    joined = [" ".join(ls) for ls in big]
    first_chapter = next(i for i, j in enumerate(joined) if "Chapitre 1 · Calculer" in j and kind[i] != "toc" and i > 0)
    pages = {}
    ptr = 0
    for h in headings:
        target = norm(h["text"])
        found = None
        for i in range(ptr, len(doc)):
            if kind[i] == "toc" or i == 0:
                continue
            if target in joined[i]:
                found = i
                break
        if found is None:
            raise SystemExit("titre introuvable dans le PDF : " + target)
        ptr = found
        if found < first_chapter:      # avant-propos : chiffres romains (la table des matières porte le numéro i)
            tocpages = [i for i, k in enumerate(kind) if k == "toc"]
            pages[h["bm"]] = roman(found - tocpages[0] + 1)
        else:
            pages[h["bm"]] = found - first_chapter + 1
    return pages, len(doc), first_chapter


for it in range(5):
    run(["node", "build.js", DOCX])
    run([sys.executable, "convert.py", DOCX])
    pages, n, first = locate()
    old = json.load(open(PAGES, encoding="utf-8")) if os.path.exists(PAGES) else {}
    print(f"passe {it + 1} : {n} pages, chapitre 1 en page physique {first + 1}")
    if pages == old:
        print("table des matières stable")
        break
    json.dump(pages, open(PAGES, "w", encoding="utf-8"), ensure_ascii=False, indent=0)
else:
    print("attention : la table des matières n'a pas convergé")
