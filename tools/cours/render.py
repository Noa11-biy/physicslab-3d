import fitz, sys, os
pdf = sys.argv[1]
out = sys.argv[2]
pages = [int(x) for x in sys.argv[3:]] if len(sys.argv) > 3 else None
os.makedirs(out, exist_ok=True)
doc = fitz.open(pdf)
print("pages :", len(doc))
for i, page in enumerate(doc, 1):
    if pages and i not in pages:
        continue
    pix = page.get_pixmap(dpi=70)
    pix.save(os.path.join(out, f"p{i:03d}.png"))
