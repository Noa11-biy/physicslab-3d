"""Découpe les captures de tools/screenshot.ps1 (fenêtre PhysicsLab de 1630 x 952 pixels sur un écran 1920 x 1080)
en img/shot_<nom>.png, et fabrique la vue 3D de la couverture. Usage : python crop_shots.py [dossier des captures]"""
import glob, os, sys
from PIL import Image

here = os.path.dirname(os.path.abspath(__file__))
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("TEMP", "."), "physicslab_shots")
os.makedirs(os.path.join(here, "img"), exist_ok=True)
for f in sorted(glob.glob(os.path.join(src, "s*.png"))):
    name = os.path.splitext(os.path.basename(f))[0]
    crop = Image.open(f).convert("RGB").crop((146, 46, 1776, 998))
    crop.save(os.path.join(here, "img", f"shot_{name}.png"), optimize=True)
    print("shot_" + name)
s04 = os.path.join(here, "img", "shot_s04.png")
if os.path.exists(s04):
    Image.open(s04).crop((440, 70, 1226, 640)).save(os.path.join(here, "img", "cover_scene.png"), optimize=True)
