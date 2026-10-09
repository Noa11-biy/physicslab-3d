import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
IMG = os.path.join(HERE, "img")
os.makedirs(IMG, exist_ok=True)

NAVY, TEAL, ORANGE, RED, PURPLE, GREEN, GREY = "#1F3A5F", "#0F8B8D", "#E07A1F", "#C0392B", "#7B4B94", "#3C8D5A", "#6B7280"
COL = {"euler": ORANGE, "symplectic": "#B5A800", "verlet": PURPLE, "rk4": GREEN}
NAME = {"euler": "Euler", "symplectic": "Euler symplectique", "verlet": "Verlet", "rk4": "RK4"}

plt.rcParams.update({
    "font.family": ["Calibri", "DejaVu Sans"],
    "font.size": 11,
    "axes.edgecolor": "#9AA3AE",
    "axes.labelcolor": "#222B36",
    "axes.titleweight": "bold",
    "axes.titlesize": 12,
    "axes.grid": True,
    "grid.color": "#E4E8EC",
    "grid.linewidth": 0.8,
    "legend.frameon": False,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
})


def read(name):
    with open(os.path.join(DATA, name), newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


# ---------------------------------------------------------------- convergence
d = read("convergence.csv")
fig, ax = plt.subplots(figsize=(6.4, 4.0))
slopes = {}
for key in ["euler", "symplectic", "verlet", "rk4"]:
    slopes[key] = -np.polyfit(np.log(d["steps"][-4:]), np.log(d[key][-4:]), 1)[0]
    ax.loglog(d["steps"], d[key], "o-", color=COL[key], label=f"{NAME[key]} (pente −{slopes[key]:.1f})".replace(".", ","), lw=2, ms=4)
    print("pente", key, round(slopes[key], 2))
ax.set_xlabel("nombre de pas (plus il est grand, plus le pas $dt$ est petit)")
ax.set_ylabel("erreur finale (espace des phases)")
ax.set_title("L'erreur diminue avec le pas, à des vitesses différentes")
ax.legend(loc="lower left")
fig.savefig(os.path.join(IMG, "fig_convergence.png"))
plt.close(fig)

# ---------------------------------------------------------------- énergie oscillateur
e = read("energy.csv")
fig, axs = plt.subplots(1, 3, figsize=(7.6, 3.0))
axs[0].semilogy(e["t"], e["euler"], color=COL["euler"], lw=2)
axs[0].set_title("Euler : elle explose", fontsize=10)
axs[0].set_ylabel("énergie / énergie de départ")
axs[1].plot(e["t"], e["symplectic"], color=COL["symplectic"], lw=0.9, label="Euler symplectique")
axs[1].plot(e["t"], e["verlet"], color=COL["verlet"], lw=1.2, label="Verlet")
axs[1].set_xlim(0, 4)
axs[1].set_title("Elle oscille (4 premières périodes)", fontsize=10)
axs[1].legend(loc="upper right", fontsize=8)
axs[2].plot(e["t"], e["verlet"], color=COL["verlet"], lw=1.2, label="Verlet")
axs[2].plot(e["t"], e["rk4"], color=COL["rk4"], lw=2.2, label="RK4")
axs[2].set_ylim(0.96, 1.01)
axs[2].set_title("Zoom sur 30 périodes : RK4 perd 0,8 %", fontsize=10)
axs[2].legend(loc="lower left", fontsize=8)
for ax in axs:
    ax.set_xlabel("temps (en périodes)")
fig.tight_layout()
fig.savefig(os.path.join(IMG, "fig_energie.png"))
plt.close(fig)
print("énergie finale (30 périodes) :", {k: float(e[k][-1]) for k in ["euler", "symplectic", "verlet", "rk4"]})

# ---------------------------------------------------------------- Kepler
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
for ax, (key, title, col) in zip(axs, [("euler", "Euler", ORANGE), ("verlet", "Verlet", PURPLE), ("rk4", "RK4", GREEN)]):
    k = read(f"kepler_{key}.csv")
    ax.plot(k["x"], k["y"], color=col, lw=0.7)
    ax.plot([0], [0], "o", color="#F2B705", ms=7)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    lim = max(abs(k["x"]).max(), abs(k["y"]).max()) * 1.05
    ax.set_xlim(-lim if key != "euler" else -lim, lim)
    ax.set_ylim(-lim, lim)
fig.suptitle("Dix tours d'une orbite elliptique, 200 pas par tour", fontsize=11, fontweight="bold", y=1.0)
fig.tight_layout()
fig.savefig(os.path.join(IMG, "fig_kepler.png"))
plt.close(fig)

# ---------------------------------------------------------------- chaos (pendule double)
c = read("double_pendulum.csv")
t, dist = c["t"], c["distance"]
mask = (dist > 3e-9) & (dist < 1e-1)
lam, lnd0 = np.polyfit(t[mask], np.log(dist[mask]), 1)
print("lambda ajusté (pendule double, RK45 serré) :", round(lam, 3), "n =", int(mask.sum()))
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.semilogy(t, dist, color=NAVY, lw=2, label="écart mesuré entre les deux pendules")
tt = np.linspace(t[mask][0], t[mask][-1] + 3, 50)
ax.semilogy(tt, np.exp(lnd0 + lam * tt), "--", color=ORANGE, lw=1.8, label=f"croissance exponentielle $e^{{\\lambda t}}$, $\\lambda \\approx$ {lam:.2f} /s")
ax.set_xlabel("temps (s)")
ax.set_ylabel("distance dans l'espace des phases")
ax.set_title("Deux pendules doubles qui diffèrent de 0,000000001 rad au départ")
ax.set_ylim(1e-10, 5)
ax.legend(loc="lower right")
fig.savefig(os.path.join(IMG, "fig_chaos.png"))
plt.close(fig)

# ---------------------------------------------------------------- rebonds (Zénon), formule exacte
g, h0, ee = 9.80665, 2.0, 0.8
t0 = math.sqrt(2 * h0 / g)
times, heights = [], []
tcur = 0.0
v = 0.0
# phase de chute puis rebonds successifs : durée du vol n = 2 t0 e^n
tt, hh = [], []
tfall = t0
x = np.linspace(0, t0, 60)
tt += list(x)
hh += list(h0 - 0.5 * g * x ** 2)
start = t0
for n in range(1, 40):
    vn = ee ** n * g * t0  # vitesse juste après le n-ième choc
    dur = 2 * vn / g
    x = np.linspace(0, dur, 80)
    tt += list(start + x)
    hh += list(vn * x - 0.5 * g * x ** 2)
    start += dur
ttot = t0 * (1 + ee) / (1 - ee)
print("durée totale de Zénon (e = 0,8, h0 = 2 m) :", round(ttot, 4), "s ; fin des rebonds simulés à", round(start, 4), "s")
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(tt, hh, color=NAVY, lw=1.8)
ax.axvline(ttot, color=ORANGE, ls="--", lw=1.5)
ax.text(ttot + 0.1, 1.2, f"arrêt à $t_0\\,(1+e)/(1-e)$\n= {ttot:.3f} s", color=ORANGE, fontsize=10)
ax.set_xlim(0, 7)
ax.set_xlabel("temps (s)")
ax.set_ylabel("hauteur (m)")
ax.set_title("Une balle ($e = 0{,}8$) lâchée de 2 m : une infinité de rebonds en 5,75 s")
fig.savefig(os.path.join(IMG, "fig_rebonds.png"))
plt.close(fig)

# ---------------------------------------------------------------- solide : axe intermédiaire (équations d'Euler, RK4)
I = np.array([1.0, 2.0, 3.0])


def rhs(w):
    return np.array([(I[1] - I[2]) / I[0] * w[1] * w[2], (I[2] - I[0]) / I[1] * w[2] * w[0], (I[0] - I[1]) / I[2] * w[0] * w[1]])


def run(w0, tend=20.0, dt=0.002):
    w = np.array(w0, float)
    ts, out = [0.0], [w.copy()]
    n = int(tend / dt)
    for k in range(1, n + 1):
        k1 = rhs(w); k2 = rhs(w + 0.5 * dt * k1); k3 = rhs(w + 0.5 * dt * k2); k4 = rhs(w + dt * k3)
        w = w + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k % 25 == 0:
            ts.append(k * dt)
            out.append(w.copy())
    return np.array(ts), np.array(out)


eps = 1e-3
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.8), sharey=True)
cases = [("axe de plus petite inertie", [1, eps, eps], GREEN), ("axe intermédiaire", [eps, 1, eps], RED), ("axe de plus grande inertie", [eps, eps, 1], TEAL)]
for ax, (title, w0, col), idx in zip(axs, cases, [0, 1, 2]):
    ts, w = run(w0)
    others = [j for j in range(3) if j != idx]
    for j, ls in zip(others, ["-", "--"]):
        ax.plot(ts, w[:, j], ls, color=col, lw=1.4, label=f"$\\omega_{j+1}$")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("temps (s)")
    ax.legend(loc="upper left", fontsize=8)
axs[0].set_ylabel("composantes perpendiculaires\nà l'axe de rotation (rad/s)")
axs[0].set_ylim(-1.3, 1.3)
fig.tight_layout()
fig.savefig(os.path.join(IMG, "fig_axes.png"))
plt.close(fig)
ts, w = run([eps, 1, eps], tend=8.0, dt=0.001)
m = (ts > 4.0) & (ts < 8.0) & (np.abs(w[:, 0]) > 0)
lam_axis = np.polyfit(ts[m], np.log(np.abs(w[m, 0])), 1)[0]
print("croissance de l'axe intermédiaire (I = 1, 2, 3 ; ω2 = 1) : λ mesuré ≈", round(lam_axis, 4), " prévu", round(math.sqrt(1 / 3), 4))

# ---------------------------------------------------------------- GPU (valeurs relevées par --gpu-test, Release, Intel UHD)
N = [1000, 2000, 4000, 8000, 16000, 32000, 64000, 100000]
gpu = [5.70, 11.67, 18.42, 28.65, 31.47, 31.54, 32.57, 30.25]  # milliards d'interactions par seconde
cpu = [0.341, 0.322, 0.323, 0.321, 0.323]
Nerr = [3, 100, 1000, 5000]
err = [1.62e-8, 2.68e-7, 4.73e-7, 1.05e-6]
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3))
axs[0].loglog(N, gpu, "o-", color=ORANGE, lw=2, ms=5, label="carte graphique (float)")
axs[0].loglog(N[:5], cpu, "o-", color=NAVY, lw=2, ms=5, label="processeur (double)")
axs[0].set_xlabel("nombre d'étoiles $N$")
axs[0].set_ylabel("milliards d'interactions / s")
axs[0].set_title("Vitesse de calcul")
axs[0].legend(fontsize=9, loc="lower right")
axs[1].loglog(Nerr, err, "o-", color=ORANGE, lw=2, ms=5, label="erreur mesurée")
xs = np.array([3, 5000])
axs[1].loglog(xs, 6e-8 * np.sqrt(xs), "--", color=GREY, lw=1.4, label="$6\\times10^{-8}\\sqrt{N}$")
axs[1].set_xlabel("nombre d'étoiles $N$")
axs[1].set_ylabel("erreur relative des accélérations")
axs[1].set_title("Précision du calcul en float")
axs[1].legend(fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(IMG, "fig_gpu.png"))
plt.close(fig)
print("figures écrites dans", IMG)
