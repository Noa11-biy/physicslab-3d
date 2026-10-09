"""Calcule les réponses numériques des exercices du cours (aucune valeur n'est écrite à la main)."""
import csv
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
pi = math.pi
g = 9.81


def hdr(s):
    print("\n==", s)


# ---------------- chapitre 1
hdr("1.3 Euler chute libre g=10 dt=0,5, 4 pas")
v = x = 0.0
rows = []
for n in range(4):
    x, v = x + v * 0.5, v + 10 * 0.5
    rows.append((n + 1, v, x))
print(rows, "exact v=20 x=20")
hdr("1.4 parachutiste m=80 k=20")
m, k = 80, 20
vl, tau = m * 10 / k, m / k
print("vlim", vl, "tau", tau, "v(tau)", vl * (1 - math.exp(-1)), "v(3tau)", vl * (1 - math.exp(-3)), "t90", tau * math.log(10))
hdr("1.5 énergie gagnée Euler m=2,g=10,t=2,dt=0,5")
m = 2
v = x = 0.0
for n in range(4):
    x, v = x + v * 0.5, v + 10 * 0.5
print("E = ", 0.5 * m * v * v - m * 10 * x, " formule", 0.5 * m * 100 * 2 * 0.5)
hdr("1.7 Euler v'=g-v/tau, g=10, tau=0,5, dt=0,1")
v = 0.0
for n in range(5):
    v = v + 0.1 * (10 - v / 0.5)
ex = 5 * (1 - math.exp(-1))
print("v Euler t=0,5:", v, "exact", ex, "erreur rel", (v - ex) / ex)
v = 0.0
for n in range(8):
    v = v + 1.2 * (10 - v / 0.5)
    print(" dt=1,2 pas", n + 1, v)
hdr("1.8 ordres (CSV)")
rows = list(csv.DictReader(open(os.path.join(HERE, "data", "convergence.csv"))))
for r in rows:
    print(r["steps"], r["euler"], r["verlet"], r["rk4"])
for key in ["euler", "verlet", "rk4"]:
    e = [float(r[key]) for r in rows]
    print(key, "ordres successifs", [round(math.log2(e[i] / e[i + 1]), 2) for i in range(len(e) - 1)])
hdr("1.9 pas RK45")
print(0.9 * 0.2 * (1e-6 / 3.2e-5) ** 0.2)

# ---------------- chapitre 2
hdr("2.3 ressort m=0,5 k=200 A=4 cm")
m, k, A = 0.5, 200, 0.04
w0 = math.sqrt(k / m)
print("w0", w0, "T", 2 * pi / w0, "f", w0 / 2 / pi, "E", 0.5 * k * A * A, "vmax", A * w0)
hdr("2.4 masse 0,2 kg allonge 5 cm")
kk = 0.2 * g / 0.05
print("k", kk, "T", 2 * pi * math.sqrt(0.2 / kk))
hdr("2.5 v à x=A/2")
print(math.sqrt(3) / 2 * A * w0)
hdr("2.6 amortissement m=1,k=10,c=2")
m, k, c = 1, 10, 2
zeta = c / (2 * math.sqrt(k * m))
w0 = math.sqrt(k / m)
gam = c / (2 * m)
wd = math.sqrt(w0 ** 2 - gam ** 2)
Td = 2 * pi / wd
print("zeta", zeta, "wd", wd, "Td", Td, "facteur par période", math.exp(-gam * Td))
hdr("2.7 Euler énergie, T=1 s, dt=0,01")
f1 = 1 + (2 * pi * 0.01) ** 2
print("par pas", f1, "par période", f1 ** 100, "10 périodes", f1 ** 1000)
hdr("2.8 RK4 énergie")
y = 2 * pi / 20
per_step = 1 - y ** 6 / 72 + y ** 8 / 576
print("par pas", per_step, "après 600 pas", per_step ** 600, "1-y^6/72", 1 - y ** 6 / 72)
y = 0.1
n = 2 * pi / 0.1
print("y=0,1 : perte par période ≈", n * y ** 6 / 72)
print("Euler 20 pas/période, 30 périodes:", (1 + y ** 2) ** 0, (1 + (2 * pi / 20) ** 2) ** 600)

# ---------------- chapitre 3
hdr("3.2 T pendule")
print("L=0,5:", 2 * pi * math.sqrt(0.5 / g), " L pour T=2 s:", g * 4 / (4 * pi ** 2), " Lune:", 2 * pi * math.sqrt(1 / 1.62))
hdr("3.4 pendule L=1,5, 60°")
L = 1.5
th = math.radians(60)
vb = math.sqrt(2 * g * L * (1 - math.cos(th)))
print("v bas", vb, "tension m=0,2", 0.2 * g * (3 - 2 * math.cos(th)))
hdr("3.5 série période à 60°")


def agm(a, b):
    for _ in range(40):
        a, b = (a + b) / 2, math.sqrt(a * b)
    return a


def K(k):
    return pi / (2 * agm(1, math.sqrt(1 - k * k)))


th0 = math.radians(60)
ratio = 2 * K(math.sin(th0 / 2)) / pi
print("exact", ratio, "série", 1 + th0 ** 2 / 16 + 11 * th0 ** 4 / 3072, "ordre 1", 1 + th0 ** 2 / 16)
for deg in (10, 30, 60, 90, 120, 170):
    t = math.radians(deg)
    print(deg, "°", round(2 * K(math.sin(t / 2)) / pi, 4))
hdr("3.6 Lyapunov")
lam = math.log(4400) / 8
print("lambda", lam, "horizon", math.log(1e8) / lam)
hdr("3.8 horizons float / double (λ=1)")
for name, d0 in (("float", 6e-8), ("double", 1.1e-16)):
    print(name, math.log(0.1 / d0) / 1.0)

# ---------------- chapitre 4
hdr("4.2 Kepler")
for nom, a in (("Mars", 1.524), ("Jupiter", 5.203)):
    print(nom, a ** 1.5)
hdr("4.3 géostationnaire")
GM = 3.986004e14
T = 86164.1
r = (GM * T ** 2 / (4 * pi ** 2)) ** (1 / 3)
print("r", r, "alt km", (r - 6378137) / 1e3)
hdr("4.4 ISS")
r = 6378e3 + 400e3
v = math.sqrt(GM / r)
print("v", v, "T s", 2 * pi * r / v, "min", 2 * pi * r / v / 60, "v_lib", math.sqrt(2) * v)
hdr("4.5 énergie massique Terre")
GMs = 1.32712440018e20
a = 1.495978707e11
print("E/m", -GMs / (2 * a), "v", math.sqrt(GMs / a))
hdr("4.6 Lagrange")
print("omega", math.sqrt(3), "T", 2 * pi / math.sqrt(3))
hdr("4.7 précession numérique")
for steps in (100, 200, 400):
    h = 2 * pi / steps
    e = 0.5
    dw = (pi / 8) * h * h * (4 + e * e) / (1 - e * e) ** 3
    print(steps, "pas/orbite : Δω =", dw, "rad =", math.degrees(dw), "°")
hdr("4.8 huit")
pos = [(0.97000436, -0.24308753), (-0.97000436, 0.24308753), (0.0, 0.0)]
vel = [(0.466203685, 0.43236573), (0.466203685, 0.43236573), (-0.93240737, -0.86473146)]
Tk = 0.5 * sum(vx * vx + vy * vy for vx, vy in vel)
U = 0.0
for i in range(3):
    for j in range(i + 1, 3):
        U -= 1.0 / math.hypot(pos[i][0] - pos[j][0], pos[i][1] - pos[j][1])
print("T", Tk, "U", U, "E", Tk + U)
hdr("4.9 amplification")
print("lambda", math.log(6.6e4) / 14)

# ---------------- chapitre 5
hdr("5.2 plan incliné")
th = math.radians(20)
print("tan20", math.tan(th), "a (mu=0,3)", g * (math.sin(th) - 0.3 * math.cos(th)))
hdr("5.3 rebond")
h0, e = 1.6, 0.5
v = math.sqrt(2 * g * h0)
print("v impact", v, "v après", e * v, "h1", e * e * h0)
hdr("5.4 Zénon")
t0 = math.sqrt(2 * h0 / g)
print("t0", t0, "total", t0 * (1 + e) / (1 - e))
print("n pour h < 1 mm:", math.log(h0 / 1e-3) / math.log(1 / e ** 2))
hdr("5.5 montée sur plan")
th = math.radians(30)
a = g * (math.sin(th) + 0.2 * math.cos(th))
print("décélération", a, "distance", 16 / (2 * a), "temps", 4 / a)
hdr("5.6 choc 1D")
m1, v1, m2, v2, e = 2, 3, 1, 0, 0.5
w1 = ((m1 - e * m2) * v1 + (1 + e) * m2 * v2) / (m1 + m2)
w2 = ((m2 - e * m1) * v2 + (1 + e) * m1 * v1) / (m1 + m2)
print(w1, w2, "p", m1 * w1 + m2 * w2, "Ec avant", 0.5 * m1 * v1 ** 2, "après", 0.5 * m1 * w1 ** 2 + 0.5 * m2 * w2 ** 2)
hdr("5.7 Hertz")
c = 4 / 5 * math.gamma(2 / 5) * math.gamma(1 / 2) / math.gamma(9 / 10)
print("constante", c)
for v0 in (1.0, 2.0):
    mr, kh = 0.5, 1e4
    dm = (5 * mr * v0 ** 2 / (4 * kh)) ** 0.4
    print("v0", v0, "delta_max", dm, "T", c * dm / v0)
hdr("5.8 bissection et Hunt-Crossley")
print("bissection", math.ceil(math.log2(0.05 / 1e-10)))
for alpha in (0.01, 0.1):
    ee = 1 / (1 + alpha * 1.0)
    print(alpha, "e", ee, "perte énergie", 1 - ee ** 2)

# ---------------- chapitre 6
hdr("6.2 inerties")
M, R = 2, 0.1
print("disque", 0.5 * M * R ** 2, "sphère", 0.4 * M * R ** 2, "L", 0.5 * M * R ** 2 * 20, "E", 0.5 * 0.5 * M * R ** 2 * 400)
hdr("6.3 patineur")
I1, I2 = 3, 1.5
w1 = 4 * pi
w2 = I1 * w1 / I2
print("w2 tr/s", w2 / 2 / pi, "E1", 0.5 * I1 * w1 ** 2, "E2", 0.5 * I2 * w2 ** 2)
hdr("6.4 toupie")
m, l, I1, I3, w3 = 1, 0.5, 1.2, 0.4, 25
print("rapide", m * 9.80665 * l / (I3 * w3))
th = 0.6
a = I1 * math.cos(th)
b = -I3 * w3
cc = m * 9.80665 * l
phi = (-b - math.sqrt(b * b - 4 * a * cc)) / (2 * a)
print("exact lente", phi, "période", 2 * pi / phi, "cos", math.cos(th))
print("spin critique", 2 * math.sqrt(I1 * m * 9.80665 * l) / I3)
hdr("6.6 axe intermédiaire")
I = (1, 2, 3)
for w2 in (1, 2):
    lam = w2 * math.sqrt((I[1] - I[0]) * (I[2] - I[1]) / (I[0] * I[2]))
    print("w2", w2, "lambda", lam, "facteur sur 3 s", math.exp(lam * 3))
hdr("6.7 corps libre symétrique")
print("Omega", (0.5 - 1) * 3 / 1)
hdr("6.9 norme du quaternion Euler")
for N in (100, 2000):
    print(N, (1 + 0.05 ** 2 * 4 / 4) ** (N / 2))

# ---------------- chapitre 7
hdr("7.2 paires")
for N in (500, 2000, 10000, 100000):
    print(N, N * (N - 1) // 2, N * N)
hdr("7.3 temps")
print("GPU", 1e10 / 3e10, "CPU", 1e10 / 0.33e9, "rapport", (1e10 / 0.33e9) / (1e10 / 3e10))
hdr("7.4 pas par image")
print("steps/frame", 14 / 1.5, "temps simulé par s", 9 * 60 * 0.004)
hdr("7.5 erreur float")
u = 2.0 ** -24
for N in (100, 1000, 5000):
    print(N, u * math.sqrt(N))
print("N pour 1e-4:", (1e-4 / u) ** 2)
hdr("7.6 tuiles")
print(math.ceil(2000 / 256), math.ceil(100000 / 256))
hdr("7.7 centrage")
print("ulp(1000)", 2.0 ** (9 - 23), "ulp(1)", 2.0 ** -23, "rapport", 2.0 ** (9 - 23) / 2.0 ** -23)
print("relatif à 0,05:", 2.0 ** (9 - 23) / 0.05, 2.0 ** -23 / 0.05)
