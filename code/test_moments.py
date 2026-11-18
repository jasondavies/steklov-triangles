"""Sanity test (numerical): rig.moments vs brute-force quadrature for random harmonic polynomials
on a random triangle."""
import sys; sys.path.insert(0, '.')
import numpy as np
from flint import acb, arb
from rig import moments, residual_sq
rng = np.random.default_rng(1)
Zc = [0j, 1 + 0j, 0.3 + 0.7j]
z0 = 0.4 + 0.2j; R = 0.9
cs = [rng.normal(size=6) + 1j * rng.normal(size=6) for _ in range(2)]
funcs = [([acb(arb(c.real), arb(c.imag)) for c in cc], acb(z0.real, z0.imag), arb(R)) for cc in cs]
M = moments(funcs, [acb(z.real, z.imag) for z in Zc])
def f(cc, z): w = (z - z0) / R; return sum(c * w ** (k + 1) for k, c in enumerate(cc))
def fp(cc, z): w = (z - z0) / R; return sum(c * (k + 1) * w ** k / R for k, c in enumerate(cc))
# area quadrature: collapsed Gauss
xg, wg = np.polynomial.legendre.leggauss(40); xg = (xg + 1) / 2; wg = wg / 2
A, B, C = Zc
Qn = np.zeros((2, 2, 2, 2))
for i, s in enumerate(xg):
    for j, t in enumerate(xg):
        z = A + s * (B - A) + s * t * (C - B)
        jac = s * abs(((B - A).conjugate() * (C - B)).imag)
        g = [np.array([fp(cc, z).real, -fp(cc, z).imag]) for cc in cs]   # grad Re f = (Re f', -Im f')
        for a in range(2):
            for b in range(2):
                Qn[a, b] += wg[i] * wg[j] * jac * np.outer(g[a], g[b])
err = max(abs(float(M['Q'][a][b][i][j].mid()) - Qn[a, b, i, j]) for a in range(2) for b in range(2) for i in range(2) for j in range(2))
print("max |Q_rig - Q_quad| =", err, " (scale", np.abs(Qn).max(), ")")
Gn = np.zeros((3, 2, 2)); mn = np.zeros((3, 2))
for e in range(3):
    a, b = Zc[e], Zc[(e + 1) % 3]
    z = a + xg * (b - a); L = abs(b - a)
    u = [f(cc, z).real for cc in cs]
    for x in range(2):
        mn[e, x] = L * np.sum(wg * u[x])
        for y in range(2):
            Gn[e, x, y] = L * np.sum(wg * u[x] * u[y])
print("max |G| err:", max(abs(float(M['G'][e][x][y].mid()) - Gn[e, x, y]) for e in range(3) for x in range(2) for y in range(2)),
      " max |m| err:", max(abs(float(M['m'][e][x].mid()) - mn[e, x]) for e in range(3) for x in range(2)))

okG = max(abs(float(M['G'][e][x][y].mid()) - Gn[e, x, y]) for e in range(3) for x in range(2) for y in range(2)) < 1e-12
okm = max(abs(float(M['m'][e][x].mid()) - mn[e, x]) for e in range(3) for x in range(2)) < 1e-12
ok = err < 1e-12 and okG and okm
import json
json.dump(dict(ok=bool(ok), max_err_Q=float(err)), open("../verification/test_moments.json", "w"), indent=1)
print("TEST_MOMENTS:", "ALL PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
