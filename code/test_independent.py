"""Independent numerical checks (brute-force quadrature, float64) of the two closed-form
upper-bound formulas used in the certified proof:
 (1) far_region.bound(p,q): Rayleigh quotient (times L) of u = cos(pi x) on T(p,q);
 (2) local.py: L_T * lambda_min(N_A, D_A) for the transplanted pair u_a o A^{-1}, recomputed
     here directly on T from the polynomial trial functions (no moment formulas)."""
import sys, math; sys.path.insert(0, '.')
import numpy as np
from flint import arb
xg, wg = np.polynomial.legendre.leggauss(60); xg = (xg + 1) / 2; wg = wg / 2
def area_quad(Z, F):
    A, B, C = Z; s = 0.0
    for i, a in enumerate(xg):
        z = A + a * (B - A) + a * xg * (C - B)
        jac = a * abs(((B - A).conjugate() * (C - B)).imag)
        s += wg[i] * np.sum(wg * jac * F(z))
    return s
def edges(Z):
    for j in range(3):
        a, b = Z[j], Z[(j + 1) % 3]
        yield a + xg * (b - a), wg * abs(b - a)
import far_region as R3
rng = np.random.default_rng(0)
worst = 0
for _ in range(20):
    p = rng.uniform(0.01, 0.5); q = rng.uniform(0.01, math.sqrt(2 * p - p * p))
    Z = [0j, 1 + 0j, complex(p, q)]
    N = area_quad(Z, lambda z: (math.pi * np.sin(math.pi * z.real)) ** 2)
    S = sum(np.sum(w * np.cos(math.pi * z.real) ** 2) for z, w in edges(Z))
    M = sum(np.sum(w * np.cos(math.pi * z.real)) for z, w in edges(Z))
    L = sum(np.sum(w) for z, w in edges(Z))
    direct = L * N / (S - M * M / L)
    formula = float(R3.bound(arb(p), arb(q))[0].mid())
    worst = max(worst, abs(direct - formula) / direct)
print("(1) cos(pi x) bound: max rel. difference formula vs quadrature =", worst)
worst1 = worst

import local as LC
c1, c2 = LC.c1, LC.c2
zc = complex(0.5, math.sqrt(3) / 6); R = 1 / math.sqrt(3)
def grad(cc, z):        # gradient of Re f, f = sum c_k ((z-zc)/R)^k
    w = (z - zc) / R; fp = sum(c * (k + 1) * w ** k / R for k, c in enumerate(cc))
    return fp.real, -fp.imag
def val(cc, z):
    w = (z - zc) / R; return sum(c * w ** (k + 1) for k, c in enumerate(cc)).real
worst = 0
for rho, psi in [(0.1, 4.0), (0.3, 3.8), (0.55, 4.5), (0.6, 0.3), (0.45, 2.0)]:
    a = 2 * math.cos(psi) / math.sqrt(3); b = 2 * math.sin(psi) / math.sqrt(3)
    A = np.array([[1, rho * a], [0, 1 + rho * b]]); Ai = np.linalg.inv(A)
    p = 0.5 + rho * math.cos(psi); q = math.sqrt(3) / 2 + rho * math.sin(psi)
    ZT = [0j, 1 + 0j, complex(p, q)]
    def pull(z):   # A^{-1} z
        x, y = z.real, z.imag
        return (Ai[0, 0] * x + Ai[0, 1] * y) + 1j * (Ai[1, 0] * x + Ai[1, 1] * y)
    cs = [c1, c2]
    Nm = np.zeros((2, 2)); Bm = np.zeros((2, 2)); mv = np.zeros(2); LT = 0
    for x in range(2):
        for y in range(2):
            def F(z, x=x, y=y):
                gx = np.array(grad(cs[x], pull(z))); gy = np.array(grad(cs[y], pull(z)))
                gx = Ai.T @ gx; gy = Ai.T @ gy        # grad(u o A^{-1}) = A^{-T} grad u
                return gx[0] * gy[0] + gx[1] * gy[1]
            Nm[x, y] = area_quad(ZT, F)
    for z, w in edges(ZT):
        LT += np.sum(w)
        vals = [val(cs[x], pull(z)) for x in range(2)]
        for x in range(2):
            mv[x] += np.sum(w * vals[x])
            for y in range(2):
                Bm[x, y] += np.sum(w * vals[x] * vals[y])
    Dm = Bm - np.outer(mv, mv) / LT
    direct = LT * min(np.linalg.eigvals(np.linalg.solve(Dm, Nm)).real)
    # formula route: X = rho*Y has a positive eigenvalue  <=>  bound < sigma* L_E ; compare values:
    Y, D = LC.Ybox(arb(rho), arb(psi))
    Ym = np.array([[float(Y[i][j].mid()) for j in range(2)] for i in range(2)])
    Dmf = np.array([[float(D[i][j].mid()) for j in range(2)] for i in range(2)])
    sig = float(LC.sig.mid())
    # N_A = sigma*(3/LT) D_A - rho Y  =>  bound = LT * lambda_min(N_A, D_A)
    NA = sig * (3 / LT) * Dmf - rho * Ym
    formula = LT * min(np.linalg.eigvals(np.linalg.solve(Dmf, NA)).real)
    worst = max(worst, abs(direct - formula) / direct, np.abs(Dm - Dmf).max())
    print("  rho=%.2f psi=%.1f: direct %.12f  formula %.12f" % (rho, psi, direct, formula))
print("(2) transplant bound: max rel. difference =", worst)
ok = worst1 < 1e-11 and worst < 1e-11
import json
json.dump(dict(ok=bool(ok), cos_rel_diff=float(worst1), transplant_rel_diff=float(worst)), open("../verification/test_independent.json", "w"), indent=1)
print("TEST_INDEPENDENT:", "ALL PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
