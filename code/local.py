"""
Part (b): RIGOROUS local maximality of the equilateral triangle, with an explicit radius.

Triangles T(p,q) = conv{(0,0),(1,0),(p,q)}, p = 1/2 + rho cos(psi), q = sqrt3/2 + rho sin(psi).
T = A(E) with A = I + rho [[0, a],[0, b]], a = 2cos(psi)/sqrt3, b = 2 sin(psi)/sqrt3.
Trial space on T: {u o A^{-1} : u in E*}, E* = EXACT first eigenspace of E (dim 2).
Exact moments of E* (energy tensor Q, edge Gram g^e, edge means m^e, eigenvalue sigma*)
are enclosed from the explicit equivariant pair (u1,u2) of equilateral.py plus Davis-Kahan
bounds (residual from equilateral.py, gap from the CR bound sigma_3(E) >= 2.2).

sigma_1(T) L_T < sigma* L_E  <=  X := sigma* (L_E/L_T) D_A - N_A has a positive eigenvalue
(with D_A positive definite).  Using the EXACT identities  sum_e m^e = 0  and  N_I = sigma* D_I,
X = rho * Y with the singularity-free expression
  Y = sigma* [ sum_e c_e g^e - 3 rho v v^T / L_T^2 ] - sum_ij H_ij Q_ij,
  H = (|det A| A^{-1}A^{-T} - I)/rho,  c_e = ((L_E/L_T) s_e - 1)/rho,  v = sum_e ((s_e-1)/rho) m^e.
We cover (rho, psi) in [0, rho1] x [0, 6.3] (6.3 > 2pi) by boxes and check, in ball arithmetic,
tr Y > 0 or det Y < 0, together with D_A > 0.  This proves sigma_1 L(T) < sigma_1 L(E) for
0 < rho <= rho1 (every direction), i.e. strict local maximality, including the conical rate.
"""
import sys, json
import numpy as np
sys.path.insert(0, '.')
from flint import acb, arb, ctx
from rig import moments, residual_sq, enc, dec

ctx.prec = 120
cert = json.load(open('../data/equilateral_pair.json'))          # committed certificate data
eq = json.load(open('../verification/equilateral.json'))         # certified output of equilateral.py
cr = json.load(open('../verification/cr_lower.json'))            # certified output of cr_lower.py
n = cert['n']; ks = cert['ks']; alpha = [float.fromhex(h) for h in cert['c_hex']]
sig = arb(float.fromhex(eq['sigma_tilde_hex'])); eta = dec(eq['eta_upper'])
s3lo = dec(cr['sigma3_lower'])
assert sig > 0          # Lemma 2.2 needs sigma~ >= 0
sq3 = arb(3).sqrt()
Zb = [acb(0), acb(1), acb(arb(1) / 2, sq3 / 2)]
zc = acb(arb(1) / 2, sq3 / 6); Rb = 1 / sq3
c1 = [0j] * n; c2 = [0j] * n
for a, k in zip(alpha, ks):
    c1[k - 1] = a * (-1j) ** k
    c2[k - 1] = ((-1j * a) if k % 3 == 1 else (1j * a)) * (-1j) ** k
cb = [[acb(arb(x.real), arb(x.imag)) for x in cc] for cc in (c1, c2)]
funcs = [(cb[0], zc, Rb), (cb[1], zc, Rb)]
M = moments(funcs, Zb)
# residuals and norms
r2 = []; u2 = []
for c in cb:
    a_, b_ = residual_sq(c, zc, Rb, Zb, sig); r2.append(a_); u2.append(b_)
# paper, Lemma 2.2: ||r_a|| <= eta ||u_a||_{L2(dE)}, since eta (equilateral.py) is the relative residual maximised over span(u1,u2)
rnorm = arb((eta * max(u2[0].sqrt().upper(), u2[1].sqrt().upper())).upper())
gap = s3lo - sig
assert gap > 0
assert s3lo > sig + eta, "sigma_3 lower bound must exceed the inclusion window (paper, Prop. 3.5)"
epsL = rnorm / gap
epsE = s3lo.sqrt() * rnorm / gap
print("||r|| <=", rnorm, " eps_L2(bdry) <=", epsL, " eps_energy <=", epsE)
sigstar = sig + arb(0, eta.upper())          # sigma* in [sig - eta, sig + eta]
def widen(x, r):
    return x + arb(0, r.upper())
# energy norms of tilde u
En = [(M['Q'][a][a][0][0] + M['Q'][a][a][1][1]).sqrt() for a in range(2)]
Bn = [u2[a].sqrt() for a in range(2)]
Q = [[[[widen(M['Q'][a][b][i][j], epsE * (En[a] + En[b] + epsE)) for j in range(2)] for i in range(2)]
      for b in range(2)] for a in range(2)]
g = [[[widen(M['G'][e][a][b], epsL * (Bn[a] + Bn[b] + epsL)) for b in range(2)] for a in range(2)] for e in range(3)]
mm = [[widen(M['m'][e][a], M['l'][e].real.sqrt() * epsL) for a in range(2)] for e in range(3)]
# sanity: equivariance (numerical consistency checks; not used for rigour)
print("Gram of (u1,u2) on dE:", [[sum(M['G'][e][a][b] for e in range(3)).mid() for b in range(2)] for a in range(2)])
print("sum_e m^e:", [sum(M['m'][e][a] for e in range(3)) for a in range(2)])

def Ybox(rho, psi):
    """rho, psi: arb balls. Returns (Y 2x2 list, D_A 2x2 list)."""
    cs, sn = psi.cos(), psi.sin()
    a = 2 * cs / sq3; b = 2 * sn / sq3
    d = 1 + rho * b
    H = [[b + rho * a * a / d, -a / d], [-a / d, -b / d]]
    s1 = ((-arb(1) / 2 + rho * cs) ** 2 + (sq3 / 2 + rho * sn) ** 2).sqrt()
    s2 = ((arb(1) / 2 + rho * cs) ** 2 + (sq3 / 2 + rho * sn) ** 2).sqrt()
    d1 = (-cs + sq3 * sn + rho) / (1 + s1)
    d2 = (cs + sq3 * sn + rho) / (1 + s2)
    LT = 1 + s1 + s2
    ce = [-(d1 + d2) / LT, (2 * d1 - d2) / LT, (2 * d2 - d1) / LT]
    v = [d1 * mm[1][x] + d2 * mm[2][x] for x in range(2)]
    se = [arb(1), s1, s2]
    Y = [[None, None], [None, None]]; D = [[None, None], [None, None]]
    sv = [sum(se[e] * mm[e][x] for e in range(3)) for x in range(2)]
    for x in range(2):
        for y in range(2):
            t = sum(ce[e] * g[e][x][y] for e in range(3)) - 3 * rho * v[x] * v[y] / LT ** 2
            Y[x][y] = sigstar * t - sum(H[i][j] * Q[x][y][i][j] for i in range(2) for j in range(2))
            D[x][y] = sum(se[e] * g[e][x][y] for e in range(3)) - sv[x] * sv[y] / LT
    return Y, D

LEAVES = []   # certified leaves (r0, r1, psi0, psi1) as exact floats, for coverage_check.py

def ball(lo, hi):
    """A ball containing the closed interval [lo, hi] (lo, hi exact binary floats)."""
    return arb(lo).union(arb(hi))

def certify_box(r0, r1, p0, p1):
    rho = ball(r0, r1); psi = ball(p0, p1)
    Y, D = Ybox(rho, psi)
    trY = Y[0][0] + Y[1][1]
    detY = Y[0][0] * Y[1][1] - Y[0][1] * Y[1][0]
    Dpd = (D[0][0] > 0) and (D[0][0] * D[1][1] - D[0][1] * D[1][0] > 0)
    return bool(Dpd and ((trY > 0) or (detY < 0)))

def check(r0, r1, p0, p1, depth=0):
    if certify_box(r0, r1, p0, p1):
        LEAVES.append((r0, r1, p0, p1))
        return 1, 0
    if depth > 14:
        return 0, 1
    rm = (r0 + r1) / 2; pm = (p0 + p1) / 2
    ok = 0; bad = 0
    for (a, b_) in [(r0, rm), (rm, r1)]:
        for (c, d_) in [(p0, pm), (pm, p1)]:
            o, bd = check(a, b_, c, d_, depth + 1); ok += o; bad += bd
    return ok, bad

if __name__ == "__main__":
    # first-order coefficient: Y at rho = 0 (eigenvalues +-c up to the traceless check)
    cone = []
    for psi in [0.0, 0.7, 1.9, 4.0]:
        Y, D = Ybox(arb(0), arb(psi))
        tr = Y[0][0] + Y[1][1]; det = Y[0][0] * Y[1][1] - Y[0][1] * Y[1][0]
        lam = (-det).sqrt()
        print("psi=%.2f: trY(0)=%s  sqrt(-detY(0))=%s" % (psi, tr.mid(), lam.mid()))
        cone.append(dict(psi=psi, tr_abs_upper=enc(abs(tr).upper()), eig_lower=enc(lam.lower()), eig_upper=enc(lam.upper())))
    gram = sum(M['G'][e][0][0] for e in range(3))
    rho1 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
    nps = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    PSI_MAX = 6.3            # > 2*pi: the psi-range [0, 6.3] covers a full period
    tot_ok = 0; tot_bad = 0
    for k in range(nps):
        o, b = check(0.0, rho1, PSI_MAX * k / nps, PSI_MAX * (k + 1) / nps)
        tot_ok += o; tot_bad += b
    print("rho1=%r: boxes verified=%d, failed=%d" % (rho1, tot_ok, tot_bad))
    json.dump(dict(rho1=rho1.hex(), psi_max=PSI_MAX.hex(), verified=tot_ok, failed=tot_bad,
                   eps_bdry_upper=enc(epsL.upper()), eps_grad_upper=enc(epsE.upper()),
                   gram_lower=enc(gram.lower()), gram_upper=enc(gram.upper()), cone=cone,
                   leaves=[[x.hex() for x in L] for L in LEAVES]),
              open('../verification/local_leaves.json', 'w'))
    print("LOCAL:", "ALL PASS" if tot_bad == 0 else "FAIL")
    sys.exit(0 if tot_bad == 0 else 1)
