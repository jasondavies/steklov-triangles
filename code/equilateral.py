"""
Paper, Section 3.1 and Prop. 3.5: certified enclosure of sigma_1 of the equilateral triangle E
(side 1, perimeter 3).  Reads ONLY the committed certificate ../data/equilateral_pair.json
(coefficients c_k as exact binary64 numbers); no floating-point linear algebra is used.

1. (u1, u2) = (Re f, Im f), f = sum_{k=1 mod 3} c_k w^k + sum_{k=2 mod 3} c_k wbar^k, w = -i sqrt3 (z - z0):
   an EXACTLY D3-equivariant pair of harmonic polynomials, whatever the c_k.
2. sigma~ = the Rayleigh quotient of u1 (ball arithmetic), rounded to binary64 (deterministic).
3. eta = max over span(u1,u2) of ||d_n f - sigma~ f|| / ||f||  (L2(dE)), from exact 2x2 Gram matrices.
By the inclusion lemma, >= 2 eigenvalues of the DtN map lie in [sigma~ - eta, sigma~ + eta]; with
sigma_3(E) > sigma~ + eta (cr_lower.py) these are sigma_1 = sigma_2.
Output (machine readable, exact 'man*2^exp' encodings): ../verification/equilateral.json
"""
import json, sys
sys.path.insert(0, '.')
from flint import acb, arb, ctx
from rig import residual_sq, moments, residual_gram, enc

ctx.prec = 200
cert = json.load(open('../data/equilateral_pair.json'))
n = cert['n']; ks = cert['ks']
alpha = [float.fromhex(h) for h in cert['c_hex']]

# complex coefficient lists in the variable (z - z0)/R, R = 1/sqrt3: coefficient of ((z-z0)/R)^k
c1 = [0j] * n; c2 = [0j] * n
for a, k in zip(alpha, ks):
    c1[k - 1] = a * (-1j) ** k
    c2[k - 1] = ((-1j * a) if k % 3 == 1 else (1j * a)) * (-1j) ** k
# (-1j)**k and the products are exact in binary64 (multiplication by +-1, +-i)

sq3 = arb(3).sqrt()
Zb = [acb(0), acb(1), acb(arb(1) / 2, sq3 / 2)]
zc = acb(arb(1) / 2, sq3 / 6)
Rb = 1 / sq3
coefb1 = [acb(arb(x.real), arb(x.imag)) for x in c1]
coefb2 = [acb(arb(x.real), arb(x.imag)) for x in c2]
funcs = [(coefb1, zc, Rb), (coefb2, zc, Rb)]
M = moments(funcs, Zb)
E11 = M['Q'][0][0][0][0] + M['Q'][0][0][1][1]
G11 = sum(M['G'][e][0][0] for e in range(3))
rq = E11 / G11
sig = arb(float(rq.mid()))           # exact binary64 shift
print("Rayleigh quotient of u1 (ball):", rq, "  sigma~ =", float(sig), "=", float(sig).hex())
print("Gram of (u1,u2):", [[sum(M['G'][e][a][b] for e in range(3)).mid() for b in range(2)] for a in range(2)])
print("boundary means:", [sum(M['m'][e][a] for e in range(3)) for a in range(2)])
Rg, Ug = residual_gram(funcs, Zb, sig)
a11, a12, a22 = Rg[0][0], Rg[0][1], Rg[1][1]
b11, b12, b22 = Ug[0][0], Ug[0][1], Ug[1][1]
detB = b11 * b22 - b12 * b12
assert detB > 0
tr = (a11 * b22 + a22 * b11 - 2 * a12 * b12) / detB
de = (a11 * a22 - a12 * a12) / detB
disc = tr * tr - 4 * de
lam_max_up = (tr + arb(max(disc.upper(), arb(0))).sqrt()).upper() / 2      # >= largest eigenvalue
eta = arb(lam_max_up).sqrt().upper()                                       # exact upper bound of eta
lo = (sig - eta).lower(); hi = (sig + eta).upper()
SigE_lo = (3 * lo).lower(); SigE_hi = (3 * hi).upper()
f = (sq3 / 4).sqrt() / 3                                                    # sqrt|E| / L(E)
area_lo = (SigE_lo * f).lower(); area_hi = (SigE_hi * f).upper()
print("eta <=", eta)
print("sigma_1(E) in [%s, %s]" % (lo, hi))
print("Sigma(E) = sigma_1 L in [%s, %s]" % (SigE_lo, SigE_hi))
print("sigma_1 sqrt|E| in [%s, %s]" % (area_lo, area_hi))
gram = sum(M['G'][e][0][0] for e in range(3))
json.dump(dict(n=n, sigma_tilde_hex=float(sig).hex(), eta_upper=enc(eta),
               sigma_lo=enc(lo), sigma_hi=enc(hi), Sigma_lo=enc(SigE_lo), Sigma_hi=enc(SigE_hi),
               area_lo=enc(area_lo), area_hi=enc(area_hi), gram_u1_lower=enc(gram.lower()), gram_u1_upper=enc(gram.upper())),
          open('../verification/equilateral.json', 'w'), indent=1)
