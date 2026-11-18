"""
Far region (paper, Lemma 5.1 and Prop. 5.2): certified by a one-dimensional check.

Parameter domain F = {(p,q): 0 < p <= 1/2, q > 0, (1-p)^2 + q^2 <= 1}; rho = |(p,q) - (1/2, sqrt3/2)|.
For 0 <= p <= 1/2, F ∩ {rho >= rho1} = {0 < q <= Q(p)} with
    Q(p) = min( sqrt(2p - p^2),  sqrt3/2 - sqrt(rho1^2 - (p - 1/2)^2) ).
Trial function u = cos(pi x) on T(p,q):
    energy        I_grad = pi^2 q J(p),  J(p) = j(p) + j(1-p),
                  j(t) = t/4 - sin(2 pi t)/(4 pi) + sin(pi t) sinc(pi t)/(4 pi)  ( = (1/t) int_0^t x sin^2(pi x) dx )
    bdry variance I_bdry = min_c int_{dT} (u - c)^2 >= 2 min_c int_0^1 (cos(pi x) - c)^2 = 1,
so Sigma(T) <= L * I_grad / I_bdry <= (1 + l1 + l2) pi^2 q J(p), which is nondecreasing in q.
On each p-interval [a,b] (b <= 1/2), Q(p) <= Qab := min(sqrt(2b - b^2), sqrt3/2 - sqrt(rho1^2 - (a-1/2)^2)),
l1 <= sqrt(b^2 + Qab^2), l2 <= sqrt((1-a)^2 + Qab^2).  We certify
    pi^2 Qab J([a,b]) (1 + sqrt(b^2+Qab^2) + sqrt((1-a)^2+Qab^2)) < S_lo
on the 16 intervals [k/32, (k+1)/32], k = 0..15 (exact binary endpoints), in Arb ball arithmetic.

`bound(p, q)` (the full closed-form Rayleigh quotient of cos(pi x), including the exact boundary
variance) is kept for the consistency checks; the proof uses only the inequality above.
"""
import sys, json
from flint import arb, ctx
ctx.prec = 64
sys.path.insert(0, '.')
from rig import enc, dec
eq = json.load(open('../verification/equilateral.json'))
S_lo = dec(eq['Sigma_lo'])              # certified lower bound for Sigma(E) (equilateral.py)
RHO1 = 0.6                 # the binary double nearest 0.6 (just below it); same value as local.py
rho1 = arb(RHO1)
pi = arb.pi()
sq3h = arb(3).sqrt() / 2

def sinc(x):
    return x.sinc()

def nsqrt(x):
    """sqrt of a ball whose true value is >= 0 (clamps the negative part of the ball)."""
    if x.lower() >= 0:
        return x.sqrt()
    U = arb(x.upper()).sqrt()
    return arb(0).union(arb(U.upper()))

def j(t):
    return t / 4 - (2 * pi * t).sin() / (4 * pi) + (pi * t).sin() * sinc(pi * t) / (4 * pi)

def bound(p, q):
    """Full closed-form Rayleigh quotient times L for u = cos(pi x) (consistency checks only)."""
    N = pi ** 2 * q * (j(p) + j(1 - p))
    l1 = nsqrt(p * p + q * q); l2 = nsqrt((1 - p) ** 2 + q * q)
    S = arb(1) / 2 + l1 * (arb(1) / 2 + sinc(2 * pi * p) / 2) + l2 * (arb(1) / 2 - (2 * pi * p).sin() / (4 * pi * (1 - p)))
    M = l1 * sinc(pi * p) - l2 * (pi * p).sin() / (pi * (1 - p))
    L = 1 + l1 + l2
    D = S - M * M / L
    return L * N / D, D

def interval_bound(a, b):
    """Certified upper bound for Sigma(T) over {a <= p <= b, 0 < q <= Q(p)} (a, b exact floats)."""
    A, B = arb(a), arb(b)
    U1 = nsqrt(2 * B - B * B)
    U2 = sq3h - nsqrt(rho1 * rho1 - (A - arb(1) / 2) ** 2)
    Qab = arb(min(U1.upper(), U2.upper()))          # an exact float >= Q(p) on [a, b]
    P = A.union(B)                                   # ball containing [a, b]
    J = j(P) + j(1 - P)
    G = pi ** 2 * Qab * J * (1 + (B * B + Qab * Qab).sqrt() + ((1 - A) ** 2 + Qab * Qab).sqrt())
    return G, Qab

if __name__ == "__main__":
    print("S_lo =", S_lo)
    worst = 0.0; ok = True; worst_arb = arb(0)
    n = 16
    for k in range(n):
        a, b = k / 32, (k + 1) / 32
        G, Qab = interval_bound(a, b)
        good = bool(G < S_lo)
        ok &= good
        worst = max(worst, float(G.upper()))
        if G.upper() > worst_arb:
            worst_arb = G.upper()
        print("p in [%-8s %-8s]  Q <= %.6f   Sigma <= %.4f  %s" % (a, b, float(Qab), float(G.upper()), "ok" if good else "FAIL"))
    # the intervals [k/32,(k+1)/32] tile [0,1/2] exactly
    tiles = (0 / 32 == 0.0) and (n / 32 == 0.5)
    print("intervals tile [0, 1/2]:", tiles)
    print("max certified bound on the far region: %.4f" % worst)
    print("FAR REGION:", "ALL PASS" if (ok and tiles) else "FAIL")
    json.dump(dict(ok=bool(ok and tiles), rho1=RHO1.hex(), intervals=n, max_bound_upper=enc(worst_arb)),
              open('../verification/far_region.json', 'w'), indent=1)
    sys.exit(0 if (ok and tiles) else 1)
