"""Consistency check (numerical): the Region-I upper bound L_T*lambda_min(N_A,D_A) from the
transplanted equilateral eigenspace must dominate the converged Ritz value sigma_1 L(T)."""
import sys, math; sys.path.insert(0, '.')
import numpy as np
import local as LC
from flint import arb
from hritz import steklov, triangle
sq3 = math.sqrt(3)
print(" rho    psi(deg)   transplant-bound   Ritz sigma1*L   cos(pi x)-bound")
import far_region as R3
bad = 0
for rho in [0.05, 0.2, 0.4, 0.6]:
    for psideg in [0, 90, 210, 240, 270]:
        psi = math.radians(psideg)
        p = 0.5 + rho * math.cos(psi); q = sq3 / 2 + rho * math.sin(psi)
        # recompute N_A, D_A directly (midpoints) from the enclosed moments
        Y, D = LC.Ybox(arb(rho), arb(psi))
        a = 2 * math.cos(psi) / sq3; b = 2 * math.sin(psi) / sq3
        A = np.array([[1, rho * a], [0, 1 + rho * b]]); Ai = np.linalg.inv(A); G = Ai @ Ai.T
        Q = np.array([[[[float(LC.Q[x][y][i][j].mid()) for j in range(2)] for i in range(2)] for y in range(2)] for x in range(2)])
        N = abs(np.linalg.det(A)) * np.einsum('ij,xyij->xy', G, Q)
        Dm = np.array([[float(D[x][y].mid()) for y in range(2)] for x in range(2)])
        s1 = math.hypot(p - 1, q); s2 = math.hypot(p, q); LT = 1 + s1 + s2
        lam = min(np.linalg.eigvals(np.linalg.solve(Dm, N)).real)
        ev, L = steklov(triangle(p, q), 30)
        r3 = float(R3.bound(arb(p), arb(q))[0].mid())
        print("%5.2f %6d      %12.6f      %12.6f    %12.6f" % (rho, psideg, LT * lam, ev[0] * L, r3))
        if not (LT * lam >= ev[0] * L - 1e-9 and r3 >= ev[0] * L - 1e-9 and LT * lam < 3.8724615957722):
            bad += 1
import json
json.dump(dict(ok=bad == 0, failures=bad), open("../verification/check_local_vs_ritz.json", "w"), indent=1)
print("CHECK_LOCAL_VS_RITZ:", "ALL PASS" if bad == 0 else "FAIL (%d)" % bad)
sys.exit(0 if bad == 0 else 1)
