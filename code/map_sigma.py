"""NUMERICAL map of sigma_1*L and sigma_1*sqrt(A) over triangles (apex (p,q), base [0,1]),
fundamental region 0<=p<=1/2, (1-p)^2+q^2<=1 (base = longest side, |v0 v2| <= |v1 v2|).
Harmonic-polynomial Ritz, degree n=30 (converged to ~1e-8 except very thin triangles)."""
import sys, json
import numpy as np
sys.path.insert(0, '.')
from hritz import steklov, triangle
out = []
for p in np.linspace(0.02, 0.5, 25):
    qmax = np.sqrt(1 - (1 - p) ** 2)
    for q in np.linspace(0.01, qmax, 30):
        ev, L = steklov(triangle(p, q), 30)
        A = q / 2
        out.append(dict(p=float(p), q=float(q), sL=float(ev[0] * L), s2L=float(ev[1] * L), sA=float(ev[0] * np.sqrt(A))))
json.dump(out, open('../verification/map_sigma.json', 'w'))
best = max(out, key=lambda d: d['sL']); bestA = max(out, key=lambda d: d['sA'])
print("max sigma1 L on grid:", best); print("max sigma1 sqrt(A):", bestA)
E = steklov(triangle(0.5, np.sqrt(3) / 2), 30)
print("equilateral sigma1 L = %.12f, sigma1 sqrt(A) = %.12f" % (E[0][0] * E[1], E[0][0] * np.sqrt(np.sqrt(3) / 4)))
ok = abs(best['p'] - 0.5) < 1e-12 and abs(best['q'] - np.sqrt(3) / 2) < 1e-12 and \
     all(d['sL'] < E[0][0] * E[1] - 1e-6 for d in out if abs(d['q'] - np.sqrt(3) / 2) > 1e-9 or abs(d['p'] - 0.5) > 1e-9)
json.dump(dict(ok=bool(ok), grid_points=len(out), argmax=[best["p"], best["q"]], max_Sigma=best["sL"]), open("../verification/map_sigma_summary.json", "w"), indent=1)
print("MAP:", "ALL PASS" if ok else "FAIL", "(grid maximum at the equilateral apex, all other grid values smaller)")
sys.exit(0 if ok else 1)
