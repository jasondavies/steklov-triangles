"""NUMERICAL ONLY (floating point, not part of the proof): thin-triangle asymptotics, the explicit
eigenvalue 2 sqrt3 of E, and the Crouzeix-Raviart discrete values quoted in the paper.
Output: ../verification/numerics_extra.json (written last)."""
import sys, math, json; sys.path.insert(0, '.')
import numpy as np, scipy.linalg as sla
from hritz import steklov, triangle
from cr_lower import assemble
out = {}

# thin triangles: Sigma ~ 2 q theta(p), theta(p) = min int Phi_p u'^2 / (2 int u^2), int u = 0
def theta(p, M=4000):
    x = np.linspace(0, 1, M + 1); h = 1 / M
    xm = 0.5 * (x[:-1] + x[1:])
    phi = (1 - xm) if p == 0 else np.where(xm <= p, xm / p, (1 - xm) / (1 - p))
    K = np.zeros((M + 1, M + 1)); Mm = np.zeros((M + 1, M + 1))
    for i in range(M):
        K[i:i + 2, i:i + 2] += phi[i] / h * np.array([[1, -1], [-1, 1]])
        Mm[i:i + 2, i:i + 2] += 2 * h / 6 * np.array([[2, 1], [1, 2]])
    return sla.eigh(K, Mm, eigvals_only=True, subset_by_index=[0, 1])[1]
th0 = theta(0.0)
print("theta(0): exact j_{1,1}^2/8 = %.6f, numeric %.6f" % (3.8317059702075125 ** 2 / 8, th0))
deg = []
for p in [0.05, 0.25, 0.5]:
    t = theta(p)
    ratios = []
    for q in [0.04, 0.02, 0.01]:
        ev, L = steklov(triangle(p, q), 40); ratios.append(ev[0] * L / (2 * q))
    print("p=%.2f  theta(p)=%.6f   Sigma/(2q) at q=0.04,0.02,0.01: %s" % (p, t, " ".join("%.6f" % r for r in ratios)))
    deg.append(dict(p=p, theta=t, ratios=ratios))
out["theta0"] = th0
out["thin"] = deg

# sigma = 2 sqrt3 on E (side 1): explicit eigenfunction u = Im((z-z0)^3) + 1/(12 sqrt3)
r = 1 / (2 * math.sqrt(3)); V3 = [complex(-0.5, -r), complex(0.5, -r), complex(0, 2 * r)]
res = 0.0
for k in range(3):
    a_, b_ = V3[k], V3[(k + 1) % 3]; nrm = -1j * (b_ - a_) / abs(b_ - a_)
    for t in np.linspace(0, 1, 11):
        z = a_ + t * (b_ - a_)
        res = max(res, abs((3 * z ** 2 * nrm).imag - 2 * math.sqrt(3) * ((z ** 3).imag + 1 / (12 * math.sqrt(3)))))
evE, LE = steklov(triangle(0.5, math.sqrt(3) / 2), 30, nev=4)
print("explicit eigenfunction: max |d_n u - 2sqrt3 u| on dE = %.2e" % res)
print("Ritz sigma_3(E) = %.15f, 2 sqrt3 = %.15f, difference %.1e" % (evE[2], 2 * math.sqrt(3), abs(evE[2] - 2 * math.sqrt(3))))
out["explicit_residual"] = res
out["sigma3_ritz"] = float(evE[2])

# Crouzeix-Raviart discrete eigenvalues sigma_{h,1} (paper, consistency checks)
cr = {}
for m in (8, 16):
    Sint, B3, _, _ = assemble(m)
    ev = sla.eigvals((2 / np.sqrt(3)) * Sint, (1.0 / (3 * m)) * B3)
    ev = np.sort(ev[np.isfinite(ev)].real); ev = ev[ev > -1e-8]
    cr[m] = float(ev[1])
    print("CR n_h=%d: sigma_{h,1} = %.6f" % (m, ev[1]))
out["cr_sigma_h1"] = {str(k): v for k, v in cr.items()}

okexp = res < 1e-12 and abs(evE[2] - 2 * math.sqrt(3)) < 1e-9
okdeg = all(abs(d['ratios'][-1] - d['theta']) / d['theta'] < 1e-3 for d in deg) and abs(th0 - 3.8317059702075125 ** 2 / 8) < 1e-6
okcr = cr[8] < cr[16] < evE[0]
out["ok"] = bool(okexp and okdeg and okcr)
json.dump(out, open('../verification/numerics_extra.json', 'w'), indent=1)
print("NUMERICS_EXTRA:", "ALL PASS" if out["ok"] else "FAIL")
sys.exit(0 if out["ok"] else 1)
