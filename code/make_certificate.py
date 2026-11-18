"""
OPTIONAL: regenerate the committed certificate data in ../data/.  NOT run by scripts/verify.py.

The certified steps read only these files, so their results do not depend on BLAS, on the CPU, or on
floating-point eigensolvers.  Regenerating with a different BLAS may give slightly different (equally
valid) data; the claims in ../claims.json must then be re-checked by scripts/verify.py.

  ../data/equilateral_pair.json  coefficients c_k (exact binary64, hex) of the D3-equivariant pair
                                 (u1, u2) on E, degree n = 20
  ../data/cr_congruence.npz      congruence matrix O (binary64) for S_h - mu B_h, n_h = 16, mu = 3
"""
import sys, json
import numpy as np
sys.path.insert(0, '.')
from hritz import poly_data

def equilateral_pair(n=20):
    s3 = np.sqrt(3) / 2
    P = np.array([[0, 0], [1, 0], [0.5, s3]])
    c = complex(0.5, s3 / 3)
    R = 1 / np.sqrt(3)
    ks = [k for k in range(1, n + 1) if k % 3 != 0]
    Z, Wt, N, L, _, _ = poly_data(P, n, center=c)
    w = -1j * (Z - c) / R
    dw = -1j / R
    Phi = np.array([(w ** k).real for k in ks]).T
    dPhi = np.array([((k * w ** (k - 1) * dw) * N).real for k in ks]).T
    K = Phi.T @ (Wt[:, None] * dPhi); K = 0.5 * (K + K.T)
    B = Phi.T @ (Wt[:, None] * Phi)
    Lc = np.linalg.cholesky(B)
    Kt = np.linalg.solve(Lc, np.linalg.solve(Lc, K).T)
    ev, V = np.linalg.eigh(0.5 * (Kt + Kt.T))
    alpha = np.linalg.solve(Lc.T, V[:, 0])
    alpha = alpha / np.max(np.abs(alpha))
    print("float Ritz sigma_1 L(E) =", 3 * ev[0])
    return dict(n=n, ks=ks, c_hex=[float(a).hex() for a in alpha])

def cr_congruence(m=16, mu=3.0):
    from cr_lower import assemble
    Sint, B3, _, _ = assemble(m)
    S = (2 / np.sqrt(3)) * Sint - mu * (1.0 / (3 * m)) * B3
    w, O = np.linalg.eigh(S)
    return O

if __name__ == "__main__":
    pair = equilateral_pair()
    json.dump(pair, open('../data/equilateral_pair.json', 'w'), indent=1)
    O = cr_congruence()
    np.savez_compressed('../data/cr_congruence.npz', O=O, m=16, mu=3.0)
    print("wrote ../data/equilateral_pair.json and ../data/cr_congruence.npz")
