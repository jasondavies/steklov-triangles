"""
Harmonic-polynomial Ritz (Galerkin) method for Steklov eigenvalues of a polygon.
NUMERICAL (float64).  Trial space: span{Re w^k, Im w^k : 1<=k<=n}, w=(z-c)/R.
Since every basis function is harmonic, the energy form is a boundary integral,
  E(u,v) = int_{dOmega} u d_n v ds,
and the Steklov boundary form is B(u,v) = int_{dOmega} u v ds; the mean-zero
constraint is handled with the variance form B0 = B - b b^T / L.
Ritz values are upper bounds for sigma_k (min-max), and converge from above.
FLOAT64 ONLY: directions of the variance Gram with relative eigenvalue < 1e-14 are dropped, so for
large n the effective trial space is smaller and values can be too high (e.g. 30-30-120 at n = 30 keeps
50 of 60 directions and is ~2e-6 too high).  Accurate values: hp_table.py; certified bounds: table_bounds.py.
"""
import numpy as np
import scipy.linalg as sla

def triangle(p, q):
    return np.array([[0.0, 0.0], [1.0, 0.0], [p, q]])

def poly_data(P, n, m=None, center=None):
    """Boundary quadrature data for a polygon P (k x 2, ccw)."""
    P = np.asarray(P, float)
    k = len(P)
    if m is None:
        m = n + 2
    xg, wg = np.polynomial.legendre.leggauss(m)
    xg = 0.5 * (xg + 1); wg = 0.5 * wg
    Z = []; Wt = []; N = []
    L = 0.0
    for i in range(k):
        a = complex(*P[i]); b = complex(*P[(i + 1) % k])
        ell = abs(b - a); L += ell
        t = (b - a) / ell
        nrm = -1j * t           # outward normal for ccw polygon
        Z.append(a + xg * (b - a)); Wt.append(wg * ell); N.append(np.full(m, nrm))
    Z = np.concatenate(Z); Wt = np.concatenate(Wt); N = np.concatenate(N)
    if center is None:
        center = complex(*P.mean(axis=0))
    R = np.max(np.abs(P[:, 0] + 1j * P[:, 1] - center))
    return Z, Wt, N, L, center, R

def basis(Z, N, n, center, R):
    w = (Z - center) / R
    cols = []; dcols = []
    for k in range(1, n + 1):
        f = w ** k
        fp = k * w ** (k - 1) / R
        dn = fp * N
        cols += [f.real, f.imag]
        dcols += [dn.real, dn.imag]
    return np.array(cols).T, np.array(dcols).T

def steklov(P, n=20, nev=4, return_vec=False, center=None):
    Z, Wt, N, L, c, R = poly_data(P, n, center=center)
    Phi, dPhi = basis(Z, N, n, c, R)
    K = Phi.T @ (Wt[:, None] * dPhi); K = 0.5 * (K + K.T)
    B = Phi.T @ (Wt[:, None] * Phi)
    b = Phi.T @ Wt
    B0 = B - np.outer(b, b) / L
    # stabilise: orthonormalise w.r.t. B0
    s, U = np.linalg.eigh(B0)
    keep = s > s.max() * 1e-14
    T = U[:, keep] / np.sqrt(s[keep])
    Kt = T.T @ K @ T
    ev, V = np.linalg.eigh(0.5 * (Kt + Kt.T))
    if return_vec:
        return ev[:nev], T @ V[:, :nev], L, (Z, Wt, N, c, R)
    return ev[:nev], L

if __name__ == "__main__":
    s3 = np.sqrt(3) / 2
    for n in [8, 12, 16, 20, 24, 30, 36]:
        ev, L = steklov(triangle(0.5, s3), n)
        print("equilateral n=%2d  sigma*L:" % n, " ".join("%.12f" % (e * L) for e in ev))
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    for n in [10, 20, 30]:
        ev, L = steklov(sq, n)
        print("square[-1,1]^2 n=%d sigma:" % n, ev)
