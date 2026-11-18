"""
RIGOROUS crude lower bound for sigma_3(E), E = equilateral triangle of side 1 (perimeter 3),
via Crouzeix-Raviart (CR) finite elements and the Liu-type projection argument:

  sigma_k >= sigma_{h,k} / (1 + C_h^2 sigma_{h,k})        (k counts the zero eigenvalue as k=0)

where sigma_{h,k} are the discrete CR Steklov eigenvalues and C_h is a constant with
  || u - Pi_h u ||_{L2(dE)} <= C_h || grad_h (u - Pi_h u) ||_{L2(E)}  for all u in H^1(E),
Pi_h = CR interpolation (edge means).  Derivation of C_h (paper, Lemmas 3.1-3.2):
  for w in H^1(K) with zero mean on the boundary edge e of K (K equilateral, side h):
  ||w||_e <= ||w - wbar_K||_e,  trace identity with the field (x-P)(w-wbar)^2 (P = vertex opposite e),
  Payne-Weinberger ||w - wbar_K||_K <= (h/pi)||grad w||_K  give
     ||w||_e^2 <= (|e|/|K|) h^2 (1/pi^2 + 1/pi) ||grad w||_K^2 = (4/sqrt3)(1/pi^2+1/pi) h ||grad w||_K^2.
  Corner elements carry two boundary edges, so C_h^2 = 2 (4/sqrt3)(1/pi^2+1/pi) h.

Discrete part: we certify that S_h - mu B_h has at most 3 negative eigenvalues (=> sigma_{h,3} >= mu)
by a congruence C = O^T (S_h - mu B_h) O, evaluated in Arb, and Gershgorin discs (Sylvester's law of
inertia; O is certified nonsingular by Gershgorin on O^T O).  The matrix O is COMMITTED certificate data
(../data/cr_congruence.npz, exact binary64 entries, produced once by make_certificate.py), so this
step does not depend on BLAS or on floating-point eigensolvers.  The matrices S_h, B_h are integer
matrices times the exact factors 2/sqrt3 and h/3.  Output: ../verification/cr_lower.json.
"""
import sys, json
import numpy as np
sys.path.insert(0, '.')
from flint import arb, arb_mat, ctx
from rig import enc
ctx.prec = 80


def assemble(m):
    """Uniform mesh of E (side 1) into m^2 equilateral cells; CR stiffness S_h = (2/sqrt3) Sint and
    boundary mass B_h = (h/3) B3 with integer matrices Sint, B3.  Returns (Sint, B3, max boundary
    edges per cell, number of cells)."""
    vid = {}
    for j in range(m + 1):
        for i in range(m + 1 - j):
            vid[(i, j)] = len(vid)
    tris = []
    for j in range(m):
        for i in range(m - j):
            tris.append((vid[(i, j)], vid[(i + 1, j)], vid[(i, j + 1)]))            # up
            if i + j < m - 1:
                tris.append((vid[(i + 1, j)], vid[(i + 1, j + 1)], vid[(i, j + 1)]))  # down
    inv = {v: k for k, v in vid.items()}
    def on_bdry_edge(a, b):
        (i1, j1), (i2, j2) = inv[a], inv[b]
        return (j1 == 0 and j2 == 0) or (i1 == 0 and i2 == 0) or (i1 + j1 == m and i2 + j2 == m)
    edges = {}
    def eid(a, b):
        key = (min(a, b), max(a, b))
        if key not in edges:
            edges[key] = len(edges)
        return edges[key]
    elem = []
    for (a, b, c) in tris:
        # local CR dof i <-> edge opposite vertex i
        elem.append(((a, b, c), (eid(b, c), eid(c, a), eid(a, b))))
    n = len(edges)
    Sint = np.zeros((n, n), dtype=np.int64)
    B3 = np.zeros((n, n), dtype=np.int64)
    nbmax = 0
    loc = np.array([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])
    for (verts, dofs) in elem:
        for x in range(3):
            for y in range(3):
                Sint[dofs[x], dofs[y]] += loc[x, y]
        nb = 0
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            if on_bdry_edge(verts[j], verts[k]):
                nb += 1
                di, dj, dk = dofs[i], dofs[j], dofs[k]
                B3[di, di] += 3
                B3[dj, dj] += 1; B3[dk, dk] += 1
                B3[dj, dk] -= 1; B3[dk, dj] -= 1
        nbmax = max(nbmax, nb)
    return Sint, B3, nbmax, len(elem)


def certify(Sint, B3, m, mu, O):
    """Ball-arithmetic inertia certificate for S_h - mu B_h with congruence matrix O (exact floats)."""
    n = Sint.shape[0]
    sq3 = arb(3).sqrt()
    mua = arb(mu)
    Sa = arb_mat(n, n)
    for x in range(n):
        for y in range(n):
            if Sint[x, y] != 0 or B3[x, y] != 0:
                Sa[x, y] = arb(2) / sq3 * int(Sint[x, y]) - mua * int(B3[x, y]) / (3 * m)
    Oa = arb_mat(O.tolist())
    C = Oa.transpose() * Sa * Oa
    G = Oa.transpose() * Oa
    ok_O = all((G[i, i] - sum(abs(G[i, j]) for j in range(n) if j != i)) > 0 for i in range(n))
    neg = 0; pos = 0
    for i in range(n):
        rad = sum(abs(C[i, j]) for j in range(n) if j != i)
        if (C[i, i] + rad) < 0:
            neg += 1
        elif (C[i, i] - rad) > 0:
            pos += 1
    return neg, pos, ok_O


if __name__ == "__main__":
    cert = np.load('../data/cr_congruence.npz')
    m = int(cert['m']); mu = float(cert['mu']); O = cert['O']
    Sint, B3, nbmax, ncells = assemble(m)
    n = Sint.shape[0]
    assert O.shape == (n, n) and O.dtype == np.float64
    print("n_h=%d  dofs=%d  cells=%d  max boundary edges per cell=%d" % (m, n, ncells, nbmax))
    neg, pos, okO = certify(Sint, B3, m, mu, O)
    print("mu=%r: certified #neg discs=%d, #pos discs=%d (n=%d), O nonsingular=%s" % (mu, neg, pos, n, okO))
    ok = okO and neg + pos == n and neg <= 3 and nbmax <= 2
    Ch2 = 2 * (4 / arb(3).sqrt()) * (1 / arb.pi() ** 2 + 1 / arb.pi()) / m
    lb = arb(mu) / (1 + Ch2 * arb(mu))
    print("=> sigma_{h,3} >= %r and sigma_3(E) >= mu/(1+C_h^2 mu) = %s   (C_h^2 = %s)" % (mu, lb, Ch2))
    json.dump(dict(ok=bool(ok), n_h=m, dofs=n, mu=mu, neg=neg, pos=pos, O_nonsingular=bool(okO),
                   max_bdry_edges_per_cell=nbmax, Ch2_upper=enc(Ch2.upper()),
                   sigma3_lower=enc(lb.lower())),
              open('../verification/cr_lower.json', 'w'), indent=1)
    print("CR LOWER:", "ALL PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
