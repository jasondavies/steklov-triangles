"""
Numerical table of the paper (Section 7), recomputed carefully.  NOT part of the proof.

For each triangle (and the regular N-gons) we run the harmonic-polynomial Ritz method with the basis
{Re w^k, Im w^k : 1 <= k <= n}, w = (z - c)/R (c, R exact binary64), in high precision:
  * the Ritz matrices (energy, boundary Gram, boundary means) are exact polynomial line integrals,
    evaluated in Arb at 400 bits (rig.moments);
  * the generalised eigenproblem (variance form, no truncation of directions) is solved in mpmath at
    `dps` digits; the retained rank (= 2n, full) and the smallest eigenvalue of the variance Gram are
    recorded;
  * convergence is checked in the degree n (several values) and in the precision (two dps values);
  * CERTIFIED UPPER BOUND: the computed Ritz vector, rounded to binary64 coefficients, is a fixed
    trial function; its Rayleigh quotient on the EXACT polygon (vertices as Arb balls) is enclosed in
    ball arithmetic, which gives a rigorous upper bound for Sigma = sigma_1 L.
The estimate is the highest-degree Ritz value; the error estimate is max(last degree change, geometric
tail estimate from the last two changes); Ritz values decrease monotonically to Sigma.
OPTIONAL and slow; NOT run by scripts/verify.py.  Output: ../data/table_estimates.json, including the
binary64 trial vectors, from which code/table_bounds.py (run by verify.py) recomputes the certified
upper bounds.
"""
import sys, json, math
sys.path.insert(0, '.')
import mpmath as mp
from flint import acb, arb, ctx
from rig import moments, enc

ctx.prec = 400

# (name, angle at (0,0), angle at (1,0)) in degrees; the paper's table
TRIANGLES = [("60-60-60", 60, 60), ("50-50-80", 50, 50), ("70-70-40", 70, 70), ("40-60-80", 40, 60),
             ("45-45-90", 45, 45), ("30-60-90", 60, 30), ("30-30-120", 30, 30)]
NGONS = (4, 5, 6)

def polygon(name):
    for nm, A, B in TRIANGLES:
        if nm == name:
            return tri_vertices(A, B)
    N = int(name.split()[1].split('-')[0])
    return ngon_vertices(N)

def tri_vertices(A, B):
    """Exact (ball) vertices (0,0), (1,0), apex for angles A at (0,0), B at (1,0), in degrees."""
    a = arb.pi() * A / 180; b = arb.pi() * B / 180; c = arb.pi() - a - b
    side = b.sin() / c.sin()
    return [acb(0), acb(1), acb(side * a.cos(), side * a.sin())]

def ngon_vertices(N):
    return [acb((2 * arb.pi() * k / N).cos(), (2 * arb.pi() * k / N).sin()) for k in range(N)]

def basis_coeffs(n):
    """Basis functions Re(w^k) and Im(w^k) = Re(-i w^k) as coefficient lists."""
    fs = []
    for k in range(1, n + 1):
        for cval in (acb(1), acb(0, -1)):
            cc = [acb(0)] * n; cc[k - 1] = cval; fs.append(cc)
    return fs

def ritz_data(Z, n, cz, R):
    """Exact (ball) Ritz data: energy matrix, boundary Gram, boundary means, perimeter."""
    fs = basis_coeffs(n)
    M = moments([(cc, cz, R) for cc in fs], Z)
    m = len(fs); ne = len(Z)
    K = [[M['Q'][i][j][0][0] + M['Q'][i][j][1][1] for j in range(m)] for i in range(m)]
    G = [[sum(M['G'][e][i][j] for e in range(ne)) for j in range(m)] for i in range(m)]
    b = [sum(M['m'][e][i] for e in range(ne)) for i in range(m)]
    L = sum(M['l'][e].real for e in range(ne))
    return K, G, b, L, m

def ritz_solve(data, dps):
    K_, G_, b_, L_, m = data
    mp.mp.dps = dps
    tomp = lambda x: mp.mpf(x.mid().str(dps + 10, radius=False))
    Lm = tomp(L_); bv = [tomp(x) for x in b_]
    K = mp.matrix(m, m); B0 = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            K[i, j] = tomp(K_[i][j])
            B0[i, j] = tomp(G_[i][j]) - bv[i] * bv[j] / Lm
    ev_b, _ = mp.eigsy(B0)
    Li = mp.inverse(mp.cholesky(B0))
    Kt = Li * K * Li.T
    Kt = (Kt + Kt.T) / 2
    ev, V = mp.eigsy(Kt)
    i0 = min(range(m), key=lambda i: ev[i])
    x = Li.T * V[:, i0]
    return ev[i0] * Lm, x, min(ev_b), m

def certified_upper(Z, n, x, cz, R):
    """Rigorous upper bound for Sigma from the fixed trial function sum x_i phi_i (binary64 coefficients)."""
    xs = [float(xi) for xi in x]
    cc = [acb(0)] * n
    for k in range(1, n + 1):
        cc[k - 1] = acb(arb(xs[2 * (k - 1)]), arb(-xs[2 * (k - 1) + 1]))     # Re((a - i b) w^k)
    M = moments([(cc, cz, R)], Z)
    ne = len(Z)
    E = M['Q'][0][0][0][0] + M['Q'][0][0][1][1]
    G = sum(M['G'][e][0][0] for e in range(ne)); m = sum(M['m'][e][0] for e in range(ne))
    L = sum(M['l'][e].real for e in range(ne))
    D = G - m * m / L
    assert D > 0
    return (L * E / D).upper()

def run(name, Z, degrees, dps_list):
    zc = sum(Z, acb(0)) / len(Z)
    cz = acb(float(zc.real.mid()), float(zc.imag.mid()))
    R = arb(float(max(abs(z - cz).mid() for z in Z)))
    rows = []; vals = {}; xlast = None
    for n in degrees:
        data = ritz_data(Z, n, cz, R)
        for dps in dps_list:
            val, x, bmin, m = ritz_solve(data, dps)
            vals[(n, dps)] = val; xlast = x
            rows.append(dict(n=n, dps=dps, value=mp.nstr(val, 20), rank=m, basis_size=m, min_gram_eig=mp.nstr(bmin, 5)))
            print("  %-16s n=%2d dps=%2d  Sigma ~ %s  (rank %d of %d; min eigenvalue of variance Gram %s)" %
                  (name, n, dps, mp.nstr(val, 16), m, m, mp.nstr(bmin, 3)), flush=True)
    nmax = degrees[-1]; dmax = dps_list[-1]
    up = certified_upper(Z, nmax, xlast, cz, R)
    est = vals[(nmax, dmax)]
    d1 = abs(vals[(degrees[-2], dmax)] - est)                 # last degree change
    d0 = abs(vals[(degrees[-3], dmax)] - vals[(degrees[-2], dmax)]) if len(degrees) >= 3 else None
    r = (d1 / d0) if (d0 and d0 > 0) else mp.mpf(0)
    tail = d1 * r / (1 - r) if r < 1 else mp.inf                # geometric tail estimate of the remaining error
    err = max(d1, tail)
    err_prec = abs(vals[(nmax, dps_list[0])] - est)
    print("  %-16s estimate %s +- %s (degree changes %s, %s; precision change %s); certified upper bound %s" %
          (name, mp.nstr(est, 16), mp.nstr(err, 2), mp.nstr(d0, 2) if d0 is not None else '-', mp.nstr(d1, 2),
           mp.nstr(err_prec, 2), arb(up).str(16, radius=False)), flush=True)
    return dict(name=name, rows=rows, estimate=mp.nstr(est, 20), err_estimate=mp.nstr(err, 3),
                err_precision=mp.nstr(err_prec, 3), certified_upper=enc(up),
                certified_upper_str=arb(up).str(16, radius=False),
                trial=dict(n=nmax, center_hex=[float(cz.real.mid()).hex(), float(cz.imag.mid()).hex()],
                           R_hex=float(R.mid()).hex(), x_hex=[float(xi).hex() for xi in xlast]))

if __name__ == "__main__":
    import time, argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--degrees', default='30,40,50')
    ap.add_argument('--only', default='', help='comma-separated polygon names (default: all)')
    ap.add_argument('--out', default='../data/table_estimates.json')
    ap.add_argument('--merge', nargs='*', help='merge partial outputs into --out (in table order)')
    args = ap.parse_args()
    names = [t[0] for t in TRIANGLES] + ["regular %d-gon" % N for N in NGONS]
    if args.merge:
        parts = [json.load(open(f)) for f in args.merge]
        polys = {p['name']: p for part in parts for p in part['polygons']}
        json.dump(dict(degrees=parts[0]['degrees'], dps=parts[0]['dps'], seconds=sum(pt['seconds'] for pt in parts),
                       polygons=[polys[n] for n in names if n in polys]), open(args.out, 'w'), indent=1)
        sys.exit(0)
    t0 = time.time()
    degrees = [int(d) for d in args.degrees.split(',')]
    dps_list = [50, 80]
    want = set(args.only.split(',')) if args.only else set(names)
    out = []
    for name in names:
        if name in want:
            out.append(run(name, polygon(name), degrees, dps_list))
    json.dump(dict(degrees=degrees, dps=dps_list, seconds=round(time.time() - t0), polygons=out),
              open(args.out, 'w'), indent=1)
