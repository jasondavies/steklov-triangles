"""
Certified upper bounds for the numerical table (paper, Section 7).  Run by scripts/verify.py.

Reads the committed trial functions in ../data/table_estimates.json (binary64 coefficients, produced
by the optional, slow hp_table.py) and encloses, in Arb ball arithmetic on the EXACT polygon, the
Rayleigh quotient times the perimeter.  By the variational principle this is a rigorous upper bound for
Sigma = sigma_1 L.  Output: ../verification/table_bounds.json.
"""
import sys, json
sys.path.insert(0, '.')
from flint import acb, arb, ctx
from rig import enc
import hp_table as H
import mpmath as mp
mp.mp.dps = 50

ctx.prec = 300
data = json.load(open('../data/table_estimates.json'))
out = []
for poly in data['polygons']:
    t = poly['trial']
    Z = H.polygon(poly['name'])
    cz = acb(float.fromhex(t['center_hex'][0]), float.fromhex(t['center_hex'][1]))
    R = arb(float.fromhex(t['R_hex']))
    x = [float.fromhex(h) for h in t['x_hex']]
    up = H.certified_upper(Z, t['n'], x, cz, R)
    # error estimate (heuristic, not rigorous): the largest of the last change in n, the geometric tail
    # (hp_table.py), a power-law tail fitted to the last two changes, and the precision change
    vals = {r['n']: mp.mpf(r['value']) for r in poly['rows'] if r['dps'] == max(rr['dps'] for rr in poly['rows'])}
    ns = sorted(vals)
    d0 = abs(vals[ns[-3]] - vals[ns[-2]]); d1 = abs(vals[ns[-2]] - vals[ns[-1]])
    alg = mp.mpf(0)
    if d0 > 0 and d1 > 0 and 0.01 <= d1 / d0 < 1:     # algebraic convergence (fast cases: geometric tail)
        f = lambda q: (mp.mpf(ns[-2]) ** -q - mp.mpf(ns[-1]) ** -q) / (mp.mpf(ns[-3]) ** -q - mp.mpf(ns[-2]) ** -q) - d1 / d0
        try:
            q = mp.findroot(f, 3)
            Cc = d1 / (mp.mpf(ns[-2]) ** -q - mp.mpf(ns[-1]) ** -q)
            alg = abs(Cc) * mp.mpf(ns[-1]) ** -q
        except Exception:
            alg = mp.inf
    err_total = max(mp.mpf(poly['err_estimate']), mp.mpf(poly['err_precision']), alg)
    print("%-16s certified upper bound %s   estimate %s +- %s" % (poly['name'], arb(up).str(15, radius=False),
                                                              poly['estimate'][:16], mp.nstr(err_total, 2)))
    out.append(dict(name=poly['name'], certified_upper=enc(up), estimate=poly['estimate'], err_estimate=poly['err_estimate'],
                    err_precision=poly['err_precision'], err_power_tail=mp.nstr(alg, 3), err_total=mp.nstr(err_total, 3),
                    basis_size=poly['rows'][-1]['basis_size'], retained_rank=poly['rows'][-1]['rank'],
                    min_gram_eig=poly['rows'][-1]['min_gram_eig'], dps=poly['rows'][-1]['dps'], n=t['n']))
json.dump(out, open('../verification/table_bounds.json', 'w'), indent=1)
print("TABLE BOUNDS: done (%d polygons)" % len(out))
