"""
Machine check that the local covering of local.py really covers the closed disc rho <= rho1
(paper, proof of Prop. 4.4).

Input: ../verification/local_leaves.json, written by local.py.  Every leaf endpoint is an exact
binary float (stored as float.hex), so float comparisons are exact and Fraction(float) is exact.

 L1  the leaves lie in [0,RHO1] x [0,PSI_MAX], have pairwise disjoint interiors and their areas sum
     EXACTLY to RHO1*PSI_MAX; hence they cover that rectangle (a finite union of closed boxes of full
     measure in a box is the whole box).  PSI_MAX > 2 pi (ball arithmetic), so the polar map
     (rho, psi) -> (1/2 + rho cos psi, sqrt3/2 + rho sin psi) sends the covered set ONTO the closed disc
     of radius RHO1 about the equilateral apex.
 L2  every leaf is re-certified here by local.certify_box.
 L3  RHO1 is the same binary double as in far_region.py, so the local and far regions meet exactly.
"""
import sys, json
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, '.')
from flint import arb

ok_all = True
def report(name, ok, extra=''):
    global ok_all
    ok_all &= bool(ok)
    print('%-4s %s %s' % (name, 'PASS' if ok else 'FAIL', extra))

def tiles(boxes, X0, X1, Y0, Y1):
    a = np.array(boxes, dtype=float)
    inside = bool(np.all(a[:, 0] >= X0) and np.all(a[:, 1] <= X1) and np.all(a[:, 2] >= Y0) and np.all(a[:, 3] <= Y1)
                  and np.all(a[:, 0] < a[:, 1]) and np.all(a[:, 2] < a[:, 3]))
    overlap = 0
    for i in range(len(a)):
        b = a[i + 1:]
        ov = (np.minimum(a[i, 1], b[:, 1]) > np.maximum(a[i, 0], b[:, 0])) & \
             (np.minimum(a[i, 3], b[:, 3]) > np.maximum(a[i, 2], b[:, 2]))
        overlap += int(ov.sum())
    area = sum((Fr(x1) - Fr(x0)) * (Fr(y1) - Fr(y0)) for (x0, x1, y0, y1) in boxes)
    return inside, overlap, area == (Fr(X1) - Fr(X0)) * (Fr(Y1) - Fr(Y0)), len(a)

LCd = json.load(open('../verification/local_leaves.json'))
RHO1 = float.fromhex(LCd['rho1']); PSI_MAX = float.fromhex(LCd['psi_max'])
import far_region
report('L3', far_region.RHO1 == RHO1, 'RHO1 = %r = %s (local.py and far_region.py)' % (RHO1, RHO1.hex()))
import local
lboxes = [[float.fromhex(x) for x in L] for L in LCd['leaves']]
inside, overlap, area_ok, n = tiles(lboxes, 0.0, RHO1, 0.0, PSI_MAX)
two_pi_ok = bool(arb(PSI_MAX) > 2 * arb.pi())
report('L1', inside and overlap == 0 and area_ok and two_pi_ok,
       '%d leaves tile [0,RHO1]x[0,%r]; overlaps=%d exact-area-sum=%s; PSI_MAX > 2pi: %s' % (n, PSI_MAX, overlap, area_ok, two_pi_ok))
recert = all(local.certify_box(*b) for b in lboxes)
report('L2', recert, '%d leaves re-certified (positive direction of 2Y - D_A/2, D_A > 0)' % n)
json.dump(dict(ok=bool(ok_all), leaves=n, overlaps=overlap, exact_area_sum=bool(area_ok), psi_max_gt_2pi=two_pi_ok,
               rho1=RHO1.hex(), psi_max=PSI_MAX.hex()), open('../verification/coverage_check.json', 'w'), indent=1)
print('COVERAGE CHECK:', 'ALL PASS' if ok_all else 'FAIL')
sys.exit(0 if ok_all else 1)
