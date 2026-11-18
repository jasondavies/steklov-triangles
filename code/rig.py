"""
Rigorous (ball arithmetic, python-flint/Arb) building blocks.

A trial function is a real harmonic polynomial u = Re f, f(z) = sum_k c_k w^k,
w = (z - z0)/R, with complex coefficients c_k (any floats; they are converted
exactly to balls, so the trial function is a well-defined explicit function).

For a triangle with vertices Z[0..2] (acb balls, ccw) we compute, exactly up to
ball rounding:
  per edge e (from Z[j] to Z[j+1]): length l_e, G_e(a,b) = int_e u_a u_b ds,
  m_e(a) = int_e u_a ds;
  Q[a][b][i][j] = int_T d_i u_a d_j u_b dx  (i,j in {x,y}).
Area integrals are converted into boundary integrals by Green's formula
  int_T F_zbar dA = (1/2i) oint F dz,
with F = conj(f_b) f_a' (giving int f_a' conj(f_b') dA) and F = zbar f_a' f_b'
(giving int f_a' f_b' dA).  Along an edge z = a + t (b-a), 0<=t<=1, every
integrand is a polynomial in t, integrated exactly.
"""
from flint import acb, arb, acb_poly, ctx, fmpz
from fractions import Fraction

def enc(x):
    """Exact encoding 'man*2^exp' of an exact (radius-0) arb, e.g. x.lower() or x.upper()."""
    assert x.rad() == 0
    m, e = x.man_exp()
    return "%s*2^%s" % (m, e)

def dec(s):
    """Inverse of enc: an exact arb."""
    m, e = s.split("*2^")
    return arb(fmpz(int(m))) * arb(2) ** int(e)

def dec_fraction(s):
    """Exact Fraction of an enc() string."""
    m, e = s.split("*2^")
    m, e = int(m), int(e)
    return Fraction(m * 2 ** e) if e >= 0 else Fraction(m, 2 ** (-e))

def to_acb(c):
    if isinstance(c, acb):
        return c
    c = complex(c)
    return acb(arb(c.real), arb(c.imag))

def conj_poly(P):
    return acb_poly([c.conjugate() for c in P.coeffs()])

def re_poly(P):
    return acb_poly([acb(c.real) for c in P.coeffs()])

def int01(P):
    """integral_0^1 P(t) dt for an acb_poly."""
    s = acb(0)
    for k, c in enumerate(P.coeffs()):
        s += c / (k + 1)
    return s

def f_on_edge(coef, z0, R, a, b, deriv=False):
    """Return acb_poly in t of f(z(t)) (or f'(z(t))), z = a + t(b-a), f = sum_k coef[k-1] w^k."""
    w_lin = acb_poly([(a - z0) / R, (b - a) / R])
    n = len(coef)
    if not deriv:
        P = acb_poly([0])
        for k in range(n, 0, -1):          # Horner: f = w*(c1 + w*(c2 + ...))
            P = (P + acb_poly([coef[k - 1]])) * w_lin
        return P
    else:
        P = acb_poly([0])
        for k in range(n, 0, -1):          # f' = (1/R) sum k c_k w^{k-1}
            P = P * w_lin + acb_poly([k * coef[k - 1] / R])
        return P

def moments(funcs, Z):
    """funcs: list of (coef(list of acb), z0(acb), R(arb)). Z: 3 acb vertices (ccw).
    Returns dict with l[e], G[e][a][b], m[e][a], Q[a][b] (2x2 list [i][j])."""
    nf = len(funcs)
    l = []; G = []; m = []
    I1 = [[acb(0)] * nf for _ in range(nf)]   # int f_a' conj(f_b') dA
    I2 = [[acb(0)] * nf for _ in range(nf)]   # int f_a' f_b' dA
    for j in range(len(Z)):
        a = Z[j]; b = Z[(j + 1) % len(Z)]
        d = b - a
        le = abs(d)
        Fs = [f_on_edge(c, z0, R, a, b) for (c, z0, R) in funcs]
        Ds = [f_on_edge(c, z0, R, a, b, deriv=True) for (c, z0, R) in funcs]
        Us = [re_poly(F) for F in Fs]
        Ge = [[le * int01(Us[x] * Us[y]).real for y in range(nf)] for x in range(nf)]
        me = [le * int01(Us[x]).real for x in range(nf)]
        zbar = acb_poly([a.conjugate(), d.conjugate()])
        for x in range(nf):
            for y in range(nf):
                I1[x][y] += int01(conj_poly(Fs[y]) * Ds[x]) * d
                I2[x][y] += int01(zbar * Ds[x] * Ds[y]) * d
        l.append(le); G.append(Ge); m.append(me)
    half_over_i = acb(0, -0.5)   # 1/(2i)
    Q = [[None] * nf for _ in range(nf)]
    for x in range(nf):
        for y in range(nf):
            J1 = I1[x][y] * half_over_i
            J2 = I2[x][y] * half_over_i
            Qxx = (J1.real + J2.real) / 2
            Qyy = (J1.real - J2.real) / 2
            Qxy = (J1.imag - J2.imag) / 2      # int d_x u_a d_y u_b
            Qyx = (-J1.imag - J2.imag) / 2     # int d_y u_a d_x u_b
            Q[x][y] = [[Qxx, Qxy], [Qyx, Qyy]]
    return dict(l=l, G=G, m=m, Q=Q)

def residual_sq(coef, z0, R, Z, sig):
    """|| d_n u - sig u ||^2_{L2(dT)} and ||u||^2_{L2(dT)} for u = Re f (exact polynomials)."""
    r2 = arb(0); u2 = arb(0)
    for j in range(len(Z)):
        a = Z[j]; b = Z[(j + 1) % len(Z)]
        d = b - a; le = abs(d)
        nrm = acb(0, -1) * d / le          # outward normal (ccw polygon)
        F = f_on_edge(coef, z0, R, a, b)
        D = f_on_edge(coef, z0, R, a, b, deriv=True)
        U = re_poly(F)
        DN = re_poly(D * acb_poly([nrm]))    # d_n Re f = Re(f' n)
        res = DN - U * acb_poly([acb(sig)])
        r2 += le * int01(res * res).real
        u2 += le * int01(U * U).real
    return r2, u2

def residual_gram(funcs, Z, sig):
    """Gram matrices <r_a, r_b> and <u_a, u_b> in L2(dT), r_a = d_n u_a - sig u_a (exact polynomials)."""
    nf = len(funcs)
    Rg = [[arb(0)] * nf for _ in range(nf)]; Ug = [[arb(0)] * nf for _ in range(nf)]
    for j in range(len(Z)):
        a = Z[j]; b = Z[(j + 1) % len(Z)]
        d = b - a; le = abs(d)
        nrm = acb(0, -1) * d / le
        Us = []; Rs = []
        for (coef, z0, R) in funcs:
            U = re_poly(f_on_edge(coef, z0, R, a, b))
            DN = re_poly(f_on_edge(coef, z0, R, a, b, deriv=True) * acb_poly([nrm]))
            Us.append(U); Rs.append(DN - U * acb_poly([acb(sig)]))
        for x in range(nf):
            for y in range(nf):
                Rg[x][y] += le * int01(Rs[x] * Rs[y]).real
                Ug[x][y] += le * int01(Us[x] * Us[y]).real
    return Rg, Ug
