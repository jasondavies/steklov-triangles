#!/usr/bin/env python3
"""Reproduce every certified inequality, every consistency check and every number in the paper.

  .venv/bin/python scripts/verify.py                 # everything
  .venv/bin/python scripts/verify.py --only certified|checks|numerics|paper|lean
  .venv/bin/python scripts/verify.py --skip-numerics # skip the float-only numerics

Each step writes a machine-readable JSON file to verification/.  Each gate below parses that JSON and
compares it with the paper's claims in claims.json, in EXACT rational arithmetic (certified endpoints
are stored as exact 'man*2^exp' numbers).  The certified steps read only committed certificate data
(data/), so they do not depend on BLAS or the CPU.  Logs go to verification/<step>.log and a summary
to verification/summary.json.  Exit status is nonzero if any gate fails.

Lean: by default lean/scripts/build_with_lean_path.sh uses the repository's own Lake environment
(lean/.lake/packages, e.g. after `cd lean && lake exe cache get`); set MATHLIB_PACKAGES to a
.lake/packages directory of an existing Mathlib v4.35.0 build to override.
"""
import argparse, hashlib, json, os, re, subprocess, sys, time
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / 'verification'
PY = sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
CLAIMS = json.loads((ROOT / 'claims.json').read_text())


def C(key):
    return Fr(CLAIMS[key])


def X(s):
    """Exact Fraction of an encoded 'man*2^exp' number."""
    m, e = s.split('*2^')
    m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** (-e))


def H(h):
    return Fr(float.fromhex(h))


def J(name):
    return json.loads((LOGS / name).read_text())


def sanitize(text):
    text = text.replace(str(ROOT), '.')
    pk = os.environ.get('MATHLIB_PACKAGES')
    if pk:
        text = text.replace(pk, '$MATHLIB_PACKAGES')
    return text.replace(str(Path.home()), '~')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(name, cmd, cwd, gate):
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=ENV)
    shown = ['python' if c == PY else str(c) for c in cmd]
    out = sanitize(proc.stdout + proc.stderr)
    msgs = []
    ok = proc.returncode == 0
    (LOGS / (name + '.log')).write_text('$ ' + ' '.join(shown) + '\n' + out)   # gates may read the log
    if ok:
        try:
            msgs = gate()
            ok = all(m[0] for m in msgs)
        except Exception as exc:          # a missing or malformed JSON output is a failure
            ok = False; msgs = [(False, 'gate error: %r' % exc)]
    (LOGS / (name + '.log')).write_text('$ ' + ' '.join(shown) + '\n' + out +
                                        ''.join('\n[gate] %s %s' % ('ok  ' if g else 'FAIL', t) for g, t in msgs) + '\n')
    print('%-22s %s  (%.0fs)' % (name, 'PASS' if ok else 'FAIL', time.time() - t0), flush=True)
    for g, t in msgs:
        if not g:
            print('    gate failed:', t)
    return {'name': name, 'ok': ok, 'seconds': round(time.time() - t0, 1), 'cmd': shown,
            'gates': [{'ok': g, 'check': t} for g, t in msgs]}


def near(x, key, tol):
    return abs(x - C(key)) <= Fr(tol)


# ------------------------------------------------------------------ gates (one per step)
def g_equilateral():
    d = J('equilateral.json')
    return [(X(d['Sigma_lo']) >= C('SigmaElo'), 'Sigma(E) lower endpoint >= %s' % CLAIMS['SigmaElo']),
            (X(d['Sigma_hi']) <= C('SigmaEhi'), 'Sigma(E) upper endpoint <= %s' % CLAIMS['SigmaEhi']),
            (X(d['area_lo']) >= C('areaElo'), 'sigma_1 sqrt|E| lower endpoint >= %s' % CLAIMS['areaElo']),
            (X(d['area_hi']) <= C('areaEhi'), 'sigma_1 sqrt|E| upper endpoint <= %s' % CLAIMS['areaEhi']),
            (X(d['eta_upper']) <= C('etaMax'), 'residual eta <= %s' % CLAIMS['etaMax']),
            (near(H(d['sigma_tilde_hex']), 'sigmaTilde', '5e-13'), 'sigma~ = %s to 12 digits' % CLAIMS['sigmaTilde'])]


def g_cr():
    d = J('cr_lower.json')
    return [(d['ok'], 'inertia certificate valid, O nonsingular, <= 2 boundary edges per cell'),
            (d['n_h'] == int(CLAIMS['crNh']) and d['dofs'] == int(CLAIMS['crDofs']), 'n_h = %s, %s dofs' % (CLAIMS['crNh'], CLAIMS['crDofs'])),
            (Fr(d['mu']) == C('crMu') and d['neg'] == int(CLAIMS['crNeg']), 'exactly %s negative eigenvalues of S_h - %s B_h' % (CLAIMS['crNeg'], CLAIMS['crMu'])),
            (d['pos'] == int(CLAIMS['crPos']), 'exactly %s positive eigenvalues' % CLAIMS['crPos']),
            (X(d['sigma4_lower']) > C('sigmaFourLo') > 0 and C('sigmaFourLo') ** 2 > 12,
             'sigma_4(E) > %s > 2 sqrt(3)' % CLAIMS['sigmaFourLo']),
            (near(X(d['Ch2_upper']), 'crChSq', '5e-7'), 'C_h^2 = %s...' % CLAIMS['crChSq'])]


def g_local():
    d = J('local_leaves.json')
    msgs = [(d['failed'] == 0, 'every local box passes'),
            (d['verified'] == int(CLAIMS['localBoxes']), '%s local boxes' % CLAIMS['localBoxes']),
            (H(d['rho1']) == Fr(float(CLAIMS['rhoOne'])), 'rho1 = double(%s)' % CLAIMS['rhoOne']),
            (H(d['psi_max']) == Fr(float(CLAIMS['psiMax'])), 'psi range [0, %s]' % CLAIMS['psiMax']),
            (X(d['eps_bdry_upper']) <= C('epsBdryMax'), 'eps_bdry <= %s' % CLAIMS['epsBdryMax']),
            (X(d['eps_grad_upper']) <= C('epsGradMax'), 'eps_grad <= %s' % CLAIMS['epsGradMax'])]
    glo, ghi = X(d['gram_lower']), X(d['gram_upper'])
    msgs.append((near(glo, 'gramE', '5e-7') and near(ghi, 'gramE', '5e-7'), 'Gram of u_1 = %s' % CLAIMS['gramE']))
    lo, hi = X(d['slope_lower']), X(d['slope_upper'])
    msgs.append((0 < C('coneLo') < lo <= hi < C('coneHi'),
                 'first-order splitting constant lies strictly between %s and %s' % (CLAIMS['coneLo'], CLAIMS['coneHi'])))
    return msgs


def g_far():
    d = J('far_region.json'); m = X(d['max_bound_upper'])
    return [(d['ok'], 'every far interval passes and the intervals tile [0,1/2]'),
            (d['intervals'] == int(CLAIMS['farIntervals']), '%s intervals' % CLAIMS['farIntervals']),
            (H(d['rho1']) == Fr(float(CLAIMS['rhoOne'])), 'same rho1 as the local region'),
            (m <= C('farBound'), 'far-region maximum <= %s' % CLAIMS['farBound']),
            (near(m, 'farMax', '5e-5'), 'far-region maximum = %s (4 decimals)' % CLAIMS['farMax']),
            (C('SigmaElo') - C('farBound') > Fr(1, 2), 'far-region deficit exceeds 1/2')]


def g_coverage():
    d = J('coverage_check.json')
    return [(d['ok'] and d['overlaps'] == 0 and d['exact_area_sum'] and d['psi_max_gt_2pi'], 'local boxes tile the disc'),
            (d['leaves'] == int(CLAIMS['localBoxes']), 'all %s local boxes re-certified' % CLAIMS['localBoxes']),
            (H(d['rho1']) == Fr(float(CLAIMS['rhoOne'])), 'same rho1')]


def g_ok(name):
    return lambda: [(J(name)['ok'], '%s reports ok' % name)]


def g_numerics_extra():
    d = J('numerics_extra.json')
    th = [t['theta'] for t in d['thin']]
    agree = max(abs(t['ratios'][-1] - t['theta']) / t['theta'] for t in d['thin'])
    return [(d['ok'], 'numerics_extra reports ok'),
            (near(Fr(d['cr_sigma_h1']['8']), 'crSigmaEight', '5e-6') and near(Fr(d['cr_sigma_h1']['16']), 'crSigmaSixteen', '5e-6'),
             'CR discrete sigma_{h,1} = %s, %s' % (CLAIMS['crSigmaEight'], CLAIMS['crSigmaSixteen'])),
            (near(Fr(d['theta0']), 'thetaZero', '5e-6'), 'theta(0) = %s' % CLAIMS['thetaZero']),
            (all(near(Fr(t), k, '5e-6') for t, k in zip(th, ['thetaA', 'thetaB', 'thetaC'])), 'theta(0.05), theta(0.25), theta(1/2)'),
            (Fr(agree) <= C('thinAgree'), 'thin-triangle agreement <= %s' % CLAIMS['thinAgree'])]


def g_table():
    d = {r['name']: r for r in J('table_bounds.json')}
    msgs = []
    for row in CLAIMS['table']:
        r = d[row['name']]
        up = X(r['certified_upper']); est = Fr(r['estimate']); err = Fr(r['err_total'])
        pest = Fr(row['estimate']); perr = Fr(row['error'])
        msgs.append((up <= Fr(row['upper']), '%s: certified upper bound <= printed %s' % (row['name'], row['upper'])))
        msgs.append((est - max(err, Fr(1, 10 ** 18)) <= up, '%s: estimate (minus its error) <= certified upper bound' % row['name']))
        msgs.append((r['retained_rank'] == r['basis_size'] and Fr(r['min_gram_eig']) > Fr(10) ** (20 - r['dps']),
                     '%s: all %d directions retained; variance Gram well conditioned at %d digits' % (row['name'], r['basis_size'], r['dps'])))
        msgs.append((abs(est - pest) + 2 * err <= perr,
                     '%s: printed %s +- %s covers rounding and twice the largest error indicator' % (row['name'], row['estimate'], row['error'])))
    return msgs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', choices=['certified', 'checks', 'numerics', 'paper', 'lean'])
    ap.add_argument('--skip-numerics', action='store_true')
    a = ap.parse_args()
    LOGS.mkdir(exist_ok=True)
    code = ROOT / 'code'
    want = lambda g: a.only in (None, g)
    steps = []
    if want('certified'):
        steps += [
            ('equilateral', [PY, 'equilateral.py'], code, g_equilateral),
            ('cr_lower', [PY, 'cr_lower.py'], code, g_cr),
            ('local', [PY, 'local.py', '0.6', '64'], code, g_local),
            ('far_region', [PY, 'far_region.py'], code, g_far),
            ('coverage_check', [PY, 'coverage_check.py'], code, g_coverage),
        ]
    if want('checks'):
        steps += [
            ('test_moments', [PY, 'test_moments.py'], code, g_ok('test_moments.json')),
            ('test_independent', [PY, 'test_independent.py'], code, g_ok('test_independent.json')),
            ('check_local_vs_ritz', [PY, 'check_local_vs_ritz.py'], code, g_ok('check_local_vs_ritz.json')),
            ('table_bounds', [PY, 'table_bounds.py'], code, g_table),
        ]
    if want('numerics') and not a.skip_numerics:
        steps += [
            ('numerics_extra', [PY, 'numerics_extra.py'], code, g_numerics_extra),
            ('map_sigma', [PY, 'map_sigma.py'], code, g_ok('map_sigma_summary.json')),
        ]
    if want('paper'):
        def g_claims_tex():
            tex = (ROOT / 'paper' / 'steklov-triangles.tex').read_text()
            abstract = tex[tex.index(r'\begin{abstract}'):tex.index(r'\end{abstract}')]
            return [(True, 'paper/claims.tex is in sync with claims.json'),
                    (CLAIMS['SigmaElo'] in abstract and CLAIMS['SigmaEhi'] in abstract,
                     'the literal enclosure in the abstract matches claims.json')]
        steps += [('claims_tex', [PY, 'scripts/make_claims_tex.py', '--check'], ROOT, g_claims_tex),
                  ('paper', ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-gg',
                             'steklov-triangles.tex'], ROOT / 'paper', lambda: [(True, 'latexmk succeeded')])]
    if want('lean'):
        def g_lean():
            o = (LOGS / 'lean.log').read_text()
            n_ax = o.count('depends on axioms')
            return [('error' not in o and 'sorryAx' not in o, 'no errors, no sorry'),
                    (n_ax > 0 and n_ax == o.count('depends on axioms: [propext, Classical.choice, Quot.sound]'),
                     'only the standard axioms (%d theorems)' % n_ax)]
        steps += [('lean', ['bash', 'scripts/build_with_lean_path.sh'], ROOT / 'lean', g_lean)]
    results = [run(*s) for s in steps]
    if want('paper'):
        log = ROOT / 'paper' / 'steklov-triangles.log'
        if log.exists():
            txt = log.read_text(errors='replace')
            bad = re.findall(r'(LaTeX Warning|Package \w+ Warning|Overfull \\[hv]box|Underfull \\[hv]box|There were undefined)', txt)
            results.append({'name': 'paper_log_clean', 'ok': not bad, 'issues': len(bad)})
            print('%-22s %s' % ('paper_log_clean', 'PASS' if not bad else 'FAIL (%d issues)' % len(bad)))
    inputs = sorted(p for p in ROOT.rglob('*') if p.is_file() and p.suffix in {'.py', '.lean', '.tex', '.md', '.toml', '.sh', '.txt', '.cff', '.json', '.npz'}
                    and not ({'.git', '.venv', 'build', 'verification', '.lake'} & set(p.parts)))
    import flint, numpy, scipy, mpmath
    summary = {'time_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'only': a.only,
               'skip_numerics': a.skip_numerics, 'python': sys.version.split()[0],
               'python_flint': flint.__version__, 'numpy': numpy.__version__, 'scipy': scipy.__version__,
               'mpmath': mpmath.__version__, 'openblas_coretype': os.environ.get('OPENBLAS_CORETYPE'),
               'mathlib_override': bool(os.environ.get('MATHLIB_PACKAGES')), 'results': results,
               'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in inputs}}
    pdf = ROOT / 'paper' / 'steklov-triangles.pdf'
    if pdf.exists():
        summary['pdf_sha256'] = sha(pdf)
    allok = all(r['ok'] for r in results)
    summary['overall'] = 'PASS' if allok else 'FAIL'
    full = a.only is None and not a.skip_numerics
    (LOGS / ('summary.json' if full else 'summary-partial.json')).write_text(json.dumps(summary, indent=1))
    print('OVERALL', summary['overall'])
    return 0 if allok else 1


if __name__ == '__main__':
    raise SystemExit(main())
