# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""
Anten-O-Matic: classic antennas.  Wire types are generated from the usual design formulas,
simulated in NEC-2 and optionally tuned in NEC; patch and dish are calculators.

Coordinates in metres.  Directional antennas (Yagi, Moxon) beam along +x with horizontal
elements along y; vertical / axial antennas (ground plane, QFH, helix) have their axis on +z.
"""
import json, math, os, time
import numpy as np

import core
import pattern as PAT
from core import C, MATERIALS, fillet, arc_points, bend_plan, develop, swr, _chk

# ------------------------------------------------------------------ type catalogue

def P(key, label, unit, default, lo=None, hi=None, tune=False, step='any', help=''):
    return dict(key=key, label=label, unit=unit, default=default, lo=lo, hi=hi, tune=tune, step=step, help=help)


def _lam(f0):
    return C / f0 * 1000  # mm


def moxon_dims(f0, wire_d_mm):
    """Cebik / AC6LA regression for the Moxon rectangle (dimensions in wavelengths)."""
    lam = _lam(f0)
    dw = min(max(wire_d_mm / lam, 1e-5), 1e-2)
    x = 0.4342945 * math.log(dw)
    A = -0.0008571428571 * x * x - 0.009571428571 * x + 0.3398571429
    B = -0.002142857143 * x * x - 0.02035714286 * x + 0.008285714286
    Cc = 0.001809523381 * x * x + 0.01780952381 * x + 0.05164285714
    D = 0.001 * x + 0.07178571429
    return A * lam, B * lam, Cc * lam, D * lam


def qfh_loop(L, ratio):
    """loop length L (mm) and width/height ratio -> (diameter, height) of a half-turn QFH loop."""
    h = L / (2 * math.sqrt(1 + (math.pi * ratio / 2) ** 2) + 2 * ratio)
    return ratio * h, h


def type_params(t, f0, wire_d_mm):
    lam = _lam(f0)
    if t == 'dipole':
        return [P('arm', 'Arm length', 'mm', round(0.24 * lam, 1), 0.15 * lam, 0.3 * lam, True),
                P('apex', 'Apex angle of the V (180 = straight dipole)', '°', 120, 60, 180, help='Inverted V: arms drooping downwards.')]
    if t == 'gp':
        return [P('vert', 'Radiator (vertical)', 'mm', round(0.24 * lam, 1), 0.15 * lam, 0.3 * lam, True),
                P('rad', 'Radials', 'mm', round(0.25 * lam, 1), 0.15 * lam, 0.35 * lam, True),
                P('nrad', 'Number of radials', '', 4, 2, 8, step=1),
                P('slope', 'Radial droop', '°', 45, 0, 60, help='0° = horizontal (~36 Ω), 45° ≈ 50 Ω.')]
    if t == 'yagi':
        n = 5
        p = [P('n', 'Number of elements', '', n, 3, 12, step=1),
             P('folded', 'Folded dipole (1 = yes)', '', 0, 0, 1, step=1, help='A folded dipole has ~4x the impedance (4:1 balun for 200 Ω).')]
        return p
    if t == 'moxon':
        A, B, Cc, D = moxon_dims(f0, wire_d_mm)
        return [P('A', 'A: width', 'mm', round(A, 1), 0.8 * A, 1.2 * A, True),
                P('B', 'B: driven element bent tail', 'mm', round(B, 1), 0.5 * B, 1.6 * B, True),
                P('C', 'C: gap between the tails', 'mm', round(Cc, 1), 0.6 * Cc, 1.5 * Cc, True),
                P('D', 'D: reflector bent tail', 'mm', round(D, 1), 0.7 * D, 1.4 * D, True)]
    if t == 'qfh':
        return [P('Lbig', 'Large loop length', 'mm', round(1.026 * lam, 1), 0.95 * lam, 1.12 * lam, True),
                P('Lsmall', 'Small loop length', 'mm', round(0.974 * lam, 1), 0.88 * lam, 1.05 * lam, True),
                P('ratio', 'Width/height ratio', '', 0.44, 0.3, 0.8, step=0.01),
                P('cone', 'Coverage ± around the zenith', '°', 80, 10, 90)]
    if t == 'helix':
        return [P('turns', 'Number of turns', '', 6, 2, 20, step=1),
                P('circ', 'Turn circumference', 'mm', round(1.0 * lam, 1), 0.75 * lam, 1.33 * lam, True),
                P('pitch', 'Pitch (turn spacing)', 'mm', round(0.24 * lam, 1), 0.12 * lam, 0.4 * lam, True),
                P('cone', 'Coverage ± around the axis', '°', 15, 0, 60)]
    if t == 'patch':
        return [P('er', 'Substrate permittivity εr', '', 1.0, 1, 12, step=0.01, help='1.0 = air (patch on standoffs), FR4 ≈ 4.3'),
                P('h', 'Substrate thickness', 'mm', round(max(1.6, 0.03 * lam), 2), 0.1, 50),
                P('t', 'Copper thickness', 'mm', 0.035, 0.01, 2)]
    if t == 'dish':
        return [P('D', 'Dish diameter', 'mm', round(10 * lam), 1 * lam, 200 * lam),
                P('fd', 'f/D ratio', '', 0.38, 0.2, 1.0, step=0.01),
                P('eta', 'Efficiency', '', 0.55, 0.2, 0.8, step=0.01),
                P('gores', 'Number of gores (template)', '', 12, 4, 36, step=1)]
    raise ValueError('unknown type %s' % t)


TYPES = {
    'dipole': ('Dipole / Inverted V', 'wire', 'Half-wave dipole, optionally with arms drooping into a V. Tunes arm length for SWR.'),
    'gp': ('Ground plane', 'wire', 'Quarter-wave vertical with drooping radials. Tunes the radiator and radials.'),
    'yagi': ('Yagi-Uda', 'wire', 'Reflector, driven element and directors on a boom. Tunes lengths and spacings for gain, front/back ratio and SWR.'),
    'moxon': ('Moxon', 'wire', 'Two-element antenna with bent tails (Cebik formulas), ~50 Ω, excellent front/back ratio.'),
    'qfh': ('QFH (quadrifilar helix)', 'wire', 'Two self-phased loops, RHCP for satellites (Meteor, NOAA). Tunes the loop lengths.'),
    'helix': ('Helix (axial mode)', 'wire', 'Helix over a reflector, circular polarisation, directional (Kraus).'),
    'patch': ('Patch (microstrip)', 'calc', 'Rectangular patch with an inset feed for 50 Ω. Computed from formulas, STEP of the board.'),
    'dish': ('Dish', 'calc', 'Dimensions, gain, feed (cantenna) and a gore template for building it.'),
}


def catalogue(f0, wire_d_mm):
    return [{'type': t, 'label': v[0], 'kind': v[1], 'desc': v[2], 'params': type_params(t, f0, wire_d_mm)}
            for t, v in TYPES.items()]

# ------------------------------------------------------------------ geometry

def yagi_init(n, lam):
    """reasonable starting point (DL6WU-like taper), lengths and positions in mm."""
    L = [0.495 * lam, 0.473 * lam] + [(0.440 - 0.008 * min(i, 6)) * lam for i in range(n - 2)]
    x = [0.0, 0.2 * lam]
    sp = [0.075, 0.18, 0.215, 0.25, 0.28, 0.30, 0.315, 0.33, 0.34, 0.35]
    for i in range(n - 2):
        x.append(x[-1] + sp[min(i, len(sp) - 1)] * lam)
    return L, x


def tuned_keys(t, p):
    """the parameter vector the NEC tuner works on: [(key, lo, hi)]"""
    if t == 'yagi':
        n, lam = int(p['n']), p['_lam']
        keys = [('L%d' % i, 0.85 * p['L%d' % i], 1.12 * p['L%d' % i]) for i in range(n)]
        keys += [('X%d' % i, p['X%d' % i] - 0.08 * lam, p['X%d' % i] + 0.08 * lam) for i in range(1, n)]
        return keys
    return [(q['key'], q['lo'], q['hi']) for q in type_params(t, p['f0'], p['wire_d'] * 1000) if q['tune']]


def build(t, p):
    """-> struct dict (metres): conductors (polylines to build), feeds, ground, notes."""
    mm = 1e-3
    wd = p['wire_d']
    lam = C / p['f0']
    g = max(0.006 * lam, 3 * wd)  # feed gap
    cond, feeds, ground, tl, sources = [], [], False, [], [(0, 0.0)]
    solids = []
    if t == 'dipole':
        a = p['arm'] * mm
        beta = math.radians((180 - p['apex']) / 2)
        e = np.array([math.cos(beta), 0, -math.sin(beta)])
        A = np.array([g / 2, 0, 0]); B = np.array([-g / 2, 0, 0])
        cond = [('Right arm', [A, A + e * a]), ('Left arm', [B, B + np.array([-e[0], 0, e[2]]) * a])]
        feeds = [(B, A)]
    elif t == 'gp':
        v, r, n, sl = p['vert'] * mm, p['rad'] * mm, int(p['nrad']), math.radians(p['slope'])
        hub, top = np.zeros(3), np.array([0, 0, g])
        cond = [('Radiator', [top, top + np.array([0, 0, v])])]
        for k in range(n):
            ph = 2 * math.pi * k / n
            cond.append(('Radial %d' % (k + 1), [hub, hub + r * np.array([math.cos(sl) * math.cos(ph), math.cos(sl) * math.sin(ph), -math.sin(sl)])]))
        feeds = [(hub, top)]
    elif t == 'yagi':
        n = int(p['n'])
        S = max(0.02 * lam, 4 * wd + 2 * p['bend_r'] * mm * 2)
        for i in range(n):
            L, x = p['L%d' % i] * mm, (p['X%d' % i] * mm if i else 0.0)
            name = 'Reflector' if i == 0 else 'Driven' if i == 1 else 'Director %d' % (i - 1)
            if i != 1:
                cond.append((name, [np.array([x, -L / 2, 0]), np.array([x, L / 2, 0])]))
            elif int(p['folded']):
                cond.append((name + ' (folded)', [np.array([x, -g / 2, 0]), np.array([x, -L / 2, 0]), np.array([x, -L / 2, S]),
                                                     np.array([x, L / 2, S]), np.array([x, L / 2, 0]), np.array([x, g / 2, 0])]))
                feeds = [(np.array([x, -g / 2, 0]), np.array([x, g / 2, 0]))]
            else:
                cond.append((name + ' left half', [np.array([x, -g / 2, 0]), np.array([x, -L / 2, 0])]))
                cond.append((name + ' right half', [np.array([x, g / 2, 0]), np.array([x, L / 2, 0])]))
                feeds = [(np.array([x, -g / 2, 0]), np.array([x, g / 2, 0]))]
    elif t == 'moxon':
        A, B, Cc, D = (p[k] * mm for k in 'ABCD')
        E = B + Cc + D
        cond = [('Driven left half', [np.array([E, -g / 2, 0]), np.array([E, -A / 2, 0]), np.array([E - B, -A / 2, 0])]),
                ('Driven right half', [np.array([E, g / 2, 0]), np.array([E, A / 2, 0]), np.array([E - B, A / 2, 0])]),
                ('Reflector', [np.array([D, -A / 2, 0]), np.array([0, -A / 2, 0]), np.array([0, A / 2, 0]), np.array([D, A / 2, 0])])]
        feeds = [(np.array([E, -g / 2, 0]), np.array([E, g / 2, 0]))]
    elif t == 'qfh':
        tw = p.get('_twist', 1)
        loops = []
        for j, (key, name) in enumerate((('Lbig', 'Large loop'), ('Lsmall', 'Small loop'))):
            dia, h = qfh_loop(p[key] * mm, p['ratio'])
            r = dia / 2
            ztop = g * 2 * j  # small loop's feed slightly above the big one's
            base = j * math.pi / 2
            def P3(rad, ang, z):
                return np.array([rad * math.cos(ang + base), rad * math.sin(ang + base), z])
            nh = 16
            # terminal A (+x side of the gap) -> radial out -> helix down -> bottom radial across -> helix up -> radial in -> terminal B
            pts = [P3(g / 2, 0, ztop), P3(r, 0, ztop)]
            pts += [P3(r, tw * math.pi * k / nh, ztop - h * k / nh) for k in range(1, nh + 1)]
            pts += [P3(r, tw * math.pi + math.pi, ztop - h)]  # across the bottom through the axis
            # second element climbs back on the opposite side: always 180 deg away from the first one
            pts += [P3(r, math.pi + tw * math.pi * (1 - k / nh), ztop - h + h * k / nh) for k in range(1, nh + 1)]
            end = pts[-1]
            pts += [P3(g / 2, math.pi, ztop)]
            loops.append((name, pts, (P3(g / 2, math.pi, ztop), P3(g / 2, 0, ztop))))
        cond = [(n_, pts) for n_, pts, _ in loops]
        feeds = [loops[0][2], loops[1][2]]
        tl = [(0, 1)]
    elif t == 'helix':
        n, circ, pitch = int(p['turns']), p['circ'] * mm, p['pitch'] * mm
        R = circ / (2 * math.pi)
        tw = p.get('_twist', 1)
        h0 = max(0.02 * lam, 3 * wd)
        nk = 24
        pts = [np.array([R, 0, h0])]
        pts += [np.array([R * math.cos(tw * 2 * math.pi * k / nk), R * math.sin(tw * 2 * math.pi * k / nk), h0 + pitch * k / nk])
                for k in range(1, n * nk + 1)]
        cond = [('Helix', pts)]
        feeds = [(np.array([R, 0, 0]), np.array([R, 0, h0]))]
        ground = True
        solids.append(('reflektor', 0.4 * lam + R))
    return {'cond': cond, 'feeds': feeds, 'ground': ground, 'tl': tl, 'solids': solids, 'gap': g}


def wires(st, p):
    """filleted wire list [(cond_index | 'feed', P, Q)] and the bend radius actually used."""
    rc = p['bend_r'] * 1e-3 + p['wire_d'] / 2
    out = [('feed', a, b) for a, b in st['feeds']]
    for ci, (name, pts) in enumerate(st['cond']):
        r, pieces = rc, None
        while pieces is None and r > 1e-6:
            pieces = fillet(pts, r, max(20.0, min(179.0, math.degrees(2.5 * p['wire_d'] / r))))  # short arcs -> sharp corners for NEC
            if pieces is None: r /= 2
        if pieces is None: pieces = fillet(pts, 0) or [('line', pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
        for pc in pieces:
            if pc[0] == 'line':
                out.append((ci, pc[1], pc[2]))
            else:
                _, t1, t2, c, ax, th = pc
                m = max(1, min(math.ceil(math.degrees(th) / 45), int(r * th / (2 * p['wire_d']))))
                ap = arc_points(t1, c, ax, th, m)
                out += [(ci, ap[j], ap[j + 1]) for j in range(m)]
    return out

# ------------------------------------------------------------------ NEC

def nec_run(st, p, freqs, rp, pat_step=None):
    import necwrap as N
    n = N.nec_create()
    try:
        feeds = []
        dd = {'wire_d': p['wire_d'], 'segs': 25}
        for tag, (k, a, b) in enumerate(wires(st, p), start=1):
            s = core.nsegs(float(np.linalg.norm(b - a)), dd, max(freqs), k == 'feed')
            _chk(N.nec_wire(n, tag, s, *a, *b, p['wire_d'] / 2, 1.0, 1.0))
            if k == 'feed': feeds.append((tag, (s + 1) // 2 if not st['ground'] else 1))
        if st['ground']:
            _chk(N.nec_geometry_complete(n, 1)); _chk(N.nec_gn_card(n, 1, 0, 0, 0, 0, 0, 0, 0))
        else:
            _chk(N.nec_geometry_complete(n, 0))
        _chk(N.nec_ek_card(n, 0))
        sig = MATERIALS.get(p.get('material', 'copper'), 5.8e7)
        if sig: _chk(N.nec_ld_card(n, 5, 0, 0, 0, sig, 0, 0))
        nf = len(freqs)
        _chk(N.nec_fr_card(n, 0, nf, freqs[0], (freqs[-1] - freqs[0]) / (nf - 1) if nf > 1 else 0))
        for i, j in st['tl']:
            _chk(N.nec_tl_card(n, feeds[i][0], feeds[i][1], feeds[j][0], feeds[j][1], p['z0'], 1e-5, 0, 0, 0, 0))
        _chk(N.nec_excitation_voltage(n, feeds[0][0], feeds[0][1], 1.0, 0.0))
        nth, nph, th0, ph0, dth, dph = rp
        _chk(N.nec_rp_card(n, 0, nth, nph, 0, 0, 0, 0, th0, ph0, dth, dph, 0, 0))
        res = []
        for fi, f in enumerate(freqs):
            z = complex(N.nec_impedance_real_indexed(n, fi, 0), N.nec_impedance_imag_indexed(n, fi, 0))
            r = {'f': f, 'z': z}
            for pol in ('total', 'rhcp', 'lhcp'):
                r[pol] = {s: N.gain(n, pol, s, fi) for s in ('min', 'mean', 'max')}
            if pat_step:
                r['pattern'] = PAT.from_nec(N, n, fi, pat_step, st['ground'])
            elif nth * nph <= 5000:
                r['grid'] = np.array([[N.nec_gain(n, fi, i, j) for j in range(nph)] for i in range(nth)])
            res.append(r)
        return res
    finally:
        N.nec_delete(n)


GRID = (19, 36, 0, 0, 10, 10)  # theta 0..180, phi 0..350
_th, _ph = np.meshgrid(np.radians(np.arange(19) * 10), np.radians(np.arange(36) * 10), indexing='ij')
BACK = np.sin(_th) * np.cos(_ph) <= -0.5


def metrics(t, st, p, freqs):
    """evaluate the figures of merit the type cares about, worst case over freqs."""
    cp = t in ('qfh', 'helix')
    if cp:
        cone = p['cone']
        rp = (int(cone // 10) + 1 if cone >= 10 else 2, 24, 0, 0, max(1, min(10, cone)), 15)
    else:
        rp = GRID
    res = nec_run(st, p, freqs, rp)
    m = {'swr': max(swr(r['z'], p['z0']) for r in res), 'z': [r['z'] for r in res]}
    if cp:
        pol = p.get('pol', 'rhcp')
        m['floor'] = min(r[pol]['min'] for r in res)
        m['mean'] = min(r[pol]['mean'] for r in res)
        m['gmax'] = min(r[pol]['max'] for r in res)
    else:
        m['gmax'] = min(r['total']['max'] for r in res)
        if t in ('yagi', 'moxon'):
            fw, fb, fr = [], [], []
            for r in res:
                G = r['grid']  # theta x phi, step 10
                fwd = G[9, 0]
                fw.append(fwd); fb.append(fwd - G[9, 18]); fr.append(fwd - G[BACK].max())  # F/B straight back, F/R within 60 deg
            m['fwd'], m['fb'], m['fr'] = min(fw), min(fb), min(fr)
        if t == 'gp':
            m['low'] = min(r['grid'][6:10, :].mean() for r in res)  # theta 60..90 = elevation 0..30
    return m


def score(t, m, p):
    over = max(0.0, m['swr'] - p['swr'])
    pen = 3 * over + 0.3 * over ** 2
    if t in ('qfh', 'helix'):
        return (m['mean'] + 0.3 * m['floor'] if t == 'qfh' else m['floor'] + 0.2 * m['gmax']) - pen
    if t in ('yagi', 'moxon'):
        return m['fwd'] + 0.08 * min(m['fb'], 30) + 0.06 * min(m['fr'], 25) - pen
    if t == 'gp':
        return 0.3 * m['low'] - m['swr'] - pen
    return 0.1 * m['gmax'] - m['swr'] - pen


_T = None


def _init(t, p):
    global _T
    _T = (t, p)


def _evalvec(x):
    t, p = _T
    q = dict(p)
    for (k, lo, hi), v in zip(tuned_keys(t, p), x):
        q[k] = float(v)
    try:
        st = build(t, q)
        m = metrics(t, st, q, [q['f0'] - q['bw'] / 2, q['f0'], q['f0'] + q['bw'] / 2])
        return score(t, m, q), m
    except Exception:
        return -1e3, None


def prepare(t, p):
    """fill defaults and derived parameters (input lengths in mm)."""
    p = dict(p)
    p['f0'], p['bw'] = float(p['f0']), float(p.get('bw') or 0)
    p['wire_d'] = float(p['wire_d'])  # metres
    p['bend_r'] = float(p.get('bend_r') or p['wire_d'] * 1000)  # mm
    # a helix is ~140 ohm (matched with a strip transformer), a QFH is a low-Z beast: do not chase 50 ohm too hard
    p['z0'] = float(p.get('z0') or (140 if t == 'helix' else 50))
    p['swr'] = float(p.get('swr') or {'qfh': 3.0, 'helix': 2.0}.get(t, 1.5))
    p['_lam'] = C / p['f0'] * 1000
    for q in type_params(t, p['f0'], p['wire_d'] * 1000):
        if p.get(q['key']) in (None, ''): p[q['key']] = q['default']
        p[q['key']] = float(p[q['key']])
    if t == 'yagi':
        n = int(p['n'])
        L, x = yagi_init(n, p['_lam'])
        for i in range(n):
            p.setdefault('L%d' % i, L[i]); p.setdefault('X%d' % i, x[i])
    return p


class Tuner:
    def __init__(self, t, p, progress=None, stop=None, log=None):
        self.t, self.p = t, prepare(t, p)
        self.progress, self.stop = progress or (lambda s: None), stop
        self.log = log or (lambda m: None)

    def _stopped(self):
        return self.stop is not None and self.stop.is_set()

    def run(self):
        t, p = self.t, self.p
        t0 = time.time()
        freqs = [p['f0'] - p['bw'] / 2, p['f0'], p['f0'] + p['bw'] / 2]
        self.log('Start: %s, %.3f MHz ± %.3f, wire Ø %.2f mm, bend radius %.2f mm'
                 % (TYPES[t][0], p['f0'], p['bw'] / 2, p['wire_d'] * 1000, p['bend_r']))
        if t in ('qfh', 'helix'):
            # winding sense decides RHCP vs LHCP: let NEC tell which one is right
            best = None
            for tw in (1, -1):
                q = dict(p, _twist=tw)
                m = metrics(t, build(t, q), q, [p['f0']])
                self.log('Winding %s: %s on axis %.2f dBic' % ('left-handed' if tw > 0 else 'right-handed', p.get('pol', 'rhcp').upper(), m['gmax']))
                if best is None or m['gmax'] > best[0]: best = (m['gmax'], tw)
            p['_twist'] = best[1]
        st = build(t, p)
        m = metrics(t, st, p, freqs)
        f = score(t, m, p)
        self.log('Starting design from formulas: %s' % describe(t, m))
        self.progress({'phase': 'tune', 'msg': 'Starting design computed', 'wires': wires_json(st, p), 'info': m,
                       '_pat': (lambda pp: (lambda: full_pattern(build(t, pp), pp)))(dict(p)),
                       'fit': f, 'hist': [f], 'beat': time.time(), 'gen': 0, 'gens': 1})
        keys = tuned_keys(t, p)
        budget = int(p.get('iters') or 0)
        if budget and keys:
            from multiprocessing import Pool
            jobs = int(p.get('jobs') or 0) or max(1, (os.cpu_count() or 2) - 1)
            lo = np.array([k[1] for k in keys]); hi = np.array([k[2] for k in keys])
            x = np.array([p[k[0]] for k in keys], float)
            rng = np.random.default_rng(int(p.get('seed') or 0) or None)
            sigma, hist, evals = 0.05, [f], 0
            rounds = max(1, budget // (jobs * 2))
            self.log('Tuning %d parameters in NEC: %d rounds of %d candidates on %d cores' % (len(keys), rounds, jobs * 2, jobs))
            with Pool(jobs, initializer=_init, initargs=(t, p)) as pool:
                for it in range(rounds):
                    if self._stopped():
                        self.log('Stopped by the user.'); break
                    cand = [np.clip(x + rng.normal(0, sigma, x.size) * (hi - lo), lo, hi) for _ in range(jobs * 2)]
                    r = pool.map(_evalvec, cand)
                    evals += len(cand)
                    j = int(np.argmax([q[0] for q in r]))
                    if r[j][0] > f:
                        f, x, m = r[j][0], cand[j], r[j][1]
                        self.log('Round %d/%d: better (%.3f) %s' % (it + 1, rounds, f, describe(t, m)))
                    else:
                        sigma = max(0.003, sigma * 0.85)
                    hist.append(f)
                    for (k, _, _), v in zip(keys, x): p[k] = float(v)
                    self.progress({'phase': 'tune', 'gen': it + 1, 'gens': rounds, 'fit': f, 'info': m, 'hist': hist,
                                   'wires': wires_json(build(t, p), p), 'beat': time.time(), 'evals': evals,
                                   '_pat': (lambda pp: (lambda: full_pattern(build(t, pp), pp)))(dict(p)),
                                   'rate': evals / max(1e-6, time.time() - t0),
                                   'msg': 'NEC tuning %d/%d' % (it + 1, rounds)})
            for (k, _, _), v in zip(keys, x): p[k] = float(v)
        self.log('Done in %.0f s' % (time.time() - t0))
        return p


def describe(t, m):
    s = 'SWR %.2f' % m['swr']
    if 'fwd' in m: s += ', forward gain %.2f dBi, F/B %.1f dB, F/R %.1f dB' % (m['fwd'], m['fb'], m['fr'])
    elif 'floor' in m: s += ', min in cone %.2f dBic, max %.2f dBic' % (m['floor'], m['gmax'])
    else: s += ', max gain %.2f dBi' % m['gmax']
    if 'low' in m: s += ', mean below 30° elevation %.2f dBi' % m['low']
    return s

# ------------------------------------------------------------------ outputs

def wires_json(st, p):
    return [{'arm': '-1' if k == 'feed' else str(k), 'feed': k == 'feed', 'p': [float(v) for v in a], 'q': [float(v) for v in b]}
            for (k, a, b) in wires(st, p)]


def full_pattern(st, p, step=5):
    return nec_run(st, p, [p['f0']], PAT.rp_for(step, st['ground']), pat_step=step)[0]['pattern']


def cut_gains(t, st, p, cut, pat=None):
    """gain along a great-circle cut: 'az' = horizontal plane, 'el0' / 'el90' = vertical plane at phi 0 / 90."""
    pol = 'total' if t not in ('qfh', 'helix') else p.get('pol', 'rhcp')
    pat = pat or full_pattern(st, p)
    if cut == 'az':
        pts = [(90, ph) for ph in range(0, 361, 5)]
    else:
        ph = 0 if cut == 'el0' else 90
        pts = [(abs(a), ph if a >= 0 else ph + 180) for a in range(-180, 181, 5)]
    return pts, [core.pat_sample(pat, pol, th, ph) for th, ph in pts]


def render(t, st, p, path, title, pat=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 5.2))
    ax = fig.add_subplot(131, projection='3d')
    cols = ['#2a78d6', '#e8710a', '#1aa36f', '#d63a4a', '#8a5cd6', '#c79a00']
    allp = []
    for k, a, b in wires(st, p):
        ax.plot(*zip(a * 1000, b * 1000), color='k' if k == 'feed' else cols[k % len(cols)], lw=2.5)
        allp += [a * 1000, b * 1000]
    allp = np.array(allp)
    c = (allp.max(0) + allp.min(0)) / 2; s = max((allp.max(0) - allp.min(0)).max() / 2, 1) * 1.1
    if st['ground']:
        gx, gy = np.meshgrid([c[0] - s, c[0] + s], [c[1] - s, c[1] + s])
        ax.plot_surface(gx, gy, 0 * gx, alpha=0.15, color='grey')
    ax.set_xlim(c[0] - s, c[0] + s); ax.set_ylim(c[1] - s, c[1] + s); ax.set_zlim(c[2] - s, c[2] + s)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]'); ax.set_zlabel('z [mm]')
    ax.set_title(title)
    unit = 'dBi' if t not in ('qfh', 'helix') else 'dBic ' + p.get('pol', 'rhcp').upper()
    cuts = [('az', 'horizontal plane'), ('el0', 'vertical plane xz')] if t in ('yagi', 'moxon', 'dipole') else \
           [('el0', 'vertical plane xz'), ('el90', 'vertical plane yz')]
    for i, (cut, lab) in enumerate(cuts):
        pts, g = cut_gains(t, st, p, cut, pat)
        g = np.array(g)
        top = np.ceil(np.nanmax(g) / 5) * 5
        floor = top - 30
        axp = fig.add_subplot(1, 3, 2 + i, projection='polar')
        if cut == 'az':
            ang = np.radians([q[1] for q in pts]); axp.set_theta_zero_location('E')
        else:
            ang = np.radians(np.arange(-180, 181, 5)); axp.set_theta_zero_location('N'); axp.set_theta_direction(-1)
        axp.plot(ang, np.clip(g, floor, None), color='#e8710a', lw=2)
        axp.set_rlim(floor, top)
        axp.set_rticks([top - 30, top - 20, top - 10, top])
        axp.set_title('%s [%s], max %.1f' % (lab, unit, np.nanmax(g)), fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=105); plt.close(fig)


def build_text(t, st, p):
    rc = p['bend_r'] + p['wire_d'] * 500
    lines = []
    if t == 'yagi':
        n = int(p['n'])
        lines.append('Elements on the boom (position from the reflector, tip-to-tip lengths, elements insulated from the boom):')
        for i in range(n):
            name = 'Reflector' if i == 0 else 'Driven' if i == 1 else 'Director %d' % (i - 1)
            lines.append('  %-12s length %7.1f mm   position %7.1f mm' % (name, p['L%d' % i], p['X%d' % i] if i else 0.0))
        lines.append('  Boom length %.0f mm (+ overhangs). The driven element is split in the middle with a %.1f mm feed gap.'
                     % (p['X%d' % (n - 1)], st['gap'] * 1000))
        if int(p['folded']): lines.append('  Folded driven element: two conductors one above the other, bends below, fed through a 4:1 balun.')
        lines.append('  Caution: a metal boom the elements pass through lengthens them electrically. The model assumes insulated elements.')
        lines.append('')
    if t == 'moxon':
        A, B, Cc, D = (p[k] for k in 'ABCD')
        lines += ['Rectangle dimensions: A = %.1f mm (width), B = %.1f mm, C = %.1f mm (gap), D = %.1f mm, E = %.1f mm (depth).'
                  % (A, B, Cc, D, B + Cc + D), 'Gap C is critical, hold it with insulating spacers.', '']
    if t == 'qfh':
        for key, name in (('Lbig', 'Large'), ('Lsmall', 'Small')):
            dia, h = qfh_loop(p[key], p['ratio'])
            lines.append('%s loop: length %.1f mm, diameter %.1f mm, height %.1f mm, half turn %s.'
                         % (name, p[key], dia, h, 'left-handed' if p.get('_twist', 1) > 0 else 'right-handed'))
        lines += ['The loops are rotated 90° and both fed in parallel at the top (one coax; the loops cross with a small height offset).',
                  'NEC picked the winding sense that gives the chosen polarisation upwards.', '']
    if t == 'helix':
        R = p['circ'] / (2 * math.pi)
        lines += ['Helix: %d turns, diameter %.1f mm (circumference %.1f mm), pitch %.1f mm, height %.0f mm, winding %s.'
                  % (int(p['turns']), 2 * R, p['circ'], p['pitch'], p['pitch'] * p['turns'], 'left-handed' if p.get('_twist', 1) > 0 else 'right-handed'),
                  'Wind it on a tube/former of %.1f mm diameter (inner turn diameter); the end goes through a connector in the reflector.' % (2 * R - p['wire_d'] * 1000),
                  'Reflector: a round plate at least %.0f mm across.' % (0.8 * p['_lam'] + 2 * R),
                  'A helix is usually 100-150 Ω: a matching strip along the first quarter turn brings it to 50 Ω.', '']
    if t == 'gp':
        lines += ['Radiator and radials meet at one point (connector centre = radiator, shell = radials).', '']
    if t == 'dipole':
        lines += ['The arms are separated by a %.1f mm feed gap. Apex angle %.0f°.' % (st['gap'] * 1000, p['apex']), '']
    lines.append('Conductors (cut length = developed length along the axis of Ø %.2f mm wire, inner bend radius %.1f mm):' % (p['wire_d'] * 1000, p['bend_r']))
    for name, pts in st['cond']:
        steps = bend_plan([q for q in pts])
        total = develop(steps, rc, False)
        nb = sum(1 for s_ in steps if s_['bend'] and s_['bend'] >= 0.5)
        lines.append('%s: cut %.1f mm%s' % (name, total, '' if nb else ' (straight)'))
        if 0 < nb <= 10:
            for i, s_ in enumerate(steps, 1):
                if s_['bend'] is None:
                    lines.append('  %d. straight %.1f mm, then the end' % (i, s_['straight']))
                elif s_['bend'] >= 0.5:
                    lines.append('  %d. straight %.1f mm, mark at %.1f mm = start of bend: bend %.0f°, rotate the plane %+.0f° %s'
                                 % (i, s_['straight'], s_['mark'], s_['bend'], s_['rot'], 'from the reference plane' if i == 1 else 'from the previous bend plane'))
        elif nb > 10:
            lines.append('  (%d small bends: wind it on a former using the dimensions above; the STEP model shows the exact shape)' % nb)
    return '\n'.join(lines)


def export_step(st, p, path, name):
    import stepout
    rc = p['bend_r'] + p['wire_d'] * 500
    bodies = []
    for ci, (cname, pts) in enumerate(st['cond']):
        P_ = [q * 1000 for q in pts]
        r, pieces = rc, None
        while pieces is None and r > 1e-3:
            pieces = fillet(P_, r)
            if pieces is None: r /= 2
        if pieces is None: continue
        safe = ''.join(c if c.isalnum() else '_' for c in cname)
        bodies.append((safe, pieces, p['wire_d'] * 500))
    for sname, R in st['solids']:
        bodies.append((sname, [('line', np.array([0, 0, -2.0]), np.array([0, 0, 0.0]))], R * 1000))
    stepout.write_step(path, bodies, name)


def save(t, p, out_root):
    st = build(t, p)
    freqs = list(np.linspace(p['f0'] - p['bw'] / 2, p['f0'] + p['bw'] / 2, 5)) if p['bw'] else [p['f0']]
    rows = []
    for f in freqs:
        m = metrics(t, st, p, [f])
        rows.append({'f': f, 'zr': m['z'][0].real, 'zi': m['z'][0].imag, 'swr': m['swr'],
                     'gmax': m['gmax'], 'fwd': m.get('fwd'), 'fb': m.get('fb'), 'floor': m.get('floor'), 'mean': m.get('mean')})
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in p.get('name') or t)[:40] or t
    out = os.path.join(out_root, '%s_%s_%gMHz_%s' % (stamp, safe, p['f0'], t))
    os.makedirs(out, exist_ok=True)
    title = '%s %.1f MHz' % (TYPES[t][0], p['f0'])
    pat = full_pattern(st, p)
    with open(os.path.join(out, 'pattern.json'), 'w') as f: json.dump(pat, f)
    render(t, st, p, os.path.join(out, 'render.png'), title, pat)
    if p.get('step', True): export_step(st, p, os.path.join(out, 'antenna.step'), safe)
    write_nec(st, p, os.path.join(out, 'antenna.nec'))
    btxt = build_text(t, st, p)
    unit = 'dBic ' + p.get('pol', 'rhcp').upper() if t in ('qfh', 'helix') else 'dBi'
    rep = ['Anten-O-Matic: %s' % title, '=' * 40, TYPES[t][2], '',
           'λ = %.1f mm, wire Ø %.2f mm, feed impedance %.0f Ω' % (p['_lam'], p['wire_d'] * 1000, p['z0']), '']
    hdr = 'MHz         Z [Ω]            SWR    max [%s]' % unit
    if t in ('yagi', 'moxon'): hdr += '   fwd     F/B [dB]'
    if t in ('qfh', 'helix'): hdr += '   min in cone    mean'
    rep.append(hdr)
    for r in rows:
        line = '%9.3f  %7.1f %+7.1fj   %5.2f   %6.2f' % (r['f'], r['zr'], r['zi'], r['swr'], r['gmax'])
        if r['fwd'] is not None: line += '   %6.2f   %6.1f' % (r['fwd'], r['fb'])
        if r['floor'] is not None: line += '   %8.2f      %6.2f' % (r['floor'], r['mean'])
        rep.append(line)
    rep += ['', 'CONSTRUCTION', '-' * 40, btxt]
    text = '\n'.join(rep) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as f: f.write(text)
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as f:
        json.dump({'type': t, 'params': {k: v for k, v in p.items()}, 'summary': rows}, f, indent=1, ensure_ascii=False, default=str)
    summ = [{'f': r['f'], 'floor': r['floor'] if r['floor'] is not None else (r['fwd'] if r['fwd'] is not None else r['gmax']),
             'mean': r['mean'] if r['mean'] is not None else (r['fb'] if r['fb'] is not None else None),
             'max': r['gmax'], 'zr': r['zr'], 'zi': r['zi'], 'swr': r['swr']} for r in rows]
    return out, text, summ, btxt, wires_json(st, p)


def write_nec(st, p, path):
    with open(path, 'w') as f:
        f.write('CM Anten-O-Matic %s %.3f MHz\nCE\n' % (p.get('_type', ''), p['f0']))
        feeds = []
        dd = {'wire_d': p['wire_d'], 'segs': 25}
        for tag, (k, a, b) in enumerate(wires(st, p), start=1):
            s = core.nsegs(float(np.linalg.norm(b - a)), dd, p['f0'] + p['bw'], k == 'feed')
            f.write('GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % (tag, s, *a, *b, p['wire_d'] / 2))
            if k == 'feed': feeds.append((tag, (s + 1) // 2 if not st['ground'] else 1))
        f.write('GE 1\nGN 1\n' if st['ground'] else 'GE 0\n')
        f.write('EK 0\n')
        for i, j in st['tl']:
            f.write('TL %d %d %d %d %g 0.00001\n' % (*feeds[i], *feeds[j], p['z0']))
        f.write('FR 0 21 0 0 %.4f %.4f\n' % (p['f0'] - p['bw'], max(p['bw'], 0.001) / 10))
        f.write('EX 0 %d %d 0 1 0\n' % feeds[0])
        f.write('RP 0 37 73 1000 0 0 5 5\nEN\n')

# ------------------------------------------------------------------ calculators

def patch_calc(p):
    c = 299792458.0
    f = p['f0'] * 1e6
    er, h = p['er'], p['h'] * 1e-3
    W = c / (2 * f) * math.sqrt(2 / (er + 1))
    eeff = (er + 1) / 2 + (er - 1) / 2 / math.sqrt(1 + 12 * h / W)
    dL = 0.412 * h * (eeff + 0.3) * (W / h + 0.264) / ((eeff - 0.258) * (W / h + 0.8))
    L = c / (2 * f * math.sqrt(eeff)) - 2 * dL
    lam0 = c / f
    # edge resistance (Balanis, slot conductance model)
    k0 = 2 * math.pi / lam0
    G1 = W / (120 * lam0) * (1 - (k0 * h) ** 2 / 24) if W < lam0 else W / (120 * lam0)
    Rin = 1 / (2 * G1)
    z0 = p['z0']
    y0 = L / math.pi * math.acos(math.sqrt(min(1.0, z0 / Rin)))
    # 50 ohm microstrip width (Hammerstad / Wheeler)
    A = z0 / 60 * math.sqrt((er + 1) / 2) + (er - 1) / (er + 1) * (0.23 + 0.11 / er)
    wh = 8 * math.exp(A) / (math.exp(2 * A) - 2)
    if wh > 2:
        B = 377 * math.pi / (2 * z0 * math.sqrt(er))
        wh = 2 / math.pi * (B - 1 - math.log(2 * B - 1) + (er - 1) / (2 * er) * (math.log(B - 1) + 0.39 - 0.61 / er))
    wf = wh * h
    gap = wf  # inset notch gap each side
    Lg, Wg = L + 6 * h + lam0 / 4, W + 6 * h + lam0 / 4
    bw = 3.77 * (er - 1) / er ** 2 * h / lam0 * W / L * 100 if er > 1 else 3.77 * h / lam0 * W / L * 100 * 0.5
    D = 10 * math.log10(max(1.0, 6.6 if er > 1.5 else 8.5))
    return dict(W=W * 1e3, L=L * 1e3, eeff=eeff, Rin=Rin, y0=y0 * 1e3, wf=wf * 1e3, gap=gap * 1e3, Lg=Lg * 1e3, Wg=Wg * 1e3,
                bw=bw, lam=lam0 * 1e3, gain=D)


def patch_save(p, out_root):
    r = patch_calc(p)
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in p.get('name') or 'patch')[:40] or 'patch'
    out = os.path.join(out_root, '%s_%s_%gMHz_patch' % (stamp, safe, p['f0']))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'pattern.json'), 'w') as f:
        json.dump(PAT.patch_pattern(r['L'], r['W'], r['eeff'], p['f0'], r['gain']), f)
    h, t = p['h'], p['t']
    # patch outline with inset feed notch and feed line (mm), x along L (feed side at x=0), y along W
    L, W, y0, wf, gp = r['L'], r['W'], r['y0'], r['wf'], r['gap']
    xl, yl = -(r['Lg'] - L) / 2, 0  # feed line runs from ground edge to the notch
    poly = [(0, -W / 2), (L, -W / 2), (L, W / 2), (0, W / 2), (0, wf / 2 + gp), (y0, wf / 2 + gp), (y0, wf / 2),
            (xl, wf / 2), (xl, -wf / 2), (y0, -wf / 2), (y0, -wf / 2 - gp), (0, -wf / 2 - gp)]
    ground = [(xl, -r['Wg'] / 2), (xl + r['Lg'], -r['Wg'] / 2), (xl + r['Lg'], r['Wg'] / 2), (xl, r['Wg'] / 2)]
    import stepout
    stepout.write_step_prisms(os.path.join(out, 'antenna.step'), [
        ('zem', ground, -t, 0.0), ('substrat', ground, 0.0, h), ('patch', poly, h, h + t)], safe)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 7))
    gx, gy = zip(*(ground + [ground[0]])); ax.fill(gx, gy, color='#cfd6dd'); ax.plot(gx, gy, color='#5d6670')
    px, py = zip(*(poly + [poly[0]])); ax.fill(px, py, color='#e8a33a'); ax.plot(px, py, color='#9a5d00')
    ax.annotate('', (0, -W / 2 - 3), (L, -W / 2 - 3), arrowprops=dict(arrowstyle='<->'))
    ax.text(L / 2, -W / 2 - 5, 'L = %.2f mm' % L, ha='center', va='top')
    ax.annotate('', (L + 3, -W / 2), (L + 3, W / 2), arrowprops=dict(arrowstyle='<->'))
    ax.text(L + 5, 0, 'W = %.2f mm' % W, rotation=90, va='center')
    ax.text(y0 / 2, wf / 2 + gp + 1, 'inset %.2f' % y0, ha='center', fontsize=8)
    ax.set_aspect('equal'); ax.set_title('Patch %.1f MHz, εr %.2f, h %.2f mm (top view, mm)' % (p['f0'], p['er'], h))
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'render.png'), dpi=110); plt.close(fig)
    text = '\n'.join([
        'Anten-O-Matic: Patch %.3f MHz' % p['f0'], '=' * 40,
        'Substrate εr = %.2f, thickness %.2f mm, copper %.3f mm; λ = %.1f mm, εeff = %.3f' % (p['er'], h, t, r['lam'], r['eeff']), '',
        'Patch:          L = %.2f mm (resonant dimension, along the feed)' % L,
        '                W = %.2f mm' % W,
        'Edge resistance: %.0f Ω -> feed inset %.2f mm for %.0f Ω' % (r['Rin'], y0, p['z0']),
        'Feed line:      width %.2f mm (%.0f Ω), gap around the inset %.2f mm' % (wf, p['z0'], gp),
        'Ground plane:   %.1f x %.1f mm (at least L,W + 6h + λ/4)' % (r['Lg'], r['Wg']),
        'Estimated bandwidth (SWR<2): ~%.1f %%, directivity ~%.0f dBi (FR4 really 1-3 dB less, losses)' % (r['bw'], r['gain']), '',
        'Notes: the formulas (transmission line model, Balanis) are typically 1-3 % off. Make it slightly longer and trim.',
        'For an air patch a microstrip inset is impractical: feed it with coax from below at the inset point (pin through the ground).',
        'Circular polarisation: trim two opposite corners (by about 5-10 % of L), feed on the axis, fine-tune by measurement.',
    ]) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as f: f.write(text)
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as f: json.dump({'type': 'patch', 'params': p, 'result': r}, f, indent=1, default=str)
    return out, text


def dish_calc(p):
    lam = _lam(p['f0'])
    D, fd = p['D'], p['fd']
    f = fd * D
    depth = D * D / (16 * f)
    gain = 10 * math.log10(p['eta'] * (math.pi * D / lam) ** 2)
    hpbw = 70 * lam / D
    psi = 4 * math.degrees(math.atan(1 / (4 * fd)))  # full angle subtended by the rim, seen from the focus
    # cantenna (circular waveguide, TE11): pick diameter so f0 sits 1.25x above TE11 cutoff, below TM01
    dwg = 1.25 * 1.841 * lam / math.pi
    lc = 1.706 * dwg
    lg = lam / math.sqrt(1 - (lam / lc) ** 2)
    return dict(lam=lam, f=f, depth=depth, gain=gain, hpbw=hpbw, psi=psi, dwg=dwg, lg=lg, probe=lam / 4)


def dish_save(p, out_root):
    r = dish_calc(p)
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in p.get('name') or 'dish')[:40] or 'dish'
    out = os.path.join(out_root, '%s_%s_%gMHz_dish' % (stamp, safe, p['f0']))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'pattern.json'), 'w') as fh:
        json.dump(PAT.dish_pattern(p['D'], p['f0'], r['gain']), fh)
    D, f = p['D'], r['f']
    R = D / 2
    # arc length along the profile from the vertex, for the gore template
    rr = np.linspace(0, R, 400)
    z = rr ** 2 / (4 * f)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(rr), np.diff(z)))])
    n = int(p['gores'])
    half = np.pi * rr / n  # half width of a gore at radius r (circumference / n / 2)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 6))
    a1.plot(np.concatenate([-rr[::-1], rr]), np.concatenate([z[::-1], z]), color='#2a78d6', lw=2)
    a1.plot([0], [f], 'o', color='#e8710a'); a1.text(0, f, '  focus f = %.1f mm' % f, va='center')
    a1.set_aspect('equal'); a1.grid(alpha=0.3)
    a1.set_title('Profile: D %.0f mm, f/D %.2f, depth %.1f mm' % (D, p['fd'], r['depth']))
    a2.fill_betweenx(s, -half, half, color='#e8a33a', alpha=0.6); a2.plot(half, s, color='#9a5d00'); a2.plot(-half, s, color='#9a5d00')
    a2.set_aspect('equal'); a2.grid(alpha=0.3)
    a2.set_title('Gore (%d pcs): width vs. length along the surface [mm]' % n)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'render.png'), dpi=105); plt.close(fig)
    # 1:1 printable SVG of one gore
    sv = [(float(x), float(y)) for x, y in zip(half, s)]
    outline = sv + [(-x, y) for x, y in reversed(sv)]
    W_, H_ = 2 * half.max() + 20, s.max() + 20
    path = 'M ' + ' L '.join('%.2f %.2f' % (x + W_ / 2, H_ - 10 - y) for x, y in outline) + ' Z'
    with open(os.path.join(out, 'gore_template_1_1.svg'), 'w', encoding='utf-8') as fh:
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="%.1fmm" height="%.1fmm" viewBox="0 0 %.1f %.1f">'
                 '<path d="%s" fill="none" stroke="black" stroke-width="0.3"/>'
                 '<text x="5" y="8" font-size="5">Gore 1:1, %d pcs, D %.0f mm, f/D %.2f</text></svg>'
                 % (W_, H_, W_, H_, path, n, D, p['fd']))
    prof = ['  r %7.1f mm   depth %7.2f mm   length along surface %7.1f mm   gore width %7.1f mm'
            % (x, x * x / (4 * f), np.interp(x, rr, s), 2 * np.interp(x, rr, half)) for x in np.linspace(0, R, 11)]
    text = '\n'.join([
        'Anten-O-Matic: Dish %.3f MHz' % p['f0'], '=' * 40,
        'Diameter D = %.0f mm, f/D = %.2f, focal length f = %.1f mm, depth %.1f mm' % (D, p['fd'], f, r['depth']),
        'λ = %.1f mm, gain %.1f dBi (efficiency %.0f %%), -3 dB beamwidth ~%.1f°' % (r['lam'], r['gain'], p['eta'] * 100, r['hpbw']),
        'The feed must "see" a full angle of %.0f° (dish rim seen from the focus).' % r['psi'], '',
        'Profile (z = r² / 4f):'] + prof + ['',
        'Cantenna feed (circular waveguide, TE11 mode):',
        '  inner diameter %.1f mm, guide wavelength %.1f mm' % (r['dwg'], r['lg']),
        '  probe (connector centre) %.1f mm long, %.1f mm from the back wall, can at least %.0f mm long' % (r['probe'], r['lg'] / 4, 0.75 * r['lg'] + r['lg'] / 4),
        '  can mouth at the focus (the phase centre is slightly inside the mouth, fine-tune by moving it).',
        '  Note: below f/D ~0.35 a cantenna is too narrow to illuminate the dish; a choke collar or a shorter f helps.',
        '',
        'Template: print gore_template_1_1.svg at 100 %% scale, %d pieces. No STEP model is generated for the dish.' % int(p['gores'])]) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as fh: fh.write(text)
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as fh: json.dump({'type': 'dish', 'params': p, 'result': r}, fh, indent=1, default=str)
    return out, text


# ------------------------------------------------------------------ open a classic antenna in the CAD

CAD_MAT = {'copper': 'Copper', 'brass': 'Brass', 'silver': 'Silver', 'aluminium': 'Aluminium', 'steel': 'Steel', 'ideal': 'Perfect conductor'}


def to_cad(t, p):
    """the classic antenna as a CAD document: every conductor a bent rod, every feed gap merged into
    the wire(s) it connects (so NEC sees the same gap), touching ends joined. Each rod gets free-form
    evolution freedom so the design can be evolved further in the CAD."""
    st = build(t, p)
    lam = p['_lam']
    tol = 1e-4
    conds = [[name, [np.asarray(q, float) * 1000 for q in pts], None] for name, pts in st['cond']]   # name, pts mm, feed_at
    same = lambda a, b: float(np.linalg.norm(a - b)) < tol
    def touch(pt):
        for i, c in enumerate(conds):
            if c[2] is not None: continue
            if same(c[1][0], pt): return i, 0
            if same(c[1][-1], pt): return i, -1
        return None, None
    arclen = lambda P: float(sum(np.linalg.norm(P[i + 1] - P[i]) for i in range(len(P) - 1)))
    for a, b in st['feeds']:
        a, b = np.asarray(a, float) * 1000, np.asarray(b, float) * 1000
        mid = (a + b) / 2
        ia, ea = touch(a)
        ib, eb = touch(b)
        if ia is not None and ib is not None and ia != ib:
            A = conds[ia][1] if ea == -1 else conds[ia][1][::-1]          # ends at a
            B = conds[ib][1] if eb == 0 else conds[ib][1][::-1]           # starts at b
            P = list(A) + [mid] + list(B)
            fa = arclen(P[:len(A) + 1]) / arclen(P)
            name = conds[ia][0].replace(' left half', '').replace(' right half', '')
            conds[ia] = [name, P, fa]; conds[ib] = None
            conds = [c for c in conds if c is not None]
        elif ia is not None and ia == ib:                                  # one conductor from a round to b
            P = conds[ia][1] if ea == -1 else conds[ia][1][::-1]          # ... ends at a
            if not same(P[0], b): P = P[::-1]
            conds[ia] = [conds[ia][0], [mid] + list(P) + [mid], 0.0]
        elif ia is not None or ib is not None:
            i, e, far = (ia, ea, b) if ia is not None else (ib, eb, a)
            P = conds[i][1] if e == 0 else conds[i][1][::-1]              # starts at the touching point
            conds[i] = [conds[i][0], [far] + list(P), 0.0]
        else:
            conds.append(['Feed', [a, b], 0.5])
    wd = p['wire_d'] * 1000
    els = []
    for k, (name, P, fa) in enumerate(conds):
        P = np.array(P)
        els.append({'id': 'e%d' % (k + 1), 'type': 'rod', 'name': name[:40], 'pos': [round(float(v), 4) for v in P[0]], 'rot': [0, 0, 0],
                    'params': {'pts': (P - P[0]).round(4).tolist(), 'bend_r': p.get('bend_r', 0)}, 'mat': CAD_MAT.get(p.get('material', 'copper'), 'Copper'),
                    'd': wd, 'feed': fa is not None, 'feed_at': 0.5 if fa is None else float(fa), 'fixed': False,
                    'evolve': {'mode': 'free', 'dev': round(lam / 50, 3), 'pts': min(8, max(3, len(P)))}})
    joins = []
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            for na, ia in (('start', 0), ('end', -1)):
                for nb, ib in (('start', 0), ('end', -1)):
                    if same(np.array(conds[i][1][ia]), np.array(conds[j][1][ib])):
                        joins.append({'a': {'el': els[i]['id'], 'pt': na}, 'b': {'el': els[j]['id'], 'pt': nb}})
    doc = {'name': str(p.get('name') or t)[:34] + '-cad', 'f0': p['f0'], 'bw': p['bw'], 'z0': p['z0'], 'swr': p['swr'], 'segs': 25,
           'ground': 'perfect' if st['ground'] else 'none', 'step': True, 'engine': 'nec', 'elements': els, 'joins': joins, 'dims': []}
    notes = []
    if len(st['feeds']) > 1: notes.append('This antenna has %d feeds: in the CAD they are driven in parallel with the same voltage.' % len(st['feeds']))
    if st['solids']: notes.append('The reflector disc is modelled as the perfect ground plane.')
    return doc, notes
