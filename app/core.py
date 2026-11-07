# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""
Anten-O-Matic core: evolved wire antennas with NEC-2 (necpp) + genetic algorithm.

All lengths in metres.  The antenna axis is +z: with a ground plane the antenna stands on
the plane (z >= 0); in free space it is centred on the feed point at the origin.

Modes
  mono   one bent arm over a ground plane, fed against the plane (one coax)
  dipole one bent arm + a copy rotated 180 deg, fed in the middle (one coax, free space)
  self   two different bent dipoles crossed at 90 deg, in parallel on one coax: QFH-style
         self-phasing, the 90 deg shift has to come from the geometry (free space)
  quad   four identical arms over a ground plane, fed 0/90/180/270 deg by an ideal
         phasing network (NASA ST5 style)
  asym   dipole with two independently evolved arms (no symmetry, so it can beam sideways)

Any mode can get npar passive (parasitic) bent wires, evolved together with the arms.
Coverage: 'cone' = theta 0..cone around +z; 'sector' = azimuth az1..az2 (deg from +x towards
+y) x elevation el1..el2 (deg above the xy plane, 90 = +z), optionally with a required
suppression of everything outside the sector (rej dB, measured 15 deg away from its edges).
"""
import json, math, os, time
import numpy as np
import pattern as PAT

C = 299.792458  # m * MHz

MATERIALS = {'copper': 5.8e7, 'brass': 1.5e7, 'aluminium': 3.5e7, 'steel': 1.4e6, 'silver': 6.3e7, 'ideal': 0.0}
MODES = ('mono', 'dipole', 'self', 'quad', 'asym')

DEFAULTS = dict(
    name='antenna', f0=436.0, bw=6.0, mode='self', pol='rhcp', cone=70.0,
    rmax=0.0, hmax=0.0, maxlen=0.0,           # 0 = automatic (from wavelength)
    wire_d=0.002, bend_r=0.0, material='copper', z0=50.0, swr=1.8, step=True,
    k=5, pop=60, gens=80, seed=0, polish=300, jobs=0, islands=4, migrate=8, patience=12, stop_after=0, maxsim=85.0,
    target=None, cov='cone', az1=-30.0, az2=30.0, el1=20.0, el2=90.0, rej=8.0, guard=20.0, npar=0, below=False,
)


def lam(d):
    return C / d['f0']


def normalize(d):
    """fill automatic values, return a new dict."""
    d = dict(DEFAULTS, **{k: v for k, v in d.items() if v is not None})
    for k in ('f0', 'bw', 'cone', 'rmax', 'hmax', 'maxlen', 'wire_d', 'bend_r', 'z0', 'swr', 'az1', 'az2', 'el1', 'el2', 'rej', 'guard', 'maxsim'):
        d[k] = float(d[k])
    for k in ('k', 'pop', 'gens', 'seed', 'polish', 'jobs', 'npar', 'islands', 'migrate', 'patience', 'stop_after'):
        d[k] = int(float(d[k]))
    if d['f0'] <= 0 or d['wire_d'] <= 0: raise ValueError('frequency and wire diameter must be positive')
    if d['bw'] < 0 or d['k'] < 1 or d['pop'] < 4 or d['gens'] < 1: raise ValueError('invalid evolution parameters')
    L = lam(d)
    if d['mode'] not in MODES: raise ValueError('unknown mode %s' % d['mode'])
    if d['pol'] not in ('rhcp', 'lhcp', 'total'): raise ValueError('pol must be rhcp, lhcp or total')
    if d['cov'] not in ('cone', 'sector', 'target'): raise ValueError('cov must be cone, sector or target')
    if d['cov'] == 'target':
        d['target'] = PAT.clean_target(d['target'] or PAT.default_target())
        d['_taz'], d['_tel'] = PAT.target_peak(d['target'])
    if d['cov'] == 'sector':
        if d['el2'] < d['el1']: d['el1'], d['el2'] = d['el2'], d['el1']
        d['el1'], d['el2'] = max(-90.0, d['el1']), min(90.0, d['el2'])
        if ground(d) and d['el1'] < 0: d['el1'] = 0.0
        if d['az2'] <= d['az1']: d['az2'] += 360.0
        if d['az2'] - d['az1'] > 360: d['az2'] = d['az1'] + 360.0
    d['npar'] = max(0, min(4, d['npar']))
    big = d['cov'] in ('sector', 'target') and d['npar'] > 0  # beams with passive elements need room
    if not d['rmax']: d['rmax'] = (0.35 if big else 0.25) * L
    if not d['hmax']: d['hmax'] = (0.3 if ground(d) else 0.5) * L * (1.4 if big else 1.0)
    if not d['maxlen']: d['maxlen'] = 1.5 * L
    d['radius'] = d['wire_d'] / 2
    if d['bend_r'] <= 0: d['bend_r'] = d['wire_d']        # inner bend radius, default 1x wire diameter
    d['rc'] = d['bend_r'] + d['radius']                    # bend radius of the wire axis
    d['step'] = str(d['step']).lower() not in ('0', 'false', 'no', '')
    d['below'] = str(d['below']).lower() not in ('0', 'false', 'no', '')  # True: suppress also below the horizon
    d['sigma'] = MATERIALS.get(d['material'], 5.8e7)
    d['r0'] = 0.0 if d['mode'] == 'mono' else max(0.02 * L, 1.5 * d['wire_d'])  # feed offset
    d['stub'] = max(0.02 * L, 3 * d['wire_d'])
    d['gap'] = max(0.02 * L, 3 * d['wire_d'])        # self: height offset of the two feeds
    d['clear'] = max(0.01 * L, d['wire_d']) + d['wire_d']  # centre-to-centre wire clearance
    d['minseg'] = max(0.02 * L, 3 * d['wire_d'])
    d['segs'] = 25
    d['jobs'] = d['jobs'] or max(1, (os.cpu_count() or 2) - 1)
    return d


def ground(d):
    return d['mode'] in ('mono', 'quad')


def sector_az_limited(d):
    return d.get('cov') == 'sector' and d['az2'] - d['az1'] < 359.0


def coverage_text(d):
    if d['cov'] == 'target':
        return 'drawn target pattern (peak towards azimuth %.0f°, elevation %.0f°, floor %.0f dB)' % (d['_taz'], d['_tel'], d['target']['floor'])
    if d['cov'] == 'cone':
        return 'cone ±%.0f° around the z axis' % d['cone']
    az = 'azimuth %+.0f° to %+.0f°' % (d['az1'], d['az2'] if d['az2'] <= 180 else d['az2'] - 360) if sector_az_limited(d) else 'full azimuth'
    s = 'sector: %s, elevation %.0f° to %.0f°' % (az, d['el1'], d['el2'])
    if d['rej'] > 0: s += ', suppress outside the sector (%.0f° from its edge%s) by %.0f dB' % (d['guard'], '' if d['below'] else ', above the horizon only', d['rej'])
    return s


def warnings(d):
    w = []
    L = lam(d)
    if d['wire_d'] > L / 30:
        w.append('The wire is very thick for this frequency (diameter > λ/30). The NEC-2 thin-wire model becomes inaccurate.')
    if d['mode'] == 'quad':
        w.append('Quad mode needs a 0/90/180/270° phasing network (hybrids). It does not work without one.')
    if sector_az_limited(d) and d['mode'] in ('dipole', 'self', 'quad') and not d['npar']:
        w.append('This type is symmetric about its axis (arm + rotated copy), so it radiates equally forward and backward '
                 'and physically cannot meet an azimuth sector. For azimuth directivity use the asymmetric dipole, a single arm '
                 'over ground, or add passive elements.')
    if d['mode'] in ('dipole', 'self', 'asym'):
        w.append('Balanced feed: put a choke (ferrites) or a balun on the coax.')
    if ground(d):
        w.append('The model assumes an infinite ground plane. A real plane should be at least ~λ/2 = %.0f mm across.' % (L / 2 * 1000))
    return w

# ---------------------------------------------------------------- geometry

def rotate(p, ang):
    c, s = math.cos(ang), math.sin(ang)
    return np.array([c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]])


def n_arms_genome(d):
    return 2 if d['mode'] in ('self', 'asym') else 1


def n_par_points(d):
    return max(2, d['k'])  # points of one passive element (k-1 bends)


def bounds(d):
    R, H = d['rmax'], d['hmax']
    zlo, zhi = (0.01 * lam(d) + d['wire_d'], H) if ground(d) else (-H / 2, H / 2)
    n = d['k'] * n_arms_genome(d) + d.get('npar', 0) * n_par_points(d)
    return np.tile([-R, -R, zlo], n), np.tile([R, R, zhi], n)


def parasitics(p, d):
    """passive elements: [(arm_id, [points])], ids 10, 11, ..."""
    base, kp = d['k'] * n_arms_genome(d), n_par_points(d)
    return [(10 + j, list(p[base + j * kp: base + (j + 1) * kp])) for j in range(d.get('npar', 0))]


def arm_paths(g, d):
    """raw polylines of every physical arm: [(arm_id, first_is_stub, [points])] + feed wires.

    ground modes: base on the plate -> top of the vertical stub -> evolved points
    free modes:   feed terminal -> evolved points (the terminal is a solder joint)"""
    p = list(g.reshape(-1, 3))
    r0, arms, feeds = d['r0'], [], []
    if ground(d):
        n = 1 if d['mode'] == 'mono' else 4
        pts = [np.array([r0, 0.0, 0.0]), np.array([r0, 0.0, d['stub']])] + p[:d['k']]
        for k in range(n):
            arms.append((k, True, [rotate(q, 2 * math.pi * k / n) for q in pts]))
        arms += [(k, False, q) for k, q in parasitics(p, d)]
        return arms, feeds
    if d['mode'] == 'asym':
        feeds.append(('f0', np.array([-r0, 0, 0.0]), np.array([r0, 0, 0.0])))
        arms.append((0, False, [np.array([r0, 0.0, 0.0])] + p[:d['k']]))
        arms.append((1, False, [np.array([-r0, 0.0, 0.0])] + p[d['k']:2 * d['k']]))
        arms += [(k, False, q) for k, q in parasitics(p, d)]
        return arms, feeds
    halves = [(p[:d['k']], 0.0)]
    if d['mode'] == 'self': halves.append((p[d['k']:2 * d['k']], d['gap']))
    for j, (pp, z) in enumerate(halves):
        base = j * math.pi / 2
        feeds.append(('f%d' % j, rotate(np.array([-r0, 0, z]), base), rotate(np.array([r0, 0, z]), base)))
        pts = [np.array([r0, 0.0, z])] + list(pp)
        for h in (0, 1):
            arms.append((j + 2 * h, False, [rotate(q, base + h * math.pi) for q in pts]))
    arms += [(k, False, q) for k, q in parasitics(p, d)]
    return arms, feeds


def fillet(pts, rc, min_angle=0.5):
    """round every inner vertex of a polyline with radius rc (wire axis).
    Vertices turning by less than min_angle degrees stay sharp (a polyline approximating a
    smooth curve, e.g. a helix, needs no extra arcs).

    Returns [('line', A, B) | ('arc', T1, T2, centre, axis, angle)] or None if a bend does
    not fit (the tangent lengths of two neighbouring bends overlap)."""
    pts = [np.asarray(q, float) for q in pts]
    n = len(pts)
    dirs, lens = [], []
    for i in range(n - 1):
        v = pts[i + 1] - pts[i]
        L = float(np.linalg.norm(v))
        if L < 1e-9: return None
        dirs.append(v / L); lens.append(L)
    tan, bends, kinks = [0.0] * n, [None] * n, [False] * n
    for i in range(1, n - 1):
        a, b = dirs[i - 1], dirs[i]
        th = math.acos(float(np.clip(a @ b, -1, 1)))
        if th < math.radians(0.5): continue
        if th < math.radians(min_angle):
            kinks[i] = True; continue
        if th > math.radians(175): return None
        tan[i] = rc * math.tan(th / 2)
        bends[i] = th
    for i in range(n - 1):
        if tan[i] + tan[i + 1] > lens[i] + 1e-9: return None
    out = []
    cur = pts[0]
    for i in range(1, n):
        if i < n - 1 and bends[i] is not None:
            a, b = dirs[i - 1], dirs[i]
            t1 = pts[i] - a * tan[i]
            t2 = pts[i] + b * tan[i]
            if np.linalg.norm(t1 - cur) > 1e-9: out.append(('line', cur, t1))
            inward = b - (a @ b) * a
            inward /= np.linalg.norm(inward)
            centre = t1 + inward * rc
            axis = np.cross(a, b); axis /= np.linalg.norm(axis)
            out.append(('arc', t1, t2, centre, axis, bends[i]))
            cur = t2
        elif i == n - 1 or kinks[i]:
            out.append(('line', cur, pts[i]))
            cur = pts[i]
    return out


def arc_points(t1, centre, axis, angle, n):
    v = t1 - centre
    return [centre + v * math.cos(angle * k / n) + np.cross(axis, v) * math.sin(angle * k / n) + axis * (axis @ v) * (1 - math.cos(angle * k / n))
            for k in range(n + 1)]


def all_wires(g, d):
    """[(arm, i, P, Q)] with rounded bends; i == 0 marks feed wires (stubs / feed gaps), arc
    pieces get half-integer i.  None if a bend radius does not fit."""
    arms, feeds = arm_paths(g, d)
    out = [(k, 0, p, q) for (k, p, q) in feeds]
    for k, stub, pts in arms:
        # a bend whose arc would be shorter than 2.5 wire diameters is modelled as a sharp corner:
        # NEC rejects segments that short next to a thick wire, and electrically it makes no difference
        pieces = fillet(pts, d['rc'], max(10.0, min(179.0, math.degrees(2.5 * d['wire_d'] / d['rc']))))
        if pieces is None: return None
        seg = 0 if stub else 1
        for pc in pieces:
            if pc[0] == 'line':
                out.append((k, seg, pc[1], pc[2])); seg += 1
            else:
                _, t1, t2, c, ax, th = pc
                m = max(1, min(math.ceil(math.degrees(th) / 45), int(d['rc'] * th / (2 * d['wire_d']))))
                ap = arc_points(t1, c, ax, th, m)
                out += [(k, seg - 0.5, ap[j], ap[j + 1]) for j in range(m)]
    return out


def seg_dist(p1, q1, p2, q2):
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    c, b = d1 @ r, d1 @ d2
    den = a * e - b * b
    s = np.clip((b * f - c * e) / den, 0, 1) if den > 1e-15 else 0.0
    t = (b * s + f) / e
    if t < 0: t, s = 0.0, np.clip(-c / a, 0, 1)
    elif t > 1: t, s = 1.0, np.clip((b - c) / a, 0, 1)
    return float(np.linalg.norm(p1 + d1 * s - p2 - d2 * t))


def geometry_ok(g, d):
    arms, _ = arm_paths(g, d)
    for k, stub, pts in arms:
        L = 0.0
        for j in range(len(pts) - 1):
            seg = float(np.linalg.norm(pts[j + 1] - pts[j]))
            if (j > 0 or not stub) and seg < d['minseg']: return False
            L += seg
        if d['maxlen'] and L > d['maxlen']: return False
        if any(math.hypot(q[0], q[1]) > d['rmax'] * 1.0001 for q in pts): return False
    wires = all_wires(g, d)
    if wires is None: return False
    n = len(wires)
    for i in range(n):
        for j in range(i + 1, n):
            k1, i1, p1, q1 = wires[i]
            k2, i2, p2, q2 = wires[j]
            if k1 == k2 and abs(i1 - i2) <= 1.5: continue  # neighbours around one bend
            if min(np.linalg.norm(u - v) for u in (p1, q1) for v in (p2, q2)) < 1e-9: continue
            if seg_dist(p1, q1, p2, q2) < d['clear']: return False
    return True


def nsegs(L, d, fmax, fixed):
    """segments per wire: ~25/lambda, but never shorter than 2 wire diameters (thin-wire limit)."""
    want = math.ceil(L / (C / fmax / d['segs']))
    cap = max(1, int(L / (2 * d['wire_d'])))
    n = max(1, min(want, cap))
    if fixed: n = 3 if cap >= 3 else 1  # odd, source sits on the middle segment
    return n

# ---------------------------------------------------------------- NEC

def _chk(r):
    if r != 0:
        from necwrap import nec_error_message
        raise RuntimeError(nec_error_message().decode(errors='replace'))


def build(n, g, d, freqs, rp):
    import necwrap as N
    wires = all_wires(g, d)
    feeds = []
    for t, (k, i, p, q) in enumerate(wires, start=1):
        s = nsegs(np.linalg.norm(q - p), d, max(freqs), i == 0)
        _chk(N.nec_wire(n, t, s, *p, *q, d['radius'], 1.0, 1.0))
        if i == 0: feeds.append((t, (s + 1) // 2))
    if ground(d):
        _chk(N.nec_geometry_complete(n, 1))
        _chk(N.nec_gn_card(n, 1, 0, 0, 0, 0, 0, 0, 0))
    else:
        _chk(N.nec_geometry_complete(n, 0))
    _chk(N.nec_ek_card(n, 0))  # extended thin-wire kernel: short segments in the bends
    if d['sigma']: _chk(N.nec_ld_card(n, 5, 0, 0, 0, d['sigma'], 0, 0))
    nf = len(freqs)
    _chk(N.nec_fr_card(n, 0, nf, freqs[0], (freqs[-1] - freqs[0]) / (nf - 1) if nf > 1 else 0))
    if d['mode'] == 'quad':
        sign = 1 if d['pol'] == 'lhcp' else -1
        for k, (t, s) in enumerate(feeds):
            ph = sign * math.pi / 2 * k
            _chk(N.nec_excitation_voltage(n, t, 1, math.cos(ph), math.sin(ph)))
    elif d['mode'] == 'mono':
        _chk(N.nec_excitation_voltage(n, feeds[0][0], 1, 1.0, 0.0))
    else:
        if d['mode'] == 'self':  # dipole B in parallel with A via a ~zero-length line
            _chk(N.nec_tl_card(n, feeds[0][0], feeds[0][1], feeds[1][0], feeds[1][1], d['z0'], 1e-5, 0, 0, 0, 0))
        _chk(N.nec_excitation_voltage(n, feeds[0][0], feeds[0][1], 1.0, 0.0))
        feeds = feeds[:1]
    nth, nph, th0, ph0, dth, dph = rp
    _chk(N.nec_rp_card(n, 0, nth, nph, 0, 0, 0, 0, th0, ph0, dth, dph, 0, 0))
    return len(feeds)


def _wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def _unit(th, ph):
    th, ph = np.radians(th), np.radians(ph)
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], -1)


def _in_sector(d, th, ph):
    el = 90.0 - th
    azc = (d['az1'] + d['az2']) / 2; half = (d['az2'] - d['az1']) / 2
    pole = (th < 1e-6) | (th > 180 - 1e-6)  # azimuth is meaningless at the zenith / nadir
    az_ok = True if not sector_az_limited(d) else (np.abs(_wrap(ph - azc)) <= half + 1e-6) | pole
    return (el >= d['el1'] - 1e-6) & (el <= d['el2'] + 1e-6) & az_ok


_MASKS = {}


def sector_masks(d, th, ph):
    """boolean masks (inside the sector, farther than GUARD deg from it) for grids of theta/phi in deg.
    The distance is the true angle on the sphere, so the zenith region is handled right."""
    key = (d['az1'], d['az2'], d['el1'], d['el2'], d['below'], d['guard'], th.shape, float(th.max()))
    if key not in _MASKS:
        ins = _in_sector(d, th, ph)
        # dense sample of the sector for the distance test
        T, P = np.meshgrid(np.arange(0, 180.1, 2.0), np.arange(0, 360, 2.0), indexing='ij')
        sv = _unit(T[_in_sector(d, T, P)], P[_in_sector(d, T, P)])
        gv = _unit(th, ph).reshape(-1, 3)
        dist = np.degrees(np.arccos(np.clip(gv @ sv.T, -1, 1))).min(1).reshape(th.shape)
        out = dist > d['guard']
        if not d['below']: out &= th <= 90 + 1e-6  # below the horizon is ground anyway: do not care
        _MASKS.clear(); _MASKS[key] = (ins, out)
    return _MASKS[key]


def sector_grid(d):
    span = min(d['el2'] - d['el1'], d['az2'] - d['az1'])
    return 5 if span < 40 else 10


def simulate_sector(g, d, freqs):
    """per-frequency floor / mean / max gain inside the sector and max total gain outside."""
    import necwrap as N
    st = sector_grid(d)
    thmax = 90 if ground(d) else 180
    nth, nph = thmax // st + 1, 360 // st
    TH, PH = np.meshgrid(np.arange(nth) * st, np.arange(nph) * st, indexing='ij')
    ins, outs = sector_masks(d, TH.astype(float), PH.astype(float))
    n = N.nec_create()
    try:
        nfeed = build(n, g, d, freqs, (nth, nph, 0, 0, st, st))
        res = []
        for fi, f in enumerate(freqs):
            z = [complex(N.nec_impedance_real_indexed(n, fi, k), N.nec_impedance_imag_indexed(n, fi, k)) for k in range(nfeed)]
            G = np.array([[N.nec_gain(n, fi, i, j) for j in range(nph)] for i in range(nth)])
            res.append({'f': f, 'z': z, 'grid': G,
                        'outmax': float(G[outs].max()) if outs.any() else -999.0,
                        'floor': float(G[ins].min()), 'mean': float(10 * np.log10(np.mean(10 ** (G[ins] / 10)))),
                        'max': float(G[ins].max())})
    finally:
        N.nec_delete(n)
    if d['pol'] != 'total':
        # NEC only gives polarised gain as statistics over one RP block: rerun on the sector itself
        el1, el2 = d['el1'], d['el2']
        th0, th1 = 90 - el2, 90 - el1
        ntt = max(2, int(math.ceil((th1 - th0) / st)) + 1)
        npp = max(2, int(math.ceil((d['az2'] - d['az1']) / st)) + 1)
        if not sector_az_limited(d): npp = 360 // st
        rp = (ntt, npp, th0, d['az1'], (th1 - th0) / (ntt - 1), (d['az2'] - d['az1']) / (npp - 1) if sector_az_limited(d) else st)
        for r, r2 in zip(res, simulate(g, d, freqs, rp)):
            r.update(floor=r2['floor'], mean=r2['mean'], max=r2['max'])
    return res


_TGRID = {}


def target_grid(d, step):
    key = (json.dumps(d['target'], sort_keys=True), step, ground(d))
    if key not in _TGRID:
        th, ph = PAT.grid_axes(step, ground(d))
        _TGRID.clear(); _TGRID[key] = (th, ph, PAT.target_grid(d['target'], th, ph))
    return _TGRID[key]


def simulate_target(g, d, freqs, step=10):
    """per-frequency RMS difference between the normalised pattern and the drawn target."""
    import necwrap as N
    th, ph, T = target_grid(d, step)
    n = N.nec_create()
    try:
        nfeed = build(n, g, d, freqs, (len(th), len(ph), 0, 0, step, step))
        res = []
        for fi, f in enumerate(freqs):
            z = [complex(N.nec_impedance_real_indexed(n, fi, k), N.nec_impedance_imag_indexed(n, fi, k)) for k in range(nfeed)]
            G = PAT.read_grid(N, n, fi, len(th), len(ph), d['pol'])
            err, _ = PAT.match(G, T, th, d['target']['floor'], ground(d))
            res.append({'f': f, 'z': z, 'err': err, 'max': float(G.max()),
                        'floor': float(G.max()), 'mean': float(10 * np.log10(np.mean(10 ** (G / 10))))})
        return res
    finally:
        N.nec_delete(n)


def full_pattern(g, d, f=None, step=5):
    """full-sphere pattern (total, RHCP, LHCP) at one frequency, for the 3D view and the plots."""
    import necwrap as N
    n = N.nec_create()
    try:
        build(n, g, d, [f or d['f0']], PAT.rp_for(step, ground(d)))
        return PAT.from_nec(N, n, 0, step, ground(d))
    finally:
        N.nec_delete(n)


def simulate(g, d, freqs, rp=None):
    import necwrap as N
    if rp is None and d.get('cov') == 'sector':
        return simulate_sector(g, d, freqs)
    if rp is None and d.get('cov') == 'target':
        return simulate_target(g, d, freqs)
    if rp is None:  # coverage cone around +z: theta 0..cone (step 10), phi 0..345 (step 15)
        rp = (int(d['cone'] // 10) + 1, 24, 0, 0, 10, 15)
    n = N.nec_create()
    try:
        nfeed = build(n, g, d, freqs, rp)
        res = []
        for fi, f in enumerate(freqs):
            z = [complex(N.nec_impedance_real_indexed(n, fi, k), N.nec_impedance_imag_indexed(n, fi, k))
                 for k in range(nfeed)]
            res.append({'f': f, 'z': z,
                        'floor': N.gain(n, d['pol'], 'min', fi),
                        'mean': N.gain(n, d['pol'], 'mean', fi),
                        'max': N.gain(n, d['pol'], 'max', fi)})
        return res
    finally:
        N.nec_delete(n)


def swr(z, z0):
    g = abs((z - z0) / (z + z0))
    return (1 + g) / (1 - g) if g < 0.999 else 999.0

# ---------------------------------------------------------------- fitness

_D = None


def _init(d):
    global _D
    _D = d


def fitness(g, d=None):
    d = d or _D
    if not geometry_ok(g, d): return -1000.0, None
    try:
        res = simulate(g, d, [d['f0'] - d['bw'] / 2, d['f0'], d['f0'] + d['bw'] / 2])
    except Exception:
        return -1000.0, None
    floor = min(r['floor'] for r in res)
    mean = min(r['mean'] for r in res)
    s = max(swr(z, d['z0']) for r in res for z in r['z'])
    over = max(0.0, s - d['swr'])
    f = max(floor, -40) + 0.25 * mean - 3 * over - 0.2 * over ** 2
    info = {'floor': floor, 'mean': mean, 'swr': s}
    if d.get('cov') == 'target':
        err = max(r['err'] for r in res)
        # mismatch loss in dB instead of the SWR excess: it sits on the same dB scale as the pattern
        # error, so a good shape with a poor match is not thrown away for a matched blob
        ml = lambda q: -10 * math.log10(max(1e-6, 1 - ((q - 1) / (q + 1)) ** 2))
        f = -err - 2.0 * max(0.0, ml(min(s, 999)) - ml(d['swr']))
        info = {'err': err, 'floor': floor, 'mean': mean, 'swr': s}
        return f, info
    if d.get('cov') == 'sector':
        sup = min(r['max'] - r['outmax'] for r in res)  # main lobe in the sector vs strongest unwanted lobe (like F/B)
        info['sup'] = sup
        if d['rej'] > 0:
            f -= 1.0 * max(0.0, d['rej'] - sup)
    return f, info


def _eval(g):
    return fitness(g)


def _eval_many(gs):
    return [fitness(g) for g in gs]

# ---------------------------------------------------------------- GA

def project(g, d):
    """pull points back into the cylinder of radius rmax."""
    p = g.reshape(-1, 3).copy()
    r = np.hypot(p[:, 0], p[:, 1])
    s = np.where(r > d['rmax'], d['rmax'] / np.maximum(r, 1e-12), 1.0)
    p[:, 0] *= s; p[:, 1] *= s
    return p.ravel()


def random_walk(rng, d, lo, hi):
    """random arm(s) as a smooth-ish random walk from the feed point (tangles less than
    uniformly scattered points, so tight boxes and thick wires still give valid shapes)."""
    L = lam(d)
    g = np.empty(lo.size)
    walks = []
    for a in range(n_arms_genome(d)):
        z = d['stub'] if ground(d) else (d['gap'] if a and d['mode'] == 'self' else 0.0)
        sgn = -1.0 if (a and d['mode'] == 'asym') else 1.0
        walks.append((a * d['k'], d['k'], np.array([sgn * d['r0'], 0.0, z]), sgn))
    kp, base = n_par_points(d), d['k'] * n_arms_genome(d)
    for j in range(d.get('npar', 0)):
        # passive elements start as gently bent ~lambda/2 wires 0.1..0.35 lambda away from the feed,
        # roughly where a reflector or director would sit
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        start = project(np.clip(u * rng.uniform(0.1, 0.35) * L + rng.normal(0, 0.12 * L, 3), lo[:3], hi[:3]), d)
        g[(base + j * kp) * 3:(base + j * kp) * 3 + 3] = start
        walks.append((base + j * kp + 1, kp - 1, start, None))
    for i0, npts, p, sgn in walks:
        p = p.copy()
        par = sgn is None
        v = rng.normal(size=3)
        if not par: v[0] = sgn * (abs(v[0]) + 0.5)
        for i in range(npts):
            v = v / np.linalg.norm(v) + rng.normal(0, 0.3 if par else 0.9, 3)
            step = (rng.uniform(0.4, 0.55) * L / max(1, npts)) if par else (rng.uniform(0.04, 0.16) * L + d['minseg'])
            step = max(step, d['minseg'] * 1.05)
            p = p + v / np.linalg.norm(v) * step
            j = (i0 + i) * 3
            p = project(np.clip(p, lo[j:j + 3], hi[j:j + 3]), d)
            g[j:j + 3] = p
    return g

def shape_sig(g, d, n=16):
    """shape signature for similarity: every wire resampled to n points spaced evenly along it."""
    out = []
    for _, _, pts in arm_paths(g, d)[0]:
        P = np.array(pts)
        s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        t = np.linspace(0, s[-1], n)
        out.append(np.stack([np.interp(t, s, P[:, k]) for k in range(3)], 1))
    return np.concatenate(out)


def similarity(a, b, d):
    """0..100 %: 100 = identical wires, 0 = shapes differ by 0.15 lambda RMS or more."""
    rms = float(np.sqrt(np.mean(np.sum((a - b) ** 2, 1))))
    return max(0.0, 1.0 - rms / (0.15 * lam(d))) * 100


def beam_seed(rng, d, lo, hi):
    """Yagi-like starting shape aimed at the sector centre: driven element across the beam,
    passive elements as reflector (behind) and directors (in front), plus noise.  Evolution then
    bends it into whatever works; without these seeds a random population hardly ever finds
    a beam."""
    L = lam(d)
    g = random_walk(rng, d, lo, hi)
    if d['cov'] == 'target':
        azc, elc = math.radians(d['_taz']), math.radians(max(-60.0, min(60.0, d['_tel'])))
    else:
        azc = math.radians((d['az1'] + d['az2']) / 2)
        elc = math.radians(min(60.0, (d['el1'] + d['el2']) / 2) if sector_az_limited(d) else (d['el1'] + d['el2']) / 2)
    u = np.array([math.cos(elc) * math.cos(azc), math.cos(elc) * math.sin(azc), math.sin(elc)])
    if ground(d):
        u = np.array([math.cos(azc), math.sin(azc), 0.0])  # over ground: vertical elements, beam along the plate
        e, half = np.array([0.0, 0.0, 1.0]), 0.24 * L
    else:
        e = np.cross(u, [0.0, 0.0, 1.0]) if abs(u[2]) < 0.95 else np.array([0.0, 1.0, 0.0])
        e /= np.linalg.norm(e); half = 0.245 * L
    k, kp = d['k'], n_par_points(d)
    pts = g.reshape(-1, 3)
    if d['mode'] == 'asym':
        for a, sgn in ((0, 1), (1, -1)):
            t0 = np.array([sgn * d['r0'], 0.0, 0.0])
            for i in range(k):
                pts[a * k + i] = t0 + sgn * e * half * (i + 1) / k * rng.uniform(0.9, 1.1)
    elif d['mode'] == 'mono':
        for i in range(k):
            pts[i] = np.array([d['r0'], 0.0, d['stub'] + (half - d['stub']) * (i + 1) / k])
    base = k * n_arms_genome(d)
    offs = [-0.2] + [0.12 + 0.18 * j for j in range(3)]
    lens = [1.04] + [0.94 - 0.03 * j for j in range(3)]
    for j in range(d['npar']):
        c = u * offs[j] * L * rng.uniform(0.85, 1.15)
        hl = half * lens[j] * rng.uniform(0.95, 1.05)
        for i in range(kp):
            t = i / (kp - 1)
            pts[base + j * kp + i] = (c + e * (0.03 * L + (hl - 0.03 * L) * t) if ground(d) else c + e * hl * (2 * t - 1))
            if ground(d): pts[base + j * kp + i][2] = 0.02 * L + 2 * hl * t * 0.98
    g = pts.ravel() + rng.normal(0, 0.015 * L, g.size)
    return project(np.clip(g, lo, hi), d)


class _Stopped(Exception):
    pass


class Evolver:
    """run() blocks; progress(dict) is called every generation; stop is a threading.Event."""

    def __init__(self, d, progress=None, stop=None, log=None):
        self.d, self.progress, self.stop = normalize(d), progress or (lambda s: None), stop
        self.log = log or (lambda m: None)
        self.best = None
        self.evals = 0
        self._busy = 0.0

    def _stopped(self):
        return self.stop is not None and self.stop.is_set()

    def _new_pool(self):
        from multiprocessing import Pool
        self.pool = Pool(self.d['jobs'], initializer=_init, initargs=(self.d,))

    def _kill_pool(self):
        try: self.pool.terminate()
        except Exception: pass

    def _map(self, pool, items, what):
        """parallel evaluation with live progress, a Stop that works mid-generation and a watchdog:
        if NEC hangs or a worker process dies, no result arrives; after a timeout the pool is
        killed and rebuilt, the shape it was stuck on scores as invalid and the rest is re-run."""
        out, n, t0, last = [], len(items), time.time(), 0.0
        try:
            return self._map_inner(items, what, out, n, last)
        finally:
            self._busy += time.time() - t0

    def _map_inner(self, items, what, out, n, last):
        while len(out) < n:
            todo = items[len(out):]
            chunk = max(1, len(todo) // (self.d['jobs'] * 4)) if not getattr(self, '_careful', 0) else 1
            # chunks are sent by hand: imap with chunksize > 1 gives a plain generator without a timeout
            it = self.pool.imap(_eval_many, [todo[i:i + chunk] for i in range(0, len(todo), chunk)])
            t_last, hung = time.time(), False
            while len(out) < n:
                if self._stopped():
                    self._kill_pool()
                    raise _Stopped()
                # an evaluation normally takes well under a second; allow generously for slow PCs and big chunks
                per = self._busy * self.d['jobs'] / self.evals if self.evals > 20 else 2.0   # s per eval per process
                limit = max(60.0, 20.0 * chunk * per)
                try:
                    r = it.next(timeout=0.5)
                except StopIteration:
                    break
                except Exception as e:
                    if not isinstance(e, __import__('multiprocessing').TimeoutError): raise
                    if time.time() - t_last > limit:
                        hung = True; break
                    now = time.time()
                    if now - last > 2:
                        last = now
                        self.progress(dict(self._last, beat=now, submsg='%s: %d/%d (waiting for NEC %.0f s)' % (what, len(out), n, now - t_last)))
                    continue
                t_last = time.time()
                out += r; self.evals += len(r)
                now = time.time()
                if now - last > 0.5 or len(out) == n:
                    last = now
                    self.progress(dict(self._last, sub=len(out) / n, submsg='%s: %d/%d' % (what, len(out), n),
                                       evals=self.evals, rate=self.evals / max(1e-6, now - self.t0), beat=now))
            if hung:
                self.hangs = getattr(self, 'hangs', 0) + 1
                self.log('Watchdog: NEC gave no result for %.0f s (hung or crashed process), restarting the worker processes. '
                         'The shape it got stuck on counts as invalid.' % (time.time() - t_last))
                self._kill_pool(); self._new_pool()
                if getattr(self, '_careful', 0):
                    out.append((-1000.0, None)); self.evals += 1   # with chunk 1 the stuck one is exactly this
                self._careful = 2  # re-run the rest one by one, so the culprit can be pinned down
            elif getattr(self, '_careful', 0):
                self._careful -= 1
        return out

    def _patfn(self, g):
        if g is None: return None
        g = np.array(g, copy=True); d = self.d
        return lambda: full_pattern(g, d)

    def _report(self, p):
        self._last = dict(p, beat=time.time(), evals=self.evals, rate=self.evals / max(1e-6, time.time() - self.t0))
        self.progress(self._last)

    def run(self):
        from multiprocessing import Pool
        d = self.d
        rng = np.random.default_rng(d['seed'] or None)
        lo, hi = bounds(d)
        t0 = self.t0 = time.time()
        self._report({'phase': 'init', 'gen': 0, 'gens': d['gens'], 'msg': 'Generating the first population'})
        self.log('Start: %s, %.3f MHz ± %.3f, %s, wire Ø %.2f mm, bend radius %.2f mm, %d cores'
                 % (MODE_LABEL[d['mode']], d['f0'], d['bw'] / 2, d['pol'].upper(), d['wire_d'] * 1000,
                    d['bend_r'] * 1000, d['jobs']))
        self.log('Space: radius %.1f mm, height %.1f mm, max %.0f mm of wire per arm, %d bends per arm'
                 % (d['rmax'] * 1000, d['hmax'] * 1000, d['maxlen'] * 1000, d['k'])
                 + (', %d passive elements' % d['npar'] if d['npar'] else ''))
        self.log('Coverage: %s' % coverage_text(d))
        for w in warnings(d):
            if 'symmetric' in w: self.log('WARNING: ' + w)
        seeded = ((d['cov'] == 'sector' and sector_az_limited(d)) or (d['cov'] == 'target' and d['_tel'] < 60)) \
            and (d['npar'] or d['mode'] in ('asym', 'mono'))

        def fresh(n, seed_share=0.0):
            """n valid new genomes (random walks, optionally partly Yagi-style seeds)."""
            out, tries = [], 0
            if seeded and seed_share:
                cand = [beam_seed(rng, d, lo, hi) for _ in range(max(4, int(n * seed_share * 3)))]
                out += [g for g in cand if geometry_ok(g, d)][:int(n * seed_share)]
            while len(out) < n:
                cand = [random_walk(rng, d, lo, hi) for _ in range(max(8, (n - len(out)) * 4))]
                out += [g for g in cand if geometry_ok(g, d)]
                tries += 1
                if tries > 100:
                    raise RuntimeError('No valid shape fits in the given space. Increase the dimensions or use a thinner wire.')
                if self._stopped(): return None
            return out[:n]

        nisl = max(1, min(d['islands'], d['pop'] // 8))   # every island needs a few individuals
        hold = {}   # (a, b) -> generations the two champions have stayed too similar
        size = d['pop'] // nisl
        pop0 = fresh(size * nisl, 0.4)
        if pop0 is None: return None
        rng.shuffle(pop0)
        if seeded: self.log('Part of the initial population are Yagi-style shapes aimed at the %s' % ('target peak' if d['cov'] == 'target' else 'sector centre'))
        self.log('Initial population: %d valid shapes on %d %s of %d' % (len(pop0), nisl, 'islands' if nisl > 1 else 'island', size))
        if nisl > 1:
            self.log('Islands evolve separately; every %d generations each sends its %d best to the neighbouring island' % (d['migrate'], 2))
        if d['patience']:
            self.log('Stagnation watchdog: an island with no improvement for %d generations renews half its population and boosts mutation' % d['patience'])
        hist = []
        self.log('Starting %d worker processes' % d['jobs'])
        self._new_pool(); pool = None
        try:
            r = self._map(pool, pop0, 'Generation 0')
            isl = []
            for j in range(nisl):
                part = slice(j * size, (j + 1) * size)
                isl.append({'pop': np.array(pop0[part]), 'fit': np.array([x[0] for x in r[part]]),
                            'info': [x[1] for x in r[part]], 'best': -1e9, 'stall': 0, 'boost': 0.0, 'restarts': 0})
            best_f, best_g, best_i, gstall = -1e9, None, None, 0
            elite = max(1, size // 15)
            for gen in range(d['gens']):
                improved = False
                for j, I in enumerate(isl):
                    I['order'] = np.argsort(-I['fit'])
                    top = I['order'][0]
                    if I['fit'][top] > I['best'] + 1e-6:
                        I['best'], I['stall'] = I['fit'][top], 0
                        I['bestg'], I['besti'] = I['pop'][top].copy(), I['info'][top]
                        I['wires'] = wires_json(I['bestg'], d) if I['besti'] else None
                    else:
                        I['stall'] += 1
                    if I['fit'][top] > best_f + 1e-6:
                        best_f, best_g, best_i, improved = I['fit'][top], I['pop'][top].copy(), I['info'][top], True
                # diversity guard: champions of two islands must not be too similar; the weaker island
                # throws out everything resembling the stronger champion and refills with new shapes
                sims = [0.0] * nisl
                hold_need = max(3, d['patience'] // 2)
                if nisl > 1 and d['maxsim'] > 0:
                    sig = [shape_sig(I['bestg'], d) if I.get('bestg') is not None else None for I in isl]
                    for a in range(nisl):
                        for b in range(a + 1, nisl):
                            if sig[a] is None or sig[b] is None: continue
                            sv = similarity(sig[a], sig[b], d)
                            sims[a], sims[b] = max(sims[a], sv), max(sims[b], sv)
                            # a fresh migrant may briefly be a clone, that is the point of migration; only a
                            # similarity that persists for a while counts as convergence
                            hold[a, b] = hold.get((a, b), 0) + 1 if sv > d['maxsim'] else 0
                            if hold[a, b] < hold_need or isl[a].get('clear_wait', 0) or isl[b].get('clear_wait', 0): continue
                            hold[a, b] = 0
                            w = b if isl[b]['best'] <= isl[a]['best'] else a
                            win = a if w == b else b
                            W = isl[w]
                            close = [q for q in range(len(W['pop'])) if similarity(shape_sig(W['pop'][q], d), sig[win], d) > d['maxsim']]
                            newg = fresh(len(close), 0.3)
                            if newg is None: break
                            ev = self._map(pool, newg, 'Island %d: new shapes' % (w + 1)) if newg else []
                            for q, gnew, rr in zip(close, newg, ev):
                                W['pop'][q], W['fit'][q], W['info'][q] = gnew, rr[0], rr[1]
                            top = int(np.argmax(W['fit']))
                            W['best'], W['bestg'], W['besti'] = W['fit'][top], W['pop'][top].copy(), W['info'][top]
                            W['wires'] = wires_json(W['bestg'], d) if W['besti'] else None
                            W['order'] = np.argsort(-W['fit'])
                            W['stall'], W['boost'] = 0, 1.0
                            W['clear_wait'] = d['migrate'] or 8  # no migrants for a while, let it find its own way
                            W['clears'] = W.get('clears', 0) + 1
                            sig[w] = shape_sig(W['bestg'], d)
                            self.log('Islands %d and %d are %.0f %% similar (limit %.0f %%): island %d dropped %d similar shapes and searches elsewhere'
                                     % (a + 1, b + 1, sv, d['maxsim'], w + 1, len(close)))
                for I in isl:
                    if I.get('clear_wait'): I['clear_wait'] -= 1
                gstall = 0 if improved else gstall + 1
                hist.append(best_f)
                self.best = (best_g, best_f, best_i)
                valid = int(sum(np.sum(I['fit'] > -999) for I in isl))
                self._report({'phase': 'ga', 'gen': gen + 1, 'gens': d['gens'], 'fit': best_f, 'info': best_i,
                              'hist': hist, 'wires': wires_json(best_g, d) if best_g is not None else None,
                              '_pat': self._patfn(best_g) if improved or gen == 0 else None,
                              'elapsed': time.time() - t0, 'valid': valid,
                              'stall': gstall,
                              'isl': [{'fit': float(I['best']), 'info': I.get('besti'), 'wires': I.get('wires'),
                                       'stall': I['stall'], 'restarts': I['restarts'], 'sim': sims[j], 'clears': I.get('clears', 0),
                                       'champ': I.get('bestg') is not None and best_g is not None and np.array_equal(I['bestg'], best_g)}
                                      for j, I in enumerate(isl)] if nisl > 1 else None})
                if best_i and 'err' in best_i:
                    self.log('Generation %d/%d: fitness %.2f, pattern error %.2f dB RMS, peak %.2f %s, SWR %.2f, valid %d/%d, %.1f evaluations/s%s'
                             % (gen + 1, d['gens'], best_f, best_i['err'], best_i['floor'], POL_LABEL[d['pol']], best_i['swr'],
                                valid, d['pop'], self.evals / max(1e-6, time.time() - t0),
                                ('; islands ' + ' | '.join('%.2f' % I['best'] if I['best'] > -999 else '-' for I in isl)) if nisl > 1 else ''))
                elif best_i:
                    self.log('Generation %d/%d: fitness %.2f, min %.2f / mean %.2f %s, SWR %.2f%s, valid %d/%d, %.1f evaluations/s%s'
                             % (gen + 1, d['gens'], best_f, best_i['floor'], best_i['mean'], POL_LABEL[d['pol']],
                                best_i['swr'], ', suppression outside sector %.1f dB' % best_i['sup'] if 'sup' in best_i else '',
                                valid, d['pop'], self.evals / max(1e-6, time.time() - t0),
                                ('; islands ' + ' | '.join('%.2f' % I['best'] if I['best'] > -999 else '-' for I in isl)) if nisl > 1 else ''))
                else:
                    self.log('Generation %d/%d: no valid shape yet' % (gen + 1, d['gens']))
                if self._stopped(): break
                if gen == d['gens'] - 1: break
                if d['stop_after'] and gstall >= d['stop_after']:
                    self.log('Early stop: %d generations without improving the best shape' % gstall)
                    break
                # migration: the best of every island travel to the next one (ring), replacing its worst
                if nisl > 1 and d['migrate'] and (gen + 1) % d['migrate'] == 0:
                    moves = [(I['pop'][I['order'][:2]].copy(), I['fit'][I['order'][:2]].copy(), [I['info'][q] for q in I['order'][:2]]) for I in isl]
                    for j, I in enumerate(isl):
                        if I.get('clear_wait'): continue  # freshly cleared island: keep it unpolluted for a while
                        g2, f2, i2 = moves[(j - 1) % nisl]
                        worst = I['order'][-2:]
                        I['pop'][worst], I['fit'][worst] = g2, f2
                        for q, inf in zip(worst, i2): I['info'][q] = inf
                        I['order'] = np.argsort(-I['fit'])
                    self.log('Migration: the best shapes moved to the neighbouring islands')
                base_sigma = 0.15 * (1 - gen / d['gens']) + 0.01
                batch, plan = [], []
                for j, I in enumerate(isl):
                    nnew = 0
                    if d['patience'] and I['stall'] >= d['patience']:
                        nnew = size // 2
                        I['stall'], I['boost'], I['restarts'] = 0, 1.0, I['restarts'] + 1
                        self.log('Island %d: %d generations without improvement, renewing %d shapes and boosting mutation (restart #%d)'
                                 % (j + 1, d['patience'], nnew, I['restarts']))
                    sigma = max(base_sigma, 0.12 * I['boost'])
                    I['boost'] *= 0.75
                    newg = fresh(nnew, 0.3) if nnew else []
                    if newg is None: break
                    fit, pop = I['fit'], I['pop']
                    kids = []
                    while len(kids) < size - elite - len(newg):
                        for attempt in range(8):  # retry kids with impossible geometry (cheap, no NEC)
                            p1 = pop[max(rng.integers(0, len(pop), 3), key=lambda i: fit[i])]
                            p2 = pop[max(rng.integers(0, len(pop), 3), key=lambda i: fit[i])]
                            if rng.random() < 0.8:  # BLX-0.3 crossover
                                a, b = np.minimum(p1, p2), np.maximum(p1, p2)
                                c = a - 0.3 * (b - a) + rng.random(a.size) * 1.6 * (b - a)
                            else:
                                c = p1.copy()
                            m = rng.random(c.size) < 0.25
                            kid = project(np.clip(c + m * rng.normal(0, sigma * (0.7 ** attempt), c.size) * (hi - lo), lo, hi), d)
                            try:
                                if geometry_ok(kid, d): break
                            except Exception:
                                pass
                        kids.append(kid)
                    plan.append((j, len(batch), len(newg) + len(kids)))
                    batch += list(newg) + kids
                if self._stopped(): break
                r = self._map(pool, batch, 'Generation %d' % (gen + 1))
                for j, start, n in plan:
                    I = isl[j]
                    keep = I['order'][:elite]
                    I['pop'] = np.vstack([I['pop'][keep], np.array(batch[start:start + n])])
                    I['fit'] = np.concatenate([I['fit'][keep], [x[0] for x in r[start:start + n]]])
                    I['info'] = [I['info'][q] for q in keep] + [x[1] for x in r[start:start + n]]
            restarts = sum(I['restarts'] for I in isl)
            if restarts: self.log('The stagnation watchdog restarted islands %d times in total' % restarts)
            # polish: parallel stochastic hill climb around the winner
            if self._stopped(): self.log('Stopped by the user, keeping the best shape so far.')
            if d['polish'] and best_g is not None and not self._stopped():
                sigma, rounds = 0.02, max(1, min(d['polish'], d['pop'] * d['gens'] // 4) // max(1, d['jobs'] * 2))
                self.log('Polishing the winner with a local search: %d rounds of %d candidates' % (rounds, d['jobs'] * 2))
                for it in range(rounds):
                    if self._stopped(): break
                    cand = [project(np.clip(best_g + rng.normal(0, sigma, best_g.size) * (hi - lo), lo, hi), d)
                            for _ in range(d['jobs'] * 2)]
                    r = self._map(pool, cand, 'Polishing %d/%d' % (it + 1, rounds))
                    j = int(np.argmax([x[0] for x in r]))
                    if r[j][0] > best_f:
                        self.log('Polishing %d/%d: improved %.3f -> %.3f' % (it + 1, rounds, best_f, r[j][0]))
                        best_f, best_g, best_i = r[j][0], cand[j], r[j][1]
                    else:
                        sigma = max(0.002, sigma * 0.85)
                    self.best = (best_g, best_f, best_i)
                    self._report({'phase': 'polish', 'gen': d['gens'], 'gens': d['gens'], 'fit': best_f,
                                   'info': best_i, 'hist': hist + [best_f], 'wires': wires_json(best_g, d),
                                   '_pat': self._patfn(best_g),
                                  'elapsed': time.time() - t0, 'msg': 'Polishing the winner (%d/%d)' % (it + 1, rounds)})
        except _Stopped:
            self.log('Stopped by the user mid-generation, keeping the best shape so far.')
        finally:
            self._kill_pool()
        self.log('Done: %d NEC evaluations in %.0f s' % (self.evals, time.time() - t0))
        return self.best

# ---------------------------------------------------------------- outputs

def wires_json(g, d):
    return [{'arm': str(k), 'feed': i == 0, 'p': [float(x) for x in p], 'q': [float(x) for x in q]}
            for (k, i, p, q) in all_wires(g, d)]


def arm_chains(g, d):
    """[(label, copies, [points], free_start)] - each distinct arm to bend.  In free modes the
    polyline starts in the middle of the feed gap, so the first bend (at the terminal) is the
    solder joint and stays sharp."""
    arms, _ = arm_paths(g, d)
    par = [('Passive element %d' % (k - 9), 1, list(p), False) for (k, _, p) in arms if k >= 10]
    if ground(d):
        n = 1 if d['mode'] == 'mono' else 4
        return [('Arm' if n == 1 else 'Arm (4x identical)', n, arms[0][2], False)] + par
    out = []
    if d['mode'] == 'asym':
        for j, lab in ((0, 'A'), (1, 'B')):
            pts = next(p for (k, _, p) in arms if k == j)
            out.append(('Arm %s' % lab, 1, [np.zeros(3)] + list(pts), True))
        return out + par
    for j, lab in [(0, 'A')] + ([(1, 'B')] if d['mode'] == 'self' else []):
        pts = next(p for (k, _, p) in arms if k == j)
        centre = np.array([0.0, 0.0, pts[0][2]])
        out.append(('Arm %s (2x identical)' % lab, 2, [centre] + list(pts), True))
    return out + par


def develop(steps, rc_mm, free_start):
    """add the real (developed) wire with bend radius rc: straight part before every bend,
    position of the bend start mark measured from the start of the wire, arc length."""
    n = len(steps)
    tan = [0.0] * n
    for i, s in enumerate(steps):
        if s['bend'] and s['bend'] >= 0.5 and not (free_start and i == 0):
            tan[i] = rc_mm * math.tan(math.radians(s['bend']) / 2)
    pos = 0.0
    for i, s in enumerate(steps):
        if free_start and i == 0:  # half of the feed gap: not part of the wire
            s.update(straight=0.0, mark=0.0, arc=0.0, sharp=True); continue
        straight = s['len'] - tan[i] - (tan[i - 1] if i else 0.0)
        pos += straight
        arc = rc_mm * math.radians(s['bend']) if tan[i] else 0.0
        s.update(straight=straight, mark=pos, arc=arc, sharp=False)
        pos += arc
    return pos


def bend_plan(pts, ref_axis=(0, 0, 1)):
    """Turn a polyline into wire-bender steps.

    Returns [{'len': mm, 'bend': deg, 'rot': deg}] where for step i the wire runs 'len' mm,
    then gets bent by 'bend' degrees (0 = straight on, 90 = right angle) in a plane rotated by
    'rot' degrees from the previous bending plane: positive = clockwise when looking along the
    wire in the direction you are bending towards the free end.  The first bending plane is
    measured from the reference plane spanned by the first wire and the antenna axis.
    """
    pts = [np.asarray(p, float) for p in pts]
    dirs = [(pts[i + 1] - pts[i]) / np.linalg.norm(pts[i + 1] - pts[i]) for i in range(len(pts) - 1)]
    ref = np.asarray(ref_axis, float)
    n_prev = np.cross(dirs[0], ref)
    if np.linalg.norm(n_prev) < 1e-6:
        n_prev = np.cross(dirs[0], [1.0, 0, 0])
    n_prev /= np.linalg.norm(n_prev)
    steps = []
    for i in range(len(dirs)):
        L = float(np.linalg.norm(pts[i + 1] - pts[i])) * 1000
        if i == len(dirs) - 1:
            steps.append({'len': L, 'bend': None, 'rot': None}); break
        a, b = dirs[i], dirs[i + 1]
        bend = math.degrees(math.acos(np.clip(a @ b, -1, 1)))
        n = np.cross(a, b)
        if np.linalg.norm(n) < 1e-9:
            steps.append({'len': L, 'bend': 0.0, 'rot': 0.0}); continue
        n /= np.linalg.norm(n)
        rot = math.degrees(math.atan2(np.cross(n_prev, n) @ a, n_prev @ n))
        steps.append({'len': L, 'bend': bend, 'rot': rot})
        n_prev = n
    return steps


def rebuild_from_plan(steps, d0, ref_axis=(0, 0, 1)):
    """inverse of bend_plan (used for self-test): returns polyline from origin."""
    def rot_about(v, axis, ang):
        axis = axis / np.linalg.norm(axis)
        return v * math.cos(ang) + np.cross(axis, v) * math.sin(ang) + axis * (axis @ v) * (1 - math.cos(ang))
    dvec = np.asarray(d0, float) / np.linalg.norm(d0)
    n = np.cross(dvec, ref_axis)
    if np.linalg.norm(n) < 1e-6: n = np.cross(dvec, [1.0, 0, 0])
    n /= np.linalg.norm(n)
    p, pts = np.zeros(3), [np.zeros(3)]
    for s in steps:
        p = p + dvec * s['len'] / 1000; pts.append(p)
        if s['bend'] is None: break
        n = rot_about(n, dvec, math.radians(s['rot']))
        dvec = rot_about(dvec, n, math.radians(s['bend']))
    return pts


def plan_text(d, g):
    lines = []
    rc = d['rc'] * 1000
    for label, copies, pts, free in arm_chains(g, d):
        steps = bend_plan(pts)
        total = develop(steps, rc, free)
        lines.append('%s: cut %.1f mm of Ø %.2f mm wire (developed length%s)'
                     % (label, total, d['wire_d'] * 1000, ', add a margin for soldering' if not label.startswith('Passive') else ''))
        if label.startswith('Passive'):
            lines.append('  the wire starts at x %.1f, y %.1f, z %.1f mm (feed centre = 0, z axis = antenna axis) and connects to nothing'
                         % tuple(pts[0] * 1000))
            lines.append('  the wire axis passes through the vertices [mm]: ' + ' -> '.join('(%.0f, %.0f, %.0f)' % tuple(q * 1000) for q in pts))
        for i, s in enumerate(steps, 1):
            ref = 'from the reference plane' if i == 1 else 'from the previous bend plane'
            if s['sharp']:
                lines.append('  %d. the wire starts at the feed terminal (%.1f mm from the centre). Angle it %.0f° away from the feed axis, '
                             'with the bend plane rotated %+.0f° from the reference plane (soldered joint, sharp).'
                             % (i, s['len'], s['bend'], s['rot']))
            elif s['bend'] is None:
                lines.append('  %d. straight %.1f mm, then the end (wire ends at %.1f mm)' % (i, s['straight'], s['mark']))
            elif s['bend'] < 0.5:
                lines.append('  %d. straight %.1f mm, no bend' % (i, s['straight']))
            else:
                lines.append('  %d. straight %.1f mm, mark at %.1f mm = start of bend: bend %.0f° (inner radius %.1f mm), '
                             'rotate the plane %+.0f° %s'
                             % (i, s['straight'], s['mark'], s['bend'], d['bend_r'] * 1000, s['rot'], ref))
        lines.append('')
    lines += [
        'How to read: "straight" is the straight piece of wire between bends; marks are measured from the start of the wire along its axis',
        '(on the straight wire before bending). Each bend starts at its mark. Bend 0° = straight, 90° = right angle.',
        'Plane rotation: positive = clockwise when looking along the wire towards its free end.',
        'The bend radius (%.1f mm inner, %.1f mm at the wire axis) is included in the electrical model too.' % (d['bend_r'] * 1000, rc),
    ]
    if d['mode'] == 'asym':
        lines += ['The reference plane for the first step is the vertical plane containing the feed (the antenna axis is vertical, the feed lies on the x axis).',
                  'Arms A and B differ: A leaves the terminal at +x, B the terminal at -x.']
    elif d['mode'] in ('dipole', 'self'):
        lines += ['The reference plane for the first step is the vertical plane containing the feed (the antenna axis is vertical).',
                  'The second arm of the pair is bent exactly the same, only rotated 180° about the antenna axis.']
        if d['mode'] == 'self':
            lines.append('Dipole B is rotated 90° and its feed point sits %.1f mm higher than A. Both dipoles are connected in parallel.'
                         % (d['gap'] * 1000))
    elif d['mode'] == 'mono':
        lines += ['The wire starts at the feed point on the ground plane; the first section is vertical (perpendicular to the plane).',
                  'Reference plane for the first bend: the vertical plane through the stub and the x axis on the plane (pick the direction).']
    else:
        lines += ['Each arm starts vertically from the plane, %.1f mm from the axis; the arms are rotated in 90° steps.' % (d['r0'] * 1000),
                  'Reference plane for the first bend: the vertical plane through the stub and the outward direction from the axis.']
    if d['npar']:
        lines += ['Passive elements are separate bent wires with no connection. Hold them on non-conductive supports (PLA, acrylic, wood)',
                  'exactly at the given position; their position relative to the arms matters as much for directivity as their shape.']
    return '\n'.join(lines)


def summary(g, d):
    freqs = list(np.linspace(d['f0'] - d['bw'] / 2, d['f0'] + d['bw'] / 2, 5))
    res = simulate(g, d, freqs)
    return [{'f': r['f'], 'floor': r['floor'], 'mean': r['mean'], 'max': r['max'],
             'sup': (r['max'] - r['outmax']) if 'outmax' in r else None, 'err': r.get('err'),
             'zr': r['z'][0].real, 'zi': r['z'][0].imag, 'swr': max(swr(z, d['z0']) for z in r['z'])} for r in res]


def pat_sample(pat, pol, th, ph):
    """gain from a full pattern grid at (theta, phi) [deg], nearest grid point; below a ground plane -99."""
    G, st = pat[pol] if pat.get(pol) is not None else pat['total'], pat['step']
    th = th % 360; ph = ph % 360
    if th > 180: th, ph = 360 - th, (ph + 180) % 360
    i = int(round(th / st))
    if i >= len(G): return -99.0
    return float(G[i][int(round(ph / st)) % len(G[0])])


def plane_cut(pat, pol, plane, step=5):
    """[(angle, gain)] around one of the three main planes (pattern.PLANES conventions)."""
    out = []
    for a in range(0, 361, step):
        r = math.radians(a)
        v = {'xy': (math.cos(r), math.sin(r), 0.0), 'xz': (math.cos(r), 0.0, math.sin(r)), 'yz': (0.0, math.cos(r), math.sin(r))}[plane]
        th = math.degrees(math.acos(max(-1, min(1, v[2])))); ph = math.degrees(math.atan2(v[1], v[0]))
        out.append((a, pat_sample(pat, pol, th, ph)))
    return out


def pattern_cut(g, d, phi, pat=None):
    pat = pat or full_pattern(g, d)
    thmax = 90 if ground(d) else 180
    return [(th, pat_sample(pat, d['pol'], th, phi)) for th in range(0, thmax + 1, 5)]


def sector_cuts(g, d, pat=None):
    """azimuth cut at the middle elevation of the sector and the vertical cut through its centre."""
    pat = pat or full_pattern(g, d)
    elm = (d['el1'] + d['el2']) / 2
    if elm > 80: elm = max(d['el1'], min(60.0, d['el2']))  # a cut right at the zenith says nothing
    elm = round(elm / 5) * 5
    azc = (d['az1'] + d['az2']) / 2
    az = [(ph, pat_sample(pat, d['pol'], 90 - elm, ph)) for ph in range(-180, 181, 5)]
    # angle in the vertical plane through the sector centre: 0 = horizon in front, 90 = zenith,
    # 180 = horizon behind, 270 = nadir
    el = []
    for e in (range(0, 181, 5) if ground(d) else range(-90, 271, 5)):
        th, ph = (90 - e, azc) if e <= 90 else (e - 90, azc + 180)
        el.append((e, pat_sample(pat, d['pol'], th, ph)))
    return {'elm': elm, 'azc': azc, 'az': az, 'el': el}


def write_nec(g, d, path):
    wires = all_wires(g, d)
    fmax = d['f0'] + d['bw']
    with open(path, 'w') as f:
        f.write('CM Anten-O-Matic: %s, mode %s, f0 %.3f MHz, pol %s, wire %.2f mm\nCE\n'
                % (d['name'], d['mode'], d['f0'], d['pol'], d['wire_d'] * 1000))
        feeds = []
        for t, (k, i, p, q) in enumerate(wires, start=1):
            s = nsegs(np.linalg.norm(q - p), d, fmax, i == 0)
            f.write('GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % (t, s, *p, *q, d['radius']))
            if i == 0: feeds.append((t, (s + 1) // 2))
        f.write('GE 1\nGN 1\n' if ground(d) else 'GE 0\n')
        f.write('EK 0\n')
        if d['sigma']: f.write('LD 5 0 0 0 %g\n' % d['sigma'])
        f.write('FR 0 21 0 0 %.4f %.4f\n' % (d['f0'] - d['bw'], d['bw'] / 10))
        if d['mode'] == 'quad':
            sign = 1 if d['pol'] == 'lhcp' else -1
            for k, (t, s) in enumerate(feeds):
                ph = sign * math.pi / 2 * k
                f.write('EX 0 %d 1 0 %.5f %.5f\n' % (t, math.cos(ph), math.sin(ph)))
        elif d['mode'] == 'mono':
            f.write('EX 0 %d 1 0 1 0\n' % feeds[0][0])
        else:
            if d['mode'] == 'self':
                f.write('TL %d %d %d %d %g 0.00001\n' % (*feeds[0], *feeds[1], d['z0']))
            f.write('EX 0 %d %d 0 1 0\n' % feeds[0])
        f.write('RP 0 37 73 1000 0 0 5 5\nEN\n')


POL_LABEL = {'rhcp': 'dBic RHCP', 'lhcp': 'dBic LHCP', 'total': 'dBi'}
MODE_LABEL = {'mono': 'single arm over a ground plane', 'dipole': 'bent dipole (one arm + copy)',
              'self': 'two dipoles, self-phasing (QFH style)', 'quad': '4 arms + phasing network (ST5)',
              'asym': 'asymmetric dipole (two different arms)'}


def render(g, d, path, cuts=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sec = d.get('cov') == 'sector'
    tgt = d.get('cov') == 'target'
    fig = plt.figure(figsize=(20 if tgt else 16 if sec else 12, 5.4))
    ax = fig.add_subplot(141 if tgt else 131 if sec else 121, projection='3d')
    cols = ['#2a78d6', '#e8710a', '#1aa36f', '#d63a4a', '#8a5cd6', '#c79a00']
    lw = max(1.5, min(6, d['wire_d'] * 1000 / (lam(d) * 1000) * 300))
    allp = []
    for (k, i, p, q) in all_wires(g, d):
        kk = int(str(k).lstrip('f')) if not str(k).startswith('f') else -1
        ax.plot(*zip(p * 1000, q * 1000), color='k' if i == 0 else ('#7a828c' if kk >= 10 else cols[kk % 4]), lw=lw, solid_capstyle='round')
        allp += [p * 1000, q * 1000]
    allp = np.array(allp)
    R = max(np.abs(allp[:, :2]).max(), 1) * 1.1
    if ground(d):
        gx, gy = np.meshgrid([-R, R], [-R, R])
        ax.plot_surface(gx, gy, 0 * gx, alpha=0.15, color='grey')
    ax.set_xlim(-R, R); ax.set_ylim(-R, R)
    z0, z1 = allp[:, 2].min(), allp[:, 2].max()
    if ground(d): z0 = 0
    ax.set_zlim(z0, max(z1, z0 + 1))
    ax.set_box_aspect((1, 1, max(0.25, (z1 - z0) / (2 * R))))
    ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]'); ax.set_zlabel('z [mm]')
    ax.set_title('%s  %.1f MHz, wire Ø %.1f mm' % (d['name'], d['f0'], d['wire_d'] * 1000))
    if tgt:
        pat = cuts
        th_t, ph_t = np.array([90.0]), np.arange(0, 361, 5.0)
        for i, (plane, title) in enumerate((('xy', 'horizontal plane XY (0° = +x)'), ('xz', 'vertical plane XZ (0° = +x, 90° = +z)'),
                                            ('yz', 'vertical plane YZ (0° = +y, 90° = +z)'))):
            axp = fig.add_subplot(1, 4, 2 + i, projection='polar')
            cut = plane_cut(pat, d['pol'], plane)
            gmax = max(pat_sample(pat, d['pol'], th, ph) for th in range(0, 181, 5) for ph in range(0, 360, 5))
            ang = np.radians([c[0] for c in cut]); gv = np.array([c[1] for c in cut]) - gmax
            fl = d['target']['floor']
            tv = PAT._cut(d['target'][plane], [c[0] for c in cut])
            axp.set_theta_zero_location('E')
            axp.plot(ang, np.clip(tv, fl, 0), color='#2a78d6', lw=1.5, ls='--', label='target')
            axp.plot(ang, np.clip(gv, fl, 0), color='#e8710a', lw=2, label='antenna')
            axp.set_rlim(fl, 0); axp.set_rticks([fl, fl * 2 / 3, fl / 3, 0])
            axp.set_title('%s [dB rel.]' % title, fontsize=9)
            if i == 0: hl = axp.get_legend_handles_labels()
        fig.legend(*hl, loc='lower center', ncol=2, fontsize=9, frameon=False)
        fig.suptitle('%s, %s; peak %.1f %s' % (MODE_LABEL[d['mode']], coverage_text(d), gmax, POL_LABEL[d['pol']]), fontsize=10, y=0.99)
        fig.tight_layout(rect=(0, 0.05, 1, 0.93)); fig.savefig(path, dpi=100); plt.close(fig)
        return
    if sec:
        sc = cuts
        unit = POL_LABEL[d['pol']]
        for i, (key, title) in enumerate((('az', 'azimuth at elevation %.0f°' % sc['elm']), ('el', 'vertical cut at azimuth %.0f°' % sc['azc']))):
            axp = fig.add_subplot(1, 3, 2 + i, projection='polar')
            ang = np.radians([c[0] for c in sc[key]]); gv = np.array([c[1] for c in sc[key]])
            top = np.ceil(np.nanmax(gv) / 5) * 5; floor = top - 30
            if key == 'az':
                axp.set_theta_zero_location('E')
                a1, a2 = np.radians(d['az1']), np.radians(d['az2'])
                if not sector_az_limited(d): a1, a2 = 0, 2 * np.pi
            else:
                axp.set_theta_zero_location('E')  # 0 = horizon in the sector direction, 90 = zenith
                a1, a2 = np.radians(d['el1']), np.radians(d['el2'])
            axp.fill_between(np.linspace(a1, a2, 60), floor, top, color='#2a78d6', alpha=0.1)
            axp.plot(ang, np.clip(gv, floor, None), color='#e8710a', lw=2)
            axp.set_rlim(floor, top); axp.set_rticks([top - 30, top - 20, top - 10, top])
            if key == 'el' and ground(d): axp.set_thetamin(0); axp.set_thetamax(180)
            axp.set_title('%s [%s]' % (title, unit), fontsize=10)
        fig.suptitle('%s, %s' % (MODE_LABEL[d['mode']], coverage_text(d)), fontsize=10, y=0.99)
        fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(path, dpi=105); plt.close(fig)
        return
    ax2 = fig.add_subplot(122)
    cuts = cuts or {phi: pattern_cut(g, d, phi) for phi in (0, 45, 90)}
    lo = 0
    for phi, cut in cuts.items():
        ax2.plot([c[0] for c in cut], [c[1] for c in cut], label='azimuth %d°' % phi)
        lo = min(lo, min(c[1] for c in cut))
    ax2.axvspan(0, d['cone'], color='#2a78d6', alpha=0.07, label='required coverage')
    ax2.set_xlabel('angle from the antenna axis [°]  (0 = z axis)')
    ax2.set_ylabel('gain [%s]' % POL_LABEL[d['pol']])
    ax2.set_ylim(max(-25, lo - 3), None)
    ax2.grid(alpha=0.3); ax2.legend(loc='lower left')
    ax2.set_title(MODE_LABEL[d['mode']])
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def export_step(g, d, path):
    import stepout
    arms, _ = arm_paths(g, d)
    bodies = []
    for k, stub, pts in arms:
        pieces = fillet([p * 1000 for p in pts], d['rc'] * 1000)
        bodies.append(('arm_%s' % k, pieces, d['radius'] * 1000))
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in d['name'])[:40] or 'antenna'
    stepout.write_step(path, bodies, safe)


def save_outputs(g, d, out_root):
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in d['name'])[:40] or 'antenna'
    out = os.path.join(out_root, '%s_%s_%gMHz_%s' % (stamp, safe, d['f0'], d['mode']))
    os.makedirs(out, exist_ok=True)
    summ = summary(g, d)
    pat = full_pattern(g, d)
    with open(os.path.join(out, 'pattern.json'), 'w') as f: json.dump(pat, f)
    cuts = pat if d['cov'] == 'target' else sector_cuts(g, d, pat) if d['cov'] == 'sector' else {phi: pattern_cut(g, d, phi, pat) for phi in (0, 45, 90)}
    render(g, d, os.path.join(out, 'render.png'), cuts)
    write_nec(g, d, os.path.join(out, 'antenna.nec'))
    if d['step']: export_step(g, d, os.path.join(out, 'antenna.step'))
    plan = plan_text(d, g)
    unit = POL_LABEL[d['pol']]
    rep = ['Anten-O-Matic: %s' % d['name'], '=' * 40,
           'Type: %s' % MODE_LABEL[d['mode']],
           'Frequency %.3f MHz ± %.3f MHz, polarisation %s, coverage %s'
           % (d['f0'], d['bw'] / 2, d['pol'].upper(), coverage_text(d)),
           'Coordinates: antenna axis = z (vertical), azimuth 0° = +x axis, 90° = +y axis, elevation 90° = zenith.',
           'Wire: Ø %.2f mm, %s; λ = %.1f mm' % (d['wire_d'] * 1000, d['material'], lam(d) * 1000), '',
           'Freq. [MHz]     min in region   mean   max [%s]   Z [Ω]          SWR (%g Ω)%s'
           % (unit, d['z0'], '   suppression outside [dB]' if d['cov'] == 'sector' else '   pattern error [dB RMS]' if d['cov'] == 'target' else '')]
    for s in summ:
        rep.append('  %10.3f      %6.2f      %6.2f   %6.2f        %7.1f %+7.1fj   %5.2f%s'
                   % (s['f'], s['floor'], s['mean'], s['max'], s['zr'], s['zi'], s['swr'],
                      '        %6.1f' % s['sup'] if s.get('sup') is not None else '        %6.2f' % s['err'] if s.get('err') is not None else ''))
    rep += ['', 'BENDING INSTRUCTIONS', '-' * 40, plan, '', 'Notes:'] + ['- ' + w for w in warnings(d)]
    rep += ['', 'Coordinates of all wires [mm] (x, y, z  ->  x, y, z):']
    for (k, i, p, q) in all_wires(g, d):
        rep.append('  %-3s %s  (%7.1f %7.1f %7.1f) -> (%7.1f %7.1f %7.1f)'
                   % (k, 'feed   ' if i == 0 else '       ', *(p * 1000), *(q * 1000)))
    text = '\n'.join(rep) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as f: f.write(text)
    steps = {}
    for lab, _, pts, free in arm_chains(g, d):
        steps[lab] = bend_plan(pts); develop(steps[lab], d['rc'] * 1000, free)
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as f:
        json.dump({'params': {k: v for k, v in d.items()}, 'genome': g.tolist(), 'summary': summ,
                   'bend_plan': steps, 'wires': wires_json(g, d)}, f, indent=1, ensure_ascii=False)
    return out, text, summ, steps
