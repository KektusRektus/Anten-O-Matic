# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""Antenna CAD: a document of elements (rods, arcs/fillets, tori, helices, plates, blocks) placed in
3D, joined at their end points, with a material per element. Turns the document into NEC wires,
simulates it, writes the usual outputs (render, report, NEC, STEP, pattern.json).

Units in the document: mm and degrees. Element pose: pos [x, y, z] + rot [rx, ry, rz] (Euler XYZ,
same convention as three.js), world = pos + Rx Ry Rz . local.
"""
import json, math, os, time
import numpy as np

import core
import pattern as PAT
import pcb as PCB
from core import C, swr, _chk

TYPES = ('rod', 'fillet', 'torus', 'helix', 'meander', 'custom', 'plate', 'block', 'pcb')
PATH_TYPES = ('rod', 'fillet', 'torus', 'helix', 'meander', 'custom')

DEFAULT_MATERIALS = [
    {'name': 'Copper', 'kind': 'conductor', 'sigma': 5.8e7, 'color': '#c77b3a', 'builtin': True},
    {'name': 'Brass', 'kind': 'conductor', 'sigma': 1.5e7, 'color': '#c9a640', 'builtin': True},
    {'name': 'Aluminium', 'kind': 'conductor', 'sigma': 3.5e7, 'color': '#a9b2bb', 'builtin': True},
    {'name': 'Steel', 'kind': 'conductor', 'sigma': 1.4e6, 'color': '#7c858e', 'builtin': True},
    {'name': 'Silver', 'kind': 'conductor', 'sigma': 6.3e7, 'color': '#d7dade', 'builtin': True},
    {'name': 'Perfect conductor', 'kind': 'conductor', 'sigma': 0, 'color': '#e8710a', 'builtin': True},
    {'name': 'FR4', 'kind': 'dielectric', 'eps': 4.4, 'tand': 0.02, 'color': '#3f8f4f', 'builtin': True},
    {'name': 'Rogers RO4003C', 'kind': 'dielectric', 'eps': 3.38, 'tand': 0.0027, 'color': '#d9c9a0', 'builtin': True},
    {'name': 'PTFE', 'kind': 'dielectric', 'eps': 2.1, 'tand': 0.0002, 'color': '#eeeeee', 'builtin': True},
    {'name': 'PLA (3D print)', 'kind': 'dielectric', 'eps': 2.7, 'tand': 0.01, 'color': '#5b9ef0', 'builtin': True},
]


# ------------------------------------------------------------------ material library (data/materials.json)

def clean_material(m):
    out = {'name': str(m.get('name', '')).strip()[:60], 'kind': m.get('kind')}
    if not out['name']: raise ValueError('material needs a name')
    if out['kind'] not in ('conductor', 'dielectric'): raise ValueError('material kind must be conductor or dielectric')
    if out['kind'] == 'conductor':
        out['sigma'] = float(m.get('sigma', 0))
        if out['sigma'] < 0: raise ValueError('conductivity must be >= 0 (0 = perfect conductor)')
    else:
        out['eps'] = float(m.get('eps', 1)); out['tand'] = float(m.get('tand', 0))
        if out['eps'] < 1 or not 0 <= out['tand'] < 1: raise ValueError('need eps_r >= 1 and 0 <= tan d < 1')
    c = str(m.get('color', '#888888'))
    out['color'] = c if len(c) == 7 and c[0] == '#' else '#888888'
    return out


def load_materials(data_dir):
    lib = [dict(m) for m in DEFAULT_MATERIALS]
    try:
        with open(os.path.join(data_dir, 'materials.json'), encoding='utf-8') as f:
            user = json.load(f)
        names = {m['name'] for m in lib}
        for m in user:
            try:
                m = clean_material(m)
                if m['name'] not in names: lib.append(m); names.add(m['name'])
            except Exception:
                pass
    except Exception:
        pass
    return lib


def save_materials(data_dir, mats):
    builtin = {m['name'] for m in DEFAULT_MATERIALS}
    user = [clean_material(m) for m in mats if m.get('name') not in builtin]
    os.makedirs(data_dir, exist_ok=True)
    with open(os.path.join(data_dir, 'materials.json'), 'w', encoding='utf-8') as f:
        json.dump(user, f, indent=1, ensure_ascii=False)
    return load_materials(data_dir)


# ------------------------------------------------------------------ geometry

def rotation(rot):
    rx, ry, rz = np.radians(np.asarray(rot, float))
    Rx = np.array([[1, 0, 0], [0, math.cos(rx), -math.sin(rx)], [0, math.sin(rx), math.cos(rx)]])
    Ry = np.array([[math.cos(ry), 0, math.sin(ry)], [0, 1, 0], [-math.sin(ry), 0, math.cos(ry)]])
    Rz = np.array([[math.cos(rz), -math.sin(rz), 0], [math.sin(rz), math.cos(rz), 0], [0, 0, 1]])
    return Rx @ Ry @ Rz


def _num(p, k, lo=None):
    v = float(p[k])
    if not math.isfinite(v) or (lo is not None and v < lo): raise ValueError('%s must be >= %g' % (k, lo))
    return v


def local_path(e):
    """polyline (N x 3, mm, local frame) of a wire element. A rod may carry explicit 'pts' (bent
    by the user or by the evolution); otherwise it is a straight rod from 0 to L along +x."""
    t, p = e['type'], e['params']
    if t == 'rod':
        if p.get('pts'): return np.array(p['pts'], float)
        return np.array([[0, 0, 0], [_num(p, 'L', 0.01), 0, 0]], float)
    if t == 'fillet':  # arc starting at the origin heading +x, turning towards +y
        R, a = _num(p, 'R', 0.01), math.radians(_num(p, 'angle', 1))
        n = max(2, int(math.ceil(math.degrees(a) / 10)))
        s = np.linspace(0, a, n + 1)
        return np.stack([R * np.sin(s), R * (1 - np.cos(s)), 0 * s], -1)
    if t == 'torus':   # loop (or open arc) around the origin in the local xy plane, starts at (R, 0, 0)
        R, a = _num(p, 'R', 0.01), math.radians(min(360.0, _num(p, 'arc', 1)))
        n = max(3, int(math.ceil(math.degrees(a) / 10)))
        s = np.linspace(0, a, n + 1)
        P = np.stack([R * np.cos(s), R * np.sin(s), 0 * s], -1)
        return P
    if t == 'helix':   # along +z from (R, 0, 0)
        R, pitch, turns = _num(p, 'R', 0.01), float(p['pitch']), _num(p, 'turns', 0.05)
        hand = -1 if p.get('hand', 'R') == 'L' else 1
        n = max(4, int(math.ceil(turns * 24)))
        s = np.linspace(0, 2 * math.pi * turns, n + 1)
        return np.stack([R * np.cos(s), hand * R * np.sin(s), pitch * s / (2 * math.pi)], -1)
    if t == 'meander':  # rectangular meander in the local xy plane: n half-periods over L, amplitude A
        L, A, n = _num(p, 'L', 0.01), float(p['A']), max(1, int(round(_num(p, 'n', 1))))
        s_ = L / n
        P = [[0, 0, 0]]
        for k in range(n):
            y0, y1 = (0, A) if k % 2 == 0 else (A, 0)
            P += [[k * s_, y1, 0], [(k + 1) * s_, y1, 0]]
        P.append([L, 0 if n % 2 == 0 else A, 0])
        P = [q for i, q in enumerate(P) if i == 0 or q != P[i - 1]]
        return np.array(P, float)
    if t == 'custom':
        return local_paths(e)[0]
    raise ValueError('not a wire element: %s' % t)


def local_paths(e):
    """all polylines of a wire element (a custom sketch can have several, a closed shape is one
    polyline whose last point equals the first)."""
    if e['type'] != 'custom': return [local_path(e)]
    out = []
    for P in e['params'].get('paths') or []:
        A = np.array([[float(q[0]), float(q[1]), float(q[2]) if len(q) > 2 else 0.0] for q in P], float)
        if len(A) >= 2 and np.isfinite(A).all(): out.append(A)
    if not out: raise ValueError('%s: the sketch has no lines' % e.get('name', 'custom element'))
    return out


def pixel_mask(e):
    """(nx, ny, mask[ny][nx] bool) of a plate cut into pixels by the evolution, or None."""
    p = e['params']
    if not p.get('mask') or not p.get('pix'): return None
    nx, ny = int(p['pix'][0]), int(p['pix'][1])
    m = np.array(p['mask'], bool).reshape(ny, nx)
    return nx, ny, m


def local_points(e):
    """named snap points in the local frame."""
    t, p = e['type'], e['params']
    if t == 'custom':
        Ps = local_paths(e)
        pts = {'start': Ps[0][0], 'end': Ps[-1][-1], 'center': np.zeros(3)}
        for k, v in (p.get('pts') or {}).items():
            pts[str(k)] = np.array([float(v[0]), float(v[1]), float(v[2]) if len(v) > 2 else 0.0])
        return pts
    if t == 'pcb': return PCB.local_points(e)
    if t in PATH_TYPES:
        P = local_path(e)
        pts = {'start': P[0], 'end': P[-1]}
        if t == 'torus': pts['center'] = np.zeros(3)
        elif t == 'helix': pts['center'] = np.array([0, 0, P[-1][2] / 2])
        else: pts['center'] = P[len(P) // 2] if len(P) > 2 else (P[0] + P[-1]) / 2
        return pts
    W, H = _num(p, 'W', 0.01), _num(p, 'H', 0.01)
    T = _num(p, 'T', 0) if t == 'block' else 0.0
    pts = {'center': np.array([0, 0, T / 2])}
    for nm, sx, sy in (('c1', -1, -1), ('c2', 1, -1), ('c3', 1, 1), ('c4', -1, 1), ('e1', 0, -1), ('e2', 1, 0), ('e3', 0, 1), ('e4', -1, 0)):
        pts[nm] = np.array([sx * W / 2, sy * H / 2, 0.0])
        if t == 'block': pts[nm + 't'] = np.array([sx * W / 2, sy * H / 2, T])
    return pts


def world(e, P):
    return np.asarray(e['pos'], float) + np.asarray(P, float) @ rotation(e['rot']).T


def world_points(e):
    return {k: world(e, v) for k, v in local_points(e).items()}


def euler_xyz(R):
    """Euler XYZ angles [deg] (three.js convention, R = Rx Ry Rz) of a rotation matrix."""
    m13 = max(-1.0, min(1.0, R[0, 2]))
    y = math.asin(m13)
    if abs(m13) < 0.9999999:
        x = math.atan2(-R[1, 2], R[2, 2]); z = math.atan2(-R[0, 1], R[0, 0])
    else:
        x = math.atan2(R[2, 1], R[1, 1]); z = 0.0
    return [math.degrees(x), math.degrees(y), math.degrees(z)]


def axis_rot(axis, ang):
    """rotation matrix about a unit axis by ang [rad] (Rodrigues)."""
    a = np.asarray(axis, float); a = a / max(np.linalg.norm(a), 1e-12)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


def element_dir(e):
    """the direction an angle dimension measures: start -> end of a wire (its normal when the wire
    is closed), the normal of a plate or block."""
    R = rotation(e['rot'])
    if e['type'] in PATH_TYPES:
        Ps = local_paths(e)
        v = Ps[-1][-1] - Ps[0][0]
        if np.linalg.norm(v) > 1e-6: return R @ (v / np.linalg.norm(v))
    return R[:, 2]


def ref_point(doc, ref):
    if ref.get('el') == 'origin': return np.zeros(3)
    for e in doc['elements']:
        if e['id'] == ref.get('el'):
            l = local_points(e).get(ref.get('pt'))
            return None if l is None else world(e, l)
    return None


def _by_id(doc, i):
    for e in doc['elements']:
        if e['id'] == i: return e
    return None


def rotate_about(e, pivot, Rd):
    """turn element e rigidly by Rd about the world point pivot."""
    R = Rd @ rotation(e['rot'])
    e['rot'] = euler_xyz(R)
    e['pos'] = list(np.asarray(pivot) + Rd @ (np.asarray(e['pos'], float) - np.asarray(pivot)))


def angle_now(doc, d):
    ea, eb = _by_id(doc, d['a']['el']), _by_id(doc, d['b']['el'])
    if ea is None or eb is None: return None
    u, v = element_dir(ea), element_dir(eb)
    return math.degrees(math.acos(max(-1.0, min(1.0, float(u @ v)))))


def shared_point(doc, ea, eb):
    for j in doc['joins']:
        for x, y in ((j['a'], j['b']), (j['b'], j['a'])):
            if x['el'] == ea['id'] and y['el'] == eb['id']: return ref_point(doc, y)
    return None


def solve(doc, locked=(), iters=400):
    """same solver as the CAD window: elements are rigid; joints (coincident points), driving
    distance / x / y / z dimensions are met by translating elements, angle dimensions by turning
    one of the two elements (about their shared joint, or its own start). Fixed and locked elements
    never move. Returns the worst remaining error [mm]."""
    free = lambda e: e is not None and not e.get('fixed') and e['id'] not in locked
    cons = [{'a': j['a'], 'b': j['b'], 'axis': 'join', 'value': 0} for j in doc['joins']]
    cons += [d for d in doc.get('dims') or [] if d.get('drive', True)]
    worst = 0.0
    for _ in range(iters):
        worst = 0.0
        for c in cons:
            if c['axis'] == 'angle':
                ea, eb = _by_id(doc, c['a']['el']), _by_id(doc, c['b']['el'])
                if ea is None or eb is None: continue
                cur = angle_now(doc, c)
                err = float(c['value']) - cur
                worst = max(worst, abs(err) * math.pi / 180 * 10)    # 10 mm lever: comparable to the mm errors
                if abs(err) < 1e-5: continue
                mov, oth, sgn = (eb, ea, 1) if free(eb) else (ea, eb, 1) if free(ea) else (None, None, 0)
                if mov is None: continue
                u, v = element_dir(oth), element_dir(mov)
                n = np.cross(u, v)
                if np.linalg.norm(n) < 1e-9:
                    n = np.cross(u, [0, 0, 1]) if abs(u[2]) < 0.9 else np.cross(u, [1, 0, 0])
                piv = shared_point(doc, mov, oth)
                if piv is None:
                    piv = world(mov, local_points(mov).get('center', np.zeros(3)))
                rotate_about(mov, piv, axis_rot(n, math.radians(err)))
                continue
            A, B = ref_point(doc, c['a']), ref_point(doc, c['b'])
            if A is None or B is None: continue
            v = B - A
            if c['axis'] == 'join': corr = v
            elif c['axis'] == 'd':
                L = float(np.linalg.norm(v))
                if L < 1e-9: continue
                corr = v * (L - float(c['value'])) / L
            else:
                k = 'xyz'.index(c['axis']); corr = np.zeros(3); corr[k] = v[k] - float(c['value'])
            err = float(np.linalg.norm(corr)); worst = max(worst, err)
            if err < 1e-6: continue
            ea = None if c['a']['el'] == 'origin' else _by_id(doc, c['a']['el'])
            eb = None if c['b']['el'] == 'origin' else _by_id(doc, c['b']['el'])
            wa, wb = (1 if free(ea) else 0), (1 if free(eb) else 0)
            if not wa and not wb: continue
            s_ = 1.0 / (wa + wb)
            if wa: ea['pos'] = list(np.asarray(ea['pos'], float) + corr * s_)
            if wb: eb['pos'] = list(np.asarray(eb['pos'], float) - corr * s_)
        if worst < 1e-5: break
    return worst


def clean_doc(doc, materials):
    if not isinstance(doc, dict): raise ValueError('CAD document missing')
    out = {'name': str(doc.get('name') or 'cad')[:40], 'f0': float(doc.get('f0', 0)), 'bw': float(doc.get('bw', 0)),
           'z0': float(doc.get('z0', 50)), 'ground': doc.get('ground', 'none'), 'segs': int(doc.get('segs', 25)),
           'swr': float(doc.get('swr', 2.0)), 'step': bool(doc.get('step', True)), 'svg': bool(doc.get('svg', True)), 'elements': [], 'dims': list(doc.get('dims') or []),
           'joins': list(doc.get('joins') or []), 'engine': doc.get('engine', 'nec'), 'ems_mesh': doc.get('ems_mesh') or 'normal', 'ems_air': float(doc.get('ems_air') or 20), 'ems_fine': float(doc.get('ems_fine') or 40), 'ems_margin': float(doc.get('ems_margin') or 0.25)}
    if out['engine'] not in ('nec', 'openems'): raise ValueError('engine must be nec or openems')
    if out['ems_mesh'] not in ('draft', 'coarse', 'normal', 'fine', 'vfine', 'custom'): out['ems_mesh'] = 'normal'
    out['ems_margin'] = min(0.5, max(0.1, out['ems_margin']))
    if out['f0'] <= 0: raise ValueError('frequency must be positive')
    if out['bw'] < 0: raise ValueError('bandwidth must be >= 0')
    if out['ground'] not in ('none', 'perfect'): raise ValueError('ground must be none or perfect')
    mats = {m['name']: m for m in materials}
    ids = set()
    for e in doc.get('elements') or []:
        t = e.get('type')
        if t not in TYPES: raise ValueError('unknown element type %s' % t)
        q = {'id': str(e.get('id')), 'type': t, 'name': str(e.get('name') or t)[:40],
             'pos': [float(v) for v in (e.get('pos') or [0, 0, 0])], 'rot': [float(v) for v in (e.get('rot') or [0, 0, 0])],
             'params': dict(e.get('params') or {}), 'mat': e.get('mat') or 'Copper',
             'd': float(e.get('d', 2.0)), 'feed': bool(e.get('feed')), 'feed_at': float(e.get('feed_at', 0.5)),
             'feed_path': int(e.get('feed_path') or 0), 'fixed': bool(e.get('fixed')), 'evolve': dict(e.get('evolve') or {})}
        if q['id'] in ids: raise ValueError('duplicate element id %s' % q['id'])
        ids.add(q['id'])
        if q['mat'] not in mats: raise ValueError('unknown material %s on %s' % (q['mat'], q['name']))
        if q['d'] <= 0: raise ValueError('wire diameter of %s must be positive' % q['name'])
        if not 0 <= q['feed_at'] <= 1: raise ValueError('feed position must be 0..1')
        if t == 'pcb':
            PCB.check(q, mats); q['feed'] = bool(q['params'].get('feed'))
        local_points(q)  # validates the parameters
        out['elements'].append(q)
    if not out['elements']: raise ValueError('the design has no elements')
    out['_mats'] = {k: mats[k] for k in {e['mat'] for e in out['elements']} | {e['params'].get('sub') for e in out['elements'] if e['type'] == 'pcb' and e['params'].get('sub')}}
    out['_lam'] = C / out['f0'] * 1000
    return out


def is_conductor(doc, e):
    return doc['_mats'][e['mat']]['kind'] == 'conductor'


def _split_at(P, t):
    """insert the point at arc-length fraction t; returns (polyline, index of that point)."""
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    L = seg.sum()
    s = t * L
    acc = 0.0
    for i, l in enumerate(seg):
        if acc + l >= s - 1e-9:
            u = (s - acc) / l if l > 0 else 0
            q = P[i] + u * (P[i + 1] - P[i])
            if u * l < 2e-3: return P, i               # [mm]: snap onto an existing corner
            if (1 - u) * l < 2e-3: return P, i + 1
            return np.vstack([P[:i + 1], q, P[i + 1:]]), i + 1
        acc += l
    return P, len(P) - 1


def _decimate(P, m):
    """drop vertices closer than m [mm] to the last kept one (ends stay): NEC-2 needs every wire
    a few diameters long, or a short wire's midpoint ends up inside its neighbour on a tight bend."""
    P = np.asarray(P, float)
    if len(P) < 3 or m <= 0: return P
    keep = [P[0]]
    for q in P[1:-1]:
        if np.linalg.norm(q - keep[-1]) >= m: keep.append(q)
    if len(keep) > 1 and np.linalg.norm(P[-1] - keep[-1]) < m: keep.pop()   # merge a short tail
    keep.append(P[-1])
    return np.array(keep)


def nec_wires(doc):
    """[(element index, kind, A, B, radius_m)] in metres; kind 'feed' marks the excited wire,
    'wire' a normal piece, 'grid' a plate mesh line. Plate meshes come with their own segment
    count (so mesh nodes coincide with segment ends and NEC connects the crossing lines)."""
    out, warn = [], []
    lam = doc['_lam']
    for k, e in enumerate(doc['elements']):
        if not is_conductor(doc, e):
            warn.append('%s is a dielectric: NEC-2 ignores it (the openEMS engine models dielectrics)' % e['name'])
            continue
        r = e['d'] / 2000
        if e['type'] in PATH_TYPES:
            paths = [_decimate(world(e, Q), 2 * e['d']) for Q in local_paths(e)]
            fp = min(max(0, e.get('feed_path', 0)), len(paths) - 1)
            for pi, P in enumerate(paths):
                if not (e['feed'] and pi == fp):
                    for j in range(len(P) - 1):
                        if np.linalg.norm(P[j + 1] - P[j]) > 1e-6: out.append((k, 'wire', P[j] / 1000, P[j + 1] / 1000, r))
                    continue
                P, i = _split_at(P, e['feed_at'])
                L_tot = float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum())
                g = min(lam / doc['segs'], L_tot / 7)   # feed segment as long as its neighbours [mm]
                c = P[i]
                # never past the neighbouring corners, and no sliver left before them (a feed gap stays a gap)
                lim = 4 * e['d'] / 2
                if 0 < i < len(P) - 1:
                    for nb in (2 * float(np.linalg.norm(P[i - 1] - c)), 2 * float(np.linalg.norm(P[i + 1] - c))):
                        if g > nb or nb - g < 2 * lim: g = nb
                else:
                    nb = float(np.linalg.norm(P[1 if i == 0 else -2] - c))
                    if g > nb or nb - g < lim: g = nb
                if i == 0 or i == len(P) - 1:          # feed at an end (e.g. monopole on the ground)
                    j = 1 if i == 0 else len(P) - 2
                    v = P[j] - c; v = v / np.linalg.norm(v)
                    a, b = c, c + v * g
                    rest = [np.vstack([b, P[1:]])] if i == 0 else [np.vstack([P[:-1], b])]
                    if i != 0: a, b = b, c
                    out.append((k, 'feed', a / 1000, b / 1000, r))
                else:
                    v1 = P[i - 1] - c; v1 /= np.linalg.norm(v1)
                    v2 = P[i + 1] - c; v2 /= np.linalg.norm(v2)
                    a, b = c + v1 * g / 2, c + v2 * g / 2
                    out.append((k, 'feed', a / 1000, b / 1000, r))
                    rest = [np.vstack([P[:i], a]), np.vstack([b, P[i + 1:]])]
                for Q in rest:
                    for j in range(len(Q) - 1):
                        if np.linalg.norm(Q[j + 1] - Q[j]) > 1e-6: out.append((k, 'wire', Q[j] / 1000, Q[j + 1] / 1000, r))
        elif e['type'] == 'plate':
            p = e['params']
            W, H = float(p['W']), float(p['H'])
            fmax = doc['f0'] + doc['bw'] / 2
            mesh = float(p.get('mesh') or 0) or C / fmax * 1000 / 20
            nx = max(2, 2 * int(math.ceil(W / mesh / 2))); ny = max(2, 2 * int(math.ceil(H / mesh / 2)))
            pm = pixel_mask(e)
            if pm is not None:            # NEC cells inside the pixels: a multiple of the pixel grid
                nx = pm[0] * max(1, int(math.ceil(nx / pm[0]))); ny = pm[1] * max(1, int(math.ceil(ny / pm[1])))
            xs, ys = plate_lines(doc, k, W, H, nx, ny)
            rg = (float(p['mesh_r']) if p.get('mesh_r') else max(W / nx, H / ny) / (2 * math.pi)) / 1000   # equal-area rule
            on = cell_mask(e, xs, ys)
            for i, x in enumerate(xs):          # one wire per mesh cell edge: NEC joins them at the nodes
                for j, y in enumerate(ys):
                    if j < len(ys) - 1 and (on is None or on(i - 1, j) or on(i, j)):
                        out.append((k, ('grid', 1), world(e, [x, y, 0]) / 1000, world(e, [x, ys[j + 1], 0]) / 1000, rg))
                    if i < len(xs) - 1 and (on is None or on(i, j - 1) or on(i, j)):
                        out.append((k, ('grid', 1), world(e, [x, y, 0]) / 1000, world(e, [xs[i + 1], y, 0]) / 1000, rg))
        elif e['type'] == 'pcb':
            out += PCB.nec_wires(doc, k, lam, world, warn)
        else:
            warn.append('%s: conducting blocks are not supported by NEC-2, use a plate' % e['name'])
    out = heal_wires(out, warn, [e['name'] for e in doc['elements']])
    if not any(w[1] == 'feed' for w in out): raise ValueError('no feed: tick "Feed" on one wire element')
    if sum(1 for w in out if w[1] == 'feed') > 1: warn.append('more than one feed: all are driven in phase with 1 V')
    if doc['ground'] == 'perfect':
        zmin = min(min(w[2][2], w[3][2]) for w in out)
        if zmin < -1e-6: raise ValueError('with a ground plane, nothing may go below z = 0')
    return out, warn


def _grid(w): return isinstance(w[1], tuple)


def heal_wires(wires, warn, names):
    """make NEC-2 happy where wires touch or cross away from a shared end (its 'segment midpoint intersects
    wire' error): a wire end lying on another wire splits that wire there, an end that nearly meets
    another end snaps onto it, and two wires crossing in their middles are split and joined at the
    crossing. Moves are at most the two radii. Wires of the same plate grid are left alone."""
    W = [list(w) for w in wires]
    told, budget = set(), 4 * len(W) + 50

    def note(i, j, what):
        a, b = names[W[i][0]], names[W[j][0]]
        key = (min(a, b), max(a, b), what)
        if key not in told:
            told.add(key)
            warn.append(('%s touches %s away from a joint: joined there' if what == 't' else
                         '%s crosses %s: joined at the crossing' if what == 'x' else
                         '%s nearly meets %s: ends snapped together') % (a, b) if a != b else
                        '%s touches itself away from a joint: joined there' % a)

    def move(p, q):   # every wire end at p moves to q
        for w in W:
            for f in (2, 3):
                if np.linalg.norm(w[f] - p) < 1e-9: w[f] = q.copy()

    def split(j, F):
        w = W[j]
        W[j] = [w[0], w[1], w[2], F.copy(), w[4]]
        W.insert(j + 1, [w[0], w[1], F.copy(), w[3], w[4]])

    def skip(i, j):
        return i == j or (_grid(W[i]) and _grid(W[j]) and W[i][0] == W[j][0])

    def step():
        A = np.array([w[2] for w in W]); B = np.array([w[3] for w in W]); R = np.array([w[4] for w in W])
        D = B - A; L2 = np.maximum((D * D).sum(1), 1e-30); L = np.sqrt(L2)
        # 1) a wire end on (or touching) another wire; the closest case first, so a true contact (a tangent)
        #    is joined before its near neighbours are considered
        best = None
        for i, w in enumerate(W):
            for p, o in ((w[2], w[3]), (w[3], w[2])):
                at = (np.linalg.norm(A - p, axis=1) < 1e-9) | (np.linalg.norm(B - p, axis=1) < 1e-9)
                Li = float(L[at].min()) if at.any() else float(np.linalg.norm(o - p))   # shortest wire hanging on this end
                t = np.clip(((p - A) * D).sum(1) / L2, 0, 1)
                F = A + D * t[:, None]; d = np.linalg.norm(p - F, axis=1)
                tol = R + w[4]
                for j in np.argsort(d):
                    if d[j] >= tol[j]: break
                    if skip(i, j): continue
                    eps = max(w[4], R[j], 1e-7)
                    if np.linalg.norm(p - A[j]) < 1e-9 or np.linalg.norm(p - B[j]) < 1e-9: continue   # already joined
                    if np.linalg.norm(o - A[j]) < 1e-9 or np.linalg.norm(o - B[j]) < 1e-9: continue   # its other end is joined there: a corner
                    if d[j] > 0.25 * Li: continue                                                         # would squash this wire
                    if t[j] * L[j] > eps and (1 - t[j]) * L[j] > eps:
                        if W[j][1] == 'feed': continue
                        if best is None or d[j] < best[0]: best = (d[j], 't', i, j, p.copy(), F[j].copy())
                        break
                    q = A[j] if t[j] * L[j] <= eps else B[j]
                    if np.linalg.norm(p - q) > 1e-9:
                        dq = float(np.linalg.norm(p - q))
                        if best is None or dq < best[0]: best = (dq, 's', i, j, p.copy(), q.copy())
                        break
        if best is not None:
            _, kind, i, j, p, q = best
            note(i, j, kind); move(p, q)
            if kind == 't': split(j, q)
            return True
        # 2) two wires crossing in their middles
        for i, w in enumerate(W):
            u = D[i]; w0 = A[i] - A
            a = L2[i]; b = (D * u).sum(1); c = L2; e = (w0 * u).sum(1)   # e = u·w0
            den = a * c - b * b
            ok = den > 1e-12 * a * c
            # closest points of the two infinite lines: s on wire i, t on wire j
            sc = np.where(ok, (b * ((D * w0).sum(1)) - c * e) / np.where(ok, den, 1), 0)
            tc = np.where(ok, (a * ((D * w0).sum(1)) - b * e) / np.where(ok, den, 1), 0)
            Pi = A[i] + np.outer(sc, u); Pj = A + D * tc[:, None]
            dist = np.linalg.norm(Pi - Pj, axis=1)
            for j in np.where(ok & (dist < R + w[4]) & (sc * L[i] > R + w[4]) & ((1 - sc) * L[i] > R + w[4]) & (tc * L > R + w[4]) & ((1 - tc) * L > R + w[4]))[0]:
                if j <= i or skip(i, j) or W[i][1] == 'feed' or W[j][1] == 'feed': continue
                q = (Pi[j] + Pj[j]) / 2
                note(i, j, 'x'); split(max(i, j), q); split(min(i, j), q); return True
        return False

    while budget > 0 and step(): budget -= 1
    W = [w for w in W if np.linalg.norm(w[3] - w[2]) > 1e-9]
    # wires that now lie exactly on each other (a tangent contact pulled together): keep one
    seen, keep = {}, []
    for w in W:
        key = tuple(sorted((tuple(np.round(w[2] * 1e7)), tuple(np.round(w[3] * 1e7)))))
        if key in seen:
            o = seen[key]
            if w[1] == 'feed' and o[1] != 'feed': keep[keep.index(o)] = w; seen[key] = w
            if o[4] < w[4]: o[4] = w[4]
            continue
        seen[key] = w; keep.append(w)
    W = keep
    # what cannot be joined: two wires running along each other closer than their radii
    A = np.array([w[2] for w in W]); B = np.array([w[3] for w in W]); R = np.array([w[4] for w in W])
    D = B - A; L = np.linalg.norm(D, axis=1); U = D / L[:, None]
    for i, w in enumerate(W):
        sin = np.linalg.norm(np.cross(U, U[i]), axis=1)
        off = np.linalg.norm(np.cross(A - A[i], U[i]), axis=1)
        s0 = ((A - A[i]) * U[i]).sum(1); s1 = ((B - A[i]) * U[i]).sum(1)
        ov = np.minimum(np.maximum(s0, s1), L[i]) - np.maximum(np.minimum(s0, s1), 0)
        for j in np.where((sin < 0.05) & (off < R + w[4]) & (ov > np.maximum(R, w[4]) * 2))[0]:
            if j <= i or skip(i, j): continue
            shared = sum(float(np.linalg.norm(p - q)) < 1e-9 for p in (A[i], B[i]) for q in (A[j], B[j]))
            if shared == 1: continue      # joined at one end and leaving at a shallow angle (a tangent): NEC-2 accepts that
            a, b = names[W[i][0]], names[W[j][0]]
            at = ' (near %.1f, %.1f, %.1f mm)' % tuple((W[i][2] + W[i][3]) / 2 * 1000)
            raise ValueError(('%s and %s run along each other closer than their wire radius%s. '
                              'NEC-2 cannot model that: move them apart by a few wire diameters, or make them one wire.' % (a, b, at)) if a != b else
                             ('%s runs back over itself closer than its wire radius%s. NEC-2 cannot model that: open the shape up.' % (a, at)))
    return [tuple(w) for w in W]


def cell_mask(e, xs, ys):
    """for a pixel plate: on(i, j) -> is NEC mesh cell (xs[i]..xs[i+1], ys[j]..ys[j+1]) metal."""
    pm = pixel_mask(e)
    if pm is None: return None
    nx, ny, m = pm
    W, H = float(e['params']['W']), float(e['params']['H'])

    def on(i, j):
        if i < 0 or j < 0 or i >= len(xs) - 1 or j >= len(ys) - 1: return False
        cx, cy = (xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2
        a = min(nx - 1, max(0, int((cx + W / 2) / W * nx))); b = min(ny - 1, max(0, int((cy + H / 2) / H * ny)))
        return bool(m[b, a])
    return on


def pixel_rects(e):
    """metal rectangles (x0, y0, x1, y1) in the plate frame: whole plate or its pixels, merged along x."""
    W, H = float(e['params']['W']), float(e['params']['H'])
    pm = pixel_mask(e)
    if pm is None: return [(-W / 2, -H / 2, W / 2, H / 2)]
    nx, ny, m = pm
    out = []
    for b in range(ny):
        a = 0
        while a < nx:
            if not m[b, a]: a += 1; continue
            a1 = a
            while a1 + 1 < nx and m[b, a1 + 1]: a1 += 1
            out.append((-W / 2 + a * W / nx, -H / 2 + b * H / ny, -W / 2 + (a1 + 1) * W / nx, -H / 2 + (b + 1) * H / ny))
            a = a1 + 1
    return out


def plate_lines(doc, k, W, H, nx, ny):
    """grid line coordinates of plate k: a uniform mesh, with the nearest lines moved onto every wire
    end point that touches the plate, so the wire connects to a mesh node (e.g. a patch feed probe)."""
    e = doc['elements'][k]
    xs, ys = list(np.linspace(-W / 2, W / 2, nx + 1)), list(np.linspace(-H / 2, H / 2, ny + 1))
    R = rotation(e['rot']); o = np.asarray(e['pos'], float)
    tol = max(1e-3, 0.02 * min(W / nx, H / ny))
    for j, o_ in enumerate(doc['elements']):
        if j == k or o_['type'] not in PATH_TYPES or not is_conductor(doc, o_): continue
        for q in [P_[i_] for P_ in (world(o_, Q) for Q in local_paths(o_)) for i_ in (0, -1)]:
            l = (q - o) @ R                      # into the plate frame
            if abs(l[2]) > tol or abs(l[0]) > W / 2 + tol or abs(l[1]) > H / 2 + tol: continue
            for lines, v in ((xs, l[0]), (ys, l[1])):
                i = int(np.argmin([abs(a - v) for a in lines]))
                if abs(lines[i] - v) < 1e-6: continue
                if 0 < i < len(lines) - 1: lines[i] = float(v)
                else: lines.insert(1 if i == 0 else len(lines) - 1, float(v))
    return sorted(xs), sorted(ys)


def _segs(L, r, fmax, segs):
    want = math.ceil(L / (C / fmax / segs))
    cap = max(1, int(L / (4 * r)))
    return max(1, min(want, cap))


def _geom_msg(doc, wires, ex):
    """NEC-2's geometry errors in words: which elements, and what to do."""
    import re
    m = str(ex)
    if 'GEOMETRY' not in m and 'INTERSECT' not in m and 'SEGMENT' not in m: return m
    tags = [int(t) for t in re.findall(r'TAG ID #(\d+)', m)]
    names = [doc['elements'][wires[t - 1][0]]['name'] for t in tags if 0 < t <= len(wires)] + ['?', '?']
    if 'INTERSECT' in m:
        a, b = names[0], names[1]
        who = ('two parts of %s' % a) if a == b else ('%s and %s' % (a, b))
        return ('%s lie on top of each other (closer than their wire radius) and could not be joined automatically: '
                'usually two parallel wires overlapping, or a wire running inside a thicker one. Move them apart by a few '
                'wire diameters, or make them one wire. A bend or small arc tighter than about 3 wire diameters does the '
                'same: draw it larger or use a thinner wire. (NEC-2: %s)' % (who, m.strip()))
    return 'NEC-2 cannot use the geometry of %s: %s' % (names[0], m.strip())


def nec_run(doc, freqs, rp, pat_step=None):
    import necwrap as N
    wires, warn = nec_wires(doc)
    n = N.nec_create()
    try:
        feeds, tags = [], {}
        fmax = max(freqs)
        for tag, (k, kind, a, b, r) in enumerate(wires, start=1):
            L = float(np.linalg.norm(b - a))
            if kind == 'feed': s = 1
            elif isinstance(kind, tuple): s = kind[1]
            else: s = _segs(L, r, fmax, doc['segs'])
            try: _chk(N.nec_wire(n, tag, s, *a, *b, r, 1.0, 1.0))
            except RuntimeError as ex: raise RuntimeError(_geom_msg(doc, wires, ex)) from None
            if kind == 'feed': feeds.append(tag)
            tags.setdefault(k, []).append(tag)
        gnd = doc['ground'] == 'perfect'
        try:
            _chk(N.nec_geometry_complete(n, 1 if gnd else 0))
        except RuntimeError as ex:
            raise RuntimeError(_geom_msg(doc, wires, ex)) from None
        if gnd: _chk(N.nec_gn_card(n, 1, 0, 0, 0, 0, 0, 0, 0))
        _chk(N.nec_ek_card(n, 0))
        for k, tl in tags.items():
            sig = doc['_mats'][doc['elements'][k]['mat']].get('sigma', 0)
            if sig:
                for tag in tl: _chk(N.nec_ld_card(n, 5, tag, 0, 0, sig, 0, 0))
        nf = len(freqs)
        _chk(N.nec_fr_card(n, 0, nf, freqs[0], (freqs[-1] - freqs[0]) / (nf - 1) if nf > 1 else 0))
        for tag in feeds: _chk(N.nec_excitation_voltage(n, tag, 1, 1.0, 0.0))
        nth, nph, th0, ph0, dth, dph = rp
        _chk(N.nec_rp_card(n, 0, nth, nph, 0, 0, 0, 0, th0, ph0, dth, dph, 0, 0))
        res = []
        for fi, f in enumerate(freqs):
            z = complex(N.nec_impedance_real_indexed(n, fi, 0), N.nec_impedance_imag_indexed(n, fi, 0))
            r = {'f': f, 'z': z}
            for pol in ('total', 'rhcp', 'lhcp'):
                r[pol] = {s: N.gain(n, pol, s, fi) for s in ('min', 'mean', 'max')}
            if pat_step: r['pattern'] = PAT.from_nec(N, n, fi, pat_step, gnd)
            res.append(r)
        return res, warn
    finally:
        N.nec_delete(n)


def wires_json(doc):
    wires, _ = nec_wires(doc)
    return [{'arm': str(k), 'feed': kind == 'feed', 'p': [float(v) for v in a], 'q': [float(v) for v in b]}
            for k, kind, a, b, r in wires]


# ------------------------------------------------------------------ outputs

def write_nec(doc, path):
    wires, _ = nec_wires(doc)
    fmax = doc['f0'] + doc['bw'] / 2
    with open(path, 'w') as f:
        f.write('CM Anten-O-Matic CAD design %s, %.3f MHz\nCE\n' % (doc['name'], doc['f0']))
        feeds, tags = [], {}
        for tag, (k, kind, a, b, r) in enumerate(wires, start=1):
            L = float(np.linalg.norm(b - a))
            s = 1 if kind == 'feed' else kind[1] if isinstance(kind, tuple) else _segs(L, r, fmax, doc['segs'])
            f.write('GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n' % (tag, s, *a, *b, r))
            if kind == 'feed': feeds.append(tag)
            tags.setdefault(k, []).append(tag)
        f.write('GE 1\nGN 1\n' if doc['ground'] == 'perfect' else 'GE 0\n')
        f.write('EK 0\n')
        for k, tl in tags.items():
            sig = doc['_mats'][doc['elements'][k]['mat']].get('sigma', 0)
            if sig:
                for tag in tl: f.write('LD 5 %d 0 0 %g\n' % (tag, sig))
        nf = 21 if doc['bw'] else 1
        f.write('FR 0 %d 0 0 %.4f %.4f\n' % (nf, doc['f0'] - doc['bw'] / 2, doc['bw'] / 20 if doc['bw'] else 0))
        for tag in feeds: f.write('EX 0 %d 1 0 1 0\n' % tag)
        f.write('RP 0 37 73 1000 0 0 5 5\nEN\n')


def step_bodies(doc):
    """tube bodies (wire elements) and prisms (plates, blocks) for the STEP export, in mm."""
    tubes, prisms = [], []
    for e in doc['elements']:
        safe = ''.join(c if c.isalnum() else '_' for c in e['name']) or e['type']
        if e['type'] in PATH_TYPES:
            for pi, Q in enumerate(local_paths(e)):
                P = world(e, Q)
                r = e['d'] / 2
                rc = float(e.get('params', {}).get('bend_r') or 0) + 1.5 * r
                pieces = None
                while pieces is None and rc > 1e-3:
                    pieces = core.fillet(P, rc)
                    if pieces is None: rc /= 2
                if pieces: tubes.append((safe + ('_%d' % (pi + 1) if pi else ''), pieces, r))
        elif e['type'] == 'pcb':
            continue                       # planar layouts go out as SVG (layout_*.svg)
        else:
            p = e['params']
            T = float(p.get('T') or 0) if e['type'] == 'block' else float(p.get('t') or 0.035)
            rects = pixel_rects(e) if e['type'] == 'plate' else [(-float(p['W']) / 2, -float(p['H']) / 2, float(p['W']) / 2, float(p['H']) / 2)]
            for ri, (x0, y0, x1, y1) in enumerate(rects):
                q = dict(e, pos=list(world(e, [(x0 + x1) / 2, (y0 + y1) / 2, 0])))
                prisms.append((safe + ('_%d' % (ri + 1) if len(rects) > 1 else ''), q, x1 - x0, y1 - y0, T))
    return tubes, prisms


def export_step(doc, path):
    """wire elements as exact tubes; plates and blocks as boxes (placed with their own pose)."""
    import stepout
    tubes, prisms = step_bodies(doc)
    w = stepout._W()

    def solids():
        out = [stepout.tube_solid(w, pieces, r, name) for name, pieces, r in tubes]
        for name, e, W, H, T in prisms:
            out.append(box_solid(w, e, W, H, T, name))
        return out
    stepout._write(path, w, solids, doc['name'], 'Anten-O-Matic CAD antenna')


def box_solid(w, e, W, H, T, name):
    """a box with an arbitrary pose (8 vertices, 6 planar faces)."""
    import stepout
    loc = [np.array([sx * W / 2, sy * H / 2, z]) for z in (0.0, T) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    V = [world(e, q) for q in loc]
    vp = [w.add("VERTEX_POINT('',%s)" % w.pt(q)) for q in V]

    def edge(i, j):
        d = V[j] - V[i]
        ln = w.add("LINE('',%s,%s)" % (w.pt(V[i]), w.add("VECTOR('',%s,%s)" % (w.dir(d), w.num(np.linalg.norm(d))))))
        return w.add("EDGE_CURVE('',%s,%s,%s,.T.)" % (vp[i], vp[j], ln))
    E = {}
    for a, b in ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)):
        E[(a, b)] = edge(a, b)

    def oe(a, b):
        return w.add("ORIENTED_EDGE('',*,*,%s,%s)" % ((E[(a, b)], '.T.') if (a, b) in E else (E[(b, a)], '.F.')))
    R = rotation(e['rot'])
    faces = []
    for ring, nrm in (((0, 3, 2, 1), -R[:, 2]), ((4, 5, 6, 7), R[:, 2]), ((0, 1, 5, 4), -R[:, 1]), ((1, 2, 6, 5), R[:, 0]),
                      ((2, 3, 7, 6), R[:, 1]), ((3, 0, 4, 7), -R[:, 0])):
        loop = w.add("EDGE_LOOP('',(%s))" % ','.join(oe(ring[i], ring[(i + 1) % 4]) for i in range(4)))
        xd = V[ring[1]] - V[ring[0]]
        faces.append(w.add("ADVANCED_FACE('',(%s),%s,.T.)" % (w.add("FACE_OUTER_BOUND('',%s,.T.)" % loop),
                                                           w.add("PLANE('',%s)" % w.ax2(V[ring[0]], nrm, xd)))))
    shell = w.add("CLOSED_SHELL('',(%s))" % ','.join(faces))
    return w.add("MANIFOLD_SOLID_BREP('%s',%s)" % (stepout._ascii(name), shell))


def render(doc, path, pat, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 5.2))
    ax = fig.add_subplot(131, projection='3d')
    wires, _ = nec_wires(doc)
    allp = []
    for k, kind, a, b in [(w[0], w[1], w[2], w[3]) for w in wires]:
        col = 'k' if kind == 'feed' else doc['_mats'][doc['elements'][k]['mat']]['color']
        ax.plot(*zip(a * 1000, b * 1000), color=col, lw=1 if isinstance(kind, tuple) else 2.5)
        allp += [a * 1000, b * 1000]
    for e in doc['elements']:
        if e['type'] == 'block':
            c = world_points(e)
            ring = [c[k] for k in ('c1', 'c2', 'c3', 'c4', 'c1')]
            ax.plot(*zip(*ring), color=doc['_mats'][e['mat']]['color'], lw=1, ls='--')
            allp += ring
    allp = np.array(allp)
    c = (allp.max(0) + allp.min(0)) / 2; s = max((allp.max(0) - allp.min(0)).max() / 2, 1) * 1.1
    if doc['ground'] == 'perfect':
        gx, gy = np.meshgrid([c[0] - s, c[0] + s], [c[1] - s, c[1] + s])
        ax.plot_surface(gx, gy, 0 * gx, alpha=0.15, color='grey')
    ax.set_xlim(c[0] - s, c[0] + s); ax.set_ylim(c[1] - s, c[1] + s); ax.set_zlim(c[2] - s, c[2] + s)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]'); ax.set_zlabel('z [mm]')
    ax.set_title('%s %.1f MHz' % (doc['name'], doc['f0']))
    for i, (plane, lab) in enumerate((('xz', 'vertical plane xz'), ('xy', 'horizontal plane xy'))):
        cut = core.plane_cut(pat, 'total', plane)
        a = [q[0] for q in cut]; g = np.array([q[1] for q in cut])
        top = np.ceil(np.nanmax(g) / 5) * 5
        axp = fig.add_subplot(1, 3, 2 + i, projection='polar')
        axp.set_theta_zero_location('E')
        axp.plot(np.radians(a), np.clip(g, top - 30, None), color='#e8710a', lw=2)
        axp.set_rlim(top - 30, top); axp.set_rticks([top - 30, top - 20, top - 10, top])
        axp.set_title('%s [dBi], max %.1f' % (lab, np.nanmax(g)), fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=105); plt.close(fig)


def build_text(doc):
    lines = ['Elements (positions in mm, rotations in degrees, Euler XYZ):']
    for e in doc['elements']:
        p = e['params']
        if (e['type'] == 'rod' and p.get('pts')) or e['type'] == 'custom':
            Ps = local_paths(e)
            L = sum(float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum()) for P in Ps)
            geo = ('bent rod, %d points' % len(Ps[0]) if e['type'] == 'rod' else 'drawn element, %d wire runs' % len(Ps)) + ', wire length %.1f' % L
        elif e['type'] == 'pcb':
            geo = PCB.describe(e)
        elif e['type'] == 'plate' and pixel_mask(e) is not None:
            nx, ny, m = pixel_mask(e)
            geo = 'W %.2f, H %.2f, cut into %dx%d pixels, %d metal: %s' % (float(p['W']), float(p['H']), nx, ny, int(m.sum()),
                  ' / '.join(''.join('#' if v else '.' for v in row) for row in m[::-1]))
        else:
            geo = ', '.join('%s %s' % (k, ('%.2f' % v) if isinstance(v, (int, float)) else v) for k, v in p.items() if k != 'pts')
        wp = world_points(e)
        ends = ''
        if e['type'] in PATH_TYPES:
            ends = '; start (%.1f, %.1f, %.1f), end (%.1f, %.1f, %.1f)' % (*wp['start'], *wp['end'])
        if e['type'] == 'pcb':
            lines.append('  %-14s %-7s %-18s %s' % (e['name'], e['type'], e['mat'], geo)); continue
        lines.append('  %-14s %-7s %-18s %s, wire %.2f mm%s%s' % (e['name'], e['type'], e['mat'], geo, e['d'],
                                                                ends, ', FEED at %.0f %%' % (e['feed_at'] * 100) if e['feed'] else ''))
    if doc['joins']:
        lines.append('Joined points (solder / weld): %d' % len(doc['joins']))
    for d in doc.get('dims') or []:
        if d.get('axis') == 'angle': lines.append('Angle between %s and %s: %.1f deg' % (_nm(doc, d['a']), _nm(doc, d['b']), d['value']))
    return '\n'.join(lines)


def _nm(doc, ref):
    for e in doc['elements']:
        if e['id'] == ref.get('el'): return e['name']
    return 'origin' if ref.get('el') == 'origin' else '?'


def save(doc, out_root, log=print):
    freqs = list(np.linspace(doc['f0'] - doc['bw'] / 2, doc['f0'] + doc['bw'] / 2, 5)) if doc['bw'] else [doc['f0']]
    res, warn = nec_run(doc, freqs, PAT.rp_for(10, doc['ground'] == 'perfect'))
    for w_ in warn: log('Note: ' + w_)
    rows = [{'f': r['f'], 'zr': r['z'].real, 'zi': r['z'].imag, 'swr': swr(r['z'], doc['z0']), 'gmax': r['total']['max'],
             'rhcp': r['rhcp']['max'], 'lhcp': r['lhcp']['max']} for r in res]
    pat = nec_run(doc, [doc['f0']], PAT.rp_for(5, doc['ground'] == 'perfect'), pat_step=5)[0][0]['pattern']
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in doc['name'])[:40] or 'cad'
    out = os.path.join(out_root, '%s_%s_%gMHz_cad' % (stamp, safe, doc['f0']))
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, 'pattern.json'), 'w') as f: json.dump(pat, f)
    render(doc, os.path.join(out, 'render.png'), pat, rows)
    write_nec(doc, os.path.join(out, 'antenna.nec'))
    if doc.get('step', True): export_step(doc, os.path.join(out, 'antenna.step'))
    if doc.get('svg', True): PCB.write_svgs(doc, out)
    btxt = build_text(doc)
    rep = ['Anten-O-Matic CAD: %s, %.3f MHz' % (doc['name'], doc['f0']), '=' * 40,
           'λ = %.1f mm, feed impedance %.0f Ω, ground: %s' % (doc['_lam'], doc['z0'], doc['ground']), '',
           'MHz         Z [Ω]            SWR    max [dBi]   RHCP [dBic]   LHCP [dBic]']
    for r in rows:
        rep.append('%9.3f  %7.1f %+7.1fj   %5.2f   %6.2f      %6.2f        %6.2f' % (r['f'], r['zr'], r['zi'], r['swr'], r['gmax'], r['rhcp'], r['lhcp']))
    if warn: rep += [''] + ['Note: ' + w_ for w_ in warn]
    rep += ['', 'CONSTRUCTION', '-' * 40, btxt]
    text = '\n'.join(rep) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as f: f.write(text)
    keep = {k: v for k, v in doc.items() if not k.startswith('_')}
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as f:
        json.dump({'type': 'cad', 'doc': keep, 'materials': list(doc['_mats'].values()), 'summary': rows}, f, indent=1, ensure_ascii=False)
    summ = [{'f': r['f'], 'floor': r['gmax'], 'mean': None, 'max': r['gmax'], 'zr': r['zr'], 'zi': r['zi'], 'swr': r['swr']} for r in rows]
    return out, text, summ, btxt, wires_json(doc), warn
