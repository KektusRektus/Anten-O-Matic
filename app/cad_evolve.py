# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""Evolution seeded from a CAD design. Every element can get one kind of freedom:
  shape   - its own parameters (length, radius, pitch, turns, size...) within +- size %, its
            position within +- dev mm and its rotation within +- rot degrees;
  free    - the element is approximated by a chain of control points (stations) spread along the
            wire, every station moves within +- dev mm and the wire between them follows; any wire
            (rod, arc, loop, helix, meander, drawn element) can become a free-form bent wire;
  pixels  - a plate is cut into nx x ny pixels that the evolution switches on and off.
Driving dimensions can get a range too (the distance or angle the solver keeps varies within it),
so the evolution also changes spacings and angles between elements. Joined points stay joined:
rigid elements are moved by the same solver as the CAD window, free-form wires that hang on them
follow, and wire ends sitting on a plate stay on it. The fitness is the chosen goal minus a
mismatch-loss penalty for SWR above the target.
"""
import copy, math, multiprocessing as mp, os, time
import numpy as np

import cad
import core
import pattern as PAT
import pcb as PCB

GRID = 10   # pattern grid step [deg] during the evolution

# numeric parameters the shape mode may change, and whether they are integers
SHAPE_KEYS = {'rod': (('L', 0),), 'fillet': (('R', 0), ('angle', 0)), 'torus': (('R', 0), ('arc', 0)),
              'helix': (('R', 0), ('pitch', 0), ('turns', 0)), 'meander': (('L', 0), ('A', 0), ('n', 1)),
              'plate': (('W', 0), ('H', 0)), 'block': (('W', 0), ('H', 0), ('T', 0)), 'custom': (), 'pcb': (('h', 0),)}


def _ml(s):
    g = (s - 1) / (s + 1)
    return -10 * math.log10(max(1e-6, 1 - g * g))


def evolve_opts(e):
    """normalised freedom of an element (also reads the 3.0 format {bends, dev, planar})."""
    ev = dict(e.get('evolve') or {})
    mode = ev.get('mode')
    dev = float(ev.get('dev') or 0)
    if mode is None:                                        # 3.0 design
        if dev <= 0: mode = 'off'
        elif e['type'] == 'rod': mode = 'free'; ev['pts'] = int(ev.get('bends') or 0) + 2
        else: mode = 'shape'; ev.setdefault('size', 0); ev['_mm'] = dev
    out = {'mode': mode, 'dev': dev, 'size': float(ev.get('size') or 0), 'rot': float(ev.get('rot') or 0),
           'pts': max(2, int(ev.get('pts') or (int(ev.get('bends') or 0) + 2))), 'planar': bool(ev.get('planar')),
           'pix': [max(1, int(v)) for v in (ev.get('pix') or [8, 8])][:2], '_mm': ev.get('_mm')}
    if mode == 'pixels' and e['type'] not in ('plate', 'pcb'): out['mode'] = 'off'
    if mode == 'free' and e['type'] not in cad.PATH_TYPES + ('pcb',): out['mode'] = 'shape'
    return out


def _arclen(P):
    return np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])


def _with_stations(P, n):
    """polyline with n stations spread evenly along it: returns (polyline, station vertex indices).
    The original vertices stay, so a helix or meander keeps its shape until the stations move."""
    s = _arclen(P)
    if s[-1] <= 0: return P, [0, len(P) - 1]
    want = np.linspace(0, s[-1], n)
    pts, idx = [P[0]], [0]
    for i in range(1, len(P)):
        for t in want[(want > s[i - 1] + 1e-9) & (want < s[i] - 1e-9)]:
            u = (t - s[i - 1]) / (s[i] - s[i - 1]); pts.append(P[i - 1] + u * (P[i] - P[i - 1])); idx.append(len(pts) - 1)
        pts.append(P[i])
        if np.min(np.abs(want - s[i])) < 1e-9: idx.append(len(pts) - 1)
    return np.array(pts), sorted(set(idx) | {0, len(pts) - 1})


class Space:
    """maps a genome in [-1, 1]^n onto a modified copy of the document."""

    def __init__(self, doc):
        self.doc = doc
        els = doc['elements']
        self.opt = [evolve_opts(e) for e in els]
        self.slots = []                      # (kind, a, b, scale)
        self.frame = {k: cad.rotation(e['rot']) for k, e in enumerate(els)}
        free = [k for k, o in enumerate(self.opt) if o['mode'] == 'free' and o['dev'] > 0 and els[k]['type'] != 'pcb']
        # ---- free-form wires: world polylines with stations
        self.poly = {}                       # k -> [(world polyline, station indices)]
        for k in free:
            e = els[k]
            out = []
            for pi, Q in enumerate(cad.local_paths(e)):
                W = cad.world(e, Q)
                if e['type'] == 'custom':
                    nl = e['params'].get('nodes') or []
                    nodes = nl[pi] if pi < len(nl) else []
                    st = sorted({0, len(W) - 1} | {int(i) for i in nodes if 0 <= int(i) < len(W)})
                    out.append((W, st))
                else:
                    out.append(_with_stations(W, self.opt[k]['pts']))
            self.poly[k] = out
        # stations of one element at the same place are one station (closed loops, drawn shapes)
        parent = {}

        def find(x):
            while parent.setdefault(x, x) != x: x = parent[x]
            return x

        def union(a, b): parent[find(a)] = find(b)
        self.st_pos = {}
        for k, out in self.poly.items():
            seen = {}
            for pi, (W, st) in enumerate(out):
                for i in st:
                    key = (k, pi, i); find(key); self.st_pos[key] = W[i].copy()
                    h = tuple(np.round(W[i], 4))
                    if h in seen: union(key, seen[h])
                    else: seen[h] = key
        # joins: a free station joined to another free station shares its gene; joined to a rigid
        # element it follows that element's point
        self.attach = {}                     # station root -> (element id, point name) of a rigid partner
        idx = {e['id']: k for k, e in enumerate(els)}

        def station_at(k, name):
            e = els[k]
            p = cad.world(e, cad.local_points(e)[name]) if name in cad.local_points(e) else None
            if p is None: return None
            best = None
            for key, q in self.st_pos.items():
                if key[0] == k:
                    d = float(np.linalg.norm(q - p))
                    if d < 1e-3 and (best is None or d < best[0]): best = (d, key)
            if best is None and name in ('start', 'end'):
                pi = 0 if name == 'start' else len(self.poly[k]) - 1
                return (k, pi, 0 if name == 'start' else len(self.poly[k][pi][0]) - 1)
            return best and best[1]
        for j in doc['joins']:
            ka, kb = idx.get(j['a']['el']), idx.get(j['b']['el'])
            if ka is None or kb is None: continue
            sa = station_at(ka, j['a']['pt']) if ka in self.poly else None
            sb = station_at(kb, j['b']['pt']) if kb in self.poly else None
            if sa and sb: union(sa, sb)
            elif sa: self.attach.setdefault(find(sa), []).append((j['b']['el'], j['b']['pt']))
            elif sb: self.attach.setdefault(find(sb), []).append((j['a']['el'], j['a']['pt']))
        # re-root attachments after all unions
        att = {}
        for r, lst in self.attach.items(): att.setdefault(find(r), []).extend(lst)
        self.attach = att
        self.root = {key: find(key) for key in self.st_pos}
        groups = {}
        for key in self.st_pos: groups.setdefault(find(key), []).append(key)
        self.groups = groups
        for r, mem in groups.items():
            if r in self.attach: continue
            dev = min(self.opt[m[0]]['dev'] for m in mem)
            k0 = mem[0][0]
            planar = any(self.opt[m[0]]['planar'] for m in mem)
            for c in range(2 if planar else 3): self.slots.append(('st', r, c, dev))
        # ---- shape mode
        for k, e in enumerate(els):
            o = self.opt[k]
            if o['mode'] != 'shape': continue
            for key, isint in SHAPE_KEYS[e['type']]:
                if key not in e['params'] or e['params'].get('pts'): continue
                v = float(e['params'][key])
                rng = (o['size'] / 100 * abs(v)) if o['size'] > 0 else (o['_mm'] or 0) if key not in ('angle', 'arc', 'turns', 'n') else 0
                if isint: rng = max(rng, 1.0) if o['size'] > 0 else 0
                if rng > 0: self.slots.append(('par', k, key, rng))
            if o['dev'] > 0:
                for c in range(2 if o['planar'] else 3): self.slots.append(('pos', k, c, o['dev']))
            if o['rot'] > 0:
                for c in ((2,) if o['planar'] else (0, 1, 2)): self.slots.append(('rot', k, c, o['rot']))
        # ---- pixel plates
        self.pix = {}
        for k, e in enumerate(els):
            o = self.opt[k]
            if o['mode'] != 'pixels' or e['type'] != 'plate': continue
            nx, ny = o['pix']
            keep = self._probe_cells(k, nx, ny)
            self.pix[k] = (nx, ny, keep)
            for b in range(ny):
                for a in range(nx):
                    if (a, b) not in keep: self.slots.append(('pix', k, b * nx + a, 1.0))
        # ---- planar (pcb) elements: outline corners, vias and the feed move; or the top is cut into pixels
        self.pcbst = {}
        for k, e in enumerate(els):
            if e['type'] != 'pcb': continue
            o = self.opt[k]; p = e['params']
            if o['mode'] == 'shape' and o['size'] > 0:
                for c in (0, 1): self.slots.append(('pcbscale', k, c, o['size'] / 100))
            if o['mode'] == 'free' and o['dev'] > 0:
                st = []
                for ri, R in enumerate(p.get('top') or []):
                    for li, L in enumerate([R.get('outer') or []] + list(R.get('holes') or [])):
                        Lp = PCB._loop(L); n = len(Lp)
                        if n < 3: continue
                        nodes = (R.get('nodes') if li == 0 else (R.get('hnodes') or [[]] * (li))[li - 1] if len(R.get('hnodes') or []) >= li else []) or []
                        nodes = sorted({int(i) % n for i in nodes})
                        if len(nodes) < 3:                 # a circle or a smooth outline: stations spread along it
                            nodes = sorted(set(nodes) | {int(round(i * n / max(3, o['pts']))) % n for i in range(max(3, o['pts']))})
                        st.append((ri, li, nodes))
                        for ni in range(len(nodes)):
                            for c in (0, 1): self.slots.append(('pcbn', k, (len(st) - 1, ni, c), o['dev']))
                self.pcbst[k] = st
                for vi in range(len(p.get('vias') or [])):
                    for c in (0, 1): self.slots.append(('pcbv', k, (vi, c), o['dev']))
                if p.get('feed'):
                    for c in (0, 1): self.slots.append(('pcbf', k, c, o['dev'] / 2))
            if o['mode'] == 'pixels':
                q = dict(e, params=dict(p, pix=o['pix'], mask=None))
                g = PCB.pixel_grid(q)
                if g is None: continue
                nx, ny, x0, y0, dx, dy, inside = g
                keep = set()
                for v in list(p.get('vias') or []) + ([p['feed']] if p.get('feed') else []):
                    a, b = int((float(v[0]) - x0) // dx), int((float(v[1]) - y0) // dy)
                    if 0 <= a < nx and 0 <= b < ny: keep.add(b * nx + a)
                self.pix[k] = (nx, ny, keep)
                for b in range(ny):
                    for a in range(nx):
                        if inside[b, a] and b * nx + a not in keep: self.slots.append(('pix', k, b * nx + a, 1.0))
        # ---- dimensions with a range
        for i, d in enumerate(doc.get('dims') or []):
            r = float(d.get('evo') or 0)
            if r > 0 and d.get('drive', True): self.slots.append(('dim', i, None, r))
        self.n = len(self.slots)
        # wire ends sitting on a plate stay on that plate (e.g. a patch feed probe)
        self.on_plate = []
        for k, out in self.poly.items():
            for pi, (W, st) in enumerate(out):
                for ii in (0, len(W) - 1):
                    for j, pl in enumerate(els):
                        if pl['type'] != 'plate': continue
                        R = cad.rotation(pl['rot']); l = (W[ii] - np.asarray(pl['pos'], float)) @ R
                        W_, H_ = float(pl['params']['W']), float(pl['params']['H'])
                        if abs(l[2]) < 1e-3 and abs(l[0]) <= W_ / 2 + 1e-3 and abs(l[1]) <= H_ / 2 + 1e-3:
                            self.on_plate.append((k, pi, ii, j, l))

    def _probe_cells(self, k, nx, ny):
        """pixels under a wire end that touches the plate stay metal (the feed probe needs them)."""
        els = self.doc['elements']; pl = els[k]
        R = cad.rotation(pl['rot']); o = np.asarray(pl['pos'], float)
        W, H = float(pl['params']['W']), float(pl['params']['H'])
        keep = set()
        for j, e in enumerate(els):
            if j == k or e['type'] not in cad.PATH_TYPES: continue
            for Q in cad.local_paths(e):
                P = cad.world(e, Q)
                for q in (P[0], P[-1]):
                    l = (q - o) @ R
                    if abs(l[2]) < 1e-3 and abs(l[0]) <= W / 2 + 1e-3 and abs(l[1]) <= H / 2 + 1e-3:
                        keep.add((min(nx - 1, max(0, int((l[0] + W / 2) / W * nx))), min(ny - 1, max(0, int((l[1] + H / 2) / H * ny)))))
        return keep

    def _apply_pcb(self, e, k, genes):
        p = e['params']
        sx = sy = 1.0
        for kind, b, v in genes:
            if kind == 'pcbscale':
                if b == 0: sx = 1 + v
                else: sy = 1 + v
        if sx != 1 or sy != 1:                  # stretch the top layer (with its vias and feed) about its centre
            bb = PCB.bbox(PCB.drawn_regions(e, 'top')); c = (bb[0] + bb[1]) / 2
            f = lambda q: [float(c[0] + (q[0] - c[0]) * sx), float(c[1] + (q[1] - c[1]) * sy)]
            for R in p.get('top') or []:
                R['outer'] = [f(q) for q in R['outer']]; R['holes'] = [[f(q) for q in H] for H in R.get('holes') or []]
            p['vias'] = [f(v) for v in p.get('vias') or []]
            if p.get('feed'): p['feed'] = f(p['feed'])
        st = self.pcbst.get(k)
        if st:
            off = {}
            for kind, b, v in genes:
                if kind == 'pcbn': off.setdefault((b[0], b[1]), [0.0, 0.0])[b[2]] = v
            for si, (ri, li, nodes) in enumerate(st):
                R = p['top'][ri]
                L = PCB._loop(R['outer'] if li == 0 else R['holes'][li - 1]); n = len(L)
                D = np.array([off.get((si, ni), [0.0, 0.0]) for ni in range(len(nodes))])
                seg = np.linalg.norm(np.roll(L, -1, axis=0) - L, axis=1); s = np.concatenate([[0], np.cumsum(seg)]); T = s[-1]
                sn = np.array([s[i] for i in nodes]); sn_ = np.concatenate([sn - T, sn, sn + T]); D_ = np.vstack([D, D, D])
                disp = np.stack([np.interp(s[:-1], sn_, D_[:, c]) for c in (0, 1)], -1)
                newL = (L + disp).round(4).tolist()
                if li == 0: R['outer'] = newL
                else: R['holes'][li - 1] = newL
        for kind, b, v in genes:
            if kind == 'pcbv': p['vias'][b[0]][b[1]] = float(p['vias'][b[0]][b[1]]) + v
            elif kind == 'pcbf': p['feed'] = list(p['feed']); p['feed'][b] = float(p['feed'][b]) + v
        p['sketch_stale'] = True

    def apply(self, x):
        doc = copy.deepcopy(self.doc)
        els = doc['elements']
        stoff, posoff, pcbg = {}, {}, {}
        for (kind, a, b, sc), v in zip(self.slots, x):
            if kind == 'st': stoff.setdefault(a, np.zeros(3))[b] = v * sc
            elif kind == 'pos': posoff.setdefault(a, np.zeros(3))[b] = v * sc
            elif kind == 'rot':
                els[a]['rot'] = list(els[a]['rot']); els[a]['rot'][b] = float(els[a]['rot'][b]) + v * sc
            elif kind == 'par':
                p = els[a]['params']; isint = dict(SHAPE_KEYS[els[a]['type']]).get(b)
                nv = float(p[b]) + v * sc
                p[b] = max(1, int(round(nv))) if isint else max(0.05, nv)
                if b == 'arc': p[b] = min(360.0, p[b])
            elif kind == 'pix':
                k = a; nx, ny, keep = self.pix[k]
                p = els[k]['params']
                if not p.get('mask') or len(p['mask']) != nx * ny: p['mask'] = [1] * (nx * ny); p['pix'] = [nx, ny]
                p['mask'][b] = 1 if v > -0.5 else 0
            elif kind in ('pcbscale', 'pcbn', 'pcbv', 'pcbf'):
                pcbg.setdefault(a, []).append((kind, b, v * sc))
            elif kind == 'dim':
                d = doc['dims'][a]; d['value'] = float(d['value']) + v * sc
                if d.get('axis') == 'angle': d['value'] = min(180.0, max(0.0, d['value']))
                elif d.get('axis') == 'd': d['value'] = max(0.0, d['value'])
        for k, genes in pcbg.items(): self._apply_pcb(els[k], k, genes)
        for k in self.pix:
            if els[k]['type'] == 'pcb': els[k]['params']['pix'] = list(self.pix[k][:2])
        for k in self.pix:                   # all pixels on when no gene touched them
            p = els[k]['params']; nx, ny, _ = self.pix[k]
            if not p.get('mask'): p['mask'] = [1] * (nx * ny); p['pix'] = [nx, ny]
        for k, off in posoff.items():
            els[k]['pos'] = [float(a + b) for a, b in zip(els[k]['pos'], self.frame[k] @ off)]
        # rigid elements meet the joints and driving dimensions again (free-form wires adapt after)
        fids = {els[k]['id'] for k in self.poly}
        joins = doc['joins']
        doc['joins'] = [j for j in joins if j['a']['el'] not in fids and j['b']['el'] not in fids]
        cad.solve(doc, locked=fids, iters=150)
        doc['joins'] = joins
        # free-form wires: station offsets, interpolated along every polyline
        newpos = {}
        for r, mem in self.groups.items():
            if r in self.attach:
                el, pt = self.attach[r][0]
                q = cad.ref_point(doc, {'el': el, 'pt': pt})
                base = self.st_pos[mem[0]]
                newpos[r] = (q - base) if q is not None else np.zeros(3)
            else:
                off = stoff.get(r, np.zeros(3))
                k0 = mem[0][0]
                if any(self.opt[m[0]]['planar'] for m in mem): off = np.array([off[0], off[1], 0.0])
                newpos[r] = self.frame[k0] @ off
        for k, out in self.poly.items():
            e = els[k]
            newP = []
            for pi, (W, st) in enumerate(out):
                s = _arclen(W)
                D = np.array([newpos[self.root[(k, pi, i)]] for i in st])
                disp = np.stack([np.interp(s, s[st], D[:, c]) for c in range(3)], -1)
                newP.append(W + disp)
            for (kk, pi, ii, j, l0) in self.on_plate:    # back onto the plate plane, inside the plate
                if kk != k: continue
                pl = els[j]; R = cad.rotation(pl['rot']); o = np.asarray(pl['pos'], float)
                l = (newP[pi][ii] - o) @ R
                l[2] = 0; l[0] = np.clip(l[0], -pl['params']['W'] / 2, pl['params']['W'] / 2); l[1] = np.clip(l[1], -pl['params']['H'] / 2, pl['params']['H'] / 2)
                newP[pi][ii] = o + R @ l
            if doc['ground'] == 'perfect':
                for P in newP: P[:, 2] = np.maximum(P[:, 2], 0)
            R = cad.rotation(e['rot']); o = np.asarray(e['pos'], float)
            loc = [((P - o) @ R).round(4) for P in newP]
            if e['type'] == 'custom':
                p = e['params']
                oldW = [W for W, _ in out]
                # named snap points follow the station they sit on
                for nm, v in list((p.get('pts') or {}).items()):
                    w0 = cad.world(self.doc['elements'][k], [v[0], v[1], v[2] if len(v) > 2 else 0])
                    for pi, W in enumerate(oldW):
                        hit = np.where(np.linalg.norm(W - w0, axis=1) < 1e-3)[0]
                        if len(hit): p['pts'][nm] = loc[pi][hit[0]].tolist(); break
                p['paths'] = [L.tolist() for L in loc]
                p['sketch_stale'] = True
            else:
                e['type'] = 'rod'
                e['params'] = {'pts': loc[0].tolist(), 'bend_r': e['params'].get('bend_r', 0)}
        return doc


def _resample(P, n):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, P[:, c]) for c in range(3)], -1)


def valid(doc):
    """cheap geometry checks before NEC: no degenerate pieces, no wire folding back on itself."""
    for e in doc['elements']:
        if e['type'] in cad.PATH_TYPES:
            try: Ps = cad.local_paths(e)
            except Exception: return False
            for P in Ps:
                if not np.isfinite(P).all(): return False
                seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
                if e['type'] == 'rod' and (seg < 0.75 * e['d']).any(): return False
                if len(P) > 2:
                    d = np.diff(P, axis=0)
                    cosang = (d[:-1] * d[1:]).sum(1) / np.maximum(seg[:-1] * seg[1:], 1e-12)
                    if (cosang < -0.94).any(): return False   # folds back on itself (> 160 deg)
        elif e['type'] == 'plate' and cad.pixel_mask(e) is not None and not cad.pixel_mask(e)[2].any():
            return False
        elif e['type'] == 'pcb':
            try:
                p = e['params']
                if not PCB.regions(e, 'top') or PCB.h_of(e) <= 0: return False
                for R in p.get('top') or []:
                    for L in [R['outer']] + list(R.get('holes') or []):
                        if PCB.self_crosses(PCB._loop(L)): return False
                for q in list(p.get('vias') or []) + ([p['feed']] if p.get('feed') else []):
                    if not PCB.metal_at(e, 'top', q[0], q[1]): return False
                    if PCB.regions(e, 'gnd') and not PCB.metal_at(e, 'gnd', q[0], q[1]): return False
            except Exception:
                return False
    return True


# ------------------------------------------------------------------ evaluation (worker processes)

_W = {}


def _init(space, goal):
    _W['space'], _W['goal'] = space, goal


def ems_run(doc, freqs, goal, stop):
    """one openEMS evaluation in a scratch folder inside the outputs folder; same result shape as NEC."""
    import ems, shutil
    sim = os.path.join(goal['simroot'], 'eval_%d' % os.getpid())
    try:
        info = ems.simulate(doc, goal['exe'], sim, log=lambda m: None, stop=stop)
        Z = ems.results_z(sim, freqs)
        f0 = doc['f0'] * 1e6
        th, ph, Eth, Eph, prad = ems.nf2ff(sim, info['faces'], f0, GRID, doc['ground'] == 'perfect')
        t, u = ems.read_probe(os.path.join(sim, 'port_ut_1')); ti, i = ems.read_probe(os.path.join(sim, 'port_it_1'))
        p_in = 0.5 * float(np.real(ems.dft(t, u, [f0])[0] * np.conj(ems.dft(ti, i, [f0])[0])))
        pat = ems.pattern_from_ff(th, ph, Eth, Eph, max(p_in, 1e-30), doc['ground'] == 'perfect')
    finally:
        shutil.rmtree(sim, ignore_errors=True)
    G = np.array(pat['total'])
    return [{'f': f, 'z': z, 'pattern': pat, 'total': {'max': float(G.max())}} for f, z in zip(freqs, Z)]


def evaluate(space, goal, x, stop=None):
    doc = space.apply(np.asarray(x, float))
    if not valid(doc): return -1e9, None
    f0, bw = doc['f0'], doc['bw']
    freqs = [f0] if not bw else [f0 - bw / 2, f0, f0 + bw / 2]
    try:
        if doc.get('engine') == 'openems': res = ems_run(doc, freqs, goal, stop)
        else: res, _ = cad.nec_run(doc, freqs, PAT.rp_for(GRID, doc['ground'] == 'perfect'), pat_step=GRID)
    except Exception as e:
        if stop is not None and stop.is_set(): raise Stopped()
        return -1e9, None
    s = max(core.swr(r['z'], doc['z0']) for r in res)
    pol = goal.get('pol', 'total')
    if goal['kind'] == 'dir':
        th, ph = 90 - goal['el'], goal['az']
        val = min(core.pat_sample(r['pattern'], pol, th, ph) for r in res)
        info = {'goal': val, 'gmax': max(r['total']['max'] for r in res)}
        f = val
    elif goal['kind'] == 'target':
        t = goal['target']
        errs = []
        for r in res:
            pat = r['pattern']
            th = np.array(pat['th']); G = np.array(pat[pol] if pat.get(pol) is not None else pat['total'])
            T = PAT.target_grid(t, th, np.array(pat['ph']))
            errs.append(PAT.match(G, T, th, t['floor'], doc['ground'] == 'perfect')[0])
        val = max(errs)
        info = {'goal': val, 'gmax': max(r['total']['max'] for r in res)}
        f = -val
    else:  # best match only
        val = s
        info = {'goal': s, 'gmax': max(r['total']['max'] for r in res)}
        f = -_ml(min(s, 999)) * 3 + 0.05 * info['gmax']
    f -= 2.0 * max(0.0, _ml(min(s, 999)) - _ml(doc['swr'])) if goal['kind'] != 'match' else 0
    info['swr'] = s
    return f, info


def _eval_chunk(xs):
    return [evaluate(_W['space'], _W['goal'], x) for x in xs]


class Stopped(Exception):
    pass


class CadEvolver:
    kind = 'cad'

    def __init__(self, doc, goal, pop=40, gens=30, jobs=0, seed=0, progress=None, stop=None, log=None):
        self.doc, self.goal = doc, goal
        self.space = Space(doc)
        self.pop, self.gens = max(8, int(pop)), max(1, int(gens))
        self.jobs = int(jobs) or max(1, (mp.cpu_count() or 2) - 1)
        self.rng = np.random.default_rng(int(seed) or None)
        self.progress, self.stop, self.log = progress or (lambda p: None), stop, log or (lambda m: None)
        self.pool = None
        self.best = None          # (doc, fitness, info)
        self.evals = 0
        self.t0 = time.time()

    def _stopped(self):
        return self.stop is not None and self.stop.is_set()

    def _pool(self):
        if self.doc.get('engine') == 'openems': return None
        return mp.get_context('spawn').Pool(self.jobs, initializer=_init, initargs=(self.space, self.goal))

    def _map(self, X):
        """evaluate genomes in chunks; a chunk that hangs (stuck NEC) is killed and scored invalid."""
        if self.doc.get('engine') == 'openems':      # openEMS is multithreaded itself: one at a time, in-process
            res = []
            for j, x in enumerate(X):
                if self._stopped(): raise Stopped()
                res.append(evaluate(self.space, self.goal, x, self.stop))
                self.evals += 1
                self.log('openEMS evaluation %d/%d: %s' % (j + 1, len(X), 'invalid' if res[-1][0] < -1e8 else 'fitness %.2f' % res[-1][0]))
            return res
        chunk = max(1, len(X) // (self.jobs * 3))
        parts = [X[i:i + chunk] for i in range(0, len(X), chunk)]
        jobs = [(i, self.pool.apply_async(_eval_chunk, (p,))) for i, p in enumerate(parts)]
        out = [None] * len(parts)
        t = time.time()
        per = None
        while jobs:
            if self._stopped(): raise Stopped()
            time.sleep(0.2)
            rest = []
            for i, j in jobs:
                if j.ready():
                    try: out[i] = j.get()
                    except Exception: out[i] = [(-1e9, None)] * len(parts[i])
                    per = per or (time.time() - t) / max(1, len(parts[i]))
                else: rest.append((i, j))
            jobs = rest
            limit = max(60.0, 20 * chunk * (per or 3.0))
            if jobs and time.time() - t > limit:
                self.log('NEC watchdog: a simulation hung for %.0f s, restarting the workers' % limit)
                self.pool.terminate(); self.pool = self._pool()
                for i, _ in jobs: out[i] = [(-1e9, None)] * len(parts[i])
                jobs = []
        res = [r for o in out for r in o]
        self.evals += len(res)
        return res

    def run(self):
        sp = self.space
        if sp.n == 0:
            raise ValueError('nothing may change: give at least one element evolution freedom (shape, free-form or pixels), or a dimension a range')
        self.log('Evolving %s: %d genes, population %d, %d generations, %d cores' % (self.doc['name'], sp.n, self.pop, self.gens, self.jobs))
        self.pool = self._pool()
        try:
            X = [np.zeros(sp.n)] + [self.rng.uniform(-0.5, 0.5, sp.n) for _ in range(self.pop - 1)]
            R = self._map(X)
            hist = []
            for gen in range(self.gens):
                fit = np.array([r[0] for r in R])
                order = np.argsort(-fit)
                b = order[0]
                if fit[b] > -1e8 and (self.best is None or fit[b] > self.best[1]):
                    self.best = (sp.apply(X[b]), float(fit[b]), R[b][1], X[b].copy())
                if self.best:
                    i = self.best[2]
                    hist.append(self.best[1])
                    valid_n = int((fit > -1e8).sum())
                    self.log('Generation %d/%d: fitness %.2f, %s %.2f, max gain %.2f dBi, SWR %.2f, valid %d/%d' % (
                        gen + 1, self.gens, self.best[1], {'dir': 'gain to target', 'target': 'pattern error', 'match': 'SWR'}[self.goal['kind']],
                        i['goal'], i['gmax'], i['swr'], valid_n, len(X)))
                    self.progress({'gen': gen + 1, 'gens': self.gens, 'msg': 'Generation %d/%d' % (gen + 1, self.gens),
                                   'info': {'gmax': i['gmax'], 'swr': i['swr'], 'goal': i['goal'], 'kind': self.goal['kind']},
                                   'hist': hist[-200:], 'wires': cad.wires_json(self.best[0]), 'evals': self.evals,
                                   '_pat': (lambda doc: (lambda: cad.nec_run(doc, [doc['f0']], PAT.rp_for(5, doc['ground'] == 'perfect'), pat_step=5)[0][0]['pattern']))(self.best[0]) if hist[-2:-1] != [self.best[1]] else None,
                                   'rate': self.evals / max(1e-6, time.time() - self.t0)})
                if gen == self.gens - 1 or self._stopped(): break
                sigma = 0.35 * (1 - gen / self.gens) + 0.05
                elite = [X[i] for i in order[:2]]
                kids = []
                while len(kids) < self.pop - len(elite):
                    for attempt in range(6):
                        p1 = X[max(self.rng.integers(0, len(X), 3), key=lambda i: fit[i])]
                        p2 = X[max(self.rng.integers(0, len(X), 3), key=lambda i: fit[i])]
                        a, bb = np.minimum(p1, p2), np.maximum(p1, p2)
                        c = a - 0.3 * (bb - a) + self.rng.random(sp.n) * 1.6 * (bb - a) if self.rng.random() < 0.8 else p1.copy()
                        m = self.rng.random(sp.n) < max(0.2, 2.0 / sp.n)
                        c = np.clip(c + m * self.rng.normal(0, sigma, sp.n), -1, 1)
                        if valid(sp.apply(c)): break
                    kids.append(c)
                Rk = self._map(kids)
                R = [R[i] for i in order[:2]] + Rk
                X = elite + kids
            return self.best
        except Stopped:
            self.log('Stopped, keeping the best shape found so far')
            return self.best
        finally:
            try:
                if self.pool: self.pool.terminate()
            except Exception: pass
