# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""Planar (PCB / patch) elements drawn in the sketch editor: a top copper layer at z = h and a ground
layer at z = 0, both made of filled areas (an outline with holes), shorting vias between them, a
coax probe feed from the ground up to the top and a substrate of thickness h in between.

params of a 'pcb' element (local frame, mm):
  top, gnd : [{'outer': [[x, y], ...], 'holes': [[[x, y], ...], ...], 'nodes': [i, ...], 'hnodes': [[i, ...]]}]
             'nodes' are the drawn corner points of the outline (the evolution moves these)
  vias     : [[x, y], ...]          via_d : via diameter
  feed     : [x, y] or None        h : substrate thickness     sub : substrate material ('' = air)
  mesh     : NEC grid cell (0 = auto)
  pix, mask: top layer cut into nx x ny pixels over its bounding box (pixel evolution)
"""
import math
import numpy as np

LAYERS = ('top', 'gnd')


def _loop(L):
    A = np.array([[float(q[0]), float(q[1])] for q in L], float)
    if len(A) > 1 and np.allclose(A[0], A[-1]): A = A[:-1]
    return A


def drawn_regions(e, layer):
    """[(outer Nx2, [holes])] as drawn (or as evolved), no pixels."""
    out = []
    for R in e['params'].get(layer) or []:
        o = _loop(R.get('outer') or [])
        if len(o) < 3 or not np.isfinite(o).all(): continue
        out.append((o, [h for h in (_loop(H) for H in R.get('holes') or []) if len(h) >= 3]))
    return out


def area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def in_poly(P, X, Y):
    """even-odd point in polygon, vectorised over X, Y (same shape)."""
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    ins = np.zeros(X.shape, bool)
    n = len(P)
    for i in range(n):
        x0, y0 = P[i]; x1, y1 = P[(i + 1) % n]
        if y0 == y1: continue
        c = ((y0 > Y) != (y1 > Y)) & (X < (x1 - x0) * (Y - y0) / (y1 - y0) + x0)
        ins ^= c
    return ins


def in_regions(regs, X, Y):
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    on = np.zeros(X.shape, bool)
    for o, holes in regs:
        m = in_poly(o, X, Y)
        for h in holes: m &= ~in_poly(h, X, Y)
        on |= m
    return on


def bbox(regs):
    if not regs: return None
    A = np.vstack([o for o, _ in regs])
    return A.min(0), A.max(0)


def pixel_grid(e):
    """(nx, ny, x0, y0, dx, dy, inside[ny, nx]) of the top layer's pixel grid, or None."""
    p = e['params']
    if not p.get('pix'): return None
    regs = drawn_regions(e, 'top')
    bb = bbox(regs)
    if bb is None: return None
    nx, ny = max(1, int(p['pix'][0])), max(1, int(p['pix'][1]))
    (x0, y0), (x1, y1) = bb
    dx, dy = (x1 - x0) / nx, (y1 - y0) / ny
    cx = x0 + (np.arange(nx) + 0.5) * dx; cy = y0 + (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(cx, cy)
    return nx, ny, x0, y0, dx, dy, in_regions(regs, X, Y)


def pixel_on(e):
    """metal pixels [ny, nx] (inside the drawn top area and switched on), or None without pixels."""
    g = pixel_grid(e)
    if g is None: return None
    nx, ny, x0, y0, dx, dy, inside = g
    m = e['params'].get('mask')
    on = np.array(m, bool).reshape(ny, nx) if m is not None and len(m) == nx * ny else np.ones((ny, nx), bool)
    return on & inside


def pixel_rects(e):
    g = pixel_grid(e); on = pixel_on(e)
    nx, ny, x0, y0, dx, dy, _ = g
    out = []
    for b in range(ny):
        a = 0
        while a < nx:
            if not on[b, a]: a += 1; continue
            a1 = a
            while a1 + 1 < nx and on[b, a1 + 1]: a1 += 1
            out.append((x0 + a * dx, y0 + b * dy, x0 + (a1 + 1) * dx, y0 + (b + 1) * dy))
            a = a1 + 1
    return out


def regions(e, layer):
    """the metal of a layer as [(outer, holes)]: the drawn areas, or the switched-on pixels of the top."""
    if layer == 'top' and e['params'].get('pix'):
        return [(np.array([[a, b], [c, b], [c, d], [a, d]], float), []) for a, b, c, d in pixel_rects(e)]
    return drawn_regions(e, layer)


def metal_at(e, layer, x, y):
    return bool(in_regions(regions(e, layer), np.array([x]), np.array([y]))[0])


def h_of(e): return float(e['params'].get('h') or 0)


def check(e, mats):
    p = e['params']
    if h_of(e) <= 0: raise ValueError('%s: the substrate thickness h must be positive' % e['name'])
    if not drawn_regions(e, 'top'): raise ValueError('%s: fill at least one closed area on the Top layer' % e['name'])
    s = p.get('sub') or ''
    if s and (s not in mats or mats[s]['kind'] != 'dielectric'):
        raise ValueError('%s: substrate %s is not a dielectric in the material library' % (e['name'], s))
    if float(p.get('via_d') or 0.8) <= 0: raise ValueError('%s: via diameter must be positive' % e['name'])
    for v in p.get('vias') or []:
        if len(v) < 2 or not all(math.isfinite(float(c)) for c in v[:2]): raise ValueError('%s: bad via' % e['name'])


def local_points(e):
    h = h_of(e)
    pts = {'center': np.array([0.0, 0.0, h / 2])}
    for layer, z, tag in (('top', h, 't'), ('gnd', 0.0, 'g')):
        bb = bbox(drawn_regions(e, layer))
        if bb is None: continue
        (x0, y0), (x1, y1) = bb
        for i, (x, y) in enumerate(((x0, y0), (x1, y0), (x1, y1), (x0, y1))):
            pts['%s%d' % (tag, i + 1)] = np.array([x, y, z])
        pts[tag + 'c'] = np.array([(x0 + x1) / 2, (y0 + y1) / 2, z])
    f = e['params'].get('feed')
    if f: pts['feed'] = np.array([float(f[0]), float(f[1]), 0.0]); pts['feedt'] = np.array([float(f[0]), float(f[1]), h])
    for i, v in enumerate(e['params'].get('vias') or []):
        pts['via%d' % (i + 1)] = np.array([float(v[0]), float(v[1]), h])
    return pts


# ------------------------------------------------------------------ NEC-2: wire grids on both layers, vertical wires for vias and the probe

def _lines(lo, hi, fixed, cell):
    """grid coordinates from lo to hi: every fixed coordinate is a line, the gaps are split evenly
    into cells no larger than `cell`; fixed lines closer than cell / 4 are merged."""
    pts = sorted(set([lo, hi] + [v for v in fixed if lo - 1e-9 <= v <= hi + 1e-9]))
    merged = [pts[0]]
    for v in pts[1:]:
        if v - merged[-1] < cell / 4:
            if v == hi: merged[-1] = hi
            continue
        merged.append(v)
    out = [merged[0]]
    for a, b in zip(merged[:-1], merged[1:]):
        n = max(1, int(math.ceil((b - a) / cell - 1e-9)))
        out += list(a + (b - a) * np.arange(1, n + 1) / n)
    return np.array(out)


def nec_wires(doc, k, lam, world, warn):
    e = doc['elements'][k]
    p = e['params']; h = h_of(e)
    top, gnd = regions(e, 'top'), regions(e, 'gnd')
    if (p.get('sub') or '') and doc['_mats'].get(p.get('sub'), {}).get('eps', 1) > 1.0001:
        warn.append('%s: NEC-2 ignores the %s substrate (the patch resonates too high); the openEMS engine models it' % (e['name'], p['sub']))
    if not gnd: warn.append('%s: no ground layer: the top layer radiates on its own' % e['name'])
    allr = top + gnd
    (x0, y0), (x1, y1) = bbox(allr)
    span = max(x1 - x0, y1 - y0, 1e-3)
    cell = float(p.get('mesh') or 0) or min(lam / 20, span / 12)
    cell = max(cell, span / 60)                                  # at most 60 cells across
    fx, fy = [], []
    for o, hs in allr:                                        # bounding edges of every area
        for P in [o] + hs:
            fx += [float(P[:, 0].min()), float(P[:, 0].max())]; fy += [float(P[:, 1].min()), float(P[:, 1].max())]
    g = pixel_grid(e)
    if g is not None:
        nx, ny, gx0, gy0, dx, dy, _ = g
        fx += list(gx0 + dx * np.arange(nx + 1)); fy += list(gy0 + dy * np.arange(ny + 1))
    pins = [(float(v[0]), float(v[1]), 'via') for v in p.get('vias') or []]
    if p.get('feed'): pins.append((float(p['feed'][0]), float(p['feed'][1]), 'feed'))
    xs = _lines(x0, x1, fx + [q[0] for q in pins], cell)
    ys = _lines(y0, y1, fy + [q[1] for q in pins], cell)
    rg = min(cell / (2 * math.pi), h / 3) / 1000                # equal-area rule, thinner than the gap to the other layer
    # where lines crowd (vias, small holes) every wire at a node stays thinner than 0.3 x the shortest
    # wire there, so no wire's midpoint ends up inside a neighbour (NEC-2 geometry check)
    gx, gy = np.diff(xs), np.diff(ys)
    nx_ = np.minimum(np.concatenate([[np.inf], gx]), np.concatenate([gx, [np.inf]]))
    ny_ = np.minimum(np.concatenate([[np.inf], gy]), np.concatenate([gy, [np.inf]]))
    gmin = np.minimum(ny_[:, None], nx_[None, :])                # [j, i] shortest wire at node, mm
    rad = lambda j0, i0, j1, i1: min(rg, 0.3 * min(gmin[j0, i0], gmin[j1, i1]) / 1000)
    out = []
    cx, cy = (xs[:-1] + xs[1:]) / 2, (ys[:-1] + ys[1:]) / 2
    X, Y = np.meshgrid(cx, cy)                                   # [j, i]
    node = {}
    for layer, regs, z in (('top', top, h), ('gnd', gnd, 0.0)):
        if not regs: continue
        on = in_regions(regs, X, Y)
        if not on.any(): continue
        nj, ni = on.shape
        isOn = lambda i, j: 0 <= i < ni and 0 <= j < nj and bool(on[j, i])
        nd = np.zeros((len(ys), len(xs)), bool)
        for i in range(len(xs)):
            for j in range(len(ys)):
                if j < len(ys) - 1 and (isOn(i - 1, j) or isOn(i, j)):
                    out.append((k, ('grid', 1), world(e, [xs[i], ys[j], z]) / 1000, world(e, [xs[i], ys[j + 1], z]) / 1000, rad(j, i, j + 1, i)))
                    nd[j, i] = nd[j + 1, i] = True
                if i < len(xs) - 1 and (isOn(i, j - 1) or isOn(i, j)):
                    out.append((k, ('grid', 1), world(e, [xs[i], ys[j], z]) / 1000, world(e, [xs[i + 1], ys[j], z]) / 1000, rad(j, i, j, i + 1)))
                    nd[j, i] = nd[j, i + 1] = True
        node[layer] = nd
    rv = float(p.get('via_d') or 0.8) / 2000
    for x, y, kind in pins:
        i = int(np.argmin(np.abs(xs - x))); j = int(np.argmin(np.abs(ys - y)))
        ok_t = 'top' in node and node['top'][j, i]
        ok_g = 'gnd' in node and node['gnd'][j, i]
        what = 'the feed probe' if kind == 'feed' else 'a via at (%.1f, %.1f)' % (x, y)
        if not ok_t: raise ValueError('%s: %s is not under top-layer copper' % (e['name'], what))
        if not ok_g:
            if kind == 'feed': raise ValueError('%s: the feed probe needs ground-layer copper under it (fill the ground layer)' % e['name'])
            raise ValueError('%s: %s has no ground-layer copper under it' % (e['name'], what))
        a, b = world(e, [xs[i], ys[j], 0.0]) / 1000, world(e, [xs[i], ys[j], h]) / 1000
        out.append((k, 'feed' if kind == 'feed' else ('grid', 1), a, b, min(rv, h / 8000, 0.3 * gmin[j, i] / 1000)))
    if not p.get('feed') and e.get('feed'):
        raise ValueError('%s: place the feed with the Probe tool in the drawing' % e['name'])
    return out


# ------------------------------------------------------------------ SVG (top view, 1 unit = 1 mm)

def svg(e, title=None):
    p = e['params']
    regs = {L: regions(e, L) for L in LAYERS}
    allr = regs['top'] + regs['gnd']
    (x0, y0), (x1, y1) = bbox(allr)
    vd = float(p.get('via_d') or 0.8)
    m = max(2.0, 0.05 * max(x1 - x0, y1 - y0))
    W, H = (x1 - x0) + 2 * m, (y1 - y0) + 2 * m
    f = lambda v: ('%.4f' % v).rstrip('0').rstrip('.')

    def path(o, holes):
        d = ''
        for P in [o] + holes:
            d += 'M' + ' L'.join('%s,%s' % (f(q[0]), f(-q[1])) for q in P) + ' Z '
        return '<path d="%s"/>' % d.strip()
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" '
           'width="%smm" height="%smm" viewBox="%s %s %s %s">' % (f(W), f(H), f(x0 - m), f(-y1 - m), f(W), f(H)),
           '<title>%s</title>' % (title or e['name']).replace('&', '&amp;').replace('<', '&lt;'),
           '<desc>Anten-O-Matic planar layout, top view, units mm, substrate %s h = %s mm. The ground layer is drawn as seen from the top.</desc>'
           % (p.get('sub') or 'air', f(h_of(e)))]
    for L, col, nm in (('gnd', '#2f7fb8', 'Ground (bottom)'), ('top', '#c77b3a', 'Top copper')):
        if not regs[L]: continue
        out.append('<g id="%s" inkscape:groupmode="layer" inkscape:label="%s" fill="%s" fill-rule="evenodd" stroke="none"%s>'
                   % (L, nm, col, ' opacity="0.55"' if L == 'gnd' else ''))
        out += [path(o, hs) for o, hs in regs[L]]
        out.append('</g>')
    if p.get('vias') or p.get('feed'):
        out.append('<g id="drill" inkscape:groupmode="layer" inkscape:label="Vias and feed (drill)" fill="#ffffff" stroke="#000000" stroke-width="%s">' % f(vd / 8))
        for v in p.get('vias') or []: out.append('<circle cx="%s" cy="%s" r="%s"/>' % (f(v[0]), f(-v[1]), f(vd / 2)))
        if p.get('feed'):
            fx, fy = p['feed']
            out.append('<circle cx="%s" cy="%s" r="%s" fill="#ffd400"/>' % (f(fx), f(-fy), f(max(vd, 1.0) / 2)))
        out.append('</g>')
    out.append('</svg>')
    return '\n'.join(out) + '\n'


def write_svgs(doc, out_dir):
    import os
    files = []
    for e in doc['elements']:
        if e['type'] != 'pcb': continue
        safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in e['name'])[:40] or 'pcb'
        fn = os.path.join(out_dir, 'layout_%s.svg' % safe)
        with open(fn, 'w', encoding='utf-8') as fh: fh.write(svg(e))
        files.append(fn)
    return files


def describe(e):
    p = e['params']
    def a(L): return sum(abs(area(o)) - sum(abs(area(h)) for h in hs) for o, hs in regions(e, L))
    s = 'planar: top %.1f mm², ground %.1f mm², %d via(s) ⌀%.2f, substrate %s h %.2f mm' % (
        a('top'), a('gnd'), len(p.get('vias') or []), float(p.get('via_d') or 0.8), p.get('sub') or 'air', h_of(e))
    if p.get('feed'): s += ', probe feed at (%.2f, %.2f)' % tuple(p['feed'])
    if p.get('pix'):
        on = pixel_on(e); s += ', top cut into %dx%d pixels, %d metal' % (on.shape[1], on.shape[0], int(on.sum()))
    return s


# ------------------------------------------------------------------ openEMS: polygons on two planes, vias as PEC lines, a probe port, a substrate box

def _mesh_corners(W, a1, a2, tol=1e-6):
    """points of a polygon that deserve their own mesh lines: the extremes of the loop and the ends of
    straight edges along the mesh axes. The vertices of curved edges (circles, arcs) are left to the
    normal mesh resolution: giving each one a line made a ring patch need millions of tiny cells."""
    pts = [W[int(np.argmin(W[:, a1]))], W[int(np.argmax(W[:, a1]))], W[int(np.argmin(W[:, a2]))], W[int(np.argmax(W[:, a2]))]]
    n = len(W); span = float(np.ptp(W[:, [a1, a2]], axis=0).max()) or 1.0
    for i in range(n):
        A, B = W[i], W[(i + 1) % n]
        d = np.abs(B - A)
        if (d[a1] < tol or d[a2] < tol) and max(d[a1], d[a2]) > 0.02 * span:
            pts += [A, B]
    return pts


def ems_parts(doc, k, world, rotation):
    """(boxes, curves, port) in the format of ems.geometry for one planar element."""
    e = doc['elements'][k]; p = e['params']; h = h_of(e)
    R = rotation(e['rot'])
    nd = int(np.argmax(np.abs(R[:, 2])))
    if not np.allclose(np.abs(R).round(6), np.round(np.abs(R)), atol=1e-6):
        raise ValueError('%s is rotated off the x/y/z axes: openEMS needs planar layers aligned with the axes' % e['name'])
    a1, a2 = (nd + 1) % 3, (nd + 2) % 3
    boxes, curves = [], []
    metal = doc['_mats'][e['mat']]
    for layer, z in (('top', h), ('gnd', 0.0)):
        regs = regions(e, layer)
        loops = [o for o, _ in regs] + [hh for _, hs in regs for hh in hs]
        for o, hs in regs:
            q0 = o[0]
            depth = sum(1 for L in loops if L is not o and abs(area(L)) > abs(area(o)) and in_poly(L, np.array([q0[0]]), np.array([q0[1]]))[0])
            nodes = None
            for kind, P, prio in [('metal', o, 10 + 2 * depth)] + [('hole', hh, 11 + 2 * depth) for hh in hs]:
                W = np.array([world(e, [q[0], q[1], z]) for q in P])
                corners = _mesh_corners(W, a1, a2)
                boxes.append({'el': k, 'lo': W.min(0), 'hi': W.max(0), 'mat': metal, 'name': '%s_%s' % (e['name'], layer),
                              'poly': {'nd': nd, 'elev': float(W[0][nd]), 'v': [(float(q[a1]), float(q[a2])) for q in W],
                                       'kind': kind, 'prio': prio}, 'fix': list(corners)})
    sub = p.get('sub') or ''
    if sub:
        regs = drawn_regions(e, 'gnd') or drawn_regions(e, 'top')
        (x0, y0), (x1, y1) = bbox(regs)
        Q = np.array([world(e, [x, y, z]) for x in (x0, x1) for y in (y0, y1) for z in (0.0, h)])
        boxes.append({'el': k, 'lo': Q.min(0), 'hi': Q.max(0), 'mat': doc['_mats'][sub], 'name': '%s_substrate' % e['name']})
    for v in p.get('vias') or []:
        curves.append({'el': k, 'pts': np.array([world(e, [float(v[0]), float(v[1]), 0.0]), world(e, [float(v[0]), float(v[1]), h])])})
    port = None
    if p.get('feed'):
        f = p['feed']
        if not metal_at(e, 'top', f[0], f[1]): raise ValueError('%s: the feed probe is not under top-layer copper' % e['name'])
        a, b = world(e, [float(f[0]), float(f[1]), 0.0]), world(e, [float(f[0]), float(f[1]), h])
        port = {'a': a, 'b': b, 'ax': nd, 'el': k}
    return boxes, curves, port


def self_crosses(P):
    """does a closed outline cross itself (two non-neighbouring edges intersect)?"""
    n = len(P)
    if n < 4: return False
    A = P; B = np.roll(P, -1, axis=0)
    for i in range(n):
        a, b = A[i], B[i]
        j = np.arange(i + 2, n if i > 0 else n - 1)
        if not len(j): continue
        c, d = A[j], B[j]
        o1 = np.sign((b[0] - a[0]) * (c[:, 1] - a[1]) - (b[1] - a[1]) * (c[:, 0] - a[0]))
        o2 = np.sign((b[0] - a[0]) * (d[:, 1] - a[1]) - (b[1] - a[1]) * (d[:, 0] - a[0]))
        o3 = np.sign((d[:, 0] - c[:, 0]) * (a[1] - c[:, 1]) - (d[:, 1] - c[:, 1]) * (a[0] - c[:, 0]))
        o4 = np.sign((d[:, 0] - c[:, 0]) * (b[1] - c[:, 1]) - (d[:, 1] - c[:, 1]) * (b[0] - c[:, 0]))
        if np.any((o1 * o2 < 0) & (o3 * o4 < 0)): return True
    return False
