# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""
Minimal STEP (ISO 10303-21, AP214) writer for wire antennas: every arm becomes a solid
tube (wire diameter) along its centre line, straight runs as cylinders, bends as tori with
the real bend radius, flat end caps.  Pure Python, no CAD kernel needed.  Units: mm.
"""
import math, time, unicodedata
import numpy as np


class _W:
    def __init__(self):
        self.lines, self.n = [], 0

    def add(self, s):
        self.n += 1
        self.lines.append('#%d=%s;' % (self.n, s))
        return '#%d' % self.n

    def pt(self, p):
        return self.add('CARTESIAN_POINT(\'\',(%s))' % ','.join(self.num(v) for v in p))

    def dir(self, v):
        v = np.asarray(v, float); v = v / np.linalg.norm(v)
        return self.add('DIRECTION(\'\',(%s))' % ','.join(self.num(c) for c in v))

    @staticmethod
    def num(x):
        x = float(x)
        if abs(x) < 1e-12: x = 0.0
        s = '%.10G' % x
        if 'E' in s:
            m, e = s.split('E'); s = (m if '.' in m else m + '.') + 'E' + e
        elif '.' not in s:
            s += '.'
        return s

    def ax2(self, o, z, x):
        z = np.asarray(z, float); z = z / np.linalg.norm(z)
        x = np.asarray(x, float); x = x - (x @ z) * z; x = x / np.linalg.norm(x)
        return self.add('AXIS2_PLACEMENT_3D(\'\',%s,%s,%s)' % (self.pt(o), self.dir(z), self.dir(x)))


def _perp(v):
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    p = np.cross(v, a)
    return p / np.linalg.norm(p)


def _rot(v, axis, ang):
    return v * math.cos(ang) + np.cross(axis, v) * math.sin(ang) + axis * (axis @ v) * (1 - math.cos(ang))


def tube_solid(w, pieces, r, name):
    """pieces: [('line', A, B) | ('arc', T1, T2, centre, axis, angle)] in mm, continuous."""
    # junction frames: centre point, tangent, seam reference (parallel transported)
    first = pieces[0]
    t0 = (first[2] - first[1]) if first[0] == 'line' else np.cross(first[4], first[1] - first[3])
    t0 = t0 / np.linalg.norm(t0)
    nvec = _perp(t0)
    juncs = [(np.asarray(first[1], float), t0, nvec)]
    for pc in pieces:
        c, t, nv = juncs[-1]
        if pc[0] == 'line':
            juncs.append((np.asarray(pc[2], float), t, nv))
        else:
            _, t1, t2, cen, ax, th = pc
            juncs.append((np.asarray(t2, float), _rot(t, ax, th), _rot(nv, ax, th)))
    # one closed circle edge with one vertex per junction
    verts, circles = [], []
    for c, t, nv in juncs:
        v = w.add('VERTEX_POINT(\'\',%s)' % w.pt(c + r * nv))
        circ = w.add('CIRCLE(\'\',%s,%s)' % (w.ax2(c, t, nv), w.num(r)))
        verts.append(v)
        circles.append(w.add('EDGE_CURVE(\'\',%s,%s,%s,.T.)' % (v, v, circ)))
    faces = []
    # start cap (normal = -t0) and end cap (normal = +t_end)
    c, t, nv = juncs[0]
    loop = w.add('EDGE_LOOP(\'\',(%s))' % w.add('ORIENTED_EDGE(\'\',*,*,%s,.F.)' % circles[0]))
    faces.append(w.add('ADVANCED_FACE(\'\',(%s),%s,.T.)' % (w.add('FACE_OUTER_BOUND(\'\',%s,.T.)' % loop),
                                                         w.add('PLANE(\'\',%s)' % w.ax2(c, -t, nv)))))
    c, t, nv = juncs[-1]
    loop = w.add('EDGE_LOOP(\'\',(%s))' % w.add('ORIENTED_EDGE(\'\',*,*,%s,.T.)' % circles[-1]))
    faces.append(w.add('ADVANCED_FACE(\'\',(%s),%s,.T.)' % (w.add('FACE_OUTER_BOUND(\'\',%s,.T.)' % loop),
                                                         w.add('PLANE(\'\',%s)' % w.ax2(c, t, nv)))))
    # side faces
    for k, pc in enumerate(pieces):
        (c0, t0_, n0), (c1, t1_, n1) = juncs[k], juncs[k + 1]
        s0, s1 = c0 + r * n0, c1 + r * n1
        if pc[0] == 'line':
            L = float(np.linalg.norm(c1 - c0))
            seam_c = w.add('LINE(\'\',%s,%s)' % (w.pt(s0), w.add('VECTOR(\'\',%s,%s)' % (w.dir(t0_), w.num(L)))))
            surf = w.add('CYLINDRICAL_SURFACE(\'\',%s,%s)' % (w.ax2(c0, t0_, n0), w.num(r)))
        else:
            _, t1, t2, cen, ax, th = pc
            cen = np.asarray(cen, float); ax = np.asarray(ax, float)
            rc = float(np.linalg.norm(np.asarray(t1) - cen))
            sc = cen + ax * ((s0 - cen) @ ax)
            seam_c = w.add('CIRCLE(\'\',%s,%s)' % (w.ax2(sc, ax, s0 - sc), w.num(np.linalg.norm(s0 - sc))))
            surf = w.add('TOROIDAL_SURFACE(\'\',%s,%s,%s)' % (w.ax2(cen, ax, np.asarray(t1) - cen), w.num(rc), w.num(r)))
        seam = w.add('EDGE_CURVE(\'\',%s,%s,%s,.T.)' % (verts[k], verts[k + 1], seam_c))
        oe = [w.add('ORIENTED_EDGE(\'\',*,*,%s,.T.)' % circles[k]),
              w.add('ORIENTED_EDGE(\'\',*,*,%s,.T.)' % seam),
              w.add('ORIENTED_EDGE(\'\',*,*,%s,.F.)' % circles[k + 1]),
              w.add('ORIENTED_EDGE(\'\',*,*,%s,.F.)' % seam)]
        loop = w.add('EDGE_LOOP(\'\',(%s))' % ','.join(oe))
        faces.append(w.add('ADVANCED_FACE(\'\',(%s),%s,.T.)' % (w.add('FACE_OUTER_BOUND(\'\',%s,.T.)' % loop), surf)))
    shell = w.add('CLOSED_SHELL(\'\',(%s))' % ','.join(faces))
    return w.add('MANIFOLD_SOLID_BREP(\'%s\',%s)' % (_ascii(name), shell))


def prism_solid(w, poly, z0, z1, name):
    """Planar polygon (x, y in mm) extruded from z0 to z1: flat faces, straight edges."""
    poly = [tuple(map(float, q)) for q in poly]
    area = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(poly, poly[1:] + poly[:1]))
    if area < 0: poly = poly[::-1]                    # counter-clockwise seen from +z
    n = len(poly)
    B = [np.array([x, y, z0]) for x, y in poly]
    T = [np.array([x, y, z1]) for x, y in poly]
    vb = [w.add("VERTEX_POINT('',%s)" % w.pt(p)) for p in B]
    vt = [w.add("VERTEX_POINT('',%s)" % w.pt(p)) for p in T]

    def edge(va, vb_, a, b):
        d = b - a
        ln = w.add("LINE('',%s,%s)" % (w.pt(a), w.add("VECTOR('',%s,%s)" % (w.dir(d), w.num(np.linalg.norm(d))))))
        return w.add("EDGE_CURVE('',%s,%s,%s,.T.)" % (va, vb_, ln))
    eb = [edge(vb[i], vb[(i + 1) % n], B[i], B[(i + 1) % n]) for i in range(n)]
    et = [edge(vt[i], vt[(i + 1) % n], T[i], T[(i + 1) % n]) for i in range(n)]
    ev = [edge(vb[i], vt[i], B[i], T[i]) for i in range(n)]

    def face(oes, o, z, x):
        loop = w.add("EDGE_LOOP('',(%s))" % ','.join(w.add("ORIENTED_EDGE('',*,*,%s,%s)" % (e, s)) for e, s in oes))
        return w.add("ADVANCED_FACE('',(%s),%s,.T.)" % (w.add("FACE_OUTER_BOUND('',%s,.T.)" % loop),
                                                     w.add("PLANE('',%s)" % w.ax2(o, z, x))))
    faces = [face([(et[i], '.T.') for i in range(n)], T[0], [0, 0, 1], [1, 0, 0]),
             face([(eb[i], '.F.') for i in reversed(range(n))], B[0], [0, 0, -1], [1, 0, 0])]
    for i in range(n):
        j = (i + 1) % n
        d = B[j] - B[i]
        out = np.array([d[1], -d[0], 0.0])
        faces.append(face([(eb[i], '.T.'), (ev[j], '.T.'), (et[i], '.F.'), (ev[i], '.F.')], B[i], out, d))
    shell = w.add("CLOSED_SHELL('',(%s))" % ','.join(faces))
    return w.add("MANIFOLD_SOLID_BREP('%s',%s)" % (_ascii(name), shell))


def _ascii(s):
    s = unicodedata.normalize('NFKD', str(s)).encode('ascii', 'ignore').decode('ascii')
    return ''.join(c if c.isalnum() or c in '-_ .' else '_' for c in s)


def _write(path, w, solids_fn, product, desc):
    product = _ascii(product)
    ctx_app = w.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
    w.add("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,%s)" % ctx_app)
    pctx = w.add("PRODUCT_CONTEXT('',%s,'mechanical')" % ctx_app)
    prod = w.add("PRODUCT('%s','%s','',(%s))" % (product, product, pctx))
    w.add("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(%s))" % prod)
    pdf = w.add("PRODUCT_DEFINITION_FORMATION('','',%s)" % prod)
    pdc = w.add("PRODUCT_DEFINITION_CONTEXT('part definition',%s,'design')" % ctx_app)
    pd = w.add("PRODUCT_DEFINITION('design','',%s,%s)" % (pdf, pdc))
    pds = w.add("PRODUCT_DEFINITION_SHAPE('','',%s)" % pd)
    u_len = w.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    u_ang = w.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
    u_sol = w.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
    unc = w.add("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),%s,'distance_accuracy_value','confusion accuracy')" % u_len)
    gctx = w.add("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((%s))"
                 "GLOBAL_UNIT_ASSIGNED_CONTEXT((%s,%s,%s))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))"
                 % (unc, u_len, u_ang, u_sol))
    origin = w.ax2([0, 0, 0], [0, 0, 1], [1, 0, 0])
    solids = solids_fn()
    rep = w.add("ADVANCED_BREP_SHAPE_REPRESENTATION('%s',(%s),%s)" % (product, ','.join([origin] + solids), gctx))
    w.add("SHAPE_DEFINITION_REPRESENTATION(%s,%s)" % (pds, rep))
    stamp = time.strftime('%Y-%m-%dT%H:%M:%S')
    head = ["ISO-10303-21;", "HEADER;", "FILE_DESCRIPTION(('%s'),'2;1');" % desc,
            "FILE_NAME('%s.step','%s',(''),(''),'Anten-O-Matic','Anten-O-Matic','');" % (product, stamp),
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));", "ENDSEC;", "DATA;"]
    with open(path, 'w', encoding='ascii', newline='\n') as f:
        f.write('\n'.join(head + w.lines + ['ENDSEC;', 'END-ISO-10303-21;']) + '\n')


def write_step(path, bodies, product='antenna'):
    """bodies: [(name, pieces_in_mm, radius_mm)]"""
    w = _W()
    _write(path, w, lambda: [tube_solid(w, pieces, r, name) for name, pieces, r in bodies], product,
           'Anten-O-Matic wire antenna')


def write_step_prisms(path, prisms, product='antenna'):
    """prisms: [(name, [(x, y), ...] in mm, z0, z1)] -- e.g. PCB layers of a patch antenna"""
    w = _W()
    _write(path, w, lambda: [prism_solid(w, poly, z0, z1, name) for name, poly, z0, z1 in prisms], product,
           'Anten-O-Matic planar antenna')
