# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""openEMS (FDTD) engine for CAD designs: models dielectrics (e.g. a patch on FR4), which NEC-2 cannot.

The app writes the CSXCAD/openEMS XML itself, runs openEMS.exe from the modules folder, reads the
port voltage/current probes (text) and the frequency-domain near fields on a box around the antenna
(VTK .vtr), and does the near-to-far-field transform here, so no Python bindings or HDF5 are needed.

Geometry rules for openEMS: plates and blocks must be axis aligned (rotations in steps of 90 deg);
wires become thin PEC curves (their radius is ignored); the feed piece must run along x, y or z.
"""
import base64, glob, math, os, re, shutil, struct, subprocess, sys, time, zlib
import xml.etree.ElementTree as ET
import numpy as np

import cad
import pattern as PAT

C0 = 299792458.0
EPS0 = 8.8541878128e-12
ETA0 = 376.730313668


# ------------------------------------------------------------------ finding openEMS

def find_exe(root):
    """openEMS executable: modules/openEMS (the official Windows zip unpacked there), or PATH."""
    names = ['openEMS.exe'] if sys.platform == 'win32' else ['openEMS', 'openEMS.sh']
    base = os.path.join(root, 'modules')
    cands = []
    for n in names:
        cands += glob.glob(os.path.join(base, 'openEMS', n)) + glob.glob(os.path.join(base, '*', n)) + \
                 glob.glob(os.path.join(base, '*', '*', n)) + glob.glob(os.path.join(base, n))
    env = os.environ.get('AOM_OPENEMS')
    if env: cands.insert(0, env)
    for c in cands:
        if os.path.isfile(c): return c
    for n in names:
        p = shutil.which(n)
        if p: return p
    return None


# ------------------------------------------------------------------ geometry -> boxes, curves, port

def _axis_box(e):
    """world-space axis-aligned box (lo, hi) of a plate / block; error if rotated off-axis."""
    R = cad.rotation(e['rot'])
    if not np.allclose(np.abs(R).round(6), np.round(np.abs(R)), atol=1e-6):
        raise ValueError('%s is rotated off the x/y/z axes: openEMS needs plates and blocks aligned with the axes' % e['name'])
    pts = np.array(list(cad.world_points(e).values()))
    return pts.min(0), pts.max(0)


def geometry(doc):
    lam = C0 / ((doc['f0'] + doc['bw'] / 2) * 1e6) * 1000
    boxes, curves, port = [], [], None
    for k, e in enumerate(doc['elements']):
        m = doc['_mats'][e['mat']]
        if e['type'] in ('plate', 'block'):
            lo, hi = _axis_box(e)
            if e['type'] == 'plate' and m['kind'] == 'dielectric':
                raise ValueError('%s: a dielectric must be a block with a thickness' % e['name'])
            if e['type'] == 'plate' and cad.pixel_mask(e) is not None:
                for (x0, y0, x1, y1) in cad.pixel_rects(e):
                    q = np.array([cad.world(e, [x0, y0, 0]), cad.world(e, [x1, y1, 0])])
                    boxes.append({'el': k, 'lo': q.min(0), 'hi': q.max(0), 'mat': m, 'name': e['name']})
                continue
            boxes.append({'el': k, 'lo': lo, 'hi': hi, 'mat': m, 'name': e['name']})
            continue
        if e['type'] == 'pcb':
            b_, c_, p_ = cad.PCB.ems_parts(doc, k, cad.world, cad.rotation)
            boxes += b_; curves += c_
            if p_ is not None:
                if port is not None: raise ValueError('openEMS engine: only one feed is supported')
                port = p_
            continue
        if m['kind'] != 'conductor': raise ValueError('%s: wires must be conductors' % e['name'])
        paths = [cad.world(e, Q) for Q in cad.local_paths(e)]
        fp = min(max(0, e.get('feed_path', 0)), len(paths) - 1)
        for pi, P in enumerate(paths):
            if not (e['feed'] and pi == fp): curves.append({'el': k, 'pts': P})
        P = paths[fp]
        if e['feed']:
            if port is not None: raise ValueError('openEMS engine: only one feed is supported')
            seg = np.diff(P, axis=0); L = np.linalg.norm(seg, axis=1)
            straight = len(P) == 2
            ax = int(np.argmax(np.abs(seg[0]))) if straight else None
            if straight and abs(seg[0][ax]) / L[0] > 0.9999 and L[0] < lam / 6:
                port = {'a': P[0].copy(), 'b': P[1].copy(), 'ax': ax, 'el': k}        # a probe: the whole rod is the port
                continue
            P2, i = cad._split_at(P, e['feed_at'])
            j = min(i, len(P2) - 2) if i < len(P2) - 1 else i - 1
            d = P2[j + 1] - P2[j]
            ax = int(np.argmax(np.abs(d)))
            if abs(d[ax]) / np.linalg.norm(d) < 0.9999:
                raise ValueError('%s: openEMS needs the wire at the feed point to run exactly along x, y or z' % e['name'])
            g = min(lam / 60, 0.2 * float(np.linalg.norm(d)) if 0 < i < len(P2) - 1 else 0.5 * float(np.linalg.norm(d)))
            c = P2[i]; u = np.zeros(3); u[ax] = np.sign(d[ax])
            if i == 0: a, b = c, c + u * g; rest = [np.vstack([b, P2[1:]])]
            elif i == len(P2) - 1: a, b = c - u * g, c; rest = [np.vstack([P2[:-1], a])]
            else: a, b = c - u * g / 2, c + u * g / 2; rest = [np.vstack([P2[:i], a]), np.vstack([b, P2[i + 1:]])]
            port = {'a': a, 'b': b, 'ax': ax, 'el': k}
            curves += [{'el': k, 'pts': Q} for Q in rest if len(Q) > 1]
    if port is None: raise ValueError('no feed: tick "Feed" on one wire element')
    return boxes, curves, port


# ------------------------------------------------------------------ mesh

def _q(v):
    v = round(float(v), 6)
    return 0.0 if v == 0 else v


def smooth_lines(fixed, res_at, ratio=1.4, hard=(), merge=1e-3):
    """rectilinear mesh lines: keep every fixed line, fill each gap with cells that grow
    geometrically (by at most `ratio`) from the small features at its ends up to the local resolution.
    Fixed lines closer than `merge` become one (a hard line, e.g. the port, wins): two edges a few
    micrometres apart would otherwise force a tiny cell and a timestep thousands of times shorter."""
    hard = set(_q(v) for v in hard)
    L = sorted(set(_q(v) for v in fixed))
    merged = [L[0]]
    for v in L[1:]:
        if v - merged[-1] > merge: merged.append(v)
        elif v in hard and merged[-1] not in hard: merged[-1] = v
    L = merged
    if len(L) < 2: return np.array(L)
    gaps = np.diff(L)
    size = [min(gaps[max(0, i - 1)], gaps[min(len(gaps) - 1, i)], res_at(L[i])) for i in range(len(L))]
    out = [L[0]]
    for i, (a, b) in enumerate(zip(L[:-1], L[1:])):
        r = min(res_at(a), res_at(b), res_at((a + b) / 2))
        gap = b - a
        sa, sb = min(size[i], r), min(size[i + 1], r)
        if gap <= r * 1.0001 and gap <= min(sa, sb) * ratio:
            out.append(b); continue
        left, right, x, y, tot = [], [], sa, sb, 0.0
        while tot < gap:
            if x <= y: left.append(x); tot += x; x = min(x * ratio, r)
            else: right.append(y); tot += y; y = min(y * ratio, r)
        cells = np.array(left + right[::-1]) * gap / tot
        out += list(a + np.cumsum(cells)[:-1]) + [b]     # land exactly on the fixed line
    out[-1] = L[-1]
    # whole nanometres: a port or probe box of zero width must sit exactly on its line, older
    # openEMS builds (e.g. the v0.37 Windows release) do not snap 1e-15 to 0 and drop the port
    return np.array(sorted(set(_q(v) for v in out)))


# mesh density presets: air cell = lambda/a, cell inside the structure = lambda/b,
# cell in a dielectric = (air cell / sqrt(eps)) / c
MESH_Q = {'draft': (10, 16, 1.0), 'coarse': (15, 25, 1.0), 'normal': (20, 40, 2.0), 'fine': (28, 60, 3.0), 'vfine': (40, 80, 4.0)}


def mesh_q(doc):
    """(air, structure, dielectric) density of the design's mesh setting; 'custom' takes the user's lambda/N values."""
    m = doc.get('ems_mesh') or 'normal'
    if m == 'custom':
        a = min(80.0, max(6.0, float(doc.get('ems_air') or 20))); f = min(160.0, max(a, float(doc.get('ems_fine') or 40)))
        return a, f, max(1.0, f / a)
    return MESH_Q.get(m, MESH_Q['normal'])


def mesh_info(doc):
    """mesh size of a design without running anything (for the UI estimate)."""
    boxes, curves, port = geometry(doc)
    lines, _, _ = build_mesh(doc, boxes, curves, port)
    n = [len(L) for L in lines]
    cells = (n[0] - 1) * (n[1] - 1) * (n[2] - 1)
    return {'lines': n, 'cells': cells, 'min_cell': min(float(np.diff(L).min()) for L in lines),
            'ram_mb': cells * 150 / 1e6}


def build_mesh(doc, boxes, curves, port):
    fmax = (doc['f0'] + doc['bw'] / 2) * 1e6
    fmin = max(1e6, (doc['f0'] - max(doc['bw'], 0.3 * doc['f0']) / 2) * 1e6)
    lam_min = C0 / fmax * 1000
    lam_max = C0 / fmin * 1000
    q_air, q_fine, q_diel = mesh_q(doc)
    res_air = lam_min / q_air
    res_fine = lam_min / q_fine        # inside the structure: metal edges and resonant lengths need it
    pts = [port['a'], port['b']]
    for b in boxes: pts += [b['lo'], b['hi']]
    for c in curves: pts += list(c['pts'])
    pts = np.array(pts)
    lo, hi = pts.min(0), pts.max(0)
    margin = lam_max * float(doc.get('ems_margin') or 0.25) + 8 * res_air   # air around the antenna + 8 PML cells
    ground = doc['ground'] == 'perfect'
    sim_lo = lo - margin; sim_hi = hi + margin
    if ground: sim_lo[2] = 0.0
    nf_off = max(0.55 * (margin - 10 * res_air), 3 * res_fine)   # near-field box: in the air, clear of the absorbing layer
    nf_lo = lo - nf_off; nf_hi = hi + nf_off
    if ground: nf_lo[2] = 0.0
    lines = []
    for ax in range(3):
        fixed = [sim_lo[ax], sim_hi[ax], nf_lo[ax], nf_hi[ax], port['a'][ax], port['b'][ax]]
        diel = []
        for b in boxes:
            fixed += [b['lo'][ax], b['hi'][ax]]
            t = b['hi'][ax] - b['lo'][ax]
            if b['mat']['kind'] == 'dielectric':
                diel.append((b['lo'][ax], b['hi'][ax], res_air / math.sqrt(b['mat']['eps'])))
                if 0 < t < res_air:   # thin substrate: at least 4 cells through it
                    fixed += list(np.linspace(b['lo'][ax], b['hi'][ax], 5))
        for c in curves: fixed += list(c['pts'][:, ax])
        for b in boxes: fixed += [q[ax] for q in b.get('fix', [])]

        for b in boxes:                 # thirds rule at metal plate edges (openEMS recommendation)
            if b['mat']['kind'] == 'conductor' and b['hi'][ax] - b['lo'][ax] > 1e-6 and (b.get('poly') or {}).get('kind') != 'hole':
                h = min(res_fine, (b['hi'][ax] - b['lo'][ax]) / 6)
                fixed += [b['lo'][ax] - h / 3, b['lo'][ax] + 2 * h / 3, b['hi'][ax] - 2 * h / 3, b['hi'][ax] + h / 3]
        slo, shi = lo[ax], hi[ax]

        def res_at(x, diel=diel, slo=slo, shi=shi):
            r = res_air if not (slo - 1e-6 <= x <= shi + 1e-6) else res_fine
            for a, b, rd in diel:
                if a - 1e-6 <= x <= b + 1e-6: r = min(r, rd / q_diel)
            return r
        lines.append(smooth_lines(fixed, res_at, hard=[port['a'][ax], port['b'][ax]], merge=max(1e-3, min(res_fine / 12, 0.1))))
    return lines, (nf_lo, nf_hi), (sim_lo, sim_hi)


# ------------------------------------------------------------------ XML

def _f(v): return '%.9g' % v


def _c(v): return '%.9g' % _q(v)       # a coordinate, on the same nanometre grid as the mesh lines


def _prop(parent, tag, name, pid, **attrs):
    e = ET.SubElement(parent, tag, dict({'ID': str(pid), 'Name': name}, **{k: str(v) for k, v in attrs.items()}))
    return e


def _box(prims, lo, hi, prio):
    b = ET.SubElement(prims, 'Box', {'Priority': str(prio)})
    ET.SubElement(b, 'P1', {'X': _c(lo[0]), 'Y': _c(lo[1]), 'Z': _c(lo[2])})
    ET.SubElement(b, 'P2', {'X': _c(hi[0]), 'Y': _c(hi[1]), 'Z': _c(hi[2])})


def build_xml(doc, path, ts=400000, end=1e-4):
    boxes, curves, port = geometry(doc)
    lines, (nf_lo, nf_hi), _ = build_mesh(doc, boxes, curves, port)
    f0 = doc['f0'] * 1e6
    fc = max(doc['bw'] * 1e6, 0.5 * f0) if doc['bw'] else 0.5 * f0
    fc = min(fc, 0.9 * f0)
    root = ET.Element('openEMS')
    fd = ET.SubElement(root, 'FDTD', {'NumberOfTimesteps': str(ts), 'endCriteria': _f(end), 'f_max': _f(f0 + fc)})
    ET.SubElement(fd, 'Excitation', {'Type': '0', 'f0': _f(f0), 'fc': _f(fc)})
    bc = 'PML_8'
    ET.SubElement(fd, 'BoundaryCond', {'xmin': bc, 'xmax': bc, 'ymin': bc, 'ymax': bc,
                                       'zmin': 'PEC' if doc['ground'] == 'perfect' else bc, 'zmax': bc})
    cs = ET.SubElement(root, 'ContinuousStructure', {'CoordSystem': '0'})
    grid = ET.SubElement(cs, 'RectilinearGrid', {'DeltaUnit': '0.001', 'CoordSystem': '0'})
    for tag, L in zip(('XLines', 'YLines', 'ZLines'), lines):
        ET.SubElement(grid, tag, {'Qty': str(len(L))}).text = ','.join(_c(v) for v in L)
    ET.SubElement(cs, 'BackgroundMaterial', {'Epsilon': '1', 'Mue': '1', 'Kappa': '0', 'Sigma': '0'})
    ET.SubElement(cs, 'ParameterSet')
    props = ET.SubElement(cs, 'Properties')
    pid = [0]

    def nid():
        pid[0] += 1; return pid[0] - 1
    # dielectrics (kappa from the loss tangent at f0) and metals
    for i, b in enumerate(boxes):
        m = b['mat']
        if b.get('poly'):            # planar layer: a polygon on its plane (holes: air with a higher priority)
            q = b['poly']
            if q['kind'] == 'hole':
                p = _prop(props, 'Material', 'hole_%d_%s' % (i, _safe(b['name'])), nid(), Isotropy=1)
                ET.SubElement(p, 'Property', {'Epsilon': '1', 'Mue': '1', 'Kappa': '0', 'Sigma': '0'})
            else:
                p = _prop(props, 'Metal', 'metal_%d_%s' % (i, _safe(b['name'])), nid())
            pg = ET.SubElement(ET.SubElement(p, 'Primitives'), 'Polygon', {'Priority': str(q['prio']), 'Elevation': _c(q['elev']), 'NormDir': str(q['nd']), 'CoordSystem': '0'})
            for (u, v) in q['v']: ET.SubElement(pg, 'Vertex', {'X1': _c(u), 'X2': _c(v)})
            continue
        if m['kind'] == 'dielectric':
            p = _prop(props, 'Material', 'diel_%d_%s' % (i, _safe(b['name'])), nid(), Isotropy=1)
            kappa = 2 * math.pi * f0 * EPS0 * m['eps'] * m.get('tand', 0)
            ET.SubElement(p, 'Property', {'Epsilon': _f(m['eps']), 'Mue': '1', 'Kappa': _f(kappa), 'Sigma': '0'})
            _box(ET.SubElement(p, 'Primitives'), b['lo'], b['hi'], 1)
        else:
            p = _prop(props, 'Metal', 'metal_%d_%s' % (i, _safe(b['name'])), nid())
            _box(ET.SubElement(p, 'Primitives'), b['lo'], b['hi'], 10)
    if curves:
        p = _prop(props, 'Metal', 'wires', nid())
        prims = ET.SubElement(p, 'Primitives')
        for c in curves:
            cu = ET.SubElement(prims, 'Curve', {'Priority': '10'})
            for q in c['pts']: ET.SubElement(cu, 'Vertex', {'X': _c(q[0]), 'Y': _c(q[1]), 'Z': _c(q[2])})
    # lumped port 1 (50 ohm resistor + soft E-field excitation + U/I probes), as openEMS' LumpedPort
    a, b, ax = np.array(port['a'], float), np.array(port['b'], float), port['ax']
    if a[ax] > b[ax]: a, b = b, a
    direction = 1.0
    p = _prop(props, 'LumpedElement', 'port_resist_1', nid(), Direction=ax, Caps=1, R=_f(doc['z0']), LEtype='0')
    _box(ET.SubElement(p, 'Primitives'), a, b, 5)
    exc = np.zeros(3); exc[ax] = -1 * direction
    p = _prop(props, 'Excitation', 'port_excite_1', nid(), Number=0, Enabled=1, Frequency='0', Delay='0', Type=0,
              Excite=','.join(_f(v) for v in exc), PropDir='0,0,0')
    _box(ET.SubElement(p, 'Primitives'), a, b, 5)
    ET.SubElement(p, 'Weight', {'X': '1', 'Y': '1', 'Z': '1'})
    us, ue = (a + b) / 2, (a + b) / 2; us[ax], ue[ax] = a[ax], b[ax]
    p = _prop(props, 'ProbeBox', 'port_ut_1', nid(), Number=0, Type=0, Weight=-1, NormDir=-1, StartTime=0, StopTime=0)
    _box(ET.SubElement(p, 'Primitives'), us, ue, 0)
    i_s, i_e = a.copy(), b.copy(); i_s[ax] = i_e[ax] = (a[ax] + b[ax]) / 2
    p = _prop(props, 'ProbeBox', 'port_it_1', nid(), Number=0, Type=1, Weight=_f(direction), NormDir=ax, StartTime=0, StopTime=0)
    _box(ET.SubElement(p, 'Primitives'), i_s, i_e, 0)
    # near-field box: frequency-domain E and H on its faces (VTK files)
    faces = []
    for k in range(3):
        for side, v in (('n', nf_lo[k]), ('p', nf_hi[k])):
            if doc['ground'] == 'perfect' and k == 2 and side == 'n': continue     # the ground plane closes the box
            lo_, hi_ = nf_lo.copy(), nf_hi.copy(); lo_[k] = hi_[k] = v
            faces.append(('xyz'[k] + side, lo_, hi_))
    for fld, dt in (('E', 10), ('H', 11)):
        for nm, lo_, hi_ in faces:
            p = _prop(props, 'DumpBox', 'nf2ff_%s_%s' % (fld, nm), nid(), Number=0, Type=0, Weight=1, NormDir=-1, StartTime=0, StopTime=0,
                      DumpType=dt, DumpMode=1, FileType=0, MultiGridLevel=0)
            ET.SubElement(p, 'FD_Samples').text = _f(f0)
            _box(ET.SubElement(p, 'Primitives'), lo_, hi_, 0)
    ET.indent(root) if hasattr(ET, 'indent') else None
    with open(path, 'wb') as f:
        f.write(b'<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n')
        f.write(ET.tostring(root))
    ncells = (len(lines[0]) - 1) * (len(lines[1]) - 1) * (len(lines[2]) - 1)
    return {'lines': lines, 'faces': [f[0] for f in faces], 'cells': ncells, 'port': port, 'fc': fc}


def _safe(s):
    return ''.join(c if c.isalnum() else '_' for c in s)[:30]


# ------------------------------------------------------------------ results

def read_probe(path):
    t, v = [], []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith('%'): continue
            a = line.split()
            t.append(float(a[0])); v.append(float(a[1]))
    return np.array(t), np.array(v)


def dft(t, v, freqs):
    dt = np.gradient(t) if len(t) > 1 else np.ones(1)
    return np.array([2 * np.sum(v * np.exp(-2j * np.pi * f * t) * dt) for f in freqs])


def read_vtr(path):
    """coordinates (x, y, z) and point data {name: array} of a VTK XML RectilinearGrid file
    (inline ascii, inline base64 binary with or without zlib, or appended raw / base64 data)."""
    raw = open(path, 'rb').read()
    at = raw.find(b'<AppendedData')
    app = None
    if at >= 0:
        enc = _enc(raw, at)
        start = raw.find(b'_', at) + 1
        stop_ = raw.rfind(b'</AppendedData>')
        app = (raw[start:stop_], enc)
        raw = raw[:at] + b'</VTKFile>'
    root = ET.fromstring(raw)
    hdr64 = root.get('header_type', 'UInt32') == 'UInt64'
    comp = root.get('compressor')
    bo = '<' if root.get('byte_order', 'LittleEndian') == 'LittleEndian' else '>'
    dt = {'Float32': 'f4', 'Float64': 'f8', 'Int32': 'i4', 'Int64': 'i8', 'UInt8': 'u1', 'UInt32': 'u4', 'UInt64': 'u8'}

    def arr(da):
        typ = np.dtype(bo + dt[da.get('type')])
        fmt = da.get('format', 'ascii')
        if fmt == 'ascii': return np.array(da.text.split(), dtype=typ)
        if fmt == 'binary': return _decode_b64(da.text.strip(), typ, hdr64, comp, bo)
        data, enc = app
        off = int(da.get('offset'))
        if enc == 'base64':
            txt = data.decode('latin1')[off:]
            return _decode_b64(txt.split()[0] if txt.strip() else txt, typ, hdr64, comp, bo, stream=True)
        return _decode_raw(data[off:], typ, hdr64, comp, bo)
    grid = root.find('RectilinearGrid')
    piece = grid.find('Piece')
    coords = [arr(da) for da in piece.find('Coordinates').findall('DataArray')]
    fields = {}
    for da in piece.find('PointData').findall('DataArray'):
        a = arr(da)
        n = int(da.get('NumberOfComponents', 1))
        fields[da.get('Name')] = a.reshape(-1, n) if n > 1 else a
    return coords, fields


def _enc(raw, at):
    m = re.search(rb'encoding="(\w+)"', raw[at:at + 200])
    return m.group(1).decode() if m else 'raw'


def _decode_b64(text, typ, hdr64, comp, bo, stream=False):
    hs = 8 if hdr64 else 4
    hc = bo + ('Q' if hdr64 else 'I')
    if comp:
        nb = struct.unpack(hc, base64.b64decode(text[:12])[:hs])[0]
        hlen = (3 + nb) * hs
        hb64 = 4 * math.ceil(hlen / 3)
        vals = struct.unpack(bo + hc[1] * (3 + nb), base64.b64decode(text[:hb64])[:hlen])
        csum = sum(vals[3:])
        body = base64.b64decode(text[hb64:hb64 + 4 * math.ceil(csum / 3)])
        out, pos = [], 0
        for sz in vals[3:]:
            out.append(zlib.decompress(body[pos:pos + sz])); pos += sz
        return np.frombuffer(b''.join(out), dtype=typ)
    n = struct.unpack(hc, base64.b64decode(text[:12])[:hs])[0]
    data = base64.b64decode(text[:4 * math.ceil((hs + n) / 3)])
    return np.frombuffer(data[hs:hs + n], dtype=typ)


def _decode_raw(data, typ, hdr64, comp, bo):
    hs = 8 if hdr64 else 4
    hc = bo + ('Q' if hdr64 else 'I')
    if comp:
        nb = struct.unpack(hc, data[:hs])[0]
        vals = struct.unpack(bo + hc[1] * (3 + nb), data[:(3 + nb) * hs])
        pos = (3 + nb) * hs
        out = []
        for sz in vals[3:]:
            out.append(zlib.decompress(data[pos:pos + sz])); pos += sz
        return np.frombuffer(b''.join(out), dtype=typ)
    n = struct.unpack(hc, data[:hs])[0]
    return np.frombuffer(data[hs:hs + n], dtype=typ)


def read_fd_dump(simdir, name, f0):
    """complex field (N, 3) and coordinates of one frequency-domain dump."""
    absf = sorted(glob.glob(os.path.join(simdir, name + '_f=*_abs.vtr')))
    argf = sorted(glob.glob(os.path.join(simdir, name + '_f=*_arg.vtr')))
    if not absf or not argf: raise RuntimeError('openEMS did not write the near-field dump %s' % name)
    c, fa = read_vtr(absf[0]); _, fp = read_vtr(argf[0])
    A = list(fa.values())[0].astype(float); P = list(fp.values())[0].astype(float)
    return c, A * np.exp(1j * P)


def nf2ff(simdir, faces, f0, step=5, ground=False):
    """far field from the E/H dumps on the near-field box (equivalence principle): returns the
    radiated power [W] and E_theta, E_phi * r [V] on the pattern grid."""
    k = 2 * math.pi * f0 / C0
    th, ph = PAT.grid_axes(step, ground)
    T, Pp = np.meshgrid(np.radians(th), np.radians(ph), indexing='ij')
    rhat = np.stack([np.sin(T) * np.cos(Pp), np.sin(T) * np.sin(Pp), np.cos(T)], -1)
    that = np.stack([np.cos(T) * np.cos(Pp), np.cos(T) * np.sin(Pp), -np.sin(T)], -1)
    phat = np.stack([-np.sin(Pp), np.cos(Pp), np.zeros_like(T)], -1)
    N = np.zeros(T.shape + (3,), complex); Lm = np.zeros_like(N)
    prad = 0.0
    for nm in faces:
        ce, E = read_fd_dump(simdir, 'nf2ff_E_' + nm, f0)
        ch, H = read_fd_dump(simdir, 'nf2ff_H_' + nm, f0)
        x, y, z = [np.asarray(v, float) for v in ce]     # the VTK coordinates are in metres
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        shape = (len(x), len(y), len(z))
        E = E.reshape(len(z), len(y), len(x), 3).transpose(2, 1, 0, 3).reshape(-1, 3)
        H = H.reshape(len(z), len(y), len(x), 3).transpose(2, 1, 0, 3).reshape(-1, 3)
        ax = 'xyz'.index(nm[0]); sgn = 1.0 if nm[1] == 'p' else -1.0
        n = np.zeros(3); n[ax] = sgn
        tang = [i for i in range(3) if i != ax]
        ws = []
        for i in tang:
            v = [x, y, z][i]
            w = np.zeros_like(v)
            if len(v) > 1:
                w[1:-1] = (v[2:] - v[:-2]) / 2; w[0] = (v[1] - v[0]) / 2; w[-1] = (v[-1] - v[-2]) / 2
            ws.append(w)
        dA = np.multiply.outer(ws[0], ws[1])
        full = [1, 1, 1]; full[tang[0]] = len(ws[0]); full[tang[1]] = len(ws[1])
        dA = dA.reshape(full).repeat(shape[ax], axis=ax).reshape(-1)
        if ground: dA = dA * np.where(Z.reshape(-1) < 0, 0, 1)
        J = np.cross(n, H); M = -np.cross(n, E)
        prad += 0.5 * float(np.real(np.sum((np.cross(E, np.conj(H)) @ n) * dA)))
        R = np.stack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)], -1)
        ph_ = np.exp(1j * k * (rhat.reshape(-1, 3) @ R.T))          # (ndir, npts)
        N += (ph_ @ (J * dA[:, None])).reshape(N.shape)
        Lm += (ph_ @ (M * dA[:, None])).reshape(N.shape)
        if ground:    # image currents of the PEC ground (z = 0): J_t -> -J_t, J_z -> J_z ; M_t -> M_t, M_z -> -M_z
            Ri = R * np.array([1, 1, -1])
            phi_ = np.exp(1j * k * (rhat.reshape(-1, 3) @ Ri.T))
            Ji = J * np.array([-1, -1, 1]); Mi = M * np.array([1, 1, -1])
            N += (phi_ @ (Ji * dA[:, None])).reshape(N.shape)
            Lm += (phi_ @ (Mi * dA[:, None])).reshape(N.shape)
    Nt, Np = (N * that).sum(-1), (N * phat).sum(-1)
    Lt, Lp = (Lm * that).sum(-1), (Lm * phat).sum(-1)
    Eth = -1j * k / (4 * math.pi) * (Lp + ETA0 * Nt)
    Eph = 1j * k / (4 * math.pi) * (Lt - ETA0 * Np)
    if ground: prad *= 1.0   # the box is open at z = 0 (PEC), its flux already is the half-space total
    return th, ph, Eth, Eph, prad


def pattern_from_ff(th, ph, Eth, Eph, p_in, ground):
    U = lambda e: np.abs(e) ** 2 / (2 * ETA0)            # radiation intensity per polarisation [W/sr]
    G = {}
    for pol, (e1) in (('total', None), ('rhcp', (Eth + 1j * Eph) / math.sqrt(2)), ('lhcp', (Eth - 1j * Eph) / math.sqrt(2))):
        u = U(Eth) + U(Eph) if pol == 'total' else U(e1)
        G[pol] = np.maximum(10 * np.log10(np.maximum(4 * math.pi * u / p_in, 1e-30)), PAT.FLOOR)
    return {'step': float(th[1] - th[0]), 'th': th.tolist(), 'ph': ph.tolist(), 'ground': bool(ground), 'approx': False,
            'total': G['total'].round(2).tolist(), 'rhcp': G['rhcp'].round(2).tolist(), 'lhcp': G['lhcp'].round(2).tolist(),
            'note': 'openEMS (FDTD) result, gain relative to the power accepted at the port'}


# ------------------------------------------------------------------ run

def _kcells(n): return '%.2f M' % (n / 1e6) if n >= 1e6 else '%.0f k' % (n / 1e3)


def _mmss(t):
    t = int(max(0, t)); return '%d:%02d' % (t // 60, t % 60) if t < 3600 else '%d:%02d:%02d' % (t // 3600, t // 60 % 60, t % 60)


class _NoField(Exception):
    pass


def _run(exe, xml, simdir, engine, threads, log, stop, status, end_db):
    """one openEMS run; every output line goes to openems.log. Raises _NoField when the
    energy is still exactly zero after a few reports (the excitation did not start)."""
    cmd = [exe, os.path.basename(xml), '--engine=' + engine] + (['--numThreads=%d' % threads] if threads else [])
    kw = {'creationflags': 0x08000000} if sys.platform == 'win32' else {}
    log('openEMS: building the FDTD operator (%s engine; no output for a while is normal here)' % engine)
    if status: status('openEMS: building the operator', 0)
    t0 = time.time()
    lf = open(os.path.join(simdir, 'openems.log'), 'a', encoding='utf-8', errors='replace')
    lf.write('> %s\n' % ' '.join(cmd))
    pr = subprocess.Popen(cmd, cwd=simdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace', bufsize=1, **kw)
    tail = []; tstep = None; best = 0.0; emax = 0.0; reports = 0; exc = ''; exc_len = 0
    try:
        for line in pr.stdout:
            line = line.rstrip()
            lf.write(line + '\n')
            if line: tail.append(line); del tail[:-30]
            if 'excitations' in line.lower() and ':' in line: exc += (' ' if exc else '') + ' '.join(line.split())
            me = re.search(r'Excitation signal length is:\s*(\d+)', line)
            if me: exc_len = int(me.group(1))
            if stop is not None and stop.is_set():
                pr.terminate(); raise RuntimeError('stopped')
            m = re.search(r'Timestep:\s*(\d+)(?:.*?Speed:\s*([\d.]+)\s*MC/s)?.*?Energy:\s*~?\s*([\d.e+-]+)\s*\(\s*-?\s*([\d.]+)\s*dB', line)
            if not m: continue
            now = time.time(); reports += 1
            if tstep is None:
                tstep = now; log('openEMS: operator ready after %s, time stepping' % _mmss(now - t0))
            try: emax = max(emax, abs(float(m.group(3))))
            except ValueError: pass
            if reports >= 3 and emax == 0.0 and int(m.group(1)) > exc_len:   # the pulse is over and nothing moved
                pr.terminate(); raise _NoField(exc or 'energy 0 after %s timesteps' % m.group(1))
            db = float(m.group(4)); best = max(best, db)
            pct = min(99.0, 100.0 * best / end_db)
            run = now - tstep
            eta = ' · ~%s left' % _mmss(run * (end_db - best) / best) if best > 3 and run > 8 else ''
            spd = ' · %s MC/s' % m.group(2) if m.group(2) else ''
            msg = 'openEMS %2.0f %% · energy -%.1f dB of -%.0f dB · step %s%s · %s%s' % (pct, db, end_db, m.group(1), spd, _mmss(now - t0), eta)
            log(msg)
            if status: status(msg, pct)
    finally:
        pr.wait(); lf.close()
    if pr.returncode != 0:
        raise RuntimeError('openEMS failed (exit code %s): %s' % (pr.returncode, ' | '.join(tail[-6:])))
    if emax == 0.0 and reports:
        raise _NoField(exc or 'energy 0')


def simulate(doc, exe, simdir, log=print, stop=None, threads=0, status=None, end_db=40.0):
    """runs openEMS; `log` gets a progress line at every openEMS report (about every 4 s),
    `status(msg, percent)` (optional) feeds the status bar."""
    os.makedirs(simdir, exist_ok=True)
    for f in glob.glob(os.path.join(simdir, '*')):
        if os.path.isfile(f): os.remove(f)
    xml = os.path.join(simdir, 'antenna_openems.xml')
    info = build_xml(doc, xml)
    log('openEMS: mesh %d x %d x %d = %s cells (%s mesh), about %.0f MB RAM' % (*(len(L) for L in info['lines']), _kcells(info['cells']), doc.get('ems_mesh') or 'normal', info['cells'] * 150 / 1e6))
    if info['cells'] > 30e6: raise ValueError('the mesh has %.0f M cells: too big (use the coarse mesh, make the design smaller or the frequency lower)' % (info['cells'] / 1e6))
    engines = ['multithreaded', 'basic']
    for k, eng in enumerate(engines):
        try:
            _run(exe, xml, simdir, eng, threads, log, stop, status, end_db)
            break
        except _NoField as e:
            if k + 1 < len(engines):
                log('openEMS: the field stays at zero with the %s engine (%s); trying the %s engine' % (eng, e, engines[k + 1]))
                for f in glob.glob(os.path.join(simdir, '*')):
                    if os.path.isfile(f) and not f.endswith(('.xml', '.log')): os.remove(f)
            else:
                raise RuntimeError('openEMS: the port excitation never reached the field (%s). The full openEMS output is in %s'
                                   % (e, os.path.join(simdir, 'openems.log')))
    return info


def results(doc, simdir, info, freqs):
    t, u = read_probe(os.path.join(simdir, 'port_ut_1'))
    ti, i = read_probe(os.path.join(simdir, 'port_it_1'))
    F = np.array(freqs) * 1e6
    U, I = dft(t, u, F), dft(ti, i, F)
    Z = U / I
    f0 = doc['f0'] * 1e6
    U0, I0 = dft(t, u, [f0])[0], dft(ti, i, [f0])[0]
    p_in = 0.5 * float(np.real(U0 * np.conj(I0)))
    ground = doc['ground'] == 'perfect'
    th, ph, Eth, Eph, prad = nf2ff(simdir, info['faces'], f0, 5, ground)
    if p_in <= 0: raise RuntimeError('openEMS: no power accepted at the port (check the feed)')
    pat = pattern_from_ff(th, ph, Eth, Eph, p_in, ground)
    eff = min(1.0, prad / p_in) if p_in > 0 else 0
    return Z, pat, eff, prad, p_in


def save(doc, out_root, exe, log=print, stop=None, status=None):
    """openEMS counterpart of cad.save: same outputs (render, report, STEP, pattern, design)."""
    import json
    from core import swr
    stamp = time.strftime('%Y%m%d-%H%M%S')
    safe = ''.join(c if c.isalnum() or c in '-_' else '_' for c in doc['name'])[:40] or 'cad'
    out = os.path.join(out_root, '%s_%s_%gMHz_openems' % (stamp, safe, doc['f0']))
    sim = os.path.join(out, 'openems')
    os.makedirs(sim, exist_ok=True)
    info = simulate(doc, exe, sim, log, stop, status=status)
    log('openEMS finished, computing impedance and the far field')
    if status: status('openEMS finished, computing impedance and the far field', 100)
    bw = doc['bw'] or 0.1 * doc['f0']
    freqs = list(np.linspace(doc['f0'] - bw / 2, doc['f0'] + bw / 2, 5))
    Z, pat, eff, prad, p_in = results(doc, sim, info, freqs)
    wide = list(np.linspace(doc['f0'] - info['fc'] / 1e6 * 0.9, doc['f0'] + info['fc'] / 1e6 * 0.9, 41))
    Zw = results_z(sim, wide)
    for f in glob.glob(os.path.join(sim, '*.vtr')): os.remove(f)     # near-field dumps: big and no longer needed
    G = np.array(pat['total'])
    gmax = float(G.max())
    pr = lambda key: float(np.array(pat[key]).max())
    rows = [{'f': f, 'zr': z.real, 'zi': z.imag, 'swr': swr(z, doc['z0']), 'gmax': gmax if abs(f - doc['f0']) < 1e-9 else None,
             'rhcp': None, 'lhcp': None} for f, z in zip(freqs, Z)]
    mid = len(rows) // 2
    rows[mid].update(gmax=gmax, rhcp=pr('rhcp'), lhcp=pr('lhcp'))
    with open(os.path.join(out, 'pattern.json'), 'w') as f: json.dump(pat, f)
    cad_rows = [dict(r, gmax=r['gmax'] if r['gmax'] is not None else float('nan')) for r in rows]
    render(doc, os.path.join(out, 'render.png'), pat, Zw, wide)
    if doc.get('step', True):
        try: cad.export_step(doc, os.path.join(out, 'antenna.step'))
        except Exception as e: log('STEP export failed: %s' % e)
    if doc.get('svg', True): cad.PCB.write_svgs(doc, out)
    btxt = cad.build_text(doc)
    best = min(zip(wide, Zw), key=lambda q: abs((q[1] - doc['z0']) / (q[1] + doc['z0'])))
    rep = ['Anten-O-Matic CAD (openEMS FDTD): %s, %.3f MHz' % (doc['name'], doc['f0']), '=' * 40,
           'λ = %.1f mm, port impedance %.0f Ω, mesh %d x %d x %d cells' % (doc['_lam'], doc['z0'], *(len(L) - 1 for L in info['lines'])),
           'radiation efficiency %.0f %% (radiated %.3g W of %.3g W accepted)' % (eff * 100, prad, p_in),
           'best match in the simulated span: %.2f MHz, SWR %.2f' % (best[0], swr(best[1], doc['z0'])), '',
           'MHz         Z [Ω]            SWR    max [dBi]   RHCP [dBic]   LHCP [dBic]']
    for r in rows:
        g = lambda v: ('%6.2f' % v) if v is not None else '     –'
        rep.append('%9.3f  %7.1f %+7.1fj   %5.2f   %s      %s        %s' % (r['f'], r['zr'], r['zi'], r['swr'], g(r['gmax']), g(r['rhcp']), g(r['lhcp'])))
    rep += ['', '(gain and pattern at the centre frequency)', '', 'CONSTRUCTION', '-' * 40, btxt,
            '', 'openEMS input and probe files are in the openems subfolder (antenna_openems.xml can be opened in AppCSXCAD).']
    text = '\n'.join(rep) + '\n'
    with open(os.path.join(out, 'report.txt'), 'w', encoding='utf-8') as f: f.write(text)
    keep = {k: v for k, v in doc.items() if not k.startswith('_')}
    with open(os.path.join(out, 'design.json'), 'w', encoding='utf-8') as f:
        json.dump({'type': 'cad', 'engine': 'openems', 'doc': keep, 'materials': list(doc['_mats'].values()), 'summary': rows}, f, indent=1, ensure_ascii=False, default=float)
    summ = [{'f': r['f'], 'floor': r['gmax'], 'mean': None, 'max': r['gmax'], 'zr': r['zr'], 'zi': r['zi'], 'swr': r['swr']} for r in rows]
    warn = ['Wires are thin conductors in openEMS (their radius is ignored).'] if any(e['type'] in cad.PATH_TYPES and not e['feed'] for e in doc['elements']) else []
    return out, text, summ, btxt, cad_wires(doc), warn


def results_z(simdir, freqs):
    t, u = read_probe(os.path.join(simdir, 'port_ut_1'))
    ti, i = read_probe(os.path.join(simdir, 'port_it_1'))
    F = np.array(freqs) * 1e6
    return dft(t, u, F) / dft(ti, i, F)


def cad_wires(doc):
    """wire list for the 3D view: wires and the outlines of plates / blocks."""
    out = []
    for k, e in enumerate(doc['elements']):
        if e['type'] in cad.PATH_TYPES:
            for Q in cad.local_paths(e):
                P = cad.world(e, Q) / 1000
                out += [{'arm': str(k), 'feed': bool(e['feed']), 'p': P[j].tolist(), 'q': P[j + 1].tolist()} for j in range(len(P) - 1)]
        elif e['type'] == 'pcb':
            for layer, z in (('top', cad.PCB.h_of(e)), ('gnd', 0.0)):
                for o, hs in cad.PCB.regions(e, layer):
                    for L in [o] + hs:
                        ring = [cad.world(e, [q[0], q[1], z]) / 1000 for q in list(L) + [L[0]]]
                        out += [{'arm': str(k), 'feed': False, 'p': ring[j].tolist(), 'q': ring[j + 1].tolist()} for j in range(len(ring) - 1)]
            f = e['params'].get('feed')
            if f:
                out.append({'arm': str(k), 'feed': True, 'p': (cad.world(e, [f[0], f[1], 0]) / 1000).tolist(), 'q': (cad.world(e, [f[0], f[1], cad.PCB.h_of(e)]) / 1000).tolist()})
        else:
            rects = cad.pixel_rects(e) if e['type'] == 'plate' else [(-e['params']['W'] / 2, -e['params']['H'] / 2, e['params']['W'] / 2, e['params']['H'] / 2)]
            for (x0, y0, x1, y1) in rects:
                ring = [cad.world(e, [x, y, 0]) / 1000 for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0))]
                out += [{'arm': str(k), 'feed': False, 'p': ring[j].tolist(), 'q': ring[j + 1].tolist()} for j in range(4)]
    return out


def render(doc, path, pat, Zw, freqs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import core
    from core import swr
    fig = plt.figure(figsize=(14, 5.2))
    ax = fig.add_subplot(131)
    s11 = [20 * math.log10(max(1e-6, abs((z - doc['z0']) / (z + doc['z0'])))) for z in Zw]
    ax.plot(freqs, s11, color='#2a78d6', lw=2)
    ax.axvline(doc['f0'], color='#e8710a', lw=1, ls='--')
    if doc['bw']: ax.axvspan(doc['f0'] - doc['bw'] / 2, doc['f0'] + doc['bw'] / 2, color='#e8710a', alpha=.12)
    ax.set_xlabel('MHz'); ax.set_ylabel('S11 [dB]'); ax.grid(alpha=.3)
    ax.set_title('%s: reflection (openEMS)' % doc['name'])
    for i, (plane, lab) in enumerate((('xz', 'vertical plane xz'), ('yz', 'vertical plane yz'))):
        cut = core.plane_cut(pat, 'total', plane)
        a = [q[0] for q in cut]; g = np.array([q[1] for q in cut])
        top = np.ceil(np.nanmax(g) / 5) * 5
        axp = fig.add_subplot(1, 3, 2 + i, projection='polar')
        axp.set_theta_zero_location('E')
        axp.plot(np.radians(a), np.clip(g, top - 30, None), color='#e8710a', lw=2)
        axp.set_rlim(top - 30, top); axp.set_rticks([top - 30, top - 20, top - 10, top])
        axp.set_title('%s [dBi], max %.1f' % (lab, np.nanmax(g)), fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=105); plt.close(fig)
