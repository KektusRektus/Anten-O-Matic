# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""Radiation patterns shared by every antenna type: full-sphere grids read from NEC, analytic
patterns for the calculators (patch, dish), the drawn three-plane target and its 3D blend.

Angles in degrees.  theta = angle from +z (0 = zenith), phi = azimuth from +x towards +y.
A pattern dict is JSON-ready:
  {'step': 5, 'th': [...], 'ph': [...], 'total': [[...]], 'rhcp': [[...]] | None, 'lhcp': ...,
   'ground': bool, 'approx': bool}
gains are dBi / dBic, rows = theta, columns = phi (phi = 0 .. 360 - step).
"""
import math
import numpy as np

FLOOR = -60.0


def grid_axes(step, ground=False):
    th = np.arange(0, (90 if ground else 180) + 1e-9, step)
    ph = np.arange(0, 360 - 1e-9, step)
    return th, ph


def rp_for(step, ground=False):
    th, ph = grid_axes(step, ground)
    return (len(th), len(ph), 0.0, 0.0, float(step), float(step))


def read_grid(N, n, fi, nth, nph, pol='total'):
    """gain grid [dB] of one NEC radiation-pattern block (needs the per-point polarised gain
    added to necpp for Anten-O-Matic)."""
    fn = N.nec_gain if pol == 'total' else N.nec_gain_point_rhcp if pol == 'rhcp' else N.nec_gain_point_lhcp
    G = np.array([[fn(n, fi, i, j) for j in range(nph)] for i in range(nth)], dtype=float)
    return np.maximum(G, FLOOR)


def from_nec(N, n, fi, step, ground, pols=('total', 'rhcp', 'lhcp')):
    th, ph = grid_axes(step, ground)
    out = {'step': step, 'th': th.tolist(), 'ph': ph.tolist(), 'ground': bool(ground), 'approx': False}
    for p in ('total', 'rhcp', 'lhcp'):
        out[p] = read_grid(N, n, fi, len(th), len(ph), p).round(2).tolist() if p in pols else None
    return out


def unit(th, ph):
    th, ph = np.radians(th), np.radians(ph)
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], -1)


def directivity(P, th):
    """directivity [dBi] of a linear power grid P[theta, phi] on a regular grid."""
    w = np.sin(np.radians(th))[:, None] * np.ones_like(P)
    w[0] = w[-1] = max(w[1, 0] / 4, 1e-6) if len(th) > 2 else w[0]   # pole rings
    mean = (P * w).sum() / w.sum() if th[-1] > 179 else (P * w).sum() / w.sum() / 2  # hemisphere grid: rest is zero
    return 10 * math.log10(max(P.max(), 1e-30) / max(mean, 1e-30))


def _finish(P, th, ph, gain_dbi, approx_note, ground=False):
    P = np.maximum(P, 1e-12)
    D = directivity(P, th)
    G = 10 * np.log10(P / P.max()) + (gain_dbi if gain_dbi is not None else D)
    return {'step': float(th[1] - th[0]), 'th': th.tolist(), 'ph': ph.tolist(), 'total': np.maximum(G, FLOOR).round(2).tolist(),
            'rhcp': None, 'lhcp': None, 'ground': ground, 'approx': True, 'note': approx_note, 'directivity': D}


def patch_pattern(L_mm, W_mm, eeff, f_mhz, gain_dbi=None, step=5):
    """cavity / two-slot model of a rectangular patch on a large ground plane, feed along x.
    Front hemisphere only; the back side gets a fixed -20 dB leak so the plot is not empty."""
    th, ph = grid_axes(step)
    T, Pp = np.meshgrid(np.radians(th), np.radians(ph), indexing='ij')
    k = 2 * math.pi * f_mhz * 1e6 / 299792458.0
    L, W = L_mm / 1000, W_mm / 1000
    X = k * W / 2 * np.sin(T) * np.sin(Pp)
    elem = np.where(np.abs(X) < 1e-9, 1.0, np.sin(X) / np.where(np.abs(X) < 1e-9, 1, X))
    af = np.cos(k * L / 2 * np.sin(T) * np.cos(Pp))
    Et = np.cos(Pp) * elem * af
    Ep = -np.cos(T) * np.sin(Pp) * elem * af
    P = (np.abs(Et) ** 2 + np.abs(Ep) ** 2)
    P = np.where(T <= math.pi / 2 + 1e-9, P, P.max() * 0.01 * np.cos(T - math.pi) ** 2)
    return _finish(P, th, ph, gain_dbi, 'cavity model of the patch (infinite ground plane), approximate')


def dish_pattern(D_mm, f_mhz, gain_dbi=None, step=2):
    """uniformly illuminated circular aperture facing +z, with a taper-like -25 dB sidelobe floor."""
    from math import pi
    th, ph = grid_axes(step)
    T = np.radians(th)[:, None] * np.ones(len(ph))[None, :]
    k = 2 * pi * f_mhz * 1e6 / 299792458.0
    u = k * D_mm / 2000 * np.sin(T)
    # 2 J1(u)/u via the series / asymptotic form (no scipy in the bundled Python)
    j = np.where(u < 1e-6, 1.0, 2 * _j1(u) / np.where(u < 1e-6, 1, u))
    P = np.where(T <= pi / 2, j ** 2 * ((1 + np.cos(T)) / 2) ** 2, 0) + 10 ** (-3.5)
    P = np.where(T <= pi / 2, P, 10 ** (-3.5) * 0.3)
    return _finish(P, th, ph, gain_dbi, 'aperture model of the dish, approximate')


def _j1(x):
    x = np.asarray(x, float)
    out = np.empty_like(x)
    s = x < 8
    xs = x[s]; y = xs * xs
    out[s] = xs * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606)))))) / \
        (144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y)))))
    xl = x[~s]; z = 8.0 / xl; y = z * z; xx = xl - 2.356194491
    p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))))
    p2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))
    out[~s] = np.sqrt(0.636619772 / xl) * (np.cos(xx) * p1 - z * np.sin(xx) * p2)
    return out

# ---------------------------------------------------------------- drawn target pattern

PLANES = ('xy', 'xz', 'yz')   # xy: angle = azimuth from +x to +y; xz: from +x up to +z; yz: from +y up to +z
NCUT = 72                       # samples per plane, 5 deg apart


def default_target():
    return {'xy': [0.0] * NCUT, 'xz': [0.0] * NCUT, 'yz': [0.0] * NCUT, 'floor': -30.0}


def clean_target(t):
    if not isinstance(t, dict): raise ValueError('target pattern missing')
    out = {'floor': float(t.get('floor', -30.0)), 'mode': 'axi' if t.get('mode') == 'axi' else 'planes'}
    if not -60 <= out['floor'] <= -3: raise ValueError('target floor must be between -60 and -3 dB')
    for p in PLANES:
        v = [float(x) for x in t.get(p, [])]
        if len(v) != NCUT: raise ValueError('target plane %s needs %d samples' % (p, NCUT))
        out[p] = [min(0.0, max(out['floor'], x)) for x in v]
    return out


def _cut(vals, a):
    """periodic linear interpolation of a 72-sample cut at angle(s) a [deg]."""
    v = np.asarray(vals, float)
    x = (np.asarray(a, float) % 360) / (360 / NCUT)
    i = np.floor(x).astype(int) % NCUT
    f = x - np.floor(x)
    return v[i] * (1 - f) + v[(i + 1) % NCUT] * f


def target_grid(t, th, ph):
    """3D target [dB relative, <= 0] blended from the three drawn plane cuts. Each direction takes
    the cuts' values at its projection on each plane, weighted by how close it is to that plane."""
    if t.get('mode') == 'axi':      # axially symmetric around z: only the elevation cut (xz plane) counts
        TH = np.meshgrid(th, ph, indexing='ij')[0]
        return np.clip(_cut(t['xz'], 90 - TH), t['floor'], 0.0)
    U = unit(*np.meshgrid(th, ph, indexing='ij'))
    x, y, z = U[..., 0], U[..., 1], U[..., 2]
    vals, ws = [], []
    for p, (a, b, nrm) in zip(PLANES, ((x, y, z), (x, z, y), (y, z, x))):
        ang = np.degrees(np.arctan2(b, a))
        dist = np.arcsin(np.clip(np.abs(nrm), 0, 1))          # angle between the direction and the plane
        inplane = np.hypot(a, b) > 1e-6
        vals.append(_cut(t[p], ang))
        ws.append(np.where(inplane, 1.0 / (dist ** 2 + 1e-4), 0.0))
    V, Wt = np.array(vals), np.array(ws)
    T = (V * Wt).sum(0) / np.maximum(Wt.sum(0), 1e-12)
    return np.clip(T, t['floor'], 0.0)


def target_peak(t):
    """direction (az, el) of the strongest blended target value."""
    th, ph = grid_axes(5)
    T = target_grid(t, th, ph)
    i, j = np.unravel_index(np.argmax(T), T.shape)
    return float(ph[j]), float(90 - th[i])


def target_directivity(t, step=5):
    th, ph = grid_axes(step)
    T = target_grid(t, th, ph)
    return directivity(10 ** (T / 10), th)


def harrington(radius_m, f_mhz):
    """normal directivity limit [dBi] of an antenna that fits in a sphere of this radius."""
    ka = 2 * math.pi * f_mhz * 1e6 / 299792458.0 * radius_m
    return 10 * math.log10(max(1.5, ka * ka + 2 * ka))


def match(G, T, th, floor, upper_only=False):
    """weighted RMS difference [dB] between a simulated grid (normalised to its own peak) and the
    target, both clipped at the target floor; solid-angle weights."""
    rel = np.clip(G - G.max(), floor, 0)
    w = np.sin(np.radians(th))[:, None] * np.ones_like(G) + 0.02
    if upper_only: w = w * (np.asarray(th)[:, None] <= 90 + 1e-6)
    err = rel - T
    return float(np.sqrt((w * err ** 2).sum() / w.sum())), err
