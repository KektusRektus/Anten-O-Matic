// Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
// SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
// Anten-O-Matic patch library (3.7): design formulas for the common microstrip patches, a sketch
// generator (the result is an ordinary, editable drawing) and the live resonance estimate of the
// drawing editor. Units: mm, MHz. Formulas: transmission-line model (rectangular), cavity model
// with fringing radius (circular, annular, shorted annular), quarter-wave shorted patch;
// every family has a correction factor fitted to openEMS runs (FR4 1.6 mm, 2.4 and 5.8 GHz).

export const C0 = 299792.458;   // mm * MHz

// ------------------------------------------------------------------ Bessel functions (Numerical Recipes approximations)
export function J0(x) {
  const ax = Math.abs(x);
  if (ax < 8) { const y = x * x;
    return (57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * -184.9052456))))) /
           (57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 0.785398164;
  return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) -
         z * Math.sin(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)))));
}
export function J1(x) {
  const ax = Math.abs(x);
  if (ax < 8) { const y = x * x;
    return x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * -30.16036606))))) /
           (144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
  const r = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)))) -
            z * Math.sin(xx) * (0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))));
  return x < 0 ? -r : r;
}
export function Y0(x) {
  if (x < 8) { const y = x * x;
    return (-2957821389.0 + y * (7062834065.0 + y * (-512359803.6 + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))))) /
           (40076544269.0 + y * (745249964.8 + y * (7189466.438 + y * (47447.26470 + y * (226.1030244 + y))))) + 0.636619772 * J0(x) * Math.log(x); }
  const z = 8 / x, y = z * z, xx = x - 0.785398164;
  return Math.sqrt(0.636619772 / x) * (Math.sin(xx) * (1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) +
         z * Math.cos(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)))));
}
export function Y1(x) {
  if (x < 8) { const y = x * x;
    return x * (-0.4900604943e13 + y * (0.1275274390e13 + y * (-0.5153438139e11 + y * (0.7349264551e9 + y * (-0.4237922726e7 + y * 0.8511937935e4))))) /
           (0.2499580570e14 + y * (0.4244419664e12 + y * (0.3733650367e10 + y * (0.2245904002e8 + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))))) +
           0.636619772 * (J1(x) * Math.log(x) - 1 / x); }
  const z = 8 / x, y = z * z, xx = x - 2.356194491;
  return Math.sqrt(0.636619772 / x) * (Math.sin(xx) * (1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)))) +
         z * Math.cos(xx) * (0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))));
}
const J1p = x => J0(x) - J1(x) / x, Y1p = x => Y0(x) - Y1(x) / x;

// first root of g(k) above k0 (scan, then bisection)
function firstRoot(g, k0, k1, n = 2000) {
  let a = k0, ga = g(a);
  for (let i = 1; i <= n; i++) {
    const b = k0 + (k1 - k0) * i / n, gb = g(b);
    if (Number.isFinite(ga) && Number.isFinite(gb) && ga * gb < 0) {
      let lo = a, hi = b, glo = ga;
      for (let j = 0; j < 60; j++) { const m = (lo + hi) / 2, gm = g(m); if (glo * gm <= 0) hi = m; else { lo = m; glo = gm; } }
      return (lo + hi) / 2;
    }
    a = b; ga = gb;
  }
  return NaN;
}
function solveFor(fn, f, lo, hi) {      // dimension x with fn(x) = f, fn falling with x
  for (let j = 0; j < 80; j++) { const m = Math.sqrt(lo * hi); if (fn(m) > f) lo = m; else hi = m; }
  return Math.sqrt(lo * hi);
}

// ------------------------------------------------------------------ resonance formulas (MHz)
// openEMS correction of each model (f_openEMS / f_formula), fitted on FR4 1.6 mm
export const CAL = {rect: 0.928, disk: 0.944, sarp: 0.873, ring: 1.0, qw: 0.824, mono: 1.055, redge: 0.44, trunc: 1.04};
// the fitted corrections drift with the electrical thickness x = h / lambda0 (thicker: lower f)
const drift = (f, h, k0, x0, s) => k0 - s * (Math.min(0.045, Math.max(0.006, h / (C0 / f))) - x0);

export function epsEff(W, h, er) { return (er + 1) / 2 + (er - 1) / 2 / Math.sqrt(1 + 12 * h / W); }
export function dL(W, h, er) { const ee = epsEff(W, h, er); return 0.412 * h * (ee + 0.3) * (W / h + 0.264) / ((ee - 0.258) * (W / h + 0.8)); }
// rectangular patch: resonance along L (width W across)
export function rectF(L, W, h, er) { const f = C0 / (2 * (L + 2 * dL(W, h, er)) * Math.sqrt(epsEff(W, h, er))); return CAL.rect * f * drift(CAL.rect * f, h, 1.026, 0.0128, 1.7); }
// quarter-wave patch shorted along one edge of width W: resonance along L
export function qwF(L, W, h, er) { const f = C0 / (4 * (L + dL(W, h, er)) * Math.sqrt(epsEff(W, h, er))); return CAL.qw * f * drift(CAL.qw * f, h, 1.144, 0.0128, 7.5); }
// circular disk: fringing (effective) radius, TMnm with chi = 1.8412 (TM11) or 3.8317 (TM01 / monopolar)
export function aEff(a, h, er) { return a * Math.sqrt(1 + 2 * h / (Math.PI * a * er) * (Math.log(Math.PI * a / (2 * h)) + 1.7726)); }
export function diskF(a, h, er, chi = 1.8412) {
  const f = chi * C0 / (2 * Math.PI * aEff(a, h, er) * Math.sqrt(er));
  return chi === 1.8412 ? CAL.disk * f * drift(CAL.disk * f, h, 1.024, 0.0128, 1.5) : CAL.mono * f * drift(CAL.mono * f, h, 1.188, 0.0128, 1.15);
}
// annular ring TM11: open inner edge (magnetic walls both sides) or shorted inner edge (SARP)
export function ringF(a, b, h, er, shorted) {
  const be = aEff(b, h, er), ae = shorted ? a : Math.max(0.05 * b, a - (aEff(b, h, er) - b));
  const g = shorted ? k => J1p(k * be) * Y1(k * ae) - J1(k * ae) * Y1p(k * be) : k => J1p(k * be) * Y1p(k * ae) - J1p(k * ae) * Y1p(k * be);
  const k = firstRoot(g, 0.05 / be, 6 / (be - ae + 1e-9) + 4 / be);
  const f = k * C0 / (2 * Math.PI * Math.sqrt(er));
  return shorted ? CAL.sarp * f * drift(CAL.sarp * f, h, 1.077, 0.0128, 4.9) : CAL.ring * f;
}

// ------------------------------------------------------------------ library
export const TYPES = [
  {id: 'rect', label: 'Rectangular patch (probe fed)', pol: ['lin'], beam: 'broadside, ~60-80° beamwidth'},
  {id: 'trunc', label: 'Square patch, truncated corners (single-feed CP)', pol: ['rhcp', 'lhcp'], beam: 'broadside, circular polarisation'},
  {id: 'disk', label: 'Circular patch TM11', pol: ['lin'], beam: 'broadside'},
  {id: 'ellipse', label: 'Elliptical patch (single-feed CP)', pol: ['rhcp', 'lhcp'], beam: 'broadside, circular polarisation'},
  {id: 'sarp', label: 'Shorted annular ring TM11 (SARP)', pol: ['lin'], beam: 'broad dome, low surface waves'},
  {id: 'qw', label: 'Quarter-wave shorted patch', pol: ['lin'], beam: 'broad, half the size of a patch'},
  {id: 'mono', label: 'Monopolar disk TM01 (conical beam)', pol: ['lin'], beam: 'null at zenith, peak towards the horizon (like a monopole)'},
];

const circle = (r, n = 72, cx = 0, cy = 0) => Array.from({length: n}, (_, i) => [cx + r * Math.cos(2 * Math.PI * i / n), cy + r * Math.sin(2 * Math.PI * i / n)]);

// substrate Q estimate (Jackson & Alexopoulos radiation Q + dielectric loss) for the CP perturbations
function qTotal(f, h, er, tand) {
  const n = Math.sqrt(er), c1 = 1 - 1 / (n * n) + 2 / (5 * n ** 4), lam = C0 / f;
  const qr = 3 / 16 * er / c1 * (lam / h);
  return 1 / (1 / qr + (tand || 0));
}

// returns {dims: [[label, text]], notes: [..], geo: {top: [{shape}], gnd: {...}, vias, feed}}
export function design(type, o) {
  const f = o.f, h = o.h, er = Math.max(1, o.er), lam = C0 / f, tand = o.tand || 0;
  const r1 = v => Math.round(v * 100) / 100;
  const out = {dims: [], notes: [], geo: {top: [], vias: [], feed: null}};
  const Q = qTotal(f, h, er, tand);
  let extent = 0;
  if (type === 'rect' || type === 'trunc') {
    let W = C0 / (2 * f) * Math.sqrt(2 / (er + 1));
    const kc = 1.4 * (o.cpk || 1);   // corner cut relative to S / sqrt(2Q) (best of the openEMS runs)
    // the CP point lies a few % below the resonance of the cut square: the side is sized for that, not centred
    if (type === 'trunc') W = solveFor(x => CAL.trunc * rectF(x, x, h, er), f, 0.05 * lam, lam);
    const L = type === 'trunc' ? W : solveFor(x => rectF(x, W, h, er), f, 0.05 * lam, lam);
    // probe: R_in = R_edge cos^2(pi x0 / L), R_edge from the radiating-slot conductance
    const Redge = (type === 'trunc' ? 0.5 : 1) * CAL.redge * 90 * er * er / (er - 1 + 1e-9) * (L / W) ** 2;
    const x0 = Redge > 50 ? L / Math.PI * Math.acos(Math.sqrt(50 / Redge)) : L * 0.05;   // from the edge
    const xf = -(L / 2 - x0);
    if (type === 'rect') {
      out.geo.top.push({poly: [[-L / 2, -W / 2], [L / 2, -W / 2], [L / 2, W / 2], [-L / 2, W / 2]]});
      out.dims.push(['L (resonant, along x)', r1(L)], ['W (along y)', r1(W)], ['probe from centre', r1(xf) + ' mm on x']);
    } else {
      const c = L / Math.sqrt(2 * Q) * kc;       // leg of each cut corner: cut area = S^2 / (2 Q)
      const s = o.pol === 'lhcp' ? 1 : -1;                  // which diagonal is cut sets the sense (checked in openEMS)
      const H = L / 2;
      const poly = s > 0
        ? [[-H, -H + c], [-H + c, -H], [H, -H], [H, H - c], [H - c, H], [-H, H]]
        : [[-H, -H], [H - c, -H], [H, -H + c], [H, H], [-H + c, H], [-H, H - c]];
      out.geo.top.push({poly});
      out.dims.push(['side', r1(L)], ['corner cut (each leg)', r1(c)], ['probe from centre', r1(xf) + ' mm on x']);
      out.notes.push(`Q ≈ ${Math.round(Q)}: the cut corners split the two orthogonal modes so they are 90° apart at f0. openEMS (FR4): axial ratio ≈ 3 dB at 5.8 GHz, ≈ 7 dB at 2.4 GHz: trim the cuts in openEMS (bigger cut = more split). The elliptical patch reaches 1-1.5 dB without trimming.`);
    }
    out.geo.feed = [xf, 0]; extent = Math.max(L, W);
  } else if (type === 'disk' || type === 'ellipse' || type === 'mono') {
    const chi = type === 'mono' ? 3.8317 : 1.8412;
    const a = solveFor(x => diskF(x, h, er, chi), f, 0.02 * lam, lam);   // (CP ellipse: the CP point sits ~2 % below its resonance, as wanted)
    if (type === 'mono') {
      out.geo.top.push({circle: a});
      // TM01 has its field maximum in the middle and needs a centre feed; shorting pins on a ring tune the match
      out.geo.feed = [0, 0];
      const rp = (o.pinR || 0.22) * a, np = 4;
      for (let i = 0; i < np; i++) out.geo.vias.push([rp * Math.cos(2 * Math.PI * i / np + Math.PI / 4), rp * Math.sin(2 * Math.PI * i / np + Math.PI / 4)]);
      out.dims.push(['radius', r1(a)], ['centre probe', '(0, 0)'], [`${np} shorting pins on radius`, r1(rp)]);
      out.notes.push('The pins lower the centre impedance of TM01 towards 50 Ω and shift f0 up a little: tune their radius.');
    } else if (type === 'disk') {
      out.geo.top.push({circle: a});
      const rf = 0.32 * a;
      out.geo.feed = [-rf, 0];
      out.dims.push(['radius', r1(a)], ['probe from centre', r1(-rf) + ' mm on x']);
    } else {
      const ar = 1 + 1.4 * 1.0887 / Q * (o.cpk || 1.0);    // Long & Shen axis ratio x 1.4 (best AR in openEMS, FR4)
      const ax = a * Math.sqrt(ar), ay = a / Math.sqrt(ar);
      const poly = Array.from({length: 72}, (_, i) => [ax * Math.cos(2 * Math.PI * i / 72), ay * Math.sin(2 * Math.PI * i / 72)]);
      out.geo.top.push({poly});
      const rf = 0.38 * a, s = o.pol === 'lhcp' ? 1 : -1, t = s * Math.PI / 4;   // sense checked in openEMS
      out.geo.feed = [-rf * Math.cos(t), -rf * Math.sin(t)];
      out.dims.push(['semi-axis x', r1(ax)], ['semi-axis y', r1(ay)], ['probe', `(${r1(out.geo.feed[0])}, ${r1(out.geo.feed[1])}) at 45°`]);
      out.notes.push(`Q ≈ ${Math.round(Q)}: axis ratio ${ar.toFixed(3)}. The probe on the 45° line excites both modes; the other diagonal gives the other sense. openEMS (FR4): axial ratio 1.1 dB at 5.8 GHz, 1.5 dB at 2.4 GHz.`);
    }
    extent = 2 * a * (type === 'ellipse' ? 1.05 : 1);
  } else if (type === 'sarp') {
    const k = 0.45;                                         // inner / outer radius
    const b = solveFor(x => ringF(k * x, x, h, er, true), f, 0.03 * lam, lam), a = k * b;
    const hole = Math.max(0.2, a - 0.5);
    out.geo.top.push({circle: b, hole});
    const n = Math.max(6, Math.round(2 * Math.PI * a / Math.max(2 * (o.via_d || 0.8), 1.6)));
    for (let i = 0; i < n; i++) out.geo.vias.push([a * Math.cos(2 * Math.PI * i / n), a * Math.sin(2 * Math.PI * i / n)]);
    const rf = a + 0.30 * (b - a);
    out.geo.feed = [rf, 0];
    out.dims.push(['outer radius', r1(b)], ['inner short radius (via ring)', r1(a)], ['hole radius', r1(hole)], [`vias`, `${n} × ⌀${o.via_d || 0.8}`], ['probe from centre', r1(rf) + ' mm on x']);
    out.notes.push('Vias short the inner edge: the ring radiates a broad dome with weak surface waves. Linear polarisation; CP needs a second probe at 90° fed in quadrature.');
    extent = 2 * b;
  } else if (type === 'qw') {
    const W = 0.5 * C0 / (2 * f) * Math.sqrt(2 / (er + 1)) * 2;
    const L = solveFor(x => qwF(x, W, h, er), f, 0.02 * lam, lam);
    out.geo.top.push({poly: [[-L / 2, -W / 2], [L / 2, -W / 2], [L / 2, W / 2], [-L / 2, W / 2]]});
    const vd = o.via_d || 0.8, n = Math.max(3, Math.floor(W / (2.5 * vd)));
    for (let i = 0; i < n; i++) out.geo.vias.push([-L / 2 + 0.6 * vd, -W / 2 + W * (i + 0.5) / n]);
    const xf = -L / 2 + 0.25 * L;
    out.geo.feed = [xf, 0];
    out.dims.push(['L (short to open edge)', r1(L)], ['W', r1(W)], [`shorting vias on the left edge`, `${n} × ⌀${vd}`], ['probe from centre', r1(xf) + ' mm on x']);
    out.notes.push('The via row is the short circuit: fewer vias make it inductive and lower f0. Move the probe towards the short for a lower impedance.');
    extent = Math.max(L, W);
  }
  const G = o.gnd || r1(extent + Math.max(12 * h, 0.4 * lam));
  out.geo.gnd = G; out.dims.push(['ground plane', `${G} × ${G}`]);
  out.extent = extent; out.Q = Q;
  out.bw = 100 * (2 - 1) / (Q * Math.sqrt(2));    // SWR < 2 bandwidth of a matched resonator: (s - 1) / (Q sqrt s)
  out.dims.push(['estimated SWR<2 bandwidth', `${out.bw.toFixed(1)} % (${Math.round(out.bw * f / 100)} MHz)`]);
  if (type === 'trunc' || type === 'ellipse') out.dims.push(['polarisation', `${(o.pol || 'rhcp').toUpperCase()} towards +z (up)`]);
  return out;
}

// ------------------------------------------------------------------ geometry -> editable drawing
export function toSketch(geo, name, pcb) {
  const S = {name, pts: {o: {x: 0, y: 0, fixed: true}}, ents: [], cons: [], feed: null, grid: true, layer: 'top', fills: [], probe: null, pcb: Object.assign({h: 1.6, sub: 'FR4', via_d: 0.8}, pcb)};
  let n = 1;
  const r3 = v => Math.round(v * 1000) / 1000;
  const pt = (x, y) => { const k = 'p' + (n++); S.pts[k] = {x: r3(x), y: r3(y)}; return k; };
  const poly = (P, layer) => { const ks = P.map(p => pt(p[0], p[1])); ks.forEach((k, i) => S.ents.push({type: 'line', a: k, b: ks[(i + 1) % ks.length], id: 'l' + (n++), layer})); };
  for (const s of geo.top) {
    if (s.poly) { poly(s.poly, 'top'); const c = s.poly.reduce((a, p) => [a[0] + p[0] / s.poly.length, a[1] + p[1] / s.poly.length], [0, 0]); S.fills.push({layer: 'top', x: r3(c[0] + 0.013), y: r3(c[1] + 0.017)}); }
    else if (s.circle) {
      S.ents.push({type: 'circle', c: 'o', r: r3(s.circle), id: 'k' + (n++), layer: 'top'});
      if (s.hole) S.ents.push({type: 'circle', c: 'o', r: r3(s.hole), id: 'k' + (n++), layer: 'top'});
      const rr = s.hole ? (s.circle + s.hole) / 2 : s.circle / 2;
      S.fills.push({layer: 'top', x: r3(rr * 0.7071 + 0.013), y: r3(rr * 0.7071)});
    }
  }
  const G = geo.gnd / 2; poly([[-G, -G], [G, -G], [G, G], [-G, G]], 'gnd'); S.fills.push({layer: 'gnd', x: r3(-G * 0.93), y: r3(-G * 0.93)});
  for (const v of geo.vias) S.ents.push({type: 'point', a: pt(v[0], v[1]), aux: true, via: true, id: 'q' + (n++)});
  if (geo.feed) { const k = pt(geo.feed[0], geo.feed[1]); S.ents.push({type: 'point', a: k, aux: true, probe: true, id: 'q' + (n++)}); S.probe = k; }
  return S;
}

// ------------------------------------------------------------------ live estimate of a drawn layout
// regions: [{outer: [{x,y}], holes: [[{x,y}]]}] of the top layer; vias: [{x,y}]; returns lines of text
export function estimate(regions, vias, h, er, viaD = 0.8) {
  const out = [];
  if (!(h > 0)) return out;
  const ghz = v => Number.isFinite(v) && v > 0 ? (v >= 1000 ? (v / 1000).toFixed(3) + ' GHz' : v.toFixed(1) + ' MHz') : '?';
  regions.forEach((R, i) => {
    const P = R.outer; if (P.length < 3) return;
    let A = 0, cx = 0, cy = 0;
    for (let j = 0; j < P.length; j++) { const p = P[j], q = P[(j + 1) % P.length], c = p.x * q.y - q.x * p.y; A += c; cx += (p.x + q.x) * c; cy += (p.y + q.y) * c; }
    A /= 2; cx /= 6 * A; cy /= 6 * A; A = Math.abs(A);
    const rs = P.map(p => Math.hypot(p.x - cx, p.y - cy)), rm = rs.reduce((s, v) => s + v, 0) / rs.length;
    const round = P.length >= 12 && Math.max(...rs) / Math.min(...rs) < 1.04 && A / (Math.PI * rm * rm) > 0.95;
    const xs = P.map(p => p.x), ys = P.map(p => p.y), dx = Math.max(...xs) - Math.min(...xs), dy = Math.max(...ys) - Math.min(...ys);
    const tag = regions.length > 1 ? `Area ${i + 1}: ` : '';
    const near = (v, r) => Math.abs(Math.hypot(v.x - cx, v.y - cy) - r);
    if (round && R.holes.length === 1) {
      const H = R.holes[0], hr = H.reduce((s, p) => s + Math.hypot(p.x - cx, p.y - cy), 0) / H.length;
      const inner = vias.filter(v => near(v, hr) < Math.max(1.5 * viaD, 0.2 * (rm - hr)));
      if (inner.length >= 3) { const ri = inner.reduce((s, v) => s + Math.hypot(v.x - cx, v.y - cy), 0) / inner.length; out.push(`${tag}shorted annular ring (SARP): TM11 ≈ ${ghz(ringF(ri, rm, h, er, true))}`); }
      else out.push(`${tag}annular ring: TM11 ≈ ${ghz(ringF(hr, rm, h, er, false))}`);
      return;
    }
    if (round && !R.holes.length) { out.push(`${tag}circular patch r ${rm.toFixed(2)}: TM11 ≈ ${ghz(diskF(rm, h, er))}, TM01 (centre feed) ≈ ${ghz(diskF(rm, h, er, 3.8317))}`); return; }
    const fill = A / (dx * dy);
    // shorting vias along one straight side: quarter-wave patch
    const sides = [['x', Math.min(...xs), dx, dy], ['x', Math.max(...xs), dx, dy], ['y', Math.min(...ys), dy, dx], ['y', Math.max(...ys), dy, dx]];
    for (const [ax, at, L, W] of sides) {
      const on = vias.filter(v => Math.abs((ax === 'x' ? v.x : v.y) - at) < Math.max(1.5 * viaD, 0.06 * L));
      if (on.length >= 2) { out.push(`${tag}quarter-wave patch shorted at ${ax} = ${at.toFixed(1)}: ≈ ${ghz(qwF(L, W, h, er))}`); return; }
    }
    const fx = rectF(dx, dy, h, er), fy = rectF(dy, dx, h, er);
    if (fill > 0.84) {
      const cp = Math.abs(fx - fy) / Math.max(fx, fy) < 0.06 && fill < 0.995 ? `; with the cut corners the CP point is ≈ ${ghz(CAL.trunc * Math.sqrt(fx * fy))}` : '';
      out.push(`${tag}${fill > 0.995 ? 'rectangle' : 'patch'} ${dx.toFixed(2)} × ${dy.toFixed(2)}: TM10 (x) ≈ ${ghz(fx)}, TM01 (y) ≈ ${ghz(fy)}${cp}`);
    } else {
      const s = Math.sqrt(A);
      out.push(`${tag}free shape: rough guess ≈ ${ghz(rectF(s, s, h, er))} (equal-area square), check in openEMS`);
    }
  });
  return out;
}
