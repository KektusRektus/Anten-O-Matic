// Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
// SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
// Anten-O-Matic drawing editor (3.2): a Fusion-style 2D sketch on grid paper.
// Create: line, 2-point / centre rectangle, circle, 3-point / centre / tangent arc, polygon, point.
// Modify: trim, extend, break, sketch fillet, offset, mirror, circular / rectangular pattern, drag.
// Constrain: coincident, collinear, concentric, midpoint, parallel, perpendicular, tangent, equal,
// symmetric, horizontal / vertical, fix, and in-picture dimensions (length, distance, point-line,
// angle, radius). Construction geometry helps but is not part of the element.
// A damped least-squares solver keeps every rule. Units: mm, y up.

import * as PL from './patchlib.js';

const $ = id => document.getElementById(id);
const D2R = Math.PI / 180;
const r3 = v => Math.round(v * 1000) / 1000;

let S = null;            // the sketch
let cv, cx, pane, inputBox = null, askBox = null, target = null;
let view = {x: 0, y: 0, z: 4};   // centre [mm] and zoom [px/mm]
let tool = 'line', sel = new Set(), hover = null, chain = null, pend = [], mouse = {x: 0, y: 0, sx: 0, sy: 0}, dragPt = null, dragEnt = null, panning = null, dimPick = [], cpick = [], auxMode = false, toolSel = [];
let seq = 1, boxSel = null;
const HIST = [];

export function emptySketch(name) { return {name: name || 'Drawn element', pts: {o: {x: 0, y: 0, fixed: true}}, ents: [], cons: [], feed: null, grid: true, layer: 'wire', fills: [], probe: null, pcb: {h: 1.6, sub: 'FR4', via_d: 0.8}}; }
const pcbOf = () => (S.pcb = Object.assign({h: 1.6, sub: 'FR4', via_d: 0.8}, S.pcb || {}));
const P = id => S.pts[id];
const E = id => S.ents.find(x => x.id === id);
const nid = p => { let k; do { k = p + (seq++); } while (S.pts[k] || S.ents.some(e => e.id === k) || S.cons.some(c => c.id === k)); return k; };
function push() { HIST.push(JSON.stringify(S)); if (HIST.length > 80) HIST.shift(); }
function undo() { if (!HIST.length) return status('Nothing to undo.'); S = JSON.parse(HIST.pop()); sel.clear(); draw(); }
function revert(msg) { if (HIST.length) S = JSON.parse(HIST.pop()); sel.clear(); status(msg); draw(); }

// ------------------------------------------------------------------ geometry helpers
const sub = (a, b) => ({x: a.x - b.x, y: a.y - b.y});
const add = (a, b, k = 1) => ({x: a.x + b.x * k, y: a.y + b.y * k});
const len = v => Math.hypot(v.x, v.y);
const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const dot = (a, b) => a.x * b.x + a.y * b.y;
const crs = (a, b) => a.x * b.y - a.y * b.x;
const unit = v => { const L = len(v) || 1e-12; return {x: v.x / L, y: v.y / L}; };
const ang = v => Math.atan2(v.y, v.x);
const wrap = a => { while (a > Math.PI) a -= 2 * Math.PI; while (a < -Math.PI) a += 2 * Math.PI; return a; };
const isLine = e => !!e && e.type === 'line';
const isCurve = e => !!e && (e.type === 'arc' || e.type === 'circle');
const RAD = e => e.type === 'circle' ? e.r : dist(P(e.a), P(e.c));
const usesPt = (e, k) => e.a === k || e.b === k || e.c === k;
function arcSpan(e) { const c = P(e.c), a = P(e.a), b = P(e.b); const t0 = ang(sub(a, c)); let t1 = ang(sub(b, c)); while (t1 <= t0 + 1e-9) t1 += 2 * Math.PI; return [t0, t1, dist(a, c)]; }
function entPoints(e, step = 10) {   // polyline of an entity [mm]
  if (e.type === 'line') return [P(e.a), P(e.b)].map(p => ({x: p.x, y: p.y}));
  if (e.type === 'arc') {
    const c = P(e.c), [t0, t1, R] = arcSpan(e), n = Math.max(2, Math.ceil((t1 - t0) / (step * D2R)));
    const out = []; for (let i = 0; i <= n; i++) { const t = t0 + (t1 - t0) * i / n; out.push({x: c.x + R * Math.cos(t), y: c.y + R * Math.sin(t)}); }
    out[0] = {...P(e.a)}; out[n] = {...P(e.b)}; return out;
  }
  if (e.type === 'circle') { const c = P(e.c), n = Math.max(36, Math.ceil(360 / step)), out = []; for (let i = 0; i <= n; i++) { const t = 2 * Math.PI * i / n; out.push({x: c.x + e.r * Math.cos(t), y: c.y + e.r * Math.sin(t)}); } return out; }
  return [];
}
function circum(a, b, c) {
  const d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)); if (Math.abs(d) < 1e-9) return null;
  const a2 = a.x * a.x + a.y * a.y, b2 = b.x * b.x + b.y * b.y, c2 = c.x * c.x + c.y * c.y;
  return {x: (a2 * (b.y - c.y) + b2 * (c.y - a.y) + c2 * (a.y - b.y)) / d, y: (a2 * (c.x - b.x) + b2 * (a.x - c.x) + c2 * (b.x - a.x)) / d};
}
const lineDir = e => sub(P(e.b), P(e.a));
function signedAngle(e1, e2) { const u = lineDir(e1), v = lineDir(e2); return Math.atan2(crs(u, v), dot(u, v)); }
function lineDist(p, e) { const A = P(e.a), d = lineDir(e), L = len(d) || 1e-12; return crs(sub(p, A), d) / L; }   // signed, to the infinite line
const orient = (a, b, c) => (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
function anchor(e) {   // a representative point of an entity (label position)
  if (isLine(e)) { const A = P(e.a), B = P(e.b); return {x: (A.x + B.x) / 2, y: (A.y + B.y) / 2}; }
  if (e.type === 'arc') { const [t0, t1, R] = arcSpan(e), c = P(e.c), t = (t0 + t1) / 2; return {x: c.x + R * Math.cos(t), y: c.y + R * Math.sin(t)}; }
  if (e.type === 'circle') { const c = P(e.c); return {x: c.x + e.r * Math.SQRT1_2, y: c.y + e.r * Math.SQRT1_2}; }
  return {...P(e.a)};
}

// ------------------------------------------------------------------ solver: damped least squares over the free point coordinates
function residuals(extra) {
  const r = [];
  for (const e of S.ents) if (e.type === 'arc') r.push(dist(P(e.b), P(e.c)) - dist(P(e.a), P(e.c)));
  for (const c of S.cons) {
    const e = c.e && E(c.e), e1 = c.e1 && E(c.e1), e2 = c.e2 && E(c.e2), a = c.a && P(c.a), b = c.b && P(c.b), k = c.k && P(c.k);
    switch (c.type) {
      case 'len': if (e) r.push(dist(P(e.a), P(e.b)) - c.value); break;
      case 'dist': if (a && b) r.push(c.mode === 'h' ? (b.x - a.x) - c.s * c.value : c.mode === 'v' ? (b.y - a.y) - c.s * c.value : dist(a, b) - c.value); break;
      case 'pldist': if (k && e) r.push(lineDist(k, e) * (c.s || 1) - c.value); break;
      case 'angle': if (e1 && e2) r.push(wrap(signedAngle(e1, e2) - c.s * c.value * D2R) * 20); break;
      case 'radius': if (e) r.push(RAD(e) - c.value); break;
      case 'horiz': case 'vert':
        if (e) r.push(c.type === 'horiz' ? P(e.b).y - P(e.a).y : P(e.b).x - P(e.a).x);
        else if (a && b) r.push(c.type === 'horiz' ? b.y - a.y : b.x - a.x); break;
      case 'on': if (k && e) r.push(isLine(e) ? lineDist(k, e) : dist(k, P(e.c)) - RAD(e)); break;
      case 'mid': if (k && e) { const A = P(e.a), B = P(e.b); r.push(k.x - (A.x + B.x) / 2, k.y - (A.y + B.y) / 2); } break;
      case 'par': case 'perp': if (e1 && e2) { const u = lineDir(e1), v = lineDir(e2), n = (len(u) * len(v)) || 1e-12; r.push((c.type === 'par' ? crs(u, v) : dot(u, v)) / n * 20); } break;
      case 'colin': if (e1 && e2) r.push(lineDist(P(e2.a), e1), lineDist(P(e2.b), e1)); break;
      case 'tan': if (e1 && e2) {
        if (isLine(e1)) r.push(lineDist(P(e2.c), e1) * (c.s || 1) - RAD(e2));
        else { const d = dist(P(e1.c), P(e2.c)); r.push(c.mode === 'int' ? d - Math.abs(RAD(e1) - RAD(e2)) : d - (RAD(e1) + RAD(e2))); }
      } break;
      case 'eq': if (e1 && e2) r.push(isLine(e1) ? dist(P(e1.a), P(e1.b)) - dist(P(e2.a), P(e2.b)) : RAD(e1) - RAD(e2)); break;
      case 'sym': if (a && b && e) { const m = {x: (a.x + b.x) / 2, y: (a.y + b.y) / 2}, u = unit(lineDir(e)); r.push(lineDist(m, e), dot(sub(b, a), u)); } break;
    }
  }
  if (extra) for (const f of extra) r.push(...f());
  return r;
}
function vars() {
  const v = [];
  for (const p of Object.values(S.pts)) if (!p.fixed) v.push([p, 'x'], [p, 'y']);
  for (const e of S.ents) if (e.type === 'circle') v.push([e, 'r']);
  return v;
}
function solveLin(A, b) {   // Gaussian elimination with partial pivoting
  const n = b.length;
  for (let i = 0; i < n; i++) {
    let m = i; for (let k = i + 1; k < n; k++) if (Math.abs(A[k][i]) > Math.abs(A[m][i])) m = k;
    [A[i], A[m]] = [A[m], A[i]]; [b[i], b[m]] = [b[m], b[i]];
    const d = A[i][i]; if (Math.abs(d) < 1e-14) continue;
    for (let k = i + 1; k < n; k++) { const f = A[k][i] / d; if (!f) continue; for (let j = i; j < n; j++) A[k][j] -= f * A[i][j]; b[k] -= f * b[i]; }
  }
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = Math.abs(A[i][i]) > 1e-14 ? s / A[i][i] : 0; }
  return x;
}
export function solve(extra, iters = 60) {
  const V = vars(); const n = V.length;
  const x0 = V.map(([o, k]) => o[k]);
  let r = residuals(extra), err = r.reduce((s, v) => s + v * v, 0), lam = 1e-3;
  if (!n || err < 1e-14) return Math.sqrt(err);
  for (let it = 0; it < iters && err > 1e-14; it++) {
    const J = r.map(() => Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const [o, k] = V[j], h = 1e-6 * Math.max(1, Math.abs(o[k])); o[k] += h;
      const r2 = residuals(extra); o[k] -= h;
      for (let i = 0; i < r.length; i++) J[i][j] = (r2[i] - r[i]) / h;
    }
    // (J^T J + lam I + w I) dx = -J^T r - w (x - x0): small pull towards where things were (minimum change)
    const w = 1e-6, A = Array.from({length: n}, () => Array(n).fill(0)), g = Array(n).fill(0);
    for (let i = 0; i < r.length; i++) for (let a = 0; a < n; a++) { const Ja = J[i][a]; if (!Ja) continue; g[a] -= Ja * r[i]; for (let b = a; b < n; b++) A[a][b] += Ja * J[i][b]; }
    for (let a = 0; a < n; a++) { for (let b = 0; b < a; b++) A[a][b] = A[b][a]; A[a][a] += lam * (1 + A[a][a]) + w; g[a] -= w * (V[a][0][V[a][1]] - x0[a]); }
    const dx = solveLin(A, g), old = V.map(([o, k]) => o[k]);
    V.forEach(([o, k], j) => { o[k] += dx[j]; });
    const r2 = residuals(extra), e2 = r2.reduce((s, v) => s + v * v, 0);
    if (e2 < err) { r = r2; err = e2; lam = Math.max(1e-9, lam / 3); } else { V.forEach(([o, k], j) => { o[k] = old[j]; }); lam *= 4; if (lam > 1e8) break; }
  }
  for (const e of S.ents) if (e.type === 'circle') e.r = Math.max(1e-3, e.r);
  return Math.sqrt(err);
}
// solve after an edit; on failure undo the edit and say why
function commit(okMsg, badMsg) {
  const err = solve();
  if (err > 1e-3) return revert(badMsg || `That conflicts with the other dimensions and rules (off by ${err.toFixed(3)}): nothing changed.`), false;
  if (okMsg !== undefined) status(okMsg);
  draw(); return true;
}

// ------------------------------------------------------------------ bookkeeping
function cleanup() {
  const has = id => S.ents.some(e => e.id === id);
  S.ents = S.ents.filter(e => !(e.type === 'line' && e.a === e.b) && !(e.type === 'arc' && e.a === e.b));
  for (let pass = 0; pass < 2; pass++) {
    S.cons = S.cons.filter(c => [c.e, c.e1, c.e2].every(id => !id || has(id)) && [c.a, c.b, c.k].every(k => !k || S.pts[k]) && !((c.type === 'dist' || c.type === 'sym' || c.type === 'horiz' || c.type === 'vert') && c.a && c.a === c.b) && !(c.e1 && c.e1 === c.e2));
    for (const k of Object.keys(S.pts)) if (k !== 'o' && !S.ents.some(e => usesPt(e, k))) delete S.pts[k];
  }
  if (S.feed && !has(S.feed.e)) S.feed = null;
  if (S.probe && !S.pts[S.probe]) S.probe = null;
  for (const id of [...sel]) if (!has(id) && !S.pts[id] && !S.cons.some(c => c.id === id)) sel.delete(id);
}
function mergePts(keep, drop) {
  if (keep === drop || !S.pts[keep] || !S.pts[drop]) return;
  if (drop === 'o') [keep, drop] = [drop, keep];
  for (const e of S.ents) for (const f of ['a', 'b', 'c']) if (e[f] === drop) e[f] = keep;
  for (const c of S.cons) for (const f of ['a', 'b', 'k']) if (c[f] === drop) c[f] = keep;
  if (S.pts[drop].fixed) S.pts[keep].fixed = true;
  delete S.pts[drop];
  cleanup();
}
function addPoint(p) { const k = nid('p'); S.pts[k] = {x: r3(p.x), y: r3(p.y)}; return k; }
function addEnt(e) { e.id = nid(e.type === 'line' ? 'l' : e.type === 'arc' ? 'a' : e.type === 'circle' ? 'k' : 'q'); if (auxMode && e.type !== 'point') e.aux = true; if (e.type !== 'point' && S.layer && S.layer !== 'wire') e.layer = S.layer; S.ents.push(e); return e; }
function addCon(c) { c.id = nid('c'); S.cons.push(c); return c; }
function pointEnt(k) { if (!S.ents.some(e => usesPt(e, k))) addEnt({type: 'point', a: k, aux: true}); return k; }
function ptNear(p, tol = 1e-4) { for (const [k, q] of Object.entries(S.pts)) if (dist(p, q) < tol) return k; return null; }
function keepFeed(fn) {   // run an edit that may replace the fed entity, keep the feed where it was
  const f = S.feed && feedPoint(); const r = fn();
  if (f && (!S.feed || !E(S.feed.e) || dist(feedPoint() || {x: 1e9, y: 1e9}, f) > 1e-3)) {
    let best = null; for (const e of S.ents) { if (e.aux || e.type === 'point') continue; const q = nearestOn(e, f); if (q && (!best || q.d < best.d)) best = {d: q.d, e, q: q.q}; }
    S.feed = best && best.d < 1e-2 ? {e: best.e.id, t: fracOn(best.e, best.q)} : null;
  }
  return r;
}
function nearestOn(e, w) {
  const pts = entPoints(e, 2); let best = null;
  for (let i = 0; i < pts.length - 1; i++) { const a = pts[i], ab = sub(pts[i + 1], a), L2 = dot(ab, ab); if (!L2) continue; const u = Math.max(0, Math.min(1, dot(sub(w, a), ab) / L2)), q = add(a, ab, u), d = dist(q, w); if (!best || d < best.d) best = {d, q}; }
  return best;
}

// ------------------------------------------------------------------ planar layers (patch antennas): closed outlines, filled areas, vias, probe
const LAYC = {top: '#e8a03a', gnd: '#3fa7d6'};
function polyArea(P) { let a = 0; for (let i = 0; i < P.length; i++) { const p = P[i], q = P[(i + 1) % P.length]; a += p.x * q.y - q.x * p.y; } return a / 2; }
function inPoly(P, w) { let c = false; for (let i = 0, j = P.length - 1; i < P.length; j = i++) { const a = P[i], b = P[j]; if ((a.y > w.y) !== (b.y > w.y) && w.x < (b.x - a.x) * (w.y - a.y) / (b.y - a.y) + a.x) c = !c; } return c; }
function layerLoops(layer) {   // every closed outline drawn on a layer: {pts, nodes, key}
  const used = S.ents.filter(e => !e.aux && e.layer === layer && (e.type === 'line' || e.type === 'arc'));
  const adj = {}; for (const e of used) for (const k of [e.a, e.b]) (adj[k] = adj[k] || []).push(e);
  const seen = new Set(), out = [];
  for (const e0 of used) {
    if (seen.has(e0.id)) continue;
    const k0 = e0.a; let k = k0, e = e0, pts = [{...P(k0)}], nodes = [0], ids = [], closed = false;
    while (e && !seen.has(e.id)) {
      seen.add(e.id); ids.push(e.id);
      let seg = entPoints(e); if (e.b === k) seg = seg.slice().reverse();
      for (let i = 1; i < seg.length; i++) pts.push(seg[i]);
      k = e.a === k ? e.b : e.a;
      if (k === k0) { closed = true; break; }
      nodes.push(pts.length - 1);
      const nb = adj[k] || []; e = nb.length === 2 ? nb.find(q => q !== e) : null;
    }
    if (closed && (adj[k0] || []).length === 2 && pts.length > 3) { pts.pop(); out.push({pts, nodes, key: ids.sort().join(',')}); }
  }
  for (const e of S.ents) if (e.type === 'circle' && !e.aux && e.layer === layer) { const pts = entPoints(e); pts.pop(); out.push({pts, nodes: [0], key: e.id}); }
  for (const L of out) L.area = Math.abs(polyArea(L.pts));
  return out;
}
function regionAt(loops, w) {   // the filled area a click is in: innermost outline around it, minus the outlines directly inside
  const around = loops.filter(L => inPoly(L.pts, w)).sort((a, b) => a.area - b.area);
  if (!around.length) return null;
  const outer = around[0];
  const inner = loops.filter(L => L !== outer && L.area < outer.area && inPoly(outer.pts, L.pts[0]));
  const holes = inner.filter(L => !inner.some(M => M !== L && M.area > L.area && inPoly(M.pts, L.pts[0])));
  return {outer, holes};
}
function layerRegions(layer) {
  const loops = layerLoops(layer), out = [], keys = new Set();
  for (const f of S.fills || []) if (f.layer === layer) { const r = regionAt(loops, f); if (r && !keys.has(r.outer.key)) { keys.add(r.outer.key); out.push(r); } }
  return out;
}
function fillAt(w) {
  const layers = S.layer === 'gnd' ? ['gnd'] : S.layer === 'top' ? ['top'] : ['top', 'gnd'];
  for (const layer of layers) {
    const loops = layerLoops(layer), r = regionAt(loops, w);
    if (!r) continue;
    push(); S.fills = S.fills || [];
    const i = S.fills.findIndex(f => f.layer === layer && (regionAt(loops, f) || {outer: {}}).outer.key === r.outer.key);
    if (i >= 0) { S.fills.splice(i, 1); status(`Area on the ${layer === 'top' ? 'Top' : 'Ground'} layer emptied.`); }
    else { S.fills.push({layer, x: r3(w.x), y: r3(w.y)}); status(`Area filled on the ${layer === 'top' ? 'Top' : 'Ground'} layer${r.holes.length ? ` (${r.holes.length} hole${r.holes.length > 1 ? 's' : ''} left open)` : ''}.`); }
    draw(); return;
  }
  status(S.layer === 'wire' ? 'Fill: click inside a closed outline drawn on the Top or Ground layer (choose the layer first, then draw).' : 'Fill: click inside a closed outline of this layer (one closed loop of lines / arcs, or a circle).');
}
const isPcb = () => (S.fills || []).length > 0 || S.ents.some(e => e.via) || !!S.probe;

// ------------------------------------------------------------------ view and drawing
const toS = p => ({x: cv.width / 2 + (p.x - view.x) * view.z, y: cv.height / 2 - (p.y - view.y) * view.z});
const toW = (sx, sy) => ({x: view.x + (sx - cv.width / 2) / view.z, y: view.y - (sy - cv.height / 2) / view.z});
function gridStep() { const want = 18 / view.z; const e = Math.pow(10, Math.floor(Math.log10(want))); for (const m of [1, 2, 5, 10]) if (m * e >= want) return m * e; return 10 * e; }
const css = k => getComputedStyle(document.documentElement).getPropertyValue(k).trim();
const GLYPH = {par: '∥', perp: '⊥', colin: '≡', tan: 'T', eq: '=', on: '◦', mid: 'M', sym: '⇆', horiz: 'H', vert: 'V'};
function strokePoly(pts) { cx.beginPath(); pts.forEach((p, i) => i ? cx.lineTo(p.x, p.y) : cx.moveTo(p.x, p.y)); cx.stroke(); }
function draw() {
  if (!cv || !S) return;
  DIMBOX = [];
  const W = cv.width, H = cv.height, x = cx;
  x.fillStyle = css('--ws-bg') || '#0f151b'; x.fillRect(0, 0, W, H);
  const g = gridStep(), tl = toW(0, 0), br = toW(W, H);
  for (const [stp, col, lw] of [[g, css('--ws-grid') || '#1d2a35', 1], [g * 5, css('--ws-grid2') || '#2a3d4c', 1]]) {
    x.strokeStyle = col; x.lineWidth = lw; x.beginPath();
    for (let gx = Math.floor(tl.x / stp) * stp; gx <= br.x; gx += stp) { const s = toS({x: gx, y: 0}).x; x.moveTo(s, 0); x.lineTo(s, H); }
    for (let gy = Math.floor(br.y / stp) * stp; gy <= tl.y; gy += stp) { const s = toS({x: 0, y: gy}).y; x.moveTo(0, s); x.lineTo(W, s); }
    x.stroke();
  }
  const o = toS({x: 0, y: 0}); x.strokeStyle = css('--ws-axis') || '#3b5468'; x.lineWidth = 1.5; x.beginPath(); x.moveTo(o.x, 0); x.lineTo(o.x, H); x.moveTo(0, o.y); x.lineTo(W, o.y); x.stroke();
  const ink = css('--ws-ink') || '#e6edf3', selc = css('--ws-sel') || '#3fc6ff', dimc = css('--ws-dim') || '#ffb020', hov = css('--ws-hover') || '#9ad8ff', auxc = '#8a7bd8';
  // filled areas (ground under top)
  for (const layer of ['gnd', 'top']) {
    const regs = layerRegions(layer); if (!regs.length) continue;
    x.fillStyle = LAYC[layer]; x.globalAlpha = layer === 'top' ? .34 : .22; x.beginPath();
    for (const r of regs) for (const L of [r.outer, ...r.holes]) { L.pts.forEach((p, i) => { const q = toS(p); i ? x.lineTo(q.x, q.y) : x.moveTo(q.x, q.y); }); x.closePath(); }
    x.fill('evenodd'); x.globalAlpha = 1;
  }
  // entities
  for (const e of S.ents) {
    if (e.type === 'point') continue;
    const pts = entPoints(e, 4).map(toS);
    x.strokeStyle = sel.has(e.id) ? selc : hover === e.id ? hov : e.aux ? auxc : e.layer ? LAYC[e.layer] : ink; x.lineWidth = sel.has(e.id) ? 3.5 : e.aux ? 1.5 : 2.5;
    if (e.aux) x.setLineDash([8, 5]);
    strokePoly(pts); x.setLineDash([]);
    if (e.type === 'arc' || e.type === 'circle') { const c = toS(P(e.c)); x.strokeStyle = 'rgba(160,180,200,.5)'; x.lineWidth = 1; x.beginPath(); x.moveTo(c.x - 5, c.y); x.lineTo(c.x + 5, c.y); x.moveTo(c.x, c.y - 5); x.lineTo(c.x, c.y + 5); x.stroke(); }
  }
  // trim / break preview
  if ((tool === 'trim' || tool === 'break') && hover && E(hover) && E(hover).type !== 'point') {
    const pl = cutPlan(E(hover), mouse);
    if (pl && pl.span) { x.strokeStyle = tool === 'trim' ? '#ff5a4f' : '#ffd400'; x.lineWidth = 4; strokePoly(pl.span.map(toS)); }
  }
  // constraint glyphs
  drawGlyphs();
  // points
  for (const [k, p] of Object.entries(S.pts)) {
    const users = S.ents.filter(e => usesPt(e, k));
    if (k !== 'o' && !users.length) continue;
    const s = toS(p), lone = users.length && users.every(e => e.type === 'point');
    x.fillStyle = k === 'o' ? '#e6e6e6' : sel.has(k) ? selc : hover === k ? hov : p.fixed ? '#ff6b6b' : lone ? auxc : '#ff8a1f';
    if (lone) { x.strokeStyle = x.fillStyle; x.lineWidth = 2; x.beginPath(); x.moveTo(s.x - 5, s.y - 5); x.lineTo(s.x + 5, s.y + 5); x.moveTo(s.x + 5, s.y - 5); x.lineTo(s.x - 5, s.y + 5); x.stroke(); }
    else x.fillRect(s.x - 3.5, s.y - 3.5, 7, 7);
  }
  // vias and the probe feed
  const vd = Math.max(5, pcbOf().via_d / 2 * view.z);
  for (const e of S.ents) if (e.via && P(e.a)) { const s = toS(P(e.a)); x.strokeStyle = sel.has(e.a) ? selc : '#e6e6e6'; x.lineWidth = 2; x.fillStyle = '#20262c'; x.beginPath(); x.arc(s.x, s.y, vd, 0, 2 * Math.PI); x.fill(); x.stroke(); x.beginPath(); x.arc(s.x, s.y, vd * .4, 0, 2 * Math.PI); x.stroke(); }
  if (S.probe && P(S.probe)) { const s = toS(P(S.probe)); x.fillStyle = '#000'; x.strokeStyle = '#ffd400'; x.lineWidth = 2.5; x.beginPath(); x.arc(s.x, s.y, 8, 0, 2 * Math.PI); x.fill(); x.stroke(); x.beginPath(); x.arc(s.x, s.y, 3, 0, 2 * Math.PI); x.fillStyle = '#ffd400'; x.fill(); x.font = '700 12px Consolas, monospace'; x.fillText('PROBE', s.x + 11, s.y - 9); }
  // feed
  if (S.feed) { const q = feedPoint(); if (q) { const s = toS(q); x.fillStyle = '#000'; x.strokeStyle = '#ffd400'; x.lineWidth = 2.5; x.beginPath(); x.arc(s.x, s.y, 8, 0, 2 * Math.PI); x.fill(); x.stroke(); x.fillStyle = '#ffd400'; x.font = '700 12px Consolas, monospace'; x.fillText('FEED', s.x + 11, s.y - 9); } }
  // dimensions
  for (const c of S.cons) drawDim(c, dimc);
  // pick highlights for multi-click tools
  for (const q of [...cpick, ...toolSel.map(id => ({e: id}))]) {
    if (q.k) { const s = toS(P(q.k)); x.strokeStyle = '#ffd400'; x.lineWidth = 2; x.strokeRect(s.x - 6, s.y - 6, 12, 12); }
    else if (q.e && E(q.e) && E(q.e).type !== 'point') { x.strokeStyle = '#ffd400'; x.lineWidth = 3; x.globalAlpha = .7; strokePoly(entPoints(E(q.e), 4).map(toS)); x.globalAlpha = 1; }
  }
  // preview of the tool in progress
  const m = snapPoint(mouse.sx, mouse.sy, null, true);
  x.strokeStyle = selc; x.lineWidth = 1.5; x.setLineDash([6, 4]);
  if (tool === 'line' && chain) { const a = toS(P(chain.last)), b = toS(m.p); x.beginPath(); x.moveTo(a.x, a.y); x.lineTo(b.x, b.y); x.stroke(); lenTip(P(chain.last), m.p); }
  if (tool === 'rect' && pend.length === 1) { const a = pend[0], b = m.p; strokePoly([a, {x: b.x, y: a.y}, b, {x: a.x, y: b.y}, a].map(toS)); lenTip(a, {x: b.x, y: a.y}); }
  if (tool === 'crect' && pend.length === 1) { const c = pend[0], d = sub(m.p, c), q = [{x: c.x - d.x, y: c.y - d.y}, {x: c.x + d.x, y: c.y - d.y}, {x: c.x + d.x, y: c.y + d.y}, {x: c.x - d.x, y: c.y + d.y}]; strokePoly([...q, q[0]].map(toS)); lenTip(q[0], q[1]); }
  if (tool === 'circle' && pend.length === 1) { const c = toS(pend[0]); x.beginPath(); x.arc(c.x, c.y, dist(pend[0], m.p) * view.z, 0, 2 * Math.PI); x.stroke(); lenTip(pend[0], m.p, 'R '); }
  if (tool === 'poly' && pend.length === 1) { strokePoly(polyPts(pend[0], m.p, polyN()).concat([m.p]).map(toS)); lenTip(pend[0], m.p, 'R '); }
  if (tool === 'arc' && pend.length === 1) { const a = toS(pend[0]), b = toS(m.p); x.beginPath(); x.moveTo(a.x, a.y); x.lineTo(b.x, b.y); x.stroke(); }
  if (tool === 'arc' && pend.length === 2) { const c = circum(pend[0], pend[1], m.p); if (c) { const R = dist(c, pend[0]), s = toS(c); const t0 = ang(sub(pend[0], c)), t1 = ang(sub(pend[1], c)); const ccw = orient(pend[0], m.p, pend[1]) > 0; x.beginPath(); x.arc(s.x, s.y, R * view.z, -t0, -t1, ccw); x.stroke(); } }
  if (tool === 'carc' && pend.length === 1) { const c = toS(pend[0]); x.beginPath(); x.arc(c.x, c.y, dist(pend[0], m.p) * view.z, 0, 2 * Math.PI); x.stroke(); lenTip(pend[0], m.p, 'R '); }
  if (tool === 'carc' && pend.length === 2) { const c = pend[0], R = dist(c, pend[1]), s = toS(c), t0 = ang(sub(pend[1], c)), t1 = ang(sub(m.p, c)); x.beginPath(); x.arc(s.x, s.y, R * view.z, -t0, -t1, orient(c, pend[1], m.p) > 0); x.stroke(); }
  if (tool === 'tarc' && pend.length === 1) { const t = tanArc(pend.k, m.p); if (t) { const s = toS(t.C), A = P(pend.k); x.beginPath(); x.arc(s.x, s.y, Math.abs(t.R) * view.z, -ang(sub(A, t.C)), -ang(sub(m.p, t.C)), t.R > 0); x.stroke(); } }
  x.setLineDash([]);
  if (m.snap || m.mid || m.on) { const s = toS(m.p); x.strokeStyle = '#1fbf7a'; x.lineWidth = 2; if (m.mid) { x.beginPath(); x.moveTo(s.x, s.y - 8); x.lineTo(s.x + 7, s.y + 5); x.lineTo(s.x - 7, s.y + 5); x.closePath(); x.stroke(); } else if (m.on) { x.beginPath(); x.arc(s.x, s.y, 6, 0, 2 * Math.PI); x.stroke(); } else x.strokeRect(s.x - 7, s.y - 7, 14, 14); }
  if (m.hv) { const s = toS(m.p); x.fillStyle = '#7fd18b'; x.font = '600 11px Consolas, monospace'; x.fillText(m.hv, s.x + 10, s.y + 16); }
  // selection box: solid = wholly inside (dragged to the right), dashed = touching (dragged to the left)
  if (boxSel) {
    const b = boxSel, cross = b.x1 < b.x0, X = Math.min(b.x0, b.x1), Y = Math.min(b.y0, b.y1), Wd = Math.abs(b.x1 - b.x0), Hd = Math.abs(b.y1 - b.y0);
    x.fillStyle = cross ? 'rgba(31,191,122,.10)' : 'rgba(63,198,255,.10)'; x.fillRect(X, Y, Wd, Hd);
    x.strokeStyle = cross ? '#1fbf7a' : '#3fc6ff'; x.lineWidth = 1.5; if (cross) x.setLineDash([6, 4]); x.strokeRect(X, Y, Wd, Hd); x.setLineDash([]);
  }
  // scale and mode
  x.fillStyle = css('--ws-ink') || '#e6edf3'; x.font = '12px Consolas, monospace'; x.fillText(`grid ${g} mm · ${(mouse.x).toFixed(2)}, ${(mouse.y).toFixed(2)} mm${auxMode ? ' · CONSTRUCTION' : ''}${S.layer && S.layer !== 'wire' ? ' · LAYER ' + (S.layer === 'top' ? 'TOP' : 'GROUND') : ''}`, 12, H - 12);
  placeInput(); placeAsk(); freqPreview();
}
// live resonance estimate of the drawn patch (Top layer areas over the substrate)
let FREQKEY = '';
function freqPreview() {
  const box = $('skfreq'); if (!box) return;
  if (!isPcb()) { if (FREQKEY) { box.innerHTML = ''; FREQKEY = ''; } return; }
  const b = pcbOf(), key = JSON.stringify([S.pts, S.ents, S.fills, b, SK.mat ? (SK.mat(b.sub) || {}).eps : 0]);
  if (key === FREQKEY) return; FREQKEY = key;
  const m = b.sub && SK.mat ? SK.mat(b.sub) : null, er = m && m.eps ? m.eps : 1;
  const regs = layerRegions('top').map(r => ({outer: r.outer.pts, holes: r.holes.map(h => h.pts)}));
  const vias = S.ents.filter(e => e.via && P(e.a)).map(e => ({...P(e.a)}));
  let lines = [];
  try { lines = PL.estimate(regs, vias, +b.h, er, +b.via_d); } catch (e) { lines = ['estimate failed: ' + e.message]; }
  box.innerHTML = lines.length ? `<b>≈ resonance</b> <span class="mut">(${b.sub || 'air'} εr ${er}, h ${b.h} mm)</span><br>` + lines.map(t => t.replace(/&/g, '&amp;').replace(/</g, '&lt;')).join('<br>') : '';
}

// ------------------------------------------------------------------ patch library / calculator
function libOpen() {
  let box = $('skpl');
  if (box) { box.remove(); return; }
  box = document.createElement('div'); box.id = 'skpl';
  const subs = $('sksub') ? [...$('sksub').options].map(o => o.value) : ['FR4'];
  const b = pcbOf(), f0 = SK.f0 ? SK.f0() : 2400;
  box.innerHTML = `<div class="lt">▌ Patch library</div>
    <label>Type <select id="pl_t">${PL.TYPES.map(t => `<option value="${t.id}">${t.label}</option>`).join('')}</select></label>
    <div class="lr"><label>f0 [MHz] <input id="pl_f" value="${f0}"></label><label>Polarisation <select id="pl_p"></select></label></div>
    <div class="lr"><label>Substrate <select id="pl_s">${subs.map(v => `<option value="${v}">${v || 'air'}</option>`).join('')}</select></label><label>h [mm] <input id="pl_h" value="${b.h}"></label><label>via Ø <input id="pl_v" value="${b.via_d}"></label></div>
    <div class="lr"><label>Ground side [mm, empty = auto] <input id="pl_g" value=""></label></div>
    <div id="pl_out"></div>
    <div class="lr"><button id="pl_draw" class="hot">Draw it</button><button id="pl_x">Close</button></div>`;
  pane.appendChild(box);
  $('pl_s').value = b.sub || '';
  const read = () => {
    const t = $('pl_t').value, T = PL.TYPES.find(q => q.id === t), sub = $('pl_s').value, m = sub && SK.mat ? SK.mat(sub) : null;
    return {t, T, sub, o: {f: num($('pl_f').value), h: num($('pl_h').value), via_d: num($('pl_v').value) || 0.8, er: m && m.eps ? m.eps : 1, tand: m ? m.tand || 0 : 0, pol: $('pl_p').value, gnd: num($('pl_g').value) || 0}};
  };
  const pols = () => { const T = PL.TYPES.find(q => q.id === $('pl_t').value), v = $('pl_p').value; $('pl_p').innerHTML = T.pol.map(p => `<option value="${p}">${p === 'lin' ? 'linear' : p.toUpperCase()}</option>`).join(''); if (T.pol.includes(v)) $('pl_p').value = v; };
  const calc = () => {
    const r = read();
    if (!(r.o.f > 0) || !(r.o.h > 0)) { $('pl_out').innerHTML = '<span class="err">Type f0 and h.</span>'; return null; }
    const d = PL.design(r.t, r.o);
    $('pl_out').innerHTML = `<div class="mut">${r.T.beam} · ${r.sub || 'air'} εr ${r.o.er}${r.o.tand ? ', tan δ ' + r.o.tand : ''}</div><table>${d.dims.map(([k, v]) => `<tr><td>${k}</td><td>${typeof v === 'number' ? v + ' mm' : v}</td></tr>`).join('')}</table>` + d.notes.map(n => `<div class="mut">${n}</div>`).join('') + '<div class="mut">Formula + openEMS correction: expect ±2-3 %; check the drawing in openEMS and trim.</div>';
    return {r, d};
  };
  pols(); calc();
  box.addEventListener('input', () => calc()); box.addEventListener('change', e => { if (e.target.id === 'pl_t') pols(); calc(); });
  box.addEventListener('keydown', e => e.stopPropagation());
  $('pl_x').onclick = () => box.remove();
  $('pl_draw').onclick = () => {
    const c = calc(); if (!c) return;
    if (S.ents.length && !confirm('Replace the current drawing with the library patch?')) return;
    push();
    const name = S.name, T = c.r.T;
    S = PL.toSketch(c.d.geo, name, {h: c.r.o.h, sub: c.r.sub, via_d: c.r.o.via_d});
    if ($('sksub')) $('sksub').value = c.r.sub; if ($('skh')) $('skh').value = c.r.o.h; if ($('skviad')) $('skviad').value = c.r.o.via_d; if ($('sklayer')) $('sklayer').value = 'top';
    sel.clear(); box.remove(); fitView();
    status(`${T.label} for ${c.r.o.f} MHz drawn: edit it like any drawing, then Insert into CAD and run openEMS.`);
  };
}
function drawGlyphs() {
  const x = cx, count = {}, shown = new Set();
  x.font = '700 12px Consolas, monospace'; x.textAlign = 'center'; x.textBaseline = 'middle';
  const put = (c, w, key) => {
    const tag = key + '|' + c.type; if (c.type === 'eq' || c.type === 'colin' || c.type === 'on') { if (shown.has(tag)) return; shown.add(tag); }   // one of a kind per curve
    const n = count[key] = (count[key] || 0) + 1, s = toS(w), gx = s.x + 12 + (n - 1) * 17, gy = s.y - 14;
    const on = sel.has(c.id), hv = hover === c.id;
    x.fillStyle = 'rgba(8,12,16,.85)'; x.fillRect(gx - 8, gy - 8, 16, 16);
    x.strokeStyle = on ? (css('--ws-sel') || '#3fc6ff') : hv ? '#9ad8ff' : '#3f7a52'; x.lineWidth = 1; x.strokeRect(gx - 8, gy - 8, 16, 16);
    x.fillStyle = on ? (css('--ws-sel') || '#3fc6ff') : '#7fd18b'; x.fillText(GLYPH[c.type], gx, gy + 1);
    DIMBOX.push({id: c.id, x0: gx - 8, x1: gx + 8, y0: gy - 8, y1: gy + 8});
  };
  for (const c of S.cons) {
    if (!GLYPH[c.type]) continue;
    if (c.type === 'on' || c.type === 'mid') { if (P(c.k)) put(c, P(c.k), c.k); continue; }
    if (c.type === 'sym') { if (P(c.a) && P(c.b)) { put(c, P(c.a), c.a); put(c, P(c.b), c.b); } continue; }
    if ((c.type === 'horiz' || c.type === 'vert') && !c.e) { const a = P(c.a), b = P(c.b); if (a && b) put(c, {x: (a.x + b.x) / 2, y: (a.y + b.y) / 2}, c.id); continue; }
    for (const id of [c.e, c.e1, c.e2]) { const e = id && E(id); if (e) put(c, anchor(e), id); }
  }
  x.textAlign = 'left'; x.textBaseline = 'alphabetic';
}
function lenTip(a, b, pre = '') { const s = toS(b); cx.fillStyle = '#9ad8ff'; cx.font = '12px Consolas, monospace'; cx.fillText(pre + dist(a, b).toFixed(2), s.x + 12, s.y - 10); }
let DIMBOX = [];
function dimGeom(c) {   // where a dimension is drawn: {pts: screen polylines, label: screen point, text}
  const off = (c.off ?? 8 / view.z * 3);
  const lineOf = (A, B, o) => { const d = sub(B, A), L = len(d) || 1, n = {x: -d.y / L, y: d.x / L}; return [{x: A.x + n.x * o, y: A.y + n.y * o}, {x: B.x + n.x * o, y: B.y + n.y * o}]; };
  if (c.type === 'len' || c.type === 'dist') {
    let A, B;
    if (c.type === 'len') { const e = E(c.e); if (!e) return null; A = P(e.a); B = P(e.b); } else { A = P(c.a); B = P(c.b); if (!A || !B) return null; }
    let a2 = A, b2 = B;
    if (c.mode === 'h') { const y = Math.max(A.y, B.y) + Math.abs(off); a2 = {x: A.x, y}; b2 = {x: B.x, y}; }
    else if (c.mode === 'v') { const xx = Math.max(A.x, B.x) + Math.abs(off); a2 = {x: xx, y: A.y}; b2 = {x: xx, y: B.y}; }
    else [a2, b2] = lineOf(A, B, off);
    const cur = c.mode === 'h' ? Math.abs(B.x - A.x) : c.mode === 'v' ? Math.abs(B.y - A.y) : dist(A, B);
    return {ext: [[A, a2], [B, b2]], line: [a2, b2], label: {x: (a2.x + b2.x) / 2, y: (a2.y + b2.y) / 2}, text: (c.mode === 'h' ? '↔ ' : c.mode === 'v' ? '↕ ' : '') + cur.toFixed(2), cur};
  }
  if (c.type === 'pldist') {
    const e = E(c.e), K = P(c.k); if (!e || !K) return null;
    const d = lineDist(K, e), n = unit({x: -lineDir(e).y, y: lineDir(e).x}), F = add(K, n, -d);
    return {line: [K, F], label: {x: (K.x + F.x) / 2, y: (K.y + F.y) / 2}, text: '⟂ ' + Math.abs(d).toFixed(2), cur: Math.abs(d)};
  }
  if (c.type === 'angle') {
    const e1 = E(c.e1), e2 = E(c.e2); if (!e1 || !e2) return null;
    // vertex: shared point or the intersection of the two lines
    const sh = [e1.a, e1.b].find(k => k === e2.a || k === e2.b);
    let V = sh ? P(sh) : null;
    const u = lineDir(e1), v = lineDir(e2);
    if (!V) { const A = P(e1.a), B = P(e2.a), d = crs(u, v); if (Math.abs(d) > 1e-9) { const t = crs(sub(B, A), v) / d; V = add(A, u, t); } else V = P(e1.b); }
    const dirFrom = e => { const k = [e.a, e.b].find(q => q !== sh); const far = sh ? P(k) : (dist(P(e.a), V) > dist(P(e.b), V) ? P(e.a) : P(e.b)); return ang(sub(far, V)); };
    const t1 = dirFrom(e1), t2 = dirFrom(e2), R = c.off ?? 30 / view.z;
    const d = wrap(t2 - t1);
    const pts = []; for (let i = 0; i <= 20; i++) { const t = t1 + d * i / 20; pts.push({x: V.x + R * Math.cos(t), y: V.y + R * Math.sin(t)}); }
    const cur = Math.abs(signedAngle(e1, e2)) / D2R;
    return {arc: pts, label: pts[10], text: '∠ ' + cur.toFixed(1) + '°', cur};
  }
  if (c.type === 'radius') {
    const e = E(c.e); if (!e) return null;
    const C = P(e.c), R = RAD(e);
    let t = Math.PI / 4; if (e.type === 'arc') { const [t0, t1] = arcSpan(e); t = (t0 + t1) / 2; }
    const Q = {x: C.x + R * Math.cos(t), y: C.y + R * Math.sin(t)};
    return {line: [C, Q], label: {x: (C.x + Q.x) / 2, y: (C.y + Q.y) / 2}, text: 'R ' + R.toFixed(2), cur: R};
  }
  return null;
}
function drawDim(c, col) {
  const g = dimGeom(c); if (!g) return;
  const x = cx, on = sel.has(c.id);
  x.strokeStyle = on ? (css('--ws-sel') || '#3fc6ff') : col; x.fillStyle = x.strokeStyle; x.lineWidth = 1.2;
  if (g.ext) for (const [a, b] of g.ext) { const s = toS(a), t = toS(b); x.beginPath(); x.moveTo(s.x, s.y); x.lineTo(t.x, t.y); x.stroke(); }
  if (g.line) { const s = toS(g.line[0]), t = toS(g.line[1]); x.beginPath(); x.moveTo(s.x, s.y); x.lineTo(t.x, t.y); x.stroke(); arrow(s, t); if (c.type !== 'radius' && c.type !== 'pldist') arrow(t, s); }
  if (g.arc) { const ps = g.arc.map(toS); strokePoly(ps); arrow(ps[1], ps[0]); arrow(ps[ps.length - 2], ps[ps.length - 1]); }
  if (inputBox && inputBox.dataset.c === c.id) return;
  const L = toS(g.label); x.font = '600 13px Consolas, monospace';
  const txt = g.text, w = x.measureText(txt).width + 10;
  x.fillStyle = css('--ws-bg') || '#0f151b'; x.fillRect(L.x - w / 2, L.y - 10, w, 20);
  x.strokeRect(L.x - w / 2, L.y - 10, w, 20);
  x.fillStyle = on ? (css('--ws-sel') || '#3fc6ff') : col; x.textAlign = 'center'; x.textBaseline = 'middle'; x.fillText(txt, L.x, L.y + 1); x.textAlign = 'left'; x.textBaseline = 'alphabetic';
  DIMBOX.push({id: c.id, dim: true, x0: L.x - w / 2, x1: L.x + w / 2, y0: L.y - 10, y1: L.y + 10});
}
function arrow(from, to) { const a = Math.atan2(to.y - from.y, to.x - from.x), x = cx; x.beginPath(); x.moveTo(to.x, to.y); x.lineTo(to.x - 9 * Math.cos(a - 0.35), to.y - 9 * Math.sin(a - 0.35)); x.lineTo(to.x - 9 * Math.cos(a + 0.35), to.y - 9 * Math.sin(a + 0.35)); x.closePath(); x.fill(); }

// ------------------------------------------------------------------ picking and snapping
function pickPoint(sx, sy, except) {
  let best = null;
  for (const [k, p] of Object.entries(S.pts)) {
    if (k === except || (k !== 'o' && !S.ents.some(e => usesPt(e, k)))) continue;
    const s = toS(p), d = Math.hypot(s.x - sx, s.y - sy);
    if (d < 9 && (!best || d < best.d)) best = {d, k};
  }
  return best && best.k;
}
function pickEnt(sx, sy, filter) {
  const w = toW(sx, sy); let best = null;
  for (const e of S.ents) {
    if (e.type === 'point' || (filter && !filter(e))) continue;
    const pts = entPoints(e, 3);
    for (let i = 0; i < pts.length - 1; i++) {
      const a = pts[i], ab = sub(pts[i + 1], a), L2 = dot(ab, ab); if (!L2) continue;
      const u = Math.max(0, Math.min(1, dot(sub(w, a), ab) / L2)), q = add(a, ab, u), d = dist(q, w) * view.z;
      if (d < 7 && (!best || d < best.d)) best = {d, id: e.id, q};
    }
  }
  return best;
}
function pickDim(sx, sy) { for (const b of DIMBOX) if (sx >= b.x0 && sx <= b.x1 && sy >= b.y0 && sy <= b.y1) return b.id; return null; }
const DRAWING = ['line', 'rect', 'crect', 'circle', 'arc', 'carc', 'tarc', 'poly', 'point'];
function snapPoint(sx, sy, from, extras) {
  const k = pickPoint(sx, sy, dragPt && dragPt.k);
  if (k) return {p: {...P(k)}, snap: k};
  if (extras && DRAWING.includes(tool)) {   // midpoints and curves, Fusion style: they become constraints
    for (const e of S.ents) if (isLine(e)) { const m = anchor(e), s = toS(m); if (Math.hypot(s.x - sx, s.y - sy) < 9) return {p: m, mid: e.id}; }
    const h = pickEnt(sx, sy, e => !chain || !(usesPt(e, chain.last)));
    if (h) { const e = E(h.id); let q = h.q; if (isCurve(e)) { const c = P(e.c); q = add(c, unit(sub(h.q, c)), RAD(e)); } else { const A = P(e.a), d = unit(lineDir(e)); q = add(A, d, dot(sub(h.q, A), d)); } return {p: q, on: h.id}; }
  }
  let p = toW(sx, sy);
  if (S.grid) { const g = gridStep() / 2; p = {x: Math.round(p.x / g) * g, y: Math.round(p.y / g) * g}; }
  let hv = null;
  const base = from || (chain && P(chain.last)) || null;
  if (base && tool === 'line') {
    const d = sub(p, base);
    if (Math.abs(d.y) * view.z < 8 && Math.abs(d.x) > 0) { p.y = base.y; hv = 'H'; }
    else if (Math.abs(d.x) * view.z < 8 && Math.abs(d.y) > 0) { p.x = base.x; hv = 'V'; }
  }
  return {p, hv};
}
function ptKey(m) {   // a point key for a snapped click; snapping onto a midpoint or a curve adds the rule
  if (m.snap) return m.snap;
  const k = addPoint(m.p);
  if (m.mid) addCon({type: 'mid', k, e: m.mid});
  else if (m.on) addCon({type: 'on', k, e: m.on});
  return k;
}

// ------------------------------------------------------------------ tools
const HINTS = {
  select: 'Select: click lines, points, dimensions or rule glyphs (Shift adds), or drag a box on empty paper: to the right = wholly inside, to the left = touching. Drag a point or a line to move it. Del removes (also a via or the probe).',
  line: 'Line: click the start, then every corner. Esc, right-click or a double-click ends the line. Snaps to points, midpoints (▲) and curves (○), H / V automatically.',
  rect: 'Rectangle: click two opposite corners.', crect: 'Centre rectangle: click the centre, then a corner.',
  arc: 'Arc (3 points): click the start, the end, then a point the arc goes through.',
  carc: 'Centre arc: click the centre, the start (sets the radius), then the end. It runs the short way round.',
  tarc: 'Tangent arc: click the end of a line or arc, then where the arc ends. It leaves tangentially.',
  circle: 'Circle: click the centre, then the radius.', poly: 'Polygon: click the centre, then a corner. Sides: the number next to the button.',
  point: 'Point: click to place a construction point (a reference for dimensions and rules).',
  dim: 'Dimension: click a line (length), two points (distance, then pick ↔ / x / y), a point and a line (distance), two lines (angle) or an arc / circle (radius). Type the value.',
  feed: 'Feed: click on a line or arc where the coax connects.',
  trim: 'Trim: click the part of a curve to cut away (red preview). It cuts back to the nearest crossings.',
  extend: 'Extend: click near the end of a line: it grows to the next curve it meets.',
  break: 'Break: click a curve to split it at the crossings either side (yellow preview).',
  fillet: 'Fillet: click the corner point of two lines, then type the radius.',
  mirror: 'Mirror: click the mirror line (any line, construction lines work best).',
  cpat: 'Circular pattern: click the centre point, then type the number of copies.',
  fill: 'Fill: click inside a closed outline of the current layer to fill it with copper (outlines inside it stay open as holes). Click again to empty it.',
  via: 'Via: click where a shorting via joins the top copper to the ground.',
  probe: 'Probe: click where the coax probe comes up from the ground to the top copper (the feed of a patch). Click the probe again to remove it.',
};
const CONS = {
  coin: ['Coincident', 'click two points (they merge), or a point and a line / arc (the point stays on it).'],
  colin: ['Collinear', 'click two lines.'], conc: ['Concentric', 'click two arcs / circles.'],
  mid: ['Midpoint', 'click a point and a line.'], par: ['Parallel', 'click two lines.'], perp: ['Perpendicular', 'click two lines.'],
  tan: ['Tangent', 'click a line and an arc / circle, or two arcs / circles.'], eq: ['Equal', 'click two lines, or two arcs / circles.'],
  sym: ['Symmetric', 'click two points, then the symmetry line.'], horiz: ['Horizontal', 'click a line, or two points.'],
  vert: ['Vertical', 'click a line, or two points.'], fix: ['Fix', 'click a point to pin it (again unpins).'],
};
export function setTool(t) {
  chain = null; pend = []; dimPick = []; cpick = []; tool = t; closeInput(); closeAsk();
  if (t !== 'mirror' && t !== 'cpat') toolSel = [];
  for (const b of document.querySelectorAll('[data-sk]')) b.classList.toggle('on', b.dataset.sk === t);
  for (const b of document.querySelectorAll('[data-skc]')) b.classList.toggle('on', 'k:' + b.dataset.skc === t);
  status(t.startsWith('k:') ? `${CONS[t.slice(2)][0]}: ${CONS[t.slice(2)][1]} Esc ends.` : HINTS[t] || '');
  draw();
}
function onDown(ev) {
  const r = cv.getBoundingClientRect(), sx = (ev.clientX - r.left) * cv.width / r.width, sy = (ev.clientY - r.top) * cv.height / r.height;
  closeInput(); closeAsk();
  if (ev.button === 1 || ev.button === 2 || (ev.button === 0 && ev.altKey)) {
    if (ev.button === 2 && tool === 'line' && chain) { chain = null; draw(); ev.preventDefault(); return; }
    panning = {sx, sy, x: view.x, y: view.y, moved: false}; cv.setPointerCapture(ev.pointerId); ev.preventDefault(); return;
  }
  const m = snapPoint(sx, sy, null, true), w = toW(sx, sy);
  if (tool === 'select') {
    const dId = pickDim(sx, sy);
    if (dId) { if (!ev.shiftKey) sel.clear(); sel.add(dId); draw(); return; }
    const k = pickPoint(sx, sy);
    if (k && k !== 'o') { dragPt = {k, pushed: false}; if (!ev.shiftKey) sel.clear(); sel.add(k); cv.setPointerCapture(ev.pointerId); draw(); return; }
    const e = pickEnt(sx, sy);
    if (!ev.shiftKey && !(e && sel.has(e.id))) sel.clear();
    if (!e) { boxSel = {x0: sx, y0: sy, x1: sx, y1: sy, base: new Set(sel)}; cv.setPointerCapture(ev.pointerId); draw(); return; }
    if (e) { if (sel.has(e.id) && ev.shiftKey) sel.delete(e.id); else { sel.add(e.id); dragEnt = {id: e.id, w0: w, pushed: false, sx, sy}; cv.setPointerCapture(ev.pointerId); } }
    draw(); return;
  }
  if (tool.startsWith('k:')) {
    const k = pickPoint(sx, sy), e = k ? null : pickEnt(sx, sy);
    if (!k && !e) { status('Click a point or a curve. Esc ends.'); return; }
    cpick.push(k ? {k} : {e: e.id});
    const res = applyCons(tool.slice(2), cpick);
    if (res !== 'need') cpick = [];
    draw(); return;
  }
  if (tool === 'line') {
    push();
    if (!chain) { chain = {last: ptKey(m), first: null}; chain.first = chain.last; }
    else {
      const b = ptKey(m); if (b === chain.last) return;
      addLine(chain.last, b, m.hv); chain.last = b;
      if (m.snap && (m.snap === chain.first || S.ents.filter(e => usesPt(e, b)).length > 1)) chain = null;   // closed or joined: done
    }
    settle(); return;
  }
  if (tool === 'point') { push(); const k = ptKey(m); pointEnt(k); settle(); return; }
  if (tool === 'fill') { fillAt(w); return; }
  if (tool === 'via') { push(); const k = ptKey(m); addEnt({type: 'point', a: k, aux: true, via: true}); settle(); status('Via placed. Vias join the top copper to the ground; drag it, dimension it, Del removes it.'); return; }
  if (tool === 'probe' && S.probe && pickPoint(sx, sy) === S.probe) { push(); S.ents = S.ents.filter(e => !(e.probe && e.a === S.probe)); S.probe = null; cleanup(); draw(); status('Probe feed removed.'); return; }
  if (tool === 'probe') { push(); const k = ptKey(m); addEnt({type: 'point', a: k, aux: true, probe: true}); const old = S.probe; S.probe = k; if (old && old !== k) S.ents = S.ents.filter(e => !(e.probe && e.a === old)); cleanup(); settle(); status('Probe feed set: the coax comes up from the ground here.'); return; }
  if (tool === 'rect' || tool === 'crect') {
    if (!pend.length) { pend = [m.p]; pend.m = m; draw(); return; }
    push();
    let q, ck = null;
    if (tool === 'rect') { const a = pend[0], b = m.p; q = [a, {x: b.x, y: a.y}, b, {x: a.x, y: b.y}]; }
    else { const c = pend[0], d = sub(m.p, c); q = [{x: c.x - d.x, y: c.y - d.y}, {x: c.x + d.x, y: c.y - d.y}, {x: c.x + d.x, y: c.y + d.y}, {x: c.x - d.x, y: c.y + d.y}]; }
    const ks = q.map((p, i) => tool === 'rect' && i === 0 && pend.m.snap ? pend.m.snap : tool === 'rect' && i === 2 && m.snap ? m.snap : addPoint(p));
    const ls = []; for (let i = 0; i < 4; i++) ls.push(addLine(ks[i], ks[(i + 1) % 4], i % 2 === 0 ? 'H' : 'V'));
    if (tool === 'crect') {
      ck = pend.m.snap || pointEnt(addPoint(pend[0]));
      const was = auxMode; auxMode = true; const d1 = addLine(ks[0], ks[2]); addLine(ks[1], ks[3]); auxMode = was;
      addCon({type: 'mid', k: ck, e: d1.id});
    }
    pend = []; settle(); return;
  }
  if (tool === 'circle') {
    if (!pend.length) { pend = [m.p]; pend.m = m; draw(); return; }
    push();
    const c = ptKey(pend.m); addEnt({type: 'circle', c, r: r3(Math.max(1e-3, dist(P(c), m.p)))}); pend = []; settle(); return;
  }
  if (tool === 'poly') {
    if (!pend.length) { pend = [m.p]; pend.m = m; draw(); return; }
    push();
    const n = polyN(), c = ptKey(pend.m), R = Math.max(1e-3, dist(P(c), m.p));
    const was = auxMode; auxMode = true; const circ = addEnt({type: 'circle', c, r: r3(R)}); auxMode = was;
    const ks = polyPts(P(c), m.p, n).map(addPoint), ls = [];
    for (let i = 0; i < n; i++) { ls.push(addLine(ks[i], ks[(i + 1) % n])); addCon({type: 'on', k: ks[i], e: circ.id}); }
    for (let i = 1; i < n; i++) addCon({type: 'eq', e1: ls[0].id, e2: ls[i].id});
    pend = []; settle(); return;
  }
  if (tool === 'arc') {
    pend.push(m.p); (pend.ms = pend.ms || []).push(m);
    if (pend.length < 3) { draw(); return; }
    const C = circum(pend[0], pend[1], pend[2]); if (!C) { status('The three points are in a line: no arc.'); pend = []; draw(); return; }
    push();
    let a = ptKey(pend.ms[0]), b = ptKey(pend.ms[1]);
    if (orient(pend[0], pend[2], pend[1]) < 0) [a, b] = [b, a];     // stored counter-clockwise from a to b
    addEnt({type: 'arc', c: addPoint(C), a, b}); pend = []; settle(); return;
  }
  if (tool === 'carc') {
    pend.push(m.p); (pend.ms = pend.ms || []).push(m);
    if (pend.length < 3) { draw(); return; }
    push();
    const c = ptKey(pend.ms[0]), C = P(c), R = dist(C, pend[1]); if (R < 1e-6) { pend = []; revert('Zero radius: no arc.'); return; }
    let a = ptKey(pend.ms[1]), b = addPoint(add(C, unit(sub(pend[2], C)), R));
    if (orient(C, P(a), P(b)) < 0) [a, b] = [b, a];
    addEnt({type: 'arc', c, a, b}); pend = []; settle(); return;
  }
  if (tool === 'tarc') {
    if (!pend.length) {
      const k = pickPoint(sx, sy); if (!k || !tanStart(k)) return status('Click the free end of a line or an arc.');
      pend = [P(k)]; pend.k = k; draw(); return;
    }
    const t = tanArc(pend.k, m.p); if (!t) return status('Move the end off the tangent line.');
    push();
    const s = pend.k, from = tanStart(s), e = ptKey(m), c = addPoint(t.C), arc = addEnt({type: 'arc', c, a: t.R > 0 ? s : e, b: t.R > 0 ? e : s});
    if (isLine(from)) addCon({type: 'tan', e1: from.id, e2: arc.id, s: Math.sign(lineDist(P(c), from)) || 1});
    else { const d = dist(P(from.c), P(c)); addCon({type: 'tan', e1: from.id, e2: arc.id, mode: Math.abs(d - Math.abs(RAD(from) - RAD(arc))) < Math.abs(d - RAD(from) - RAD(arc)) ? 'int' : 'ext'}); }
    pend = [P(e)]; pend.k = e; settle(); return;    // chain on from the new end
  }
  if (tool === 'feed') {
    const e = pickEnt(sx, sy, q => !q.aux); if (!e) return status('Click on a line or an arc (not a construction line).');
    push(); S.feed = {e: e.id, t: fracOn(E(e.id), e.q)}; status('Feed point set. Insert the element into the CAD when you are done.'); draw(); return;
  }
  if (tool === 'trim' || tool === 'break') {
    const h = pickEnt(sx, sy); if (!h) return status('Click on a curve.');
    push(); const ok = keepFeed(() => cutAt(E(h.id), w, tool === 'trim'));
    if (!ok) { HIST.pop(); return draw(); }
    commit(tool === 'trim' ? 'Trimmed.' : 'Split.'); return;
  }
  if (tool === 'extend') { const h = pickEnt(sx, sy, isLine); if (!h) return status('Click near the end of a line.'); extendAt(E(h.id), w); return; }
  if (tool === 'fillet') {
    const k = pickPoint(sx, sy); if (!k) return status('Click the corner point where two lines meet.');
    const ls = S.ents.filter(e => isLine(e) && (e.a === k || e.b === k));
    if (ls.length !== 2 || S.ents.some(e => !isLine(e) && e.type !== 'point' && usesPt(e, k))) return status('A fillet needs a corner of exactly two lines.');
    const L = Math.min(...ls.map(e => len(lineDir(e))));
    ask(P(k), r3(L / 4), 'radius', v => filletAt(k, v)); return;
  }
  if (tool === 'mirror') {
    const h = pickEnt(sx, sy, isLine); if (!h) return status('Click the mirror line.');
    mirrorAbout(E(h.id)); return;
  }
  if (tool === 'cpat') {
    const at = {...m.p};
    ask(at, 6, 'copies', v => circPattern(at, m.snap, Math.round(v))); return;
  }
  if (tool === 'dim') {
    const dId = pickDim(sx, sy); if (dId) { const c = S.cons.find(c => c.id === dId); if (dimGeom(c)) openInput(c); return; }
    const k = pickPoint(sx, sy), e = k ? null : pickEnt(sx, sy);
    if (!k && !e && dimPick.length === 1 && dimPick[0].e) dimPick.push({e: dimPick[0].e});   // a line, then empty paper: its length
    else if (k) dimPick.push({k}); else if (e) dimPick.push({e: e.id}); else { dimPick = []; status('Click a line, a point, an arc or a circle.'); return; }
    const ents = dimPick.filter(q => q.e).map(q => E(q.e)), pts = dimPick.filter(q => q.k).map(q => q.k);
    let c = null;
    if (dimPick.length === 1 && ents.length === 1 && ents[0].type !== 'line') c = {type: 'radius', e: ents[0].id};
    else if (dimPick.length === 2 && ents.length === 2 && ents[0].id === ents[1].id) c = {type: 'len', e: ents[0].id};
    else if (dimPick.length === 2 && ents.length === 2 && ents.every(isLine)) { const sa = signedAngle(ents[0], ents[1]); c = {type: 'angle', e1: ents[0].id, e2: ents[1].id, s: Math.sign(sa) || 1, value: Math.abs(sa) / D2R}; }
    else if (dimPick.length === 2 && pts.length === 2) c = {type: 'dist', a: pts[0], b: pts[1], mode: 'a', s: 1};
    else if (dimPick.length === 2 && ents.length === 1 && pts.length === 1 && isLine(ents[0])) { const d = lineDist(P(pts[0]), ents[0]); c = {type: 'pldist', k: pts[0], e: ents[0].id, s: Math.sign(d) || 1, value: Math.abs(d)}; }
    else if (dimPick.length === 1 && ents.length === 1) { status('Line picked: click a second line for the angle, a point for the distance, or the same line / empty paper for its length.'); draw(); return; }
    else if (dimPick.length === 1) { status('Now a second point or a line.'); draw(); return; }
    else { dimPick = []; status('That combination cannot be dimensioned.'); return; }
    if (c.type === 'len') c.value = dist(P(ents[0].a), P(ents[0].b));
    if (c.type === 'dist') c.value = dist(P(c.a), P(c.b));
    if (c.type === 'radius') c.value = RAD(ents[0]);
    push(); c.value = r3(c.value); addCon(c); dimPick = [];
    draw(); openInput(c); return;
  }
}
function settle() { const err = solve(); if (err > 1e-3) status(`Some rules cannot all be met (off by ${err.toFixed(3)}).`); draw(); }
function addLine(a, b, hv) {
  const e = addEnt({type: 'line', a, b});
  if (hv === 'H') addCon({type: 'horiz', e: e.id}); else if (hv === 'V') addCon({type: 'vert', e: e.id});
  return e;
}
const polyN = () => Math.max(3, Math.min(64, Math.round(+(($('skpolyn') || {}).value) || 6)));
function polyPts(C, V, n) { const R = dist(C, V), t0 = ang(sub(V, C)), out = []; for (let i = 0; i < n; i++) { const t = t0 + 2 * Math.PI * i / n; out.push({x: C.x + R * Math.cos(t), y: C.y + R * Math.sin(t)}); } return out; }
function tanStart(k) {   // the single line / arc ending at point k, if any
  const us = S.ents.filter(e => (isLine(e) || e.type === 'arc') && (e.a === k || e.b === k));
  return us.length === 1 ? us[0] : null;
}
function tanArc(k, Eq) {   // arc from point k leaving tangentially, ending at Eq -> {C, R (signed: + = counter-clockwise)}
  const e = tanStart(k); if (!e) return null;
  const Sp = P(k); let d;
  if (isLine(e)) d = unit(sub(Sp, P(e.a === k ? e.b : e.a)));
  else { const r = sub(Sp, P(e.c)), t = unit({x: -r.y, y: r.x}); d = e.b === k ? t : {x: -t.x, y: -t.y}; }
  const n = {x: -d.y, y: d.x}, v = sub(Eq, Sp), den = 2 * dot(v, n);
  if (Math.abs(den) < 1e-9 || len(v) < 1e-6) return null;
  const R = dot(v, v) / den; return {C: add(Sp, n, R), R};
}

// ------------------------------------------------------------------ constraints
function applyCons(type, picks) {
  const pts = picks.filter(q => q.k).map(q => q.k), ents = picks.filter(q => q.e).map(q => E(q.e)).filter(Boolean);
  const need = m => { status(`${CONS[type][0]}: ${m}`); return 'need'; };
  const bad = m => { status(`${CONS[type][0]}: ${m}`); return 'bad'; };
  const n = picks.length;
  const go = (fn, msg) => { push(); fn(); return commit(`${CONS[type][0]}${msg ? ': ' + msg : '.'}`, `${CONS[type][0]} conflicts with the other dimensions and rules: nothing changed. Remove a rule first (select its glyph, Del).`) ? 'ok' : 'conflict'; };
  switch (type) {
    case 'fix': if (pts.length) { push(); for (const k of pts) if (k !== 'o') S.pts[k].fixed = !S.pts[k].fixed; draw(); status('Fixed / released.'); return 'ok'; } return n ? bad('click a point.') : need('click a point.');
    case 'horiz': case 'vert':
      if (ents.length === 1 && n === 1 && isLine(ents[0])) return go(() => { S.cons = S.cons.filter(c => !((c.type === 'horiz' || c.type === 'vert') && c.e === ents[0].id)); addCon({type, e: ents[0].id}); });
      if (pts.length === 1 && n === 1) return need('now the second point.');
      if (pts.length === 2) return go(() => addCon({type, a: pts[0], b: pts[1]}));
      return bad('pick a line or two points.');
    case 'coin':
      if (n === 1) return need('now the second point or a curve.');
      if (pts.length === 2) return go(() => mergePts(pts[0], pts[1]), 'points merged');
      if (pts.length === 1 && ents.length === 1 && ents[0].type !== 'point') { if (usesPt(ents[0], pts[0])) return bad('that point is already on it.'); return go(() => addCon({type: 'on', k: pts[0], e: ents[0].id})); }
      return bad('pick two points, or a point and a curve.');
    case 'mid':
      if (n === 1) return need('now the other one (point / line).');
      if (pts.length === 1 && ents.length === 1 && isLine(ents[0])) return go(() => addCon({type: 'mid', k: pts[0], e: ents[0].id}));
      return bad('pick a point and a line.');
    case 'colin': case 'par': case 'perp':
      if (n === 1) return ents.length && isLine(ents[0]) ? need('now the second line.') : bad('pick lines.');
      if (ents.length === 2 && ents.every(isLine)) return go(() => addCon({type, e1: ents[0].id, e2: ents[1].id}));
      return bad('pick two lines.');
    case 'conc':
      if (n === 1) return ents.length && isCurve(ents[0]) ? need('now the second arc / circle.') : bad('pick arcs or circles.');
      if (ents.length === 2 && ents.every(isCurve)) return go(() => mergePts(ents[0].c, ents[1].c), 'centres merged');
      return bad('pick two arcs / circles.');
    case 'eq':
      if (n === 1) return ents.length ? need('now the second one.') : bad('pick lines or arcs / circles.');
      if (ents.length === 2 && (ents.every(isLine) || ents.every(isCurve))) return go(() => addCon({type: 'eq', e1: ents[0].id, e2: ents[1].id}));
      return bad('pick two lines, or two arcs / circles.');
    case 'tan': {
      if (n === 1) return ents.length ? need('now the second one.') : bad('pick a line or an arc / circle.');
      if (ents.length !== 2) return bad('pick curves, not points.');
      let [a, b] = ents; if (isCurve(a) && isLine(b)) [a, b] = [b, a];
      if (isLine(a) && isCurve(b)) return go(() => addCon({type: 'tan', e1: a.id, e2: b.id, s: Math.sign(lineDist(P(b.c), a)) || 1}));
      if (isCurve(a) && isCurve(b)) { const d = dist(P(a.c), P(b.c)); return go(() => addCon({type: 'tan', e1: a.id, e2: b.id, mode: Math.abs(d - Math.abs(RAD(a) - RAD(b))) < Math.abs(d - RAD(a) - RAD(b)) ? 'int' : 'ext'})); }
      return bad('two lines cannot be tangent: use Collinear.');
    }
    case 'sym':
      if (pts.length < 2 && !ents.length) return need(pts.length ? 'now the second point.' : 'click the first point.');
      if (pts.length === 2 && ents.length === 1 && isLine(ents[0]) && n === 3) return go(() => addCon({type: 'sym', a: pts[0], b: pts[1], e: ents[0].id}));
      if (pts.length === 2 && n === 2) return need('now the symmetry line.');
      return bad('pick two points, then a line.');
  }
  return 'bad';
}
function consFromSel(type) {   // a constraint button with the right things already selected applies at once
  const picks = [...sel].map(id => S.pts[id] ? {k: id} : E(id) ? {e: id} : null).filter(Boolean);
  if (picks.length) { const r = applyCons(type, picks); if (r === 'ok') { sel.clear(); draw(); return; } if (r === 'conflict') return; }
  setTool('k:' + type);
}

// ------------------------------------------------------------------ modify: trim, break, extend
function curveOf(e) {
  if (isLine(e)) return {L: true, A: P(e.a), B: P(e.b)};
  const C = P(e.c), R = RAD(e);
  if (e.type === 'circle') return {C, R, t0: 0, sp: 2 * Math.PI, full: true};
  const [t0, t1] = arcSpan(e); return {C, R, t0, sp: t1 - t0};
}
function paramOn(cu, p) {
  if (cu.L) { const d = sub(cu.B, cu.A); return dot(sub(p, cu.A), d) / (dot(d, d) || 1e-12); }
  let u = ang(sub(p, cu.C)) - cu.t0; while (u < 0) u += 2 * Math.PI; while (u >= 2 * Math.PI) u -= 2 * Math.PI; return u;
}
function inRange(cu, u, tol = 1e-7) { return cu.L ? u >= -tol && u <= 1 + tol : cu.full || u <= cu.sp + tol || u >= 2 * Math.PI - tol; }
function intersect(c1, c2) {   // of the infinite lines / full circles -> world points
  if (c1.L && c2.L) { const d1 = sub(c1.B, c1.A), d2 = sub(c2.B, c2.A), den = crs(d1, d2); if (Math.abs(den) < 1e-12) return []; return [add(c1.A, d1, crs(sub(c2.A, c1.A), d2) / den)]; }
  if (c1.L || c2.L) {
    const l = c1.L ? c1 : c2, k = c1.L ? c2 : c1, d = sub(l.B, l.A), f = sub(l.A, k.C), a = dot(d, d), b = 2 * dot(f, d), cc = dot(f, f) - k.R * k.R;
    let D = b * b - 4 * a * cc; if (D < 0) { if (D > -1e-9 * a * k.R * k.R) D = 0; else return []; }
    const s = Math.sqrt(D); return (s ? [(-b - s) / (2 * a), (-b + s) / (2 * a)] : [-b / (2 * a)]).map(t => add(l.A, d, t));
  }
  const d = dist(c1.C, c2.C); if (d < 1e-12 || d > c1.R + c2.R + 1e-9 || d < Math.abs(c1.R - c2.R) - 1e-9) return [];
  const a = (c1.R * c1.R - c2.R * c2.R + d * d) / (2 * d), h = Math.sqrt(Math.max(0, c1.R * c1.R - a * a)), u = unit(sub(c2.C, c1.C)), m = add(c1.C, u, a);
  return h < 1e-12 ? [m] : [{x: m.x - u.y * h, y: m.y + u.x * h}, {x: m.x + u.y * h, y: m.y - u.x * h}];
}
function crossings(e) {   // where other curves cross e: [{u, p, o}]
  const cu = curveOf(e), out = [];
  for (const o of S.ents) {
    if (o === e || o.type === 'point') continue;
    const co = curveOf(o);
    for (const p of intersect(cu, co)) { const u = paramOn(cu, p); if (inRange(cu, u) && inRange(co, paramOn(co, p))) out.push({u, p, o}); }
  }
  for (const o of S.ents) if (o.type === 'point' && !usesPt(e, o.a)) {   // construction points lying on e cut it too
    const p = P(o.a), on = isLine(e) ? Math.abs(lineDist(p, e)) : Math.abs(dist(p, cu.C) - cu.R);
    if (on < 1e-6) { const u = paramOn(cu, p); if (inRange(cu, u)) out.push({u, p, o}); }
  }
  return out;
}
function spanPts(cu, u0, u1) {
  if (cu.L) return [add(cu.A, sub(cu.B, cu.A), u0), add(cu.A, sub(cu.B, cu.A), u1)];
  if (u1 < u0) u1 += 2 * Math.PI;
  const n = Math.max(2, Math.ceil((u1 - u0) / (3 * D2R))), out = [];
  for (let i = 0; i <= n; i++) { const t = cu.t0 + u0 + (u1 - u0) * i / n; out.push({x: cu.C.x + cu.R * Math.cos(t), y: cu.C.y + cu.R * Math.sin(t)}); }
  return out;
}
function cutPlan(e, w) {   // the crossings either side of the click, and the span between them
  const cu = curveOf(e), end = cu.L ? 1 : cu.sp, eps = cu.L ? 1e-6 : 1e-6;
  let uc = paramOn(cu, w);
  if (!cu.full && !cu.L && uc > cu.sp) uc = uc - cu.sp < 2 * Math.PI - uc ? cu.sp : 0;
  if (cu.L) uc = Math.max(0, Math.min(1, uc));
  const cs = crossings(e).filter(c => cu.full || (c.u > eps && c.u < end - eps));
  if (cu.full) {
    if (cs.length < 2) return {cu, all: true, span: spanPts(cu, 0, 2 * Math.PI)};
    const sorted = cs.sort((a, b) => a.u - b.u);
    const lo = sorted.filter(c => c.u < uc).pop() || sorted[sorted.length - 1], hi = sorted.find(c => c.u > uc) || sorted[0];
    return {cu, lo, hi, span: spanPts(cu, lo.u, hi.u)};
  }
  const lo = cs.filter(c => c.u < uc).sort((a, b) => b.u - a.u)[0], hi = cs.filter(c => c.u > uc).sort((a, b) => a.u - b.u)[0];
  return {cu, lo, hi, all: !lo && !hi, span: spanPts(cu, lo ? lo.u : 0, hi ? hi.u : end)};
}
function ptAt(c) {   // a point at a crossing: reuse one that is there, else a new one that stays on the other curve
  const k0 = ptNear(c.p, 1e-5); if (k0) return k0;
  const k = addPoint(c.p); S.pts[k].x = c.p.x; S.pts[k].y = c.p.y;
  if (c.o && c.o.type !== 'point') addCon({type: 'on', k, e: c.o.id});
  return k;
}
function copyRules(from, to) {   // a split piece keeps the direction rules of its parent
  for (const c of [...S.cons]) {
    if ((c.type === 'horiz' || c.type === 'vert') && c.e === from.id) addCon({type: c.type, e: to.id});
    if ((c.type === 'par' || c.type === 'perp') && (c.e1 === from.id || c.e2 === from.id)) addCon({type: c.type, e1: c.e1 === from.id ? to.id : c.e1, e2: c.e2 === from.id ? to.id : c.e2});
  }
  if (isLine(from)) addCon({type: 'colin', e1: from.id, e2: to.id});
}
function lengthRulesToPoints(e) {   // before cutting a line, its length dimension becomes a distance between its old ends
  for (const c of S.cons) if (c.type === 'len' && c.e === e.id) { c.type = 'dist'; c.a = e.a; c.b = e.b; c.mode = 'a'; c.s = 1; delete c.e; }
  S.cons = S.cons.filter(c => !((c.type === 'mid' || c.type === 'eq') && (c.e === e.id || c.e1 === e.id || c.e2 === e.id)));
}
function cutAt(e, w, remove) {
  const pl = cutPlan(e, w);
  if (pl.all) { if (!remove) { status('Nothing crosses it: nothing to split.'); return false; } S.ents = S.ents.filter(q => q !== e); cleanup(); return true; }
  if (pl.cu.full) {
    const kl = ptAt(pl.lo), kh = ptAt(pl.hi);
    const r = e.r; delete e.r; e.type = 'arc';
    e.a = kh; e.b = kl;   // what stays runs counter-clockwise from the crossing after the click to the one before it
    if (!remove) { const n = addEnt({type: 'arc', c: e.c, a: kl, b: kh, aux: e.aux}); if (!e.aux) delete n.aux; }
    cleanup(); return true;
  }
  lengthRulesToPoints(e);
  const oldB = e.b, kl = pl.lo && ptAt(pl.lo), kh = pl.hi && ptAt(pl.hi);
  const piece = (a, b) => { const n = addEnt({type: e.type, a, b, c: e.c, aux: e.aux}); if (!e.c) delete n.c; if (!e.aux) delete n.aux; copyRules(e, n); return n; };
  if (remove) {
    if (kl && kh) { e.b = kl; piece(kh, oldB); }
    else if (kl) e.b = kl; else e.a = kh;
  } else {
    if (kl && kh) { e.b = kl; piece(kl, kh); piece(kh, oldB); }
    else { const k = kl || kh; e.b = k; piece(k, oldB); }
  }
  cleanup(); return true;
}
function extendAt(e, w) {
  const A = P(e.a), B = P(e.b), endB = dist(w, B) < dist(w, A), k = endB ? e.b : e.a;
  if (S.ents.some(q => q !== e && q.type !== 'point' && usesPt(q, k))) return status('That end is joined to something else: disconnect it first.');
  const cu = curveOf(e); let best = null;
  for (const o of S.ents) {
    if (o === e || o.type === 'point') continue;
    const co = curveOf(o);
    for (const p of intersect(cu, co)) { const u = paramOn(cu, p); if (!inRange(co, paramOn(co, p))) continue; const ahead = endB ? u - 1 : -u; if (ahead > 1e-7 && (!best || ahead < best.ahead)) best = {ahead, p, o}; }
  }
  if (!best) return status('Nothing to extend to in that direction.');
  push(); S.pts[k].x = best.p.x; S.pts[k].y = best.p.y;
  for (const c of S.cons) if (c.type === 'len' && c.e === e.id) c.value = r3(len(lineDir(e)));
  addCon({type: 'on', k, e: best.o.id}); commit('Extended.');
}

// ------------------------------------------------------------------ modify: fillet, offset, mirror, patterns, construction
function filletAt(k, R) {
  if (!(R > 0)) return status('Type a positive radius.');
  const ls = S.ents.filter(e => isLine(e) && (e.a === k || e.b === k)); if (ls.length !== 2) return;
  const V = {...P(k)}, far = e => P(e.a === k ? e.b : e.a), u1 = unit(sub(far(ls[0]), V)), u2 = unit(sub(far(ls[1]), V));
  const th = Math.acos(Math.max(-1, Math.min(1, dot(u1, u2))));
  if (th < 1e-3 || th > Math.PI - 1e-3) return status('The lines are in line: nothing to round.');
  const t = R / Math.tan(th / 2);
  if (t >= dist(far(ls[0]), V) - 1e-6 || t >= dist(far(ls[1]), V) - 1e-6) return status(`Radius ${r3(R)} is too large for these lines (max ${r3(Math.min(dist(far(ls[0]), V), dist(far(ls[1]), V)) * Math.tan(th / 2))}).`);
  push();
  keepFeed(() => {
    const bis = unit(add(u1, u2)), C = add(V, bis, R / Math.sin(th / 2));
    const k1 = addPoint(add(V, u1, t)), k2 = addPoint(add(V, u2, t)), kc = addPoint(C);
    for (const c of S.cons) if (c.type === 'len' && (c.e === ls[0].id || c.e === ls[1].id)) { const e = E(c.e); c.type = 'dist'; c.a = e.a === k ? e.b : e.a; c.b = k; c.mode = 'a'; c.s = 1; delete c.e; }
    ls[0][ls[0].a === k ? 'a' : 'b'] = k1; ls[1][ls[1].a === k ? 'a' : 'b'] = k2;
    let a = k1, b = k2; if (orient(C, P(a), P(b)) < 0) [a, b] = [b, a];
    const arc = addEnt({type: 'arc', c: kc, a, b}); if (ls[0].aux && ls[1].aux) arc.aux = true; else delete arc.aux;
    for (const l of ls) addCon({type: 'tan', e1: l.id, e2: arc.id, s: Math.sign(lineDist(C, l)) || 1});
    addCon({type: 'radius', e: arc.id, value: r3(R)});
    // the old corner stays as a construction point on both lines, so dimensions to it keep the shape
    pointEnt(k); for (const l of ls) addCon({type: 'on', k, e: l.id});
  });
  commit(`Fillet R ${r3(R)}.`);
}
function selEnts() { return [...sel].map(E).filter(e => e && e.type !== 'point'); }
function offsetSel(d) {
  const es = selEnts(); if (!es.length) return status('Select the curves to offset first.');
  if (!d) return status('Type a non-zero distance (negative goes inwards).');
  push();
  const all = es.flatMap(e => entPoints(e, 10)), G = {x: all.reduce((s, p) => s + p.x, 0) / all.length, y: all.reduce((s, p) => s + p.y, 0) / all.length};
  const ends = {};   // original point -> [{e, field, pos, n (line normal) }]
  const made = [];
  for (const e of es) {
    if (isLine(e)) {
      const A = P(e.a), B = P(e.b), n0 = unit({x: -(B.y - A.y), y: B.x - A.x}), mid = anchor(e);
      const sgn = len(sub(mid, G)) < 1e-9 ? 1 : (dot(n0, sub(mid, G)) >= 0 ? 1 : -1), n = {x: n0.x * sgn * d, y: n0.y * sgn * d};
      const ne = {type: 'line', src: e, pa: add(A, n), pb: add(B, n), n};
      (ends[e.a] = ends[e.a] || []).push({ne, f: 'a'}); (ends[e.b] = ends[e.b] || []).push({ne, f: 'b'}); made.push(ne);
    } else {
      const C = P(e.c), R = RAD(e), out = dist(C, G) <= dist(anchor(e), G) ? 1 : -1, R2 = R + out * d;
      if (R2 <= 1e-6) { revert(`The offset would shrink ${e.type === 'circle' ? 'a circle' : 'an arc'} to nothing.`); return; }
      const ne = {type: e.type, src: e, R2};
      if (e.type === 'arc') { ne.pa = add(C, unit(sub(P(e.a), C)), R2); ne.pb = add(C, unit(sub(P(e.b), C)), R2); (ends[e.a] = ends[e.a] || []).push({ne, f: 'a'}); (ends[e.b] = ends[e.b] || []).push({ne, f: 'b'}); }
      made.push(ne);
    }
  }
  // joints between offset pieces: lines meet at their intersection, anything else at the average
  for (const list of Object.values(ends)) {
    let pos;
    if (list.length === 2 && list.every(q => q.ne.type === 'line')) {
      const [p, q] = list, c1 = {L: true, A: p.ne.pa, B: p.ne.pb}, c2 = {L: true, A: q.ne.pa, B: q.ne.pb}, X = intersect(c1, c2);
      pos = X.length ? X[0] : p.ne['p' + p.f];
    } else pos = list.reduce((s, q) => add(s, q.ne['p' + q.f], 1 / list.length), {x: 0, y: 0});
    const k = addPoint(pos); S.pts[k].x = pos.x; S.pts[k].y = pos.y;
    for (const q of list) q.ne['k' + q.f] = k;
  }
  const newSel = [];
  for (const ne of made) {
    const e = ne.src;
    if (ne.type === 'line') { const n = addEnt({type: 'line', a: ne.ka, b: ne.kb, aux: e.aux}); if (!e.aux) delete n.aux; addCon({type: 'par', e1: e.id, e2: n.id}); const dd = lineDist(P(e.a), n); addCon({type: 'pldist', k: e.a, e: n.id, s: Math.sign(dd) || 1, value: r3(Math.abs(dd))}); newSel.push(n.id); }
    else if (ne.type === 'arc') { const n = addEnt({type: 'arc', c: e.c, a: ne.ka, b: ne.kb, aux: e.aux}); if (!e.aux) delete n.aux; newSel.push(n.id); }
    else { const n = addEnt({type: 'circle', c: e.c, r: r3(ne.R2), aux: e.aux}); if (!e.aux) delete n.aux; newSel.push(n.id); }
  }
  sel = new Set(newSel); commit(`Offset ${r3(d)} mm.`);
}
function reflect(p, A, u) { const v = sub(p, A), t = dot(v, u), f = add(A, u, t); return {x: 2 * f.x - p.x, y: 2 * f.y - p.y}; }
function mirrorAbout(axis) {
  const es = toolSel.map(E).filter(e => e && e !== axis);
  if (!es.length) return status('Select what to mirror first (tool Select), then Mirror.');
  push();
  const A = P(axis.a), u = unit(lineDir(axis)), map = {};
  const mk = k => { if (map[k]) return map[k]; const p = P(k); if (Math.abs(lineDist(p, axis)) < 1e-6) return (map[k] = k); const q = reflect(p, A, u), n = addPoint(q); S.pts[n].x = q.x; S.pts[n].y = q.y; map[k] = n; return n; };
  const pairs = new Map(), made = [];
  for (const e of es) {
    const n = {type: e.type, aux: e.aux};
    if (e.type === 'line') { n.a = mk(e.a); n.b = mk(e.b); }
    else if (e.type === 'arc') { n.c = mk(e.c); n.a = mk(e.b); n.b = mk(e.a); }   // mirroring flips the turning sense
    else if (e.type === 'circle') { n.c = mk(e.c); n.r = e.r; }
    else if (e.type === 'point') { n.a = mk(e.a); }
    if (!e.aux) delete n.aux;
    const ne = addEnt(n); made.push(ne.id);
    if (e.type === 'circle') addCon({type: 'eq', e1: e.id, e2: ne.id});
  }
  for (const [k, n] of Object.entries(map)) if (k !== n) addCon({type: 'sym', a: k, b: n, e: axis.id});
  toolSel = []; sel = new Set(made); setTool('select'); commit('Mirrored. The copies stay symmetric to the line.');
}
function copyEnts(es, fn) {   // copies with every point moved by fn; points fn leaves in place are shared
  const map = {}, made = [];
  const mk = k => { if (map[k]) return map[k]; const p = P(k), q = fn(p); if (dist(p, q) < 1e-7) return (map[k] = k); const n = addPoint(q); S.pts[n].x = q.x; S.pts[n].y = q.y; return (map[k] = n); };
  for (const e of es) { const n = {type: e.type}; for (const f of ['a', 'b', 'c']) if (e[f]) n[f] = mk(e[f]); if (e.r) n.r = e.r; if (e.aux) n.aux = true; const ne = addEnt(n); if (!e.aux) delete ne.aux; made.push(ne.id); if (e.type === 'circle') addCon({type: 'eq', e1: e.id, e2: ne.id}); }
  return made;
}
function circPattern(C, ck, n) {
  const es = toolSel.map(E).filter(Boolean);
  if (!es.length) return status('Select what to repeat first, then Circular pattern.');
  if (!(n >= 2 && n <= 72)) return status('Type 2 to 72 copies (the original counts).');
  push(); const made = [...toolSel];
  for (let i = 1; i < n; i++) { const t = 2 * Math.PI * i / n, co = Math.cos(t), si = Math.sin(t); made.push(...copyEnts(es, p => ({x: C.x + (p.x - C.x) * co - (p.y - C.y) * si, y: C.y + (p.x - C.x) * si + (p.y - C.y) * co}))); }
  toolSel = []; sel = new Set(made); setTool('select'); commit(`${n} × around the centre.`);
}
function rectPattern(txt) {
  const es = selEnts(); if (!es.length) return status('Select what to repeat first.');
  const v = String(txt).replace(/,/g, ' ').split(/\s+/).filter(Boolean).map(Number);
  let nx, ny, dx, dy;
  if (v.length === 3) { nx = v[0]; ny = 1; dx = v[1]; dy = v[2]; } else if (v.length === 4) [nx, ny, dx, dy] = v; else return status('Type "columns rows dx dy" (mm), for example 3 2 20 15, or "count dx dy" for one row.');
  if (![nx, ny, dx, dy].every(Number.isFinite) || nx < 1 || ny < 1 || nx * ny > 200 || nx * ny < 2) return status('Use 2 to 200 copies in all.');
  push(); const made = es.map(e => e.id);
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) if (i || j) made.push(...copyEnts(es, p => ({x: p.x + i * dx, y: p.y + j * dy})));
  sel = new Set(made); commit(`${nx} × ${ny} pattern.`);
}
function toggleAux() {
  const es = selEnts();
  if (!es.length) { auxMode = !auxMode; $('skaux') && $('skaux').classList.toggle('on', auxMode); status(auxMode ? 'Construction mode: new geometry is construction (dashed). It guides, it is not part of the element.' : 'Normal geometry again.'); draw(); return; }
  push(); const to = !es.every(e => e.aux);
  for (const e of es) { if (to) e.aux = true; else delete e.aux; }
  if (S.feed && E(S.feed.e) && E(S.feed.e).aux) S.feed = null;
  status(to ? 'Made construction.' : 'Made normal geometry.'); draw();
}

// ------------------------------------------------------------------ feed, mouse, keys
function fracOn(e, q) {
  const pts = entPoints(e, 3); let tot = 0, at = 0, best = 1e18;
  const segs = []; for (let i = 0; i < pts.length - 1; i++) segs.push(dist(pts[i], pts[i + 1]));
  const T = segs.reduce((a, b) => a + b, 0);
  for (let i = 0; i < pts.length - 1; i++) {
    const a = pts[i], ab = sub(pts[i + 1], a), L2 = dot(ab, ab);
    const u = L2 ? Math.max(0, Math.min(1, dot(sub(q, a), ab) / L2)) : 0, d = dist(add(a, ab, u), q);
    if (d < best) { best = d; at = tot + u * segs[i]; }
    tot += segs[i];
  }
  return T ? at / T : 0.5;
}
function feedPoint() {
  const e = E(S.feed.e); if (!e) return null;
  const pts = entPoints(e, 3), seg = []; let T = 0; for (let i = 0; i < pts.length - 1; i++) { seg.push(dist(pts[i], pts[i + 1])); T += seg[i]; }
  let s = S.feed.t * T, i = 0; while (i < seg.length - 1 && s > seg[i]) { s -= seg[i]; i++; }
  const u = seg[i] ? s / seg[i] : 0; return add(pts[i], sub(pts[i + 1], pts[i]), u);
}
function onMove(ev) {
  const r = cv.getBoundingClientRect(), sx = (ev.clientX - r.left) * cv.width / r.width, sy = (ev.clientY - r.top) * cv.height / r.height;
  mouse = Object.assign(toW(sx, sy), {sx, sy});
  if (panning) { view.x = panning.x - (sx - panning.sx) / view.z; view.y = panning.y + (sy - panning.sy) / view.z; panning.moved = true; draw(); return; }
  if (boxSel) { boxSel.x1 = sx; boxSel.y1 = sy; sel = new Set([...boxSel.base, ...boxHits(boxSel)]); draw(); return; }
  if (dragPt) {
    if (!dragPt.pushed) { push(); dragPt.pushed = true; }
    const m = snapPoint(sx, sy, null), p = m.p, q = S.pts[dragPt.k];
    const wasFixed = q.fixed; q.fixed = false;
    solve([() => [(q.x - p.x) * 30, (q.y - p.y) * 30]], 25);
    q.fixed = wasFixed; draw(); return;
  }
  if (dragEnt) {
    if (!dragEnt.pushed) { if (Math.hypot(sx - dragEnt.sx, sy - dragEnt.sy) < 4) return; push(); dragEnt.pushed = true; const e = E(dragEnt.id); dragEnt.orig = {}; for (const f of ['a', 'b', 'c']) if (e[f] && !P(e[f]).fixed) dragEnt.orig[e[f]] = {...P(e[f])}; }
    const w = toW(sx, sy), d = sub(w, dragEnt.w0), ks = Object.keys(dragEnt.orig);
    solve([() => ks.flatMap(k => [(P(k).x - dragEnt.orig[k].x - d.x) * 30, (P(k).y - dragEnt.orig[k].y - d.y) * 30])], 25);
    draw(); return;
  }
  const k = pickPoint(sx, sy), e = k ? null : pickEnt(sx, sy), d = pickDim(sx, sy);
  hover = d || k || (e && e.id) || null;
  cv.style.cursor = tool === 'select' ? (hover ? 'pointer' : 'default') : 'crosshair';
  draw();
}
// box selection, CAD style: dragged left -> right selects what lies wholly inside, right -> left whatever it touches
function boxHits(b) {
  if (Math.abs(b.x1 - b.x0) < 4 && Math.abs(b.y1 - b.y0) < 4) return [];
  const A = toW(Math.min(b.x0, b.x1), Math.max(b.y0, b.y1)), B = toW(Math.max(b.x0, b.x1), Math.min(b.y0, b.y1));
  const inside = p => p.x >= A.x && p.x <= B.x && p.y >= A.y && p.y <= B.y, cross = b.x1 < b.x0;
  const segHits = (p, q) => {
    if (inside(p) || inside(q)) return true;
    const E = [[{x: A.x, y: A.y}, {x: B.x, y: A.y}], [{x: B.x, y: A.y}, {x: B.x, y: B.y}], [{x: B.x, y: B.y}, {x: A.x, y: B.y}], [{x: A.x, y: B.y}, {x: A.x, y: A.y}]];
    return E.some(([c, d]) => orient(p, q, c) * orient(p, q, d) < 0 && orient(c, d, p) * orient(c, d, q) < 0);
  };
  const out = [];
  for (const e of S.ents) {
    if (e.type === 'point') { if (P(e.a) && inside(P(e.a))) out.push(e.a); continue; }
    const pts = entPoints(e, 4);
    const hit = cross ? pts.some((p, i) => i < pts.length - 1 && segHits(p, pts[i + 1])) : pts.every(inside);
    if (hit) { out.push(e.id); for (const k of [e.a, e.b]) if (k && P(k) && inside(P(k))) out.push(k); }
  }
  for (const [k, p] of Object.entries(S.pts)) if (k !== 'o' && inside(p) && S.ents.some(e => usesPt(e, k)) && !cross) out.push(k);
  return out;
}
function onUp(ev) {
  if (boxSel) { const b = boxSel; boxSel = null; const hits = boxHits(b); sel = new Set([...b.base, ...hits]); if (hits.length) status(`${hits.length} item${hits.length > 1 ? 's' : ''} selected (${b.x1 < b.x0 ? 'touching' : 'wholly inside'} the box). Del removes, drag moves.`); draw(); return; }
  if (panning) { const p = panning; panning = null; if (!p.moved && ev.button === 2) { chain = null; draw(); } return; }
  if (dragPt || dragEnt) {
    const moved = (dragPt && dragPt.pushed) || (dragEnt && dragEnt.pushed); dragPt = null; dragEnt = null;
    if (moved) { const err = solve(); if (err > 1e-3) status(`The dimensions cannot all be met (off by ${err.toFixed(3)}). Remove one or change a value.`); }
    draw();
  }
}
function onDbl() { if (tool === 'line') { chain = null; draw(); } else if (tool === 'tarc') { pend = []; draw(); } else if (tool === 'select') { const id = [...sel].find(k => S.cons.some(c => c.id === k)); const c = id && S.cons.find(c => c.id === id); if (c && dimGeom(c)) openInput(c); } }
function onWheel(ev) {
  ev.preventDefault();
  const r = cv.getBoundingClientRect(), sx = (ev.clientX - r.left) * cv.width / r.width, sy = (ev.clientY - r.top) * cv.height / r.height;
  const before = toW(sx, sy); view.z = Math.max(0.02, Math.min(400, view.z * Math.exp(-ev.deltaY * 0.0015)));
  const after = toW(sx, sy); view.x += before.x - after.x; view.y += before.y - after.y; draw();
}
function startSelTool(t) {   // mirror / circular pattern: remember the selection, then pick the axis / centre
  const es = selEnts(); if (!es.length) return status('Select what to ' + (t === 'mirror' ? 'mirror' : 'repeat') + ' first (tool Select, Shift adds), then press the button.');
  toolSel = es.map(e => e.id); const keep = toolSel; setTool(t); toolSel = keep; draw();
}
function onKey(ev) {
  if (!pane || !pane.classList.contains('on') || /INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName)) return;
  const k = ev.key.toLowerCase();
  if ((ev.ctrlKey || ev.metaKey) && k === 'z') { undo(); ev.preventDefault(); return; }
  if ((ev.ctrlKey || ev.metaKey) && k === 'a') { setTool('select'); sel = new Set([...S.ents.filter(e => e.type !== 'point').map(e => e.id), ...Object.keys(S.pts).filter(q => q !== 'o' && S.ents.some(e => usesPt(e, q))), ...S.cons.map(c => c.id)]); status(`Everything selected (${sel.size} items). Del removes it, Ctrl+Z brings it back.`); draw(); ev.preventDefault(); return; }
  if (ev.ctrlKey || ev.metaKey || ev.altKey) return;
  if (ev.key === 'Escape') { if (chain || pend.length || dimPick.length || cpick.length) { chain = null; pend = []; dimPick = []; cpick = []; draw(); } else setTool('select'); return; }
  if (ev.key === 'Delete' || ev.key === 'Backspace') { delSel(); ev.preventDefault(); return; }
  const map = {l: 'line', a: 'arc', c: 'circle', r: 'rect', d: 'dim', f: 'feed', s: 'select', t: 'trim', e: 'extend', p: 'point', b: 'fill'};
  if (map[k]) setTool(map[k]);
  else if (k === 'h') consFromSel('horiz'); else if (k === 'v') consFromSel('vert');
  else if (k === 'x') toggleAux();
  else if (k === 'o') offsetAsk();
  else if (k === 'm') startSelTool('mirror');
  else if (k === 'g') { S.grid = !S.grid; $('skgrid').checked = S.grid; status('Grid snap ' + (S.grid ? 'on' : 'off')); }
  else if (k === 'z' || k === 'home') fitView();
}
function offsetAsk() { if (!selEnts().length) return status('Select the curves to offset first, then Offset.'); const es = selEnts(); ask(anchor(es[0]), 5, 'offset mm', offsetSel); }
function delSel() {
  if (!sel.size) return;
  push();
  const pts = new Set([...sel].filter(k => S.pts[k] && k !== 'o'));
  // a via or the probe sitting on a corner of other geometry: delete only the via / probe, keep the corner
  for (const k of [...pts]) {
    const users = S.ents.filter(e => usesPt(e, k) && !sel.has(e.id)), marks = users.filter(e => e.via || e.probe);
    if (marks.length && users.length > marks.length) { S.ents = S.ents.filter(e => !marks.includes(e)); if (S.probe === k) S.probe = null; pts.delete(k); }
  }
  if ([...pts].includes(S.probe)) S.probe = null;
  S.ents = S.ents.filter(e => !sel.has(e.id) && !pts.has(e.a) && !pts.has(e.b) && !pts.has(e.c));
  S.cons = S.cons.filter(c => !sel.has(c.id));
  for (const k of pts) delete S.pts[k];
  cleanup(); sel.clear(); solve(); draw();
}

// in-picture value input for dimensions
function openInput(c) {
  closeInput();
  const g = dimGeom(c); if (!g) return;
  inputBox = document.createElement('div'); inputBox.className = 'skin'; inputBox.dataset.c = c.id;
  const modes = c.type === 'dist' ? `<span class="axb">${[['a', '↔'], ['h', 'x'], ['v', 'y']].map(([m, t]) => `<button data-m="${m}" class="${(c.mode || 'a') === m ? 'on' : ''}">${t}</button>`).join('')}</span>` : '';
  inputBox.innerHTML = `<input value="${r3(g.cur)}" spellcheck="false">${modes}<button data-x="1" title="delete this dimension">×</button>`;
  pane.appendChild(inputBox);
  const inp = inputBox.querySelector('input'); inp.focus(); inp.select(); setTimeout(() => { if (document.activeElement !== inp) { inp.focus(); inp.select(); } }, 0);
  inp.onkeydown = ev => { ev.stopPropagation(); if (ev.key === 'Enter') { const v = num(inp.value); if (v === null || v < 0) return status('Type a positive number (mm or °; λ/4 works too).'); push(); c.value = v; const err = solve(); closeInput();
      if (err > 1e-3) { S = JSON.parse(HIST.pop()); status(`${r3(v)} cannot be met together with the other dimensions and rules, so the drawing stays as it was. Remove a conflicting rule (select it, Del) first.`); }
      else status(''); draw(); } else if (ev.key === 'Escape') { closeInput(); draw(); } };
  for (const b of inputBox.querySelectorAll('[data-m]')) b.onmousedown = ev => { ev.preventDefault(); push(); const A = P(c.a), B = P(c.b); c.mode = b.dataset.m; c.s = c.mode === 'h' ? Math.sign(B.x - A.x) || 1 : c.mode === 'v' ? Math.sign(B.y - A.y) || 1 : 1; c.value = r3(dimGeom(c).cur); openInput(c); draw(); };
  inputBox.querySelector('[data-x]').onmousedown = ev => { ev.preventDefault(); push(); S.cons = S.cons.filter(q => q !== c); closeInput(); draw(); };
  placeInput();
}
function placeInput() {
  if (!inputBox) return;
  const c = S.cons.find(q => q.id === inputBox.dataset.c); if (!c) return closeInput();
  const g = dimGeom(c); if (!g) return;
  const L = toS(g.label), r = cv.getBoundingClientRect(), k = r.width / cv.width;
  inputBox.style.left = (cv.offsetLeft + L.x * k) + 'px'; inputBox.style.top = (cv.offsetTop + L.y * k) + 'px';
}
function closeInput() { if (inputBox) { inputBox.remove(); inputBox = null; } }
// a small in-picture prompt for tool values (fillet radius, offset, copies)
function ask(at, init, label, ok) {
  closeAsk();
  askBox = document.createElement('div'); askBox.className = 'skin'; askBox.at = at;
  askBox.innerHTML = `<input value="${init}" spellcheck="false" title="${label}"><button data-x="1" title="cancel">×</button>`;
  pane.appendChild(askBox); placeAsk();
  const inp = askBox.querySelector('input'); inp.focus(); inp.select(); setTimeout(() => { inp.focus(); inp.select(); }, 0);   // after the click's own focus change
  inp.onkeydown = ev => { ev.stopPropagation(); if (ev.key === 'Enter') { const raw = inp.value, v = label === 'pattern' ? raw : num(raw); closeAsk(); if (v === null) return status('Type a number.'); ok(v); } else if (ev.key === 'Escape') { closeAsk(); draw(); } };
  askBox.querySelector('[data-x]').onmousedown = ev => { ev.preventDefault(); closeAsk(); };
  status(`Type the ${label}, Enter.`);
}
function placeAsk() { if (!askBox) return; const L = toS(askBox.at), r = cv.getBoundingClientRect(), k = r.width / cv.width; askBox.style.left = (cv.offsetLeft + L.x * k + 60) + 'px'; askBox.style.top = (cv.offsetTop + L.y * k - 24) + 'px'; }
function closeAsk() { if (askBox) { askBox.remove(); askBox = null; } }
function num(s) {
  const t = String(s).trim().replace(/,/g, '.').replace(/mm|°/g, '').replace(/λ|lambda/gi, '(' + (SK.lam ? SK.lam() : 1000) + ')');
  if (!t || !/^[\d\s.+\-*/()eE]+$/.test(t)) return null;
  try { const v = Function('"use strict";return (' + t + ')')(); return Number.isFinite(v) ? v : null; } catch (e) { return null; }
}
function status(t) { const h = $('skhint'); if (h) h.textContent = t; }
function fitView() {
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
  for (const e of S.ents) for (const p of e.type === 'point' ? [P(e.a)] : entPoints(e)) { x0 = Math.min(x0, p.x); x1 = Math.max(x1, p.x); y0 = Math.min(y0, p.y); y1 = Math.max(y1, p.y); }
  if (x0 > x1 || (x1 - x0 < 1e-6 && y1 - y0 < 1e-6)) { const L = SK.lam ? SK.lam() : 300; x0 = -L / 4; x1 = L / 4; y0 = -L / 6; y1 = L / 6; }
  view.x = (x0 + x1) / 2; view.y = (y0 + y1) / 2; view.z = Math.min(400, 0.8 * Math.min(cv.width / Math.max(1e-3, x1 - x0), cv.height / Math.max(1e-3, y1 - y0)));
  draw();
}

// ------------------------------------------------------------------ export to a CAD element (construction geometry stays behind)
function pcbElement() {
  const lp = L => L.pts.map(p => [r3(p.x), r3(p.y)]);
  const regs = layer => layerRegions(layer).map(r => ({outer: lp(r.outer), nodes: r.outer.nodes, holes: r.holes.map(lp), hnodes: r.holes.map(h => h.nodes)}));
  const top = regs('top'), gnd = regs('gnd');
  if (!top.length) return {error: 'A patch needs at least one filled area on the Top layer (layer Top, draw a closed outline, Fill).'};
  const b = pcbOf(), vias = S.ents.filter(e => e.via && P(e.a)).map(e => [r3(P(e.a).x), r3(P(e.a).y)]);
  const feed = S.probe && P(S.probe) ? [r3(P(S.probe).x), r3(P(S.probe).y)] : null;
  const pts = {};
  return {kind: 'pcb', name: S.name, sketch: JSON.parse(JSON.stringify(S)), params: {top, gnd, vias, feed, h: +b.h, sub: b.sub, via_d: +b.via_d, mesh: 0}, feed};
}
export function toElement() {
  if (isPcb()) return pcbElement();
  const used = S.ents.filter(e => !e.aux && (e.type === 'line' || e.type === 'arc'));
  // graph over point ids: lines and arcs are edges
  const adj = {}; for (const e of used) for (const k of [e.a, e.b]) (adj[k] = adj[k] || []).push(e);
  const done = new Set(), paths = [], nodes = [], entOf = [];
  const walk = (start, first) => {
    const pts = [], nd = [], eo = []; let k = start, e = first;
    pts.push({...P(k)}); nd.push(0);
    while (e && !done.has(e.id)) {
      done.add(e.id);
      let seg = entPoints(e);
      if (e.b === k) seg = seg.slice().reverse();
      const startIdx = pts.length - 1;
      for (let i = 1; i < seg.length; i++) pts.push(seg[i]);
      eo.push({e: e.id, i0: startIdx, i1: pts.length - 1, rev: e.b === k});
      k = e.a === k ? e.b : e.a; nd.push(pts.length - 1);
      const next = (adj[k] || []).filter(q => !done.has(q.id));
      e = (adj[k] || []).length === 2 ? next[0] : null;
    }
    return {pts, nd, eo};
  };
  for (const k of Object.keys(adj)) if (adj[k].length !== 2) for (const e of adj[k]) if (!done.has(e.id)) { const w = walk(k, e); paths.push(w.pts); nodes.push(w.nd); entOf.push(w.eo); }
  for (const e of used) if (!done.has(e.id)) { const w = walk(e.a, e); paths.push(w.pts); nodes.push(w.nd); entOf.push(w.eo); }
  for (const e of S.ents) if (e.type === 'circle' && !e.aux) { const pts = entPoints(e); paths.push(pts); nodes.push([0, pts.length - 1]); entOf.push([{e: e.id, i0: 0, i1: pts.length - 1, rev: false}]); }
  if (!paths.length) return null;
  const pts = {};
  for (const [k, p] of Object.entries(S.pts)) if (k !== 'o' && adj[k]) pts['p' + k.replace(/^p/, '')] = [r3(p.x), r3(p.y), 0];
  let feed = null;
  if (S.feed) {
    paths.forEach((P_, pi) => {
      for (const o of entOf[pi]) if (o.e === S.feed.e) {
        const seg = []; let T = 0; for (let i = 0; i < P_.length - 1; i++) { seg.push(dist(P_[i], P_[i + 1])); T += seg[i]; }
        let s0 = 0; for (let i = 0; i < o.i0; i++) s0 += seg[i];
        let se = 0; for (let i = o.i0; i < o.i1; i++) se += seg[i];
        const t = o.rev ? 1 - S.feed.t : S.feed.t;
        feed = {path: pi, at: T ? (s0 + t * se) / T : 0.5};
      }
    });
  }
  return {name: S.name, sketch: JSON.parse(JSON.stringify(S)), paths: paths.map(P_ => P_.map(p => [r3(p.x), r3(p.y), 0])), nodes, pts, feed};
}

// ------------------------------------------------------------------ public
export const SK = {lam: null, onPlace: null, onClose: null, mat: null, f0: null};
export function init() {
  pane = $('psketch'); cv = $('skcanvas'); cx = cv.getContext('2d');
  const fitCanvas = () => { const r = cv.getBoundingClientRect(); if (!r.width) return; const k = Math.min(2, devicePixelRatio || 1); cv.width = Math.round(r.width * k); cv.height = Math.round(r.height * k); draw(); };
  new ResizeObserver(fitCanvas).observe(cv);
  cv.addEventListener('pointerdown', onDown); cv.addEventListener('pointermove', onMove); cv.addEventListener('pointerup', onUp);
  cv.addEventListener('dblclick', onDbl); cv.addEventListener('wheel', onWheel, {passive: false}); cv.addEventListener('contextmenu', e => e.preventDefault());
  window.addEventListener('keydown', onKey);
  for (const b of document.querySelectorAll('[data-sk]')) b.onclick = () => setTool(b.dataset.sk);
  for (const b of document.querySelectorAll('[data-skc]')) b.onclick = () => consFromSel(b.dataset.skc);
  const on = (id, f) => { const b = $(id); if (b) b.onclick = f; };
  on('skdel', delSel); on('skundo', undo); on('skfit', fitView); on('sklib', libOpen);
  on('skaux', toggleAux); on('skoffset', offsetAsk); on('skmirror', () => startSelTool('mirror')); on('skcpat', () => startSelTool('cpat'));
  on('skrpat', () => { const es = selEnts(); if (!es.length) return status('Select what to repeat first, then Rect. pattern.'); ask(anchor(es[0]), '3 1 20 0', 'pattern', rectPattern); status('Pattern: columns rows dx dy (mm), Enter. Example: 3 2 20 15.'); });
  $('skgrid').onchange = () => { S.grid = $('skgrid').checked; };
  if ($('sklayer')) $('sklayer').onchange = () => {
    const v = $('sklayer').value, es = selEnts();
    if (es.length) { push(); for (const e of es) { if (v === 'wire') delete e.layer; else e.layer = v; } status(`${es.length} curve${es.length > 1 ? 's' : ''} moved to the ${v === 'top' ? 'Top' : v === 'gnd' ? 'Ground' : 'Wire'} layer. New curves go there too.`); }
    else status(v === 'wire' ? 'Layer Wire: curves become wires of a drawn element.' : `Layer ${v === 'top' ? 'Top' : 'Ground'}: draw closed outlines, then Fill them. ${v === 'top' ? 'The top copper sits on the substrate at height h.' : 'The ground plane is at z = 0.'}`);
    S.layer = v; draw();
  };
  const bnum = (id, k) => { if ($(id)) $(id).onchange = () => { const v = num($(id).value); if (v > 0) { pcbOf()[k] = v; draw(); } else $(id).value = pcbOf()[k]; }; };
  bnum('skh', 'h'); bnum('skviad', 'via_d');
  if ($('sksub')) $('sksub').onchange = () => { pcbOf().sub = $('sksub').value; };
  if ($('sksvg')) $('sksvg').onclick = async () => {
    const r = toElement(); if (!r || r.kind !== 'pcb') return status(r && r.error ? r.error : 'SVG export is for patch layouts: fill an area on the Top layer first.');
    try {
      const q = await fetch('/api/pcb_svg', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({element: {id: 'svg', type: 'pcb', name: r.name, mat: 'Copper', params: r.params}})});
      const j = await q.json(); status(q.ok ? 'SVG saved: ' + j.path : 'SVG export failed: ' + j.error);
    } catch (e) { status('SVG export failed: ' + e); }
  };
  $('skname').onchange = () => { S.name = $('skname').value.trim() || S.name; };
  $('skclear').onclick = () => { if (!S.ents.length || confirm('Clear the drawing?')) { push(); const n = S.name; S = emptySketch(n); sel.clear(); draw(); } };
  $('skplace').onclick = () => {
    const r = toElement(); if (!r) return status('Draw at least one normal (not construction) line first.');
    if (r.error) return status(r.error);
    if (r.kind === 'pcb' && !r.feed) status('No probe feed yet: the patch is placed without one (tool Probe adds it).');
    else if (!r.feed) status('No feed point yet: the element is placed without one (tool Feed adds it).');
    SK.onPlace && SK.onPlace(r, target);
  };
  $('skcancel').onclick = () => SK.onClose && SK.onClose();
}
export function open(elem) {
  target = elem || null;
  S = elem && elem.params && elem.params.sketch ? JSON.parse(JSON.stringify(elem.params.sketch)) : emptySketch(elem ? elem.name : 'Drawn element ' + (Math.floor(Math.random() * 90) + 10));
  if (elem) S.name = elem.name;
  HIST.length = 0; sel.clear(); chain = null; pend = []; dimPick = []; cpick = []; toolSel = []; seq = 1;
  $('skname').value = S.name; $('skgrid').checked = S.grid !== false;
  S.fills = S.fills || []; S.layer = S.layer || 'wire'; const b = pcbOf();
  if ($('sklayer')) $('sklayer').value = S.layer; if ($('skh')) $('skh').value = b.h; if ($('skviad')) $('skviad').value = b.via_d;
  if ($('sksub')) { if ($('sksub').options.length > 1 && ![...$('sksub').options].some(o => o.value === b.sub)) b.sub = ''; $('sksub').value = b.sub; }
  $('skplace').textContent = elem ? 'Update element' : 'Insert into CAD';
  $('sktitle').textContent = elem ? `Editing ${elem.name}` : 'New drawn element';
  if (elem && elem.params.sketch_stale) status('The evolution changed this element: updating it from the drawing replaces the evolved shape.');
  setTimeout(() => { fitView(); setTool(S.ents.length ? 'select' : 'line'); }, 30);
}
export function show() { if (!S) open(null); else setTimeout(() => { const r = cv.getBoundingClientRect(); if (r.width && cv.width !== Math.round(r.width * Math.min(2, devicePixelRatio || 1))) { cv.width = Math.round(r.width * Math.min(2, devicePixelRatio || 1)); cv.height = Math.round(r.height * Math.min(2, devicePixelRatio || 1)); } draw(); }, 0); }
export const _test = {layerRegions, fillAt, get: () => S, solve, toElement, setTool, view: () => view, toS, set: s => { S = s; sel.clear(); draw(); }, sel: ids => { sel = new Set(ids); draw(); },
  applyCons, filletAt, offsetSel, rectPattern, cutAt: (id, w, rm) => { push(); const r = cutAt(E(id), w, rm); solve(); draw(); return r; }, extendAt: (id, w) => extendAt(E(id), w),
  mirror: (ids, axis) => { toolSel = ids; mirrorAbout(E(axis)); }, circ: (ids, C, n) => { toolSel = ids; circPattern(C, null, n); }, toggleAux, consFromSel};
