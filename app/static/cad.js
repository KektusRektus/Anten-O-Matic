// Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
// SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
// Anten-O-Matic antenna CAD (3.1): elements in 3D, a gumball that can sit on any pivot point and
// move a whole selection, numeric input right at the gumball, dimensions edited in the view
// (distance, x / y / z, angle between elements), a fillet operation on a joint, drawn custom
// elements (sketch.js), undo. The geometry rules mirror app/cad.py exactly.
import * as THREE from 'three';
import { OrbitControls } from '/static/OrbitControls.js';
import { TransformControls } from '/static/TransformControls.js';

const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
const r3 = v => Math.round(v * 1000) / 1000;
const D2R = Math.PI / 180;
const PATH = ['rod', 'fillet', 'torus', 'helix', 'meander', 'custom'];
const LABEL = {rod: 'Rod', fillet: 'Fillet', torus: 'Torus / loop', helix: 'Helix', meander: 'Meander line', custom: 'Drawn element', plate: 'Plate', block: 'Block', pcb: 'Patch / PCB'};
const PARAMS = {
  rod: [['L', 'Length', 'mm']],
  fillet: [['R', 'Radius', 'mm'], ['angle', 'Angle', '°']],
  torus: [['R', 'Radius', 'mm'], ['arc', 'Arc (360 = closed)', '°']],
  helix: [['R', 'Radius', 'mm'], ['pitch', 'Pitch', 'mm'], ['turns', 'Turns', ''], ['hand', 'Hand', '', ['R', 'L']]],
  meander: [['L', 'Length x', 'mm'], ['A', 'Amplitude y', 'mm'], ['n', 'Half-periods', '']],
  custom: [],
  plate: [['W', 'Width x', 'mm'], ['H', 'Height y', 'mm'], ['mesh', 'NEC mesh (0 = auto)', 'mm'], ['t', 'Thickness', 'mm']],
  block: [['W', 'Width x', 'mm'], ['H', 'Height y', 'mm'], ['T', 'Thickness z', 'mm']],
  pcb: [['h', 'Substrate h', 'mm'], ['via_d', 'Via Ø', 'mm'], ['mesh', 'NEC mesh (0 = auto)', 'mm']],
};
const SHAPE_KEYS = {rod: ['L'], fillet: ['R', 'angle'], torus: ['R', 'arc'], helix: ['R', 'pitch', 'turns'], meander: ['L', 'A', 'n'], plate: ['W', 'H'], block: ['W', 'H', 'T'], custom: [], pcb: ['h', 'top width / height']};

export const CAD = {doc: null, mats: [], lam: () => 1000, wire: () => 2, onChange: null, editSketch: null};
let DOC = null, SEL = [], seq = 1, dseq = 1;
let renderer, scene, camera, orbit, tMove, tRot, root, ptsGroup, dimGroup, view, labels, dirty = true, active = false;
let pickMode = null, pickA = null, snapHit = null, dragging = false, PIVOT = null, PROXY = null, drag = null, lastDrag = null;
const UNDO = [], REDO = [];

// ------------------------------------------------------------------ geometry (same as cad.py)

function localPath(e) {
  const p = e.params, t = e.type;
  if (t === 'rod') return p.pts && p.pts.length > 1 ? p.pts.map(q => q.slice()) : [[0, 0, 0], [+p.L, 0, 0]];
  if (t === 'fillet') {
    const R = +p.R, a = +p.angle * D2R, n = Math.max(2, Math.ceil(+p.angle / 10)), out = [];
    for (let i = 0; i <= n; i++) { const s = a * i / n; out.push([R * Math.sin(s), R * (1 - Math.cos(s)), 0]); }
    return out;
  }
  if (t === 'torus') {
    const R = +p.R, A = Math.min(360, +p.arc), n = Math.max(3, Math.ceil(A / 10)), out = [];
    for (let i = 0; i <= n; i++) { const s = A * D2R * i / n; out.push([R * Math.cos(s), R * Math.sin(s), 0]); }
    return out;
  }
  if (t === 'helix') {
    const R = +p.R, pitch = +p.pitch, turns = +p.turns, h = p.hand === 'L' ? -1 : 1, n = Math.max(4, Math.ceil(turns * 24)), out = [];
    for (let i = 0; i <= n; i++) { const s = 2 * Math.PI * turns * i / n; out.push([R * Math.cos(s), h * R * Math.sin(s), pitch * s / (2 * Math.PI)]); }
    return out;
  }
  if (t === 'meander') {
    const L = +p.L, A = +p.A, n = Math.max(1, Math.round(+p.n)), st = L / n, P = [[0, 0, 0]];
    for (let k = 0; k < n; k++) { const y1 = k % 2 === 0 ? A : 0; P.push([k * st, y1, 0], [(k + 1) * st, y1, 0]); }
    P.push([L, n % 2 === 0 ? 0 : A, 0]);
    return P.filter((q, i) => i === 0 || q.some((v, k) => v !== P[i - 1][k]));
  }
  if (t === 'custom') return localPaths(e)[0];
  return null;
}
function localPaths(e) {
  if (e.type !== 'custom') return [localPath(e)];
  const out = (e.params.paths || []).filter(P => P.length > 1).map(P => P.map(q => [+q[0], +q[1], q.length > 2 ? +q[2] : 0]));
  return out.length ? out : [[[0, 0, 0], [1, 0, 0]]];
}
// planar (pcb) layers, same as pcb.py
const loopOf = L => { const A = (L || []).map(q => [+q[0], +q[1]]); if (A.length > 1 && A[0][0] === A[A.length - 1][0] && A[0][1] === A[A.length - 1][1]) A.pop(); return A; };
function inLoop(P, x, y) { let c = false; for (let i = 0, j = P.length - 1; i < P.length; j = i++) { const a = P[i], b = P[j]; if ((a[1] > y) !== (b[1] > y) && x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]) c = !c; } return c; }
function drawnRegions(e, layer) { return (e.params[layer] || []).map(R => ({outer: loopOf(R.outer), holes: (R.holes || []).map(loopOf).filter(h => h.length >= 3)})).filter(r => r.outer.length >= 3); }
const inRegions = (regs, x, y) => regs.some(r => inLoop(r.outer, x, y) && !r.holes.some(h => inLoop(h, x, y)));
function boxOf(regs) { if (!regs.length) return null; let x0 = 1e18, y0 = 1e18, x1 = -1e18, y1 = -1e18; for (const r of regs) for (const [x, y] of r.outer) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); } return [x0, y0, x1, y1]; }
function pcbRegions(e, layer) {
  const p = e.params, regs = drawnRegions(e, layer);
  if (layer !== 'top' || !p.pix) return regs;
  const bb = boxOf(regs); if (!bb) return regs;
  const [nx, ny] = p.pix, dx = (bb[2] - bb[0]) / nx, dy = (bb[3] - bb[1]) / ny, out = [];
  for (let b = 0; b < ny; b++) for (let a = 0; a < nx; a++) {
    const cx = bb[0] + (a + .5) * dx, cy = bb[1] + (b + .5) * dy;
    if ((p.mask && p.mask.length === nx * ny && !p.mask[b * nx + a]) || !inRegions(regs, cx, cy)) continue;
    const x0 = bb[0] + a * dx, y0 = bb[1] + b * dy; out.push({outer: [[x0, y0], [x0 + dx, y0], [x0 + dx, y0 + dy], [x0, y0 + dy]], holes: []});
  }
  return out;
}
function pcbPoints(e) {
  const p = e.params, h = +p.h || 0, o = {center: [0, 0, h / 2]};
  for (const [layer, z, t] of [['top', h, 't'], ['gnd', 0, 'g']]) {
    const bb = boxOf(drawnRegions(e, layer)); if (!bb) continue;
    [[bb[0], bb[1]], [bb[2], bb[1]], [bb[2], bb[3]], [bb[0], bb[3]]].forEach(([x, y], i) => { o[t + (i + 1)] = [x, y, z]; });
    o[t + 'c'] = [(bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2, z];
  }
  if (p.feed) { o.feed = [+p.feed[0], +p.feed[1], 0]; o.feedt = [+p.feed[0], +p.feed[1], h]; }
  (p.vias || []).forEach((v, i) => { o['via' + (i + 1)] = [+v[0], +v[1], h]; });
  return o;
}
function localPoints(e) {
  if (e.type === 'pcb') return pcbPoints(e);
  if (e.type === 'custom') {
    const Ps = localPaths(e), o = {start: Ps[0][0], end: Ps[Ps.length - 1][Ps[Ps.length - 1].length - 1], center: [0, 0, 0]};
    for (const [k, v] of Object.entries(e.params.pts || {})) o[k] = [+v[0], +v[1], v.length > 2 ? +v[2] : 0];
    return o;
  }
  if (PATH.includes(e.type)) {
    const P = localPath(e), o = {start: P[0], end: P[P.length - 1]};
    if (e.type === 'torus') o.center = [0, 0, 0];
    else if (e.type === 'helix') o.center = [0, 0, P[P.length - 1][2] / 2];
    else o.center = P.length > 2 ? P[Math.floor(P.length / 2)] : P[0].map((v, k) => (v + P[1][k]) / 2);
    return o;
  }
  const W = +e.params.W, H = +e.params.H, T = e.type === 'block' ? +e.params.T : 0, o = {center: [0, 0, T / 2]};
  for (const [nm, sx, sy] of [['c1', -1, -1], ['c2', 1, -1], ['c3', 1, 1], ['c4', -1, 1], ['e1', 0, -1], ['e2', 1, 0], ['e3', 0, 1], ['e4', -1, 0]]) {
    o[nm] = [sx * W / 2, sy * H / 2, 0];
    if (e.type === 'block') o[nm + 't'] = [sx * W / 2, sy * H / 2, T];
  }
  return o;
}
const quatOf = e => new THREE.Quaternion().setFromEuler(new THREE.Euler(e.rot[0] * D2R, e.rot[1] * D2R, e.rot[2] * D2R, 'XYZ'));
function mtx(e) { return new THREE.Matrix4().compose(new THREE.Vector3(...e.pos), quatOf(e), new THREE.Vector3(1, 1, 1)); }
function setMtx(e, m) {
  const p = new THREE.Vector3(), q = new THREE.Quaternion(), s = new THREE.Vector3(); m.decompose(p, q, s);
  const eu = new THREE.Euler().setFromQuaternion(q, 'XYZ');
  e.pos = [p.x, p.y, p.z].map(v => Math.round(v * 1e5) / 1e5); e.rot = [eu.x, eu.y, eu.z].map(v => Math.round(v / D2R * 1e4) / 1e4);
}
function worldPt(e, name) {
  if (e === 'origin') return new THREE.Vector3();
  const l = localPoints(e)[name]; if (!l) return null;
  return new THREE.Vector3(...l).applyMatrix4(mtx(e));
}
const worldPaths = e => { const m = mtx(e); return localPaths(e).map(P => P.map(q => new THREE.Vector3(...q).applyMatrix4(m))); };
const elById = id => DOC.elements.find(e => e.id === id);
const ptOf = ref => ref.el === 'origin' ? new THREE.Vector3() : (elById(ref.el) ? worldPt(elById(ref.el), ref.pt) : null);
const matOf = name => CAD.mats.find(m => m.name === name) || {kind: 'conductor', color: '#888888'};
const sameRef = (a, b) => a.el === b.el && a.pt === b.pt;
function elDir(e) {   // what an angle dimension measures (same as cad.element_dir)
  const q = quatOf(e);
  if (PATH.includes(e.type)) {
    const Ps = localPaths(e), a = Ps[0][0], b = Ps[Ps.length - 1][Ps[Ps.length - 1].length - 1];
    const v = new THREE.Vector3(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    if (v.length() > 1e-6) return v.normalize().applyQuaternion(q);
  }
  return new THREE.Vector3(0, 0, 1).applyQuaternion(q);
}
function angleNow(d) {
  const a = elById(d.a.el), b = elById(d.b.el); if (!a || !b) return null;
  return Math.acos(Math.max(-1, Math.min(1, elDir(a).dot(elDir(b))))) / D2R;
}
function dimNow(d) {
  if (d.axis === 'angle') return angleNow(d);
  const A = ptOf(d.a), B = ptOf(d.b); if (!A || !B) return null;
  return d.axis === 'd' ? A.distanceTo(B) : B.getComponent('xyz'.indexOf(d.axis)) - A.getComponent('xyz'.indexOf(d.axis));
}

// ------------------------------------------------------------------ document, undo

export function newDoc() {
  return {name: 'cad', ground: 'none', z0: 50, swr: 2, step: true, engine: 'nec', elements: [], joins: [], dims: []};
}
function defaults(t) {
  const L = CAD.lam(), d = CAD.wire();
  const P = {rod: {L: r3(L / 4)}, fillet: {R: r3(Math.max(3 * d, L / 40)), angle: 90, bend_r: 0}, torus: {R: r3(L / (2 * Math.PI)), arc: 360},
             helix: {R: r3(L / (2 * Math.PI)), pitch: r3(L / 4), turns: 3, hand: 'R'}, meander: {L: r3(L / 6), A: r3(L / 30), n: 8}, plate: {W: r3(L / 2), H: r3(L / 2), mesh: 0, t: 1},
             block: {W: r3(L / 4), H: r3(L / 4), T: 1.6}, custom: {paths: [], pts: {}}, pcb: {top: [], gnd: [], vias: [], feed: null, h: 1.6, sub: 'FR4', via_d: 0.8, mesh: 0}}[t];
  const mat = t === 'block' ? (CAD.mats.find(m => m.kind === 'dielectric') || {}).name || 'FR4' : 'Copper';
  return {id: 'e' + (seq++), type: t, name: LABEL[t].split(' ')[0] + ' ' + (DOC.elements.filter(e => e.type === t).length + 1),
          pos: [0, 0, 0], rot: [0, 0, 0], params: P, mat, d, feed: false, feed_at: 0.5, fixed: false, evolve: {mode: 'off', dev: r3(L / 40), size: 15, rot: 0, pts: 4}};
}
function snapshot() { UNDO.push(JSON.stringify(DOC)); if (UNDO.length > 60) UNDO.shift(); REDO.length = 0; }
export function undo() { if (!UNDO.length) return hint('Nothing to undo.'); REDO.push(JSON.stringify(DOC)); restore(UNDO.pop()); hint('Undone.'); }
export function redo() { if (!REDO.length) return hint('Nothing to redo.'); UNDO.push(JSON.stringify(DOC)); restore(REDO.pop()); hint('Redone.'); }
function restore(json) { DOC = JSON.parse(json); CAD.doc = DOC; SEL = SEL.filter(id => elById(id)); PIVOT = null; rebuild(); panel(); CAD.onChange && CAD.onChange(DOC); }

export function addElement(t) {
  snapshot();
  const e = defaults(t);
  const s = SEL.length ? elById(SEL[SEL.length - 1]) : null;
  if (s && PATH.includes(s.type) && PATH.includes(t) && t !== 'custom') {        // chain new wires onto the end of the selected one
    const end = worldPt(s, 'end'); e.pos = [r3(end.x), r3(end.y), r3(end.z)]; e.rot = s.rot.slice();
    DOC.elements.push(e); DOC.joins.push({a: {el: s.id, pt: 'end'}, b: {el: e.id, pt: 'start'}});
  } else DOC.elements.push(e);
  select([e.id]); changed();
  return e;
}
export function template(name) {
  const L = CAD.lam(), d = CAD.wire(), r = r3;
  const el = (type, o) => { const e = Object.assign(defaults(type), o); e.evolve = Object.assign(defaults(type).evolve, o.evolve || {}); return e; };
  snapshot();
  DOC = newDoc(); seq = 1; dseq = 1; CAD.doc = DOC;
  const free = (dev, pts) => ({mode: 'free', dev: r(dev), pts});
  const shape = (size, dev = 0, rot = 0, planar = false) => ({mode: 'shape', size, dev: r(dev), rot, planar});
  if (name === 'dipole') DOC.elements = [el('rod', {name: 'Dipole', pos: [r(-0.24 * L), 0, 0], params: {L: r(0.48 * L)}, feed: true, evolve: free(L / 20, 5)})];
  else if (name === 'mono') { DOC.ground = 'perfect'; DOC.elements = [el('rod', {name: 'Monopole', rot: [0, -90, 0], params: {L: r(0.24 * L)}, feed: true, feed_at: 0, evolve: free(L / 20, 4)})]; }
  else if (name === 'gp') {
    DOC.elements = [el('rod', {name: 'Radiator', rot: [0, -90, 0], params: {L: r(0.24 * L)}, feed: true, feed_at: 0, evolve: shape(15)})];
    for (let k = 0; k < 4; k++) {
      const e = el('rod', {name: 'Radial ' + (k + 1), rot: [0, 30, 90 * k], params: {L: r(0.25 * L)}, evolve: shape(15, 0, 15)}); DOC.elements.push(e);
      DOC.joins.push({a: {el: DOC.elements[0].id, pt: 'start'}, b: {el: e.id, pt: 'start'}});
    }
  }
  else if (name === 'yagi') {
    DOC.elements = [
      el('rod', {name: 'Reflector', pos: [r(-0.2 * L), r(-0.25 * L), 0], rot: [0, 0, 90], params: {L: r(0.5 * L)}, evolve: shape(8)}),
      el('rod', {name: 'Driven', pos: [0, r(-0.235 * L), 0], rot: [0, 0, 90], params: {L: r(0.47 * L)}, feed: true, evolve: shape(8)}),
      el('rod', {name: 'Director', pos: [r(0.15 * L), r(-0.22 * L), 0], rot: [0, 0, 90], params: {L: r(0.44 * L)}, evolve: shape(8)})];
    const [a, b, c] = DOC.elements.map(e => e.id);
    DOC.dims = [{id: 'd' + (dseq++), a: {el: a, pt: 'center'}, b: {el: b, pt: 'center'}, axis: 'x', value: r(0.2 * L), drive: true, evo: r(0.06 * L)},
                {id: 'd' + (dseq++), a: {el: b, pt: 'center'}, b: {el: c, pt: 'center'}, axis: 'x', value: r(0.15 * L), drive: true, evo: r(0.06 * L)},
                {id: 'd' + (dseq++), a: {el: a, pt: 'center'}, b: {el: b, pt: 'center'}, axis: 'y', value: 0, drive: true},
                {id: 'd' + (dseq++), a: {el: b, pt: 'center'}, b: {el: c, pt: 'center'}, axis: 'y', value: 0, drive: true}];
  }
  else if (name === 'patch') {     // air-dielectric patch over a ground plate, probe fed, four parasitic meanders around it
    const P = 0.415 * L, h = 0.04 * L, g = 0.065 * L, mA = 0.03 * L, mL = 0.3 * L;
    DOC.elements = [
      el('plate', {name: 'Ground', params: {W: r(0.95 * L), H: r(0.95 * L), mesh: 0, t: 1}, fixed: true}),
      el('plate', {name: 'Patch', pos: [0, 0, r(h)], params: {W: r(P), H: r(0.9 * P), mesh: 0, t: 1}, evolve: shape(6, 0.03 * L, 0, true)}),
      el('rod', {name: 'Probe', pos: [r(0.3 * P), 0, 0], rot: [0, -90, 0], params: {L: r(h)}, feed: true, feed_at: 0, evolve: {mode: 'free', dev: r(0.06 * P), pts: 2, planar: false}}),
      el('meander', {name: 'Meander +x', pos: [r(P / 2 + g), r(-mL / 2), r(h)], rot: [0, 0, 90], params: {L: r(mL), A: r(mA), n: 8}, evolve: shape(15, 0.02 * L, 0, true)}),
      el('meander', {name: 'Meander -x', pos: [r(-P / 2 - g), r(mL / 2), r(h)], rot: [0, 0, -90], params: {L: r(mL), A: r(mA), n: 8}, evolve: shape(15, 0.02 * L, 0, true)}),
      el('meander', {name: 'Meander +y', pos: [r(mL / 2), r(0.45 * P + g), r(h)], rot: [0, 0, 180], params: {L: r(mL), A: r(mA), n: 8}, evolve: shape(15, 0.02 * L, 0, true)}),
      el('meander', {name: 'Meander -y', pos: [r(-mL / 2), r(-0.45 * P - g), r(h)], params: {L: r(mL), A: r(mA), n: 8}, evolve: shape(15, 0.02 * L, 0, true)})];
    DOC.swr = 2.5;
  }
  else if (name === 'pixel') {     // pixel patch: the evolution switches the cells of the patch on and off
    const P = 0.42 * L, h = 0.04 * L;
    DOC.elements = [
      el('plate', {name: 'Ground', params: {W: r(0.9 * L), H: r(0.9 * L), mesh: 0, t: 1}, fixed: true}),
      el('plate', {name: 'Pixel patch', pos: [0, 0, r(h)], params: {W: r(P), H: r(P), mesh: 0, t: 1}, evolve: {mode: 'pixels', pix: [8, 8]}}),
      el('rod', {name: 'Probe', pos: [r(0.25 * P), 0, 0], rot: [0, -90, 0], params: {L: r(h)}, feed: true, feed_at: 0})];
    DOC.swr = 2.5;
  }
  else if (name === 'fr4patch') {   // classic inset-free probe-fed patch on 1.6 mm FR4 (transmission-line model), for openEMS
    const f = 299792.458 / L, er = 4.4, hs = 1.6, c = 299792458e3 / (f * 1e6);
    const W = c / 2 * Math.sqrt(2 / (er + 1)), ee = (er + 1) / 2 + (er - 1) / 2 / Math.sqrt(1 + 12 * hs / W);
    const dL = 0.412 * hs * (ee + 0.3) * (W / hs + 0.264) / ((ee - 0.258) * (W / hs + 0.8)), Lp = 0.955 * (c / (2 * Math.sqrt(ee)) - 2 * dL);   // 0.955: transmission-line model vs openEMS
    const fr4 = (CAD.mats.find(m => m.name === 'FR4') || {}).name || 'FR4';
    DOC.elements = [
      el('plate', {name: 'Ground', params: {W: r(Lp + 12 * hs), H: r(W + 12 * hs), mesh: 0, t: 0.035}, fixed: true}),
      el('block', {name: 'FR4 board', params: {W: r(Lp + 12 * hs), H: r(W + 12 * hs), T: hs}, mat: fr4, fixed: true}),
      el('plate', {name: 'Patch', pos: [0, 0, hs], params: {W: r(Lp), H: r(W), mesh: 0, t: 0.035}, evolve: shape(4, 0, 0, true)}),
      el('rod', {name: 'Probe', pos: [r(-Lp * 0.17), 0, 0], rot: [0, -90, 0], params: {L: hs}, feed: true, feed_at: 0, d: 1.2, evolve: {mode: 'free', dev: r(Lp * 0.08), pts: 2}})];
    DOC.engine = 'openems'; DOC.swr = 2;
  }
  else if (name === 'loop') DOC.elements = [el('torus', {name: 'Loop', params: {R: r(L / (2 * Math.PI)), arc: 360}, feed: true, feed_at: 0, evolve: free(L / 30, 8)})];
  else if (name === 'helix') { DOC.ground = 'perfect'; DOC.elements = [
      el('rod', {name: 'Riser', pos: [r(L / (2 * Math.PI) * 1.05), 0, 0], rot: [0, -90, 0], params: {L: r(L / 40)}, feed: true, feed_at: 0}),
      el('helix', {name: 'Helix', pos: [0, 0, r(L / 40)], params: {R: r(L / (2 * Math.PI) * 1.05), pitch: r(0.23 * L), turns: 5, hand: 'R'}, evolve: shape(12)})];
    DOC.joins.push({a: {el: DOC.elements[0].id, pt: 'end'}, b: {el: DOC.elements[1].id, pt: 'start'}}); DOC.z0 = 140; }
  for (const e of DOC.elements) if (!e.d) e.d = d;
  SEL = []; PIVOT = null; solve(); rebuild(); fit(); panel(); CAD.onChange && CAD.onChange(DOC);
}
export function setDoc(doc) {
  DOC = doc || newDoc();
  for (const e of DOC.elements) { e.evolve = e.evolve || {mode: 'off'}; seq = Math.max(seq, (parseInt(String(e.id).slice(1)) || 0) + 1); }
  for (const d of DOC.dims) dseq = Math.max(dseq, (parseInt(String(d.id).slice(1)) || 0) + 1);
  CAD.doc = DOC; SEL = []; PIVOT = null; UNDO.length = 0; REDO.length = 0; rebuild(); fit(); panel();
}
export function getDoc() { return DOC; }
function changed(noSolve) {
  if (!noSolve) solve();
  rebuild(); panel(); CAD.onChange && CAD.onChange(DOC);
}

// ------------------------------------------------------------------ constraint solver (same as cad.solve)

export function solve(lock) {
  const locked = new Set(lock instanceof Set ? lock : lock ? [lock] : []);
  const free = e => e && !e.fixed && !locked.has(e.id);
  const cons = [];
  for (const j of DOC.joins) cons.push({a: j.a, b: j.b, axis: 'join', value: 0});
  for (const d of DOC.dims) if (d.drive !== false) cons.push(d);
  let worst = 0;
  for (let it = 0; it < 400; it++) {
    worst = 0;
    for (const c of cons) {
      if (c.axis === 'angle') {
        const ea = elById(c.a.el), eb = elById(c.b.el); if (!ea || !eb) continue;
        const err = c.value - angleNow(c); worst = Math.max(worst, Math.abs(err) * D2R * 10); if (Math.abs(err) < 1e-5) continue;
        const [mov, oth] = free(eb) ? [eb, ea] : free(ea) ? [ea, eb] : [null, null]; if (!mov) continue;
        const u = elDir(oth), v = elDir(mov); let n = u.clone().cross(v);
        if (n.length() < 1e-9) n = Math.abs(u.z) < 0.9 ? u.clone().cross(new THREE.Vector3(0, 0, 1)) : u.clone().cross(new THREE.Vector3(1, 0, 0));
        let piv = null;
        for (const j of DOC.joins) for (const [x, y] of [[j.a, j.b], [j.b, j.a]]) if (!piv && x.el === mov.id && y.el === oth.id) piv = ptOf(y);
        if (!piv) piv = worldPt(mov, 'center') || new THREE.Vector3(...mov.pos);
        const R = new THREE.Matrix4().makeRotationAxis(n.normalize(), err * D2R);
        const T = new THREE.Matrix4().makeTranslation(piv.x, piv.y, piv.z).multiply(R).multiply(new THREE.Matrix4().makeTranslation(-piv.x, -piv.y, -piv.z));
        setMtx(mov, T.multiply(mtx(mov)));
        continue;
      }
      const A = ptOf(c.a), B = ptOf(c.b); if (!A || !B) continue;
      const v = B.clone().sub(A); let corr;
      if (c.axis === 'join') corr = v;
      else if (c.axis === 'd') { const L = v.length(); if (L < 1e-9) continue; corr = v.multiplyScalar((L - c.value) / L); }
      else { const k = 'xyz'.indexOf(c.axis); corr = new THREE.Vector3(); corr.setComponent(k, v.getComponent(k) - c.value); }
      const err = corr.length(); worst = Math.max(worst, err); if (err < 1e-6) continue;
      const ea = c.a.el === 'origin' ? null : elById(c.a.el), eb = c.b.el === 'origin' ? null : elById(c.b.el);
      const wa = free(ea) ? 1 : 0, wb = free(eb) ? 1 : 0; if (!wa && !wb) continue;
      const s = 1 / (wa + wb);
      if (wa) ea.pos = ea.pos.map((p, k) => p + corr.getComponent(k) * s * wa);
      if (wb) eb.pos = eb.pos.map((p, k) => p - corr.getComponent(k) * s * wb);
    }
    if (worst < 1e-5) break;
  }
  for (const e of DOC.elements) e.pos = e.pos.map(v => Math.round(v * 1e5) / 1e5);
  CAD.conflict = worst > 1e-3 ? worst : 0;
  return worst;
}

// ------------------------------------------------------------------ 3D

const tubeMat = {}; const M = (c, o) => { const k = c + JSON.stringify(o || {}); return tubeMat[k] || (tubeMat[k] = new THREE.MeshStandardMaterial(Object.assign({color: c, roughness: .4, metalness: .35}, o || {}))); };
const CYL = new THREE.CylinderGeometry(1, 1, 1, 12), SPH = new THREE.SphereGeometry(1, 12, 8);
const SELCOL = 0x3fc6ff;

// ------------------------------------------------------------------ optional interface style for the 3D views (one accent on near-black, flat phosphor lines)
export function uiPalette() {
  const acc = new THREE.Color(getComputedStyle(document.documentElement).getPropertyValue('--acc').trim() || '#ffb000'), k = new THREE.Color(0, 0, 0);
  const mix = (a, b, t) => a.clone().lerp(b, t);
  return {acc, ink: mix(acc, new THREE.Color(1, 1, 1), .55), muted: mix(k, acc, .5), line: mix(k, acc, .26), dim: mix(k, acc, .12), bg: mix(new THREE.Color('#040506'), acc, .03)};
}
// replace the materials of a built group: wires flat accent, sheets translucent, edges accent, feed and selection bright
export function stylize(g, sel) {
  const P = uiPalette();
  g.traverse(o => {
    if (!o.material || o.userData.keepStyle) return;
    const old = o.material, c = old.color ? old.color.getHex() : 0;
    const feed = c === 0x111111 || c === 0x222222, sheet = o.geometry && !(o.geometry.type === 'CylinderGeometry' || o.geometry.type === 'SphereGeometry');
    let m;
    if (o.isLineSegments || o.isLine) m = new THREE.LineBasicMaterial({color: sel ? P.ink : P.acc, transparent: true, opacity: old.opacity < 1 ? Math.max(.25, old.opacity) : 1});
    else if (feed) m = new THREE.MeshBasicMaterial({color: 0xffffff});
    else if (sheet) m = new THREE.MeshBasicMaterial({color: sel ? P.ink : P.acc, transparent: true, opacity: old.transparent && old.opacity < .5 ? .05 : .11, side: THREE.DoubleSide, depthWrite: false});
    else m = new THREE.MeshBasicMaterial({color: sel ? P.ink : P.acc});
    if (old.userData && old.userData.own) old.dispose();
    m.userData = {own: true}; o.material = m;
  });
  return g;
}

function pixelRects(e) {
  const W = +e.params.W, H = +e.params.H, p = e.params;
  if (!p.mask || !p.pix) return [[-W / 2, -H / 2, W / 2, H / 2]];
  const [nx, ny] = p.pix, out = [];
  for (let b = 0; b < ny; b++) for (let a = 0; a < nx; a++) if (p.mask[b * nx + a]) out.push([-W / 2 + a * W / nx, -H / 2 + b * H / ny, -W / 2 + (a + 1) * W / nx, -H / 2 + (b + 1) * H / ny]);
  return out;
}
function elementMesh(e) {
  const g = new THREE.Group(); g.userData.id = e.id;
  const m = matOf(e.mat), sel = SEL.includes(e.id);
  const col = new THREE.Color(m.color || '#888');
  const mt = M(col.getHex(), sel ? {emissive: SELCOL, emissiveIntensity: .55} : null);
  if (PATH.includes(e.type)) {
    const r = Math.max(e.d / 2, sceneSize() / 400);
    localPaths(e).forEach((Q, pi) => {
      const P = Q.map(q => new THREE.Vector3(...q));
      for (let i = 0; i < P.length - 1; i++) {
        const a = P[i], b = P[i + 1], len = a.distanceTo(b); if (len <= 0) continue;
        const c = new THREE.Mesh(CYL, mt); c.scale.set(r, len, r); c.position.copy(a).add(b).multiplyScalar(.5);
        c.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), b.clone().sub(a).normalize()); c.userData = {id: e.id, path: pi, seg: i}; g.add(c);
        const s = new THREE.Mesh(SPH, mt); s.scale.setScalar(r); s.position.copy(b); s.userData.id = e.id; g.add(s);
      }
      const s0 = new THREE.Mesh(SPH, mt); s0.scale.setScalar(r); s0.position.copy(P[0]); s0.userData.id = e.id; g.add(s0);
      if (e.feed && pi === Math.min(e.feed_path || 0, localPaths(e).length - 1)) {   // feed marker: a black bead at the feed point
        const q = alongPath(Q, e.feed_at);
        const f = new THREE.Mesh(SPH, M(0x111111)); f.scale.setScalar(r * 2.6); f.position.set(...q); f.userData.id = e.id; g.add(f);
      }
    });
  } else if (e.type === 'pcb') {
    const p = e.params, h = +p.h || 0.1, cu = M(col.getHex(), Object.assign({side: THREE.DoubleSide}, sel ? {emissive: SELCOL, emissiveIntensity: .55} : {}));
    const gcol = M(0x8f99a3, Object.assign({side: THREE.DoubleSide}, sel ? {emissive: SELCOL, emissiveIntensity: .35} : {}));
    const shape = r => { const s = new THREE.Shape(r.outer.map(([x, y]) => new THREE.Vector2(x, y))); for (const H of r.holes) s.holes.push(new THREE.Path(H.map(([x, y]) => new THREE.Vector2(x, y)))); return s; };
    for (const [layer, z, mt_] of [['top', h, cu], ['gnd', 0, gcol]]) {
      const regs = pcbRegions(e, layer); if (!regs.length) continue;
      const geo = new THREE.ShapeGeometry(regs.map(shape));
      const m_ = new THREE.Mesh(geo, mt_); m_.position.z = z + (layer === 'top' ? 0.01 : -0.01); m_.userData.id = e.id; m_.userData.own = true; g.add(m_);
      const ed = new THREE.LineSegments(new THREE.EdgesGeometry(geo), new THREE.LineBasicMaterial({color: sel ? SELCOL : layer === 'top' ? 0x5a3412 : 0x3a434c}));
      ed.position.z = m_.position.z; ed.userData.own = true; g.add(ed);
    }
    const sm = CAD.mats.find(x => x.name === p.sub);
    const bb = boxOf(drawnRegions(e, 'gnd')) || boxOf(drawnRegions(e, 'top'));
    if (sm && bb) {
      const b = new THREE.Mesh(new THREE.BoxGeometry(bb[2] - bb[0], bb[3] - bb[1], h), new THREE.MeshStandardMaterial({color: new THREE.Color(sm.color || '#3f8f4f'), transparent: true, opacity: sel ? .45 : .3, roughness: .8, metalness: 0, depthWrite: false}));
      b.position.set((bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2, h / 2); b.userData.id = e.id; b.userData.own = true; b.material.userData = {own: true}; g.add(b);
    }
    const rv = Math.max((+p.via_d || 0.8) / 2, 0.05);
    for (const v of p.vias || []) { const c = new THREE.Mesh(CYL, M(0xd7dade)); c.scale.set(rv, h, rv); c.rotation.x = Math.PI / 2; c.position.set(+v[0], +v[1], h / 2); c.userData.id = e.id; g.add(c); }
    if (p.feed) { const c = new THREE.Mesh(CYL, M(0x111111)); c.scale.set(rv * 1.6, h * 1.6, rv * 1.6); c.rotation.x = Math.PI / 2; c.position.set(+p.feed[0], +p.feed[1], h / 2 - h * 0.3); c.userData.id = e.id; g.add(c); }
  } else {
    const W = +e.params.W, H = +e.params.H, T = e.type === 'block' ? +e.params.T : Math.max(+e.params.t || 0.5, sceneSize() / 600);
    const diel = m.kind === 'dielectric';
    const rects = e.type === 'plate' ? pixelRects(e) : [[-W / 2, -H / 2, W / 2, H / 2]];
    const mat = diel ? new THREE.MeshStandardMaterial({color: col, transparent: true, opacity: sel ? .55 : .35, roughness: .8, metalness: 0, depthWrite: false})
                     : M(col.getHex(), Object.assign({side: THREE.DoubleSide}, sel ? {emissive: SELCOL, emissiveIntensity: .55} : {}));
    for (const [x0, y0, x1, y1] of rects) {
      const b = new THREE.Mesh(new THREE.BoxGeometry(x1 - x0, y1 - y0, T), mat);
      b.position.set((x0 + x1) / 2, (y0 + y1) / 2, e.type === 'block' ? T / 2 : 0); b.userData.id = e.id; b.userData.own = true; g.add(b);
    }
    const ob = new THREE.BoxGeometry(W, H, T);
    const ed = new THREE.LineSegments(new THREE.EdgesGeometry(ob), new THREE.LineBasicMaterial({color: sel ? SELCOL : 0x6b7480})); ob.dispose();
    ed.position.z = e.type === 'block' ? T / 2 : 0; ed.userData.own = true; g.add(ed);
    if (diel) g.children[0].material.userData = {own: true};
  }
  g.position.set(...e.pos); g.rotation.set(e.rot[0] * D2R, e.rot[1] * D2R, e.rot[2] * D2R, 'XYZ');
  if (CAD.uiStyle) stylize(g, sel);
  return g;
}
function alongPath(Q, f) {
  const L = []; let tot = 0; for (let i = 0; i < Q.length - 1; i++) { L.push(Math.hypot(...Q[i + 1].map((v, k) => v - Q[i][k]))); tot += L[i]; }
  let s = f * tot, i = 0; while (i < L.length - 1 && s > L[i]) { s -= L[i]; i++; }
  const u = L[i] ? Math.min(1, s / L[i]) : 0;
  return Q[i].map((v, k) => v + u * (Q[i + 1][k] - v));
}
function fractionAt(e, wp) {   // nearest position along the element's wires to a world point -> {path, f}
  const inv = mtx(e).invert(), p = wp.clone().applyMatrix4(inv);
  let best = null;
  localPaths(e).forEach((Q, pi) => {
    let tot = 0; const acc = [0]; for (let i = 0; i < Q.length - 1; i++) { tot += Math.hypot(...Q[i + 1].map((v, k) => v - Q[i][k])); acc.push(tot); }
    for (let i = 0; i < Q.length - 1; i++) {
      const a = new THREE.Vector3(...Q[i]), b = new THREE.Vector3(...Q[i + 1]), ab = b.clone().sub(a), L2 = ab.lengthSq(); if (!L2) continue;
      const u = Math.max(0, Math.min(1, p.clone().sub(a).dot(ab) / L2)), d = a.clone().addScaledVector(ab, u).distanceTo(p);
      if (!best || d < best.d) best = {d, path: pi, f: tot ? (acc[i] + u * Math.sqrt(L2)) / tot : 0};
    }
  });
  return best;
}
function centreAt(e, wp) {   // the point on the element's wire axis (centre of its cross-section) nearest to a world point
  if (!PATH.includes(e.type)) return wp.clone();
  const inv = mtx(e).invert(), p = wp.clone().applyMatrix4(inv);
  let best = null;
  for (const Q of localPaths(e)) for (let i = 0; i < Q.length - 1; i++) {
    const a = new THREE.Vector3(...Q[i]), b = new THREE.Vector3(...Q[i + 1]), ab = b.clone().sub(a), L2 = ab.lengthSq(); if (!L2) continue;
    const u = Math.max(0, Math.min(1, p.clone().sub(a).dot(ab) / L2)), q = a.clone().addScaledVector(ab, u), d = q.distanceTo(p);
    if (!best || d < best.d) best = {d, q};
  }
  return best ? best.q.applyMatrix4(mtx(e)) : wp.clone();
}
function pivotFrom(ev, h) {   // a node within 30 px wins, otherwise the wire axis under the cursor
  const p = hitPoint(ev, 30);
  if (p) return {at: p.p.clone(), node: true};
  h = h || hitElement(ev); if (!h) return null;
  const e = elById(h.id); return {at: e ? centreAt(e, h.point) : h.point.clone(), node: false};
}
let _size = 0;
function sceneSize() {
  if (_size) return _size;
  let s = CAD.lam() / 4;
  for (const e of DOC.elements) for (const p of Object.values(localPoints(e))) s = Math.max(s, new THREE.Vector3(...p).applyMatrix4(mtx(e)).length());
  return (_size = s);
}
function disposeGroup(g) { g.traverse(o => { if (o.geometry && o.geometry !== CYL && o.geometry !== SPH) o.geometry.dispose(); if (o.material && (o.userData.own || (o.material.userData || {}).own)) o.material.dispose(); }); g.clear(); }

function rebuild() {
  if (!scene) return;
  _size = 0;
  disposeGroup(root);
  for (const e of DOC.elements) root.add(elementMesh(e));
  rebuildOverlays();
  placeGumball();
  GROUND.visible = DOC.ground === 'perfect';
  GROUND.scale.setScalar(sceneSize() * 1.4);
  dirty = true;
}
function rebuildOverlays() {   // snap points, dimension lines, labels (element meshes stay)
  disposeGroup(ptsGroup); disposeGroup(dimGroup);
  const ps = sceneSize() / 120;
  const showAll = pickMode === 'dimA' || pickMode === 'dimB' || pickMode === 'pivot';
  for (const e of DOC.elements) for (const [nm, l] of Object.entries(localPoints(e))) {
    const w = new THREE.Vector3(...l).applyMatrix4(mtx(e));
    const joined = DOC.joins.some(j => (j.a.el === e.id && j.a.pt === nm) || (j.b.el === e.id && j.b.pt === nm));
    const m = new THREE.Mesh(SPH, M(joined ? 0x1fbf7a : nm === 'center' ? 0x9b6cf0 : 0xff8a1f, {depthTest: false}));
    m.scale.setScalar(ps * (joined ? 1.2 : 1)); m.position.copy(w); m.renderOrder = 5; m.userData = {el: e.id, pt: nm};
    m.visible = showAll || SEL.includes(e.id) || joined;
    ptsGroup.add(m);
  }
  const o = new THREE.Mesh(SPH, M(0xe6e6e6, {depthTest: false})); o.scale.setScalar(ps * 1.4); o.userData = {el: 'origin', pt: 'origin'}; o.renderOrder = 5; ptsGroup.add(o);
  if (PIVOT) { const m = new THREE.Mesh(new THREE.OctahedronGeometry(ps * 2.2), new THREE.MeshBasicMaterial({color: 0xffd400, depthTest: false, wireframe: true})); m.position.copy(PIVOT); m.renderOrder = 8; m.userData.own = true; dimGroup.add(m); }
  for (const d of DOC.dims) {
    if (d.axis === 'angle') {
      const a = elById(d.a.el), b = elById(d.b.el); if (!a || !b) continue;
      const pa = worldPt(a, 'center'), pb = worldPt(b, 'center');
      const ln = new THREE.Line(new THREE.BufferGeometry().setFromPoints([pa, pb]), new THREE.LineDashedMaterial({color: 0xffb020, dashSize: sceneSize() / 60, gapSize: sceneSize() / 60, depthTest: false}));
      ln.computeLineDistances(); ln.renderOrder = 6; ln.userData.own = true; dimGroup.add(ln);
      continue;
    }
    const A = ptOf(d.a), B = ptOf(d.b); if (!A || !B) continue;
    const pts = [A];
    if (d.axis !== 'd') { const C = A.clone(); C.setComponent('xyz'.indexOf(d.axis), B.getComponent('xyz'.indexOf(d.axis))); pts.push(C); }
    pts.push(B);
    const ln = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), new THREE.LineDashedMaterial({color: 0x3fc6ff, dashSize: sceneSize() / 40, gapSize: sceneSize() / 80, depthTest: false}));
    ln.computeLineDistances(); ln.renderOrder = 6; ln.userData.own = true; dimGroup.add(ln);
  }
  if (snapHit) { const m = new THREE.Mesh(SPH, M(0x1fbf7a, {depthTest: false, transparent: true, opacity: .6})); m.scale.setScalar(ps * 2.4); m.position.copy(snapHit.p); m.renderOrder = 7; dimGroup.add(m); }
  buildLabels();
  dirty = true;
}

// ------------------------------------------------------------------ labels in the view (dimensions, parameters), edited in place

let LBL = [];   // {el: div, at: Vector3 (world), kind, ...}
function buildLabels() {
  if (!labels) return;
  labels.innerHTML = ''; LBL = [];
  const add = (html, at, cls, data) => { const d = document.createElement('div'); d.className = 'clab ' + cls; d.innerHTML = html; labels.appendChild(d); LBL.push(Object.assign({el: d, at}, data)); return d; };
  for (const d of DOC.dims) {
    const cur = dimNow(d); if (cur === null) continue;
    let at;
    if (d.axis === 'angle') at = worldPt(elById(d.a.el), 'center').lerp(worldPt(elById(d.b.el), 'center'), .5);
    else { const A = ptOf(d.a), B = ptOf(d.b), C = A.clone(); if (d.axis !== 'd') C.setComponent('xyz'.indexOf(d.axis), B.getComponent('xyz'.indexOf(d.axis))); at = d.axis === 'd' ? A.clone().lerp(B, .5) : C.clone().lerp(B, .5); }
    const bad = d.drive !== false && Math.abs(cur - d.value) > 0.01;
    const txt = d.axis === 'angle' ? `∠ ${cur.toFixed(1)}°` : `${d.axis === 'd' ? '' : d.axis + ' '}${cur.toFixed(2)}`;
    const div = add(`${txt}${d.evo ? ` <i>±${+d.evo}</i>` : ''}${d.drive === false ? ' <i>ref</i>' : ''}`, at, 'dim' + (bad ? ' bad' : '') + (d.axis === 'angle' ? ' ang' : ''), {kind: 'dim', dim: d.id});
    div.title = 'Double-click: type a value. Right-click: options.';
    div.ondblclick = ev => { ev.stopPropagation(); editDim(d, div); };
    div.oncontextmenu = ev => { ev.preventDefault(); ev.stopPropagation(); dimMenu(d, ev.clientX, ev.clientY); };
  }
  // parameters of the selected element, Fusion-style: double-click to type
  const e = SEL.length === 1 ? elById(SEL[0]) : null;
  if (e && !dragging) {
    const m = mtx(e), W = v => new THREE.Vector3(...v).applyMatrix4(m), p = e.params;
    const SPOTS = {
      rod: p.pts ? [] : [['L', W([+p.L / 2, 0, 0])]],
      fillet: () => { const P = localPath(e) || []; return [['R', W([0, +p.R, 0])]].concat(P.length ? [['angle', W(P[Math.floor(P.length / 2)])]] : []); },
      torus: [['R', W([+p.R / 2, 0, 0])], ['arc', W([-(+p.R), 0, 0])]],
      helix: [['R', W([+p.R / 2, 0, 0])], ['pitch', W([+p.R * 1.15, 0, +p.pitch / 2])], ['turns', W([-(+p.R) * 1.15, 0, +p.pitch * +p.turns / 2])]],
      meander: [['L', W([+p.L / 2, -(+p.A) * 0.6 - sceneSize() / 40, 0])], ['A', W([0, +p.A / 2, 0])], ['n', W([+p.L, +p.A * 1.4, 0])]],
      plate: [['W', W([0, -(+p.H) / 2, 0])], ['H', W([+p.W / 2, 0, 0])]],
      block: [['W', W([0, -(+p.H) / 2, 0])], ['H', W([+p.W / 2, 0, 0])], ['T', W([+p.W / 2, +p.H / 2, +p.T / 2])]],
      custom: [],
      pcb: () => { const bb = boxOf(drawnRegions(e, 'top')); return bb ? [['h', W([bb[2], bb[1], +p.h / 2])]] : []; },
    }[e.type];
    const spots = (typeof SPOTS === 'function' ? SPOTS() : SPOTS) || [];
    for (const [k, at] of spots) {
      const def = PARAMS[e.type].find(q => q[0] === k); if (!def) continue;
      const v = p[k];
      const div = add(`${k === 'turns' ? 'turns ' : k === 'n' ? 'n ' : k === 'angle' || k === 'arc' ? '' : k + ' '}${typeof v === 'number' ? +(+v).toFixed(2) : esc(v)}${def[2] === '°' ? '°' : ''}`, at, 'par', {kind: 'par', key: k});
      div.title = def[1] + ': double-click to type a new value';
      div.ondblclick = ev => { ev.stopPropagation(); editParam(e, k, div); };
    }
  }
  placeLabels();
}
function placeLabels() {
  if (!labels || !LBL.length) return;
  const w = view.clientWidth, h = view.clientHeight;
  for (const L of LBL) {
    const p = L.at.clone().project(camera);
    const vis = p.z < 1 && p.z > -1;
    L.el.style.display = vis ? '' : 'none';
    if (vis) L.el.style.transform = `translate(${((p.x + 1) / 2 * w).toFixed(1)}px, ${((1 - p.y) / 2 * h).toFixed(1)}px) translate(-50%, -50%)`;
  }
}
function inlineInput(div, value, onOk, extra = '') {
  div.classList.add('edit');
  div.innerHTML = `<input type="text" value="${value}" spellcheck="false">${extra}`;
  const inp = div.querySelector('input'); inp.focus(); inp.select();
  const done = ok => { if (div._done) return; div._done = true; if (ok) { const v = evalNum(inp.value); if (v !== null) onOk(v, div); else { hint('Not a number: ' + inp.value); rebuildOverlays(); } } else rebuildOverlays(); };
  inp.onkeydown = ev => { ev.stopPropagation(); if (ev.key === 'Enter') done(true); else if (ev.key === 'Escape') done(false); };
  inp.onblur = () => setTimeout(() => done(true), 120);
  return inp;
}
function evalNum(s) {   // numbers and simple arithmetic like "lambda/4" or "12.5*2"
  const t = String(s).trim().replace(/,/g, '.').replace(/mm|°/g, '').replace(/λ|lambda|lam/gi, '(' + CAD.lam() + ')');
  if (!/^[\d\s.+\-*/()eE]+$/.test(t) || !t) return null;
  try { const v = Function('"use strict";return (' + t + ')')(); return Number.isFinite(v) ? v : null; } catch (e) { return null; }
}
function editDim(d, div) {
  const ax = d.axis === 'angle' ? '' : `<span class="axb">${['d', 'x', 'y', 'z'].map(a => `<button data-ax="${a}" class="${a === d.axis ? 'on' : ''}">${a === 'd' ? '↔' : a}</button>`).join('')}</span>`;
  inlineInput(div, r3(dimNow(d)), v => { snapshot(); d.value = d.axis === 'angle' ? Math.max(0, Math.min(180, v)) : v; d.drive = true; changed(); }, ax);
  for (const b of div.querySelectorAll('[data-ax]')) b.onmousedown = ev => { ev.preventDefault(); snapshot(); d.axis = b.dataset.ax; d.value = r3(dimNow(d)); div._done = true; changed(); setTimeout(() => { const L = LBL.find(x => x.dim === d.id); if (L) editDim(d, L.el); }, 0); };
}
function editParam(e, k, div) {
  const def = PARAMS[e.type].find(q => q[0] === k);
  if (def[3]) { snapshot(); e.params[k] = def[3][(def[3].indexOf(e.params[k]) + 1) % def[3].length]; solve(e.id); changed(true); return; }
  inlineInput(div, e.params[k], v => { snapshot(); e.params[k] = k === 'n' ? Math.max(1, Math.round(v)) : Math.max(k === 'mesh' ? 0 : 0.001, v); solve(e.id); changed(true); });
}

// ------------------------------------------------------------------ scene, gumball with a pivot

let GROUND, GRID, AX;
export function init() {
  view = $('cadview');
  renderer = new THREE.WebGLRenderer({antialias: true});
  renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
  const bg = () => getComputedStyle(document.documentElement).getPropertyValue('--view').trim() || '#0d1014';
  renderer.setClearColor(new THREE.Color(bg()), 1);
  CAD.refreshTheme = () => {
    const P = uiPalette();
    renderer.setClearColor(CAD.uiStyle ? P.bg : new THREE.Color(bg()), 1);
    if (GRID) { const cols = GRID.geometry.attributes.color; if (cols) { const a = CAD.uiStyle ? P.line : new THREE.Color(0x3a424b), b = CAD.uiStyle ? P.muted : new THREE.Color(0x56606b); const pos = GRID.geometry.attributes.position; for (let i = 0; i < cols.count; i += 2) { const ctr = (Math.abs(pos.getX(i)) < 1e-9 && Math.abs(pos.getX(i + 1)) < 1e-9) || (Math.abs(pos.getZ(i)) < 1e-9 && Math.abs(pos.getZ(i + 1)) < 1e-9); const c = ctr ? b : a; cols.setXYZ(i, c.r, c.g, c.b); cols.setXYZ(i + 1, c.r, c.g, c.b); } cols.needsUpdate = true; } }
    if (GROUND) GROUND.material.color.set(CAD.uiStyle ? P.dim : new THREE.Color(0x9aa3ad));
    if (DOC) rebuild();
    dirty = true;
  };
  view.prepend(renderer.domElement);
  labels = document.createElement('div'); labels.id = 'cadlabels'; view.appendChild(labels);
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(40, 1, 0.1, 1e6); camera.up.set(0, 0, 1); camera.position.set(300, -400, 260);
  scene.add(new THREE.AmbientLight(0xffffff, 1.6));
  const sun = new THREE.DirectionalLight(0xffffff, 2.2); sun.position.set(1, -1, 2); scene.add(sun);
  root = new THREE.Group(); ptsGroup = new THREE.Group(); dimGroup = new THREE.Group(); scene.add(root, ptsGroup, dimGroup);
  GRID = new THREE.GridHelper(1, 20, 0x56606b, 0x3a424b); GRID.rotation.x = Math.PI / 2; GRID.material.transparent = true; GRID.material.opacity = .5; scene.add(GRID);
  AX = new THREE.AxesHelper(1); scene.add(AX);
  GROUND = new THREE.Mesh(new THREE.CircleGeometry(1, 64), new THREE.MeshStandardMaterial({color: 0x9aa3ad, transparent: true, opacity: .16, side: THREE.DoubleSide, depthWrite: false}));
  scene.add(GROUND);
  if (CAD.uiStyle) setTimeout(() => CAD.refreshTheme(), 0);
  orbit = new OrbitControls(camera, renderer.domElement); orbit.enableDamping = true; orbit.addEventListener('change', () => { dirty = true; });
  PROXY = new THREE.Object3D(); scene.add(PROXY);
  tMove = new TransformControls(camera, renderer.domElement); tMove.setMode('translate'); tMove.setSize(0.9);
  tRot = new TransformControls(camera, renderer.domElement); tRot.setMode('rotate'); tRot.setSize(1.25);
  for (const [t, other] of [[tMove, tRot], [tRot, tMove]]) {
    scene.add(t);
    t.addEventListener('change', () => { dirty = true; });
    t.addEventListener('mouseDown', () => startDrag(t));
    t.addEventListener('dragging-changed', ev => { orbit.enabled = !ev.value; other.enabled = !ev.value; if (!ev.value) endDrag(); });
    t.addEventListener('objectChange', () => onDrag(t));
  }
  gumball('both');
  new ResizeObserver(resize).observe(view);
  renderer.domElement.addEventListener('pointerdown', ev => { downAt = [ev.clientX, ev.clientY]; hideMenu(); });
  renderer.domElement.addEventListener('pointerup', onClick);
  let ghost = null, gT = 0;   // live preview of where the pivot will land while picking it
  renderer.domElement.addEventListener('pointermove', ev => {
    if (pickMode !== 'pivot' || ev.buttons) { if (ghost && ghost.visible) { ghost.visible = false; dirty = true; } return; }
    const now = performance.now(); if (now - gT < 30) return; gT = now;
    const pv = pivotFrom(ev);
    if (!ghost) { ghost = new THREE.Mesh(new THREE.OctahedronGeometry(1), new THREE.MeshBasicMaterial({color: 0xffd400, depthTest: false, wireframe: true, transparent: true, opacity: .6})); ghost.renderOrder = 9; scene.add(ghost); }
    ghost.visible = !!pv;
    if (pv) { ghost.position.copy(pv.at); ghost.scale.setScalar(sceneSize() / (pv.node ? 30 : 45)); ghost.material.color.setHex(pv.node ? 0x3fc6ff : 0xffd400); }
    dirty = true;
  });
  renderer.domElement.addEventListener('dblclick', onDbl);
  renderer.domElement.addEventListener('contextmenu', onContext);
  window.addEventListener('keydown', onKey);
  document.addEventListener('pointerdown', ev => { if (MENU && !MENU.contains(ev.target)) hideMenu(); });
  (function loop() {
    requestAnimationFrame(loop);
    if (!active) return;
    if (orbit.update()) dirty = true;
    if (!dirty || !view.clientWidth) return;
    resize(); dirty = false; renderer.render(scene, camera); placeLabels(); placeHud();
  })();
  bindPanel();
}
export function setActive(a) { active = a; if (a) { resize(); dirty = true; } }
function resize() {
  const w = view.clientWidth, h = view.clientHeight; if (!w || !h) return;
  const c = renderer.domElement;
  if (c.width !== Math.round(w * renderer.getPixelRatio()) || c.height !== Math.round(h * renderer.getPixelRatio())) renderer.setSize(w, h);
  camera.aspect = w / h; camera.updateProjectionMatrix(); dirty = true;
}
export function fit() {
  if (!scene) return;
  _size = 0; const s = sceneSize();
  const box = new THREE.Box3(); for (const e of DOC.elements) for (const p of Object.values(localPoints(e))) box.expandByPoint(new THREE.Vector3(...p).applyMatrix4(mtx(e)));
  const c = box.isEmpty() ? new THREE.Vector3() : box.getCenter(new THREE.Vector3());
  const size = box.isEmpty() ? s : Math.max(s / 4, ...box.getSize(new THREE.Vector3()).toArray());
  camera.near = size / 300; camera.far = size * 300; camera.updateProjectionMatrix();
  orbit.target.copy(c); camera.position.copy(c).add(new THREE.Vector3(size * 1.1, -size * 1.4, size * .95)); orbit.update();
  GRID.scale.setScalar(Math.max(size, CAD.lam() / 4) * 3); AX.scale.setScalar(Math.max(size, 10) * .6);
  dirty = true;
}
export function viewFrom(dir) {   // standard views: top, front, right, iso
  const t = orbit.target.clone(), d = camera.position.distanceTo(t);
  const v = {top: [0, -0.0001, 1], front: [0, -1, 0], right: [1, 0, 0], iso: [1.1, -1.4, .95]}[dir]; if (!v) return;
  camera.position.copy(t).add(new THREE.Vector3(...v).normalize().multiplyScalar(d)); orbit.update(); dirty = true;
}
export function gumball(mode) {
  tMove.visible = tMove.enabled = mode !== 'rotate';
  tRot.visible = tRot.enabled = mode !== 'move';
  for (const b of document.querySelectorAll('[data-gb]')) b.classList.toggle('on', b.dataset.gb === mode);
  CAD.gb = mode; placeGumball(); dirty = true;
}
function snapSettings() {
  const g = Number($('cadgrid').value) || 0, a = Number($('cadang').value) || 0;
  tMove.setTranslationSnap(g > 0 ? g : null); tRot.setRotationSnap(a > 0 ? a * D2R : null);
  const sp = $('cadspace').value; tMove.setSpace(sp); tRot.setSpace(sp); placeGumball();
}
function selCenter() {
  const box = new THREE.Box3();
  for (const id of SEL) { const e = elById(id); if (e) box.expandByPoint(worldPt(e, 'center')); }
  return SEL.length === 1 ? new THREE.Vector3(...elById(SEL[0]).pos) : box.getCenter(new THREE.Vector3());
}
function placeGumball() {
  if (!tMove) return;
  const sel = SEL.map(elById).filter(Boolean);
  if (!sel.length || dragging) { if (!dragging) { tMove.detach(); tRot.detach(); } return; }
  PROXY.position.copy(PIVOT || selCenter());
  PROXY.quaternion.copy(sel.length === 1 && $('cadspace') && $('cadspace').value === 'local' ? quatOf(sel[0]) : new THREE.Quaternion());
  PROXY.updateMatrixWorld(true);
  tMove.attach(PROXY); tRot.attach(PROXY);
  dirty = true;
}
function startDrag(t) {
  snapshot();
  dragging = true;
  drag = {t, axis: t.axis, mode: t.mode, P0: PROXY.matrixWorld.clone(), pos0: PROXY.position.clone(), q0: PROXY.quaternion.clone(),
          M0: new Map(SEL.map(id => [id, mtx(elById(id))])), others: new Map(DOC.elements.filter(e => !SEL.includes(e.id)).map(e => [e.id, e.pos.slice()]))};
  lastDrag = null;
}
function applyDelta(D) {   // D: world matrix change of the gumball since the drag started
  for (const [id, M0] of drag.M0) setMtx(elById(id), D.clone().multiply(M0));
  for (const [id, p] of drag.others) elById(id).pos = p.slice();   // neighbours start from where they were
}
function onDrag(t) {
  if (!drag) return;
  PROXY.updateMatrixWorld(true);
  const D = PROXY.matrixWorld.clone().multiply(drag.P0.clone().invert());
  applyDelta(D);
  snapHit = null;
  if (t === tMove && $('cadsnap').checked) {         // snap my points to the nearest point of another element
    const rad = sceneSize() / 15; let best = null;
    for (const id of SEL) {
      const e = elById(id), mine = Object.keys(localPoints(e)).map(n => [n, worldPt(e, n)]);
      for (const o of DOC.elements) {
        if (SEL.includes(o.id)) continue;
        for (const on of Object.keys(localPoints(o))) {
          const q = worldPt(o, on);
          for (const [mn, mp] of mine) { const d = mp.distanceTo(q); if (d < rad && (!best || d < best.d)) best = {d, id, mn, el: o.id, pt: on, p: q, delta: q.clone().sub(mp)}; }
        }
      }
    }
    if (best) { PROXY.position.add(best.delta); PROXY.updateMatrixWorld(true); applyDelta(PROXY.matrixWorld.clone().multiply(drag.P0.clone().invert())); snapHit = best; }
  }
  solve(new Set(SEL));                     // joined neighbours follow
  for (const c of root.children) { const o = elById(c.userData.id); c.position.set(...o.pos); c.rotation.set(o.rot[0] * D2R, o.rot[1] * D2R, o.rot[2] * D2R, 'XYZ'); }
  if (PIVOT) PIVOT.copy(PROXY.position);
  lastDrag = dragState();
  rebuildOverlays(); propsLive();
}
function dragState() {   // the change since the drag started, in gumball axes
  const q0i = drag.q0.clone().invert();
  const dp = PROXY.position.clone().sub(drag.pos0).applyQuaternion(t_space() === 'local' ? q0i : new THREE.Quaternion());
  const dq = PROXY.quaternion.clone().multiply(q0i);
  let ang = 2 * Math.acos(Math.min(1, Math.abs(dq.w))) / D2R; if (dq.w < 0) ang = -ang;
  const ax = new THREE.Vector3(dq.x, dq.y, dq.z); const sgn = ax.length() > 1e-9 ? Math.sign(ax.dot(axisVec(drag.axis)) || 1) : 1;
  return {mode: drag.mode, axis: drag.axis, dp, ang: Math.abs(ang) * sgn};
}
const t_space = () => ($('cadspace') ? $('cadspace').value : 'world');
function axisVec(a) {   // gumball axis in world coordinates at the drag start
  const v = new THREE.Vector3(a === 'X' ? 1 : 0, a === 'Y' ? 1 : 0, a === 'Z' ? 1 : 0);
  return t_space() === 'local' && drag ? v.applyQuaternion(drag.q0) : v;
}
function endDrag() {
  dragging = false;
  if (!drag) return;
  if (snapHit && $('cadautojoin').checked) {
    const j = {a: {el: snapHit.el, pt: snapHit.pt}, b: {el: snapHit.id, pt: snapHit.mn}};
    if (!DOC.joins.some(x => sameRef(x.a, j.a) && sameRef(x.b, j.b) || sameRef(x.a, j.b) && sameRef(x.b, j.a))) DOC.joins.push(j);
    hint(`Joined ${elById(snapHit.id).name} (${snapHit.mn}) to ${elById(snapHit.el).name} (${snapHit.pt}).`);
  }
  snapHit = null;
  solve(new Set(SEL));
  rebuild(); panel(); CAD.onChange && CAD.onChange(DOC);
  if (lastDrag && (lastDrag.mode === 'translate' ? lastDrag.dp.length() > 1e-9 : Math.abs(lastDrag.ang) > 1e-6)) showHud(); else drag = null;
}

// numeric box at the gumball after a drag: type the exact distance or angle (Fusion style)
let HUD = null;
function showHud() {
  hideHud();
  const s = lastDrag, single = /^[XYZ]$/.test(s.axis || '');
  HUD = document.createElement('div'); HUD.className = 'cadhud';
  if (s.mode === 'translate') {
    const ax = single ? [s.axis] : ['X', 'Y', 'Z'];
    HUD.innerHTML = ax.map(a => `<label>Δ${a.toLowerCase()} <input data-k="${a}" value="${r3(s.dp.getComponent('XYZ'.indexOf(a)))}"></label>`).join('') + '<span>mm · Enter applies</span>';
  } else if (single) HUD.innerHTML = `<label>∠${s.axis.toLowerCase()} <input data-k="A" value="${r3(s.ang)}"></label><span>° · Enter applies</span>`;
  else { drag = null; HUD = null; return; }
  view.appendChild(HUD);
  for (const inp of HUD.querySelectorAll('input')) {
    inp.onkeydown = ev => { ev.stopPropagation(); if (ev.key === 'Enter') applyHud(); else if (ev.key === 'Escape') hideHud(); };
    inp.onfocus = () => inp.select();
  }
  placeHud(); dirty = true;
}
function hideHud() { if (HUD) HUD.remove(); HUD = null; }
function placeHud() {
  if (!HUD || !PROXY) return;
  const p = PROXY.position.clone().project(camera), w = view.clientWidth, h = view.clientHeight;
  HUD.style.left = Math.min(w - 230, Math.max(4, (p.x + 1) / 2 * w + 40)) + 'px'; HUD.style.top = Math.min(h - 60, Math.max(4, (1 - p.y) / 2 * h + 30)) + 'px';
}
function applyHud() {
  if (!drag || !HUD) return hideHud();
  const s = lastDrag, val = k => { const i = HUD.querySelector(`[data-k="${k}"]`); return i ? evalNum(i.value) : null; };
  let D;
  if (s.mode === 'translate') {
    const v = new THREE.Vector3(...['X', 'Y', 'Z'].map(a => { const x = val(a); return x === null ? (HUD.querySelector(`[data-k="${a}"]`) ? 0 : (/^[XYZ]$/.test(s.axis) ? 0 : s.dp.getComponent('XYZ'.indexOf(a)))) : x; }));
    if (t_space() === 'local') v.applyQuaternion(drag.q0);
    D = new THREE.Matrix4().makeTranslation(v.x, v.y, v.z);
  } else {
    const a = val('A'); if (a === null) return hint('Not a number.');
    const c = drag.pos0;
    D = new THREE.Matrix4().makeTranslation(c.x, c.y, c.z).multiply(new THREE.Matrix4().makeRotationAxis(axisVec(s.axis).normalize(), a * D2R)).multiply(new THREE.Matrix4().makeTranslation(-c.x, -c.y, -c.z));
  }
  applyDelta(D);
  if (PIVOT) PIVOT.copy(drag.pos0).applyMatrix4(D);
  solve(new Set(SEL));
  hideHud(); drag = null;
  rebuild(); panel(); CAD.onChange && CAD.onChange(DOC);
  hint('Moved exactly as typed.');
}

// ------------------------------------------------------------------ picking, menus, keys

let downAt = [0, 0];
function raycast(ev, objs) {
  const r = renderer.domElement.getBoundingClientRect();
  const v = new THREE.Vector2((ev.clientX - r.left) / r.width * 2 - 1, -(ev.clientY - r.top) / r.height * 2 + 1);
  const rc = new THREE.Raycaster(); rc.setFromCamera(v, camera);
  rc.params.Line = {threshold: sceneSize() / 100};
  return rc.intersectObjects(objs, true);
}
function hitElement(ev) {
  const h = raycast(ev, root.children)[0]; if (!h) return null;
  let o = h.object; while (o && !o.userData.id) o = o.parent;
  return o ? {id: o.userData.id, point: h.point} : null;
}
function hitPoint(ev, rad = 14) {
  const pts = ptsGroup.children.filter(m => m.visible);
  const r = renderer.domElement.getBoundingClientRect();
  let best = null;   // screen-space nearest within 14 px: much easier to hit than the tiny spheres
  for (const m of pts) {
    const p = m.position.clone().project(camera), x = r.left + (p.x + 1) / 2 * r.width, y = r.top + (1 - p.y) / 2 * r.height, d = Math.hypot(x - ev.clientX, y - ev.clientY);
    if (d < rad && (!best || d < best.d)) best = {d, ref: {el: m.userData.el, pt: m.userData.pt}, p: m.position.clone()};
  }
  return best;
}
function onClick(ev) {
  if (ev.button !== 0) return;
  if (Math.hypot(ev.clientX - downAt[0], ev.clientY - downAt[1]) > 4 || dragging || tMove.dragging || tRot.dragging) return;
  if (!pickMode && ((tMove.enabled && tMove.axis) || (tRot.enabled && tRot.axis))) return;   // a click on the gumball
  if (pickMode === 'pivot') {
    const pv = pivotFrom(ev);
    if (!pv) return hint('Click a point (dots), a spot on an element, or the origin. Esc cancels.');
    PIVOT = pv.at; pickMode = null;
    hint((pv.node ? 'Pivot snapped to the node.' : 'Pivot set on the wire axis.') + ' The gumball turns and moves the selection around it. Shift+P puts it back.');
    rebuildOverlays(); placeGumball(); return;
  }
  if (pickMode === 'dimA' || pickMode === 'dimB') {
    const p = hitPoint(ev);
    if (!p) return hint('Click one of the points (dots, the light one is the origin).');
    if (pickMode === 'dimA') { pickA = p.ref; pickMode = 'dimB'; hint('Now the second point.'); return; }
    if (sameRef(p.ref, pickA)) return;
    snapshot();
    const d = {id: 'd' + (dseq++), a: pickA, b: p.ref, axis: 'd', value: 0, drive: false};
    d.value = r3(dimNow(d)); DOC.dims.push(d);
    pickMode = null; pickA = null; changed(true);
    hint('Type the distance and Enter (↔ / x / y / z choose what is measured). Esc keeps it as a reference.');
    const L = LBL.find(x => x.dim === d.id); if (L) editDim(d, L.el);
    return;
  }
  if (pickMode === 'angA' || pickMode === 'angB') {
    const h = hitElement(ev); if (!h) return hint('Click an element.');
    if (pickMode === 'angA') { pickA = h.id; pickMode = 'angB'; hint('Now the second element.'); select([h.id]); return; }
    if (h.id === pickA) return;
    pickMode = null; addAngle(pickA, h.id); pickA = null;
    return;
  }
  const h = hitElement(ev);
  if (!h) { if (!ev.shiftKey && !ev.ctrlKey) select([]); return; }
  if (ev.shiftKey || ev.ctrlKey) select(SEL.includes(h.id) ? SEL.filter(x => x !== h.id) : [...SEL, h.id]);
  else select([h.id]);
}
function addAngle(a, b) {
  snapshot();
  const d = {id: 'd' + (dseq++), a: {el: a}, b: {el: b}, axis: 'angle', value: 0, drive: false};
  d.value = r3(angleNow(d)); DOC.dims.push(d);
  select([b]);
  hint('Type the angle between the two elements and Enter (Esc keeps it as a reference).');
  const L = LBL.find(x => x.dim === d.id); if (L) editDim(d, L.el);
}
function onDbl(ev) {
  const h = hitElement(ev); if (!h) return;
  const e = elById(h.id);
  if ((e.type === 'custom' || e.type === 'pcb') && CAD.editSketch) CAD.editSketch(e);
}
export function select(ids) {
  const same = ids.length === SEL.length && ids.every(i => SEL.includes(i));
  SEL = ids; if (!same) { PIVOT = null; hideHud(); drag = null; }
  rebuild(); panel();
}

let MENU = null;
function hideMenu() { if (MENU) MENU.remove(); MENU = null; }
function showMenu(x, y, items) {
  hideMenu();
  MENU = document.createElement('div'); MENU.className = 'cadmenu';
  MENU.innerHTML = items.map((it, i) => it === '-' ? '<hr>' : `<button data-i="${i}"${it.off ? ' disabled' : ''}>${esc(it.label)}${it.key ? `<kbd>${it.key}</kbd>` : ''}</button>`).join('');
  document.body.appendChild(MENU);
  const r = MENU.getBoundingClientRect();
  MENU.style.left = Math.min(x, innerWidth - r.width - 6) + 'px'; MENU.style.top = Math.min(y, innerHeight - r.height - 6) + 'px';
  for (const b of MENU.querySelectorAll('button')) b.onclick = () => { const it = items[+b.dataset.i]; hideMenu(); it.do(); };
}
function onContext(ev) {
  ev.preventDefault();
  if (Math.hypot(ev.clientX - downAt[0], ev.clientY - downAt[1]) > 4) return;   // that was a pan
  const h = hitElement(ev);
  if (h && !SEL.includes(h.id)) select([h.id]);
  const e = h ? elById(h.id) : null, two = SEL.length === 2, isPath = e && PATH.includes(e.type);
  const items = [];
  if (e) {
    const pv = pivotFrom(ev, h); items.push({label: pv && pv.node ? 'Set pivot on this node' : 'Set pivot here', key: 'P', do: () => { PIVOT = pv.at; rebuildOverlays(); placeGumball(); hint('Pivot set on the element.'); }});
    if (PIVOT) items.push({label: 'Pivot back to the element', key: '⇧P', do: () => { PIVOT = null; rebuildOverlays(); placeGumball(); }});
    if (isPath && matOf(e.mat).kind === 'conductor') items.push({label: 'Feed point here', do: () => { const f = fractionAt(e, h.point); snapshot(); e.feed = true; e.feed_at = r3(f.f); e.feed_path = f.path; changed(true); }});
    if (e.type === 'rod') items.push({label: 'Split here', do: () => split(fractionAt(e, h.point).f)});
    if (e.type === 'custom' || e.type === 'pcb') items.push({label: 'Edit drawing', do: () => CAD.editSketch && CAD.editSketch(e)});
    if (e.type === 'pcb') items.push({label: 'Export layout as SVG', do: () => exportSvg(e)});
    items.push('-');
  }
  if (two) {
    items.push({label: 'Fillet the joint…', do: () => filletPrompt(ev.clientX, ev.clientY)});
    items.push({label: 'Join nearest points', key: 'J', do: joinSel});
    items.push({label: 'Angle between these two', do: () => addAngle(SEL[0], SEL[1])});
  }
  if (SEL.length) {
    items.push({label: 'Disconnect', do: disconnect});
    items.push({label: 'Duplicate', key: 'Ctrl D', do: duplicate});
    items.push({label: e && e.fixed ? 'Unfix' : 'Fix in place', do: () => { snapshot(); for (const id of SEL) elById(id).fixed = !(e && e.fixed); changed(true); }});
    items.push({label: 'Delete', key: 'Del', do: delSel});
  } else {
    items.push({label: 'Add dimension', key: 'D', do: () => startPick('dimA')}, {label: 'Add angle', key: 'A', do: () => startPick('angA')}, {label: 'Fit view', key: 'F', do: fit});
  }
  showMenu(ev.clientX, ev.clientY, items);
}
function dimMenu(d, x, y) {
  showMenu(x, y, [
    {label: 'Type a value', do: () => { const L = LBL.find(l => l.dim === d.id); if (L) editDim(d, L.el); }},
    {label: d.drive === false ? 'Make driving (the solver keeps it)' : 'Make reference (only shows the value)', do: () => { snapshot(); d.drive = d.drive === false; if (d.drive) d.value = r3(dimNow(d)); changed(); }},
    {label: d.evo ? `Evolution range ±${d.evo} (change)` : 'Let the evolution vary it…', do: () => {
      const L = LBL.find(l => l.dim === d.id); if (!L) return;
      inlineInput(L.el, d.evo || r3(d.axis === 'angle' ? 10 : CAD.lam() / 20), v => { snapshot(); d.evo = Math.max(0, v); d.drive = true; changed(); });
      hint(`Range ± ${d.axis === 'angle' ? '°' : 'mm'} the evolution may move this dimension (0 = fixed).`); }},
    '-', {label: 'Delete dimension', do: () => { snapshot(); DOC.dims = DOC.dims.filter(x => x !== d); changed(); }},
  ]);
}
function startPick(m) {
  hideHud(); pickMode = m; pickA = null; rebuildOverlays();
  hint({dimA: 'Dimension: pick the first point (dots, the light one is the origin). Esc cancels.', angA: 'Angle: click the first element. Esc cancels.',
        pivot: 'Pivot: click a point or a spot on an element. Esc cancels.'}[m]);
}
function onKey(ev) {
  if (!active || /INPUT|SELECT|TEXTAREA/.test(document.activeElement.tagName)) return;
  const k = ev.key.toLowerCase();
  if ((ev.ctrlKey || ev.metaKey) && k === 'z') { ev.shiftKey ? redo() : undo(); ev.preventDefault(); return; }
  if ((ev.ctrlKey || ev.metaKey) && k === 'y') { redo(); ev.preventDefault(); return; }
  if ((ev.ctrlKey || ev.metaKey) && k === 'd') { duplicate(); ev.preventDefault(); return; }
  if ((ev.ctrlKey || ev.metaKey) && k === 'a') { selectAll(); ev.preventDefault(); return; }
  if (ev.ctrlKey || ev.metaKey || ev.altKey) return;
  if (ev.key === 'Delete' || ev.key === 'Backspace') { delSel(); ev.preventDefault(); }
  else if (ev.key === 'Escape') { if (pickMode) { pickMode = null; rebuildOverlays(); hint(''); } else if (HUD) { hideHud(); drag = null; } else if (PIVOT) { PIVOT = null; rebuildOverlays(); placeGumball(); } else select([]); hideMenu(); }
  else if (k === 'w') gumball('move');
  else if (k === 'e') gumball('rotate');
  else if (k === 'q') gumball('both');
  else if (k === 'f') fit();
  else if (k === 'd') startPick('dimA');
  else if (k === 'a') startPick('angA');
  else if (k === 'j') joinSel();
  else if (k === 'p') { if (ev.shiftKey) { PIVOT = null; rebuildOverlays(); placeGumball(); hint('Pivot back on the element.'); } else if (SEL.length) startPick('pivot'); else hint('Select something first, then P sets its pivot.'); }
  else if (k === 'g') { $('cadsnap').checked = !$('cadsnap').checked; hint('Point snapping ' + ($('cadsnap').checked ? 'on' : 'off') + '.'); }
  else if (k === 'l') { $('cadspace').value = $('cadspace').value === 'world' ? 'local' : 'world'; snapSettings(); hint('Gumball axes: ' + ($('cadspace').value === 'world' ? 'world' : 'element')); }
  else if (k === '1') viewFrom('top'); else if (k === '2') viewFrom('front'); else if (k === '3') viewFrom('right'); else if (k === '0') viewFrom('iso');
}
function hint(t) { const h = $('cadhint'); if (h) h.textContent = t; }
export function selectAll() { select(DOC.elements.map(e => e.id)); hint(DOC.elements.length ? `All ${DOC.elements.length} elements selected. Del removes them (Ctrl+Z brings them back).` : 'The design is empty.'); }
export function clearAll() {   // throw the whole design away (undo keeps it)
  if (!DOC.elements.length) return hint('The design is already empty.');
  if (!confirm('Throw away the whole design (all elements, joints and dimensions)? Ctrl+Z brings it back.')) return;
  snapshot(); DOC.elements = []; DOC.joins = []; DOC.dims = []; PIVOT = null; select([]); changed(); hint('Design cleared. Ctrl+Z brings it back.');
}

// ------------------------------------------------------------------ editing operations

function delSel() {
  if (!SEL.length) return;
  snapshot();
  DOC.elements = DOC.elements.filter(e => !SEL.includes(e.id));
  const ok = r => r.el === 'origin' || DOC.elements.some(e => e.id === r.el);
  DOC.joins = DOC.joins.filter(j => ok(j.a) && ok(j.b)); DOC.dims = DOC.dims.filter(d => ok(d.a) && ok(d.b));
  SEL = []; PIVOT = null; changed();
}
function duplicate() {
  if (!SEL.length) return;
  snapshot();
  const out = [];
  for (const id of SEL) {
    const e = JSON.parse(JSON.stringify(elById(id))); e.id = 'e' + (seq++); e.name += ' copy'; e.feed = false; e.fixed = false;
    e.pos = e.pos.map(v => r3(v + sceneSize() / 10)); DOC.elements.push(e); out.push(e.id);
  }
  select(out); changed(true);
}
function split(f) {
  const e = elById(SEL[0]);
  if (!e || e.type !== 'rod') return hint('Split works on rods: select one rod.');
  snapshot();
  f = Math.min(0.95, Math.max(0.05, f ?? 0.5));
  const P = localPath(e); const L = []; let tot = 0;
  for (let i = 0; i < P.length - 1; i++) { L.push(Math.hypot(...P[i + 1].map((v, k) => v - P[i][k]))); tot += L[i]; }
  let s = f * tot, i = 0; while (i < L.length - 1 && s > L[i]) { s -= L[i]; i++; }
  const u = L[i] ? s / L[i] : 0, q = P[i].map((v, k) => v + u * (P[i + 1][k] - v));
  const A = [...P.slice(0, i + 1), q], B = [q, ...P.slice(i + 1)];
  const e2 = JSON.parse(JSON.stringify(e)); e2.id = 'e' + (seq++); e2.name = e.name + ' b'; e.name = e.name + ' a';
  const norm = pts => pts.map(p => p.map((v, k) => r3(v - q[k])));
  if (A.length === 2 && A[0].every(v => v === 0) && A[1][1] === 0 && A[1][2] === 0) e.params = {L: r3(A[1][0])}; else e.params = {pts: A.map(p => p.map(r3))};
  const qw = new THREE.Vector3(...q).applyMatrix4(mtx(e));
  e2.pos = [qw.x, qw.y, qw.z].map(r3); e2.params = {pts: norm(B)};
  if (B.length === 2 && Math.abs(B[1][1] - q[1]) < 1e-9 && Math.abs(B[1][2] - q[2]) < 1e-9 && B[1][0] > q[0]) e2.params = {L: r3(B[1][0] - q[0])};
  e2.feed = e.feed && e.feed_at > f; e.feed = e.feed && e.feed_at <= f;
  if (e2.feed) e2.feed_at = r3((e.feed_at - f) / (1 - f)); if (e.feed) e.feed_at = r3(e.feed_at / f);
  DOC.elements.splice(DOC.elements.indexOf(e) + 1, 0, e2);
  for (const j of DOC.joins) for (const r of [j.a, j.b]) if (r.el === e.id && r.pt === 'end') r.el = e2.id;
  DOC.joins.push({a: {el: e.id, pt: 'end'}, b: {el: e2.id, pt: 'start'}});
  select([e.id, e2.id]); changed();
}
function joinSel() {
  if (SEL.length !== 2) return hint('Select two elements (Shift+click) to join their nearest points.');
  snapshot();
  const [a, b] = SEL.map(elById); let best = null;
  for (const [na, pa] of Object.entries(localPoints(a))) for (const [nb, pb] of Object.entries(localPoints(b))) {
    const d = new THREE.Vector3(...pa).applyMatrix4(mtx(a)).distanceTo(new THREE.Vector3(...pb).applyMatrix4(mtx(b)));
    if ((na === 'center') !== (nb === 'center') && PATH.includes(a.type) && PATH.includes(b.type)) continue;
    if (!best || d < best.d) best = {d, na, nb};
  }
  DOC.joins.push({a: {el: a.id, pt: best.na}, b: {el: b.id, pt: best.nb}});
  solve(a.id); hint(`Joined ${a.name} (${best.na}) and ${b.name} (${best.nb}).`); changed(true);
}
function disconnect() {
  snapshot();
  const n = DOC.joins.length;
  DOC.joins = DOC.joins.filter(j => !(SEL.includes(j.a.el) || SEL.includes(j.b.el)));
  hint(`${n - DOC.joins.length} joint(s) removed.`); changed(true);
}
function setRodWorld(e, Pw) {   // write a world polyline back into a rod, keeping its pose where possible
  if (Pw.length === 2 && !e.params.pts) {
    const d = Pw[1].clone().sub(Pw[0]);
    e.pos = Pw[0].toArray().map(v => Math.round(v * 1e5) / 1e5); e.params = Object.assign({}, e.params, {L: r3(d.length())});
    return;
  }
  const inv = mtx(e).invert();
  e.params = {pts: Pw.map(p => p.clone().applyMatrix4(inv).toArray().map(v => Math.round(v * 1e4) / 1e4)), bend_r: e.params.bend_r || 0};
}
function filletPrompt(x, y) {
  if (SEL.length !== 2) return hint('Select the two joined rods (Shift+click), then Fillet.');
  const def = r3(Math.max(3 * CAD.wire(), CAD.lam() / 40));
  hideMenu();
  MENU = document.createElement('div'); MENU.className = 'cadmenu cadask';
  MENU.innerHTML = `<label>Fillet radius [mm] <input value="${def}"></label><span>Enter applies</span>`;
  document.body.appendChild(MENU); MENU.style.left = x + 'px'; MENU.style.top = y + 'px';
  const inp = MENU.querySelector('input'); inp.focus(); inp.select();
  inp.onkeydown = ev => { ev.stopPropagation(); if (ev.key === 'Enter') { const v = evalNum(inp.value); hideMenu(); if (v) fillet(v); } else if (ev.key === 'Escape') hideMenu(); };
}
export function fillet(R, ids) {
  const [a, b] = (ids || SEL).map(elById);
  if (!a || !b) return hint('Select the two joined rods.');
  const j = DOC.joins.find(j => (j.a.el === a.id && j.b.el === b.id) || (j.a.el === b.id && j.b.el === a.id));
  if (!j) return hint(`${a.name} and ${b.name} are not joined: join them first, then fillet the joint.`);
  const pa = j.a.el === a.id ? j.a.pt : j.b.pt, pb = j.a.el === a.id ? j.b.pt : j.a.pt;
  if (a.type !== 'rod' || b.type !== 'rod' || !['start', 'end'].includes(pa) || !['start', 'end'].includes(pb)) return hint('Fillet works on a joint between the ends of two rods.');
  let PA = worldPaths(a)[0], PB = worldPaths(b)[0];
  if (pa === 'start') PA = PA.slice().reverse();                 // A ends at the joint
  if (pb === 'end') PB = PB.slice().reverse();                   // B starts at the joint
  const J = PA[PA.length - 1], dA = J.clone().sub(PA[PA.length - 2]).normalize(), dB = PB[1].clone().sub(J).normalize();
  const th = Math.acos(Math.max(-1, Math.min(1, dA.dot(dB))));
  if (th < 1 * D2R) return hint('The rods are in line: there is no corner to round.');
  if (th > 175 * D2R) return hint('The rods fold back on each other: no fillet possible.');
  const t = R * Math.tan(th / 2), la = J.distanceTo(PA[PA.length - 2]), lb = J.distanceTo(PB[1]);
  if (t > la * 0.98 || t > lb * 0.98) return hint(`Radius too big: the fillet needs ${r3(t)} mm of straight wire on each side, there is ${r3(Math.min(la, lb))} mm.`);
  snapshot();
  const T1 = J.clone().addScaledVector(dA, -t), T2 = J.clone().addScaledVector(dB, t);
  const yv = dB.clone().addScaledVector(dA, -dA.dot(dB)).normalize(), zv = dA.clone().cross(yv);
  const F = defaults('fillet');
  F.name = 'Fillet ' + a.name.split(' ')[0]; F.mat = a.mat; F.d = a.d; F.params = {R: r3(R), angle: r3(th / D2R), bend_r: 0};
  setMtx(F, new THREE.Matrix4().makeBasis(dA, yv, zv).setPosition(T1));
  F.fil = {a: a.id, b: b.id};
  PA[PA.length - 1] = T1; PB[0] = T2;
  if (pa === 'start') PA.reverse(); if (pb === 'end') PB.reverse();
  setRodWorld(a, PA); setRodWorld(b, PB);
  DOC.joins = DOC.joins.filter(x => x !== j);
  DOC.elements.splice(DOC.elements.indexOf(a) + 1, 0, F);
  DOC.joins.push({a: {el: a.id, pt: pa}, b: {el: F.id, pt: 'start'}}, {a: {el: F.id, pt: 'end'}, b: {el: b.id, pt: pb}});
  hint(`Rounded the joint of ${a.name} and ${b.name} with R ${r3(R)} mm.`);
  select([F.id]); changed(true);
  return F;
}
function unfillet(F) {   // back to the sharp corner
  const ja = DOC.joins.find(j => (j.b.el === F.id && j.b.pt === 'start') || (j.a.el === F.id && j.a.pt === 'start'));
  const jb = DOC.joins.find(j => (j.b.el === F.id && j.b.pt === 'end') || (j.a.el === F.id && j.a.pt === 'end'));
  if (!ja || !jb) return null;
  const ra = ja.a.el === F.id ? ja.b : ja.a, rb = jb.a.el === F.id ? jb.b : jb.a;
  const a = elById(ra.el), b = elById(rb.el); if (!a || !b || a.type !== 'rod' || b.type !== 'rod') return null;
  const q = quatOf(F), x = new THREE.Vector3(1, 0, 0).applyQuaternion(q), t = F.params.R * Math.tan(F.params.angle * D2R / 2);
  const J = new THREE.Vector3(...F.pos).addScaledVector(x, t);
  let PA = worldPaths(a)[0], PB = worldPaths(b)[0];
  if (ra.pt === 'start') PA[0] = J.clone(); else PA[PA.length - 1] = J.clone();
  if (rb.pt === 'start') PB[0] = J.clone(); else PB[PB.length - 1] = J.clone();
  setRodWorld(a, PA); setRodWorld(b, PB);
  DOC.elements = DOC.elements.filter(e => e !== F);
  DOC.joins = DOC.joins.filter(j => j !== ja && j !== jb);
  DOC.joins.push({a: ra, b: rb});
  return [a.id, b.id];
}

// sketch -> custom element (called by the sketch editor)
async function exportSvg(e) {
  try {
    const q = await fetch('/api/pcb_svg', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({element: e})});
    const j = await q.json(); hint(q.ok ? 'SVG saved: ' + j.path : 'SVG export failed: ' + j.error);
  } catch (ex) { hint('SVG export failed: ' + ex); }
}
export function placeCustom(res, target) {
  snapshot();
  let e = target ? elById(target.id) : null;
  const kind = res.kind === 'pcb' ? 'pcb' : 'custom';
  if (e && e.type !== kind) { e.type = kind; e.evolve = Object.assign({}, e.evolve, {mode: 'off'}); }
  if (!e) { e = defaults(kind); e.name = res.name || e.name; e.pos = [0, 0, 0]; DOC.elements.push(e); }
  if (kind === 'pcb') { e.params = Object.assign({}, res.params, {sketch: res.sketch}); e.feed = !!res.feed; }
  else {
  e.params = {paths: res.paths, nodes: res.nodes, pts: res.pts, sketch: res.sketch};
  if (res.feed) { e.feed = true; e.feed_path = res.feed.path; e.feed_at = res.feed.at; } else e.feed = false;
  }
  if (res.name) e.name = res.name;
  const ok = r => r.el !== e.id || r.pt in localPoints(e);
  DOC.joins = DOC.joins.filter(j => ok(j.a) && ok(j.b)); DOC.dims = DOC.dims.filter(d => d.axis === 'angle' || (ok(d.a) && ok(d.b)));
  select([e.id]); changed();
  return e;
}

// ------------------------------------------------------------------ browser (element tree) and inspector, floating in the view

function numIn(id, v, step = 'any') { return `<input id="${id}" type="text" inputmode="decimal" value="${v}" data-step="${step}">`; }
export function panel() {
  if (!DOC || !$('cadels')) return;
  $('cadels').innerHTML = DOC.elements.length ? DOC.elements.map(e => {
    const j = DOC.joins.filter(x => x.a.el === e.id || x.b.el === e.id).length, ev = e.evolve || {};
    const tags = [e.feed ? '<i class="tf">feed</i>' : '', j ? `<i>${j}⌁</i>` : '', e.fixed ? '<i>fix</i>' : '', ev.mode && ev.mode !== 'off' ? `<i class="te">${ev.mode}</i>` : ''].join('');
    return `<div class="cadel${SEL.includes(e.id) ? ' on' : ''}" data-id="${e.id}"><span class="sw" style="background:${matOf(e.mat).color}"></span><b>${esc(e.name)}</b><span class="tags">${tags}</span></div>`;
  }).join('') : '<div class="hint">Empty. Create an element in the toolbar, or start from a template.</div>';
  for (const d of document.querySelectorAll('.cadel')) {
    d.onclick = ev => { const id = d.dataset.id; select(ev.shiftKey || ev.ctrlKey ? (SEL.includes(id) ? SEL.filter(x => x !== id) : [...SEL, id]) : [id]); };
    d.ondblclick = () => { const e = elById(d.dataset.id); if (e && (e.type === 'custom' || e.type === 'pcb') && CAD.editSketch) CAD.editSketch(e); };
  }
  const nm = r => r.el === 'origin' ? 'origin' : (elById(r.el) || {name: '?'}).name + (r.pt ? '.' + r.pt : '');
  $('caddims').innerHTML = DOC.dims.length ? DOC.dims.map(d => `<div class="caddim" data-d="${d.id}"><span>${d.axis === 'angle' ? '∠ ' : d.axis === 'd' ? '↔ ' : d.axis + ' '}${esc(nm(d.a))} → ${esc(nm(d.b))}</span><b>${r3(d.value)}${d.axis === 'angle' ? '°' : ''}</b>${d.evo ? `<i>±${d.evo}</i>` : ''}${d.drive === false ? '<i>ref</i>' : ''}</div>`).join('')
    : '<div class="hint">D = distance, A = angle. Values are typed right in the view.</div>';
  if (CAD.conflict) $('caddims').insertAdjacentHTML('afterbegin', `<div class="warn">Constraints cannot all be met (off by ${CAD.conflict.toFixed(2)} mm). Unfix an element or remove a dimension.</div>`);
  for (const r of document.querySelectorAll('.caddim')) r.oncontextmenu = r.onclick = ev => { ev.preventDefault(); const d = DOC.dims.find(x => x.id === r.dataset.d); if (d) dimMenu(d, ev.clientX, ev.clientY); };
  inspector();
  if ($('cadground')) { $('cadground').value = DOC.ground; $('cadz0').value = DOC.z0; $('cadswr').value = DOC.swr; $('cadengine').value = DOC.engine || 'nec'; if ($('cadmesh')) $('cadmesh').value = DOC.ems_mesh || 'normal'; }
  const u = $('cadundo'), rd = $('cadredo'); if (u) { u.disabled = !UNDO.length; rd.disabled = !REDO.length; }
}
function inspector() {
  const box = $('cadprops'); if (!box) return;
  const e = SEL.length === 1 ? elById(SEL[0]) : null;
  $('cadinsp').classList.toggle('nosel', !SEL.length);
  if (!SEL.length) { box.innerHTML = '<div class="hint">Click an element. Shift+click adds more. Right-click for the commands, double-click a value in the view to type it.</div>'; return; }
  const xform = `<div class="sub">Transform${PIVOT ? ' <i>pivot set</i>' : ''}</div>
    <div class="row3">${['x', 'y', 'z'].map(a => `<label>Δ${a}${numIn('cx_m' + a, 0)}</label>`).join('')}</div>
    <div class="row3">${['x', 'y', 'z'].map(a => `<label>∠${a}${numIn('cx_r' + a, 0)}</label>`).join('')}</div>
    <div class="ops"><button id="cx_apply">Apply</button><button id="cx_pivot">Pivot ▸ pick</button><button id="cx_pivr"${PIVOT ? '' : ' disabled'}>Pivot reset</button></div>`;
  if (!e) {
    box.innerHTML = `<div class="hint">${SEL.length} elements selected: the gumball moves them together.</div>${xform}`;
    bindXform(); return;
  }
  const isPath = PATH.includes(e.type), m = matOf(e.mat), ev = Object.assign({mode: 'off', dev: 0, size: 15, rot: 0, pts: 4, pix: [8, 8]}, e.evolve || {});
  if (!e.evolve || !e.evolve.mode) { if ((ev.dev || 0) > 0) ev.mode = e.type === 'rod' ? 'free' : 'shape'; if (ev.bends) ev.pts = ev.bends + 2; }
  let h = `<div class="r2"><input id="cp_name" value="${esc(e.name)}"><span class="badge">${LABEL[e.type]}</span></div>
    <div class="row"><label>Material<select id="cp_mat">${CAD.mats.map(x => `<option${x.name === e.mat ? ' selected' : ''}>${esc(x.name)}</option>`).join('')}</select></label>
    ${isPath || e.type === 'plate' ? `<label>${isPath ? 'Wire' : 'Grid wire'} Ø [mm]${numIn('cp_d', e.d)}</label>` : '<span></span>'}</div>`;
  if (e.type === 'rod' && e.params.pts) h += `<div class="hint">Bent rod, ${e.params.pts.length} points. <button class="link" id="cp_straight">Make straight</button></div>`;
  if (e.type === 'custom') h += `<div class="hint">${e.params.sketch_stale ? 'Shape changed by the evolution. ' : ''}<button class="link" id="cp_sketch">Edit drawing</button> (or double-click it)</div>`;
  if (e.type === 'pcb') h += `<div class="hint">${e.params.sketch_stale ? 'Shape changed by the evolution. ' : ''}<button class="link" id="cp_sketch">Edit drawing</button> · <button class="link" id="cp_svg">Export SVG</button>${e.params.feed ? ' · <button class="link" id="cp_noprobe">Remove probe feed</button>' : ''}</div>
    <div class="row"><label>Substrate<select id="cp_sub"><option value="">air (none)</option>${CAD.mats.filter(x => x.kind === 'dielectric').map(x => `<option${x.name === e.params.sub ? ' selected' : ''}>${esc(x.name)}</option>`).join('')}</select></label><span></span></div>
    ${(e.params.sub && !(CAD.mats.find(x => x.name === e.params.sub) || {}).eps) ? '' : ''}`;
  if (e.type === 'fillet') h += `<div class="hint"><button class="link" id="cp_unfil">Remove fillet (sharp corner)</button></div>`;
  const ps = PARAMS[e.type].filter(([k]) => !(e.type === 'rod' && e.params.pts && k === 'L'));
  if (ps.length) h += '<div class="row">' + ps.map(([k, lab, u, opts]) => `<label>${lab}${u ? ' [' + u + ']' : ''}${opts ? `<select id="cp_p_${k}">${opts.map(o => `<option${e.params[k] === o ? ' selected' : ''}>${o}</option>`).join('')}</select>` : numIn('cp_p_' + k, e.params[k] ?? 0)}</label>`).join('') + '</div>';
  h += `<div class="sub">Pose</div><div class="row3">${[0, 1, 2].map(k => `<label>${'xyz'[k]} [mm]${numIn('cp_pos' + k, e.pos[k])}</label>`).join('')}</div>
    <div class="row3">${[0, 1, 2].map(k => `<label>r${'xyz'[k]} [°]${numIn('cp_rot' + k, e.rot[k])}</label>`).join('')}</div>
    <label class="check"><input type="checkbox" id="cp_fixed"${e.fixed ? ' checked' : ''}> Fixed (the solver never moves it)</label>`;
  if (isPath && m.kind === 'conductor') {
    const np = localPaths(e).length;
    h += `<label class="check"><input type="checkbox" id="cp_feed"${e.feed ? ' checked' : ''}> Feed point on this wire</label>
      ${e.feed ? `<div class="row"><label>Along the wire [%]${numIn('cp_feed_at', r3(e.feed_at * 100), 1)}</label>${np > 1 ? `<label>Wire run (1–${np})${numIn('cp_feed_path', (e.feed_path || 0) + 1, 1)}</label>` : '<span></span>'}</div>` : ''}`;
  }
  if (m.kind === 'dielectric') h += `<div class="hint">Dielectric: NEC-2 ignores it, the openEMS engine models it.</div>`;
  const modes = [['off', 'off'], ['shape', 'shape: own sizes, position, rotation']];
  if (isPath) modes.push(['free', 'free form: control points move']);
  if (e.type === 'pcb') modes.push(['free', 'free form: outline corners, vias and feed move']);
  if (e.type === 'plate' || e.type === 'pcb') modes.push(['pixels', 'pixels: cells on / off']);
  h += `<div class="sub">Evolution freedom</div>
    <label>Mode<select id="cp_emode">${modes.map(([v, t]) => `<option value="${v}"${ev.mode === v ? ' selected' : ''}>${t}</option>`).join('')}</select></label>`;
  if (ev.mode === 'shape') h += `<div class="row3"><label>Sizes ± [%]${numIn('cp_esize', ev.size || 0)}</label><label>Move ± [mm]${numIn('cp_edev', ev.dev || 0)}</label><label>Turn ± [°]${numIn('cp_erot', ev.rot || 0)}</label></div>
      <label class="check"><input type="checkbox" id="cp_eplanar"${ev.planar ? ' checked' : ''}> Stay in its own plane</label>
      <div class="hint">Changes ${SHAPE_KEYS[e.type].join(', ') || 'its pose'} within the range; joints and driving dimensions are kept by the solver.</div>`;
  if (ev.mode === 'free' && e.type === 'pcb') h += `<div class="row"><label>Move ± [mm]${numIn('cp_edev', ev.dev || 0)}</label><label>Stations on round outlines${numIn('cp_epts', ev.pts || 4, 1)}</label></div>
      <div class="hint">Every drawn corner of the top copper moves within the range (the edges between them follow), and so do the vias and the probe feed. The ground stays as drawn.</div>`;
  if (ev.mode === 'free' && e.type !== 'pcb') h += `<div class="row"><label>Control points${numIn('cp_epts', ev.pts || 4, 1)}</label><label>Move ± [mm]${numIn('cp_edev', ev.dev || 0)}</label></div>
      <label class="check"><input type="checkbox" id="cp_eplanar"${ev.planar ? ' checked' : ''}> Stay in its own plane</label>
      <div class="hint">The wire becomes a chain of ${e.type === 'custom' ? 'its drawn corner points' : (ev.pts || 4) + ' control points'} that the evolution moves; the wire between them follows, so a ${LABEL[e.type].toLowerCase()} keeps its character. Joined points stay joined.</div>`;
  if (ev.mode === 'pixels') h += `<div class="row"><label>Pixels x${numIn('cp_epx', ev.pix[0], 1)}</label><label>Pixels y${numIn('cp_epy', ev.pix[1], 1)}</label></div>
      <div class="hint">The evolution switches each pixel on or off (a "pixel antenna"). Pixels under a probe stay metal.${e.type === 'pcb' ? ' The grid covers the top copper; only cells inside the drawn area can be metal.' : ''}</div>`;
  h += xform;
  box.innerHTML = h;
  const set = f => ev2 => { snapshot(); f(ev2.target); solve(e.id); changed(true); };
  const num = t => evalNum(t.value);
  $('cp_name').onchange = set(t => { e.name = t.value.trim() || e.name; });
  $('cp_mat').onchange = set(t => { e.mat = t.value; if (matOf(e.mat).kind !== 'conductor') e.feed = false; });
  if ($('cp_d')) $('cp_d').onchange = set(t => { e.d = Math.max(0.01, num(t) || e.d); });
  for (const [k, , , opts] of PARAMS[e.type]) { const el = $('cp_p_' + k); if (el) el.onchange = set(t => { e.params[k] = opts ? t.value : (num(t) ?? e.params[k]); }); }
  for (const k of [0, 1, 2]) { $('cp_pos' + k).onchange = set(t => { e.pos[k] = num(t) ?? e.pos[k]; }); $('cp_rot' + k).onchange = set(t => { e.rot[k] = num(t) ?? e.rot[k]; }); }
  $('cp_fixed').onchange = set(t => { e.fixed = t.checked; });
  if ($('cp_feed')) $('cp_feed').onchange = set(t => { e.feed = t.checked; });
  if ($('cp_feed_at')) $('cp_feed_at').onchange = set(t => { e.feed_at = Math.min(1, Math.max(0, (num(t) || 0) / 100)); });
  if ($('cp_feed_path')) $('cp_feed_path').onchange = set(t => { e.feed_path = Math.max(0, Math.min(localPaths(e).length - 1, Math.round(num(t) || 1) - 1)); });
  if ($('cp_straight')) $('cp_straight').onclick = set(() => { const P = e.params.pts; e.params = {L: r3(Math.hypot(...P[P.length - 1].map((v, k) => v - P[0][k])))}; });
  if ($('cp_sketch')) $('cp_sketch').onclick = () => CAD.editSketch && CAD.editSketch(e);
  if ($('cp_svg')) $('cp_svg').onclick = () => exportSvg(e);
  if ($('cp_noprobe')) $('cp_noprobe').onclick = set(() => { e.params.feed = null; e.feed = false; const sk = e.params.sketch; if (sk && sk.probe) { sk.ents = (sk.ents || []).filter(q => !(q.probe && q.a === sk.probe)); sk.probe = null; } });
  if ($('cp_sub')) $('cp_sub').onchange = set(t => { e.params.sub = t.value; });
  if ($('cp_unfil')) $('cp_unfil').onclick = () => { snapshot(); const r = unfillet(e); if (r) { select(r); changed(true); } else hint('This fillet is not between two rods any more.'); };
  const evo = f => set(t => { e.evolve = Object.assign({}, ev, e.evolve || {}); f(t, e.evolve); delete e.evolve.bends; });
  $('cp_emode').onchange = evo((t, o) => { o.mode = t.value; if (o.mode !== 'off' && !(o.dev > 0) && o.mode !== 'pixels') o.dev = r3(CAD.lam() / 40); });
  if ($('cp_esize')) $('cp_esize').onchange = evo((t, o) => { o.size = Math.max(0, Math.min(90, num(t) || 0)); });
  if ($('cp_edev')) $('cp_edev').onchange = evo((t, o) => { o.dev = Math.max(0, num(t) || 0); });
  if ($('cp_erot')) $('cp_erot').onchange = evo((t, o) => { o.rot = Math.max(0, Math.min(180, num(t) || 0)); });
  if ($('cp_epts')) $('cp_epts').onchange = evo((t, o) => { o.pts = Math.max(2, Math.min(24, Math.round(num(t) || 4))); });
  if ($('cp_epx')) $('cp_epx').onchange = evo((t, o) => { o.pix = [Math.max(1, Math.min(24, Math.round(num(t) || 8))), o.pix[1]]; });
  if ($('cp_epy')) $('cp_epy').onchange = evo((t, o) => { o.pix = [o.pix[0], Math.max(1, Math.min(24, Math.round(num(t) || 8)))]; });
  if ($('cp_eplanar')) $('cp_eplanar').onchange = evo((t, o) => { o.planar = t.checked; });
  bindXform();
}
function bindXform() {
  $('cx_pivot').onclick = () => startPick('pivot');
  $('cx_pivr').onclick = () => { PIVOT = null; rebuildOverlays(); placeGumball(); panel(); };
  $('cx_apply').onclick = () => {
    const g = k => evalNum($(k).value) || 0;
    const mv = new THREE.Vector3(g('cx_mx'), g('cx_my'), g('cx_mz')), rv = [g('cx_rx'), g('cx_ry'), g('cx_rz')];
    if (!mv.length() && !rv.some(Boolean)) return;
    snapshot();
    const c = PIVOT ? PIVOT.clone() : selCenter();
    let D = new THREE.Matrix4().makeTranslation(c.x, c.y, c.z);
    for (const [k, a] of [[0, 'x'], [1, 'y'], [2, 'z']]) if (rv[k]) D.multiply(new THREE.Matrix4().makeRotationAxis(new THREE.Vector3(k === 0 ? 1 : 0, k === 1 ? 1 : 0, k === 2 ? 1 : 0), rv[k] * D2R));
    D.multiply(new THREE.Matrix4().makeTranslation(-c.x, -c.y, -c.z));
    D.premultiply(new THREE.Matrix4().makeTranslation(mv.x, mv.y, mv.z));
    for (const id of SEL) { const e = elById(id); setMtx(e, D.clone().multiply(mtx(e))); }
    if (PIVOT) PIVOT.applyMatrix4(D);
    solve(new Set(SEL)); changed(true);
    hint('Transformed' + (PIVOT ? ' about the pivot.' : ' about the selection.'));
  };
}
function propsLive() {
  const e = SEL.length === 1 ? elById(SEL[0]) : null; if (!e) return;
  for (const k of [0, 1, 2]) { const p = $('cp_pos' + k), r = $('cp_rot' + k); if (p) p.value = e.pos[k]; if (r) r.value = e.rot[k]; }
}

function bindPanel() {
  for (const b of document.querySelectorAll('[data-add]')) b.onclick = () => addElement(b.dataset.add);
  for (const b of document.querySelectorAll('[data-gb]')) b.onclick = () => gumball(b.dataset.gb);
  for (const b of document.querySelectorAll('[data-view]')) b.onclick = () => viewFrom(b.dataset.view);
  const on = (id, f) => { const x = $(id); if (x) x.onclick = f; };
  on('caddel', delSel); on('caddup', duplicate); on('cadsplit', () => split(0.5)); on('cadjoin', joinSel); on('caddisc', disconnect);
  on('caddim', () => startPick('dimA')); on('cadangle', () => startPick('angA')); on('cadpivot', () => SEL.length ? startPick('pivot') : hint('Select something first.'));
  on('cadfillet', ev => filletPrompt(ev.clientX, ev.clientY + 20)); on('cadfit', fit); on('cadundo', undo); on('cadredo', redo);
  on('cadsketch', () => CAD.editSketch && CAD.editSketch(null));
  on('cadall', selectAll); on('cadclear', clearAll);
  $('cadtpl').onchange = () => { const v = $('cadtpl').value; $('cadtpl').value = ''; if (v && (!DOC.elements.length || confirm('Replace the current design with this template?'))) template(v); };
  for (const id of ['cadgrid', 'cadang', 'cadspace']) $(id).onchange = snapSettings;
  snapSettings();
  $('cadground').onchange = () => { DOC.ground = $('cadground').value; changed(true); };
  $('cadengine').onchange = () => { DOC.engine = $('cadengine').value; changed(true); };
  $('cadz0').onchange = () => { DOC.z0 = Number($('cadz0').value) || 50; changed(true); };
  $('cadswr').onchange = () => { DOC.swr = Number($('cadswr').value) || 2; changed(true); };
  for (const b of document.querySelectorAll('.fpanel .fhead')) b.onclick = () => b.parentElement.classList.toggle('min');
}

export function _screen(v) {   // test helper: screen position of a world point
  const r = renderer.domElement.getBoundingClientRect(), p = new THREE.Vector3(...v).project(camera);
  return [r.left + (p.x + 1) / 2 * r.width, r.top + (1 - p.y) / 2 * r.height];
}
export const _test = {selCenter: () => selCenter(), pivot: () => PIVOT, setPivot: v => { PIVOT = v ? new THREE.Vector3(...v) : null; placeGumball(); }, labels: () => LBL.map(l => l.el.textContent), split, joinSel};
