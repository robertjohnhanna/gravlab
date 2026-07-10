// sources.js — user-facing mass sources built on the verified physics core.
// Each source carries a world position + orientation and type-specific params.
// buildSource() precomputes caches (rotation matrix, ring stacks); sourceField()
// returns { g, phi } — acceleration [m/s²] and potential [J/kg] at a world point.
import * as P from './physics.js';

// Named reference masses [kg] — pick one to set a source's mass, then fine-tune.
// (The gravitational analogue of a magnet-grade table: it seeds one number.)
export const BODIES = {
  'Mountain (1e15)': 1e15,
  'Asteroid (1e18)': 1e18,
  'Ceres (9.4e20)':  9.4e20,
  'Moon (7.3e22)':   7.342e22,
  'Mars (6.4e23)':   6.417e23,
  'Earth (6.0e24)':  5.972e24,
  'Neptune (1.0e26)':1.024e26,
  'Jupiter (1.9e27)':1.898e27,
  'Sun (2.0e30)':    1.989e30,
};
// Real mean diameters [km] for the same named bodies, so choosing a preset sets
// a physically-consistent size as well as a mass (for sphere / shell bodies).
export const BODY_DIA = {
  'Mountain (1e15)': 10,
  'Asteroid (1e18)': 100,
  'Ceres (9.4e20)':  940,
  'Moon (7.3e22)':   3474,
  'Mars (6.4e23)':   6779,
  'Earth (6.0e24)':  12742,
  'Neptune (1.0e26)':49244,
  'Jupiter (1.9e27)':139820,
  'Sun (2.0e30)':    1391000,
};

let _id = 1;
export const nextId = () => _id++;

// Default parameter sets. Lengths are in kilometres in the UI, mass in kg;
// lengths convert to metres on build. Everything here is in UI units.
export function defaultSource(type) {
  const base = {
    id: nextId(), type, name: '', visible: true,
    pos: [0, 0, 0],            // km
    rot: [0, 0, 0],            // yaw,pitch,roll in degrees
    vel: [0, 0, 0],            // km/s — initial velocity for the N-body simulation
    mass: 5.972e24,            // kg
    color: pickColor(),
  };
  switch (type) {
    case 'point':    return { ...base, name: 'Point mass', mass: 5e23 };
    case 'sphere':   return { ...base, name: 'Planet',     dia: 12000 };
    case 'shell':    return { ...base, name: 'Shell',      dia: 16000, mass: 3e24 };
    case 'ring':     return { ...base, name: 'Ring',       dia: 30000, mass: 2e24 };
    case 'disc':     return { ...base, name: 'Disc',       dia: 28000, thick: 3000, mass: 3e24 };
    case 'cylinder': return { ...base, name: 'Cylinder',   dia: 14000, len: 26000, mass: 3e24 };
    case 'rod':      return { ...base, name: 'Rod',        len: 40000, mass: 1.5e24 };
    case 'box':      return { ...base, name: 'Slab',       size: [14000, 14000, 28000], mass: 3e24 };
    default: throw new Error('unknown source type ' + type);
  }
}

const PALETTE = ['#e9b44c', '#4aa3ff', '#37b36b', '#e0682f', '#9b7bf0', '#e05f9e', '#20c4c4', '#c7ccd6'];
let _cix = 0;
function pickColor() { return PALETTE[(_cix++) % PALETTE.length]; }

const km = (v) => v * 1000;      // km -> m
const deg = (v) => v * Math.PI / 180;

// Ring decomposition of a uniform disc of radius R (m) and mass (kg), placed at
// axial height z (m): nr coaxial rings whose masses reproduce the σ = M/πR²
// areal density.  Building block for both discs and cylinders.
function discRings(R, mass, nr, z) {
  const dr = R / nr, out = [];
  let wsum = 0;
  for (let i = 0; i < nr; i++) wsum += (i + 0.5);   // Σ r_i ∝ Σ(i+½)
  for (let i = 0; i < nr; i++) {
    const a = dr * (i + 0.5);
    out.push({ a, z, m: mass * (i + 0.5) / wsum });
  }
  return out;
}

// Build world-space caches. Call whenever pos/rot/params change.
export function buildSource(s) {
  const R = P.eulerToMatrix(deg(s.rot[0]), deg(s.rot[1]), deg(s.rot[2]));
  s._R = R;
  s._origin = [km(s.pos[0]), km(s.pos[1]), km(s.pos[2])];
  s._rings = null;     // ring stacks (ring / disc / cylinder), LOCAL frame

  if (s.type === 'ring') {
    s._rings = [{ a: km(s.dia) / 2, z: 0, m: s.mass }];
  } else if (s.type === 'disc') {
    // A disc with a settable thickness: rings over the radius, stacked over the
    // thickness (a single sheet when the thickness is ~0).
    const Rr = km(s.dia) / 2, L = km(s.thick || 0), nr = 22;
    s._rings = [];
    if (L > 1) {
      const nz = 4;
      for (let k = 0; k < nz; k++) {
        const z = -L / 2 + L * (k + 0.5) / nz;
        for (const r of discRings(Rr, s.mass / nz, nr, z)) s._rings.push(r);
      }
    } else {
      for (const r of discRings(Rr, s.mass, nr, 0)) s._rings.push(r);
    }
  } else if (s.type === 'cylinder') {
    const Rr = km(s.dia) / 2, L = km(s.len), nz = 6, nr = 12;
    s._rings = [];
    for (let k = 0; k < nz; k++) {
      const z = -L / 2 + L * (k + 0.5) / nz;
      for (const r of discRings(Rr, s.mass / nz, nr, z)) s._rings.push(r);
    }
  }
  return s;
}

// Field of a single source at world point Q -> { g:[3], phi:number }.
// Pass wantPhi = false to skip the potential — the field-only callers (heatmap
// gradient sampling, the force/torque volume integral) are the hot paths, and
// for ring stacks Φ costs a second elliptic-integral pass per ring.
export function sourceField(s, Q, wantPhi = true) {
  if (!s.visible) return { g: [0, 0, 0], phi: 0 };
  let g = [0, 0, 0], phi = 0;

  if (s.type === 'point') {
    g = P.pointField(s.mass, s._origin, Q);
    if (wantPhi) phi = P.pointPot(s.mass, s._origin, Q);
  } else if (s.type === 'sphere') {
    const Rr = km(s.dia) / 2;
    g = P.sphereField(s.mass, Rr, s._origin, Q);
    if (wantPhi) phi = P.spherePot(s.mass, Rr, s._origin, Q);
  } else if (s.type === 'shell') {
    const Rr = km(s.dia) / 2;
    g = P.shellField(s.mass, Rr, s._origin, Q);
    if (wantPhi) phi = P.shellPot(s.mass, Rr, s._origin, Q);
  } else if (s.type === 'rod') {
    const qLoc = P.matTVec(s._R, P.vsub(Q, s._origin));
    const L = km(s.len), lam = s.mass / L;
    g = P.matVec(s._R, P.rodField(lam, L, qLoc[0], qLoc[1], qLoc[2]));
    if (wantPhi) phi = P.rodPot(lam, L, qLoc[0], qLoc[1], qLoc[2]);
  } else if (s.type === 'box') {
    const qLoc = P.matTVec(s._R, P.vsub(Q, s._origin));
    const half = [km(s.size[0]) / 2, km(s.size[1]) / 2, km(s.size[2]) / 2];
    const rhod = s.mass / (8 * half[0] * half[1] * half[2]);
    g = P.matVec(s._R, P.cuboidField(qLoc, half, rhod));
    if (wantPhi) phi = P.cuboidPot(qLoc, half, rhod);
  } else if (s._rings) {
    const q = P.matTVec(s._R, P.vsub(Q, s._origin));
    let gx = 0, gy = 0, gz = 0;
    for (const r of s._rings) {
      const gl = P.ringField(r.a, r.m, q[0], q[1], q[2] - r.z);
      gx += gl[0]; gy += gl[1]; gz += gl[2];
      if (wantPhi) phi += P.ringPot(r.a, r.m, q[0], q[1], q[2] - r.z);
    }
    g = P.matVec(s._R, [gx, gy, gz]);
  }
  return { g, phi };
}

// Total mass [kg] of a source.
export function massOf(s) { return s.mass; }

// ---------------------------------------------------------------------------
// Scene: a collection of sources.  Provides total field, potential and the
// test-particle acceleration function.
// ---------------------------------------------------------------------------
export class Scene {
  constructor() { this.sources = []; }
  add(s) { this.sources.push(buildSource(s)); return s; }
  remove(id) { this.sources = this.sources.filter((s) => s.id !== id); }
  get(id) { return this.sources.find((s) => s.id === id); }
  rebuild() { this.sources.forEach(buildSource); }

  // Total gravitational field g [m/s²] at world point Q.
  g(Q) {
    let a = [0, 0, 0];
    for (const s of this.sources) a = P.vadd(a, sourceField(s, Q, false).g);
    return a;
  }
  // Total gravitational potential Φ [J/kg] at world point Q.
  potential(Q) {
    let p = 0;
    for (const s of this.sources) p += sourceField(s, Q).phi;
    return p;
  }
  // Both g and Φ in one pass over the sources — use this where both are needed
  // (grid sampling, probe): separate g() + potential() calls double the work.
  sample(Q) {
    let a = [0, 0, 0], p = 0;
    for (const s of this.sources) { const f = sourceField(s, Q); a = P.vadd(a, f.g); p += f.phi; }
    return { g: a, phi: p };
  }
  // Exact net force [N] and torque [N·m] on a source from all *other* sources.
  forceTorque(target) { return forceOn(this, target); }
}

// True if world point q lies inside the solid body of source s.
// (Point masses, rods and rings have no interior volume and never contain.)
export function bodyContains(s, q) {
  if (s.type === 'sphere' || s.type === 'shell')
    return P.vlen(P.vsub(q, s._origin)) < km(s.dia) / 2;
  if (s.type === 'box') {
    const l = P.matTVec(s._R, P.vsub(q, s._origin));
    return Math.abs(l[0]) < km(s.size[0]) / 2 && Math.abs(l[1]) < km(s.size[1]) / 2 && Math.abs(l[2]) < km(s.size[2]) / 2;
  }
  if (s.type === 'cylinder' || s.type === 'disc') {
    const l = P.matTVec(s._R, P.vsub(q, s._origin));
    const hl = s.type === 'cylinder' ? km(s.len) / 2 : km(s.thick || 0) / 2;
    return Math.hypot(l[0], l[1]) < km(s.dia) / 2 && Math.abs(l[2]) < hl;
  }
  return false;
}

// Oriented bounding box for a body (spheres/cylinders use their enclosing box —
// conservative, the safe direction for refusing an overlapping-force query).
function obbOf(s) {
  let he;
  if (s.type === 'box') he = [km(s.size[0]) / 2, km(s.size[1]) / 2, km(s.size[2]) / 2];
  else if (s.type === 'sphere' || s.type === 'shell') { const r = km(s.dia) / 2; he = [r, r, r]; }
  else if (s.type === 'cylinder') { const r = km(s.dia) / 2; he = [r, r, km(s.len) / 2]; }
  else if (s.type === 'disc') { const r = km(s.dia) / 2; he = [r, r, km(s.thick || 0) / 2]; }
  else if (s.type === 'ring') { const r = km(s.dia) / 2; he = [r, r, r * 0.02]; }
  else if (s.type === 'rod') { he = [0, 0, km(s.len) / 2]; }
  else return null;                                  // point mass — no body
  const R = s._R;
  const ax = [[R[0][0], R[1][0], R[2][0]], [R[0][1], R[1][1], R[2][1]], [R[0][2], R[1][2], R[2][2]]];
  return { c: s._origin, he, ax };
}
// Exact oriented-box overlap via the Separating Axis Theorem.
function bodiesOverlap(sa, sb) {
  const a = obbOf(sa), b = obbOf(sb);
  if (!a || !b) return false;
  const T = P.vsub(b.c, a.c);
  const axes = [a.ax[0], a.ax[1], a.ax[2], b.ax[0], b.ax[1], b.ax[2]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    const cr = P.vcross(a.ax[i], b.ax[j]);
    if (P.vlen(cr) > 1e-9) axes.push(P.vnorm(cr));
  }
  for (const L of axes) {
    let ra = 0, rb = 0;
    for (let i = 0; i < 3; i++) { ra += a.he[i] * Math.abs(P.vdot(a.ax[i], L)); rb += b.he[i] * Math.abs(P.vdot(b.ax[i], L)); }
    if (Math.abs(P.vdot(T, L)) > ra + rb + 1e-6) return false;   // a separating axis exists
  }
  return true;
}

// Exact net force [N] and torque [N·m] on `target` from all OTHER sources.
// The force density is ρ·g_ext integrated over the target body (g_ext = the
// field of the other sources only). A spherically-symmetric body (point / sphere
// / shell) responds exactly as a point mass at its centre — Newton's theorem —
// so those need no integration and feel zero self-torque. Extended bodies (rod,
// ring, disc, cylinder, box) are integrated over their mass, which also yields
// the real gravity-gradient torque. Valid whenever the bodies don't overlap.
export function forceOn(scene, target) {
  const others = scene.sources.filter((s) => s !== target && s.visible);
  if (!others.length) return { F: [0, 0, 0], tau: [0, 0, 0], valid: true, hasExternal: false };
  for (const o of others) if (bodiesOverlap(target, o)) return { F: [0, 0, 0], tau: [0, 0, 0], valid: false, hasExternal: true };
  const c = target._origin, R = target._R;
  const gExt = (q) => { let a = [0, 0, 0]; for (const s of others) a = P.vadd(a, sourceField(s, q, false).g); return a; };
  let F = [0, 0, 0], tau = [0, 0, 0], valid = true, Fabs = 0, tauAbs = 0;
  const local = (loc) => P.vadd(c, P.matVec(R, loc));          // local (m) -> world
  const inOther = (q) => { for (const s of others) if (bodyContains(s, q)) return true; return false; };
  const add = (r, dm) => {
    const dF = P.vscale(gExt(r), dm);
    F = P.vadd(F, dF); Fabs += P.vlen(dF);
    const t = P.vcross(P.vsub(r, c), dF); tau = P.vadd(tau, t); tauAbs += P.vlen(t);
  };

  if (target.type === 'point' || target.type === 'sphere' || target.type === 'shell') {
    // Newton's theorem: net force is exactly M·g_ext(centre); no self-torque.
    F = P.vscale(gExt(c), target.mass); Fabs = P.vlen(F);
    if (inOther(c)) valid = false;
  } else if (target.type === 'rod') {
    const L = km(target.len), n = 60, dm = target.mass / n;
    for (let i = 0; i < n; i++) {
      const r = local([0, 0, -L / 2 + L * (i + 0.5) / n]);
      if (inOther(r)) valid = false;
      add(r, dm);
    }
  } else if (target.type === 'ring') {
    const a = km(target.dia) / 2, n = 72, dm = target.mass / n;
    for (let i = 0; i < n; i++) {
      const th = 2 * Math.PI * (i + 0.5) / n;
      const r = local([a * Math.cos(th), a * Math.sin(th), 0]);
      if (inOther(r)) valid = false;
      add(r, dm);
    }
  } else if (target.type === 'disc') {
    const Rr = km(target.dia) / 2, L = km(target.thick || 0), nz = L > 1 ? 3 : 1, nr = 12, nth = 28;
    let wsum = 0; for (let ir = 0; ir < nr; ir++) wsum += ir + 0.5;
    for (let iz = 0; iz < nz; iz++) {
      const z = nz === 1 ? 0 : -L / 2 + L * (iz + 0.5) / nz;
      for (let ir = 0; ir < nr; ir++) {
        const rr = Rr * (ir + 0.5) / nr, ringM = (target.mass / nz) * (ir + 0.5) / wsum, dm = ringM / nth;
        for (let it = 0; it < nth; it++) {
          const th = 2 * Math.PI * (it + 0.5) / nth;
          const r = local([rr * Math.cos(th), rr * Math.sin(th), z]);
          if (inOther(r)) valid = false;
          add(r, dm);
        }
      }
    }
  } else if (target.type === 'cylinder') {
    const Rr = km(target.dia) / 2, L = km(target.len), nz = 6, nr = 8, nth = 24;
    let wsum = 0; for (let ir = 0; ir < nr; ir++) wsum += ir + 0.5;
    for (let iz = 0; iz < nz; iz++) {
      const z = -L / 2 + L * (iz + 0.5) / nz;
      for (let ir = 0; ir < nr; ir++) {
        const rr = Rr * (ir + 0.5) / nr, ringM = (target.mass / nz) * (ir + 0.5) / wsum, dm = ringM / nth;
        for (let it = 0; it < nth; it++) {
          const th = 2 * Math.PI * (it + 0.5) / nth;
          const r = local([rr * Math.cos(th), rr * Math.sin(th), z]);
          if (inOther(r)) valid = false;
          add(r, dm);
        }
      }
    }
  } else if (target.type === 'box') {
    const a = km(target.size[0]) / 2, b = km(target.size[1]) / 2, cc = km(target.size[2]) / 2;
    const nx = 6, ny = 6, nz = 6, dm = target.mass / (nx * ny * nz);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) for (let k = 0; k < nz; k++) {
      const r = local([
        -a + (i + 0.5) * (2 * a / nx),
        -b + (j + 0.5) * (2 * b / ny),
        -cc + (k + 0.5) * (2 * cc / nz),
      ]);
      if (inOther(r)) valid = false;
      add(r, dm);
    }
  }
  return { F, tau, valid, hasExternal: true, Fabs, tauAbs };
}

// Approximate physical half-extent of a source [m] — its bounding-sphere radius
// about the centre. Used for click hit-testing, the selection ring, and Fit.
export function sourceExtent(s) {
  switch (s.type) {
    case 'box':      return 0.5 * Math.hypot(km(s.size[0]), km(s.size[1]), km(s.size[2]));
    case 'sphere':
    case 'shell':    return km(s.dia) / 2;
    case 'ring':     return km(s.dia) / 2;
    case 'disc':     return Math.max(km(s.dia) / 2, km(s.thick || 0) / 2);
    case 'cylinder': return Math.max(km(s.dia) / 2, km(s.len) / 2);
    case 'rod':      return km(s.len) / 2;
    default:         return km(2000);   // point mass — small marker (~2000 km)
  }
}
