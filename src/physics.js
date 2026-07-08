// physics.js — core Newtonian gravitation, SI units throughout.
// Length: metres, mass: kg, field g: m/s², potential Φ: J/kg (= m²/s²).
//
// All field routines return the gravitational acceleration g [m/s²] at a point
// given in metres — the field a unit test mass would feel, always pointing
// TOWARD the source mass (gravity only attracts). A companion routine returns
// the scalar potential Φ. Sources compose linearly (superposition), so the
// total field/potential of a scene is the sum over its members.

export const G   = 6.67430e-11;     // gravitational constant [m³·kg⁻¹·s⁻²]
export const G0  = 9.80665;         // standard gravity        [m/s²]  (unit ref)

// ---------------------------------------------------------------------------
// Minimal 3-vector helpers.  Vectors are plain [x, y, z] arrays.
// ---------------------------------------------------------------------------
export const vadd   = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const vsub   = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const vscale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const vdot   = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const vcross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
export const vlen  = (a) => Math.hypot(a[0], a[1], a[2]);
export const vnorm = (a) => { const l = vlen(a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
export const vzero = () => [0, 0, 0];

// ---------------------------------------------------------------------------
// Rotations.  Orientation is stored as intrinsic Z-Y-X Euler angles (radians):
// yaw (about z), pitch (about y), roll (about x).  These helpers convert between
// world and a body's local frame.  (Gravity of a symmetric body is orientation-
// independent, but rods, rings, discs, cylinders and boxes are not.)
// ---------------------------------------------------------------------------
export function eulerToMatrix(yaw, pitch, roll) {
  const cy = Math.cos(yaw),   sy = Math.sin(yaw);
  const cp = Math.cos(pitch), sp = Math.sin(pitch);
  const cr = Math.cos(roll),  sr = Math.sin(roll);
  // R = Rz(yaw) · Ry(pitch) · Rx(roll)
  return [
    [cy * cp,  cy * sp * sr - sy * cr,  cy * sp * cr + sy * sr],
    [sy * cp,  sy * sp * sr + cy * cr,  sy * sp * cr - cy * sr],
    [-sp,      cp * sr,                 cp * cr],
  ];
}
export function matVec(M, v) {
  return [
    M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
    M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
    M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
  ];
}
export function matTVec(M, v) { // transpose (inverse for rotations) times vector
  return [
    M[0][0] * v[0] + M[1][0] * v[1] + M[2][0] * v[2],
    M[0][1] * v[0] + M[1][1] * v[1] + M[2][1] * v[2],
    M[0][2] * v[0] + M[1][2] * v[1] + M[2][2] * v[2],
  ];
}

// ---------------------------------------------------------------------------
// Point mass — the fundamental source.
//   g(r) = −G M (r − r₀) / |r − r₀|³      [m/s²]  (points toward the mass)
//   Φ(r) = −G M / |r − r₀|                [J/kg]
// ---------------------------------------------------------------------------
export function pointField(M, r0, Q) {
  const r = vsub(Q, r0);
  const rl = vlen(r);
  if (rl < 1e-12) return [0, 0, 0];
  return vscale(r, -G * M / (rl * rl * rl));
}
export function pointPot(M, r0, Q) {
  const rl = vlen(vsub(Q, r0));
  return rl < 1e-12 ? -Infinity : -G * M / rl;
}

// ---------------------------------------------------------------------------
// Uniform solid sphere — exact everywhere via Newton's shell theorem.
//   Outside (r ≥ R):  a point mass at the centre.
//   Inside  (r < R):  g = −G M r / R³  (linear; only the enclosed mass pulls),
//                     Φ = −G M (3R² − r²) / (2R³).
// The two forms agree at the surface, so g and Φ are continuous there.
// ---------------------------------------------------------------------------
export function sphereField(M, R, r0, Q) {
  const r = vsub(Q, r0);
  const rl = vlen(r);
  if (rl < 1e-12) return [0, 0, 0];
  if (rl >= R) return vscale(r, -G * M / (rl * rl * rl));
  return vscale(r, -G * M / (R * R * R));            // inside: ∝ r
}
export function spherePot(M, R, r0, Q) {
  const rl = vlen(vsub(Q, r0));
  if (rl >= R) return -G * M / rl;
  return -G * M * (3 * R * R - rl * rl) / (2 * R * R * R);
}

// ---------------------------------------------------------------------------
// Thin spherical shell — exact everywhere.
//   Outside (r ≥ R):  a point mass at the centre.
//   Inside  (r < R):  g = 0  (the shell exerts no net field within it),
//                     Φ = −G M / R  (constant — a flat potential well floor).
// ---------------------------------------------------------------------------
export function shellField(M, R, r0, Q) {
  const r = vsub(Q, r0);
  const rl = vlen(r);
  if (rl < R) return [0, 0, 0];
  if (rl < 1e-12) return [0, 0, 0];
  return vscale(r, -G * M / (rl * rl * rl));
}
export function shellPot(M, R, r0, Q) {
  const rl = vlen(vsub(Q, r0));
  return rl >= R ? -G * M / rl : -G * M / R;
}

// ---------------------------------------------------------------------------
// Complete elliptic integrals K(m), E(m) (parameter m = k²) via the
// arithmetic–geometric mean.  Fast and accurate for 0 ≤ m < 1.
// ---------------------------------------------------------------------------
export function ellipKE(m) {
  if (m < 0) m = 0;
  if (m > 1 - 1e-12) m = 1 - 1e-12;
  let a = 1, b = Math.sqrt(1 - m), c = Math.sqrt(m);
  let s = 0.5 * c * c;      // n = 0 term
  let twoPow = 1;
  for (let n = 1; n <= 16; n++) {
    const a1 = (a + b) / 2, b1 = Math.sqrt(a * b), c1 = (a - b) / 2;
    a = a1; b = b1; c = c1;
    twoPow *= 2;
    s += (twoPow / 2) * c * c;
    if (c < 1e-15) break;
  }
  const K = Math.PI / (2 * a);
  return [K, K * (1 - s)];
}

// ---------------------------------------------------------------------------
// Uniform circular ring (hoop) of mass M, radius a, lying in the local z = 0
// plane centred at the origin — evaluated at a local point (x, y, z).
// Closed form via the complete elliptic integrals K, E (the gravitational
// analogue of the current-loop field).  One evaluation replaces summing dozens
// of point masses around the hoop: much faster and exact.
//
//   D² = (ρ+a)² + z²,   α² = (ρ−a)² + z²,   m = 4aρ / D²
//   g_z = −(2 G M z)/(π D α²) · E(m)
//   g_ρ = −∂Φ/∂ρ   from Φ = −(2 G M/π) K(m)/D   (see docs/PHYSICS.md)
// ---------------------------------------------------------------------------
export function ringField(a, M, x, y, z) {
  const rho = Math.hypot(x, y);
  const D2 = (rho + a) * (rho + a) + z * z;
  const D = Math.sqrt(D2) || 1e-30;
  let al2 = (rho - a) * (rho - a) + z * z;
  if (al2 < 1e-24) al2 = 1e-24;                 // guard exactly on the hoop
  const m = 4 * a * rho / D2;
  const [K, E] = ellipKE(m);
  const gz = -(2 * G * M * z) / (Math.PI * D * al2) * E;
  let grho = 0;
  if (rho > 1e-12) {
    const om = al2 / D2;                          // 1 − m
    const Kp = (E - om * K) / (2 * m * om);       // dK/dm
    const m_u = 4 * a * (a * a + z * z - rho * rho) / (D2 * D2);
    const D_u = (rho + a) / D;
    grho = (2 * G * M / Math.PI) * (Kp * m_u / D - K * D_u / D2);
  }
  const bx = rho > 1e-12 ? grho * x / rho : 0;
  const by = rho > 1e-12 ? grho * y / rho : 0;
  return [bx, by, gz];
}
export function ringPot(a, M, x, y, z) {
  const rho = Math.hypot(x, y);
  const D2 = (rho + a) * (rho + a) + z * z;
  const D = Math.sqrt(D2) || 1e-30;
  const m = 4 * a * rho / D2;
  const [K] = ellipKE(m);
  return -(2 * G * M / Math.PI) * K / D;
}

// ---------------------------------------------------------------------------
// Uniform straight rod — exact closed form.
//
// A thin rod of linear density λ [kg/m] and length L lies along the local z
// axis from −L/2 to +L/2.  At a local point with cylindrical radius ρ and axial
// coordinate z (distances r_A, r_B to the two ends):
//
//   g_ρ = −(Gλ/ρ) [ (z+L/2)/r_A − (z−L/2)/r_B ]
//   g_z =  Gλ [ 1/r_A − 1/r_B ]
//   Φ   = −Gλ ln[ (z+L/2 + r_A) / (z−L/2 + r_B) ]
//
// (r_A is the distance to the −L/2 end, r_B to the +L/2 end.)  The infinite-rod
// limit is g_ρ → −2Gλ/ρ.
// ---------------------------------------------------------------------------
export function rodField(lam, L, x, y, z) {
  const rho = Math.hypot(x, y);
  const rA = Math.hypot(rho, z + L / 2) || 1e-30;
  const rB = Math.hypot(rho, z - L / 2) || 1e-30;
  const gz = G * lam * (1 / rA - 1 / rB);
  let grho = 0;
  if (rho > 1e-12) grho = -(G * lam / rho) * ((z + L / 2) / rA - (z - L / 2) / rB);
  const gx = rho > 1e-12 ? grho * x / rho : 0;
  const gy = rho > 1e-12 ? grho * y / rho : 0;
  return [gx, gy, gz];
}
export function rodPot(lam, L, x, y, z) {
  const rho = Math.hypot(x, y);
  const rA = Math.hypot(rho, z + L / 2);
  const rB = Math.hypot(rho, z - L / 2);
  const num = (z + L / 2) + rA, den = (z - L / 2) + rB;
  return -G * lam * Math.log((num / den) || 1e-300);
}

// ---------------------------------------------------------------------------
// Uniform rectangular prism (box) of density ρ_d [kg/m³], half-extents
// (a, b, c), centred at the local origin — exact closed form (Nagy 1966/2000,
// the standard result in gravimetry / geophysics).  Evaluated at local point p.
//
// With the eight corners X_i = pₓ±a, Y_j = p_y±b, Z_k = p_z±c and r the distance
// to each corner, the field components are the signed corner sum
//   g_x = G ρ_d Σ (−1)^{i+j+k} [ Y ln(Z+r) + Z ln(Y+r) − X·atan2(YZ, Xr) ]
// and cyclic permutations for g_y, g_z.  Exact everywhere outside the body.
// ---------------------------------------------------------------------------
export function cuboidField(p, half, rhod) {
  const X = [p[0] + half[0], p[0] - half[0]];
  const Y = [p[1] + half[1], p[1] - half[1]];
  const Z = [p[2] + half[2], p[2] - half[2]];
  let gx = 0, gy = 0, gz = 0;
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) {
    const s = ((i + j + k) & 1) ? -1 : 1;
    const x = X[i], y = Y[j], z = Z[k];
    const r = Math.hypot(x, y, z) || 1e-30;
    gx += s * (y * Math.log(z + r) + z * Math.log(y + r) - x * Math.atan2(y * z, x * r));
    gy += s * (z * Math.log(x + r) + x * Math.log(z + r) - y * Math.atan2(z * x, y * r));
    gz += s * (x * Math.log(y + r) + y * Math.log(x + r) - z * Math.atan2(x * y, z * r));
  }
  const kf = G * rhod;
  return [kf * gx, kf * gy, kf * gz];
}
export function cuboidPot(p, half, rhod) {
  const X = [p[0] + half[0], p[0] - half[0]];
  const Y = [p[1] + half[1], p[1] - half[1]];
  const Z = [p[2] + half[2], p[2] - half[2]];
  let s0 = 0;
  for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) {
    const s = ((i + j + k) & 1) ? -1 : 1;
    const x = X[i], y = Y[j], z = Z[k];
    const r = Math.hypot(x, y, z) || 1e-30;
    s0 += s * (
      x * y * Math.log(z + r) + y * z * Math.log(x + r) + z * x * Math.log(y + r)
      - 0.5 * x * x * Math.atan2(y * z, x * r)
      - 0.5 * y * y * Math.atan2(z * x, y * r)
      - 0.5 * z * z * Math.atan2(x * y, z * r)
    );
  }
  return -G * rhod * s0;   // Φ = −G ρ ∫dV/r
}

// ---------------------------------------------------------------------------
// Velocity-Verlet (leapfrog) step for a test mass under acceleration a(x).
// Symplectic and time-reversible: it conserves orbital energy over long runs,
// so closed orbits stay closed — the gravitational counterpart of the Boris
// pusher.  accelFn(x) returns the acceleration g [m/s²] at position x.
// ---------------------------------------------------------------------------
export function verletStep(x, v, dt, accelFn, a0) {
  const a = a0 || accelFn(x);
  const xNew = [
    x[0] + v[0] * dt + 0.5 * a[0] * dt * dt,
    x[1] + v[1] * dt + 0.5 * a[1] * dt * dt,
    x[2] + v[2] * dt + 0.5 * a[2] * dt * dt,
  ];
  const a1 = accelFn(xNew);
  const vNew = [
    v[0] + 0.5 * (a[0] + a1[0]) * dt,
    v[1] + 0.5 * (a[1] + a1[1]) * dt,
    v[2] + 0.5 * (a[2] + a1[2]) * dt,
  ];
  return { x: xNew, v: vNew, a: a1 };
}

// ---------------------------------------------------------------------------
// Full mutual N-body gravity between bodies, each { x:[3], v:[3], mass, rad }.
// Every body pulls on every other: a_i = −G Σ_j m_j (x_i−x_j)/|x_i−x_j|³.
// The pair law is EXACT (Newtonian 1/r²) while the bodies don't touch
// (d ≥ r_i+r_j) — which, by the shell theorem, is exact for spheres. Once they
// overlap it crosses over to the uniform-sphere interior law (∝ d) so a direct
// hit stays finite instead of exploding. The pair force is equal-and-opposite,
// so total momentum is conserved to machine precision.
// ---------------------------------------------------------------------------
export function nbodyAccel(bodies) {
  const n = bodies.length;
  const a = new Array(n);
  for (let i = 0; i < n; i++) a[i] = [0, 0, 0];
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    const bi = bodies[i], bj = bodies[j];
    const dx = bi.x[0] - bj.x[0], dy = bi.x[1] - bj.x[1], dz = bi.x[2] - bj.x[2];
    const d = Math.hypot(dx, dy, dz) || 1e-30;
    const R = Math.max((bi.rad || 0) + (bj.rad || 0), 1e3);   // contact scale (min 1 km guard)
    const inv = d >= R ? 1 / (d * d * d) : 1 / (R * R * R);   // exact outside; bounded within contact
    const gi = G * bj.mass * inv, gj = G * bi.mass * inv;
    a[i][0] -= gi * dx; a[i][1] -= gi * dy; a[i][2] -= gi * dz;
    a[j][0] += gj * dx; a[j][1] += gj * dy; a[j][2] += gj * dz;
  }
  return a;
}
// One velocity-Verlet step of the whole system (mutates each body's x and v).
export function nbodyStep(bodies, dt) {
  const a0 = nbodyAccel(bodies);
  for (let i = 0; i < bodies.length; i++) {
    const b = bodies[i], a = a0[i];
    b.x[0] += b.v[0] * dt + 0.5 * a[0] * dt * dt;
    b.x[1] += b.v[1] * dt + 0.5 * a[1] * dt * dt;
    b.x[2] += b.v[2] * dt + 0.5 * a[2] * dt * dt;
  }
  const a1 = nbodyAccel(bodies);
  for (let i = 0; i < bodies.length; i++) {
    const b = bodies[i];
    b.v[0] += 0.5 * (a0[i][0] + a1[i][0]) * dt;
    b.v[1] += 0.5 * (a0[i][1] + a1[i][1]) * dt;
    b.v[2] += 0.5 * (a0[i][2] + a1[i][2]) * dt;
  }
}
