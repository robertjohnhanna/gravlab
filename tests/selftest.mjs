// Physics verification. Each field routine is checked against an independent,
// analytically-known limit and against Maxwell/Poisson constraints in free space
// (∇·g = 0 and ∇×g = 0 outside matter), and every potential is checked to satisfy
// g = −∇Φ. Run with `npm test`.
import {
  G, pointField, pointPot, sphereField, spherePot, shellField, shellPot,
  ringField, ringPot, rodField, rodPot, cuboidField, cuboidPot, verletStep,
  nbodyStep, vadd, vsub, vscale, vlen, vcross,
} from '../src/physics.js';
import { Scene, defaultSource, sourceField, forceOn, massOf, BODIES, BODY_DIA } from '../src/sources.js';

let passed = 0, failed = 0;
function check(name, cond, extra = '') {
  if (cond) { passed++; console.log(`  ok  ${name}`); }
  else { failed++; console.log(`FAIL  ${name}  ${extra}`); }
}
const rel = (a, b) => Math.abs(a - b) / (Math.abs(b) + 1e-300);

// Numerical divergence & curl of a vector field fn(x)->[gx,gy,gz] at point p.
function divergence(fn, p, h) {
  let d = 0;
  for (let i = 0; i < 3; i++) {
    const pp = p.slice(); pp[i] += h; const pm = p.slice(); pm[i] -= h;
    d += (fn(pp)[i] - fn(pm)[i]) / (2 * h);
  }
  return d;
}
function curl(fn, p, h) {
  const d = (i, j) => { const pp = p.slice(); pp[j] += h; const pm = p.slice(); pm[j] -= h; return (fn(pp)[i] - fn(pm)[i]) / (2 * h); };
  return [d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)];
}
// g should equal −∇Φ.  Compare over several points.
function checkGradient(name, gfn, pfn, pts, tol) {
  let worst = 0;
  for (const p of pts) {
    const h = 1e-4 * (vlen(p) || 1);
    const gnum = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      const pp = p.slice(); pp[i] += h; const pm = p.slice(); pm[i] -= h;
      gnum[i] = -(pfn(pp) - pfn(pm)) / (2 * h);
    }
    const g = gfn(p);
    worst = Math.max(worst, vlen(vsub(g, gnum)) / (vlen(g) + 1e-300));
  }
  check(name, worst < tol, `worst rel=${worst.toExponential(2)}`);
}

console.log('\n== Point mass ==');
{
  const M = 5e24, r0 = [0, 0, 0], Q = [3e6, 0, 0];
  const g = pointField(M, r0, Q);
  check('|g| = GM/r²', rel(vlen(g), G * M / (3e6 * 3e6)) < 1e-12, vlen(g));
  check('points toward the mass (−x)', g[0] < 0 && Math.abs(g[1]) < 1e-20, g.toString());
  check('Φ = −GM/r', rel(pointPot(M, r0, Q), -G * M / 3e6) < 1e-12, pointPot(M, r0, Q));
}

console.log('\n== Uniform solid sphere (shell theorem) ==');
{
  const M = 6e24, R = 6e6, r0 = [0, 0, 0];
  const out = sphereField(M, R, r0, [1e7, 0, 0]);
  check('outside = point mass', rel(vlen(out), G * M / 1e14) < 1e-12, vlen(out));
  const inR = 3e6;
  const ins = sphereField(M, R, r0, [inR, 0, 0]);
  check('inside g ∝ r  (g = GM r/R³)', rel(vlen(ins), G * M * inR / (R ** 3)) < 1e-12, vlen(ins));
  const eIn = sphereField(M, R, r0, [R - 1, 0, 0]), eOut = sphereField(M, R, r0, [R + 1, 0, 0]);
  check('g continuous across surface', rel(vlen(eIn), vlen(eOut)) < 1e-5, `${vlen(eIn)} vs ${vlen(eOut)}`);
  check('Φ continuous across surface', rel(spherePot(M, R, r0, [R - 1, 0, 0]), spherePot(M, R, r0, [R + 1, 0, 0])) < 1e-5, '');
  checkGradient('sphere g = −∇Φ (in & out)', (p) => sphereField(M, R, r0, p), (p) => spherePot(M, R, r0, p),
    [[3e6, 0, 0], [1e7, 2e6, 0], [0, 0, 4e6]], 1e-4);
}

console.log('\n== Thin spherical shell ==');
{
  const M = 4e24, R = 7e6, r0 = [0, 0, 0];
  check('field zero inside', vlen(shellField(M, R, r0, [2e6, 1e6, 0])) === 0, '');
  check('outside = point mass', rel(vlen(shellField(M, R, r0, [2e7, 0, 0])), G * M / 4e14) < 1e-12, '');
  check('Φ inside is constant −GM/R', rel(shellPot(M, R, r0, [1e6, 0, 0]), -G * M / R) < 1e-12, '');
  check('Φ continuous at surface', rel(shellPot(M, R, r0, [R - 1, 0, 0]), shellPot(M, R, r0, [R + 1, 0, 0])) < 1e-5, '');
}

console.log('\n== Uniform ring (elliptic integrals) ==');
{
  const a = 1e7, M = 2e24;
  // on-axis closed form  g_z = −GM z /(a²+z²)^{3/2}
  for (const z of [3e6, 1e7, 3e7]) {
    const g = ringField(a, M, 0, 0, z);
    const exp = -G * M * z / (a * a + z * z) ** 1.5;
    check(`ring on-axis g_z (z=${z})`, rel(g[2], exp) < 1e-6, `${g[2]} vs ${exp}`);
    check(`ring on-axis transverse ≈ 0 (z=${z})`, Math.abs(g[0]) < 1e-15 && Math.abs(g[1]) < 1e-15, g.toString());
  }
  // off-axis vs brute-force point-mass integration around the hoop
  const brute = (x, y, z, N = 20000) => {
    let gx = 0, gy = 0, gz = 0;
    for (let i = 0; i < N; i++) { const ph = 2 * Math.PI * (i + 0.5) / N; const px = a * Math.cos(ph), py = a * Math.sin(ph); const dx = x - px, dy = y - py, dz = z, r = Math.hypot(dx, dy, dz), r3 = r ** 3, dm = M / N; gx += -G * dm * dx / r3; gy += -G * dm * dy / r3; gz += -G * dm * dz / r3; }
    return [gx, gy, gz];
  };
  for (const Q of [[5e6, 0, 3e6], [1.4e7, 2e6, -4e6], [3e6, 3e6, 8e6]]) {
    const gc = ringField(a, M, ...Q), gb = brute(...Q);
    check(`ring off-axis matches integration @${Q}`, vlen(vsub(gc, gb)) / vlen(gb) < 2e-3, `c=${gc.map((v) => v.toExponential(2))} b=${gb.map((v) => v.toExponential(2))}`);
  }
  checkGradient('ring g = −∇Φ', (p) => ringField(a, M, ...p), (p) => ringPot(a, M, ...p),
    [[5e6, 0, 3e6], [1.4e7, 2e6, -4e6], [3e6, 3e6, 8e6]], 2e-3);
  // Just off the axis the transverse field must vanish ∝ ρ — the naive dK/dm
  // closed form loses all significance there (E−K is pure roundoff at m→0).
  const gnx = ringField(a, M, 1e-6, 0, 5e6);
  check('ring transverse field ∝ ρ near axis (no cancellation noise)',
    Math.abs(gnx[0]) < 1e-8 * Math.abs(gnx[2]), `gx/gz=${(gnx[0] / gnx[2]).toExponential(2)}`);
}

console.log('\n== Uniform rod (line mass) ==');
{
  const lam = 2e19, L = 4e7;                       // λ [kg/m]
  const d = 5e6;
  // long-rod (near the middle, far from ends) → infinite-rod limit 2Gλ/d.
  // Use a rod far longer than the standoff so the finite-length correction
  // [(L/2)/√(d²+(L/2)²)] is negligible.
  const gInf = rodField(lam, 4e8, 1e6, 0, 0);
  check('mid-rod radial → 2Gλ/d', rel(Math.abs(gInf[0]), 2 * G * lam / 1e6) < 2e-3, `${Math.abs(gInf[0])} vs ${2 * G * lam / 1e6}`);
  const g = rodField(lam, L, d, 0, 0);
  check('mid-rod axial ≈ 0 by symmetry', Math.abs(g[2]) < 1e-20, g[2]);
  // vs brute-force integration
  const brute = (x, y, z, N = 40000) => { let gx = 0, gz = 0; const dm = lam * L / N; for (let i = 0; i < N; i++) { const zp = -L / 2 + L * (i + 0.5) / N; const dx = x, dz = z - zp, r = Math.hypot(dx, y, dz), r3 = r ** 3; gx += -G * dm * dx / r3; gz += -G * dm * dz / r3; } return [gx, gz]; };
  for (const Q of [[3e6, 0, 1e7], [8e6, 0, -1.5e7], [2e6, 0, 3e7]]) {
    const gc = rodField(lam, L, Q[0], 0, Q[2]), gb = brute(Q[0], 0, Q[2]);
    check(`rod off-axis matches integration @${Q}`, Math.abs(gc[0] - gb[0]) / Math.abs(gb[0]) < 1e-3 && Math.abs(gc[2] - gb[1]) / (Math.abs(gb[1]) + 1e-30) < 2e-3, '');
  }
  checkGradient('rod g = −∇Φ', (p) => rodField(lam, L, ...p), (p) => rodPot(lam, L, ...p),
    [[3e6, 0, 1e7], [8e6, 1e6, -1.5e7], [2e6, 2e6, 3e7]], 2e-3);
  // On the axis beyond the ends the naive ln[(z+L/2+rA)/(z−L/2+rB)] form
  // cancels to 0/0 and flips the sign of Φ; the stable form must not.
  const bruteP = (x, y, z, N = 40000) => { let p = 0; const dm = lam * L / N; for (let i = 0; i < N; i++) { const zp = -L / 2 + L * (i + 0.5) / N; p += -G * dm / Math.hypot(x, y, z - zp); } return p; };
  for (const z of [-3e7, 3e7]) {
    const pOn = rodPot(lam, L, 0, 0, z);
    check(`rod Φ finite & negative on axis beyond end (z=${z})`, isFinite(pOn) && pOn < 0, pOn);
    check(`rod Φ on-axis matches integration (z=${z})`, rel(pOn, bruteP(0, 0, z)) < 1e-3, `${pOn} vs ${bruteP(0, 0, z)}`);
    check(`rod Φ continuous across the axis (z=${z})`, rel(pOn, rodPot(lam, L, 1, 0, z)) < 1e-6, '');
  }
}

console.log('\n== Uniform box (Nagy closed form) ==');
{
  const half = [4e6, 6e6, 3e6], rhod = 4000;
  const M = 8 * half[0] * half[1] * half[2] * rhod;
  const fn = (p) => cuboidField(p, half, rhod);
  // far field → point mass
  const Rf = 3e8, gf = fn([Rf, 0, 0]);
  check('far field → point mass GM/r²', rel(vlen(gf), G * M / (Rf * Rf)) < 1e-3, `${vlen(gf)} vs ${G * M / (Rf * Rf)}`);
  check('far field points inward', gf[0] < 0, gf[0]);
  // Laplace: ∇·g = 0 and ∇×g = 0 outside the body
  let maxDiv = 0, maxCurl = 0;
  const pts = [[8e6, 0, 0], [7e6, 8e6, 0], [1e7, -5e6, 4e6], [0, 1e7, -5e6]];
  const scale = vlen(fn([8e6, 0, 0])) / 4e6;
  for (const p of pts) { maxDiv = Math.max(maxDiv, Math.abs(divergence(fn, p, 1e3))); maxCurl = Math.max(maxCurl, vlen(curl(fn, p, 1e3))); }
  check('∇·g ≈ 0 outside box', maxDiv / scale < 1e-4, `div/ref=${(maxDiv / scale).toExponential(2)}`);
  check('∇×g ≈ 0 outside box', maxCurl / scale < 1e-4, `curl/ref=${(maxCurl / scale).toExponential(2)}`);
  checkGradient('box g = −∇Φ', fn, (p) => cuboidPot(p, half, rhod), pts, 1e-3);
  // On a line extending a box edge two corner terms hit 0·ln(0): the true limit
  // is 0, but unguarded IEEE arithmetic yields NaN and poisons the whole sum.
  const ge = fn([half[0], half[1], -8e6]);
  const gnear = fn([half[0] + 0.5, half[1] + 0.5, -8e6]);
  check('box field finite on edge-extension line', ge.every(isFinite), ge.toString());
  check('box edge-line field matches neighbour', vlen(vsub(ge, gnear)) / vlen(gnear) < 1e-3, '');
  check('box Φ finite on edge-extension line', isFinite(cuboidPot([half[0], half[1], -8e6], half, rhod)), '');
}

console.log('\n== Disc & cylinder (ring stacks) ==');
{
  // Thin disc on-axis vs analytic  g_z = −2πGσ(1 − z/√(z²+R²))
  const disc = defaultSource('disc'); disc.dia = 30000; disc.mass = 3e24; disc.thick = 0;
  const scD = new Scene(); scD.add(disc);
  const Rm = 15000 * 1000, sigma = disc.mass / (Math.PI * Rm * Rm);
  for (const z of [5e6, 2e7]) {
    const g = scD.g([0, 0, z]);
    const exp = -2 * Math.PI * G * sigma * (1 - z / Math.sqrt(z * z + Rm * Rm));
    check(`disc on-axis g_z (z=${z})`, rel(g[2], exp) < 2e-2, `${g[2]} vs ${exp}`);
  }
  // A thick disc far away still reduces to a point mass.
  const td = defaultSource('disc'); td.dia = 20000; td.thick = 8000; td.mass = 2e24;
  const scT = new Scene(); scT.add(td);
  check('thick disc far field → point mass', rel(vlen(scT.g([5e8, 0, 0])), G * td.mass / (5e8) ** 2) < 5e-3, '');
  // Cylinder far field → point mass
  const cyl = defaultSource('cylinder'); cyl.dia = 14000; cyl.len = 26000; cyl.mass = 3e24;
  const scC = new Scene(); scC.add(cyl);
  const gf = scC.g([4e8, 0, 0]);
  check('cylinder far field → point mass', rel(vlen(gf), G * cyl.mass / (4e8) ** 2) < 5e-3, `${vlen(gf)} vs ${G * cyl.mass / (4e8) ** 2}`);
  // The single-pass sample() and the field-only fast path must agree exactly
  // with the separate g()/potential() calls.
  const Q = [3e7, 1e7, 5e6];
  const smp = scC.sample(Q);
  check('Scene.sample g equals Scene.g', vlen(vsub(smp.g, scC.g(Q))) === 0, '');
  check('Scene.sample Φ equals Scene.potential', smp.phi === scC.potential(Q), '');
}

console.log('\n== Exact force / torque ==');
{
  // Two spheres on the x-axis → point-mass attraction GM₁M₂/d².
  const A = defaultSource('sphere'); A.dia = 8000; A.mass = 5e24; A.pos = [-20000, 0, 0];
  const B = defaultSource('sphere'); B.dia = 8000; B.mass = 3e24; B.pos = [20000, 0, 0];
  const sc = new Scene(); sc.add(A); sc.add(B);
  const fA = forceOn(sc, A), fB = forceOn(sc, B);
  const d = 40000 * 1000;
  const expF = G * A.mass * B.mass / (d * d);
  check('sphere–sphere = GM₁M₂/d²', rel(vlen(fA.F), expF) < 1e-9, `${vlen(fA.F)} vs ${expF}`);
  check('force is attractive (A pulled +x toward B)', fA.F[0] > 0 && Math.abs(fA.F[1]) < 1e-6 * vlen(fA.F), fA.F.toString());
  const net = vlen(vadd(fA.F, fB.F));
  check("Newton's third law (F_AB = −F_BA)", net / vlen(fA.F) < 1e-9, net);
  check('sphere feels zero self-torque', vlen(fA.tau) === 0, '');

  const lone = new Scene(); lone.add(defaultSource('sphere'));
  check('lone body feels no force', forceOn(lone, lone.sources[0]).hasExternal === false, '');

  // Overlapping bodies → force refused, not a bogus value.
  const ov = new Scene();
  const o1 = defaultSource('sphere'); o1.dia = 20000; o1.pos = [-4000, 0, 0]; ov.add(o1);
  const o2 = defaultSource('sphere'); o2.dia = 20000; o2.pos = [4000, 0, 0]; ov.add(o2);
  check('overlapping bodies → force refused', forceOn(ov, o2).valid === false, '');

  // Gravity-gradient torque: a rod tilted 45° near a point mass feels a real
  // restoring torque toward radial alignment (the tidal-locking mechanism).
  const gg = new Scene();
  const planet = defaultSource('point'); planet.mass = 6e24; planet.pos = [0, 0, 0]; gg.add(planet);
  const rod = defaultSource('rod'); rod.len = 30000; rod.mass = 1e20; rod.pos = [60000, 0, 0]; rod.rot = [0, 45, 0]; gg.add(rod);
  const fg = forceOn(gg, rod);
  check('elongated body in a gradient feels nonzero torque', vlen(fg.tau) > 0, vlen(fg.tau));
  // A rod pointed straight at the mass (radial) is at a torque equilibrium.
  const rod2 = defaultSource('rod'); rod2.len = 30000; rod2.mass = 1e20; rod2.pos = [60000, 0, 0]; rod2.rot = [0, 90, 0];
  const gg2 = new Scene(); gg2.add(planet); gg2.add(rod2);
  const fg2 = forceOn(gg2, rod2);
  check('radial rod is at torque equilibrium (τ ≈ 0)', vlen(fg2.tau) / (vlen(fg2.F) * 30000e3) < 1e-6, vlen(fg2.tau));
}

console.log('\n== Orbit integrator (velocity Verlet) ==');
{
  // Circular orbit of a test mass about a point mass: check the radius stays
  // constant, the period matches Kepler T = 2π√(r³/GM), and energy is conserved.
  const M = 6e24, r = 2e7;
  const vc = Math.sqrt(G * M / r);                 // circular speed
  const accel = (x) => { const rl = vlen(x); return vscale(x, -G * M / (rl ** 3)); };
  let x = [r, 0, 0], v = [0, vc, 0];
  const Tkep = 2 * Math.PI * Math.sqrt(r ** 3 / (G * M));
  const dt = Tkep / 4000;
  let rMin = Infinity, rMax = 0, eMin = Infinity, eMax = -Infinity;
  const nSteps = 4000;                             // one full period
  let a0 = accel(x);
  for (let n = 0; n < nSteps; n++) {
    const s = verletStep(x, v, dt, accel, a0); x = s.x; v = s.v; a0 = s.a;
    const rl = vlen(x); rMin = Math.min(rMin, rl); rMax = Math.max(rMax, rl);
    const e = 0.5 * vlen(v) ** 2 - G * M / rl; eMin = Math.min(eMin, e); eMax = Math.max(eMax, e);
  }
  check('circular orbit keeps constant radius', (rMax - rMin) / r < 1e-3, `Δr/r=${((rMax - rMin) / r).toExponential(2)}`);
  check('orbit closes after one Kepler period', vlen(vsub(x, [r, 0, 0])) / r < 5e-3, `miss=${(vlen(vsub(x, [r, 0, 0])) / r).toExponential(2)}`);
  check('energy conserved', (eMax - eMin) / Math.abs(eMin) < 1e-4, `ΔE/E=${((eMax - eMin) / Math.abs(eMin)).toExponential(2)}`);
}

console.log('\n== Orbit energetics (escape & vis-viva) ==');
{
  const M = 6e24, r = 2e7;
  const vesc = Math.sqrt(2 * G * M / r);
  const eps = 0.5 * vesc * vesc - G * M / r;
  check('escape speed → ε ≈ 0 (parabolic)', Math.abs(eps) / (G * M / r) < 1e-12, eps);

  // Sub-circular tangential launch → bound ellipse; the measured semi-major
  // axis (r_peri+r_apo)/2 must equal −GM/2ε (vis-viva / energy relation).
  const accel = (x) => { const rl = vlen(x); return vscale(x, -G * M / (rl ** 3)); };
  let x = [r, 0, 0], v = [0, 0.8 * Math.sqrt(G * M / r), 0];
  const eps0 = 0.5 * vlen(v) ** 2 - G * M / vlen(x);
  const aSemi = -G * M / (2 * eps0);
  const Tell = 2 * Math.PI * Math.sqrt(aSemi ** 3 / (G * M));
  const dt = Tell / 6000; let rmin = Infinity, rmax = 0, a0 = accel(x);
  for (let n = 0; n < 6000; n++) { const s = verletStep(x, v, dt, accel, a0); x = s.x; v = s.v; a0 = s.a; const rl = vlen(x); rmin = Math.min(rmin, rl); rmax = Math.max(rmax, rl); }
  check('vis-viva: semi-major axis = −GM/2ε', rel((rmin + rmax) / 2, aSemi) < 5e-3, `${(rmin + rmax) / 2} vs ${aSemi}`);
}

console.log('\n== Extended-body force reduces to point mass ==');
{
  const box = defaultSource('box'); box.size = [10000, 12000, 8000]; box.mass = 3e24; box.pos = [0, 0, 0];
  const sph = defaultSource('sphere'); sph.dia = 6000; sph.mass = 5e24; sph.pos = [600000, 0, 0];   // far away
  const sc = new Scene(); sc.add(box); sc.add(sph);
  const fb = forceOn(sc, box), fs = forceOn(sc, sph);
  const d = 600000 * 1000, expF = G * box.mass * sph.mass / (d * d);
  check('box←sphere far force ≈ GM₁M₂/d²', rel(vlen(fb.F), expF) < 5e-3, `${vlen(fb.F)} vs ${expF}`);
  check('third law holds for the extended body', vlen(vadd(fb.F, fs.F)) / vlen(fb.F) < 5e-3, '');
}

console.log('\n== Named-body presets are self-consistent ==');
{
  for (const key of ['Earth (6.0e24)', 'Jupiter (1.9e27)', 'Sun (2.0e30)']) {
    const s = defaultSource('sphere'); s.mass = BODIES[key]; s.dia = BODY_DIA[key];
    const sc = new Scene(); sc.add(s);
    const R = BODY_DIA[key] * 1000 / 2;
    check(`${key}: surface g = GM/R²`, rel(vlen(sc.g([R, 0, 0])), G * BODIES[key] / (R * R)) < 1e-9, '');
  }
  const e = defaultSource('sphere'); e.mass = BODIES['Earth (6.0e24)']; e.dia = BODY_DIA['Earth (6.0e24)'];
  const sc = new Scene(); sc.add(e);
  const g = vlen(sc.g([BODY_DIA['Earth (6.0e24)'] * 1000 / 2, 0, 0]));
  check('Earth preset surface g ≈ 9.8 m/s²', Math.abs(g - 9.82) < 0.1, g);
}

console.log('\n== N-body mutual gravity ==');
{
  const mom = (bs) => bs.reduce((s, b) => vadd(s, vscale(b.v, b.mass)), [0, 0, 0]);
  const com = (bs) => { let M = 0, c = [0, 0, 0]; for (const b of bs) { M += b.mass; c = vadd(c, vscale(b.x, b.mass)); } return vscale(c, 1 / M); };
  const energy = (bs) => {
    let e = 0;
    for (const b of bs) e += 0.5 * b.mass * vlen(b.v) ** 2;
    for (let i = 0; i < bs.length; i++) for (let j = i + 1; j < bs.length; j++) e -= G * bs[i].mass * bs[j].mass / vlen(vsub(bs[i].x, bs[j].x));
    return e;
  };
  // Equal-mass circular binary about the barycentre.
  const M = 5e24, d = 1e8, vrel = Math.sqrt(G * 2 * M / d), ve = vrel / 2;
  const bodies = [
    { x: [-d / 2, 0, 0], v: [0, ve, 0], mass: M, rad: 5e6 },
    { x: [d / 2, 0, 0], v: [0, -ve, 0], mass: M, rad: 5e6 },
  ];
  const T = 2 * Math.PI * Math.sqrt(d ** 3 / (G * 2 * M)), dt = T / 8000;
  const e0 = energy(bodies);
  let sMin = Infinity, sMax = 0;
  for (let n = 0; n < 8000; n++) { nbodyStep(bodies, dt); const s = vlen(vsub(bodies[0].x, bodies[1].x)); sMin = Math.min(sMin, s); sMax = Math.max(sMax, s); }
  check('binary separation constant (circular)', (sMax - sMin) / d < 3e-3, `Δ/d=${((sMax - sMin) / d).toExponential(2)}`);
  check('binary closes after one Kepler period', vlen(vsub(bodies[0].x, [-d / 2, 0, 0])) / d < 5e-3, '');
  check('barycentre stationary', vlen(com(bodies)) / d < 1e-6, '');
  check('total momentum conserved (≈0)', vlen(mom(bodies)) / (M * ve) < 1e-9, '');
  check('total energy conserved', Math.abs(energy(bodies) - e0) / Math.abs(e0) < 1e-4, '');

  // Asymmetric 3-body with net drift: momentum & energy must still be conserved.
  const tb = [
    { x: [0, 0, 0], v: [0, 0, 0], mass: 6e24, rad: 6e6 },
    { x: [1.2e8, 0, 0], v: [0, 1800, 0], mass: 7e22, rad: 1e6 },
    { x: [0, 9e7, 0], v: [-1500, 0, 0], mass: 2e23, rad: 2e6 },
  ];
  const p0 = mom(tb), en0 = energy(tb);
  for (let n = 0; n < 6000; n++) nbodyStep(tb, 20);
  check('3-body momentum conserved', vlen(vsub(mom(tb), p0)) / (vlen(p0) + 1e-30) < 1e-9, '');
  check('3-body energy conserved', Math.abs(energy(tb) - en0) / Math.abs(en0) < 2e-3, `ΔE/E=${(Math.abs(energy(tb) - en0) / Math.abs(en0)).toExponential(2)}`);

  // A light body barely perturbs a heavy one; a comparable one moves it a lot.
  const recoil = (mProj) => {
    const bs = [{ x: [0, 0, 0], v: [0, 0, 0], mass: 6e24, rad: 6e6 }, { x: [5e7, 0, 0], v: [0, 3000, 0], mass: mProj, rad: 1e5 }];
    for (let n = 0; n < 3000; n++) nbodyStep(bs, 30);
    return vlen(bs[0].v);            // how fast the "planet" ends up moving
  };
  check('heavier flyby moves the planet more', recoil(6e23) > 50 * recoil(6e21), '');
}

console.log(`\n${passed} passed, ${failed} failed\n`);
process.exit(failed ? 1 : 0);
