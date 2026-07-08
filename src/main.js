// main.js — app glue: scene management, UI panels, interaction, orbit sim.
import * as P from './physics.js';
import { Scene, defaultSource, buildSource, sourceExtent, massOf, BODIES, BODY_DIA } from './sources.js';
import { Renderer, View, viridis } from './render.js';

const scene = new Scene();
const view = new View();
const canvas = document.getElementById('view');
const renderer = new Renderer(canvas, scene, view);

let selectedId = null;
let probe = null;             // last probed world point
let aimDir = [1, 0];          // launch aim: in-plane unit vector [u, v] (screen: +u right, +v up)
const particles = [];         // test masses: { x, v, trail:[], color, alive }
let simRunning = false;
const AIM_LEN = 48;           // px — the amber aim pointer length
// Draw a line with a filled arrowhead from (x0,y0) to (x1,y1).
function drawArrow(ctx, x0, y0, x1, y1, color, w, head) {
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = w; ctx.lineCap = 'round';
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
  const a = Math.atan2(y1 - y0, x1 - x0);
  ctx.beginPath(); ctx.moveTo(x1, y1);
  ctx.lineTo(x1 - head * Math.cos(a - 0.5), y1 - head * Math.sin(a - 0.5));
  ctx.lineTo(x1 - head * Math.cos(a + 0.5), y1 - head * Math.sin(a + 0.5));
  ctx.closePath(); ctx.fill(); ctx.lineCap = 'butt';
}

// ---- field unit display ------------------------------------------------
// g in m/s²; 1 Gal = 1e-2 m/s² (the gravimetry unit), so 1 m/s² = 1e5 mGal.
const UNITS = { 'm/s²': 1, mGal: 1e5, 'µGal': 1e8, 'g₀': 1 / P.G0 };
let fieldUnit = 'm/s²';
function fmtField(ms2) {
  const val = ms2 * UNITS[fieldUnit];
  const a = Math.abs(val);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return `${val.toExponential(2)} ${fieldUnit}`;
  const digits = a >= 100 ? 1 : a >= 1 ? 2 : 3;
  return `${val.toFixed(digits)} ${fieldUnit}`;
}
// Auto-scaled SI magnitude, e.g. 12.3 MN, 4.56 µJ/kg — spans planets to grains.
const SI_PREFIX = [
  [1e24, 'Y'], [1e21, 'Z'], [1e18, 'E'], [1e15, 'P'], [1e12, 'T'], [1e9, 'G'],
  [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p'], [1e-15, 'f'],
];
function fmtMag(x, unit) {
  const a = Math.abs(x);
  if (a === 0) return `0 ${unit}`;
  for (const [scale, pre] of SI_PREFIX) {
    if (a >= scale || scale === 1e-15) { const v = x / scale; return `${v.toFixed(Math.abs(v) >= 100 ? 0 : Math.abs(v) >= 10 ? 1 : 2)} ${pre}${unit}`; }
  }
}

// ---- canvas sizing -----------------------------------------------------
let lastW = 0, lastH = 0, lastDpr = 0;
function resize() {
  const wrap = canvas.parentElement;
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = wrap.clientWidth, h = wrap.clientHeight;
  if (w === lastW && h === lastH && dpr === lastDpr) return false;
  lastW = w; lastH = h; lastDpr = dpr;
  view.W = w; view.H = h;
  canvas.width = w * dpr; canvas.height = h * dpr;
  canvas.style.width = w + 'px'; canvas.style.height = h + 'px';
  renderer.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return true;
}
let resizeTimer = null;
window.addEventListener('resize', () => {
  if (!resize()) return;
  requestFrame();
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => { if (renderer.grid) invalidateField(); }, 120);
});

// ---- draw loop ---------------------------------------------------------
let gridDirty = false, layersDirty = false, frameQueued = false;
function tick() {
  frameQueued = false;
  if (simRunning) simStep();                       // advance the N-body system first (moves bodies)
  if (gridDirty) { renderer.computeGrid(simRunning ? 96 : 150, simRunning); renderer.renderField(); gridDirty = layersDirty = false; }
  else if (layersDirty) { renderer.renderField(); layersDirty = false; }
  draw();
  if (simRunning) requestFrame();
}
function requestFrame() { if (!frameQueued) { frameQueued = true; requestAnimationFrame(tick); } }
function draw() {
  renderer.clear();
  renderer.blitField();
  renderer.drawSources(selectedId);
  drawParticles();
  drawProbe();
  drawLegend();
}
function requestDraw() { requestFrame(); }
function invalidateField() { gridDirty = true; updateForceTile(); requestFrame(); }
function invalidateLayers() { layersDirty = true; requestFrame(); }

// ---- probe overlay -----------------------------------------------------
// Keep the probe pin (and its shooter reach) inside the viewport, so panning or
// zooming never loses it — it slides along the edge instead of vanishing.
function clampProbe() {
  if (!probe) return;
  const m = (AIM_LEN + 14) / view.scale;         // world-space margin for pin + shooter
  const cu = view.center[0], cv = view.center[1], hu = view.spanU / 2, hv = view.spanV / 2;
  let u = probe[view.uAxis], v = probe[view.vAxis];
  u = hu > m ? Math.min(cu + hu - m, Math.max(cu - hu + m, u)) : cu;
  v = hv > m ? Math.min(cv + hv - m, Math.max(cv - hv + m, v)) : cv;
  probe = view.worldFromUV(u, v);
}
function drawProbe() {
  if (!probe) return;
  clampProbe();
  const ctx = renderer.ctx;
  const s = view.toScreen(probe);
  const g = scene.g(probe);
  const phi = scene.potential(probe);
  const mag = P.vlen(g);
  // Amber aim pointer with a draggable tip — points the projectile launch and
  // stays where you turn it (drag the tip to aim).
  const tip = [s[0] + aimDir[0] * AIM_LEN, s[1] - aimDir[1] * AIM_LEN];
  drawArrow(ctx, s[0], s[1], tip[0], tip[1], '#e9b44c', 2.5, 7);
  ctx.fillStyle = '#e9b44c'; ctx.beginPath(); ctx.arc(tip[0], tip[1], 4, 0, 7); ctx.fill();
  ctx.strokeStyle = 'rgba(0,0,0,0.5)'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(tip[0], tip[1], 4, 0, 7); ctx.stroke();
  // probe pin: ring + centre dot
  ctx.strokeStyle = '#e9b44c'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(s[0], s[1], 9, 0, 7); ctx.stroke();
  ctx.fillStyle = '#e9b44c'; ctx.beginPath(); ctx.arc(s[0], s[1], 2.5, 0, 7); ctx.fill();
  const u = UNITS[fieldUnit];
  document.getElementById('probeReadout').innerHTML =
    `<div><b>|g|</b> ${fmtField(mag)}</div>` +
    `<div>${(g[0] * u).toExponential(2)}, ${(g[1] * u).toExponential(2)}, ${(g[2] * u).toExponential(2)}</div>` +
    `<div><b>Φ</b> ${isFinite(phi) ? '−' + fmtMag(Math.abs(phi), 'J/kg') : '−∞'}</div>`;
}

// ---- legend ------------------------------------------------------------
function drawLegend() {
  const ctx = renderer.ctx;
  if (!renderer.grid || !renderer.opts.heatmap) return;
  const w = Math.min(150, view.W - 24), h = 9;
  const pad = 8, panelW = w + pad * 2, panelH = h + 30;
  const px = Math.max(6, view.W - panelW - 10);
  const py = Math.max(6, view.H - panelH - 10);
  const x = px + pad, gy = py + 18;
  ctx.fillStyle = 'rgba(7,8,13,0.7)';
  ctx.beginPath();
  if (ctx.roundRect) ctx.roundRect(px, py, panelW, panelH, 6);
  else ctx.rect(px, py, panelW, panelH);
  ctx.fill();
  ctx.fillStyle = '#cdd3dd'; ctx.font = '10px system-ui'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  ctx.fillText('|g| (log)', x, py + 12);
  const grad = ctx.createLinearGradient(x, 0, x + w, 0);
  for (let i = 0; i <= 10; i++) grad.addColorStop(i / 10, viridisCss(i / 10));
  ctx.fillStyle = grad; ctx.fillRect(x, gy, w, h);
  ctx.strokeStyle = 'rgba(255,255,255,0.4)'; ctx.strokeRect(x, gy, w, h);
  ctx.fillStyle = '#cdd3dd'; ctx.textBaseline = 'top';
  ctx.fillText(fmtFieldShort(Math.pow(10, renderer.range.min)), x, gy + h + 2);
  ctx.textAlign = 'right';
  ctx.fillText(fmtFieldShort(Math.pow(10, renderer.range.max)) + ' ' + fieldUnit, x + w, gy + h + 2);
}
function fmtFieldShort(ms2) { const v = ms2 * UNITS[fieldUnit]; const a = Math.abs(v); return a >= 100 || (a < 1e-2 && a > 0) ? v.toExponential(1) : a >= 1 ? v.toFixed(1) : v.toPrecision(2); }
function viridisCss(t) { const c = viridis(t); return `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`; }

// ---- N-body simulation -------------------------------------------------
// Every visible body (placed source + launched projectile) moves under the
// mutual gravity of all the others — exact Newtonian point-mass interaction
// (exact for spheres by the shell theorem). Placed bodies respond too: a rock
// flung past a planet tugs the planet, and a captured body orbits the shared
// barycentre continuously. Positions live on s._origin / p.x; velocities on
// s._v / p.v. Reset (below) restores the placed layout.
let frameTime = 200;            // sim-seconds advanced per animation frame
let baseDt = 4;                 // s per sub-step cap (adaptively shrunk)
let simInited = false;
let fieldSkip = 0;              // throttles the heatmap recompute during playback
let timeRate = 1;              // playback time-rate multiplier (rate slider)
// Show/hide the Play↔Pause label and the rate slider (which only appears while running).
function setRunningUI(running) {
  document.getElementById('pauseSim').textContent = running ? 'Pause' : 'Play';
  document.getElementById('rateRow').style.display = running ? 'flex' : 'none';
}

// Physical radius [m] of a body — sets the contact scale for close encounters
// and the projectile's on-screen size. Rock density ~3000 kg/m³ for projectiles.
function projRadius(mass) { return Math.cbrt(3 * Math.max(mass, 0) / (4 * Math.PI * 3000)); }
function radiusOf(s) {
  switch (s.type) {
    case 'sphere': case 'shell': case 'ring': case 'disc': return s.dia * 500;   // (km/2)→m
    case 'cylinder': return Math.max(s.dia * 500, s.len * 500);
    case 'box': return 0.5 * Math.hypot(s.size[0], s.size[1], s.size[2]) * 1000;
    default: return 0;                                                            // point / rod
  }
}
function drawParticles() {
  const ctx = renderer.ctx;
  for (const p of particles) {
    ctx.strokeStyle = p.color; ctx.lineWidth = 1.6; ctx.beginPath();
    for (let i = 0; i < p.trail.length; i++) {
      const s = view.toScreen(p.trail[i]);
      i ? ctx.lineTo(s[0], s[1]) : ctx.moveTo(s[0], s[1]);
    }
    ctx.stroke();
    const s = view.toScreen(p.x);
    const rr = Math.max(3, (p.rad || 0) * view.scale);   // draw at its true physical size
    ctx.fillStyle = p.alive ? p.color : '#8a8f98';
    ctx.beginPath(); ctx.arc(s[0], s[1], rr, 0, 7); ctx.fill();
    ctx.strokeStyle = p.alive ? 'rgba(0,0,0,0.5)' : '#fff'; ctx.lineWidth = 1; ctx.stroke();
  }
}
// Assemble the live body list; x and v are shared by reference with the
// source / projectile they belong to, so integrating mutates them in place.
function simBodies() {
  const bodies = [];
  for (const s of scene.sources) if (s.visible) bodies.push({ x: s._origin, v: s._v, mass: s.mass, rad: radiusOf(s) });
  for (const p of particles) if (p.alive) bodies.push({ x: p.x, v: p.v, mass: p.mass, rad: p.rad });
  return bodies;
}
function comOf(bodies) {
  let M = 0, c = [0, 0, 0];
  for (const b of bodies) { M += b.mass; c[0] += b.mass * b.x[0]; c[1] += b.mass * b.x[1]; c[2] += b.mass * b.x[2]; }
  return M > 0 ? [c[0] / M, c[1] / M, c[2] / M] : [0, 0, 0];
}
// Shortest resolved timescale across all pairs (encounter time & orbital time).
function adaptiveDt(bodies) {
  let dtm = baseDt;
  for (let i = 0; i < bodies.length; i++) for (let j = i + 1; j < bodies.length; j++) {
    const bi = bodies[i], bj = bodies[j];
    const d = P.vlen(P.vsub(bi.x, bj.x)) || 1;
    const vr = P.vlen(P.vsub(bi.v, bj.v)) || 1;
    const tOrb = Math.sqrt(d * d * d / (P.G * (bi.mass + bj.mass) + 1e-30));
    dtm = Math.min(dtm, 0.02 * d / vr, 0.05 * tOrb);
  }
  return dtm;
}
function simStep() {
  const bodies = simBodies();
  if (!bodies.length) return;
  const trailGap = view.spanU * 0.004;
  const projs = particles.filter((p) => p.alive);
  const srcs = scene.sources.filter((s) => s.visible);
  const ft = frameTime * timeRate;                 // user-adjustable playback pace
  let tacc = 0, sub = 0;
  while (tacc < ft && sub < 4000) {
    let dt = Math.min(adaptiveDt(bodies), ft - tacc);
    if (!(dt > 0)) break;
    P.nbodyStep(bodies, dt);
    tacc += dt; sub++;
    for (const p of projs) {
      if (!p.alive) continue;
      const last = p.trail[p.trail.length - 1];
      if (!last || P.vlen(P.vsub(p.x, last)) > trailGap) p.trail.push(p.x.slice());
      for (const s of srcs) if (P.vlen(P.vsub(p.x, s._origin)) < radiusOf(s) + p.rad) { absorb(s, p); break; }   // impact
    }
  }
  const com = comOf(bodies);
  for (const p of projs) if (P.vlen(P.vsub(p.x, com)) > view.spanU * 12) p.alive = false;   // escaped the view
  for (const p of particles) if (p.trail.length > 5000) p.trail.splice(0, p.trail.length - 5000);
  updateParticleReadout();
  // Bodies moved → refresh the field, but only every other frame: the heatmap is
  // expensive (per-pixel × every source), while the body overlay redraws each
  // frame regardless, so motion stays smooth and the field keeps up closely.
  if ((fieldSkip = (fieldSkip + 1) & 1) === 0) gridDirty = true;
}
// Perfectly inelastic impact: the projectile is absorbed by the body it hits.
// Momentum is conserved — the body recoils to the combined centre-of-mass
// velocity and gains the projectile's mass (kinetic energy is lost to the
// impact, as in a real collision). The projectile stops (is removed).
function absorb(s, p) {
  const M = s.mass, m = p.mass, tot = M + m || 1;
  s._v = [(M * s._v[0] + m * p.v[0]) / tot, (M * s._v[1] + m * p.v[1]) / tot, (M * s._v[2] + m * p.v[2]) / tot];
  s.mass = tot;
  if (s._rings) { const f = tot / (M || 1); for (const r of s._rings) r.m *= f; }   // keep the field in sync without moving _origin
  p.alive = false;
}
// Snapshot the placed layout into the runtime state (velocities from each body's
// vel field; masses so an accreting impact can be undone). Called the first time
// the sim runs after a reset.
function initSim() {
  scene.sources.forEach(buildSource);                       // _origin ← pos
  for (const s of scene.sources) { s._v = (s.vel || [0, 0, 0]).map((c) => c * 1000); s._mass0 = s.mass; }   // km/s → m/s
  calibrateTime();
  simInited = true;
}
function startSim() {
  if (!simInited) initSim();
  setRunningUI(true);
  if (!simRunning) { simRunning = true; requestFrame(); }
}
function resetSim() {
  simRunning = false; simInited = false; particles.length = 0;
  for (const s of scene.sources) { if (s._mass0 != null) s.mass = s._mass0; }   // undo accretion
  scene.sources.forEach(buildSource);                       // restore placed positions (& ring masses)
  for (const s of scene.sources) s._v = [0, 0, 0];
  setRunningUI(false);
  document.getElementById('partReadout').textContent = 'Launch a body, or press Play';
  invalidateField();
}
function updateParticleReadout() {
  const el = document.getElementById('partReadout');
  const p = particles[particles.length - 1];
  if (p) {
    const speed = P.vlen(p.v);
    const eps = 0.5 * speed * speed + scene.potential(p.x);   // specific energy vs the placed sources
    el.innerHTML =
      `<div><b>v</b> ${(speed / 1000).toFixed(2)} km/s</div>` +
      `<div><b>KE</b> ${fmtMag(0.5 * p.mass * speed * speed, 'J')}</div>` +
      `<div><b>ε</b> ${(eps >= 0 ? '+' : '−') + fmtMag(Math.abs(eps), 'J/kg')}</div>` +
      `<div>${p.alive ? (eps < 0 ? 'bound orbit' : 'unbound (escape)') : 'captured / lost'}</div>`;
    return;
  }
  // No projectile: report the selected body's motion during an N-body run.
  const s = scene.get(selectedId);
  if (simInited && s && s._v) {
    const speed = P.vlen(s._v);
    el.innerHTML =
      `<div><b>${s.name}</b></div>` +
      `<div><b>v</b> ${(speed / 1000).toFixed(3)} km/s</div>` +
      `<div><b>KE</b> ${fmtMag(0.5 * s.mass * speed * speed, 'J')}</div>`;
  }
}

// ---------------------------------------------------------------------------
// UI construction
// ---------------------------------------------------------------------------
// 'mass' → mass field with the named-body preset dropdown (only for whole
// celestial bodies: point / sphere / shell). 'massonly' → a plain kg field for
// extended shapes (a rod or ring is not "Jupiter").
const paramDefs = {
  point:    [['mass']],
  sphere:   [['mass'], ['Dia (km)', 'dia', 100, 1500000, 500]],
  shell:    [['mass'], ['Dia (km)', 'dia', 100, 1500000, 500]],
  ring:     [['massonly'], ['Dia (km)', 'dia', 1000, 300000, 500]],
  disc:     [['massonly'], ['Dia (km)', 'dia', 1000, 300000, 500], ['Thickness (km)', 'thick', 0, 120000, 250]],
  cylinder: [['massonly'], ['Dia (km)', 'dia', 1000, 200000, 500], ['Len (km)', 'len', 1000, 300000, 500]],
  rod:      [['massonly'], ['Len (km)', 'len', 1000, 400000, 500]],
  box:      [['massonly'], ['W (km)', 'size.0', 1000, 200000, 500], ['H (km)', 'size.1', 1000, 200000, 500], ['L (km)', 'size.2', 1000, 200000, 500]],
};
const commonDefs = [
  ['X (km)', 'pos.0', -250000, 250000, 500],
  ['Y (km)', 'pos.1', -250000, 250000, 500],
  ['Z (km)', 'pos.2', -250000, 250000, 500],
  ['Yaw °', 'rot.0', -180, 180, 1],
  ['Pitch °', 'rot.1', -180, 180, 1],
  ['Roll °', 'rot.2', -180, 180, 1],
  ['Vx (km/s)', 'vel.0', -60, 60, 0.05],
  ['Vy (km/s)', 'vel.1', -60, 60, 0.05],
  ['Vz (km/s)', 'vel.2', -60, 60, 0.05],
];

function getPath(o, path) { const k = path.split('.'); let v = o; for (const p of k) v = v[isNaN(p) ? p : +p]; return v; }
function setPath(o, path, val) { const k = path.split('.'); let v = o; for (let i = 0; i < k.length - 1; i++) v = v[isNaN(k[i]) ? k[i] : +k[i]]; const last = k[k.length - 1]; v[isNaN(last) ? last : +last] = val; }

function buildInspector() {
  const el = document.getElementById('inspector');
  el.innerHTML = '';
  const s = scene.get(selectedId);
  if (!s) { el.innerHTML = '<p class="hint">Select an object</p>'; updateForceTile(); return; }

  const addRow = (label, path, min, max, step) => {
    const row = document.createElement('label'); row.className = 'row';
    const val = getPath(s, path);
    const m = label.match(/^(.*?)\s*\((.+)\)$/);
    row.innerHTML = m ? `<span>${m[1]} <span class="lc">(${m[2]})</span></span>` : `<span>${label}</span>`;
    const input = document.createElement('input');
    input.type = 'number'; input.value = val; input.step = step; input.min = min; input.max = max;
    input.addEventListener('input', () => {
      let n = parseFloat(input.value); if (isNaN(n)) return;
      setPath(s, path, n); buildSource(s); invalidateField();
    });
    if (min !== undefined) {
      const rng = document.createElement('input');
      rng.type = 'range'; rng.min = min; rng.max = max; rng.step = step; rng.value = val;
      rng.addEventListener('input', () => { setPath(s, path, parseFloat(rng.value)); input.value = rng.value; buildSource(s); invalidateField(); });
      input.addEventListener('input', () => { rng.value = input.value; });
      row.appendChild(input); row.appendChild(rng);
    } else { row.appendChild(input); }
    el.appendChild(row);
  };

  for (const def of (paramDefs[s.type] || [])) {
    if (def[0] === 'mass') { addMassRow(el, s, true); continue; }
    if (def[0] === 'massonly') { addMassRow(el, s, false); continue; }
    addRow(...def);
  }
  const hr = document.createElement('hr'); el.appendChild(hr);
  for (const def of commonDefs) addRow(...def);
  updateForceTile();
}

// Mass: a named-body dropdown that seeds the value (and renames the object),
// plus a free kg field. The dropdown selection persists on the source (s.body).
const cleanBodyName = (key) => key.replace(/\s*\(.*\)$/, '');
function addMassRow(el, s, withPreset) {
  let sel = null;
  const massInput = document.createElement('input');
  if (withPreset) {
    const row1 = document.createElement('label'); row1.className = 'row';
    row1.innerHTML = '<span>Body</span>';
    sel = document.createElement('select'); sel.id = 'bodySel';
    sel.innerHTML = '<option value="">— preset —</option>';
    for (const name of Object.keys(BODIES)) sel.innerHTML += `<option value="${name}">${name}</option>`;
    sel.value = s.body || '';                       // persist the chosen preset
    sel.addEventListener('change', () => {
      if (!sel.value) { s.body = ''; return; }
      s.body = sel.value; s.mass = BODIES[sel.value]; s.name = cleanBodyName(sel.value);
      // sphere/shell are actual celestial bodies — give them the real diameter too
      if ((s.type === 'sphere' || s.type === 'shell') && BODY_DIA[sel.value]) s.dia = BODY_DIA[sel.value];
      buildSource(s); buildList(); buildInspector(); fitView();   // refresh Dia field + frame it
    });
    row1.appendChild(sel); el.appendChild(row1);
  }

  const row2 = document.createElement('label'); row2.className = 'row';
  row2.innerHTML = '<span>Mass <span class="lc">(kg)</span></span>';
  massInput.type = 'number'; massInput.value = s.mass; massInput.step = '1e23'; massInput.min = 0;
  massInput.addEventListener('input', () => {
    const n = parseFloat(massInput.value); if (isNaN(n) || n < 0) return;
    s.mass = n; s.body = ''; if (sel) sel.value = '';   // no longer a named preset
    buildSource(s); invalidateField();
  });
  row2.appendChild(massInput); el.appendChild(row2);
}

// Force/torque data tile (right panel), kept current on any scene change.
function dirArrow(vec) {
  const c = view.planeComps(vec), ip = Math.hypot(c.u, c.v);
  if (ip === 0 && c.n === 0) return '';
  if (Math.abs(c.n) > ip) return c.n > 0 ? '⊙' : '⊗';
  const arr = ['→', '↗', '↑', '↖', '←', '↙', '↓', '↘'];
  return arr[((Math.round(Math.atan2(c.v, c.u) / (Math.PI / 4)) % 8) + 8) % 8];
}
function updateForceTile() {
  const el = document.getElementById('forceReadout');
  const s = scene.get(selectedId);
  if (!s) { el.textContent = 'Select an object'; return; }
  const ft = scene.forceTorque(s);
  if (!ft || !ft.hasExternal) { el.innerHTML = '<div>Net force needs</div><div class="hint">a second mass</div>'; return; }
  if (!ft.valid) { el.innerHTML = '<div>Bodies overlap</div><div class="hint">separate to read force</div>'; return; }
  const Fm = P.vlen(ft.F), tauM = P.vlen(ft.tau);
  const fZero = ft.Fabs > 0 && Fm / ft.Fabs < 1e-7;
  const tZero = ft.tauAbs > 0 && tauM / ft.tauAbs < 1e-7;
  el.innerHTML =
    `<div><b>F</b> ${fZero ? '≈ 0 N' : fmtMag(Fm, 'N') + ' ' + dirArrow(ft.F)}</div>` +
    `<div><b>τ</b> ${tZero ? '≈ 0 N·m' : fmtMag(tauM, 'N·m')}</div>`;
}

function buildList() {
  const el = document.getElementById('objlist'); el.innerHTML = '';
  for (const s of scene.sources) {
    const row = document.createElement('div');
    row.className = 'obj' + (s.id === selectedId ? ' sel' : '');
    row.innerHTML = `<span class="dot" style="background:${s.color}"></span><span class="onm">${s.name}</span>` +
      `<button class="vis" title="Show/hide">${s.visible ? '👁' : '∅'}</button>` +
      `<button class="del" title="Delete">✕</button>`;
    row.querySelector('.onm').addEventListener('click', () => { selectedId = s.id; buildList(); buildInspector(); requestDraw(); });
    row.querySelector('.vis').addEventListener('click', (e) => { e.stopPropagation(); s.visible = !s.visible; buildList(); invalidateField(); });
    row.querySelector('.del').addEventListener('click', (e) => {
      e.stopPropagation(); scene.remove(s.id);
      if (selectedId === s.id) selectedId = null;
      buildList(); buildInspector(); invalidateField();
    });
    el.appendChild(row);
  }
}

// ---- add-object buttons ------------------------------------------------
document.querySelectorAll('[data-add]').forEach((btn) => {
  btn.addEventListener('click', () => {
    const s = defaultSource(btn.dataset.add);
    scene.add(s); selectedId = s.id; buildList(); buildInspector(); invalidateField();
  });
});

// ---- layer toggles -----------------------------------------------------
function bindToggle(id, key) {
  const el = document.getElementById(id);
  el.checked = renderer.opts[key];
  el.addEventListener('change', () => { renderer.opts[key] = el.checked; invalidateLayers(); });
}
bindToggle('tglHeat', 'heatmap'); bindToggle('tglLines', 'lines');
bindToggle('tglEqui', 'equipotential'); bindToggle('tglVec', 'vectors'); bindToggle('tglGrid', 'grid');

// ---- view controls -----------------------------------------------------
const planes = {
  'XZ · side': [0, 2, 1],
  'XY · top':  [0, 1, 2],
  'YZ · front':[1, 2, 0],
};
const planeSel = document.getElementById('planeSel');
for (const name of Object.keys(planes)) planeSel.innerHTML += `<option value="${name}">${name.toUpperCase()}</option>`;
planeSel.addEventListener('change', () => {
  const [u, v, n] = planes[planeSel.value]; view.uAxis = u; view.vAxis = v; view.nAxis = n;
  document.getElementById('axU').textContent = view.axisLabel(u);
  document.getElementById('axV').textContent = view.axisLabel(v);
  invalidateField();
});
const sliceInput = document.getElementById('sliceInput');
sliceInput.addEventListener('input', () => {
  const v = parseFloat(sliceInput.value);
  if (!isFinite(v)) return;
  view.slice = v * 1000; invalidateField();       // km -> m
});
const unitSel = document.getElementById('unitSel');
for (const u of Object.keys(UNITS)) unitSel.innerHTML += `<option>${u}</option>`;
unitSel.value = fieldUnit;
unitSel.addEventListener('change', () => { fieldUnit = unitSel.value; requestDraw(); });

// ---- snap, pan & zoom --------------------------------------------------
let snap = false, snapStep = 1000;               // km
const maybeSnap = (v) => snap ? Math.round(v / snapStep) * snapStep : v;
const snapChk = document.getElementById('snapChk');
snapChk.addEventListener('change', () => { snap = snapChk.checked; });
document.getElementById('snapStep').addEventListener('change', (e) => { snapStep = Math.max(100, parseFloat(e.target.value) || 1000); });

const MAX_SPAN = 2e13;   // widest zoom-out [m] — fits the outer solar system
function zoomBy(factor, sx, sy) {
  const before = (sx !== undefined) ? view.toWorld(sx, sy) : null;
  view.spanU = Math.max(2e6, Math.min(MAX_SPAN, view.spanU * factor));
  if (before) {
    const after = view.toWorld(sx, sy);
    view.center[0] += before[view.uAxis] - after[view.uAxis];
    view.center[1] += before[view.vAxis] - after[view.vAxis];
  }
  invalidateField();
}
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  const r = canvas.getBoundingClientRect();
  zoomBy(Math.exp(e.deltaY * 0.0012), e.clientX - r.left, e.clientY - r.top);
}, { passive: false });

document.getElementById('zoomIn').addEventListener('click', () => zoomBy(1 / 1.3, view.W / 2, view.H / 2));
document.getElementById('zoomOut').addEventListener('click', () => zoomBy(1.3, view.W / 2, view.H / 2));
document.getElementById('zoomReset').addEventListener('click', () => { view.spanU = 8e7; view.center = [0, 0]; invalidateField(); });
document.getElementById('zoomFit').addEventListener('click', fitView);
function fitView() {
  const vis = scene.sources.filter((s) => s.visible);
  if (!vis.length) return;
  let uMin = Infinity, uMax = -Infinity, vMin = Infinity, vMax = -Infinity;
  for (const s of vis) {
    const e = sourceExtent(s), o = s._origin;
    uMin = Math.min(uMin, o[view.uAxis] - e); uMax = Math.max(uMax, o[view.uAxis] + e);
    vMin = Math.min(vMin, o[view.vAxis] - e); vMax = Math.max(vMax, o[view.vAxis] + e);
  }
  view.center = [(uMin + uMax) / 2, (vMin + vMax) / 2];
  const span = Math.max(uMax - uMin, (vMax - vMin) * view.W / view.H, 4e6) * 1.6;
  view.spanU = Math.max(4e6, Math.min(MAX_SPAN, span));
  invalidateField();
}

let dragMode = null, dragStart = null, dragObjStart = null;
const pointers = new Map();
let pinch = null;
const localXY = (e) => { const r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
const nearProbe = (sx, sy) => { if (!probe) return false; const p = view.toScreen(probe); return Math.hypot(p[0] - sx, p[1] - sy) < 20; };
const aimTipScreen = () => { if (!probe) return null; const p = view.toScreen(probe); return [p[0] + aimDir[0] * AIM_LEN, p[1] - aimDir[1] * AIM_LEN]; };
const nearAimTip = (sx, sy) => { const t = aimTipScreen(); return t && Math.hypot(t[0] - sx, t[1] - sy) < 16; };
// Point the aim from the probe pin toward a screen position (unit vector, +v up).
function setAimTo(sx, sy) {
  const p = view.toScreen(probe);
  const du = sx - p[0], dv = -(sy - p[1]), L = Math.hypot(du, dv) || 1;
  aimDir = [du / L, dv / L];
}

canvas.addEventListener('pointerdown', (e) => {
  canvas.setPointerCapture(e.pointerId);
  const [sx, sy] = localXY(e); pointers.set(e.pointerId, [sx, sy]);
  if (pointers.size === 2) {
    const p = [...pointers.values()];
    pinch = { dist: Math.hypot(p[0][0] - p[1][0], p[0][1] - p[1][1]) || 1, span: view.spanU };
    dragMode = null; return;
  }
  if (nearAimTip(sx, sy)) { dragMode = 'aim'; setAimTo(sx, sy); canvas.style.cursor = 'grabbing'; requestDraw(); return; }
  if (nearProbe(sx, sy)) { dragMode = 'probe'; canvas.style.cursor = 'grabbing'; return; }
  const hit = pickSource(sx, sy);
  if (hit) {
    if (hit.id !== selectedId) { selectedId = hit.id; buildList(); buildInspector(); requestDraw(); }
    dragMode = 'obj'; dragStart = [sx, sy]; dragObjStart = hit.pos.slice();
    canvas.style.cursor = 'grabbing';
  } else {
    if (selectedId !== null) { selectedId = null; buildList(); buildInspector(); requestDraw(); }
    dragMode = 'pan'; dragStart = [sx, sy, view.center[0], view.center[1]];
    canvas.style.cursor = 'grabbing';
  }
});
canvas.addEventListener('pointermove', (e) => {
  const [sx, sy] = localXY(e);
  if (pointers.has(e.pointerId)) pointers.set(e.pointerId, [sx, sy]);
  if (pinch && pointers.size >= 2) {
    const p = [...pointers.values()];
    const dist = Math.hypot(p[0][0] - p[1][0], p[0][1] - p[1][1]) || 1;
    const mx = (p[0][0] + p[1][0]) / 2, my = (p[0][1] + p[1][1]) / 2;
    const before = view.toWorld(mx, my);
    view.spanU = Math.max(2e6, Math.min(MAX_SPAN, pinch.span * pinch.dist / dist));
    const after = view.toWorld(mx, my);
    view.center[0] += before[view.uAxis] - after[view.uAxis];
    view.center[1] += before[view.vAxis] - after[view.vAxis];
    invalidateField(); return;
  }
  if (dragMode === 'aim') {
    setAimTo(sx, sy); requestDraw();
  } else if (dragMode === 'probe') {
    probe = view.toWorld(sx, sy); requestDraw();
  } else if (dragMode === 'pan') {
    view.center[0] = dragStart[2] - (sx - dragStart[0]) / view.scale;
    view.center[1] = dragStart[3] + (sy - dragStart[1]) / view.scale;
    invalidateField();
  } else if (dragMode === 'obj') {
    const s = scene.get(selectedId); if (!s) return;
    const du = (sx - dragStart[0]) / view.scale / 1000;     // km
    const dv = -(sy - dragStart[1]) / view.scale / 1000;
    s.pos[view.uAxis] = maybeSnap(dragObjStart[view.uAxis] + du);
    s.pos[view.vAxis] = maybeSnap(dragObjStart[view.vAxis] + dv);
    buildSource(s); invalidateField();
  } else if (e.pointerType === 'mouse') {
    canvas.style.cursor = nearAimTip(sx, sy) || nearProbe(sx, sy) || pickSource(sx, sy) ? 'grab' : 'crosshair';
  }
});
function endPointer(e) {
  pointers.delete(e.pointerId);
  if (pointers.size < 2) pinch = null;
  if (pointers.size === 0) {
    if (dragMode === 'obj') buildInspector();
    dragMode = null; canvas.style.cursor = 'crosshair';
  }
}
canvas.addEventListener('pointerup', endPointer);
canvas.addEventListener('pointercancel', endPointer);

// keyboard: arrow keys nudge the selected object, Delete/Backspace removes it
window.addEventListener('keydown', (e) => {
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA')) return;
  const s = scene.get(selectedId); if (!s) return;
  if (e.key === 'Delete' || e.key === 'Backspace') {
    scene.remove(s.id); selectedId = null; buildList(); buildInspector(); invalidateField();
    e.preventDefault(); return;
  }
  const step = (snap ? snapStep : 1000) * (e.shiftKey ? 5 : 1);   // km
  const move = { ArrowLeft: [view.uAxis, -step], ArrowRight: [view.uAxis, step],
                 ArrowUp: [view.vAxis, step], ArrowDown: [view.vAxis, -step] }[e.key];
  if (!move) return;
  s.pos[move[0]] += move[1];
  buildSource(s); buildInspector(); invalidateField();
  e.preventDefault();
});

document.getElementById('clearAll').addEventListener('click', () => {
  scene.sources = []; selectedId = null; particles.length = 0; simRunning = false; simInited = false;
  setRunningUI(false);
  document.getElementById('partReadout').textContent = 'Launch a body, or press Play';
  buildList(); buildInspector(); invalidateField();
});

// Hit-test in screen space: front-most (last-added) object under the cursor.
function pickSource(sx, sy) {
  for (let i = scene.sources.length - 1; i >= 0; i--) {
    const s = scene.sources[i]; if (!s.visible) continue;
    const p = view.toScreen(s._origin);
    const r = Math.max(14, sourceExtent(s) * view.scale);
    if (Math.hypot(p[0] - sx, p[1] - sy) < r) return s;
  }
  return null;
}

// ---- test-mass launcher ------------------------------------------------
// Set the animation cadence from the scene's characteristic free-fall time so
// an orbit takes a few hundred frames regardless of the masses/scale in play.
function calibrateTime() {
  let Mtot = 0; for (const s of scene.sources) if (s.visible) Mtot += massOf(s);
  const r = view.spanU / 3;
  const tau = Mtot > 0 ? 2 * Math.PI * Math.sqrt(r ** 3 / (P.G * Mtot)) : 1e4;
  frameTime = tau / 260;
  baseDt = frameTime / 30;
}
function launchParticle() {
  if (!simInited) initSim();
  const speed = (Math.abs(parseFloat(document.getElementById('pSpeed').value)) || 5) * 1000;   // km/s -> m/s
  const mass = Math.max(0, parseFloat(document.getElementById('pMass').value)) || 1;            // kg
  // Launch from the field-probe pin, along the red shooter aim (which stays
  // where the user turned it). Falls back to the view edge when unset.
  const pos = probe ? probe.slice() : view.worldFromUV(view.center[0] - view.spanU * 0.42, view.center[1]);
  const dir = [0, 0, 0]; dir[view.uAxis] = aimDir[0]; dir[view.vAxis] = aimDir[1];
  const vel = dir.map((c) => c * speed);
  particles.push({ x: pos, v: vel, mass, rad: projRadius(mass), trail: [pos.slice()], color: '#ffd27a', alive: true });
  startSim();
}
document.getElementById('launch').addEventListener('click', launchParticle);
document.getElementById('clearParts').addEventListener('click', resetSim);
// Launch-speed bar: keep the slider and the number box in sync.
{
  const rng = document.getElementById('pSpeedRange'), num = document.getElementById('pSpeed');
  rng.addEventListener('input', () => { num.value = rng.value; });
  num.addEventListener('input', () => { const n = parseFloat(num.value); if (!isNaN(n)) rng.value = n; });
}
document.getElementById('pauseSim').addEventListener('click', () => {
  if (simRunning) { simRunning = false; setRunningUI(false); }
  else { startSim(); }
});
// Playback time-rate slider (log scale, 0.05×…20×, 1× centred). Shown only
// while running (setRunningUI); the chosen rate persists across plays.
const rateSlider = document.getElementById('rateSlider');
rateSlider.addEventListener('input', () => {
  timeRate = 0.05 * Math.pow(400, parseFloat(rateSlider.value) / 100);
  const txt = timeRate < 1 ? timeRate.toFixed(2) : timeRate < 10 ? timeRate.toFixed(1) : timeRate.toFixed(0);
  document.getElementById('rateVal').textContent = txt + '×';
});

// ---- presets / scenarios ----------------------------------------------
const presets = {
  'Solar system': () => {
    // The Sun + 8 planets on circular, coplanar orbits (real masses, distances
    // and orbital speeds). Press Play — inner planets whirl, outer ones crawl,
    // and the Sun wobbles about the barycentre (mostly Jupiter's doing).
    scene.sources = [];
    const GM = P.G * 1.989e30;                                  // Sun's GM [m³/s²]
    const sun = defaultSource('sphere'); sun.name = 'Sun'; sun.mass = 1.989e30; sun.dia = 1391000; sun.pos = [0, 0, 0]; sun.color = '#ffcf3a';
    // name, mass [kg], orbital radius [km], diameter [km], colour
    const planets = [
      ['Mercury', 3.301e23, 57.9e6, 4879, '#b0a58e'],
      ['Venus', 4.867e24, 108.2e6, 12104, '#d8b47a'],
      ['Earth', 5.972e24, 149.6e6, 12742, '#4aa3ff'],
      ['Mars', 6.417e23, 227.9e6, 6779, '#e0682f'],
      ['Jupiter', 1.898e27, 778.5e6, 139820, '#d9a066'],
      ['Saturn', 5.683e26, 1434e6, 116460, '#e6d3a3'],
      ['Uranus', 8.681e25, 2871e6, 50724, '#7fd8e0'],
      ['Neptune', 1.024e26, 4495e6, 49244, '#3f66e0'],
    ];
    let pz = 0;                                                 // total planet z-momentum [kg·m/s]
    const bodies = [];
    for (const [name, mass, aKm, dia, color] of planets) {
      const v = Math.sqrt(GM / (aKm * 1000));                   // circular-orbit speed [m/s]
      const s = defaultSource('sphere');
      s.name = name; s.mass = mass; s.dia = dia; s.pos = [aKm, 0, 0]; s.vel = [0, 0, v / 1000]; s.color = color;
      pz += mass * v; bodies.push(s);
    }
    sun.vel = [0, 0, -(pz / 1.989e30) / 1000];                  // km/s — balance total momentum
    scene.add(sun); for (const s of bodies) scene.add(s);
    view.spanU = 1.1e13; view.center = [0, 0];
    return sun.id;
  },
  'Binary system': () => {
    // Two stars on circular orbits about their common barycentre. Press Play.
    scene.sources = [];
    const a = defaultSource('sphere'); a.name = 'Star A'; a.mass = 4e24; a.dia = 12000; a.pos = [-30000, 0, 0];
    const b = defaultSource('sphere'); b.name = 'Star B'; b.mass = 2e24; b.dia = 9000; b.pos = [30000, 0, 0];
    a.vel = [0, 0, 0.861]; b.vel = [0, 0, -1.721];       // km/s: m₁v₁ = m₂v₂, opposite, ⟂ separation
    scene.add(a); scene.add(b); view.spanU = 1.8e8;
  },
  'Planet + moon (orbit)': () => {
    // Moon on a near-circular orbit; the planet recoils about the barycentre. Play.
    scene.sources = [];
    const p = defaultSource('sphere'); p.name = 'Planet'; p.mass = 6e24; p.dia = 13000; p.pos = [0, 0, 0];
    const m = defaultSource('sphere'); m.name = 'Moon'; m.mass = 7e22; m.dia = 4000; m.pos = [90000, 0, 0];
    m.vel = [0, 0, 2.12]; p.vel = [0, 0, -2.12 * 7e22 / 6e24];   // km/s, momentum-balanced
    scene.add(p); scene.add(m); view.spanU = 3e8;
    document.getElementById('pSpeed').value = 2.1;
  },
  'Spherical shell (no field inside)': () => {
    scene.sources = [];
    const s = defaultSource('shell'); s.name = 'Shell'; s.mass = 6e24; s.dia = 60000; s.pos = [0, 0, 0];
    scene.add(s); view.spanU = 1.4e8;
    renderer.opts.equipotential = true; document.getElementById('tglEqui').checked = true;
  },
  'Ring world': () => {
    scene.sources = [];
    const r = defaultSource('ring'); r.name = 'Ring'; r.mass = 5e24; r.dia = 80000; r.pos = [0, 0, 0];
    scene.add(r); view.spanU = 1.6e8;
    planeSel.value = 'XZ · side'; planeSel.dispatchEvent(new Event('change'));
  },
  'Tidal rod near a planet': () => {
    scene.sources = [];
    const p = defaultSource('point'); p.name = 'Planet'; p.mass = 6e24; p.pos = [0, 0, 0];
    const rod = defaultSource('rod'); rod.name = 'Rod'; rod.len = 40000; rod.mass = 1e21; rod.pos = [90000, 0, 0]; rod.rot = [0, 45, 0];
    scene.add(p); scene.add(rod); view.spanU = 3e8;
    return rod.id;
  },
  'Slingshot fly-by': () => {
    scene.sources = [];
    const p = defaultSource('sphere'); p.name = 'Planet'; p.mass = 8e24; p.dia = 12000; p.pos = [10000, -18000, 0];
    scene.add(p); view.spanU = 1.6e8;
    document.getElementById('pSpeed').value = 6;
  },
};
const presetSel = document.getElementById('presetSel');
presetSel.innerHTML = '<option value="">SCENARIO…</option>';
for (const name of Object.keys(presets)) presetSel.innerHTML += `<option value="${name}">${name.toUpperCase()}</option>`;
presetSel.addEventListener('change', () => {
  if (!presets[presetSel.value]) return;
  particles.length = 0; simRunning = false; simInited = false;
  setRunningUI(false);
  document.getElementById('partReadout').textContent = 'Launch a body, or press Play';
  const wantSel = presets[presetSel.value]();
  scene.rebuild();
  selectedId = wantSel || (scene.sources[0] ? scene.sources[0].id : null);
  presetSel.value = ''; view.center = [0, 0]; buildList(); buildInspector(); invalidateField();
});

// ---- init --------------------------------------------------------------
function init() {
  resize();
  const s = defaultSource('sphere'); s.name = 'Planet'; s.mass = 6e24; s.dia = 13000;
  scene.add(s); view.spanU = 8e7;
  selectedId = s.id;
  snap = true; snapStep = 1000;
  snapChk.checked = true; document.getElementById('snapStep').value = 1000;
  document.getElementById('axU').textContent = view.axisLabel(view.uAxis);
  document.getElementById('axV').textContent = view.axisLabel(view.vAxis);
  probe = view.worldFromUV(view.center[0] + view.spanU * 0.28, view.center[1]);
  buildList(); buildInspector();
  invalidateField();
}
init();
