// render.js — 2D canvas visualisation of a slice through the 3D scene.
// The field layers (heatmap, field lines, equipotentials, quiver, grid) are
// expensive, so they are rendered once to an offscreen canvas whenever the
// scene or view changes; per-frame drawing (during the orbit animation) just
// blits that layer and adds lightweight overlays (bodies, test masses, probe).
import * as P from './physics.js';
import { sourceExtent } from './sources.js';

// Compact viridis colour ramp (control points), interpolated in sRGB.
const VIRIDIS = [
  [68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142],
  [31, 158, 137], [53, 183, 121], [110, 206, 88], [181, 222, 43], [253, 231, 37],
];
function viridis(t) {
  t = Math.max(0, Math.min(1, t));
  const x = t * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x));
  const f = x - i, a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

export class View {
  constructor() {
    this.uAxis = 0; this.vAxis = 2; this.nAxis = 1; // XZ plane, slice along Y
    this.slice = 0;
    this.center = [0, 0];
    this.spanU = 8e7;              // metres across the view (default ~80,000 km)
    this.W = 800; this.H = 600;
  }
  get scale() { return this.W / this.spanU; }
  toScreen(w) {
    return [this.W / 2 + (w[this.uAxis] - this.center[0]) * this.scale,
            this.H / 2 - (w[this.vAxis] - this.center[1]) * this.scale];
  }
  toWorld(sx, sy) {
    const u = this.center[0] + (sx - this.W / 2) / this.scale;
    const v = this.center[1] - (sy - this.H / 2) / this.scale;
    return this.worldFromUV(u, v);
  }
  worldFromUV(u, v) {
    const w = [0, 0, 0];
    w[this.uAxis] = u; w[this.vAxis] = v; w[this.nAxis] = this.slice;
    return w;
  }
  planeComps(v) { return { u: v[this.uAxis], v: v[this.vAxis], n: v[this.nAxis] }; }
  axisLabel(i) { return ['X', 'Y', 'Z'][i]; }
  get spanV() { return this.spanU * this.H / this.W; }
}

export class Renderer {
  constructor(canvas, scene, view) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.scene = scene;
    this.view = view;
    this.grid = null;
    this.heat = null;               // offscreen grid-resolution heatmap
    this.field = document.createElement('canvas'); // offscreen full-res field layer
    this.opts = { heatmap: true, lines: false, equipotential: false, vectors: true, grid: true, gridStep: 34 };
    this.range = { min: -6, max: 0 };
  }

  // Sample |g|, in-plane g-components and the potential Φ over a grid.
  // `smooth` eases the colour range toward the new extrema instead of snapping,
  // so during playback the heatmap doesn't flash as the field extrema wobble.
  computeGrid(cols = 150, smooth = false) {
    const v = this.view;
    const rows = Math.max(2, Math.round(cols * v.H / v.W));
    const mag = new Float32Array(cols * rows);
    const gu = new Float32Array(cols * rows);
    const gv = new Float32Array(cols * rows);
    const phi = new Float32Array(cols * rows);
    let lo = Infinity, hi = -Infinity;
    for (let j = 0; j < rows; j++) {
      for (let i = 0; i < cols; i++) {
        const sx = (i + 0.5) / cols * v.W;
        const sy = (j + 0.5) / rows * v.H;
        const w = v.toWorld(sx, sy);
        const g = this.scene.g(w);
        const m = P.vlen(g);
        const idx = j * cols + i;
        mag[idx] = m; gu[idx] = g[v.uAxis]; gv[idx] = g[v.vAxis];
        phi[idx] = this.scene.potential(w);
        if (m > 0) { const l = Math.log10(m); if (l < lo) lo = l; if (l > hi) hi = l; }
      }
    }
    if (!isFinite(lo)) { lo = -6; hi = 0; }
    lo = Math.max(lo, hi - 6);                    // clamp to ~6 decades
    this.grid = { cols, rows, mag, gu, gv, phi };
    // Ease the colour range during playback (smooth), snap otherwise.
    if (smooth && this._rangeSet) {
      const a = 0.08;
      this.range.min += (lo - this.range.min) * a;
      this.range.max += (hi - this.range.max) * a;
    } else {
      this.range = { min: lo, max: hi };
      this._rangeSet = true;
    }
    const rlo = this.range.min, rspan = (this.range.max - this.range.min) || 1;
    // build heatmap image at grid resolution
    const off = this.heat || (this.heat = document.createElement('canvas'));
    off.width = cols; off.height = rows;
    const octx = off.getContext('2d');
    const img = octx.createImageData(cols, rows);
    for (let k = 0; k < cols * rows; k++) {
      const m = mag[k];
      const t = m > 0 ? (Math.log10(m) - rlo) / rspan : 0;
      const c = viridis(t);
      img.data[k * 4] = c[0]; img.data[k * 4 + 1] = c[1]; img.data[k * 4 + 2] = c[2]; img.data[k * 4 + 3] = 255;
    }
    octx.putImageData(img, 0, 0);
  }

  // Bilinear in-plane field from the cached grid, in world (u,v) coordinates.
  sampleField(u, v) {
    const g = this.grid, view = this.view;
    const fx = (view.W / 2 + (u - view.center[0]) * view.scale) / view.W * g.cols - 0.5;
    const fy = (view.H / 2 - (v - view.center[1]) * view.scale) / view.H * g.rows - 0.5;
    const x0 = Math.floor(fx), y0 = Math.floor(fy);
    if (x0 < 0 || y0 < 0 || x0 >= g.cols - 1 || y0 >= g.rows - 1) return null;
    const tx = fx - x0, ty = fy - y0;
    const at = (xx, yy, arr) => arr[yy * g.cols + xx];
    const lerp2 = (arr) =>
      (at(x0, y0, arr) * (1 - tx) + at(x0 + 1, y0, arr) * tx) * (1 - ty) +
      (at(x0, y0 + 1, arr) * (1 - tx) + at(x0 + 1, y0 + 1, arr) * tx) * ty;
    return [lerp2(g.gu), lerp2(g.gv)];
  }

  // Render all field layers to the offscreen field canvas.
  renderField() {
    const v = this.view;
    const f = this.field;
    f.width = v.W; f.height = v.H;
    const ctx = f.getContext('2d');
    ctx.clearRect(0, 0, v.W, v.H);
    if (this.opts.heatmap && this.heat) {
      ctx.imageSmoothingEnabled = true; ctx.globalAlpha = 0.95;
      ctx.drawImage(this.heat, 0, 0, v.W, v.H);
      ctx.globalAlpha = 1;
    }
    if (this.opts.grid) this.drawGrid(ctx);
    if (this.opts.equipotential) this.drawEquipotential(ctx);
    if (this.opts.lines) this.drawFieldLines(ctx);
    if (this.opts.vectors) this.drawVectors(ctx);
  }

  // Trace an in-plane field line from a seed. Gravity is a sink field, so the
  // forward direction (dir = +1) flows INTO the masses.
  fieldLine(u0, v0, dir) {
    const v = this.view, pts = [];
    let u = u0, w = v0;
    const ds = v.spanU / 300;
    const f = (uu, vv) => {
      const s = this.sampleField(uu, vv);
      if (!s) return null;
      const L = Math.hypot(s[0], s[1]);
      if (L < 1e-16) return null;
      return [dir * s[0] / L, dir * s[1] / L];
    };
    for (let step = 0; step < 320; step++) {
      pts.push([u, w]);
      const k1 = f(u, w); if (!k1) break;
      const k2 = f(u + k1[0] * ds / 2, w + k1[1] * ds / 2); if (!k2) break;
      const k3 = f(u + k2[0] * ds / 2, w + k2[1] * ds / 2); if (!k3) break;
      const k4 = f(u + k3[0] * ds, w + k3[1] * ds); if (!k4) break;
      u += ds * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
      w += ds * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
    }
    return pts;
  }

  drawFieldLines(ctx) {
    const v = this.view;
    ctx.lineWidth = 1.1;
    ctx.strokeStyle = 'rgba(226,232,255,0.5)';
    const seeds = 16;
    for (let a = 0; a < seeds; a++) {
      for (let b = 0; b < seeds; b++) {
        const u0 = v.center[0] - v.spanU / 2 + v.spanU * (a + 0.5) / seeds;
        const v0 = v.center[1] - v.spanV / 2 + v.spanV * (b + 0.5) / seeds;
        for (const dir of [1, -1]) {
          const line = this.fieldLine(u0, v0, dir);
          if (line.length < 6) continue;
          ctx.beginPath();
          for (let i = 0; i < line.length; i++) {
            const s = v.toScreen(v.worldFromUV(line[i][0], line[i][1]));
            if (i === 0) ctx.moveTo(s[0], s[1]); else ctx.lineTo(s[0], s[1]);
          }
          ctx.stroke();
          if (dir === 1 && line.length > 12) this.arrowAt(ctx, line);
        }
      }
    }
  }

  arrowAt(ctx, line) {
    const i = Math.min(9, line.length - 2), v = this.view;
    const a = v.toScreen(v.worldFromUV(line[i][0], line[i][1]));
    const b = v.toScreen(v.worldFromUV(line[i + 1][0], line[i + 1][1]));
    const ang = Math.atan2(b[1] - a[1], b[0] - a[0]), sz = 5;
    ctx.fillStyle = 'rgba(226,232,255,0.65)';
    ctx.beginPath();
    ctx.moveTo(b[0], b[1]);
    ctx.lineTo(b[0] - sz * Math.cos(ang - 0.4), b[1] - sz * Math.sin(ang - 0.4));
    ctx.lineTo(b[0] - sz * Math.cos(ang + 0.4), b[1] - sz * Math.sin(ang + 0.4));
    ctx.closePath(); ctx.fill();
  }

  // Equipotential contours via marching squares over the Φ grid. Levels are
  // spaced in −1/Φ (so the near-body wells, where Φ→−∞, don't crowd out the
  // gentle far-field surfaces).
  drawEquipotential(ctx) {
    const g = this.grid, v = this.view;
    if (!g) return;
    let pmin = Infinity, pmax = -Infinity;
    for (let k = 0; k < g.phi.length; k++) {
      const p = g.phi[k];
      if (isFinite(p) && p < 0) { if (p < pmin) pmin = p; if (p > pmax) pmax = p; }
    }
    if (!isFinite(pmin) || pmin === pmax) return;
    const N = 11, levels = [];
    // uniform in the transformed variable w = −1/Φ ∈ (0, …)
    const wmin = -1 / pmax, wmax = -1 / pmin;
    for (let i = 1; i < N; i++) {
      const w = wmin + (wmax - wmin) * i / N;
      levels.push(-1 / w);
    }
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(120,240,210,0.34)';
    const toS = (ci, cj) => v.toScreen(v.worldFromUV(
      v.center[0] - v.spanU / 2 + v.spanU * (ci + 0.5) / g.cols,
      v.center[1] + v.spanV / 2 - v.spanV * (cj + 0.5) / g.rows));
    ctx.beginPath();
    for (const L of levels) {
      for (let j = 0; j < g.rows - 1; j++) {
        for (let i = 0; i < g.cols - 1; i++) {
          const p00 = g.phi[j * g.cols + i], p10 = g.phi[j * g.cols + i + 1];
          const p11 = g.phi[(j + 1) * g.cols + i + 1], p01 = g.phi[(j + 1) * g.cols + i];
          if (!isFinite(p00 + p10 + p11 + p01)) continue;
          const seg = msCell(p00, p10, p11, p01, L, i, j);
          for (const [ax, ay, bx, by] of seg) {
            const A = toS(ax, ay), B = toS(bx, by);
            ctx.moveTo(A[0], A[1]); ctx.lineTo(B[0], B[1]);
          }
        }
      }
    }
    ctx.stroke();
  }

  drawVectors(ctx) {
    // Quiver plot: arrow LENGTH and brightness both encode |g| (log-scaled);
    // direction is the field direction (toward mass). Thin white glyphs with a
    // dark underlay read cleanly over the viridis heatmap.
    const v = this.view, step = this.opts.gridStep;
    const span = (this.range.max - this.range.min) || 1;
    const maxLen = step * 0.9;
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    const glyph = (cx, cy, ex, ey, a, hs) => {
      ctx.beginPath();
      ctx.moveTo(cx, cy); ctx.lineTo(ex, ey);
      ctx.lineTo(ex - hs * Math.cos(a - 0.42), ey - hs * Math.sin(a - 0.42));
      ctx.moveTo(ex, ey); ctx.lineTo(ex - hs * Math.cos(a + 0.42), ey - hs * Math.sin(a + 0.42));
      ctx.stroke();
    };
    for (let sy = step / 2; sy < v.H; sy += step) {
      for (let sx = step / 2; sx < v.W; sx += step) {
        const wpt = v.toWorld(sx, sy);
        const s = this.sampleField(wpt[v.uAxis], wpt[v.vAxis]);
        if (!s) continue;
        const m = Math.hypot(s[0], s[1]);
        if (m < 1e-16) continue;
        const t = Math.min(1, Math.max(0, (Math.log10(m) - this.range.min) / span));
        const len = maxLen * (0.22 + 0.78 * t);
        const ux = s[0] / m, uy = -s[1] / m;
        const cx = sx - ux * len / 2, cy = sy - uy * len / 2;
        const ex = cx + ux * len, ey = cy + uy * len;
        const a = Math.atan2(ey - cy, ex - cx);
        const hs = Math.min(5.5, 2.4 + len * 0.2);
        ctx.strokeStyle = 'rgba(0,0,0,0.4)'; ctx.lineWidth = 2.3;
        glyph(cx, cy, ex, ey, a, hs);
        ctx.strokeStyle = `rgba(255,255,255,${0.45 + 0.5 * t})`; ctx.lineWidth = 1.2;
        glyph(cx, cy, ex, ey, a, hs);
      }
    }
  }

  drawGrid(ctx) {
    const v = this.view;
    ctx.strokeStyle = 'rgba(255,255,255,0.06)'; ctx.lineWidth = 1;
    const stepM = niceStep(60 / v.scale);
    const u0 = v.center[0] - v.spanU / 2, u1 = v.center[0] + v.spanU / 2;
    const w0 = v.center[1] - v.spanV / 2, w1 = v.center[1] + v.spanV / 2;
    ctx.beginPath();
    for (let u = Math.ceil(u0 / stepM) * stepM; u <= u1; u += stepM) {
      const s0 = v.toScreen(v.worldFromUV(u, w0)), s1 = v.toScreen(v.worldFromUV(u, w1));
      ctx.moveTo(s0[0], s0[1]); ctx.lineTo(s1[0], s1[1]);
    }
    for (let w = Math.ceil(w0 / stepM) * stepM; w <= w1; w += stepM) {
      const s0 = v.toScreen(v.worldFromUV(u0, w)), s1 = v.toScreen(v.worldFromUV(u1, w));
      ctx.moveTo(s0[0], s0[1]); ctx.lineTo(s1[0], s1[1]);
    }
    ctx.stroke();
    ctx.strokeStyle = 'rgba(255,255,255,0.16)';
    const o = v.toScreen([0, 0, 0]);
    ctx.beginPath(); ctx.moveTo(0, o[1]); ctx.lineTo(v.W, o[1]); ctx.moveTo(o[0], 0); ctx.lineTo(o[0], v.H); ctx.stroke();
  }

  // ---- overlays (drawn to visible ctx each frame) ----
  blitField() { this.ctx.drawImage(this.field, 0, 0, this.view.W, this.view.H); }

  drawSelection(s) {
    const ctx = this.ctx;
    const c = this.view.toScreen(s._origin);
    const r = Math.max(13, sourceExtent(s) * this.view.scale) + 6;
    ctx.save();
    ctx.strokeStyle = '#e9b44c'; ctx.lineWidth = 2.5;
    ctx.shadowColor = '#e9b44c'; ctx.shadowBlur = 9;
    ctx.beginPath(); ctx.arc(c[0], c[1], r, 0, 7); ctx.stroke();
    ctx.restore();
  }

  drawSources(selectedId) {
    const ctx = this.ctx;
    for (const s of this.scene.sources) {
      if (!s.visible) continue;
      const sel = s.id === selectedId;
      if (sel) this.drawSelection(s);
      ctx.save();
      if (s.type === 'point') this.drawPoint(s, sel);
      else if (s.type === 'sphere') this.drawSphere(s, sel);
      else if (s.type === 'shell') this.drawShell(s, sel);
      else if (s.type === 'ring') this.drawRing(s, sel);
      else if (s.type === 'disc') this.drawDisc(s, sel);
      else if (s.type === 'cylinder') this.drawCylinder(s, sel);
      else if (s.type === 'rod') this.drawRod(s, sel);
      else if (s.type === 'box') this.drawBox(s, sel);
      ctx.restore();
    }
  }

  localToScreen(s, local) {
    const w = P.vadd(s._origin, P.matVec(s._R, local.map((c) => c * 1000)));  // local km -> m
    return this.view.toScreen(w);
  }

  // Project a local circle (radius km, at local-z = zKm) to a screen polyline.
  ringPoly(s, radiusKm, n = 64, zKm = 0) {
    const pts = [];
    for (let i = 0; i <= n; i++) {
      const th = 2 * Math.PI * i / n;
      pts.push(this.localToScreen(s, [radiusKm * Math.cos(th), radiusKm * Math.sin(th), zKm]));
    }
    return pts;
  }

  bodyFill(ctx, c, rpx, col) {
    const grad = ctx.createRadialGradient(c[0] - rpx * 0.3, c[1] - rpx * 0.3, rpx * 0.1, c[0], c[1], rpx);
    grad.addColorStop(0, shade(col, 1.35)); grad.addColorStop(0.6, col); grad.addColorStop(1, shade(col, 0.55));
    ctx.fillStyle = grad;
  }

  drawPoint(s, sel) {
    const ctx = this.ctx, p = this.view.toScreen(s._origin);
    ctx.fillStyle = s.color; ctx.shadowColor = s.color; ctx.shadowBlur = 10;
    ctx.beginPath(); ctx.arc(p[0], p[1], 5, 0, 7); ctx.fill();
    ctx.shadowBlur = 0;
    ctx.strokeStyle = sel ? '#fff' : 'rgba(0,0,0,0.5)'; ctx.lineWidth = sel ? 2 : 1; ctx.stroke();
  }

  drawSphere(s, sel) {
    const ctx = this.ctx, c = this.view.toScreen(s._origin);
    const rpx = (s.dia / 2 * 1000) * this.view.scale;
    this.bodyFill(ctx, c, Math.max(rpx, 2), s.color);
    ctx.globalAlpha = 0.95; ctx.beginPath(); ctx.arc(c[0], c[1], Math.max(rpx, 2), 0, 7); ctx.fill(); ctx.globalAlpha = 1;
    ctx.strokeStyle = sel ? '#fff' : 'rgba(0,0,0,0.5)'; ctx.lineWidth = sel ? 2 : 1; ctx.stroke();
  }

  drawShell(s, sel) {
    const ctx = this.ctx, c = this.view.toScreen(s._origin);
    const rpx = (s.dia / 2 * 1000) * this.view.scale;
    ctx.strokeStyle = s.color; ctx.lineWidth = sel ? 3 : 2;
    ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.arc(c[0], c[1], Math.max(rpx, 2), 0, 7); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawRing(s, sel) {
    const ctx = this.ctx;
    const poly = this.ringPoly(s, s.dia / 2);
    ctx.strokeStyle = sel ? '#fff' : s.color; ctx.lineWidth = 3;
    ctx.beginPath(); poly.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.stroke();
  }

  drawDisc(s, sel) {
    const ctx = this.ctx;
    const r = s.dia / 2, t = (s.thick || 0) / 2;
    const trace = (poly) => { ctx.beginPath(); poly.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.closePath(); };
    if (t > 0) {
      // Two projected faces offset by the thickness give a 3-D puck: flat when
      // seen face-on, a clear slab when seen edge-on.
      const back = this.ringPoly(s, r, 64, -t), front = this.ringPoly(s, r, 64, t);
      trace(back); ctx.fillStyle = shade(s.color, 0.5); ctx.globalAlpha = 0.6; ctx.fill(); ctx.globalAlpha = 1;
      ctx.strokeStyle = shade(s.color, 0.6); ctx.lineWidth = 1; ctx.stroke();
      trace(front); ctx.fillStyle = s.color; ctx.globalAlpha = 0.5; ctx.fill(); ctx.globalAlpha = 1;
      ctx.strokeStyle = sel ? '#fff' : shade(s.color, 0.75); ctx.lineWidth = sel ? 2 : 1.5; ctx.stroke();
    } else {
      trace(this.ringPoly(s, r));
      ctx.fillStyle = s.color; ctx.globalAlpha = 0.5; ctx.fill(); ctx.globalAlpha = 1;
      ctx.strokeStyle = sel ? '#fff' : shade(s.color, 0.7); ctx.lineWidth = sel ? 2 : 1.5; ctx.stroke();
    }
  }

  drawCylinder(s, sel) {
    const ctx = this.ctx;
    const r = s.dia / 2, L = s.len;
    const c = [[-r, 0, -L / 2], [r, 0, -L / 2], [r, 0, L / 2], [-r, 0, L / 2]].map((p) => this.localToScreen(s, p));
    ctx.beginPath(); ctx.moveTo(c[0][0], c[0][1]);
    for (let i = 1; i < 4; i++) ctx.lineTo(c[i][0], c[i][1]);
    ctx.closePath();
    const grad = ctx.createLinearGradient(c[0][0], c[0][1], c[1][0], c[1][1]);
    grad.addColorStop(0, shade(s.color, 0.6)); grad.addColorStop(0.5, shade(s.color, 1.2)); grad.addColorStop(1, shade(s.color, 0.6));
    ctx.fillStyle = grad; ctx.globalAlpha = 0.9; ctx.fill(); ctx.globalAlpha = 1;
    ctx.strokeStyle = sel ? '#fff' : 'rgba(0,0,0,0.5)'; ctx.lineWidth = sel ? 2 : 1; ctx.stroke();
  }

  drawRod(s, sel) {
    const ctx = this.ctx, L = s.len;
    const a = this.localToScreen(s, [0, 0, -L / 2]), b = this.localToScreen(s, [0, 0, L / 2]);
    ctx.strokeStyle = sel ? '#fff' : s.color; ctx.lineWidth = 5; ctx.lineCap = 'round';
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
    ctx.lineCap = 'butt';
    ctx.fillStyle = shade(s.color, 1.3);
    for (const p of [a, b]) { ctx.beginPath(); ctx.arc(p[0], p[1], 3.5, 0, 7); ctx.fill(); }
  }

  drawBox(s, sel) {
    const ctx = this.ctx;
    const [w, h, d] = s.size;
    const corners = [
      [-w / 2, -h / 2, -d / 2], [w / 2, -h / 2, -d / 2], [w / 2, h / 2, -d / 2], [-w / 2, h / 2, -d / 2],
      [-w / 2, -h / 2, d / 2], [w / 2, -h / 2, d / 2], [w / 2, h / 2, d / 2], [-w / 2, h / 2, d / 2],
    ].map((c) => this.localToScreen(s, c));
    const poly = (idx, fill) => {
      ctx.beginPath();
      idx.forEach((i, k) => { const p = corners[i]; k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]); });
      ctx.closePath(); ctx.fillStyle = fill; ctx.fill();
    };
    ctx.globalAlpha = 0.9;
    poly([0, 1, 2, 3], shade(s.color, 0.8));
    poly([4, 5, 6, 7], shade(s.color, 1.15));
    ctx.globalAlpha = 1;
    ctx.strokeStyle = sel ? '#fff' : 'rgba(0,0,0,0.5)'; ctx.lineWidth = sel ? 2 : 1;
    const edges = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];
    ctx.beginPath();
    for (const [i, j] of edges) { ctx.moveTo(corners[i][0], corners[i][1]); ctx.lineTo(corners[j][0], corners[j][1]); }
    ctx.stroke();
  }

  clear() { this.ctx.fillStyle = '#07080d'; this.ctx.fillRect(0, 0, this.view.W, this.view.H); }
}

// Marching-squares: return line segments (in fractional cell coords) where the
// contour of value `L` crosses the cell with corners p00 (i,j) p10 (i+1,j)
// p11 (i+1,j+1) p01 (i,j+1). Coordinates returned are grid indices.
function msCell(p00, p10, p11, p01, L, i, j) {
  const b = (p00 > L ? 1 : 0) | (p10 > L ? 2 : 0) | (p11 > L ? 4 : 0) | (p01 > L ? 8 : 0);
  if (b === 0 || b === 15) return [];
  const lerp = (a, bb) => (L - a) / ((bb - a) || 1e-30);
  const top    = () => [i + lerp(p00, p10), j];
  const right  = () => [i + 1, j + lerp(p10, p11)];
  const bottom = () => [i + lerp(p01, p11), j + 1];
  const left   = () => [i, j + lerp(p00, p01)];
  const seg = (A, B) => [[A[0], A[1], B[0], B[1]]];
  switch (b) {
    case 1: case 14: return seg(left(), top());
    case 2: case 13: return seg(top(), right());
    case 3: case 12: return seg(left(), right());
    case 4: case 11: return seg(right(), bottom());
    case 6: case 9:  return seg(top(), bottom());
    case 7: case 8:  return seg(left(), bottom());
    case 5:  return [...seg(left(), top()), ...seg(right(), bottom())];
    case 10: return [...seg(top(), right()), ...seg(left(), bottom())];
  }
  return [];
}

function niceStep(x) {
  const e = Math.pow(10, Math.floor(Math.log10(x)));
  const m = x / e;
  const nm = m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10;
  return nm * e;
}

// Lighten (f>1) or darken (f<1) a #rrggbb colour.
function shade(hex, f) {
  const n = parseInt(hex.slice(1), 16);
  const r = Math.min(255, (n >> 16) * f), g = Math.min(255, ((n >> 8) & 255) * f), b = Math.min(255, (n & 255) * f);
  return `rgb(${r | 0},${g | 0},${b | 0})`;
}

export { viridis };
