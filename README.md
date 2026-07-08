# 🪐 Gravity — gravitational field lab

An interactive, **physics-accurate** Newtonian gravity simulator that runs
entirely in the browser. Build a scene from point masses, planets, hollow shells,
rings, discs, cylinders, rods and boxes; see the field as a heat-map, field
lines, equipotential contours and a vector quiver; drag a probe to read the exact
**g** vector and potential Φ anywhere; read the exact tidal force and torque
between bodies; and launch a body into a **full mutual N-body** simulation where
every mass pulls on every other — slingshot a rock past a planet and watch it tug
the planet back, or let two bodies orbit their shared barycentre.

No installation, no accounts, no telescope — just open the page. Works on desktop
and touch devices.

👉 **[Live demo](https://robertjohnhanna.github.io/gravlab/)**

## Why it's trustworthy

The field engine is not a toy approximation. It uses the same closed-form
solutions as professional tools (geophysics/gravimetry), and every routine is
checked against an independent analytical limit:

- **Point masses & planets** — exact `g = −GM r̂/r²`; uniform spheres are exact
  *inside and out* via Newton's shell theorem (linear interior field).
- **Spherical shells** — exact: a point mass outside, **exactly zero field**
  inside.
- **Rings** — exact closed form via complete elliptic integrals K and E.
- **Rods** — exact finite line-mass closed form (validated against the
  infinite-rod limit 2Gλ/d).
- **Boxes** — exact uniform-prism field via Nagy's closed form, validated against
  the Laplace constraints (∇·g = 0, ∇×g = 0) and the point-mass far field.
- **Discs & cylinders** — exact rings summed over a stack.
- **Orbits** — a symplectic velocity-Verlet integrator, validated against the
  Kepler period and energy conservation.

Run the checks yourself:

```bash
npm test        # 64 physics assertions, all from first principles
```

Full derivations, formulas and references: **[docs/PHYSICS.md](docs/PHYSICS.md)**.

## What you can do

| Feature | Notes |
|---|---|
| Masses | **point**, **sphere** (planet), **shell**, **ring**, **disc**, **cylinder**, **rod** and **box**; set mass (in kg, or seed from a named body — Moon … Sun), size, position and full 3-axis orientation |
| Field visualisation | log-scaled \|g\| heat-map, gravitational **field lines**, **equipotential** contours (Φ), a scientific vector quiver, and a reference grid |
| Slice control | view the XZ / XY / YZ plane at any offset through the 3-D scene |
| Field probe | drag the ⊕ pin to read the full 3-D **g** vector and potential Φ, in m/s² / mGal / µGal / g₀ |
| Force & torque | **exact** net gravitational force and **gravity-gradient (tidal) torque** on the selected body (real volume integral — no far-field approximation), valid at any separation |
| N-body lab | launch a body (or set bodies' velocities) into a full mutual-gravity simulation — every mass moves, momentum/energy conserved; live speed, specific orbital energy and bound/escape classification |
| Scenarios | one-click presets: **solar system** (Sun + 8 planets, real masses/distances/speeds), binary system, planet + moon, spherical shell, ring world, tidal rod, slingshot fly-by |
| Collisions | a launched body that strikes a placed body is absorbed in a perfectly inelastic impact — momentum conserved, the struck body recoils and accretes the mass |

## Using it

The **canvas** is in the centre. Data read-outs and the orbit/layer controls are
on the left; the mass palette, object list and parameter inspector are on the
right. (On a phone the canvas becomes the hero and the panels stack below.)

**Interaction**

- **Drag a body** to move it in the view plane; **drag empty space** to pan.
- **Scroll** or **pinch** to zoom (centred on the cursor / pinch point); the
  on-canvas buttons do **＋ / − / Fit / reset**.
- **Drag the ⊕ probe pin** anywhere to read **g** and Φ at that point.
- With a body selected: **arrow keys** nudge it (Shift = ×5, snaps to the grid
  when **Snap** is on), **Delete** removes it.
- Add masses from the palette; tune every parameter in the inspector; toggle
  visibility or delete from the object list.
- **Launch** a body from the ⊕ probe (drag the red tip to aim), or give bodies a
  velocity in the inspector and press **Play** — every mass then moves under
  mutual gravity. **Reset** restores the placed layout.

## Running locally

It's plain static files (ES modules) — serve the folder over HTTP:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

(Opening `index.html` via `file://` won't work: browsers block ES-module imports
on that protocol.)

## Hosting on GitHub Pages

1. Push to GitHub (this repo).
2. **Settings → Pages → Build and deployment → Source: Deploy from a branch**.
3. Pick the branch and the `/ (root)` folder, then save.
4. The site appears at `https://<user>.github.io/gravlab/` within a minute.

There is no build step; the files are served as-is.

## Project layout

```
index.html          page markup + panel layout
styles.css          styling (incl. responsive / mobile layout)
src/physics.js      SI-unit gravity engine — vectors, sphere/shell/ring/rod/box, Verlet
src/sources.js      user objects + Scene (fields, potential, exact force/torque)
src/render.js       canvas visualisation (heat-map, field lines, equipotentials, quiver)
src/main.js         UI, interaction, scene state, orbit simulation
tests/selftest.mjs  physics verification (npm test)
docs/PHYSICS.md     derivations, formulas, references
```

## Accuracy & limits

The honest list is in [docs/PHYSICS.md](docs/PHYSICS.md#approximations--caveats).
In short: fields, potentials, forces and torques are exact for the idealised
bodies (forces use the real volume integral, not a far-field approximation);
discs/cylinders are exact rings summed over a finite stack; the N-body
simulation treats bodies as point masses at their centres (exact for spheres);
motion is Newtonian; and the display is a
2-D slice of a fully 3-D calculation.

## License

MIT.
