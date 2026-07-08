# Physics notes

Everything here is Newtonian gravitation computed from first principles in **SI
units** (metres, kilograms, seconds). The field **g** is the gravitational
acceleration [m/s²] a test mass feels — it always points *toward* mass — and the
scalar potential **Φ** [J/kg] satisfies **g = −∇Φ**. Both obey superposition, so
the total field/potential of a scene is the exact sum of its members'. The engine
lives in [`src/physics.js`](../src/physics.js) and every routine below is checked
against an independent analytical limit in
[`tests/selftest.mjs`](../tests/selftest.mjs) (run `npm test`).

## Constants

| symbol | value | meaning |
|---|---|---|
| G | 6.674×10⁻¹¹ m³·kg⁻¹·s⁻² | gravitational constant |
| g₀ | 9.80665 m/s² | standard gravity (unit reference only) |

## Point mass

The fundamental source:

```
g(r) = −G M (r − r₀)/|r − r₀|³        Φ(r) = −G M/|r − r₀|
```

## Uniform solid sphere — exact everywhere (shell theorem)

Newton's shell theorem makes a spherically-symmetric body one of the few whose
field is closed-form *inside and out*:

- **Outside** (r ≥ R): identical to a point mass at the centre.
- **Inside** (r < R): only the enclosed mass pulls, so `g = −G M r/R³` (linear in
  r) and `Φ = −G M(3R² − r²)/(2R³)`.

The two forms agree at the surface, so **g** and **Φ** are continuous there.

## Thin spherical shell — exact everywhere

- **Outside** (r ≥ R): a point mass at the centre.
- **Inside** (r < R): the field is **exactly zero** and the potential is the
  constant `Φ = −G M/R` — a flat-bottomed well. (A body inside a shell is
  weightless with respect to the shell, no matter where it sits.)

## Uniform ring (hoop) — closed form via elliptic integrals

For a ring of mass M and radius a in the local z = 0 plane, with cylindrical
radius ρ and axial coordinate z, and `D² = (ρ+a)²+z²`, `α² = (ρ−a)²+z²`,
`m = 4aρ/D²`:

```
Φ   = −(2GM/π) K(m)/D
g_z = −(2GM z)/(π D α²) · E(m)
g_ρ = −∂Φ/∂ρ    (a closed combination of K, E and their m-derivative)
```

with K, E the complete elliptic integrals (evaluated by the arithmetic–geometric
mean). One evaluation replaces summing dozens of point masses around the hoop —
much faster *and* exact. Verified against the on-axis formula
`g_z = −GM z/(a²+z²)^{3/2}` and, off-axis, against direct numerical integration
around the ring, and by confirming **g = −∇Φ**.

## Uniform disc and cylinder — ring stacks

A **disc** is decomposed into coaxial rings reproducing the areal density
σ = M/πR²; a **cylinder** is a stack of such discs along its axis. Each ring uses
the exact elliptic form above, so the only approximation is the finite number of
rings. Verified against the analytic on-axis disc field
`g_z = −2πGσ(1 − z/√(z²+R²))` and the far-field point-mass limit.

## Uniform rod (line mass) — closed form

A rod of linear density λ and length L along the local z axis, at cylindrical
radius ρ with end-distances r_A, r_B:

```
g_ρ = −(Gλ/ρ)[(z+L/2)/r_A − (z−L/2)/r_B]
g_z =  Gλ[1/r_A − 1/r_B]
Φ   = −Gλ ln[(z+L/2 + r_A)/(z−L/2 + r_B)]
```

Verified against the infinite-rod limit `g_ρ → −2Gλ/ρ`, against direct
integration, and by **g = −∇Φ**.

## Uniform rectangular prism (box) — Nagy's closed form

The gravitational field and potential of a uniform-density cuboid are known in
closed form (Nagy 1966; Nagy, Papp & Benedek 2000 — the standard result in
gravimetry). For half-extents (a, b, c) and density ρ, with corners
X = pₓ±a, Y = p_y±b, Z = p_z±c and r the distance to each corner:

```
g_x = G ρ Σ (−1)^{i+j+k} [ Y ln(Z+r) + Z ln(Y+r) − X·atan2(YZ, X r) ]   (+ cyclic)
```

Verified against the far-field point-mass limit `GM/r²`, the Laplace constraints
∇·g = 0 and ∇×g = 0 outside the body, and **g = −∇Φ**.

## Forces and torques — exact, no far-field approximation

The net force and torque on the selected body from all *other* sources is the
force density **ρ g_ext** integrated over the body, using **g_ext** (the field of
the other sources only) — never a point-mass/far-field approximation. It is valid
at any separation, as long as the bodies don't interpenetrate.

- **Spherically-symmetric bodies (point / sphere / shell):** by Newton's theorem
  the net force is *exactly* `F = M·g_ext(centre)` and the self-torque is zero —
  no integration needed.
- **Extended bodies (rod / ring / disc / cylinder / box):** integrated over the
  body's mass, `F = Σ dm·g_ext`, `τ = Σ (r − r₀) × dm·g_ext`. The torque is the
  real **gravity-gradient (tidal) torque** — the mechanism that tidally locks
  moons: an elongated body is twisted toward radial alignment with the field.

Verified against the point-mass attraction `GM₁M₂/d²` between two spheres,
Newton's third law (|F_AB + F_BA| ≈ 0 to machine precision), zero self-torque on
a sphere, a nonzero gravity-gradient torque on a tilted rod, and a torque
equilibrium for a radially-aligned rod. Overlapping bodies are **refused** rather
than given a meaningless value.

## Dynamics — full mutual N-body

When the simulation runs, **every** body (each placed source and each launched
body) moves under the gravity of all the others — there is no fixed background.
The pairwise interaction is the exact Newtonian point-mass law

```
a_i = −G Σ_j m_j (x_i − x_j)/|x_i − x_j|³
```

which, by the shell theorem, is **exact for spheres** as long as they don't
overlap. (Once two bodies touch, d < r_i+r_j, the pair crosses over to the
uniform-sphere interior law ∝ d so a direct hit stays finite instead of blowing
up.) Because the pair force is equal-and-opposite, **total momentum is conserved
to machine precision** — a moon-sized body flung past a planet tugs the planet
back, and a captured body and its primary orbit their common **barycentre**
continuously. A projectile's own mass never changes *its own* trajectory (the
equivalence principle) — only its pull on everything else.

Integration uses the **velocity-Verlet (leapfrog)** scheme applied to the whole
system: symplectic and time-reversible, so orbital energy is conserved over long
runs and closed orbits stay closed. The sub-step is **adaptive** — shortened in
proportion to the shortest pairwise encounter/orbital time — so a fast, close
passage stays resolved. Verified against a circular binary's constant separation
and Kepler period `T = 2π√(d³/G(m₁+m₂))`, a stationary barycentre, and momentum
and energy conservation for 2- and 3-body systems.

The specific orbital energy `ε = ½v² + Φ` classifies a body against the placed
sources: ε < 0 is a **bound** (elliptical) orbit, ε ≥ 0 is **unbound**
(parabolic/hyperbolic escape).

## Approximations & caveats

- The N-body interaction treats each body as a **point mass at its centre of
  mass** (exact for spheres; the far-field limit for extended bodies). Rotational
  dynamics and time-varying tidal torques during an orbit are not integrated —
  the exact tidal torque is reported for the *instantaneous* configuration in the
  Force &amp; torque panel.
- The field/heatmap is drawn from the **placed sources**; a launched body's own
  field is not added to the map (its gravity still acts on the dynamics).
- **Discs and cylinders** are exact rings summed over a finite stack; the field
  converges to the true solid as the stack density rises.
- Motion is **Newtonian** (non-relativistic; no gravitomagnetism, no GR
  precession).
- The canvas shows a **2-D slice** through the 3-D scene. Field lines trace the
  in-plane projection of **g**; the probe reports the true out-of-plane component
  separately. The |g| heatmap is **log-scaled** and clamped to ~6 decades so a
  single hot cell next to a mass doesn't wash out the rest.
- Idealised bodies with sharp edges (box, rod, ring) have integrable field
  singularities exactly on an edge/wire; points off the edges are finite and
  accurate.
