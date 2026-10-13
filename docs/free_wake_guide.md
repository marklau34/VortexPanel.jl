# Free wake: choosing the settings

A practical guide to `FreeWakeSim`: what each setting does, how sensitive the lift is to it, and what to pick. Unless noted, the numbers are from MESA's reference wing (NACA 0015, rectangular, aspect ratio 3.88) at α = 5° on the L1 mesh, started impulsively with a free wake. Most come from `validation/freewake/run_freewake.jl`; the sweeps of the conversion distance and of long trailing-edge panels were separate runs of the same case.

## Quick start

```julia
using VortexPanel, LinearAlgebra
g  = build_geometry(mesh)
dm = build_dofmap(g)
dt = VortexPanel.Wake.te_panel_size(g) / norm(V∞)  # one row per step, as long as a trailing-edge panel
sim = FreeWakeSim(g, dm, V∞; opts=FreeWakeOptions(dt=dt))
run!(sim, nsteps)                                   # sim.history holds the loads of every step
```

The defaults are good for most runs. The ones you are most likely to change are `dt` and, at small `dt`, `merge_rows`.

| setting | default | what it controls |
|---|---|---|
| `dt` | required | the time step; each step sheds one panel row about \|V∞\| dt long |
| `te_len` | 0.3 \|V∞\| dt | length of the implicit trailing-edge panel (the Kutta condition) |
| `max_rows` | `:auto` | free panel rows kept before conversion to particles |
| `core` / `overlap` | `:auto` / 1.5 | the particles' core radius |
| `merge_rows` | 1 | rows merged into one particle pair per strip |
| `advection` | `:free` | `:freestream` gives a prescribed wake, with no roll-up |
| `delete_beyond` | `Inf` | removes particles beyond this many spans downstream |

## How long to run

Lift approaches its steady value slowly. The starting vortex left behind at the start keeps inducing downwash, and that fades only like (span / wake length)². With the default settings (Kutta–Joukowski lift):

| wake length | lift, relative to the steady rigid wake |
|---|---|
| 1 span | −6% |
| 2 spans | −1.4% |
| 5 spans | −0.4% |

For a steady answer:
- Run to about 5 spans.
- Fit CL(L) = CL∞ + a/L² to the last half of the run.

Extrapolated this way, the free wake settles 0.24–0.35% below the rigid wake, at L1 and L2. That difference is roll-up; the rest of the gap at 5 spans is the slow approach.

## Core size

**Rule:** a particle's core is `overlap` × half the longest edge of the element it came from (`core = :auto`, `overlap = 1.5`).

Each converted quad becomes two particles, about half the quad's longer side apart. The rule therefore keeps neighbouring cores overlapping in every direction, which a smooth vortex sheet needs.
- **Too small** leaves gaps, and the wake velocities get noisy.
- **Too large** smears the tip vortices.

The rule follows the mesh, so nothing needs tuning when you refine it or change dt.

| `overlap` | lift change | peak wake velocity |
|---|---|---|
| 1.0 | −0.03% | 1.028 \|V∞\| |
| 1.5 (default) | — | 1.020 \|V∞\| |
| 2.5 | +0.01–0.02% | 1.011 \|V∞\| |

- **Pick:** keep `:auto`. Raise `overlap` (2–2.5) if roll-up looks noisy; lower it (toward 1) if a tip vortex looks too diffuse.
- **The core matters for lift only when particles are close to the wing** (next section).
- **Alternatives:** `core = 0.1` gives every particle the same fixed core, like MESA's `r_vortex`. Use it only to match MESA. A fixed core doesn't follow the mesh, and MESA's default (1e-5) is effectively no smoothing at all. `core = :scaled` is the earlier area-based rule; it leaves spanwise gaps at small dt.

## When to convert: `max_rows`

Particles created close to the wing push the lift up. Their smoothed velocity field underpredicts the near wake's downwash on the wing, so the wing sees a slightly larger effective angle of attack. Lift change against the same run without conversion:

| free rows kept | particles created (cores behind the trailing edge) | lift change, extrapolated | at 1 span |
|---|---|---|---|
| 0 (convert immediately) | 0.24 | +1.9% | +4.1% |
| 1 | 0.55 | +0.6% | +2.0% |
| 4 | 1.5 | +0.09% | +0.5% |
| 10 (`:auto`) | 3.3 | +0.02% | +0.15% |

- **Pick:** keep `max_rows = :auto`. It keeps as many rows as put new particles three cores behind the trailing edge, whatever dt is (10 rows at the default dt, 38 at dt/4, 77 at dt/8).
- **A fixed number** works like MESA's `maxwakelen`. `FreeWakeSim` warns when it would create particles within two cores of the wing, and says how many rows avoid that.
- **Smaller cores reduce the error of immediate conversion** but don't remove it. With the earlier area-based cores, halving the core cut it from +2.3% to +1.3%, and about a quarter of the core to +1.0%, with a noisier wake (peak velocity 1.07 |V∞|). Keeping a few rows is the better fix.

## Trailing-edge panel length: `te_len`

The trailing-edge panel carries the Kutta condition and is solved implicitly. Everything behind it couples in explicitly.
- **Steady flow:** the solution does not depend on te_len. The coupling reproduces the implicit rigid-wake solution to rounding, for panels from 0.25 to 4 trailing-edge panels long.
- **Start-up:** te_len sets where each step's newly shed vorticity sits, so it matters during the start-up transient. Kutta–Joukowski lift (prescribed wake) at two wake lengths:

| te_len | at 1 span | at 2 spans |
|---|---|---|
| 0.1 \|V∞\| dt | 0.313 | 0.336 |
| 0.3 \|V∞\| dt (default) | 0.324 | 0.338 |
| 1.0 \|V∞\| dt | 0.331 | 0.339 |

The default, 0.3 \|V∞\| dt, is about where Katz & Plotkin place the newest shed vortex. It scales with dt, which is what time accuracy needs.

**Don't use a long panel instead of free rows.** With immediate conversion, a long panel does keep particles away from the wing, but each version is worse than keeping free rows (+0.02% extrapolated, +0.15% at 1 span). Lift change against the default run without conversion:

| te_len (with `max_rows = 0`) | lift change, extrapolated | at 1 span | peak wake velocity |
|---|---|---|---|
| 0.3 \|V∞\| dt | +1.9% | +4.1% | 1.029 |
| 1 core | +0.36% | +4.3% | 1.044 |
| 3 cores (0.8 chord) | +0.20% | +4.9% | 1.055 |

The panel is rigid, so the near wake can't roll up along it. It also holds one strength along its whole length, which puts each step's shed vorticity at its far end. And the panel is too long for a single smoothed particle per triangle to represent it well.

**Pick:** leave te_len at its default.

**MESA comparison:** MESA's lift depends strongly on its trailing-edge panel. With a prescribed wake at 2 spans on this wing, it goes from 0.281 to 0.312 for te_len from 0.015 to 0.117 m, against ours at 0.335 to 0.337. To compare, give both codes the same te_len, dt, row count and a fixed core (`core = r_vortex`).

## Small time steps and merging

A small dt makes each row thin. Lift is unaffected as long as the conversion distance stays fixed (`max_rows = :auto`). Cost grows quickly, though: there are more steps, more rows, and particle sums over all pairs.

`merge_rows = :auto` merges each strip's converted rows into one particle pair once they span about the strip's width. Narrow tip strips keep their fine particles. Total vorticity is conserved exactly.

| run (2 spans) | lift change | particles | time |
|---|---|---|---|
| dt (one trailing-edge panel) | — | 3,188 | 3 s |
| dt / 4 | −0.03% | 12,688 | 56 s |
| dt / 4, merged | +0.03% | 3,176 | 31 s |
| dt / 8 | +0.03% | 25,304 | 268 s |
| dt / 8, merged | +0.08% | 5,054 | 133 s |

**Pick:** use `merge_rows = :auto` whenever dt is below about one trailing-edge panel, or for long runs. It is off by default (validated results use 1).

## Recipes

- **Steady lift:**
  - dt = one trailing-edge panel, defaults otherwise;
  - run about 5 spans;
  - extrapolate in 1/L².
- **Start-up transients:** a smaller dt, with `te_len` at its default (it follows dt), `max_rows = :auto` and `merge_rows = :auto`.
- **Matching MESA:** `te_len` = MESA's `te_len`, `max_rows` = `maxwakelen`, `core` = `r_vortex`, the same dt.
- **Long wakes:** `merge_rows = :auto`, and `delete_beyond` (in spans). Deleting particles truncates the wake. The error falls like (span / D)² and is negligible beyond about 10 spans: for the rigid wake, 50 spans instead of 20 changed CL by under 1e-4.
- **Visualization:** `run!(sim, n; series=FreeWakeSeries(dir; every=4))`, or `validation/examples/free_wake.jl`. Open the four `.pvd` files in ParaView.

## Limitations

- **No vortex stretching.** Particle strengths never change, so the total vorticity drifts: 0.7% (L1) and 4% (L2) of the peak circulation × span after 5 spans. The impulse force route depends on that drift and is unreliable in free mode. The pressure and Kutta–Joukowski routes are fine.
- **Cost.** Particle velocities are summed over all pairs. L2 to 5 spans (33k particles) takes about 7 minutes on 16 threads. L3, or wakes much longer than 5 spans, need merging, deletion or a faster summation method.
