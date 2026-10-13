# VortexPanel.jl

Galerkin linear-vortex panel solver. Unknowns are nodal doublet strengths on flat triangles. Each triangle carries a constant vortex sheet γ = n×∇μ. Sources are prescribed. The Kutta condition is applied at the trailing-edge nodes, and the wake is rigid.

The design is in `Plan Galerkin Linear-Vortex Panel Solver (Rigid Wake Phase).md`. The reference solver is [MESA.jl](../MESA.jl), a Morino code that is used unmodified from the validation scripts.

## Status

Steps 1–5 are done, with validation stages 1–8. They cover:

- **`VortexPanel.Geometry`** (step 1): turns surface meshes of quads and triangles into a closed, outward, classified triangulation. It welds nodes, splits quads, detects T-junctions, gaps and non-matching trailing edges, and caps open tips without adding nodes. It also detects Kutta edges and reports mesh quality.
- **`VortexPanel.MeshTools`** (step 1): generates closed test bodies and seeded, reproducible abused meshes. The abuse covers jitter, spanwise skew, T-junctions, gaps, slivers and non-matching trailing edges. Jitter moves each node by up to δ times its clearance, its smallest height over an edge of an incident element, so elements skew without collapsing.
- **MESA baseline** (step 1): a rigid-wake baseline on the reference wing family.
- **`VortexPanel.DOFMap`** (step 2): the constraint matrix C, with μ_slots = C q, and the nodal Kutta rows K, with μ_w = K q.
- **`VortexPanel.Kernels`** (step 3): the velocity of a unit-strength source triangle, w = (1/4π)∫(x−y)/|x−y|³ dA, the one kernel of the method (sources σw, vortex sheets γ×w).
- **`VortexPanel.Elements`, `Assembly`, `Post`** (step 4): the flat-triangle element and quadrature tiers, the Galerkin assembly with the bordered dense solve, and surface velocity, Cp and leak. Validated on the sphere (validation stage 3) and on ellipsoids (stage 4) in uniform flow, on clean and abused meshes.
- **`VortexPanel.Wake`, lifting assembly and forces** (step 5): rigid wake strips with their tip filaments, the wake's columns in the Galerkin system, and both force routes (pressure integration; Kutta–Joukowski with the Trefftz plane). Validated on the reference wing (stage 5).
- **Comparison with MESA** (stage 6): both codes on the same meshes, with an exact two-dimensional referee for the section. The zero-tip Kutta condition (`DOFMapOptions(tip_kutta=:zero)`) became the default here.
- **Isoparametric quads** (stage 7): the element the design doc defers. `MixedMesh` keeps the input's quads whole, with bilinear μ and numerically integrated inner integrals, and is compared with the split triangles on the cube-sphere and the wing.
- **Test weights** (stage 8): the hat functions (default), median dual cells and one-rings, on the wing with clean and abused meshes.
- **`VortexPanel.FreeWake`** (phase 2): the design doc's first deferred item, a free wake shed into vortex particles.
  - The time stepping mirrors MESA's: an implicit trailing-edge panel, and free panel rows advected by forward Euler.
  - The conversion is the design doc's: local to each wake triangle, one particle per triangle (two per quad).
  - The wake couples into the body through its velocity alone, with the matrix factored once.

## Layout

```
src/geometry/    Geometry: build_geometry(mesh; opts) -> SurfaceGeometry, print_report, write_vtk, add_caps
src/dofmap/      DOFMap: build_dofmap(g; opts) -> DOFMapping (C, K, slots), reduce_system, bordered_rows, wake_rows
src/kernels/     Kernels: SourceTriangle, source_velocity, sheet_velocity, segment_velocity (vortex filament)
src/elements/    Elements: FlatTriangle (hat basis, γ_b), quadrature rules and tiers (TierRules); IsoQuad (bilinear),
                 project, quad_field (its inner integrals), rules on the parameter square (SquareTiers)
src/wake/        Wake: build_wake(g, dm, V∞; opts) -> RigidWake (strips, tip filaments), WakeOptions, te_panel, write_vtk
src/assembly/    Assembly: assemble(g, dm, els, σ, V∞; opts, wake) -> A, b (tiers incl. contact_rule); bordered_system,
                 solve_galerkin, factorize_galerkin; explicit_rhs (known vortex elements in b); consistency;
                 MixedMesh(g, dm) and assemble(g, dm, mesh::MixedMesh, V∞; opts, wake)
src/post/        Post: surface_velocity (singularity sums), gradient_velocity, pressure_coefficient, nodal_average,
                 body_velocity; pressure_force, kutta_joukowski_force, trefftz_drag, spanwise_loading; surface_pressure
src/freewake/    FreeWake: FreeWakeSim(g, dm, V∞; opts) with step!/run!, FreeWakeOptions, convert_oldest!,
                 write_free_wake_vtk (the Rosenhead–Moore Blob kernel is in Kernels)
src/meshtools/   MeshTools: icosphere, cube_sphere, uv_sphere, abuse, skew_wing, ...
test/            Pkg.test() (needs the local PAGE release tree)
validation/      reference wing, MESA baseline, geometry report, sphere, ellipsoid, wing, comparison, element and
                 test-weight studies, and examples (own environment); common/exact_flows.jl: the exact ellipsoid flow;
                 common/closed_body.jl: the closed-body study tables; common/wing_case.jl: solve_wing, a wing with its
                 wake and both force routes, and its midspan section; common/panel2d.jl: the 2-D linear-vortex panel
                 method; common/svgplot.jl: dependency-free SVG plots
```

## Environments

The local checkouts are pinned with `[sources]`:

- **Package (root `Project.toml`)**: MeshLibrary is pinned to `../MeshLibrary.jl`. The registered 0.4.5 predates commit 1df0407, which keeps element sets through `merge_duplicate_nodes`; without that commit, PAGE's `trailing_edge_up/low` sets are lost. Drop the pin once a newer MeshLibrary is registered.
- **`test/`**: PAGE is pinned to `../../_tools_release/PAGE`. The registered PAGE 0.6.8 tags the trailing edge differently.
- **`validation/`**: MESA, PAGE and MeshLibrary are pinned to the same checkouts the `aero` environment uses.

## Running

```
julia --project -e "using Pkg; Pkg.test()"
julia --project=validation validation/common/test_loads.jl                   # post-processing unit tests
julia -t auto --project=validation validation/morino_baseline/run_clean.jl   # MESA baseline: levels 1-3, open and capped tips
julia -t auto --project=validation validation/morino_baseline/smoke_abuse.jl # MESA on each mesh abuse
julia --project=validation validation/geometry/report_geometry.jl            # geometry and DOF map tables, and VTK
julia --project=validation validation/kernels/compare_ase.jl                 # source kernel against AeroSingularityElements
julia -t auto --project=validation validation/sphere/run_sphere.jl           # validation stage 3: sphere in uniform flow
julia -t auto --project=validation validation/ellipsoid/run_ellipsoid.jl     # validation stage 4: ellipsoids in uniform flow
julia -t auto --project=validation validation/wing/run_wing.jl               # validation stage 5: the wing with a rigid wake
julia -t auto --project=validation validation/comparison/run_comparison.jl   # validation stage 6: against MESA (~25 min)
julia -t auto --project=validation validation/elements/run_elements.jl       # validation stage 7: triangles against quads (~2 h)
julia -t auto --project=validation validation/test_weights/run_test_weights.jl # validation stage 8: test weights (~20 min)
julia -t auto --project=validation validation/freewake/run_freewake.jl     # phase 2: free wake and particles, against MESA (~25 min)
julia -t auto --project=validation validation/examples/abused_wing.jl level=2 abuse=all   # example: abused wing, Cp in VTK
```

`abused_wing.jl` solves the clean and an abused reference wing and writes both for ParaView: colour the `*_tris.vtu` files by `Cp`. The `abuse` argument is one of `all`, `jitter`, `slivers`, `gaps`, `tjunctions`, `tjunctions_midspan`, `skew`, `te_misaligned` or `clean`; `count=N` sets how many quads the T-junction abuses split (default 4). The script header lists every exported field.

Hanging nodes cost nothing measurable. At L2, 64 quads split within |η| ≤ 0.2 of midspan leave 128 hanging nodes and change CL by −0.001% (4 and 16 splits: under 1e-4%). The column defect barely moves, from 6.6e-5 to 7.0e-5. MESA fails on all three meshes (`calc_CHTLS_pan!`).

Outputs go to `validation/output/`, which is not tracked:

- `morino/summary.csv`, which records the git SHAs of the checkouts;
- spanwise loading, Cp and trailing-edge μ CSVs, and VTK;
- `geometry/*.vtu` with a geometry report CSV;
- `sphere/` convergence (clean and jittered), abused-mesh, quadrature-sensitivity and test-weight CSVs, and VTK of the finest icosphere and the fully abused icosphere 3 (Cp, μ, leak, their exact values and the μ error);
- `ellipsoid/` the exact-solution checks, convergence, abused-mesh and quadrature-sensitivity CSVs, Cp along the cut y = 0 of the triaxial body, and VTK of the finest triaxial, prolate and thin oblate meshes and the fully abused triaxial one;
- `wing/` a summary CSV of every case, the spanwise loading per level (Γ, sectional lift by both routes), chordwise Cp at η = 0, 0.5, 0.8 and 0.95, and VTK of the body (Cp, μ) and the wake (μ_w, tip filaments);
- `examples/` the VTK of `abused_wing.jl`: body (Cp, Cp at the nodes, μ, quality), classified edges, wake (μ_w) and tip filaments;
- `comparison/` CSV tables and SVG plots of every section of the stage-6 study (refinement, refinement by direction, the two-dimensional limit, Cp and spanwise loading, skew, trailing-edge misalignment, zero incidence, robustness), and `provenance.txt` with the codes' git commits;
- `elements/` CSV tables and SVG plots of the stage-7 study (cube-sphere, refinement, sectional accuracy with each aspect ratio's values, Γ around midspan, midspan Cp, skew with both tip conditions, trailing-edge misalignment, zero incidence, abused meshes);
- `test_weights/` CSV tables and SVG plots of the stage-8 study (refinement, sectional accuracy, midspan Cp, abused meshes, zero incidence);
- `freewake/` CSV tables and SVG plots of the free-wake study (consistency, impulsive start, conversion, free wake and its sensitivity, MESA, abused meshes), and VTK of the final free wakes (free rows, particles, trailing-edge panel).

## Geometry conventions (used by the later steps)

- **Triangles and half-edges:** triangles are outward, with n ∝ (x₂-x₁)×(x₃-x₁). Half-edge `h = 3(t-1)+k` runs from `tris[t][k]` to `tris[t][mod1(k+1,3)]`.
- **Trailing-edge direction τ:** each trailing-edge chain runs along τ, chosen so that (wake_dir × τ)·up_ref ≥ 0. With the defaults (+x, +z), τ runs along +y on the reference wing.
- **Upper and lower sides:** the upper triangle owns the half-edge along τ, the same rule as FLOWPanel's `calc_shedding`. `kutta.corner_side[t][k]` labels every triangle corner at a trailing-edge node as Upper, Lower or Cap.
- **Tip openings:** `:cap` (default) adds ladder triangles on the existing nodes.
  - `:zip` adds no panels. It records the two sides as a `ChainPair` whose μ is tied in the constraint matrix C, and marks the body not closed.
  - `:open` is for debugging.
- **Pathology records:** T-junctions, gaps, non-matching trailing edges (with the union of their nodes along τ) and geometric touch pairs are all recorded for the DOF map and the touching quadrature tier.

## DOF map conventions

All topology decisions about the unknowns are made here (design doc §4 and §7).

- **Slots.** A slot is a nodal doublet value: the coefficient of one hat function.
  - Every body node has a slot.
  - A trailing-edge node whose corners lie on both an upper and a lower surface gets a second slot for its lower copy (`lower_slot[node]`).
  - `corner_slot[t][k]` is the slot of each triangle corner. The Galerkin system is assembled over slots: test functions are rows, and trial functions are columns, followed by one column per wake node.
- **Tip Kutta condition:** `DOFMapOptions(tip_kutta=:zero)` (default) ties each tip's lower trailing-edge copy to its upper one (constraint kind `:tip`). So μ_w = 0 at the tips, the Kutta condition at the tip corner: there is no tip filament in the wake and no jump at the cap seam. `:free` treats the trailing edge's end nodes like the others and leaves μ_w there to the solution (stage 6 compares the two).
- **Cap corners at the tip trailing-edge node** use the upper copy (`DOFMapOptions(cap_te_copy=:upper)`). This matters only with free tips: the jump μ_u − μ_l then runs as a filament along the cap's seam with the lower surface.
- **Constraints.** Each constrained slot is an affine combination of other slots, resolved to the unknowns q:
  - **T-junction hanging nodes** are linear along the long edge.
  - **Stitched gaps:** the nodes of the side with more nodes are linear along the coarser side, at their projection.
  - **Zipped openings:** the nodes of one side are linear along the other side, at equal normalised arclength. Trailing-edge nodes keep both copies, so zipping a tip does not tie the Kutta jump.
- **Row checks.** Rows of C sum to one, so a linear μ is reproduced exactly through T-junctions and coincident stitches. A gap of width w adds at most |∇μ|·w.
- **Kutta rows.** Each wake trailing-edge node k has a row giving μ_w,k = μ_u − μ_l.
  - On a non-matching trailing edge, these nodes are the union of both chains' nodes, and μ_u and μ_l are interpolated along their own chains.
  - Rows of K sum to zero, so a uniform μ sheds nothing.
  - `wake_rows(dm, te_of_wake_node)` gives the C rows of the wake nodes once the strips exist.
- **Reduced system.** It is (Cᵀ A [C; W]) q = Cᵀ b (`reduce_system`). Every body adds a bordered row: `bordered_rows` gives e and c = Cᵀ∫φ dA. A uniform μ on a body sheds no sheet and no wake, and the body carries no filaments, so the row is needed whether the body is closed, zipped or open.

## Kernel

`source_velocity(SourceTriangle(a, b, c), x; side)` gives w = (1/4π)∫_T (x−y)/|x−y|³ dA for a unit source spread over a flat triangle.

- **Closed form.** The tangential part is (1/4π)Σₑ νₑ ln((rᵢ+rⱼ+dₑ)/(rᵢ+rⱼ−dₑ)), with outward in-plane edge normals νₑ. Its denominator is formed without cancellation, so points 1e-8 diameters from an edge stay accurate. The normal part is Ω/4π, from the signed solid angle (Van Oosterom–Strackee), positive on the side n points to.
- **On the triangle's own plane:** the normal part is the one-sided limit chosen by `side`:
  - `+1` on the n side, `−1` on the other, `0` for the principal value;
  - this gives ±½ inside, ±¼ on an edge, ±α/4π at a vertex and 0 outside.
  - The tangential part is log-singular on the edges and at the vertices, so it must not be evaluated there.
- **Far field.** Beyond `FAR_DEFAULT` = 10 diameters the kernel switches to a monopole + quadrupole + octupole expansion. These use the triangle's exact moments, so the error is O((diam/r)⁴), and the jump at the switch is below 5e-6 relative.
- **Number types.** The kernel is generic: ForwardDiff works through it (field and vertex derivatives), and Float32 stays Float32. It runs in about 40 ns per evaluation (25 ns in the far field) with no allocations.
- **Validation.** Tests compare against two independent references: a semi-analytic polar quadrature that covers far, near, in-plane, edge- and vertex-adjacent points, slivers and 175° triangles, and adaptive cubature. They also check divergence- and curl-free fields, Gauss's law and the ½ jump.
- **AeroSingularityElements' `vel_srcpan`** (MESA's routine) is not reused. It agrees to ~1e-13 at regular points, but it zeroes an edge's log term within an absolute 1e-8 of the edge (100% wrong at 1e-9 diameters there), has no one-sided limits and no far field; `validation/kernels/compare_ase.jl` shows this.

## Galerkin system

Solving a problem takes six calls, and a lifting body adds its wake:

```julia
g    = build_geometry(mesh)                  # closed, classified triangulation
dm   = build_dofmap(g)                       # slots, C, Kutta rows
els  = build_elements(g)                     # flat triangles with hat functions
wake = build_wake(g, dm, V∞)                 # rigid strips from the trailing edge (lifting bodies)
σ    = source_strengths(els, V∞)             # σ = −V∞·n
A, b, stats = assemble(g, dm, els, σ, V∞; wake)            # A μ = b states ∫ψᵢ n·u dA = 0
sol  = solve_galerkin(dm, A, b; W=wake_rows(dm, wake.te))  # q, μ = C q, λ, cond; μ_w = dm.K * sol.q
```

`validation/common/wing_case.jl` wraps these with both force routes (`solve_wing`).

With `els = MixedMesh(g, dm)` instead of `build_elements(g)`, the input's quads stay whole as isoparametric quads (stage 7). Then `assemble(g, dm, els, V∞; wake)` takes no σ, because the quads integrate their sources with their sheets, and `surface_pressure(els, sol.μ, V∞)` gives Cp at the quads' Gauss points. `solve_wing(mesh; element=:quad)` does both.

**Formulation.**
- The doublet μ is piecewise linear, and each triangle carries the vortex sheet γ = n × ∇μ.
- A unit nodal value at corner b of triangle T gives γ_b = −e_b/(2A_T), where e_b is the edge opposite b, counterclockwise about the outward normal.
- The boundary condition is ∫ψᵢ n·(V∞ + u_σ + u_γ) dA = 0, with velocities taken on the fluid side.
- The test functions ψ are the hat functions by default (`test=:onering` or `:dual` select the alternatives). Stage 8 compares them: the dual cell is a sound alternative, while the one-ring is unstable on thin elements.
- The exact solution makes μ the outer perturbation potential.

**Quadrature tiers.** These apply on the test triangle; the source triangle's integral is the analytic kernel.
- *Far:* beyond `far_ratio = 3` diameters, a 3-point rule.
- *Near:* a 7-point rule on sub-triangles, split while they are within 1.5 of their size of the source triangle.
- *Touching (self, shared edge, shared vertex):* collapsed-coordinate rules graded toward the singular entity. These are hp-graded, with 6 levels of ratio 0.15 and the Gauss order rising from 3 by one per level, so the log-singular integrals converge exponentially.
  - *Wide corners:* the vertex rule's angular coordinate runs linearly along the opposite edge, so it under-samples the middle of a wide corner. Its error is 5e-9 at 90°, 2e-6 at 120° and 4e-4 at 150°.
  - A shared corner wider than 90° is therefore split by rays into equal angles of at most `corner_angle = 75°`, each piece carrying the vertex rule (3e-8 at 150°).
  - Stretched meshes, such as a scaled sphere with panels of aspect ratio 6, have many corners near 150°.
- *Contact:* triangles that touch beyond their shared nodes. This covers part of an edge at a T-junction, a hanging node, a gap, and a sliver's short edge, all within `contact_rel = 0.1` of the smaller diameter.
  - The test triangle is split at the source triangle's corners that lie on its edges. Each piece is graded toward its own contact: an edge when exactly two of its corners and the middle of the edge between them are on the source, otherwise a corner. A piece with contacts at several corners takes the self rule, graded toward all three edges.
  - Against adaptive cubature, a T-junction pair block is accurate to 3–5e-8; the shared corner's vertex rule alone gives 7–8e-3.
  - Across a slit of 1% or 0.1% of the edge length it is accurate to 3–4e-8. The near rule degrades as the slit narrows.
  - Near contacts are integrated as accurately as the regular tiers: on the ellipsoids, the column defects do not change with `contact_rel` from 0.1 to 0.5 (up to 34,000 contact pairs).
  - Clean meshes have no contact pairs, so their matrices are unchanged.

**Diagnostics** (`consistency(dm, A, b)` and `solve_galerkin`):
- **Row sums of A** vanish exactly.
- **Column sums of A·C** vanish for exact integration. An unknown's basis function is continuous through T-junctions and stitches, so its sheet has no net flux through a closed surface. Their relative defect therefore measures the quadrature. At a gap the flux through the slit adds O(ε log ε): about 3e-2 for slits of 1% of the edge length.
- **λ** absorbs the incompatibility of the right-hand side and shrinks as the quadrature is refined.
  - On icosphere 3 with jitter and T-junctions it is 6e-10 with the defaults and 8e-12 with the reference quadrature. A gap's slit flux shows up as about 1e-6.
  - On centrally symmetric meshes, such as the plain icosphere and cube-sphere, λ vanishes to rounding whatever the quadrature. The quadrature defects are even under x → −x while σ and V∞·n are odd, so they cancel.
  - λ is a single number and can be small by cancellation. The column-sum defect is the per-unknown check.
- **Leak:** the normal velocity at centroids, O(h) because facet normals are piecewise constant.

## Sphere in uniform flow (validation stage 3)

The exact solution is μ = V∞·x/2 and Cp = 1 − (9/4)sin²θ. The runs below use a generic V∞ direction, from `validation/sphere/run_sphere.jl`:

| mesh | nodes | μ L2 | μ max | Cp L2 (∇μ) | Cp L2 (sums) | leak rms | cond₁ |
|---|---|---|---|---|---|---|---|
| icosphere 4 | 2562 | 1.4e-4 | 4.9e-4 | 4.6e-3 | 1.9e-2 | 9.9e-3 | 2.6e3 |
| cube-sphere 24 (split quads) | 3458 | 1.4e-4 | 1.0e-3 | 1.2e-2 | 1.8e-2 | 9.9e-3 | 3.5e3 |
| icosphere 4, jittered | 2562 | 4.4e-4 | 2.0e-3 | 1.6e-2 | 2.3e-2 | 1.1e-2 | 2.6e3 |
| cube-sphere 24, jittered | 3458 | 3.9e-4 | 1.8e-3 | 1.6e-2 | 2.1e-2 | 1.0e-2 | 3.5e3 |

- **Orders of convergence:** μ converges at second order: 1.8–2.3 on clean meshes and 1.9–2.1 on jittered ones. Cp and leak converge at first order. The two velocity routes agree to O(h).
- **Quadrature error:** the default quadrature changes μ by 1–2e-5, on clean and abused meshes alike, against a discretization error (max) of 2e-3 to 2e-2. The far tier dominates that residual; the touching and contact rules are saturated.
- **Cost:** the 2,562-node assembly takes about 1 s on 16 threads.

**Abused meshes.** The design doc requires every stage to run on an abused mesh as well. The meshes are abused with `MeshTools.abuse` (seed 7), keeping moved nodes on the sphere. Results on icosphere 3 (642 nodes):

| abuse | μ L2 | μ max | Cp L2 (∇μ) | λ | column defect | Δμ from quadrature |
|---|---|---|---|---|---|---|
| clean | 5.8e-4 | 1.9e-3 | 9.8e-3 | 1e-19 | 3.5e-6 | 1.9e-5 |
| jitter 0.3 | 1.7e-3 | 6.8e-3 | 2.9e-2 | 5e-10 | 7.0e-6 | 1.2e-5 |
| 10 slivers at 5% | 2.2e-3 | 2.1e-2 | 5.5e-2 | 7e-10 | 4.3e-6 | 1.9e-5 |
| 6 gaps of 1% | 6.2e-4 | 1.9e-3 | 9.9e-3 | 6e-7 | 3.3e-2 (slit flux) | 1.9e-5 |
| 12 T-junctions | 5.8e-4 | 1.9e-3 | 1.0e-2 | 5e-11 | 4.7e-6 | 1.9e-5 |
| all at once | 2.9e-3 | 1.9e-2 | 6.1e-2 | 9e-7 | 2.2e-2 (slit flux) | 1.2e-5 |

- **Robustness:** every variant solves as one closed body, with cond₁ between 640 and 650. The cube-sphere 12 results are similar.
- **T-junctions and gaps** cost no accuracy. The constraints reproduce the linear μ, and the contact tier integrates the touching pairs.
- **Before the contact tier,** T-junction pairs got the shared corner's vertex rule. The column defect was then 2.3e-3, and refining the rules only brought it to 1e-3.
- **Jitter** raises the error constant about 3× but keeps the order.
- **Slivers** cost accuracy locally, for two reasons. Needle facets on a curved surface have tilted normals (6° against 0.6°), and short edges amplify nodal errors in ∇μ. The change of μ falls about 10× from a sliver's neighbours to two rings away. The local error converges under refinement, and the conditioning is unchanged.

## Ellipsoids in uniform flow (validation stage 4)

**Exact solution** (Lamb, *Hydrodynamics* §114), in `validation/common/exact_flows.jl`:
- The depolarization factors are αᵢ = abc∫₀^∞ dλ/((aᵢ²+λ)Δ(λ)). They come from Carlson's R_D, computed to rounding, and sum to 2.
- On the surface the perturbation potential stays linear: μ = Σ αᵢ/(2−αᵢ) V∞ᵢ xᵢ.
- The surface velocity is the tangential part of the constant vector k, with kᵢ = 2V∞ᵢ/(2−αᵢ). So Cp = 1 − |k − (k·n)n|²/|V∞|² follows the curvature.
- The tests check these formulas against the exterior potential, computed from the ellipsoidal coordinate. The two agree on the surface, and the exterior potential meets the Neumann condition and Laplace's equation to finite-difference accuracy and decays like a dipole.

**Bodies and meshes.** All bodies use the generic V∞ of the sphere study:
- a triaxial ellipsoid 2 : 1 : 0.5;
- a slender prolate spheroid 6 : 1 : 1;
- thin oblate spheroids 1 : 1 : 0.2 and 1 : 1 : 0.1.

The meshes are scaled unit spheres, so the prolate carries panels of aspect ratio up to 6 with corners near 150°, and the thinner oblate has a rim of curvature radius 0.01. From `validation/ellipsoid/run_ellipsoid.jl` (μ relative to its largest exact value):

| body, mesh | nodes | μ L2 | μ max | Cp L2 (∇μ) | Cp max (∇μ) | Cp L2 (sums) | leak rms | cond₁ |
|---|---|---|---|---|---|---|---|---|
| triaxial, icosphere 4 | 2562 | 1.4e-4 | 7.6e-4 | 2.2e-3 | 2.7e-2 | 8.7e-3 | 6.9e-3 | 2.7e3 |
| triaxial, icosphere 4 jittered | 2562 | 3.7e-4 | 3.6e-3 | 6.5e-3 | 1.4e-1 | 1.0e-2 | 7.3e-3 | 2.7e3 |
| triaxial, cube-sphere 24 | 3458 | 1.6e-4 | 2.4e-3 | 1.0e-2 | 7.0e-2 | 9.8e-3 | 7.8e-3 | 3.6e3 |
| triaxial, latitude-longitude 48×96 | 4514 | 1.9e-4 | 2.4e-3 | 9.5e-3 | 7.0e-2 | 9.7e-3 | 7.8e-3 | 4.8e3 |
| prolate 6:1:1, icosphere 4 | 2562 | 1.5e-4 | 1.2e-3 | 1.8e-3 | 3.9e-2 | 6.9e-3 | 5.9e-3 | 2.7e3 |
| oblate 1:1:0.2, icosphere 4 | 2562 | 1.4e-4 | 8.8e-4 | 2.9e-3 | 4.1e-2 | 1.4e-2 | 8.1e-3 | 2.8e3 |
| oblate 1:1:0.1, icosphere 4 | 2562 | 1.4e-4 | 1.5e-3 | 5.5e-3 | 9.3e-2 | 2.9e-2 | 9.7e-3 | 2.9e3 |

- **Orders of convergence:** μ is second order on every body and mesh family, 2.0–2.2 between the two finest levels, including the pole fans of the latitude-longitude spheres and the jittered meshes. Cp and leak are first order.
- **Cp along the cut y = 0** of the finest triaxial mesh matches the exact value at the centroids to 2.7e-2 (gradient) and 3.7e-2 (sums); `cp_cut_y0.csv` has the curves.
- **Quadrature error:** on icosphere 3 the default quadrature changes μ by 1.6–2.8e-5, against a discretization error (max) of 2.7–5.5e-3. The column defects are 2.6–6.8e-6, as on the sphere, and the reference quadrature reaches 4e-9 to 4e-8.
- **Wide-corner split:** without it (`corner_angle = 180`), the column defect of the prolate and the thin oblate is 10× larger (3.9e-5 and 6.5e-5), and even the reference quadrature stalls at about 2e-6.
- **Cost:** a 2,562-node assembly takes 1–1.5 s; 4,514 nodes take 2.8 s.

**Abused meshes** (triaxial body, icosphere 3, seed 7, moved nodes put back on the surface along the ray from the centre):

| abuse | μ L2 | μ max | Cp L2 (∇μ) | λ | column defect | Δμ from quadrature |
|---|---|---|---|---|---|---|
| clean | 5.8e-4 | 2.7e-3 | 4.6e-3 | 1e-19 | 4.8e-6 | 1.6e-5 |
| jitter 0.3 | 1.5e-3 | 1.3e-2 | 1.4e-2 | 1e-10 | 4.7e-6 | 1.7e-5 |
| 10 slivers at 5% | 2.0e-3 | 1.9e-2 | 1.7e-2 | 1e-10 | 3.8e-6 | 1.9e-5 |
| 6 gaps of 1% | 7.3e-4 | 3.7e-3 | 4.6e-3 | 4e-7 | 5.3e-2 (slit flux) | 1.6e-5 |
| 12 T-junctions | 5.8e-4 | 2.7e-3 | 4.7e-3 | 5e-11 | 4.8e-6 | 1.7e-5 |
| all at once | 2.8e-3 | 3.4e-2 | 2.4e-2 | 3e-7 | 2.3e-2 (slit flux) | 1.7e-5 |

- **Robustness:** every variant solves as one closed body, with cond₁ between 670 and 680. The cube-sphere 12 results are similar.
- **T-junctions** cost nothing, gaps 27% in μ L2, jitter 2.6× and slivers 3.5×.
- **All at once:** a sliver on the stretched, jittered mesh leaves a facet tilted 25° from the surface. The largest μ error grows 12× and the L2 error 5×.
- **Jitter bound:** jitter was first bounded by each node's shortest edge. On these stretched meshes that collapsed obtuse triangles to 176°, leaving facets tilted 31° with Cp errors of order one there. It is now bounded by the node's clearance; the sphere's jittered results above improved by about 15%.
- **Cp max error on jittered meshes** (0.14 at 2,562 nodes) comes from the gradient route on irregular triangles, as on the sphere. Its L2 converges at first order.

## Rigid wake and lifting bodies

**Wake strips** (`build_wake`, design doc §7):
- Every pair of neighbouring wake trailing-edge nodes spans a strip of two triangles per row. The rows sit at 0, ℓ₁, ℓ₁(1 + r), … from the trailing edge until 20 spans (`WakeOptions`):
  - the first row ℓ₁ defaults to the trailing-edge panels' size;
  - the growth ratio r defaults to 1.2;
  - the first row runs along the freestream or the trailing-edge bisector, and the rest along the freestream;
  - the strips end there (`:truncate`) or continue to 1000 times further (`:semi_infinite`).
- Every node of a strip carries the μ_w of its trailing-edge node (a row of K), so the wake sheds trailing vorticity only. The system has a column per wake node, folded back onto the unknowns by `wake_rows(dm, wake.te)`.
- **Orientation:** wake triangles face the upper side (wake_dir × τ). The trailing-edge filaments of the upper surface, the lower surface and the wake then cancel under μ_w = μ_u − μ_l, so no trailing-edge filament is modelled.
- **Tip filaments:** the strips' outer sides carry filaments of strength μ_w at the tip, the tip vortices, and nothing crosses the far end (no starting vortex). With the default zero-tip Kutta condition their strength is zero; they carry load only with `tip_kutta=:free`. With γ = n × ∇μ, a doublet layer is bounded by a filament running *clockwise* about its normal: where μ drops across a thin band, γ points against the counterclockwise tangent. So for a lifting wing the right tip filament runs downstream. The first implementation had them backwards, which gave the wing nearly two-dimensional loading.
- **Quadrature:** wake triangles share no mesh node with the body, so the contact tier finds the first row's trailing-edge contact geometrically. A filament starting at a test triangle's corner (the tip node) takes the vertex rule; its 1/r singularity there is integrable.

**No filaments on the body.** Wherever μ jumps on the body, or the layer ends at a free edge, the doublet layer equals the sheets plus a filament whose strength is the jump. On a panel that holds the filament along an edge, its swirl Γ/(2πd) is entirely normal to the panel. So the Galerkin row of a node on that edge has an integrand ~1/d, with d the distance from the edge. (Nodes off the edge are fine, since their hat functions vanish on it.) The body is therefore modelled by its sheets alone:
- **Free edges** (open or zipped tips): ∫₀ dd/d diverges, so no finite discrete equation exists. The filament has to be left out, and the tip vorticity then leaks.
- **Interior jump lines** (the cap seam below): the faces on either side see equal and opposite 1/d terms, because the swirl goes in through one face and out through the other. A seam node's row is then a finite principal value, but only when both sides are integrated together with matched singular rules. The pairwise assembly doesn't do that, so the seam filament is left out.
- **MESA doesn't run into this.** Morino's potential-based formulation sees a doublet panel through its solid angle, which stays bounded even at the panel's edges, and it collocates at centroids. So jumps in μ are harmless there, and MESA has them across every panel edge.

The consequences:
- **Cap seam:** with the default zero tip, μ_u = μ_l at the tip trailing-edge node, so there is no jump along the cap's seam. With `tip_kutta=:free` the cap's corners there take the upper copy, and the residual jump μ_u − μ_l along the seam with the lower surface is left out. It converges away:
  - choosing the lower copy instead changes CL by 0.4% at L1 and 0.2% at L2;
  - the lift at α = 0 is −5e-4 at L1 and −2e-4 at L2.
- **Zipped and open tips** leave their boundary vorticity out, so they are not physical models. Open tips lose the tip vortex (CL 0.47, tip loading 0.6), and zipped tips leave a large net flux through the tip (column defect 0.4–0.6). Caps are required, as the design doc's locked decisions say.

**Tip loading.** At a free tip the loading goes to zero, in inviscid theory and in real flow:
- A nonzero circulation at the tip would leave a concentrated vortex at the trailing edge's corner. The velocity there would be infinite, breaking the Kutta condition, and so would the induced drag.
- So Γ falls like √(distance from the tip); the elliptic Γ ∝ √(1 − (2y/b)²) shows the shape. The shed vorticity, dΓ/dy ~ 1/√d, concentrates near the tips and rolls up downstream into the tip vortex.
- In real flow the pressure still equalizes around the tip. The tip vortex sits over the upper surface there and often adds a local suction bump before the load drops.
- None of this holds where the tip isn't free: at a root against a wall, fuselage or symmetry plane, at an endplate, or at a winglet that carries the load on.

**The tip filament is a discretization artifact in both codes, of different kinds:**
- **Ours was avoidable, and the default now avoids it.** μ_w is linear between nodes and the tip node sits on the tip, so the model can go to zero there. With `tip_kutta=:free` the nodal Kutta condition leaves μ_u − μ_l free at the corner, and the solve settles on a nonzero value. Its share of the peak falls 0.24, 0.17, 0.13 at L1–L3, about ×0.73 per level (roughly h^½). The default `tip_kutta=:zero` sets μ_u = μ_l at the tip node, the Kutta condition at the corner, which removes both the tip filament and the cap-seam jump. The cost is that a linear ramp from zero under-represents the √d profile on the last strip.
- **MESA's is structural, not a mistake.** Its constant-strength strips form a vortex lattice, with a filament at every strip boundary. The tip filament is the last of these and carries the whole tip-strip strength. That strength represents the loading a finite distance in from the tip (~√(strip width)), and it shrinks as the strip narrows.
- Both converge with refinement. Our two tip conditions differ in CL by 1.8%, 0.7% and 0.3% at L1–L3, and reach the same lift (stage 6).

**Forces** (`Post`, in units of q∞; doc §9):
- **Pressure:** `pressure_force` integrates the gradient route's Cp, which is constant on each triangle, so the integral is exact.
- **Kutta–Joukowski:** `kutta_joukowski_force` gives F = (2/|V∞|²) V∞ × ∫μ_w dl along the trailing edge; the trapezoidal sum is exact for linear μ_w.
- **Trefftz drag:** `trefftz_drag` puts strips of mean μ_w between the trailing-edge nodes and evaluates their trailing vortices' downwash at the strip middles. It is the vector form of the validated `loads.jl` formula, with an elliptic loading giving CL²/(πAR). With free tips, its midpoint downwash keeps the drag finite despite a nonzero μ_w at the tips, which the exact Trefftz integral would make infinite.
- **Spanwise loading:** `spanwise_loading` gives the per-strip lift by both routes.

## Reference wing with a rigid wake (validation stage 5)

Results from `validation/wing/run_wing.jl`, on MESA's reference wing (NACA 0015, rectangular, AR 3.88) at α = 5°, with capped tips, the zero-tip Kutta condition (the default since stage 6) and a wake along the freestream to 20 spans. Stage 6 has the same tables for free tips, this stage's first default.

| level | unknowns | CL | CL (KJ) | CDi (pressure) | CDi (Trefftz) | Cm c/4 | cond₁ | time |
|---|---|---|---|---|---|---|---|---|
| L1 | 700 | 0.3383 | 0.3418 | 0.01072 | 0.00895 | +0.0012 | 1.1e4 | 0.6 s |
| L2 | 2693 | 0.3389 | 0.3415 | 0.01011 | 0.00927 | +0.0007 | 4.4e4 | 2.3 s |
| L3 | 10555 | 0.3395 | 0.3411 | 0.00981 | 0.00943 | +0.0002 | 1.7e5 | 18 s |
| Richardson | | — | — | 0.00953 | 0.00960 | | | |

- **Lift converges, but uniform refinement can't extrapolate it.** CL rises slowly (0.3383, 0.3389, 0.3395) while CL_KJ falls (0.3418, 0.3415, 0.3411). The chordwise error (first order, from above) and the tip region's (from below) nearly cancel at these levels, so the observed order is about 0. Refining the two directions separately (stage 6) gives 0.3407 by the circulation route, the same as with free tips.
- **Induced drag** extrapolates to 0.00953 (pressure) and 0.00960 (Trefftz), 0.7% apart.
- **Strip by strip,** the two lift routes agree within 1% at midspan. At the tip strip they differ by design, since the cap's sharp edges carry inviscid side-edge suction.
- **Tip loading** is zero at the tip node by construction. With free tips its share of the peak falls 0.24, 0.17, 0.13, as the doc asks to check.
- **Kutta condition:** the upper and lower Cp of the midspan trailing-edge triangles differ by 0.027 at L1 and 0.019 at L2.
- **Wake geometry doesn't matter at this resolution:** the bisector first row, 50 spans, semi-infinite strips, MESA's first row (0.045) and growth 1.1 or 1.5 all change CL by at most 8e-5 at L2.
- **Quadrature** is converged: the reference rules change CL by 1e-5 at L1. At L3 the column defect reaches 2e-3 at the tip edges, where the cosine spacing makes the tip strip 0.0014 m wide against the perpendicular cap. Widening the contact tolerance to 0.25 cuts that 9×, at twice the assembly cost, and changes CL by 3e-6, so the default stays.
- **Zero incidence:** lift is 5e-6 at L1 and −6e-6 at L2. With free tips it is −5e-4 and −2e-4, from the caps' up-down asymmetry. The pressure drag, a discretization error, is 1.2e-3 and 3.8e-4.
- **Conditioning** grows 4× per level and is about 13× MESA's (Morino's second-kind Dirichlet formulation conditions better). A dense LU is unaffected at these sizes.

**Abused meshes** (the MESA smoke set, seed 1):
- At L1 every case changes CL by at most 0.91%. The largest is opposite spanwise lean of 0.3; misaligned trailing edges give 0.59%, independent jitter 0.47%, jitter 0.30%, and T-junctions, gaps and slivers under 0.08%.
- At L2 every case is within 0.36%.
- Non-matching trailing edges shed the wake from the union of both edges' nodes and change CL by under 0.11% (L1) and 0.02% (L2).
- These are larger than with free tips (0.64% and 0.26%). The skew patterns flip quad diagonals next to the trailing edge, and the zero tip responds more to those flips (stage 6).

**Next to MESA** (stage 6 explains the gaps):
- At L3 with capped tips MESA gives CL 0.3298 and CL_KJ 0.3343. That is 2.9% below our pressure lift, and its two routes still differ by 1.4% with the KJ value rising under refinement.
- MESA's Cm about c/4 is +0.005 and ours +0.0002.
- Its CDi (0.0095 pressure, 0.0091 Trefftz) is close to ours.

## Comparison with the Morino code (validation stage 6)

`validation/comparison/run_comparison.jl` solves the same meshes with both codes: the reference wing at α = 5°, the tips closed with our caps (MESA gets them through `add_caps`), and a rigid wake along the freestream to 20 spans. An exact two-dimensional referee judges the section. `validation/common/panel2d.jl` is a linear-vortex panel method (Katz & Plotkin §11.4), second order against the exact Joukowski lift. It gives the NACA 0015 section at 5° cl = 0.61701 and cm about c/4 = −0.00917.

**Refinement**, with (nspan, nchord) = (17, 19)·2^(level−1). "Galerkin" runs the defaults, with the zero-tip Kutta condition; "free tip" is `DOFMapOptions(tip_kutta=:free)`, the default until this stage:

| level | solver | CL | CL (KJ) | CDi (pressure) | CDi (Trefftz) | Cm c/4 |
|---|---|---|---|---|---|---|
| L1 | Galerkin | 0.3383 | 0.3418 | 0.01072 | 0.00895 | +0.0012 |
| | Galerkin, free tip | 0.3442 | 0.3480 | 0.01087 | 0.00942 | −0.0010 |
| | MESA | 0.3287 | 0.3320 | 0.01009 | 0.00853 | +0.0059 |
| L2 | Galerkin | 0.3389 | 0.3415 | 0.01011 | 0.00927 | +0.0007 |
| | Galerkin, free tip | 0.3413 | 0.3440 | 0.01016 | 0.00953 | −0.0002 |
| | MESA | 0.3284 | 0.3327 | 0.00979 | 0.00885 | +0.0058 |
| L3 | Galerkin | 0.3395 | 0.3411 | 0.00981 | 0.00943 | +0.0002 |
| | Galerkin, free tip | 0.3405 | 0.3422 | 0.00983 | 0.00959 | −0.0002 |
| | MESA | 0.3298 | 0.3343 | 0.00948 | 0.00910 | +0.0053 |

Uniform refinement extrapolates only the free tip's lift, to CL 0.3402 (p = 1.9). Ours with the zero tip rises too slowly to extrapolate (stage 5). MESA's lift is not monotone, because its chordwise and spanwise errors have opposite signs (below), so Richardson extrapolation fails for it too.

**The section explains both gaps.** Take the midspan section of wings of aspect ratio 15.5, 31 and 62 (68 strips each) and extrapolate it quadratically in 1/AR to the infinite wing. Its errors against the exact section are:

| panels per surface | Galerkin cl | Galerkin cm | Galerkin cl, one diagonal | Galerkin cm, one diagonal | MESA cl | MESA cm |
|---|---|---|---|---|---|---|
| 19 | +1.28% | −0.0022 | +1.57% | −0.0026 | −5.47% | +0.0055 |
| 38 | +0.66% | −0.0010 | +0.78% | −0.0011 | −3.75% | +0.0045 |
| 76 | +0.35% | −0.0004 | +0.40% | −0.0005 | −2.53% | +0.0033 |
| observed order | 0.95 | 1.2 | 1.0 | 1.1–1.2 | 0.55 | 0.3–0.4 |

- **The "one diagonal" columns are our method's own error.** The default split mirrors its diagonals at midspan, exactly where the section is measured, and that dents the circulation there (found in stage 7). Splitting every quad along the same diagonal (`GeometryOptions(split=:d13)`) removes the dent, and the section then matches the isoparametric quads' to 1e-5. The default split's figures understate the error by about a fifth.
- **The codes bracket the exact section.** Ours converges at first order from above, MESA at order 0.5 from below.
- **Ours is 6× closer in cl and in cm at 76 panels per surface** (one diagonal: 6.4× and 6.3×). MESA's cl error there is what ours has at about 12 panels.
- **For scale,** the 2-D method itself on PAGE's chordwise spacing is off by −1.02%, −0.48% and −0.22%.
- **The deficit is in MESA's circulation.** Its Kutta–Joukowski lift has the same deficit as its pressure lift, so the problem is not in its pressure post-processing.
- **The loading shows it too.** MESA's circulation loading at L2 is 3.6% below ours at midspan and 1.1% below at η = 0.8. At η = 0.95 it is 4% above ours, where our zero tip pulls the loading down.

**Refinement by direction.** Starting from the L2 mesh, the chord is refined at 34 strips (19 to 152 panels per surface) and the span at 38 panels (17 to 136 strips). Richardson extrapolation uses the last three meshes, and only where the observed order is at least 0.3:

| quantity | solver | chordwise limit (order) | spanwise limit (order) | both directions |
|---|---|---|---|---|
| CL (KJ) | Galerkin | 0.3406 (1.1) | 0.3416 (3.8) | **0.3407** |
| | Galerkin, free tip | 0.3407 (1.2) | 0.3441 (1.1) | **0.3407** |
| | MESA | 0.3446 (0.54) | 0.3296 (1.5) | **0.3415** |
| CL | Galerkin | 0.3381 (1.2) | creeps (below) | |
| | Galerkin, free tip | 0.3381 (1.2) | creeps | |
| | MESA | 0.3403 (0.46) | 0.3255 (1.4) | 0.3374 |
| Cm c/4 | Galerkin | +0.0010 (1.5) | creeps | |
| | Galerkin, free tip | +0.0010 (1.2) | creeps | |
| | MESA | +0.0043 at 152 panels (order 0.25) | +0.0064 (1.1) | |

"Both directions" assumes the errors of the two directions add: Q(∞, ∞) ≈ Q(34, ∞) + Q(∞, 38) − Q(34, 38).

- **Both codes converge to nearly the same lift.** By the circulation route, extrapolated in both directions, both our tip conditions give 0.3407 and MESA gives 0.3415, 0.24% apart. MESA's slow chordwise convergence makes its value the less certain one.
- **The gap at the refinement levels is chordwise, mostly MESA's.** Refining the chord from 19 to 152 panels raises MESA's lift by 3.0% and lowers ours by 0.55% (2.0% with free tips). MESA's Cm falls from +0.0065 to +0.0043 toward ours (+0.0010 in the limit), too slowly to extrapolate.
- **Why MESA's sequence isn't monotone:** refining its span from 17 to 136 strips lowers its lift by 2.1%, against its chordwise gain.
- **Cm at L3:** the codes differ by 0.0052 on the wing; their sectional cm errors at 76 panels differ by 0.0038. In span, our Cm falls 0.0010, 0.0008 and 0.0007 per doubling of the strips, while MESA's rises by 0.0006, 0.0004 and 0.0002.
- **Our pressure lift creeps under spanwise refinement:** +0.0017, +0.0011 and +0.0009 per doubling of the strips, while the circulation lift has converged to 1e-4 from 34 strips on. Free tips creep the same way.
  - All of the creep sits in the outer 5% of the span. Pressure lift minus circulation lift there grows from +0.0023 at 17 strips to +0.0051 at 136, as much as the overall difference shrinks.
  - The cause is the cross-flow around the sharp edge between the wing surface and the flat cap. Its suction is an integrable singularity, which the cosine spanwise spacing keeps resolving. The pressure lift therefore converges, but slowly. The extra lift acts near the trailing edge, where the tip vortex forms.
  - For wings with sharp tips the circulation route is the reliable one. MESA's pressure lift converges in span (order 1.4).
- **Induced drag:**
  - Our Trefftz drag converges in both directions, to 0.00960.
  - With free tips it grows as the strips near the tip narrow: 0.0091 to 0.0101 from 17 to 136 strips. With a nonzero μ_w at the tip, the exact Trefftz integral diverges logarithmically, and the discrete one follows it. MESA's tip strip carries the same kind of filament, and its Trefftz drag creeps too.
  - The pressure drags reach 0.00956 (ours) and 0.00907 (MESA).

**Chordwise Cp** at L2 (`cp_L2_eta*.svg`):
- The codes agree within 0.018 between 10% and 90% chord at all four stations.
- Suction peaks agree within 0.025 (−1.223 against MESA's −1.214 at midspan).
- The largest difference is on the last panel: ours 0.24 upper and 0.26 lower, against MESA's 0.18.

**Skew sensitivity** uses the patterns of MESA's dev study, though both codes run the stock rigid wake here; the dev study used `solve_decoupled`. Changes are from the clean mesh of the same level, over 4 patterns × 4 amplitudes (σ = 0.05–0.3 of the station spacing; for jitter, the mean |Δ| over 5 seeds):

| level | mean (largest) \|ΔCL\|: Galerkin | Galerkin, clean split | MESA | mean \|ΔCDi\|, pressure: Galerkin | MESA | mean \|ΔCm\|: Galerkin | MESA |
|---|---|---|---|---|---|---|---|
| L1 | 0.45% (0.92%) | 0.062% (0.22%) | 0.24% (0.65%) | 0.91% | 9.6% | 0.00056 | 0.00018 |
| L2 | 0.18% (0.36%) | 0.035% (0.13%) | 0.19% (0.94%) | 0.39% | 2.5% | 0.00022 | 0.00024 |

- **Lift, at L2:** ours is as sensitive as MESA's on average (0.18% against 0.19%) and 2.6× less in the worst case. With the clean split it is 5.3× less on average and 7.5× in the worst case.
- **Lift, at L1:** ours is 1.8× more sensitive on average, but 3.9× less with the clean split.
- **Lift, with refinement:** our sensitivity falls 2.4× per level, MESA's 1.3×. Independent jitter is MESA's weak spot: 0.94% at L2 against our 0.19%.
- **Drag:** our pressure drag is 6–11× less sensitive. Opposite lean at L1 changes MESA's by up to 42%.
- **Most of our lift sensitivity is triangulation.**
  - The patterns move the rows next to the trailing edge (and next to the leading edge), and the shortest-diagonal split then flips some of their diagonals: 34–36 quads at L1.
  - Keeping the clean mesh's diagonals (`GeometryOptions(split=:pattern, split_pattern=...)`, the "clean split" column) isolates the skew itself: 0.062% at L1 and 0.035% at L2.
  - The zero tip makes the flips cost about twice what they cost with free tips: 0.45% against 0.24% at L1, and 0.18% against 0.09% at L2, while the clean-split sensitivity is the same. The flips in the interior strips cause it; the few in the tip strips partly offset it.
  - Flipping the trailing-edge row's diagonals everywhere on the clean mesh changes CL by −0.26%, −0.10% and −0.03% at L1–L3 (free tips: +0.41%, +0.16%, +0.08%). That is an ordinary discretization difference that converges away.
- **So the design doc's "much less sensitive" holds with a fixed triangulation (5× at L2), for drag, for independent jitter and in the worst case.** With the default split, the diagonal flips alone match MESA's average skew sensitivity at L2. A split that doesn't flip under small perturbations would remove them. The isoparametric quads of stage 7 have nothing to flip and are as insensitive as the clean split (0.070% at L1, 0.030% at L2), with either tip condition.

**Trailing-edge misalignment.** The rows aft of 80% chord are offset by ±δ/2 with a shared trailing edge, or the upper rows by δ with a non-matching one.
- **Shared trailing edge:** MESA is the less sensitive here. Its lift changes by −0.11% (L1) and −0.07% (L2) at δ = 0.3, and by about 0 at 0.5.
  - Ours changes by −0.59% and −0.84% at L1, and by −0.24% and −0.33% at L2.
  - In the worst case that is 8× MESA's change at L1 and 4.7× at L2. Ours falls 2.5× per level, MESA's 1.5×. MESA's change passes through zero near δ = 0.5, so a ratio at a single δ there means little.
  - The cause is again diagonal flips. The shift skews two rows, the first row of the shifted band and the trailing-edge row, and the default split flips one diagonal per strip and surface there: 34 of 646 quads at L1, 68 of 2584 at L2.
  - With the clean split ours changes by at most 0.26% at L1 and 0.09% at L2. At L2 that is at parity with MESA: smaller at δ ≤ 0.3, larger at 0.5.
- **Non-matching trailing edge:** MESA cannot represent one, because its trailing-edge split needs matching facets. It changes our lift by at most 0.11% at L1 and 0.04% at L2.

**Zero incidence (L2):**
- **Ours (zero tip):** −6e-6 on the clean mesh, and at most 8e-5 on the skewed and non-matching ones.
- **Free tips:** −2e-4 to −4e-4 on every mesh, from the caps' up-down asymmetry.
- **MESA:** 1e-7 on the clean and leaning meshes, but 4.4e-3 under independent jitter, 1.3% of its lift at 5°.

**Robustness**, on MESA's smoke set plus all abuses at once, capped:
- **Our solver** runs every mesh, within 0.91% (L1) and 0.36% (L2) of the clean lift.
- **MESA fails on the T-junction mesh,** with a singular least-squares system in its surface gradient (`calc_CHTLS_pan!`). It does run the combined mesh, which also has T-junctions.
- **MESA cannot represent the non-matching trailing edges.**
- **Jitter of 0.2** moves MESA's lift by 4.5% (L1) and 1.6% (L2), against our 0.30% and 0.11%.

**Cost and conditioning at matched resolution.** At L3 we have 10,555 unknowns against MESA's 10,636 panels, on 16 threads:

| solver | assembly, L1 / L2 / L3 | factorization and solve | cond₁ | cond₁, rows scaled by the diagonal |
|---|---|---|---|---|
| Galerkin | 0.53 / 2.1 / 13.6 s | 0.02 / 0.23 / 5.5 s | 1.1e4 / 4.4e4 / 1.7e5 | 8.4e3 / 2.8e4 / 1.1e5 |
| MESA | 0.01 / 0.14 / 3.2 s | 0.01 / 0.09 / 3.1 s | 8.1e2 / 3.0e3 / 1.3e4 | the same |

- **Assembly:** ours costs 4.2× MESA's at L3 and 15× at L2. The Galerkin double integrals are the cost, and the near field shrinks relative to the far field as the mesh refines.
  - The cost also grows with panel aspect ratio: the 76-panel section at AR 62 takes 154 s, against MESA's 6 s.
- **Conditioning:** both condition numbers grow 4× per level. Ours is 13× MESA's, and still 8.5× with the rows scaled by the diagonal, so it is not a scaling artifact.
  - Our hypersingular first-kind operator conditions worse than Morino's second-kind one.
  - A dense LU doesn't care at these sizes. An iterative solver would need a preconditioner.
- **At equal sectional accuracy,** MESA's cl error at 76 panels per surface is ours at about 12.

**The tip Kutta condition** (`DOFMapOptions(tip_kutta)`): the default `:zero` sets μ_u = μ_l at the tip node, and `:free` leaves μ_w there to the solution. This stage made `:zero` the default.
- **Both converge to the same lift.** The chordwise limits agree to 1e-4: 0.3381 for both in pressure, and 0.3406 against 0.3407 by Kutta–Joukowski. The circulation lift extrapolated in both directions is 0.3407 for both.
- **The zero tip approaches from below,** the free tip from above. Their difference shrinks about 2.4× per level: 0.0060, 0.0024, 0.0010.
- **The zero tip removes three defects of the free one:**
  - the up-down asymmetry at α = 0 (6e-6 against −2e-4 at L2);
  - the tip filament;
  - the spanwise creep of the Trefftz drag.
- **Its cost:** the diagonal flips the default split makes next to the trailing edge cost about twice as much lift (see the skew sensitivity above). Without flips the cost is gone: the isoparametric quads of stage 7 are equally sensitive with either tip condition.
- **At a fixed chordwise resolution the two never meet.** However fine the strips, they keep different circulation lifts: 0.3416 (zero) and 0.3441 (free) at 38 panels. The tip treatment acts over a region set by the chordwise spacing.

## Isoparametric quads (validation stage 7)

The design doc's second element type keeps the input's quadrilaterals whole. `IsoQuad` maps the parameter square bilinearly, x(ξ, η) = Σ N_a x_a, and carries μ on the same bilinear functions. `validation/elements/run_elements.jl` compares it with the split triangles, the default, on the cube-sphere and on the reference wing. Both elements use the same nodes and DOF map, so they have the same unknowns.

**The element.**
- The edges are straight, so neighbouring quads, and quads next to triangles, meet exactly. A quad whose corners are not coplanar is a hyperbolic paraboloid.
- The sheet's density per unit parameter area is a polynomial: γ dA = (μ_ξ a₂ − μ_η a₁) dξ dη, with the tangents a₁ = ∂x/∂ξ and a₂ = ∂x/∂η. Nothing is divided by the Jacobian.
- `MixedMesh(g, dm)` turns every welded input quad into an `IsoQuad`. Its corner slots come from its two triangles, so trailing-edge copies, constraints and Kutta rows carry over unchanged. Caps, input triangles and the wake stay flat triangles.

**Inner integrals** (`quad_field`). A curved sheet of varying strength has no closed form, so they are numerical:
- *Far* (at least a diameter away): tensor Gauss, with 3 to 8 points per direction as the distance shrinks.
- *Near:* polar coordinates about the closest point P, with sinh substitutions in angle and radius. The parameters are first mapped by the Cholesky factor of the metric at P, so that distances and angles are physical to first order. Without that map, a sliver's sectors crowd its whole angular variation into a few degrees.
- *On the element* (self pairs): principal values, with each ray's 1/ρ term subtracted and added back analytically (Guiggiani). `side = ±1` gives the one-sided limits, adding half the sheet's jump γ × n and the source's σn.
- *Closest point* (`project`): Newton from the nearest node of a 5×5 parameter grid, then the four edges. Each step is held inside [−½, 3/2]² and halved until the distance falls. Starting from the centre alone misses on trapezoids. A step that leaves the square must be cut short, not abandoned: abandoning it returned an edge point 5 mm away on sheared leading-edge slivers, where Cp then reached −33.
- Against adaptive cubature they agree to 1e-7–1e-10, from far away down to 1e-6 diameters. The test shapes are squares, trapezoids, warped quads, slivers and a sheared sliver with 5.7° corners. A near evaluation costs 15–21 μs, against 40 ns for the triangle's closed form.

**Outer rules** (`SquareTiers`): the triangle tiers carried over to the parameter square.
- *Tiers:* tensor Gauss far (2×2) and near (3×3, refined by a quadtree near the source), and rules graded toward a shared side, a shared corner (two triangles from it) and all four sides (the self rule: four triangles from the centre). The contact rule for T-junctions, gaps and the wake cuts the square at the source's corners, as on triangles.
- *Wide corners:* a shear makes a quad's corner wider than its 90° in parameters, which the vertex rule's angular coordinate under-samples. Such a corner is split by physical angle into pieces of at most `corner_angle` (`vertex_rect!`). On a 150° corner with unequal sides the stored rule is off by 4e-5, and more levels don't help; the split rule is off by 1e-7.
- *Levels:* `AssemblyOptions(quad_touch_levels=4)` grades over 4 levels where the triangles use 6. On the cube-sphere, 6 levels change μ by only 2.5e-6, at 2.6× the cost.
- *Pair blocks against adaptive cubature:* 2.7e-6 for a shared side (the neighbour folded 25° and warped), 1.9e-7 for a shared 150° corner (1e-4 with the stored rule) and 3.8e-7 for a warped quad's self block.

**The two-dimensional check.** Split every quad of a prismatic wing along the same diagonal. For a spanwise-uniform μ, each triangle then carries the same constant sheet as its bilinear quad. The triangles' hat functions, summed along the span, give the quads' chordwise test functions. So in the two-dimensional limit the two elements give the same equations, by entirely different inner integrals and rules.
- At the midspan of wings of aspect ratio 15.5–62 (68 strips) they agree to 2e-5 in cl at every aspect ratio.
- The sections extrapolated to the infinite wing agree to 6e-5 (19 panels per surface) and 5e-6 (38).

**The default split dents the circulation at its mirror line.** Its tie-break mirrors the diagonals at y = 0, so the quads on either side of midspan are split opposite ways. The triangles' test functions there differ from a uniform split's.
- At the midspan trailing-edge node of the aspect-ratio-62 wing (19 panels per surface), Γ is 0.26% below the quads'. One station away it is 0.07% below, and two stations away under 0.02%. One diagonal throughout stays within 0.004% of the quads.
- The dent grows with the strips' width and shrinks as the chord is refined: 0.30%, 0.44% and 0.60% at aspect ratios 15.5, 31 and 62 with 9 panels per surface.
- The sectional studies of stages 6 and 8 measured exactly there; their text now accounts for it.

**Cube-sphere** (`cube_sphere(24)`, 3,458 nodes, with the generic V∞ of stage 3). μ is relative to its largest exact value. Cp is measured at the element centres (the triangles' centroids and the quads' parameter centres) and in L2 over the elements:

| mesh | element | μ L2 | μ max | Cp L2 | Cp max | Cp L2 over the elements | column defect | assembly |
|---|---|---|---|---|---|---|---|---|
| clean | triangles | 1.4e-4 | 1.0e-3 | 1.2e-2 | 3.0e-2 | 2.6e-2 | 4.4e-6 | 1.3 s |
| | quads | 8.9e-5 | 2.8e-4 | 5.7e-4 | 2.4e-3 | 3.0e-2 | 1.8e-6 | 24 s |
| jitter 0.3 | triangles | 3.9e-4 | 1.8e-3 | 1.6e-2 | 0.10 | 2.8e-2 | 8.0e-6 | 1.3 s |
| | quads | 3.8e-4 | 1.6e-3 | 1.1e-2 | 0.079 | 3.4e-2 | 4.4e-6 | 26 s |

- **μ is second order for both.** In L2 the orders are 2.0–2.1 between n = 12 and 24 on clean meshes and 1.9–2.0 jittered. On clean meshes the quads' error is 0.64× the triangles' in L2 and 0.27× at its largest. Jittered, the two are within 4–8% of each other.
- **Cp at the quads' centres is second order on clean meshes,** 20× below the triangles' at their centroids.
  - A nearly flat quad's normal at its centre matches the sphere's to second order. Away from the centre it is off by O(h), as a triangle's is away from its centroid.
  - Over the elements (the quads' 3×3 Gauss points, 7 points per triangle), both are first order, and the quads' error is 15% higher.
  - Jittered quads are not symmetric about their centres, so their Cp there is first order too, 0.7× the triangles'.
- **The quadrature is consistent:** column defects stay below 1e-5, and λ on the jittered meshes is at most 1.4e-9 (3e-11 at n = 24).

**Wing refinement** (clean meshes, α = 5°, the defaults otherwise):

| level | element | CL | CL (KJ) | CDi (pressure) | CDi (Trefftz) | Cm c/4 | column defect | cond₁ | assembly |
|---|---|---|---|---|---|---|---|---|---|
| L1 | triangles | 0.3383 | 0.3418 | 0.01072 | 0.00895 | +0.0012 | 3.6e-5 | 1.1e4 | 0.4 s |
| | quads | 0.3384 | 0.3419 | 0.01079 | 0.00895 | +0.0012 | 1.4e-5 | 8.9e3 | 7.7 s |
| L2 | triangles | 0.3389 | 0.3415 | 0.01011 | 0.00927 | +0.0007 | 6.6e-5 | 4.4e4 | 2.1 s |
| | quads | 0.3388 | 0.3414 | 0.01013 | 0.00926 | +0.0008 | 1.7e-5 | 3.5e4 | 30 s |
| L3 | triangles | 0.3395 | 0.3411 | 0.00981 | 0.00943 | +0.0002 | 2.3e-3 | 1.7e5 | 13 s |
| | quads | 0.3394 | 0.3411 | 0.00981 | 0.00943 | +0.0002 | 7.3e-4 | 1.4e5 | 121 s |

- **The elements agree on the wing.**
  - Lift differs by 0.03%, 0.02% and 0.01% at L1–L3, by either route.
  - The induced drags extrapolate to the same values: 0.00953 against 0.00952 (pressure), and 0.00960 against 0.00961 (Trefftz).
  - Only the pressure drag at L1 differs by more (+0.7%), and the gap shrinks to −0.03% at L3.
- **Midspan Cp** (L2) agrees within 0.001 between 10% and 90% chord. The suction peaks are −1.223 and −1.224. The trailing-edge panels are 0.240/0.261 (upper/lower) for the triangles against 0.227/0.250 for the quads.
- **The quads' column defect is lower:** 3× lower at L3, where the triangles' tip-edge pairs dominate (stage 5).

**Sectional accuracy** (as in stage 6: the midspan section of wings of aspect ratio 15.5, 31 and 62 with 68 strips, extrapolated in 1/AR):

| panels per surface | triangles cl | triangles cm | one diagonal cl | one diagonal cm | quads cl | quads cm |
|---|---|---|---|---|---|---|
| 19 | +1.28% | −0.0022 | +1.57% | −0.0026 | +1.57% | −0.0026 |
| 38 | +0.66% | −0.0010 | +0.78% | −0.0011 | +0.78% | −0.0011 |
| observed order | 0.95 | 1.2 | 1.01 | 1.17 | 1.01 | 1.18 |

- **The quads match the one-diagonal triangles,** as the two-dimensional check requires.
- **The default split's figures are lower only because they are measured on its mirror line.**
- **One diagonal throughout raises the triangles' column defect** to 8e-4 at the trailing edge, against the default split's 5e-5 (at the tips). The section is unaffected.
- **Stage 6's comparison with MESA is corrected:** with one diagonal at 76 panels per surface the error is +0.40% (cm −0.0005). This was measured separately with the triangles only; the quads would take about 40 minutes there.

**Skew sensitivity** (stage 6's patterns and amplitudes). Changes are from the clean mesh of the same level and element; for jitter, the mean over 5 seeds:

| level | mean (largest) \|ΔCL\|: triangles | triangles, clean split | quads | MESA (stage 6) | mean \|ΔCDi\|, pressure: triangles | quads | mean \|ΔCm\|: triangles | quads |
|---|---|---|---|---|---|---|---|---|
| L1 | 0.45% (0.92%) | 0.062% (0.22%) | 0.070% (0.33%) | 0.24% (0.65%) | 0.91% | 0.11% | 0.00056 | 0.00009 |
| L2 | 0.18% (0.36%) | 0.035% (0.13%) | 0.030% (0.14%) | 0.19% (0.94%) | 0.39% | 0.09% | 0.00022 | 0.00004 |

- **The quads are as insensitive as the triangles that keep the clean mesh's diagonals.**
  - Against the default split they are 6× less sensitive in lift, 4–8× in drag and 6× in moment.
  - Against MESA they are 3.4× less sensitive on average at L1 and 6.3× at L2, and 6.6× in the worst case at L2.
  - With the quads, the design doc's "much less sensitive" holds without fixing the triangulation.
- **Their response is smooth.** It grows with the square of the skew, as an ordinary change of discretization does: for lean at L1, 0.009%, 0.036%, 0.14% and 0.32% at σ = 0.05–0.3. The default split's flips already jump at the smallest skew: −0.26% at σ = 0.05.
- **The tip condition no longer matters.**
  - With free tips the quads' sensitivity is the same as with the zero tip: 0.070% (0.32%) at L1 and 0.030% (0.14%) at L2.
  - The default split's is 0.24% at L1 and 0.087% at L2 with free tips, half its zero-tip figures, as stage 6 found.
  - So the zero tip's extra sensitivity came entirely from the diagonal flips.

**Trailing-edge misalignment** (stage 6's offsets of the rows aft of 80% chord, δ = ±0.1 to ±0.5 of the station spacing):
- **Shared trailing edge:**
  - The quads change lift by at most 0.23% at L1 and 0.09% at L2. That matches the clean split (0.26% and 0.09%), against the default split's 0.84% and 0.33%.
  - Their change is symmetric in ±δ and quadratic.
  - At δ = 0.3 it is 0.08% at L1 and 0.03% at L2, below MESA's 0.11% and 0.07% (stage 6).
- **Non-matching trailing edge:** at most 0.06% at L1 and 0.02% at L2, against the triangles' 0.11% and 0.04%.

**Zero incidence (L2):** the quads' lift is under 1.5e-6 on the clean, misaligned and leaning meshes, against the triangles' 8e-6. Under independent jitter it is −4e-5, against the triangles' +8e-5.

**Abused meshes** (MESA's smoke set, plus all abuses at once, seed 1): every mesh runs.
- Mean (largest) |ΔCL| from the clean mesh is 0.096% (0.33%) at L1 and 0.032% (0.14%) at L2. The triangles' is 0.26% (0.91%) and 0.10% (0.36%).
- The quads' largest changes come from opposite lean of 0.3 (+0.33% and +0.14%) and jitter of 0.2 (+0.25% at L1, +0.07% at L2).
- T-junctions, gaps, slivers, sliver stations and non-matching trailing edges change their lift by under 0.06%.
- On the gap meshes the column defect is the slit flux, as in stage 3: 4.5e-2 for the quads, 3.2e-2 for the triangles.

**Distorted quads integrate less accurately.** The quads' column defect at L1 and L2 reaches:
- 5.5e-3 and 4.6e-3 with slivers. The moved nodes leave quads with one side 2% of the opposite one, so their corner Jacobians differ 50-fold.
- 2.4e-3 and 2.8e-3 with independent jitter of 0.3, which shears corners to 5° and 175°.

The triangles' defects on these meshes are 2e-5–5e-5 and 5e-4. The self blocks of such quads converge slowly as the rules are refined, to about 1e-3. Sliver stations (quads of aspect ratio 50 and more) give 1.5e-3 at L1 and 7e-4 at L2, about the triangles' 8e-4.

Lift is unaffected: on these meshes the quads change it less than the triangles do. Splitting quads whose corner Jacobians differ by more than some factor would remove the first kind.

**Cost and conditioning** (the wing refinement above, 16 threads):

| element | assembly, L1 / L2 / L3 | quadrature points, L1 / L2 / L3 | factorization and solve | cond₁ | cond₁, rows scaled by the diagonal |
|---|---|---|---|---|---|
| triangles | 0.41 / 2.1 / 13 s | 48M / 270M / 2.2G | 0.01 / 0.27 / 5.2 s | 1.1e4 / 4.4e4 / 1.7e5 | 8.4e3 / 2.8e4 / 1.1e5 |
| quads | 7.7 / 30 / 121 s | 16M / 94M / 0.82G | 0.03 / 0.25 / 5.0 s | 8.9e3 / 3.5e4 / 1.4e5 | 6.1e3 / 2.1e4 / 8.2e4 |

- **Assembly costs 19×, 14× and 9× the triangles'** at L1–L3.
  - The quads use a third of the quadrature points, but each point costs about 60× as much, since its inner integrals are numerical.
  - The ratio falls under refinement because the far field's share grows, and a far quad costs little (a 2×2 outer rule and low-order inner Gauss).
  - On the cube-sphere the ratio is 19×, and on the sectional wings about 10×.
- **The quads condition slightly better:** 0.79× the triangles' cond₁, and 0.72–0.76× with the rows scaled. The unknowns are the same, so the solve costs the same.

**The split triangles remain the default.** On the wing the two elements give the same lift to 0.03% and the same extrapolated induced drag. They give the same section too, once the triangles' mirror line is set aside. The triangles cost a tenth to a nineteenth of the quads' assembly. The quads are the better discretization in three ways:
- **No diagonals.** Lift and drag don't jump when a perturbation flips a split. The quads are 6× less sensitive to skew than the default split, as insensitive as a frozen triangulation, and the tip condition stops mattering.
- **No mirror line** at the plane of symmetry.
- **Smooth curved bodies:** a smaller μ error (0.64× in L2 and 0.27× at its largest on the cube-sphere), and second-order Cp at the quads' centres.

Their price is the assembly cost and a larger quadrature error on badly distorted quads.

A split that never flips under small perturbations would give the triangles the same skew insensitivity at their own cost. Keeping a mesh family's diagonals (`split=:pattern`) does that in these studies, but no general rule of that kind is built yet. Any mirror-symmetric rule would also keep the mirror line.

## Free wake with vortex particles (phase 2)

This is the design doc's first deferred item. Its time stepping mirrors MESA's:
- an implicit trailing-edge panel carries the Kutta condition;
- the free panel rows behind it keep their strength once shed and move with the flow by forward Euler;
- beyond `max_rows` rows the oldest row is converted into vortex particles;
- there is no vortex stretching yet.

The conversion itself follows the design doc, not MESA. It is local to each wake triangle: a triangle's particle comes from its own corner doublets and corners only, and nothing is converted across the interface edge (details below).

`validation/freewake/run_freewake.jl` validates it on the reference wing at α = 5°, started impulsively.

```julia
sim = FreeWakeSim(g, dm, V∞; opts=FreeWakeOptions(dt=0.06, max_rows=10))  # assembles and factors once
run!(sim, 200)                            # sim.history: the loads of every step
write_free_wake_vtk(sim, "out/wing")      # the free rows and the particles, for ParaView
```

**Coupling.**
- **The trailing-edge panel is implicit.** `te_panel` builds one row of wake strips of length te_len (default 0.3 |V∞| dt), closed by a filament along its downstream edge (`WakeOptions(rows, ending=:closed)`). It is assembled with the body, and the bordered matrix is factored once (`factorize_galerkin`). Its doublets follow the Kutta rows, as for the rigid wake.
- **Everything behind it is explicit.** The free rows' sheets and filaments and the particles enter only the right-hand side, through their velocity (`explicit_rhs`). A step costs one right-hand side and one back-substitution with the stored factors.
- **The coupling is exact.** Take a rigid wake, keep its first row as the trailing-edge panel, and freeze the other rows as explicit elements at the implicit solution. The solve gives that solution back to 1e-14 at L1 and L2, for panels 0.25 to 4 trailing-edge panels long.
- **The right-hand side must be integrated with the assembly's own tiers** (contact and near rules). Through a fixed 7-point rule per test triangle, the lift came out 1.6% low with a trailing-edge panel a quarter of a panel long: the first free row sits that close to the trailing edge.

**The free rows.** The rows are stored from the newest, on the trailing-edge panel's downstream edge, to the oldest, each with the μ_w it was shed with. Panels interpolate μ linearly between rows, so shed vorticity is a sheet, not a filament.
- **Filaments** sit only where the layer is cut:
  - the interface along the newest row, with that row's μ, which cancels the trailing-edge panel's downstream filament in steady flow;
  - the tips;
  - the far end, with strength μ − μ_line, which is the starting vortex until the first conversion.
- **One step**, in MESA's order:
  1. solve, then compute the loads;
  2. compute the velocities at every free node and particle;
  3. advance them by forward Euler;
  4. shed a new row on the trailing-edge panel's edge, with this step's μ_w;
  5. convert the oldest row while more than `max_rows` rows remain.
- **The first step** sheds and advects a row before shedding the next one, as MESA does. Otherwise the starting vortex would be released a step late.

**Conversion** of the oldest row:
1. A particle at each triangle's centroid takes α = ∮μ dl = γA, two per quad.
2. Two particles at the 2-point Gauss positions of each tip side and of the far edge take that edge's strength, μ − μ_line. That is the starting vortex at the first conversion and zero afterwards.
3. μ_line then takes the row's upstream doublets, so nothing is created at the interface. With `max_rows = 0` this leaves MESA's leftover line by itself.

A conversion leaves the total vorticity and its first moment (the impulse) unchanged, to rounding.

Each triangle's particle uses only that triangle's corners and corner doublets, as the design doc specifies. MESA's conversion is not local:
- it makes one particle per strip, at the quad's centroid, from the panel's edge circulations;
- it splits each side edge's circulation half and half with the neighbouring strip, so it uses the neighbour's strength;
- it folds the starting vortex into the first particles;
- it leaves a vortex line at the interface.

MESA's constant-strength panels carry all their vorticity on their edges, so its conversion has to share those edges between neighbours. The linear doublet here carries the vorticity as sheets inside the triangles, and there are no filaments between triangles to share.

**Velocities at wake points.** The body uses its exact kernels. Every wake element is a Rosenhead–Moore blob:
- the trailing-edge panel's and free rows' sheets are represented by the particles a conversion would create;
- their filaments by 2-point Gauss blobs;
- the particles by themselves.

Cores are 1.5 × the element's size, and at least |V∞| dt. The exact panel kernels are singular at the wake's own nodes; MESA uses lines with Rankine cores there. With this representation a conversion changes no wake velocity, only the body's right-hand side.

**Loads** (each step):
- **Pressure:** Cp = 1 − |u|²/V∞² − 2μ̇/V∞², with μ̇ by backward difference.
- **Kutta–Joukowski** from μ_w; during a transient this is only the circulatory part.
- **Impulse:** F = −(2/V∞²) dI/dt, with I = ½Σ x × α over all vorticity: body sheets, trailing-edge panel, free rows and particles.
- **Check:** the total vorticity Ω.

**Impulsive start with a prescribed wake** (rows moved with V∞, no conversion, L1, dt = trailing-edge panel size):

| wake length (spans) | CL (pressure) | CL (KJ) | CL (impulse) | steady CL (KJ), closed wake of the same rows | gap | closed wake − open wake |
|---|---|---|---|---|---|---|
| 0.26 | 0.2715 | 0.2363 | 0.2698 | 0.3006 | −6.4e-2 | −4.1e-2 |
| 1.0 | 0.3253 | 0.3238 | 0.3284 | 0.3364 | −1.3e-2 | −5.5e-3 |
| 2.0 | 0.3353 | 0.3382 | 0.3389 | 0.3403 | −2.1e-3 | −1.6e-3 |
| 3.0 | 0.3372 | 0.3406 | 0.3407 | 0.3411 | −5.5e-4 | −7.4e-4 |
| 5.0 | 0.3379 | 0.3415 | 0.3415 | 0.3416 | −9.4e-5 | −2.8e-4 |

- **Convergence to the steady state:**
  - The marched lift converges onto the steady solve on the same rows, faster than that solve converges to the open rigid wake (whose gap falls as (b/L)²).
  - At 5 spans the lift is 0.1% below the rigid wake's (CL 0.3383, CL_KJ 0.3418), and the impulse route matches Kutta–Joukowski to 1e-4.
- **Kelvin's theorem** holds to rounding: |Ω_total| ≤ 1.6e-15.
- **Time step:** at a fixed ratio te_len/(|V∞| dt), halving or doubling dt changes CL at 2 spans by under 0.06%.
- **The trailing-edge panel's length matters only in the transient.**
  - At 1 span, te_len = 0.1 |V∞| dt gives 0.313 and |V∞| dt gives 0.331.
  - At 2 spans the gap has narrowed to 0.336 against 0.339.
  - In steady flow the length doesn't matter at all (the exact-coupling check above).
  - 0.3 |V∞| dt is about where Katz & Plotkin place the newest shed vortex.

**Conversion** (prescribed wake, L1, 5 spans). Lift change from never converting:

| `max_rows` | 0 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| ΔCL (KJ), at 5 spans | +1.55% | +0.74% | +0.44% | +0.21% | +0.08% | +0.02% |
| largest \|ΔCL (KJ)\| over the run | 4.4e-2 | 1.7e-2 | 9.1e-3 | 3.6e-3 | 1.1e-3 | 2.5e-4 |

Particles close to the body stand for the sheet less well than its panels, so converting further downstream costs less. MESA's example uses 10 rows; at 10 the change is under 0.1%.

**Free wake** (forward Euler, `max_rows` = 10, to 5 spans):

| level | steps | time | particles | CL (pressure) at 5 spans | CL (KJ) at 5 spans | extrapolated in 1/L²: CL, CL (KJ) | rigid wake: CL, CL (KJ) |
|---|---|---|---|---|---|---|---|
| L1 | 232 | 18 s | 8,470 | 0.3367 | 0.3404 | 0.3372, 0.3410 | 0.3383, 0.3418 |
| L2 | 463 | 7 min | 32,684 | 0.3372 | 0.3400 | 0.3377, 0.3406 | 0.3389, 0.3415 |

- **Roll-up lowers the lift slightly.** Free against prescribed at 5 spans (L1) is −0.41% by pressure and −0.37% by Kutta–Joukowski. Extrapolated, the free wake sits 0.24–0.35% below the rigid wake at both levels.
- **The settings hardly matter** (L1, at 5 spans):
  - particle cores 1.0× and 2.5× the element size: ±0.05%;
  - `max_rows` 4 and 20: +0.13% and +0.01%;
  - dt doubled or halved: under 0.01%.
- **It stays well behaved.** Velocities stay below 1.03 |V∞|.
- **Without stretching the total vorticity drifts.** The far edge of the panel rows stretches while the particles' strengths stay fixed. By 5 spans |Ω_total| is 0.7% (L1) and 4% (L2) of the largest μ_w times the span. The impulse route depends on that drift, so in free mode it is unreliable: +1.7% against Kutta–Joukowski at L2. The pressure and Kutta–Joukowski routes are unaffected.
- **Cost.** At the end of the L2 run a step takes 1.8 s on 16 threads, with 32,684 particles: direct particle sums, the body's velocity at every particle, and every particle in the right-hand side. L3, or longer wakes, would need particle deletion (`delete_beyond`) or a treecode.

**Against MESA's free wake** (L1, the same mesh, dt, te_len = 0.045, `max_rows` 10 and a single particle core 1.5 × the trailing-edge panel size in both):

| | ours CL (pressure) | ours CL (KJ) | MESA CL |
|---|---|---|---|
| prescribed wake, 5 spans | 0.3380 | 0.3416 | 0.2991 |
| free wake, 5 spans | 0.3370 | 0.3406 | 0.2984 |
| rigid wake (stages 5 and 6) | 0.3383 | 0.3418 | 0.3287 |

| te_len (prescribed wake, 2 spans) | 0.0146 | 0.0292 | 0.0585 | 0.117 |
|---|---|---|---|---|
| ours CL | 0.3352 | 0.3358 | 0.3362 | 0.3365 |
| MESA CL | 0.2814 | 0.2909 | 0.3021 | 0.3123 |

- **MESA's free-wake lift is 9% below its own rigid wake, and it depends on the trailing-edge panel's length.** It rises 11% from te_len 0.0146 to 0.117.
  - This is the error its dev notes describe (dev/wake_coupling_findings.md): the trailing-edge panel enters its matrix through the potential, while the shed wake enters through velocity.
  - Ours varies by 0.4% at 2 spans, the transient described above, and agrees with our rigid wake.
- **Roll-up costs both codes about as much lift:** −0.32% for ours and −0.24% for MESA, free against prescribed.

**Abused meshes** (MESA's smoke set and all abuses at once, L1, free wake to 5 spans):
- Every mesh runs, with velocities below 1.03 |V∞|.
- The lift changes match the rigid wake's of stage 5 to within 0.04%. Opposite lean of 0.3 gives −0.93%, misaligned trailing edges −0.60%, independent jitter −0.49%, jitter 0.30%, and T-junctions, gaps and slivers under 0.12%.

## Test weights on the wing (validation stage 8)

`validation/test_weights/run_test_weights.jl` solves the reference wing at α = 5°, with the default options otherwise, using the three test weights of `AssemblyOptions(test=...)`. Only the weights change; the quadrature points are the same:
- **hat** (default): the nodal hat functions φᵢ, which makes the method a Galerkin method;
- **dual cell:** one on the node's median dual cell, the part of each triangle where the node's barycentric coordinate is the largest;
- **one-ring:** one on every triangle around the node.

**Refinement** on clean meshes:

| level | weight | CL | CL (KJ) | CDi (pressure) | CDi (Trefftz) | Cm c/4 | cond₁ |
|---|---|---|---|---|---|---|---|
| L1 | hat | 0.3383 | 0.3418 | 0.01072 | 0.00895 | +0.0012 | 1.1e4 |
| | dual cell | 0.3372 | 0.3412 | 0.01138 | 0.00892 | +0.0012 | 1.4e4 |
| | one-ring | 0.3499 | 0.3429 | 0.00570 | 0.00903 | +0.0010 | 3.3e4 |
| L2 | hat | 0.3389 | 0.3415 | 0.01011 | 0.00927 | +0.0007 | 4.4e4 |
| | dual cell | 0.3385 | 0.3413 | 0.01032 | 0.00925 | +0.0007 | 5.5e4 |
| | one-ring | 0.0141 | 0.3184 | −0.520 | 0.01724 | +0.457 | 2.9e6 |
| L3 | hat | 0.3395 | 0.3411 | 0.00981 | 0.00943 | +0.0002 | 1.7e5 |
| | dual cell | 0.3393 | 0.3410 | 0.00993 | 0.00942 | +0.0001 | 2.1e5 |
| | one-ring | −9.63 | 0.3521 | −5.31 | 0.01745 | +1.51 | 3.1e7 |

**Sectional accuracy:** the midspan section of wings of aspect ratio 15.5, 31 and 62 (68 strips), extrapolated to the infinite wing and compared with the exact section, as in stage 6:

| panels per surface | hat cl | hat cl (KJ) | hat cm | dual cl | dual cl (KJ) | dual cm |
|---|---|---|---|---|---|---|
| 19 | +1.28% | +1.21% | −0.0022 | +0.81% | +1.22% | −0.0024 |
| 38 | +0.66% | +0.61% | −0.0010 | +0.46% | +0.63% | −0.0011 |
| 76 | +0.35% | +0.32% | −0.0004 | +0.25% | +0.33% | −0.0005 |

The one-ring weight is already unstable at 19 panels (cond₁ 5e7), so its finer sections are skipped.

These sections use the default split, which mirrors its diagonals at midspan, where they are measured, and dents the circulation there (stage 7). Both weights are measured on the same meshes at the same place, so the comparison is like for like. The hat's own errors, measured with one diagonal throughout, are +1.57%, +0.78% and +0.40% in cl (stage 6).

**The hat weight stays the default.** It is the most self-consistent of the three: in the section its pressure and circulation lift agree within 0.07%, and it conditions best.

**The dual cell is a sound alternative.**
- **It solves for the same circulation.** Its Kutta–Joukowski lift is within 0.2% of the hat's on every mesh, and its sectional circulation error matches the hat's at every resolution (first order).
- **Its pressure route differs.** Its pressure lift sits below its own circulation lift by 0.4%, 0.17% and 0.08% at 19, 38 and 76 panels. That lands closer to the exact section (+0.25% against +0.35% at 76 panels), but the solution is no better: the circulation errors are equal. Its pressure drag is 6%, 2% and 1% above the hat's at L1–L3, and extrapolates to 0.00971 against 0.00953.
- **Pressure:** at midspan (L2) its Cp is within 0.006 of the hat's between 10% and 90% chord; the suction peaks are −1.214 and −1.223.
- **Conditioning:** about 1.25× the hat's.
- **Abused meshes** (MESA's smoke set and all at once): mean (largest) |ΔCL| 0.31% (1.07%) at L1 and 0.13% (0.43%) at L2, against the hat's 0.26% (0.91%) and 0.10% (0.36%). That is 20–30% more sensitive.
- **Zero incidence (L2):** 6e-6 on the clean mesh and at most 7e-5 on the skewed and non-matching ones (hat: 8e-5). With all abuses at once, both reach 2e-4.

**The one-ring weight is unstable and should not be used.**
- **The system degenerates.** Its rows integrate over overlapping one-rings, so on thin elements neighbouring rows nearly coincide. cond₁ is 3.3e4 on the clean L1 wing, 2.9e6 at L2 and 3.1e7 at L3, 1.7e6 on the L1 sliver-station mesh and 5e7 on the stretched sections.
- **Its lift is then meaningless:** CL 0.014 at L2 and −9.6 at L3, and a 2,199% change on the L1 sliver stations.
- **The quadrature is not to blame:** its column defects stay at 7e-5 to 5e-3, within 2× of the hat's on the same meshes.
- **Even where it still solves** (clean L1), its pressure is off: the suction peak is 37% above the hat's and the pressure drag half of it, while its circulation lift is within 0.3%.
- On the sphere (stage 3) it stayed stable but was 4–5× less accurate than the hat weight.
