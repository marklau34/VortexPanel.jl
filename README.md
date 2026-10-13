# VortexPanel.jl

Galerkin linear-vortex panel solver. Unknowns are nodal doublet strengths on flat triangles. Each triangle carries a constant vortex sheet γ = n×∇μ. Sources are prescribed. The Kutta condition is applied at the trailing-edge nodes, and the wake is rigid.

The design is in `Plan Galerkin Linear-Vortex Panel Solver (Rigid Wake Phase).md`. The reference solver is [MESA.jl](../MESA.jl), a Morino code that is used unmodified from the validation scripts.

## Status

Steps 1–4 are done, through validation stage 4 (ellipsoids in uniform flow). They cover:

- **`VortexPanel.Geometry`** (step 1): turns surface meshes of quads and triangles into a closed, outward, classified triangulation. It welds nodes, splits quads, detects T-junctions, gaps and non-matching trailing edges, and caps open tips without adding nodes. It also detects Kutta edges and reports mesh quality.
- **`VortexPanel.MeshTools`** (step 1): generates closed test bodies and seeded, reproducible abused meshes. The abuse covers jitter, spanwise skew, T-junctions, gaps, slivers and non-matching trailing edges. Jitter moves each node by up to δ times its clearance, its smallest height over an edge of an incident element, so elements skew without collapsing.
- **MESA baseline** (step 1): a rigid-wake baseline on the reference wing family.
- **`VortexPanel.DOFMap`** (step 2): the constraint matrix C, with μ_slots = C q, and the nodal Kutta rows K, with μ_w = K q.
- **`VortexPanel.Kernels`** (step 3): the velocity of a unit-strength source triangle, w = (1/4π)∫(x−y)/|x−y|³ dA, the one kernel of the method (sources σw, vortex sheets γ×w).
- **`VortexPanel.Elements`, `Assembly`, `Post`** (step 4): the flat-triangle element and quadrature tiers, the Galerkin assembly with the bordered dense solve, and surface velocity, Cp and leak. Validated on the sphere (validation stage 3) and on ellipsoids (stage 4) in uniform flow, on clean and abused meshes.

Not yet built: wake geometry and the lifting solve, forces (pressure, Kutta–Joukowski, Trefftz), the isoparametric quad element.

## Layout

```
src/geometry/    Geometry: build_geometry(mesh; opts) -> SurfaceGeometry, print_report, write_vtk, add_caps
src/dofmap/      DOFMap: build_dofmap(g; opts) -> DOFMapping (C, K, slots), reduce_system, bordered_rows, wake_rows
src/kernels/     Kernels: SourceTriangle, source_velocity, sheet_velocity
src/elements/    Elements: FlatTriangle (hat basis, γ_b), quadrature rules and tiers (TierRules)
src/assembly/    Assembly: assemble(g, dm, els, σ, V∞; opts) -> A, b (tiers incl. contact_rule); solve_galerkin; consistency
src/post/        Post: surface_velocity (singularity sums), gradient_velocity, pressure_coefficient, nodal_average
src/meshtools/   MeshTools: icosphere, cube_sphere, uv_sphere, abuse, skew_wing, ...
test/            Pkg.test() (needs the local PAGE release tree)
validation/      reference wing, MESA baseline, geometry report, sphere and ellipsoid studies (own environment);
                 common/exact_flows.jl: the exact ellipsoid flow; common/closed_body.jl: the closed-body study tables
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
```

Outputs go to `validation/output/`, which is not tracked:

- `morino/summary.csv`, which records the git SHAs of the checkouts;
- spanwise loading, Cp and trailing-edge μ CSVs, and VTK;
- `geometry/*.vtu` with a geometry report CSV;
- `sphere/` convergence (clean and jittered), abused-mesh, quadrature-sensitivity and test-weight CSVs, and VTK of the finest icosphere and the fully abused icosphere 3 (Cp, μ, leak, their exact values and the μ error);
- `ellipsoid/` the exact-solution checks, convergence, abused-mesh and quadrature-sensitivity CSVs, Cp along the cut y = 0 of the triaxial body, and VTK of the finest triaxial, prolate and thin oblate meshes and the fully abused triaxial one.

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
- **Cap corners at the tip trailing-edge node** use the upper copy (`DOFMapOptions(cap_te_copy=:upper)`). The jump μ_u − μ_l then runs as a filament along the cap's seam with the lower surface.
- **Constraints.** Each constrained slot is an affine combination of other slots, resolved to the unknowns q:
  - **T-junction hanging nodes** are linear along the long edge.
  - **Stitched gaps:** the nodes of the side with more nodes are linear along the coarser side, at their projection.
  - **Zipped openings:** the nodes of one side are linear along the other side, at equal normalised arclength. Trailing-edge nodes keep both copies, so zipping a tip does not tie the Kutta jump.
- **Row checks.** Rows of C sum to one, so a linear μ is reproduced exactly through T-junctions and coincident stitches. A gap of width w adds at most |∇μ|·w.
- **Kutta rows.** Each wake trailing-edge node k has a row giving μ_w,k = μ_u − μ_l.
  - On a non-matching trailing edge, these nodes are the union of both chains' nodes, and μ_u and μ_l are interpolated along their own chains.
  - Rows of K sum to zero, so a uniform μ sheds nothing.
  - `wake_rows(dm, te_of_wake_node)` gives the C rows of the wake nodes once the strips exist.
- **Reduced system.** It is (Cᵀ A [C; W]) q = Cᵀ b (`reduce_system`). Each closed body adds a bordered row: `bordered_rows` gives e and c = Cᵀ∫φ dA.

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

Solving a problem takes six calls:

```julia
g   = build_geometry(mesh)                   # closed, classified triangulation
dm  = build_dofmap(g)                        # slots, C, Kutta rows
els = build_elements(g)                      # flat triangles with hat functions
σ   = source_strengths(els, V∞)              # σ = −V∞·n
A, b, stats = assemble(g, dm, els, σ, V∞)    # A μ = b states ∫ψᵢ n·u dA = 0
sol = solve_galerkin(dm, A, b)               # q, μ = C q, λ, cond
```

**Formulation.**
- The doublet μ is piecewise linear, and each triangle carries the vortex sheet γ = n × ∇μ.
- A unit nodal value at corner b of triangle T gives γ_b = −e_b/(2A_T), where e_b is the edge opposite b, counterclockwise about the outward normal.
- The boundary condition is ∫ψᵢ n·(V∞ + u_σ + u_γ) dA = 0, with velocities taken on the fluid side.
- The test functions ψ are the hat functions by default (`test=:onering` or `:dual` select the alternatives).
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
