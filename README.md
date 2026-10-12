# VortexPanel.jl

Galerkin linear-vortex panel solver. Unknowns are nodal doublet strengths on flat triangles. Each triangle carries a constant vortex sheet γ = n×∇μ. Sources are prescribed. The Kutta condition is applied at the trailing-edge nodes, and the wake is rigid.

The design is in `Plan Galerkin Linear-Vortex Panel Solver (Rigid Wake Phase).md`. The reference solver is [MESA.jl](../MESA.jl), a Morino code that is used unmodified from the validation scripts.

## Status

Steps 1 and 2 are done. They cover:

- **`VortexPanel.Geometry`** (step 1): turns surface meshes of quads and triangles into a closed, outward, classified triangulation. It welds nodes, splits quads, detects T-junctions, gaps and non-matching trailing edges, and caps open tips without adding nodes. It also detects Kutta edges and reports mesh quality.
- **`VortexPanel.MeshTools`** (step 1): generates closed test bodies and seeded, reproducible abused meshes. The abuse covers jitter, spanwise skew, T-junctions, gaps, slivers and non-matching trailing edges.
- **MESA baseline** (step 1): a rigid-wake baseline on the reference wing family.
- **`VortexPanel.DOFMap`** (step 2): the constraint matrix C, with μ_slots = C q, and the nodal Kutta rows K, with μ_w = K q.

Not yet built: kernels, elements, wake geometry, assembly and post-processing.

## Layout

```
src/geometry/    Geometry: build_geometry(mesh; opts) -> SurfaceGeometry, print_report, write_vtk, add_caps
src/dofmap/      DOFMap: build_dofmap(g; opts) -> DOFMapping (C, K, slots), reduce_system, bordered_rows, wake_rows
src/meshtools/   MeshTools: icosphere, cube_sphere, uv_sphere, abuse, skew_wing, ...
test/            Pkg.test() (needs the local PAGE release tree)
validation/      reference wing, MESA baseline, geometry report (own environment)
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
```

Outputs go to `validation/output/`, which is not tracked:

- `morino/summary.csv`, which records the git SHAs of the checkouts;
- spanwise loading, Cp and trailing-edge μ CSVs, and VTK;
- `geometry/*.vtu` with a geometry report CSV.

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
