"""
    GeometryOptions(; kwargs...)

Tolerances and switches of the geometry pipeline. Tolerances are relative to the local edge length.

Welding and quad splitting
- `weld_rel = 1e-3`: nodes closer than `weld_rel` times their shortest incident edge are merged.
- `split = :auto`: quad diagonal. `:auto` takes the shorter diagonal, or the one with the smaller fold when the quad
  folds by more than `split_fold_deg` along either diagonal. Also `:shortest`, `:min_fold`, `:d13`, `:d24` (fixed
  diagonals, `:d13` as `MeshLibrary.split_quads`) and `:pattern`, which takes the diagonal of input element e from
  `split_pattern[e]` (1: nodes 1-3, 2: nodes 2-4), e.g. the pattern of the unperturbed mesh in a skew study.
- `split_tie_rel = 1e-6`: diagonals within this relative length are tied. Ties are broken by a rule that is invariant
  to the element's node order and mirror symmetric about the plane `x[split_mirror[1]] = split_mirror[2]` (set
  `nothing` to skip it), then by the diagonal through the lowest node id.

Pathologies
- `tj_rel = 1e-3`: a boundary node within `tj_rel` times an edge's length of its interior hangs on it (T-junction);
  `tj_tmin` excludes the ends of the edge.
- `gap_rel = 0.05`, `gap_angle_deg = 10`, `gap_overlap = 0.25`: two boundary edges form part of a gap when they are
  within `gap_angle_deg` of parallel, overlap by `gap_overlap` of the shorter one, and are everywhere within `gap_rel`
  times the shorter one's length of each other (Hausdorff distance over the overlap).
- `openings = :cap`: an opening (a loop enclosing an area, such as an open wing tip) is capped with a ladder of
  triangles on its own nodes (`:cap`), or its two sides are tied in the constraint matrix without panels (`:zip`), or
  it is left open (`:open`, for debugging).

Trailing edge
- `te_angle_deg = 30`: a convex edge whose faces meet at an included angle of at most `te_angle_deg`, pointing along
  `wake_dir`, is a Kutta edge (the meaning of PAGE's `TrailingEdgeParams.te_angle`).
- `kutta_source = :dihedral`: Kutta edges from the dihedral test (`:dihedral`), from the input edge tags named in
  `kutta_tags` (`:tags`), or the union of both (`:both`). The tags are always cross-checked against the dihedral test.
- `kutta_include`, `kutta_exclude`: input node pairs forced to be, or not to be, Kutta edges.
- `wake_dir = +x`, `up_ref = +z`: the trailing edge direction τ is chosen with (wake_dir × τ)·up_ref ≥ 0, so that
  "upper" is the side facing `up_ref`.

Quality
- `sliver_angle_deg = 5`, `sliver_aspect = 20`: triangles with a smaller angle or a larger aspect ratio are slivers.
"""
Base.@kwdef struct GeometryOptions
    weld_rel::Float64 = 1e-3
    split::Symbol = :auto
    split_pattern::Vector{Int8} = Int8[]
    split_fold_deg::Float64 = 1.0
    split_tie_rel::Float64 = 1e-6
    split_mirror::Union{Nothing,Tuple{Int,Float64}} = (2, 0.0)
    tj_rel::Float64 = 1e-3
    tj_tmin::Float64 = 1e-3
    gap_rel::Float64 = 0.05
    gap_angle_deg::Float64 = 10.0
    gap_overlap::Float64 = 0.25
    openings::Symbol = :cap
    te_angle_deg::Float64 = 30.0
    kutta_source::Symbol = :dihedral
    kutta_tags::Vector{String} = ["trailing_edge", "trailing_edge_1", "trailing_edge_2"]
    kutta_include::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    kutta_exclude::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    wake_dir::SVector{3,Float64} = SA[1.0, 0.0, 0.0]
    up_ref::SVector{3,Float64} = SA[0.0, 0.0, 1.0]
    sliver_angle_deg::Float64 = 5.0
    sliver_aspect::Float64 = 20.0
end

function check_options(opts::GeometryOptions)
    opts.split in (:auto, :shortest, :min_fold, :d13, :d24, :pattern) || throw(ArgumentError("Unknown split rule $(opts.split)"))
    opts.openings in (:cap, :zip, :open) || throw(ArgumentError("Unknown openings mode $(opts.openings)"))
    opts.kutta_source in (:dihedral, :tags, :both) || throw(ArgumentError("Unknown kutta_source $(opts.kutta_source)"))
    return opts
end
