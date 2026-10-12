"""Kind of a triangle: part of the input surface, or added to close an opening"""
@enum TriKind::UInt8 BodyTri = 1 CapTri = 2

"""Topological class of an edge"""
@enum EdgeClass::UInt8 begin
    Interior = 1  # shared by two triangles, μ continuous
    Kutta = 2     # trailing edge shared by two body triangles: μ duplicated, wake attached
    KuttaPair = 3 # boundary edge on one side of a non-matching trailing edge
    Gap = 4       # boundary edge on one side of a stitched small gap
    TJLong = 5    # long edge of a T-junction
    TJShort = 6   # short edge along a T-junction's long edge
    Open = 7      # boundary of an opening left open or zipped: carries a filament
    CapSeam = 8   # shared by a cap triangle and a body triangle
end

"""Side of the trailing edge a triangle corner at a trailing edge node belongs to"""
@enum TESide::UInt8 NoTE = 0 TEUpper = 1 TELower = 2 TECap = 3

"""A repair or diagnostic event of the geometry pipeline, with where it happened and its size"""
struct RepairEvent{T}
    kind::Symbol
    location::SVector{3,T}
    magnitude::T
    detail::String
end

"""An edge shared by more than two faces: intersecting or non-manifold surfaces are out of scope"""
struct NonManifoldEdgeError{T} <: Exception
    nodes::Tuple{Int,Int}
    coords::Tuple{SVector{3,T},SVector{3,T}}
    nfaces::Int
end

function Base.showerror(io::IO, e::NonManifoldEdgeError)
    print(io, "NonManifoldEdgeError: edge $(e.nodes) from $(e.coords[1]) to $(e.coords[2]) is shared by $(e.nfaces) faces. ",
        "Intersecting or non-manifold surfaces are not supported.")
end

"""Hanging node `hanging` on the edge from `a` to `b` of another element, at x ≈ (1 - t) x_a + t x_b"""
struct TJunction{T}
    hanging::Int
    a::Int
    b::Int
    t::T
end

"""
Two nearly coincident boundary chains: a stitched small gap (`kind = :gap`), the two sides of a non-matching trailing edge
(`:kutta`), or the two halves of a zipped opening (`:zip`). `proj_a[k] = (s, t)` places node `chain_a[k]` on segment
`s` of `chain_b` at parameter `t` (and the same for `proj_b`). `width` is the largest distance between the chains.
"""
struct ChainPair{T}
    kind::Symbol
    chain_a::Vector{Int}
    chain_b::Vector{Int}
    proj_a::Vector{Tuple{Int,T}}
    proj_b::Vector{Tuple{Int,T}}
    width::T
end

"""
A closed loop of boundary edges enclosing an area, such as an uncapped wing tip. `chain_upper` and `chain_lower` run
from the leading edge corner `le` to the trailing edge corner `te`; `mode` is `:cap`, `:zip` or `:open`, and
`cap_tris` the cap triangles. `planarity` is the largest distance of the loop from its best fit plane over the loop
diameter.
"""
struct Opening{T}
    loop::Vector{Int}
    le::Int
    te::Int
    chain_upper::Vector{Int}
    chain_lower::Vector{Int}
    mode::Symbol
    cap_tris::UnitRange{Int}
    planarity::T
end

"""
Kutta edges ordered into a chain along the direction τ, with arclength `s` at the nodes. Segment k joins `nodes[k]` and
`nodes[k+1]`, is edge `edges[k]` and separates `tri_upper[k]` (which owns the half-edge along τ) from `tri_lower[k]`.
`end_type` is `:cap`, `:open`, `:pair` (continues as a non-matching trailing edge), `:branch` or `:dangling` for the
first and last node, and `end_cap_tri` a cap triangle at each end node (0 if none).
"""
struct KuttaChain{T}
    nodes::Vector{Int}
    s::Vector{T}
    edges::Vector{Int}
    tri_upper::Vector{Int}
    tri_lower::Vector{Int}
    closed::Bool
    end_type::NTuple{2,Symbol}
    end_cap_tri::NTuple{2,Int}
end

"""
A non-matching trailing edge: the upper and lower boundary chains along τ with the triangle of each segment, and the
union of their nodes along τ (the wake's trailing edge nodes). `union_side` is +1 for nodes of the upper chain only,
-1 for the lower chain only and 0 for nodes on both; `union_upper[k] = (s, t)` locates union node k on segment s of the
upper chain (and `union_lower` on the lower chain), and `s` is its arclength along the upper chain.
"""
struct KuttaPairChain{T}
    pair::Int
    upper::Vector{Int}
    lower::Vector{Int}
    tri_upper::Vector{Int}
    tri_lower::Vector{Int}
    union_nodes::Vector{Int}
    union_side::Vector{Int8}
    union_upper::Vector{Tuple{Int,T}}
    union_lower::Vector{Tuple{Int,T}}
    s::Vector{T}
end

"""
Trailing edge data: Kutta edge chains, non-matching trailing edge pairs, the trailing edge nodes, and for every
triangle corner at a trailing edge node the side it belongs to (`corner_side[t][k]` for node `tris[t][k]`).
"""
struct KuttaData{T}
    chains::Vector{KuttaChain{T}}
    pairs::Vector{ChainPair{T}}
    pair_chains::Vector{KuttaPairChain{T}}
    te_nodes::Vector{Int}
    corner_side::Vector{SVector{3,TESide}}
    tag_mismatch::Int
end

"""Unique edges: end nodes (a < b), their half-edges (second is 0 on a boundary edge) and their class"""
struct EdgeTable
    nodes::Vector{SVector{2,Int}}
    he::Vector{SVector{2,Int}}
    class::Vector{EdgeClass}
end

"""
A connected body. `closed` is false when an opening was left open or zipped, in which case a uniform μ is not a null
vector of the system and the body needs no bordered constraint row. `euler` is V - E + F, with T-junction, gap and
non-matching trailing edge edges counted as they are: 2 for a conforming closed body.
"""
struct Body{T}
    tris::Vector{Int}
    nodes::Vector{Int}
    closed::Bool
    euler::Int
    volume::T
end

"""Per triangle quality (angles in degrees) and summary counts and magnitudes"""
struct QualityReport{T}
    min_angle::Vector{T}
    max_angle::Vector{T}
    aspect::Vector{T}
    skew::Vector{T}
    area::Vector{T}
    sliver::BitVector
    counts::Dict{Symbol,Int}
    stats::Dict{Symbol,T}
end

"""
A closed, consistently oriented triangulation with classified edges, as consumed by the DOF map, wake and assembly.

- `x`: welded node coordinates. `tris`: triangles with outward normal ∝ (x₂-x₁)×(x₃-x₁).
- `tri_kind`, `tri_parent` (input element, 0 for caps), `tri_slot` (0 for an input triangle, 1 or 2 for the halves of
  a quad), `tri_body`.
- `quads`, `quad_parent`, `quad_diag`: the welded, outward parent quads and the diagonal they were split along
  (1: nodes 1-3, 2: nodes 2-4), kept for the isoparametric quad element.
- Half-edge h = 3(t-1)+k runs from `tris[t][k]` to `tris[t][mod1(k+1,3)]`; `twin[h]` is its opposite half-edge (0 on a
  boundary) and `he_edge[h]` its edge.
- `input_node`: welded node of every input node.
- `tjunctions`, `gaps`, `openings`, `kutta`: the topology records the constraint matrix C is built from.
- `touch_pairs`: triangles that touch without sharing an edge: along a T-junction's long edge, or facing each other
  across a gap or a non-matching trailing edge (they may share an end node), for the contact quadrature tier.
"""
struct SurfaceGeometry{T}
    x::Vector{SVector{3,T}}
    tris::Vector{SVector{3,Int}}
    tri_kind::Vector{TriKind}
    tri_parent::Vector{Int}
    tri_slot::Vector{Int8}
    tri_body::Vector{Int}
    quads::Vector{SVector{4,Int}}
    quad_parent::Vector{Int}
    quad_diag::Vector{Int8}
    twin::Vector{Int}
    he_edge::Vector{Int}
    edges::EdgeTable
    input_node::Vector{Int}
    bodies::Vector{Body{T}}
    tjunctions::Vector{TJunction{T}}
    gaps::Vector{ChainPair{T}}
    openings::Vector{Opening{T}}
    kutta::KuttaData{T}
    touch_pairs::Vector{Tuple{Int,Int}}
    quality::QualityReport{T}
    log::Vector{RepairEvent{T}}
    edge_tags::Dict{String,Vector{Tuple{Int,Int}}}
    face_tags::Dict{String,Vector{Int}}
    options::Any
end

nnodes(g::SurfaceGeometry) = length(g.x)
ntris(g::SurfaceGeometry) = length(g.tris)
nedges(g::SurfaceGeometry) = length(g.edges.nodes)

function Base.show(io::IO, g::SurfaceGeometry{T}) where {T}
    q = g.quality.counts
    print(io, "SurfaceGeometry{$T}: $(nnodes(g)) nodes, $(ntris(g)) triangles ($(count(==(CapTri), g.tri_kind)) cap), ",
        "$(nedges(g)) edges, $(length(g.bodies)) bod$(length(g.bodies) == 1 ? "y" : "ies"), ",
        "$(get(q, :kutta_edges, 0)) Kutta edges, $(length(g.tjunctions)) T-junctions, $(length(g.gaps)) gaps, ",
        "$(length(g.openings)) openings")
end
