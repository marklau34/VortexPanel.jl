"""
    WakeOptions(; kwargs...)

Geometry of the rigid wake (doc §7).
- `direction = :freestream`: the strips leave the trailing edge along the freestream; `:bisector` takes the first row
  along the trailing edge bisector and the rows after it along the freestream.
- `first_row`: length of the first row; `nothing` takes the size of the trailing edge panels (the median height of the
  upper and lower trailing edge triangles over their trailing edge).
- `growth = 1.2`: ratio of successive row lengths.
- `length = 20.0`: the wake's length along the freestream in spans (the largest distance between two wake trailing
  edge nodes); the last row ends at or beyond it.
- `ending = :truncate`: the strips end there and the starting vortex is left out (no filament crosses the far end);
  `:semi_infinite` adds a last row reaching `far_factor` times further, standing for strips that run to infinity.
"""
Base.@kwdef struct WakeOptions
    direction::Symbol = :freestream
    first_row::Union{Nothing,Float64} = nothing
    growth::Float64 = 1.2
    length::Float64 = 20.0
    ending::Symbol = :truncate
    far_factor::Float64 = 1e3
end

function check_options(o::WakeOptions)
    o.direction in (:freestream, :bisector) || throw(ArgumentError("direction must be :freestream or :bisector, not $(o.direction)"))
    o.ending in (:truncate, :semi_infinite) || throw(ArgumentError("ending must be :truncate or :semi_infinite, not $(o.ending)"))
    o.growth >= 1 && o.length > 0 && o.far_factor > 1 || throw(ArgumentError("growth must be ≥ 1, length > 0 and far_factor > 1"))
    return o
end

"""
    RigidWake

Rigid wake strips (doc §7). Node i sits at `x[i]` in row `row[i]` (0 on the trailing edge, at distance `rows[row]`
from it) and carries the doublet strength of the trailing edge node its strip starts from, `te[i]` (a row of the
Kutta matrix K). Every node of a strip carries the same μ_w, so the wake sheds trailing vorticity only. The triangles
`tris` face the upper side (their normal is wake_dir × τ), so the trailing edge filaments of the upper surface, the
lower surface and the wake cancel under μ_w = μ_u - μ_l. `filaments` are the directed edges along the outer sides of
open trailing edge lines, where the layer ends and a filament of strength μ_w remains (the tip vortices). A doublet
layer of strength μ is bounded by a filament of strength μ running clockwise about its normal (where μ drops to zero
across a thin band, γ = n × ∇μ points against the counterclockwise tangent), so the filaments run upstream along the
first trailing edge node's side and downstream along the last's, continuing the trailing vorticity of the wake sheet
(for a lifting wing, the right tip vortex runs downstream). No filament crosses the far end (no spurious starting
vortex).
"""
struct RigidWake{T}
    x::Vector{SVector{3,T}}
    te::Vector{Int}
    row::Vector{Int}
    tris::Vector{SVector{3,Int}}
    filaments::Vector{Tuple{Int,Int}}
    rows::Vector{T}
    span::T
    options::WakeOptions
end

nwake(w::RigidWake) = length(w.x)

function Base.show(io::IO, w::RigidWake{T}) where {T}
    print(io, "RigidWake{$T}: $(nwake(w)) nodes in $(length(w.rows)) rows to $(round(w.rows[end] / w.span, sigdigits=3)) spans, ",
        "$(length(w.tris)) triangles, $(length(w.filaments)) tip filament segments")
end

"""
Size of the trailing edge panels: the median height of the trailing edge triangles (upper and lower, of matching and
non-matching trailing edges) over their trailing edge
"""
function te_panel_size(g::SurfaceGeometry{T}) where {T}
    h = T[]
    height(t, a, b) = (c = only(v for v in g.tris[t] if v != a && v != b);
    norm(cross(g.x[b] - g.x[a], g.x[c] - g.x[a])) / norm(g.x[b] - g.x[a]))
    for ch in g.kutta.chains, (s, e) in enumerate(ch.edges)
        a, b = g.edges.nodes[e]
        push!(h, height(ch.tri_upper[s], a, b), height(ch.tri_lower[s], a, b))
    end
    for pch in g.kutta.pair_chains
        for s in eachindex(pch.tri_upper)
            push!(h, height(pch.tri_upper[s], pch.upper[s], pch.upper[s+1]))
        end
        for s in eachindex(pch.tri_lower)
            push!(h, height(pch.tri_lower[s], pch.lower[s], pch.lower[s+1]))
        end
    end
    isempty(h) && throw(ArgumentError("no trailing edge panels"))
    sort!(h)
    return h[(length(h)+1)÷2]
end

"""
    build_wake(g, dm, V∞; opts=WakeOptions())

Rigid wake strips from the wake's trailing edge nodes (`dm.wake_te`): rows at distances 0, ℓ₁, ℓ₁(1 + r), ... from
the trailing edge (ℓ₁ the first row, r the growth ratio) until the wake's length, the first row along the freestream or
the trailing edge bisector and the others along the freestream. Every pair of neighbouring trailing edge nodes along a
trailing edge line spans a strip of two triangles per row; closed lines wrap around.
"""
function build_wake(g::SurfaceGeometry{T}, dm::DOFMapping{T}, V∞::SVector{3}; opts::WakeOptions=WakeOptions()) where {T}
    check_options(opts)
    wte = dm.wake_te
    n = nwake_te(wte)
    n == 0 && throw(ArgumentError("the geometry has no trailing edge to shed a wake from"))
    d = SVector{3,T}(normalize(V∞))
    span = maximum(norm(wte.x[i] - wte.x[j]) for i in 1:n, j in 1:n)
    ℓ1 = T(isnothing(opts.first_row) ? te_panel_size(g) : opts.first_row)
    rows = T[0, ℓ1]
    while rows[end] < opts.length * span
        push!(rows, rows[end] + (rows[end] - rows[end-1]) * opts.growth)
    end
    opts.ending == :semi_infinite && push!(rows, rows[end] * opts.far_factor)
    nrow = length(rows) - 1
    first_dir = opts.direction == :bisector ? wte.bisector : fill(d, n)
    any(iszero, first_dir) && throw(ArgumentError("a trailing edge node has no bisector; use direction = :freestream"))

    id(k, j) = j * n + k # wake node of trailing edge node k in row j
    x = Vector{SVector{3,T}}(undef, n * (nrow + 1))
    for j in 0:nrow, k in 1:n
        x[id(k, j)] = j == 0 ? wte.x[k] : wte.x[k] + rows[2] * first_dir[k] + (rows[j+1] - rows[2]) * d
    end
    te = repeat(1:n, nrow + 1)
    row = repeat(0:nrow, inner=n)

    tris = SVector{3,Int}[]
    filaments = Tuple{Int,Int}[]
    for (line, closed) in zip(wte.lines, wte.closed)
        m = length(line)
        segments = [(line[i], line[i+1]) for i in 1:m-1]
        closed && push!(segments, (line[m], line[1]))
        for (A, B) in segments, j in 0:nrow-1 # A → B along τ: the trailing edge is crossed from B to A, as by the lower surface
            push!(tris, SVector(id(A, j), id(A, j + 1), id(B, j)), SVector(id(A, j + 1), id(B, j + 1), id(B, j)))
        end
        if !closed # the outer sides, clockwise about the normal: upstream along the first node's, downstream along the last's
            append!(filaments, [(id(line[1], j + 1), id(line[1], j)) for j in 0:nrow-1])
            append!(filaments, [(id(line[m], j), id(line[m], j + 1)) for j in 0:nrow-1])
        end
    end
    return RigidWake(x, te, row, tris, filaments, rows, span, opts)
end

"""Flat triangle elements of the wake, in triangle order"""
wake_elements(w::RigidWake) = [FlatTriangle(w.x[t[1]], w.x[t[2]], w.x[t[3]]) for t in w.tris]

"""Doublet strength of every wake node from those of the wake's trailing edge nodes, μ_w = K q"""
wake_values(w::RigidWake, μ_w::AbstractVector) = μ_w[w.te]
