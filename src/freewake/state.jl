"""
    FreeWakeOptions(; dt, kwargs...)

Settings of a free wake run (`FreeWakeSim`):
- `dt`: the time step. Every step sheds one free panel row, about |V∞| dt long.
- `te_len`: length of the implicit trailing edge panel (`te_panel`); `nothing` takes 0.3 |V∞| dt, about where Katz &
  Plotkin put the newest shed vortex. `te_direction = :freestream` or `:bisector`.
- `max_rows = :auto`: free panel rows kept before the oldest is converted into particles (MESA's `maxwakelen`; 0
  converts every new row at once). `:auto` keeps as many as put the new particles three cores behind the trailing
  edge (the widest strip's first particles), a distance that does not change with dt: about 10 rows at a dt of one
  trailing edge panel. Particles created within two cores of the trailing edge cost lift accuracy, so `FreeWakeSim`
  warns when a given `max_rows` is too small for its cores.
- `merge_rows = 1`: converted rows merged into one particle pair per strip. `:auto` merges, strip by strip, as many as
  make the merged quad about as long as the strip is wide (round(width / |V∞| dt)), so a small dt does not multiply the
  particles: narrow strips (the tips of a cosine spacing) keep their fine particles, wide ones are coarsened. The
  fine particles of a batch act as they are until it is complete; then each pair takes the vorticity of the fine
  triangles on its side of the merged quad's diagonal, at their strength-weighted centre (the total vorticity exactly).
- `core = :auto`: how the cores of the particles are chosen, and of the blobs that stand for wake panels and filaments
  where velocities are needed at wake points (`CoreRule`): `:auto` takes `overlap` (1.5) times half the element's
  longest edge, so neighbouring particles overlap however thin the rows; `:scaled` takes `core_factor` times the
  element's size (√(2A) of a triangle, half a filament segment), at least `core_floor` (`nothing`: |V∞| dt); a number
  is one core for every particle (MESA's `r_vortex`).
- `advection = :free`: wake nodes and particles move with the local velocity (forward Euler); `:freestream` moves them
  with V∞ (a prescribed wake).
- `delete_beyond = Inf`: particles further than this many spans downstream of the trailing edge are removed.
"""
Base.@kwdef struct FreeWakeOptions
    dt::Float64
    te_len::Union{Nothing,Float64} = nothing
    te_direction::Symbol = :freestream
    max_rows::Union{Int,Symbol} = :auto
    merge_rows::Union{Int,Symbol} = 1
    core::Union{Symbol,Float64} = :auto
    overlap::Float64 = 1.5
    core_factor::Float64 = 1.5
    core_floor::Union{Nothing,Float64} = nothing
    advection::Symbol = :free
    delete_beyond::Float64 = Inf
end

function check_options(o::FreeWakeOptions)
    o.dt > 0 || throw(ArgumentError("dt must be positive"))
    (o.max_rows === :auto || (o.max_rows isa Integer && o.max_rows >= 0)) || throw(ArgumentError("max_rows must be :auto or at least 0"))
    (o.merge_rows === :auto || (o.merge_rows isa Integer && o.merge_rows >= 1)) || throw(ArgumentError("merge_rows must be :auto or at least 1"))
    o.advection in (:free, :freestream) || throw(ArgumentError("advection must be :free or :freestream, not $(o.advection)"))
    o.te_direction in (:freestream, :bisector) || throw(ArgumentError("te_direction must be :freestream or :bisector"))
    (o.core in (:auto, :scaled) || (o.core isa Real && o.core > 0)) ||
        throw(ArgumentError("core must be :auto, :scaled or a positive number, not $(o.core)"))
    o.overlap > 0 && o.core_factor > 0 || throw(ArgumentError("overlap and core_factor must be positive"))
    return o
end

"""
    CoreRule(opts, V∞)

The core radius of each particle, and of each blob that stands for a wake panel or filament where velocities are needed
at wake points: called with a triangle's corners, `rule(p, q, r)`, or a filament segment's ends, `rule(a, b)`.

`:auto` (the default) takes `overlap` times half the element's longest edge. A conversion puts two particles in each
quad, about half the quad's longer side apart in that direction and closer in the other, so with `overlap` ≥ 1 the
cores of neighbours overlap in every direction, which a smooth sheet needs (smaller cores leave gaps and a noisy
field; much larger ones smear the tip vortices). It follows the mesh and the time step without tuning: rows thinner
than the strips are wide (a small dt) keep the cores the strips' width needs, where an area-based core would shrink
and leave spanwise gaps. Filaments take `overlap` times half a segment, the spacing of their 2-point blobs, which
matches the sheet beside them. `:scaled` is the earlier rule, `core_factor` √(2A) or half a segment, at least
`core_floor`; a number gives every particle that core.
"""
struct CoreRule{T}
    mode::Symbol
    overlap::T
    factor::T
    floor::T
    value::T
end

function CoreRule(o::FreeWakeOptions, V∞::SVector{3,T}) where {T}
    floor = T(something(o.core_floor, norm(V∞) * o.dt))
    o.core isa Real && return CoreRule(:fixed, T(o.overlap), T(o.core_factor), floor, T(o.core))
    return CoreRule(o.core, T(o.overlap), T(o.core_factor), floor, zero(T))
end

function (c::CoreRule)(p::SVector{3}, q::SVector{3}, r::SVector{3})
    c.mode === :fixed && return c.value
    c.mode === :auto && return c.overlap * max(norm(q - p), norm(r - q), norm(p - r)) / 2
    return max(c.factor * sqrt(norm(cross(q - p, r - p))), c.floor)
end

function (c::CoreRule)(a::SVector{3}, b::SVector{3})
    c.mode === :fixed && return c.value
    c.mode === :auto && return c.overlap * norm(b - a) / 2
    return max(c.factor * norm(b - a) / 2, c.floor)
end

"""
Strips of the free wake: the pairs (A, B) of neighbouring wake trailing edge nodes along τ (closed lines wrap around),
and the end nodes of open lines, whose sides are the wake's tips (the first node's runs upstream, the last's downstream,
as in `build_wake`)
"""
struct Topology
    nte::Int
    segments::Vector{Tuple{Int,Int}}
    first_nodes::Vector{Int}
    last_nodes::Vector{Int}
end

function Topology(dm::DOFMapping)
    wte = dm.wake_te
    segs, firsts, lasts = Tuple{Int,Int}[], Int[], Int[]
    for (line, closed) in zip(wte.lines, wte.closed)
        m = length(line)
        append!(segs, [(line[i], line[i+1]) for i in 1:m-1])
        if closed
            push!(segs, (line[m], line[1]))
        else
            push!(firsts, line[1])
            push!(lasts, line[m])
        end
    end
    return Topology(nwake_te(wte), segs, firsts, lasts)
end

"""
    FreeWakeState

The free wake: node rows `x` from the newest (R0, on the trailing edge panel's downstream edge) to the oldest, the
doublets `μ` each row was shed with (frozen), `μ_line` (per trailing edge node) the doublets of the last converted row's
upstream nodes (zero before the first conversion), the particles, and the number of conversions so far. Panel row i
lies between node rows i and i + 1. With merging (`merge_rows`), `tag` marks each particle 0 when final or with the
strip whose batch it awaits, `group` the side of the merged quad's diagonal it goes to, and per strip `batch_rows`,
`batch_length` and `batch_width` describe the batch so far.
"""
mutable struct FreeWakeState{T}
    x::Vector{Vector{SVector{3,T}}}
    μ::Vector{Vector{T}}
    μ_line::Vector{T}
    particles::Vector{Blob{T}}
    conversions::Int
    tag::Vector{Int}
    group::Vector{Int8}
    batch_rows::Vector{Int}
    batch_length::Vector{T}
    batch_width::Vector{T}
end

FreeWakeState{T}(nte::Integer, nstrips::Integer=0) where {T} = FreeWakeState(Vector{SVector{3,T}}[], Vector{T}[], zeros(T, nte), Blob{T}[], 0,
    Int[], Int8[], zeros(Int, nstrips), zeros(T, nstrips), zeros(T, nstrips))

"""Keep the particles (and their merge tags) where `keep` is true"""
function keep_particles!(st::FreeWakeState, keep::AbstractVector{Bool})
    st.particles = st.particles[keep]
    st.tag = st.tag[keep]
    st.group = st.group[keep]
    return st
end

"""
Vortex elements of known strength: triangles (their corners and corner doublets; the ones with an area also as
`FlatTriangle`s for the right hand side) and filament segments with their end circulations
"""
struct WakeElements{T}
    corners::Vector{NTuple{3,SVector{3,T}}}
    corner_μ::Vector{SVector{3,T}}
    sheets::Vector{FlatTriangle{T}}
    sheet_μ::Vector{SVector{3,T}}
    segments::Vector{Segment{T}}
    segment_Γ::Vector{SVector{2,T}}
end

WakeElements{T}() where {T} =
    WakeElements(NTuple{3,SVector{3,T}}[], SVector{3,T}[], FlatTriangle{T}[], SVector{3,T}[], Segment{T}[], SVector{2,T}[])

"""Whether the triangle p, q, r has an area (relative to its longest edge)"""
nondegenerate(p, q, r) = norm(cross(q - p, r - p)) > 1e-12 * max(sum(abs2, q - p), sum(abs2, r - p), sum(abs2, r - q))

"""Add a sheet triangle; a collapsed one keeps its vorticity in the blobs but is left out of the right hand side"""
function add_sheet!(e::WakeElements, p, q, r, mp, mq, mr)
    μ = SVector(mp, mq, mr)
    push!(e.corners, (p, q, r))
    push!(e.corner_μ, μ)
    if nondegenerate(p, q, r)
        push!(e.sheets, FlatTriangle(p, q, r))
        push!(e.sheet_μ, μ)
    end
    return e
end

add_segment!(e::WakeElements, a, b, Γa, Γb) = (push!(e.segments, Segment(a, b)); push!(e.segment_Γ, SVector(Γa, Γb)); e)

"""Append the sheets (two triangles per strip, as `build_wake`) and tip sides of the panel row between node rows i and i + 1"""
function panel_row!(e::WakeElements, st::FreeWakeState, topo::Topology, i::Integer)
    x0, x1, μ0, μ1 = st.x[i], st.x[i+1], st.μ[i], st.μ[i+1]
    for (A, B) in topo.segments
        add_sheet!(e, x0[A], x1[A], x0[B], μ0[A], μ1[A], μ0[B])
        add_sheet!(e, x1[A], x1[B], x0[B], μ1[A], μ1[B], μ0[B])
    end
    for k in topo.first_nodes # upstream
        add_segment!(e, x1[k], x0[k], μ1[k], μ0[k])
    end
    for k in topo.last_nodes # downstream
        add_segment!(e, x0[k], x1[k], μ0[k], μ1[k])
    end
    return e
end

"""Append the far end filament along node row i (B → A, strength μ - μ_line), unless it carries nothing"""
function far_end!(e::WakeElements, st::FreeWakeState, topo::Topology, i::Integer)
    x, μ = st.x[i], st.μ[i]
    for (A, B) in topo.segments
        ΓA, ΓB = μ[A] - st.μ_line[A], μ[B] - st.μ_line[B]
        iszero(ΓA) && iszero(ΓB) && continue
        add_segment!(e, x[B], x[A], ΓB, ΓA)
    end
    return e
end

"""
    vortex_elements(st, topo)

The free wake's vortex elements: the sheets and tip sides of every panel row, the interface filament along R0 (A → B,
strength μ(R0)), which cancels the trailing edge panel's downstream filament in steady flow, and the far end filament
along the oldest row (B → A, strength μ - μ_line: the starting vortex until the first conversion, nothing after it).
With a single node row (`max_rows = 0`) the two filaments leave MESA's line of strength μ_line along R0.
"""
function vortex_elements(st::FreeWakeState{T}, topo::Topology) where {T}
    e = WakeElements{T}()
    nr = length(st.x)
    for i in 1:nr-1
        panel_row!(e, st, topo, i)
    end
    if nr >= 1
        for (A, B) in topo.segments
            add_segment!(e, st.x[1][A], st.x[1][B], st.μ[1][A], st.μ[1][B])
        end
        far_end!(e, st, topo, nr)
    end
    return e
end

"""Vortex elements of a rigid wake (the trailing edge panel) with the doublets μ_w of its trailing edge nodes"""
function vortex_elements(w::RigidWake{T}, μ_w::AbstractVector) where {T}
    e = WakeElements{T}()
    μn = μ_w[w.te]
    for t in w.tris
        add_sheet!(e, w.x[t[1]], w.x[t[2]], w.x[t[3]], μn[t[1]], μn[t[2]], μn[t[3]])
    end
    for (a, b) in w.filaments
        add_segment!(e, w.x[a], w.x[b], μn[a], μn[b])
    end
    return e
end

ExplicitSources(e::WakeElements, blobs) = ExplicitSources(e.sheets, e.sheet_μ, e.segments, e.segment_Γ, blobs)

"""
Vorticity of a triangle's sheet, γA = ∮μ dl = Σ_edges ½(μ_a + μ_b)(x_b - x_a): computed without the area, so it stays
finite on a collapsed triangle
"""
sheet_vorticity(p, q, r, μ::SVector{3}) = ((μ[1] + μ[2]) * (q - p) + (μ[2] + μ[3]) * (r - q) + (μ[3] + μ[1]) * (p - r)) / 2

"""
    free_wake_blobs(e, core)

Blob representation of vortex elements, for velocities at wake points: a particle at each triangle's centroid carrying
its sheet's vorticity (the particle a conversion creates), and two at the 2-point Gauss positions of each filament
(`segment_blobs`), with the cores `core(p, q, r)` and `core(a, b)` of a `CoreRule`. The particles and the 2-point blobs
carry the elements' total vorticity and its first moment exactly.
"""
function free_wake_blobs(e::WakeElements{T}, core) where {T}
    blobs = Blob{T}[]
    for ((p, q, r), μ) in zip(e.corners, e.corner_μ)
        α = sheet_vorticity(p, q, r, μ)
        push!(blobs, Blob((p + q + r) / 3, α, T(core(p, q, r))))
    end
    for (f, Γ) in zip(e.segments, e.segment_Γ)
        append!(blobs, segment_blobs(f.a, f.b, Γ[1], Γ[2], core(f.a, f.b)))
    end
    return blobs
end

"""
    convert_oldest!(st, topo, core, merge=nothing)

Convert the oldest panel row into particles: one at each triangle's centroid with its sheet's vorticity, two
(2-point Gauss) for each tip side and for the far edge (strength μ - μ_line: the starting vortex at the first
conversion, nothing later); then μ_line takes the row's upstream doublets and its downstream node row is dropped. The
total vorticity and its first moment are unchanged, and the remaining wake's far edge carries nothing.

With `merge` (rows per particle pair for each strip, `merge_rows`), a strip's sheet particles join its batch, tagged
with the side of the merged quad's diagonal their triangle's centroid falls on; when the batch holds `merge[s]` rows
its particles are merged into one pair (`merge_strip!`). Returns the number of particles created.
"""
function convert_oldest!(st::FreeWakeState{T}, topo::Topology, core, merge=nothing) where {T}
    nr = length(st.x)
    nr >= 2 || return 0
    e = WakeElements{T}()
    panel_row!(e, st, topo, nr - 1)
    far_end!(e, st, topo, nr)
    blobs = free_wake_blobs(e, core) # the strips' two sheet particles first, in strip order, then the filaments'
    tag, group = zeros(Int, length(blobs)), zeros(Int8, length(blobs))
    if !isnothing(merge)
        up, down = st.x[nr-1], st.x[nr]
        for (s, (A, B)) in enumerate(topo.segments)
            k = merge[s]
            k > 1 || continue
            m = k - st.batch_rows[s] # this row's place in the batch from its upstream edge (the first converted is the last)
            for (h, (u, v)) in enumerate(((1 / 3, 1 / 3), (2 / 3, 2 / 3))) # the two triangles' centroids in their quad
                tag[2s-2+h] = s
                group[2s-2+h] = u + (m - 1 + v) / k < 1 ? 1 : 2 # the merged quad's diagonal runs as the fine ones'
            end
            st.batch_rows[s] += 1
            st.batch_length[s] += (norm(down[A] - up[A]) + norm(down[B] - up[B])) / 2
            st.batch_width[s] = norm(up[B] - up[A])
        end
    end
    append!(st.particles, blobs)
    append!(st.tag, tag)
    append!(st.group, group)
    st.μ_line = copy(st.μ[nr-1])
    pop!(st.x)
    pop!(st.μ)
    st.conversions += 1
    if !isnothing(merge)
        for s in eachindex(topo.segments)
            merge[s] > 1 && st.batch_rows[s] >= merge[s] && merge_strip!(st, s, core)
        end
    end
    return length(blobs)
end

"""
    merge_strip!(st, s, core)

Merge strip s's batch of particles into one pair, one for each side of the merged quad's diagonal: the sum of their
strengths (the total vorticity exactly) at their strength-weighted centre (the impulse to within the batch's spread),
with the core of the merged quad (`core` of a triangle as wide as the strip and as long as the batch)
"""
function merge_strip!(st::FreeWakeState{T}, s::Integer, core) where {T}
    idx = findall(==(s), st.tag)
    o = zero(SVector{3,T})
    δ = T(core(o, o + SVector{3,T}(st.batch_width[s], 0, 0), o + SVector{3,T}(st.batch_width[s], st.batch_length[s], 0)))
    merged = Blob{T}[]
    for g in (1, 2)
        sel = [i for i in idx if st.group[i] == g]
        isempty(sel) && continue
        α = sum(st.particles[i].α for i in sel)
        w = [norm(st.particles[i].α) for i in sel]
        x = sum(w) > 0 ? sum(w[j] * st.particles[sel[j]].x for j in eachindex(sel)) / sum(w) : sum(st.particles[i].x for i in sel) / length(sel)
        push!(merged, Blob(x, α, δ))
    end
    keep_particles!(st, st.tag .!= s)
    append!(st.particles, merged)
    append!(st.tag, zeros(Int, length(merged)))
    append!(st.group, zeros(Int8, length(merged)))
    st.batch_rows[s] = 0
    st.batch_length[s] = zero(T)
    return st
end

"""
    vorticity_moments(blobs)

Total vorticity Ω = Σα and impulse I = ½Σ x × α of a set of particles (or of the blob representation of any vortex
elements, which carries both exactly)
"""
function vorticity_moments(blobs::AbstractVector{Blob{T}}) where {T}
    Ω = sum((b.α for b in blobs); init=zero(SVector{3,T}))
    I = sum((cross(b.x, b.α) for b in blobs); init=zero(SVector{3,T})) / 2
    return Ω, I
end
