"""
    AssemblyOptions(; kwargs...)

Galerkin assembly settings (doc §6). Test weights:
- `test = :hat`: the nodal hat functions (Galerkin); `:onering`: one on the node's triangles; `:dual`: one on its
  median dual cell (each quadrature point counts for the corner with the largest barycentric coordinate).

Quadrature tiers, on the test triangle (the source triangle's integral is the analytic kernel):
- far, centroid distance above `far_ratio` times the larger diameter: `far_points`-point rule;
- near: `near_points`-point rule on sub-triangles of the test triangle, split in four while closer to the source
  triangle than `near_ratio` times their size, up to `near_depth` levels;
- touching (self, shared edge, shared vertex): rules graded toward the singular entity with `touch_levels` levels of
  ratio `touch_ratio` and Gauss orders growing from `touch_order` by one per level away from the singularity;
- contact: triangles that touch beyond their shared nodes, within `contact_rel` times the smaller diameter (across a
  T-junction, a gap or a sliver): the test triangle is split at the source triangle's corners on its edges and each
  piece is graded toward its own contact (`contact_rule`).
`kernel_far` is the kernel's own far-field switch (in source triangle diameters).
"""
Base.@kwdef struct AssemblyOptions
    test::Symbol = :hat
    far_ratio::Float64 = 3.0
    far_points::Int = 3
    near_points::Int = 7
    near_ratio::Float64 = 1.5
    near_depth::Int = 3
    touch_levels::Int = 6
    touch_ratio::Float64 = 0.15
    touch_order::Int = 3
    contact_rel::Float64 = 0.1
    kernel_far::Float64 = FAR_DEFAULT
    threaded::Bool = true
end

function check_options(o::AssemblyOptions)
    o.test in (:hat, :onering, :dual) || throw(ArgumentError("test must be :hat, :onering or :dual, not $(o.test)"))
    return o
end

tier_rules(o::AssemblyOptions) = TierRules(; far_points=o.far_points, near_points=o.near_points, levels=o.touch_levels,
    ratio=o.touch_ratio, order=o.touch_order)

"""Test function values at barycentric coordinates λ"""
@inline function test_weights(test::Symbol, λ::SVector{3,Float64})
    test === :hat && return λ
    test === :onering && return SVector(1.0, 1.0, 1.0)
    m = maximum(λ)
    is = SVector(abs(λ[1] - m) < 1e-12, abs(λ[2] - m) < 1e-12, abs(λ[3] - m) < 1e-12)
    return is ./ count(is)
end

"""∫ψ_a dA / A of the test functions over a triangle"""
test_integrals(test::Symbol) = test === :onering ? SVector(1.0, 1.0, 1.0) : SVector(1 / 3, 1 / 3, 1 / 3)

"""
    near_rule(T, S, base, ratio, depth)

Rule on test element T for a near source element S: T is split in four (at the edge midpoints) recursively while a
sub-triangle's centroid is closer to S than `ratio` times its size, up to `depth` levels; the leaves carry `base`.
"""
function near_rule(T::FlatTriangle, S::FlatTriangle, base::BaryRule, ratio::Real, depth::Integer)
    λs, ws = SVector{3,Float64}[], Float64[]
    near_pieces!(λs, ws, T, S, base, ratio, depth, CORNERS, 1.0, 0)
    return BaryRule(λs, ws)
end

"""Barycentric coordinates of a triangle's corners"""
const CORNERS = (SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))

"""The four children of the piece with barycentric corners b, split at its edge midpoints"""
function quarters(b)
    m12, m23, m31 = (b[1] + b[2]) / 2, (b[2] + b[3]) / 2, (b[3] + b[1]) / 2
    return ((b[1], m12, m31), (m12, b[2], m23), (m31, m23, b[3]), (m12, m23, m31))
end

"""Append rule r mapped onto the piece of the test triangle with barycentric corners b and area fraction frac"""
function append_piece!(λs, ws, r::BaryRule, b, frac)
    for (μ, w) in zip(r.λ, r.w)
        push!(λs, μ[1] * b[1] + μ[2] * b[2] + μ[3] * b[3])
        push!(ws, w * frac)
    end
    return nothing
end

"""Append the near rule (`near_rule`) on the piece of T with barycentric corners b and area fraction frac"""
function near_pieces!(λs, ws, T::FlatTriangle, S::FlatTriangle, base::BaryRule, ratio::Real, depth::Integer, b, frac, level)
    p1, p2, p3 = point(T, b[1]), point(T, b[2]), point(T, b[3])
    size = max(norm(p2 - p1), norm(p3 - p2), norm(p1 - p3))
    if level < depth && distance_to(S, (p1 + p2 + p3) / 3) < ratio * size
        for sub in quarters(b)
            near_pieces!(λs, ws, T, S, base, ratio, depth, sub, frac / 4, level + 1)
        end
    else
        append_piece!(λs, ws, base, b, frac)
    end
    return nothing
end

"""Distance from p to the segment from a to c, and the parameter of the closest point"""
function segment_projection(p::SVector{3}, a::SVector{3}, c::SVector{3})
    e = c - a
    t = clamp(dot(p - a, e) / dot(e, e), 0.0, 1.0)
    return norm(p - (a + t * e)), t
end

"""Parameters along the edges (corner k to k+1) of the piece with corners p at which a corner of S lies, away from the ends"""
function edge_splits(p, S::FlatTriangle, tol::Real)
    return ntuple(3) do k
        a, c = p[k], p[mod1(k + 1, 3)]
        L = norm(c - a)
        ts = Float64[]
        for q in S.x
            d, t = segment_projection(q, a, c)
            (d <= tol && t * L > tol && (1 - t) * L > tol) && push!(ts, t)
        end
        sort!(ts)
    end
end

"""Whether a corner of S lies on an edge of the triangle with corners p, away from the edge's ends"""
function splits_an_edge(p, S::FlatTriangle, tol::Real)
    for k in 1:3
        a, c = p[k], p[mod1(k + 1, 3)]
        L = norm(c - a)
        for q in S.x
            d, t = segment_projection(q, a, c)
            (d <= tol && t * L > tol && (1 - t) * L > tol) && return true
        end
    end
    return false
end

"""Largest distance from the centroid to a corner"""
reach(el::FlatTriangle) = max(norm(el.x[1] - el.centroid), norm(el.x[2] - el.centroid), norm(el.x[3] - el.centroid))

"""
    contact_rule(T, S, shared, rules, opts)

Rule on test triangle T for a source triangle S that touches it beyond the corners they share (`shared[k]`: T's
corner k is a node of S): along part of an edge across a T-junction, at a hanging node, across a gap or past a sliver.
The singular set is then not a whole edge or a corner of T, so T is split at the corners of S that lie on its edges,
and every piece carries the rule graded toward its own contact with S: an edge (`rules.edge`) when both ends and the
middle of one of its edges are on S, a corner (`rules.vertex`) when one corner is, the near rule when none is, and a
split in four when several separate corners are. A point is on S within `opts.contact_rel` times the smaller
diameter, so near contacts (a gap's width) are graded too. Returns nothing when S touches T at shared corners only.
"""
function contact_rule(T::FlatTriangle, S::FlatTriangle, shared::NTuple{3,Bool}, rules::TierRules, opts::AssemblyOptions)
    tol = opts.contact_rel * min(T.diam, S.diam)
    norm(S.centroid - T.centroid) > reach(T) + reach(S) + tol && return nothing
    extra = any(k -> !shared[k] && distance_to(S, T.x[k]) <= tol, 1:3) || splits_an_edge(T.x, S, tol)
    extra || return nothing
    λs, ws = SVector{3,Float64}[], Float64[]
    contact_pieces!(λs, ws, T, S, rules, opts, tol, CORNERS, 1.0, 0)
    return BaryRule(λs, ws)
end

function contact_pieces!(λs, ws, T::FlatTriangle, S::FlatTriangle, rules::TierRules, opts::AssemblyOptions, tol::Real, b, frac, level)
    p = (point(T, b[1]), point(T, b[2]), point(T, b[3]))
    if level < 4
        splits = edge_splits(p, S, tol)
        k = findfirst(!isempty, splits)
        if !isnothing(k) # fan from the opposite corner through the split points
            ts = [0.0; splits[k]; 1.0]
            bk, bk1, opp = b[k], b[mod1(k + 1, 3)], b[mod1(k + 2, 3)]
            for i in 1:length(ts)-1
                piece = ((1 - ts[i]) * bk + ts[i] * bk1, (1 - ts[i+1]) * bk + ts[i+1] * bk1, opp)
                contact_pieces!(λs, ws, T, S, rules, opts, tol, piece, frac * (ts[i+1] - ts[i]), level + 1)
            end
            return nothing
        end
    end
    on = map(x -> distance_to(S, x) <= tol, p)
    for k in 1:3
        k1 = mod1(k + 1, 3)
        on[k] && on[k1] && distance_to(S, (p[k] + p[k1]) / 2) <= tol && return append_piece!(λs, ws, rules.edge[k], b, frac)
    end
    nc = count(on)
    if nc == 1
        append_piece!(λs, ws, rules.vertex[findfirst(on)], b, frac)
    elseif nc == 0 || level >= 4
        near_pieces!(λs, ws, T, S, rules.near, opts.near_ratio, opts.near_depth, b, frac, 0)
    else # separate contacts at several corners
        for sub in quarters(b)
            contact_pieces!(λs, ws, T, S, rules, opts, tol, sub, frac / 4, level + 1)
        end
    end
    return nothing
end

"""
Contribution of source element S to the rows of test element T with the rule: the 3×3 block
∫ψ_a n_T·(γ_b × w_S) dA and the source part of the right hand side -∫ψ_a σ_S n_T·w_S dA (one kernel evaluation per
point). `side` selects the limit of the kernel on S's own plane (+1 for the self pair: the fluid side).
"""
@inline function pair_block(T::FlatTriangle{F}, S::FlatTriangle{F}, rule::BaryRule, σS::Real, test::Symbol, side::Integer, kernel_far::Real) where {F}
    blk = zero(SMatrix{3,3,F,9})
    rb = zero(SVector{3,F})
    for (λ, wq) in zip(rule.λ, rule.w)
        w = source_velocity(S.src, point(T, λ); side=side, far=kernel_far)
        W = wq * T.area
        v = cross(w, T.n) # n·(γ × w) = γ·(w × n)
        c = SVector(dot(v, S.γ[1]), dot(v, S.γ[2]), dot(v, S.γ[3]))
        ψ = W * test_weights(test, λ)
        blk += ψ * c'
        rb -= (σS * dot(w, T.n)) * ψ
    end
    return blk, rb
end

const TIER_NAMES = (:self, :edge, :vertex, :contact, :near, :far)

"""
Quadrature tier (index into `TIER_NAMES`) and rule for test triangle t (corner nodes `tri`) and source triangle s.
`touch_geo` holds the geometry's touch pairs, which get a deeper near rule should they escape the contact test.
"""
function pair_rule(g::SurfaceGeometry, els, t::Integer, s::Integer, tri, rules::TierRules, opts::AssemblyOptions, touch_geo)
    T, S = els[t], els[s]
    s == t && return 1, rules.self
    norm(S.centroid - T.centroid) > opts.far_ratio * max(T.diam, S.diam) && return 6, rules.far
    gs = g.tris[s]
    shared = (tri[1] in gs, tri[2] in gs, tri[3] in gs)
    contact = contact_rule(T, S, shared, rules, opts)
    isnothing(contact) || return 4, contact
    ns = count(shared)
    if ns == 2 # the shared edge, from corner k to k + 1
        return 2, rules.edge[shared[1] && shared[2] ? 1 : shared[2] && shared[3] ? 2 : 3]
    elseif ns == 1
        return 3, rules.vertex[findfirst(shared)]
    end
    depth = (t, s) in touch_geo ? opts.near_depth + 4 : opts.near_depth
    return 5, near_rule(T, S, rules.near, opts.near_ratio, depth)
end

"""
    source_strengths(g, V∞; V_body=zero)

Prescribed source strengths σ = -(V∞ - V_body)·n of every triangle (V_body = 0 for a body at rest; the hook for body
motion).
"""
source_strengths(els::AbstractVector{<:FlatTriangle}, V∞::SVector{3}; V_body::SVector{3}=zero(V∞)) = [-dot(V∞ - V_body, el.n) for el in els]

"""
    assemble(g, dm, els, σ, V∞; opts=AssemblyOptions(), external=nothing)

Galerkin system over the slots (doc §6): A (nslot × nslot), with A_ij = ∫ψ_i n·u_γ[φ_j] dA, and
b_i = -∫ψ_i n·(V∞ + u_σ + u_ext) dA, the velocities taken on the fluid side, so that A μ = b states ∫ψ_i n·u dA = 0.
The freestream part of b is exact; the source part uses the same kernel evaluations as A. `external(x)`, if given,
is a velocity field (e.g. particles) added at the test triangles' points (7-point rule). Returns A, b and statistics
(pairs and points per tier, time).
"""
function assemble(g::SurfaceGeometry, dm::DOFMapping, els::AbstractVector{<:FlatTriangle}, σ::AbstractVector, V∞::SVector{3};
    opts::AssemblyOptions=AssemblyOptions(), external=nothing)
    check_options(opts)
    nt, ns = length(els), dm.nslot
    rules = tier_rules(opts)
    ext_rule = triangle_rule(7)
    touch_geo = Set{Tuple{Int,Int}}()
    for (a, b) in g.touch_pairs
        push!(touch_geo, (a, b), (b, a))
    end
    F = eltype(els[1].n)
    At = zeros(F, ns, ns) # transposed: column i holds row i of A
    b = zeros(F, ns)
    pairs = zeros(Int, length(TIER_NAMES))
    points = zeros(Int, length(TIER_NAMES))
    lk = ReentrantLock()
    nchunk = opts.threaded ? min(nt, 8 * Threads.nthreads()) : 1
    chunks = [collect(c) for c in Iterators.partition(1:nt, cld(nt, nchunk))]
    time = @elapsed begin
        tasks = map(chunks) do chunk
            task = () -> begin
                buf = zeros(F, ns, 3)
                rows = zeros(F, 3)
                my_pairs, my_points = zeros(Int, length(TIER_NAMES)), zeros(Int, length(TIER_NAMES))
                for t in chunk
                    T = els[t]
                    tri = g.tris[t]
                    fill!(buf, 0)
                    fill!(rows, 0)
                    for s in 1:nt
                        S = els[s]
                        tier, rule = pair_rule(g, els, t, s, tri, rules, opts, touch_geo)
                        blk, rb = pair_block(T, S, rule, σ[s], opts.test, s == t ? 1 : 0, opts.kernel_far)
                        cs = dm.corner_slot[s]
                        for bb in 1:3, a in 1:3
                            buf[cs[bb], a] += blk[a, bb]
                        end
                        rows .+= rb
                        my_pairs[tier] += 1
                        my_points[tier] += length(rule)
                    end
                    # freestream (exact) and external velocity
                    rows .-= dot(V∞, T.n) * T.area .* test_integrals(opts.test)
                    if !isnothing(external)
                        for (λ, wq) in zip(ext_rule.λ, ext_rule.w)
                            rows .-= (wq * T.area * dot(external(point(T, λ)), T.n)) .* test_weights(opts.test, λ)
                        end
                    end
                    ct = dm.corner_slot[t]
                    lock(lk) do
                        for a in 1:3
                            @views At[:, ct[a]] .+= buf[:, a]
                            b[ct[a]] += rows[a]
                        end
                    end
                end
                lock(lk) do
                    pairs .+= my_pairs
                    points .+= my_points
                end
            end
            opts.threaded ? Threads.@spawn(task()) : (task(); nothing)
        end
        opts.threaded && foreach(wait, tasks)
    end
    stats = (pairs=NamedTuple{TIER_NAMES}(Tuple(pairs)), points=NamedTuple{TIER_NAMES}(Tuple(points)), time=time)
    return permutedims(At), b, stats
end
