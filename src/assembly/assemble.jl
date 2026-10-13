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
  ratio `touch_ratio` and Gauss orders growing from `touch_order` by one per level away from the singularity. A shared
  corner wider than 90° is split by rays from it into equal angles of at most `corner_angle` degrees, each piece
  carrying the vertex rule (`vertex_pieces!`);
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
    corner_angle::Float64 = 75.0
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
function near_rule(T::FlatTriangle, S, base::BaryRule, ratio::Real, depth::Integer)
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

"""
    vertex_pieces!(λs, ws, p, k, base, b, frac, corner_angle)

Append the vertex rule `base` (graded toward its corner 1) toward corner k of the piece with physical corners p,
barycentric corners b and area fraction frac. The rule's angular coordinate runs linearly along the opposite edge,
which leaves the middle of a wide corner under-sampled (5e-4 relative error at 150°, 2e-6 at 120°, 5e-9 at 90°), so
a corner wider than 90° is split by rays from it into equal angles of at most `corner_angle` degrees. The rays divide
the opposite edge in the ratio of the areas, |p_i - p_k| sin φ : |p_j - p_k| sin(θ - φ).
"""
function vertex_pieces!(λs, ws, p, k::Integer, base::BaryRule, b, frac, corner_angle::Real)
    i, j = mod1(k + 1, 3), mod1(k + 2, 3)
    u, v = p[i] - p[k], p[j] - p[k]
    θ = acos(clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0))
    if θ <= π / 2
        return append_piece!(λs, ws, base, (b[k], b[i], b[j]), frac)
    end
    m = ceil(Int, rad2deg(θ) / corner_angle)
    ts = [0.0; [(φ = θ * q / m; norm(u) * sin(φ) / (norm(u) * sin(φ) + norm(v) * sin(θ - φ))) for q in 1:m-1]; 1.0]
    for q in 1:m
        b1, b2 = (1 - ts[q]) * b[i] + ts[q] * b[j], (1 - ts[q+1]) * b[i] + ts[q+1] * b[j]
        append_piece!(λs, ws, base, (b[k], b1, b2), frac * (ts[q+1] - ts[q]))
    end
    return nothing
end

"""Vertex rule toward corner k of T (the stored rotation unless the corner is wider than 90°)"""
function vertex_rule(T::FlatTriangle, k::Integer, rules::TierRules, opts::AssemblyOptions)
    u, v = T.x[mod1(k + 1, 3)] - T.x[k], T.x[mod1(k + 2, 3)] - T.x[k]
    dot(u, v) >= 0 && return rules.vertex[k]
    λs, ws = SVector{3,Float64}[], Float64[]
    vertex_pieces!(λs, ws, T.x, k, rules.vertex[1], CORNERS, 1.0, opts.corner_angle)
    return BaryRule(λs, ws)
end

"""
Append the near rule (`near_rule`) on the piece of T with barycentric corners b and area fraction frac; the source S is
anything with a `distance_to` (a triangle or a filament segment)
"""
function near_pieces!(λs, ws, T::FlatTriangle, S, base::BaryRule, ratio::Real, depth::Integer, b, frac, level)
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

"""A straight filament from a to b (a wake tip filament) as a source of the assembly"""
struct Segment{T}
    a::SVector{3,T}
    b::SVector{3,T}
end

Elements.distance_to(s::Segment, p::SVector{3}) = first(segment_projection(p, s.a, s.b))

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
and every piece carries the rule graded toward its own contact with S: an edge (`rules.edge`) when exactly two corners
and the middle of the edge between them are on S, a corner (`vertex_pieces!`) when one corner is, the near rule when
none is, and otherwise (two corners not joined by a contact edge, or all three, where an edge rule would miss the
singularity at the third) the self rule, graded toward all three edges and their ends. A point is on S within
`opts.contact_rel` times the smaller diameter, so near contacts (a gap's width) are graded too. Returns nothing when S
touches T at shared corners only.
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
    nc = count(on)
    if nc == 2 # the edge from corner k to k + 1 joins the two corners on S
        k = !on[3] ? 1 : !on[1] ? 2 : 3
        distance_to(S, (p[k] + p[mod1(k + 1, 3)]) / 2) <= tol && return append_piece!(λs, ws, rules.edge[k], b, frac)
    end
    if nc == 1
        vertex_pieces!(λs, ws, p, findfirst(on), rules.vertex[1], b, frac, opts.corner_angle)
    elseif nc == 0
        near_pieces!(λs, ws, T, S, rules.near, opts.near_ratio, opts.near_depth, b, frac, 0)
    else # contacts at several corners, not along one edge
        append_piece!(λs, ws, rules.self, b, frac)
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
        return 3, vertex_rule(T, findfirst(shared), rules, opts)
    end
    depth = (t, s) in touch_geo ? opts.near_depth + 4 : opts.near_depth
    return 5, near_rule(T, S, rules.near, opts.near_ratio, depth)
end

"""
Tier and rule for body test triangle T and wake triangle S. They share no mesh node, so the wake's first row, which
meets the body along the trailing edge, is found by the contact test (an edge or a corner of T on S).
"""
function wake_pair_rule(T::FlatTriangle, S::FlatTriangle, rules::TierRules, opts::AssemblyOptions)
    norm(S.centroid - T.centroid) > opts.far_ratio * max(T.diam, S.diam) && return 6, rules.far
    contact = contact_rule(T, S, (false, false, false), rules, opts)
    isnothing(contact) || return 4, contact
    return 5, near_rule(T, S, rules.near, opts.near_ratio, opts.near_depth)
end

"""
Tier and rule for body test triangle T and filament segment f: the vertex rule where the filament starts at a corner
of T (a tip filament at the trailing edge's end node, whose 1/r singularity there is integrable), else far or near
"""
function segment_rule(T::FlatTriangle, f::Segment, rules::TierRules, opts::AssemblyOptions)
    L = norm(f.b - f.a)
    norm((f.a + f.b) / 2 - T.centroid) > opts.far_ratio * max(T.diam, L) && return 6, rules.far
    for k in 1:3
        min(norm(T.x[k] - f.a), norm(T.x[k] - f.b)) <= 1e-12 * T.diam && return 3, vertex_rule(T, k, rules, opts)
    end
    return 5, near_rule(T, f, rules.near, opts.near_ratio, opts.near_depth)
end

"""
Contribution of filament segment f, of linearly varying strength, to the rows of test element T with the rule:
∫ψ_a n_T·v_a dA and ∫ψ_a n_T·v_b dA, the parts per unit strength at the segment's two ends (`segment_velocity`)
"""
@inline function segment_block(T::FlatTriangle{F}, f::Segment{F}, rule::BaryRule, test::Symbol) where {F}
    ca, cb = zero(SVector{3,F}), zero(SVector{3,F})
    for (λ, wq) in zip(rule.λ, rule.w)
        va, vb = segment_velocity(f.a, f.b, point(T, λ))
        ψ = (wq * T.area) * test_weights(test, λ)
        ca += dot(va, T.n) * ψ
        cb += dot(vb, T.n) * ψ
    end
    return ca, cb
end

"""The wake as sources of the assembly: its triangles (elements and wake nodes) and tip filament segments"""
function wake_sources(w::RigidWake, ::Type{F}) where {F}
    return (elements=[FlatTriangle(SVector{3,F}.(w.x[t])...) for t in w.tris], tris=w.tris,
        segments=[Segment(SVector{3,F}(w.x[a]), SVector{3,F}(w.x[b])) for (a, b) in w.filaments], filaments=w.filaments)
end

"""
Add the wake's columns (after the `ns` slots) of test triangle T's rows to `buf`: the sheets of the wake triangles and
the tip filaments, per unit doublet strength at each wake node
"""
function add_wake!(buf, pairs, points, T::FlatTriangle{F}, ns::Integer, w, rules::TierRules, opts::AssemblyOptions) where {F}
    for (s, S) in enumerate(w.elements)
        tier, rule = wake_pair_rule(T, S, rules, opts)
        blk, _ = pair_block(T, S, rule, zero(F), opts.test, 0, opts.kernel_far)
        c = w.tris[s]
        for bb in 1:3, a in 1:3
            buf[ns+c[bb], a] += blk[a, bb]
        end
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    for (f, seg) in enumerate(w.segments)
        tier, rule = segment_rule(T, seg, rules, opts)
        ca, cb = segment_block(T, seg, rule, opts.test)
        ia, ib = w.filaments[f]
        for a in 1:3
            buf[ns+ia, a] += ca[a]
            buf[ns+ib, a] += cb[a]
        end
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    return nothing
end

"""
    source_strengths(g, V∞; V_body=zero)

Prescribed source strengths σ = -(V∞ - V_body)·n of every triangle (V_body = 0 for a body at rest; the hook for body
motion).
"""
source_strengths(els::AbstractVector{<:FlatTriangle}, V∞::SVector{3}; V_body::SVector{3}=zero(V∞)) = [-dot(V∞ - V_body, el.n) for el in els]

"""
    assemble(g, dm, els, σ, V∞; opts=AssemblyOptions(), external=nothing, wake=nothing)

Galerkin system over the slots (doc §6): A_ij = ∫ψ_i n·u[φ_j] dA and b_i = -∫ψ_i n·(V∞ + u_σ + u_ext) dA, the
velocities taken on the fluid side, so that A μ = b states ∫ψ_i n·u dA = 0. The rows are the body's slots (the wake is
never tested); the columns are the slots (their vortex sheets) followed, with a `wake`, by the wake's nodes (their
sheets on the wake triangles and the tip filaments), which `wake_rows(dm, wake.te)` folds back onto the unknowns
(`solve_galerkin(...; W)`). The freestream part of b is exact; the source part uses the same kernel evaluations as A.
`external(x)`, if given, is a velocity field (e.g. particles) added at the test triangles' points (7-point rule).
Returns A, b and statistics (pairs and points per tier, time).
"""
function assemble(g::SurfaceGeometry, dm::DOFMapping, els::AbstractVector{<:FlatTriangle}, σ::AbstractVector, V∞::SVector{3};
    opts::AssemblyOptions=AssemblyOptions(), external=nothing, wake::Union{Nothing,RigidWake}=nothing)
    check_options(opts)
    nt, ns = length(els), dm.nslot
    rules = tier_rules(opts)
    ext_rule = triangle_rule(7)
    touch_geo = Set{Tuple{Int,Int}}()
    for (a, b) in g.touch_pairs
        push!(touch_geo, (a, b), (b, a))
    end
    F = eltype(els[1].n)
    nw = isnothing(wake) ? 0 : nwake(wake)
    wsrc = isnothing(wake) ? nothing : wake_sources(wake, F)
    At = zeros(F, ns + nw, ns) # transposed: column i holds row i of A
    b = zeros(F, ns)
    pairs = zeros(Int, length(TIER_NAMES))
    points = zeros(Int, length(TIER_NAMES))
    lk = ReentrantLock()
    nchunk = opts.threaded ? min(nt, 8 * Threads.nthreads()) : 1
    chunks = [collect(c) for c in Iterators.partition(1:nt, cld(nt, nchunk))]
    time = @elapsed begin
        tasks = map(chunks) do chunk
            task = () -> begin
                buf = zeros(F, ns + nw, 3)
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
                    isnothing(wsrc) || add_wake!(buf, my_pairs, my_points, T, ns, wsrc, rules, opts)
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
