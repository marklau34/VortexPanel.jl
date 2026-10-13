# Galerkin assembly over mixed meshes: the input's quadrilaterals kept whole as isoparametric quads, every other face
# (tip caps, input triangles) a flat triangle. Validation stage 7 compares it with the split triangles of `assemble`.

"""
    MixedMesh(g, dm)

The elements of geometry g with its input quadrilaterals kept whole: an `IsoQuad` for every welded quad (corners in its
outward order), a `FlatTriangle` for every other triangle (caps, input triangles), each with its corner nodes and DOF
slots. A quad's slots come from its two triangles, so the trailing edge copies follow theirs. `quad_tris` holds the
geometry's two triangles of each quad and `tri_ids` the geometry triangle of each triangle; `element_of` maps every
geometry triangle to its element (quads first, then triangles).
"""
struct MixedMesh{T}
    quads::Vector{IsoQuad{T}}
    quad_nodes::Vector{SVector{4,Int}}
    quad_slots::Vector{SVector{4,Int}}
    quad_tris::Vector{NTuple{2,Int}}
    tris::Vector{FlatTriangle{T}}
    tri_nodes::Vector{SVector{3,Int}}
    tri_slots::Vector{SVector{3,Int}}
    tri_ids::Vector{Int}
    element_of::Vector{Int}
end

function MixedMesh(g::SurfaceGeometry{T}, dm::DOFMapping) where {T}
    halves = Dict{Int,Vector{Int}}()
    for t in eachindex(g.tris)
        g.tri_slot[t] > 0 && push!(get!(halves, g.tri_parent[t], Int[]), t)
    end
    quads, qnodes, qslots, qtris = IsoQuad{T}[], SVector{4,Int}[], SVector{4,Int}[], NTuple{2,Int}[]
    element_of = zeros(Int, length(g.tris))
    for (k, q) in enumerate(g.quads)
        ts = get(halves, g.quad_parent[k], Int[])
        length(ts) == 2 || continue
        function slot(v)
            for t in ts
                i = findfirst(==(v), g.tris[t])
                isnothing(i) || return dm.corner_slot[t][i]
            end
            error("node $v of quad $k is on neither of its triangles")
        end
        push!(quads, IsoQuad(g.x[q[1]], g.x[q[2]], g.x[q[3]], g.x[q[4]]))
        push!(qnodes, SVector{4,Int}(q[1], q[2], q[3], q[4]))
        push!(qslots, SVector(slot(q[1]), slot(q[2]), slot(q[3]), slot(q[4])))
        push!(qtris, (ts[1], ts[2]))
        element_of[ts] .= length(quads)
    end
    tris, tnodes, tslots, tids = FlatTriangle{T}[], SVector{3,Int}[], SVector{3,Int}[], Int[]
    for t in eachindex(g.tris)
        element_of[t] == 0 || continue
        tri = g.tris[t]
        push!(tris, FlatTriangle(g.x[tri[1]], g.x[tri[2]], g.x[tri[3]]))
        push!(tnodes, SVector{3,Int}(tri[1], tri[2], tri[3]))
        push!(tslots, SVector{3,Int}(dm.corner_slot[t]...))
        push!(tids, t)
        element_of[t] = length(quads) + length(tris)
    end
    return MixedMesh(quads, qnodes, qslots, qtris, tris, tnodes, tslots, tids, element_of)
end

Base.show(io::IO, m::MixedMesh) = print(io, "MixedMesh($(length(m.quads)) quads, $(length(m.tris)) triangles)")

"""Test function values at the parameters p of a quad (see `test_weights`): the quarter of the square around a corner is its dual cell"""
@inline function square_test_weights(test::Symbol, p::SVector{2})
    test === :hat && return shape(p[1], p[2])
    test === :onering && return SVector(1.0, 1.0, 1.0, 1.0)
    k = p[1] < 0.5 ? (p[2] < 0.5 ? 1 : 4) : (p[2] < 0.5 ? 2 : 3)
    return SVector(ntuple(i -> Float64(i == k), 4))
end

"""The parameter-square tier rules of the assembly options (tensor Gauss far and near, graded touching rules)"""
square_tiers(o::AssemblyOptions) = SquareTiers(; far_points=o.far_points <= 3 ? 2 : 3, near_points=o.near_points <= 7 ? 3 : 4,
    levels=o.quad_touch_levels, ratio=o.touch_ratio, order=o.touch_order)

"""Append the square rule r mapped onto the parameter rectangle [lo, hi]"""
function append_rect!(ps, ws, r::SquareRule, lo::SVector{2}, hi::SVector{2})
    d = hi - lo
    for (q, w) in zip(r.p, r.w)
        push!(ps, lo + d .* q)
        push!(ws, w * d[1] * d[2])
    end
    return nothing
end

"""Append `base` on the rectangle [lo, hi] of T's parameters, split in four while close to S (`square_near_rule`)"""
function near_rects!(ps, ws, T::IsoQuad, S, base::SquareRule, ratio::Real, depth::Integer, lo::SVector{2}, hi::SVector{2}, level::Integer)
    corners = (point(T, lo), point(T, SVector(hi[1], lo[2])), point(T, hi), point(T, SVector(lo[1], hi[2])))
    extent = max(norm(corners[3] - corners[1]), norm(corners[4] - corners[2]))
    if level < depth && distance_to(S, point(T, (lo + hi) / 2)) < ratio * extent
        m = (lo + hi) / 2
        for (a, b) in ((lo, m), (SVector(m[1], lo[2]), SVector(hi[1], m[2])), (SVector(lo[1], m[2]), SVector(m[1], hi[2])), (m, hi))
            near_rects!(ps, ws, T, S, base, ratio, depth, a, b, level + 1)
        end
    else
        append_rect!(ps, ws, base, lo, hi)
    end
    return nothing
end

"""Point where the ray from c along d leaves the parameter rectangle [lo, hi] (c inside or on its boundary, d pointing in)"""
function ray_exit(c::SVector{2}, d::SVector{2}, lo::SVector{2}, hi::SVector{2})
    t = Inf
    for i in 1:2
        d[i] > 1e-14 && (t = min(t, (hi[i] - c[i]) / d[i]))
        d[i] < -1e-14 && (t = min(t, (lo[i] - c[i]) / d[i]))
    end
    return c + t * d
end

"""
    vertex_rect!(ps, ws, T, k, lo, hi, base, corner_angle)

Append the vertex rule toward corner k of the parameter rectangle [lo, hi] of quad T (corners counterclockwise from
lo): the `vertex_graded_rule` `base` on triangles from the corner. A shear can make the corner far wider than its 90° in
parameters, which the rule's angular coordinate under-samples (as on triangles, `vertex_pieces!`), so the corner is
split into equal physical angles of at most `corner_angle` degrees: the rays are mapped back to the parameters
through the tangents at the corner and run to the far sides, the far corner splitting a piece that straddles it.
"""
function vertex_rect!(ps, ws, T::IsoQuad, k::Integer, lo::SVector{2}, hi::SVector{2}, base::BaryRule, corner_angle::Real)
    rc = (lo, SVector(hi[1], lo[2]), hi, SVector(lo[1], hi[2]))
    C, A, O, B = rc[k], rc[mod1(k + 1, 4)], rc[mod1(k + 2, 4)], rc[mod1(k + 3, 4)]
    a1, a2 = tangents(T, C)
    phys(δ) = δ[1] * a1 + δ[2] * a2
    u, v = phys(A - C), phys(B - C)
    θ = acos(clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0))
    m = θ <= π / 2 ? 1 : ceil(Int, rad2deg(θ) / corner_angle)
    û = u / norm(u)
    ŵ = normalize(v - dot(v, û) * û)
    G = SMatrix{2,2}(dot(a1, a1), dot(a1, a2), dot(a1, a2), dot(a2, a2))
    points = map(0:m) do i
        i == 0 && return A
        i == m && return B
        d = cos(i * θ / m) * û + sin(i * θ / m) * ŵ
        ray_exit(C, G \ SVector(dot(a1, d), dot(a2, d)), lo, hi)
    end
    # position along the far boundary A → O → B: s in [0, 1] on A-O, in [1, 2] on O-B
    tol = 1e-12 * norm(hi - lo)
    along(P, X, Y) = clamp(dot(P - X, Y - X) / dot(Y - X, Y - X), 0.0, 1.0)
    s(P) = abs(cross2(P - A, O - A)) / norm(O - A) <= tol ? along(P, A, O) : 1 + along(P, O, B)
    tri!(P, Q) = Elements.append_triangle!(ps, ws, base, C, P, Q, abs(cross2(P - C, Q - C)) / 2)
    for i in 1:m
        P, Q = points[i], points[i+1]
        if s(P) < 1 - 1e-9 && s(Q) > 1 + 1e-9 # the piece straddles the far corner
            tri!(P, O)
            tri!(O, Q)
        else
            tri!(P, Q)
        end
    end
    return nothing
end

"""The vertex rule toward corner k of quad T: the stored one, or split by `vertex_rect!` when the corner is wider than 90°"""
function vertex_square(T::IsoQuad, k::Integer, tiers::SquareTiers, opts::AssemblyOptions)
    u, v = T.x[mod1(k + 1, 4)] - T.x[k], T.x[mod1(k + 3, 4)] - T.x[k]
    dot(u, v) >= 0 && return tiers.vertex[k]
    ps, ws = SVector{2,Float64}[], Float64[]
    vertex_rect!(ps, ws, T, k, SVector(0.0, 0.0), SVector(1.0, 1.0), tiers.vertex_base, opts.corner_angle)
    return SquareRule(ps, ws)
end

"""z component of the cross product of two parameter vectors"""
cross2(a::SVector{2}, b::SVector{2}) = a[1] * b[2] - a[2] * b[1]

"""Parameters along each side k of T (corner k to k + 1) at which a corner of S lies within tol, away from the side's ends"""
function side_splits(T::IsoQuad, S, tol::Real)
    return ntuple(4) do k
        a, c = T.x[k], T.x[mod1(k + 1, 4)]
        L = norm(c - a)
        ts = Float64[]
        for q in S.x
            d, t = segment_projection(q, a, c)
            (d <= tol && t * L > tol && (1 - t) * L > tol) && push!(ts, t)
        end
        ts
    end
end

"""Sorted cuts with near-duplicates (closer than 1e-9) merged"""
function cuts(v::Vector{Float64})
    sort!(v)
    out = [v[1]]
    for x in v[2:end]
        x - out[end] > 1e-9 && push!(out, x)
    end
    return out
end

"""
    square_contact_rule(T, S, shared, tiers, opts)

Rule on test quad T for a source S that touches it beyond the corners they share (`shared[k]`: T's corner k is a node
of S): across a T-junction, a gap or a sliver, or a wake triangle along the trailing edge. The parameter square is cut
by lines through the corners of S that lie on T's sides, and every rectangle carries the rule graded toward its own
contact with S, as for triangles (`contact_rule`): a side when its two corners and its middle are on S, a corner when
one corner is, the near rule when none is, and otherwise the self rule. A point is on S within `opts.contact_rel` times
the smaller diameter. Returns nothing when S touches T at shared corners only.
"""
function square_contact_rule(T::IsoQuad, S, shared::NTuple{4,Bool}, tiers::SquareTiers, opts::AssemblyOptions)
    tol = opts.contact_rel * min(T.diam, S.diam)
    norm(S.centroid - T.centroid) > reach(T) + reach(S) + tol && return nothing
    splits = side_splits(T, S, tol)
    extra = any(k -> !shared[k] && distance_to(S, T.x[k]) <= tol, 1:4) || any(!isempty, splits)
    extra || return nothing
    ξc = cuts([0.0; splits[1]; 1 .- splits[3]; 1.0]) # side 1 runs along +ξ at η = 0, side 3 along -ξ at η = 1
    ηc = cuts([0.0; splits[2]; 1 .- splits[4]; 1.0]) # side 2 along +η at ξ = 1, side 4 along -η at ξ = 0
    ps, ws = SVector{2,Float64}[], Float64[]
    for i in 1:length(ξc)-1, j in 1:length(ηc)-1
        lo, hi = SVector(ξc[i], ηc[j]), SVector(ξc[i+1], ηc[j+1])
        c = (lo, SVector(hi[1], lo[2]), hi, SVector(lo[1], hi[2]))
        on = map(p -> distance_to(S, point(T, p)) <= tol, c)
        nc = count(on)
        k = findfirst(k -> on[k] && on[mod1(k + 1, 4)], 1:4)
        if nc == 2 && !isnothing(k) && distance_to(S, point(T, (c[k] + c[mod1(k + 1, 4)]) / 2)) <= tol
            append_rect!(ps, ws, tiers.edge[k], lo, hi)
        elseif nc == 1
            vertex_rect!(ps, ws, T, findfirst(on), lo, hi, tiers.vertex_base, opts.corner_angle)
        elseif nc == 0
            near_rects!(ps, ws, T, S, tiers.near, opts.near_ratio, opts.near_depth + 2, lo, hi, 0)
        else
            append_rect!(ps, ws, tiers.self, lo, hi)
        end
    end
    return SquareRule(ps, ws)
end

"""
Quadrature tier and rule on test quad T (corner nodes `tn`) for a body source S (corner nodes `sn`), as `pair_rule`
does for triangles; `touch` marks a geometric touch pair, which gets a deeper near rule should it escape the contact
test
"""
function quad_pair_rule(T::IsoQuad, S, tn::SVector{4,Int}, sn, self::Bool, tiers::SquareTiers, opts::AssemblyOptions, touch::Bool)
    self && return 1, tiers.self
    norm(S.centroid - T.centroid) > opts.far_ratio * max(T.diam, S.diam) && return 6, tiers.far
    shared = ntuple(k -> tn[k] in sn, 4)
    contact = square_contact_rule(T, S, shared, tiers, opts)
    isnothing(contact) || return 4, contact
    ns = count(shared)
    if ns == 2
        k = findfirst(k -> shared[k] && shared[mod1(k + 1, 4)], 1:4)
        return 2, isnothing(k) ? tiers.self : tiers.edge[k]
    elseif ns == 1
        return 3, vertex_square(T, findfirst(shared), tiers, opts)
    elseif ns > 2
        return 2, tiers.self
    end
    return 5, square_near_rule(T, S, tiers.near, opts.near_ratio, touch ? opts.near_depth + 4 : opts.near_depth)
end

"""Tier and rule on test triangle T (corner nodes `tn`) for a body source S, as `pair_rule`"""
function tri_pair_rule(T::FlatTriangle, S, tn::SVector{3,Int}, sn, self::Bool, rules::TierRules, opts::AssemblyOptions, touch::Bool)
    self && return 1, rules.self
    norm(S.centroid - T.centroid) > opts.far_ratio * max(T.diam, S.diam) && return 6, rules.far
    shared = (tn[1] in sn, tn[2] in sn, tn[3] in sn)
    contact = contact_rule(T, S, shared, rules, opts)
    isnothing(contact) || return 4, contact
    ns = count(shared)
    if ns == 2
        return 2, rules.edge[shared[1] && shared[2] ? 1 : shared[2] && shared[3] ? 2 : 3]
    elseif ns == 1
        return 3, vertex_rule(T, findfirst(shared), rules, opts)
    end
    return 5, near_rule(T, S, rules.near, opts.near_ratio, touch ? opts.near_depth + 4 : opts.near_depth)
end

"""Tier and rule on test quad T for a wake triangle S (no shared nodes: the trailing edge contact is geometric)"""
function quad_wake_rule(T::IsoQuad, S::FlatTriangle, tiers::SquareTiers, opts::AssemblyOptions)
    norm(S.centroid - T.centroid) > opts.far_ratio * max(T.diam, S.diam) && return 6, tiers.far
    contact = square_contact_rule(T, S, (false, false, false, false), tiers, opts)
    isnothing(contact) || return 4, contact
    return 5, square_near_rule(T, S, tiers.near, opts.near_ratio, opts.near_depth)
end

"""Tier and rule on test quad T for a filament segment f: the vertex rule where f starts at a corner of T"""
function quad_segment_rule(T::IsoQuad, f::Segment, tiers::SquareTiers, opts::AssemblyOptions)
    L = norm(f.b - f.a)
    norm((f.a + f.b) / 2 - T.centroid) > opts.far_ratio * max(T.diam, L) && return 6, tiers.far
    for k in 1:4
        min(norm(T.x[k] - f.a), norm(T.x[k] - f.b)) <= 1e-12 * T.diam && return 3, vertex_square(T, k, tiers, opts)
    end
    return 5, square_near_rule(T, f, tiers.near, opts.near_ratio, opts.near_depth)
end

"""
Contribution of source quad S to the rows of test quad T with the rule: the 4×4 block ∫ψ_a n·v_b dA and the source
part of the right hand side -∫ψ_a n·u_σ dA (`quad_field`, which takes σ = -V∞·n on S); `side` +1 for the self pair,
whose test points lie on S at known parameters, so their principal values need no projection (`polar_field`)
"""
function block(T::IsoQuad{F}, S::IsoQuad{F}, rule::SquareRule, V∞::SVector{3,F}, test::Symbol, side::Integer) where {F}
    blk, rb = zero(SMatrix{4,4,F,16}), zero(SVector{4,F})
    for (p, wq) in zip(rule.p, rule.w)
        nJ = area_normal(T, p)
        J = norm(nJ)
        c, s = side != 0 ? Elements.polar_field(S, point(T, p), nJ / J, V∞, p, zero(F), side) : quad_field(S, point(T, p), nJ / J, V∞)
        ψ = (wq * J) * square_test_weights(test, p)
        blk += ψ * c'
        rb -= s * ψ
    end
    return blk, rb
end

"""Contribution of source triangle S (source strength σS) to the rows of test quad T: 4×3 block and right hand side"""
function block(T::IsoQuad{F}, S::FlatTriangle{F}, rule::SquareRule, σS::Real, test::Symbol, kernel_far::Real) where {F}
    blk, rb = zero(SMatrix{4,3,F,12}), zero(SVector{4,F})
    for (p, wq) in zip(rule.p, rule.w)
        nJ = area_normal(T, p)
        J = norm(nJ)
        nx = nJ / J
        w = source_velocity(S.src, point(T, p); far=kernel_far)
        v = cross(w, nx)
        ψ = (wq * J) * square_test_weights(test, p)
        blk += ψ * SVector(dot(v, S.γ[1]), dot(v, S.γ[2]), dot(v, S.γ[3]))'
        rb -= (σS * dot(w, nx)) * ψ
    end
    return blk, rb
end

"""Contribution of source quad S to the rows of test triangle T: 3×4 block and right hand side"""
function block(T::FlatTriangle{F}, S::IsoQuad{F}, rule::BaryRule, V∞::SVector{3,F}, test::Symbol) where {F}
    blk, rb = zero(SMatrix{3,4,F,12}), zero(SVector{3,F})
    for (λ, wq) in zip(rule.λ, rule.w)
        c, s = quad_field(S, point(T, λ), T.n, V∞)
        ψ = (wq * T.area) * test_weights(test, λ)
        blk += ψ * c'
        rb -= s * ψ
    end
    return blk, rb
end

"""Contribution of filament segment f to the rows of test quad T: the parts per unit strength at f's two ends"""
function segment_block(T::IsoQuad{F}, f::Segment{F}, rule::SquareRule, test::Symbol) where {F}
    ca, cb = zero(SVector{4,F}), zero(SVector{4,F})
    for (p, wq) in zip(rule.p, rule.w)
        nJ = area_normal(T, p)
        J = norm(nJ)
        va, vb = segment_velocity(f.a, f.b, point(T, p))
        ψ = (wq * J) * square_test_weights(test, p)
        ca += dot(va, nJ / J) * ψ
        cb += dot(vb, nJ / J) * ψ
    end
    return ca, cb
end

"""
    assemble(g, dm, mesh::MixedMesh, V∞; opts=AssemblyOptions(), external=nothing, wake=nothing)

The Galerkin system of `assemble` over a mixed mesh of isoparametric quads and flat triangles: the same rows (slots)
and columns (slots, then wake nodes), the quads' inner integrals by `quad_field` and their test integrals by the
parameter-square tiers (`SquareTiers`), the triangles' as for split meshes. The freestream part of the right hand side
is integrated exactly (3×3 Gauss on a quad), the external velocity at the same points. Returns A, b and statistics.
"""
function assemble(g::SurfaceGeometry, dm::DOFMapping, mesh::MixedMesh{F}, V∞::SVector{3}; opts::AssemblyOptions=AssemblyOptions(),
    external=nothing, wake::Union{Nothing,RigidWake}=nothing) where {F}
    check_options(opts)
    nq, nt, ns = length(mesh.quads), length(mesh.tris), dm.nslot
    V = SVector{3,F}(V∞)
    rules, tiers = tier_rules(opts), square_tiers(opts)
    g3 = square_gauss(3)
    ext_rule = triangle_rule(7)
    σt = [-dot(V, T.n) for T in mesh.tris]
    touch = Set{Tuple{Int,Int}}()
    for (a, b) in g.touch_pairs
        ea, eb = mesh.element_of[a], mesh.element_of[b]
        ea == eb || push!(touch, (ea, eb), (eb, ea))
    end
    nw = isnothing(wake) ? 0 : nwake(wake)
    wsrc = isnothing(wake) ? nothing : wake_sources(wake, F)
    At = zeros(F, ns + nw, ns)
    b = zeros(F, ns)
    pairs, points = zeros(Int, length(TIER_NAMES)), zeros(Int, length(TIER_NAMES))
    lk = ReentrantLock()
    ne = nq + nt
    nchunk = opts.threaded ? min(ne, 8 * Threads.nthreads()) : 1
    chunks = [collect(c) for c in Iterators.partition(1:ne, cld(ne, nchunk))]
    time = @elapsed begin
        tasks = map(chunks) do chunk
            task = () -> begin
                buf = zeros(F, ns + nw, 4)
                my_pairs, my_points = zeros(Int, length(TIER_NAMES)), zeros(Int, length(TIER_NAMES))
                for e in chunk
                    fill!(buf, 0)
                    if e <= nq
                        rows = test_quad!(buf, my_pairs, my_points, mesh, e, V, σt, tiers, opts, touch, ns, wsrc, external, g3)
                        slots = mesh.quad_slots[e]
                    else
                        rows = test_tri!(buf, my_pairs, my_points, mesh, e - nq, V, σt, rules, opts, touch, ns, wsrc, external, ext_rule)
                        slots = mesh.tri_slots[e-nq]
                    end
                    lock(lk) do
                        for a in eachindex(slots)
                            @views At[:, slots[a]] .+= buf[:, a]
                            b[slots[a]] += rows[a]
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

"""The rows of test quad e into buf (columns 1-4), returning its right hand side"""
function test_quad!(buf, pairs, points, mesh::MixedMesh{F}, e, V, σt, tiers, opts, touch, ns, wsrc, external, g3) where {F}
    T, tn = mesh.quads[e], mesh.quad_nodes[e]
    nq = length(mesh.quads)
    rows = zero(SVector{4,F})
    for (s, S) in enumerate(mesh.quads)
        tier, rule = quad_pair_rule(T, S, tn, mesh.quad_nodes[s], s == e, tiers, opts, (e, s) in touch)
        blk, rb = block(T, S, rule, V, opts.test, s == e ? 1 : 0)
        cs = mesh.quad_slots[s]
        for bb in 1:4, a in 1:4
            buf[cs[bb], a] += blk[a, bb]
        end
        rows += rb
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    for (s, S) in enumerate(mesh.tris)
        tier, rule = quad_pair_rule(T, S, tn, mesh.tri_nodes[s], false, tiers, opts, (e, nq + s) in touch)
        blk, rb = block(T, S, rule, σt[s], opts.test, opts.kernel_far)
        cs = mesh.tri_slots[s]
        for bb in 1:3, a in 1:4
            buf[cs[bb], a] += blk[a, bb]
        end
        rows += rb
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    if !isnothing(wsrc)
        for (s, S) in enumerate(wsrc.elements)
            tier, rule = quad_wake_rule(T, S, tiers, opts)
            blk, _ = block(T, S, rule, zero(F), opts.test, opts.kernel_far)
            c = wsrc.tris[s]
            for bb in 1:3, a in 1:4
                buf[ns+c[bb], a] += blk[a, bb]
            end
            pairs[tier] += 1
            points[tier] += length(rule)
        end
        for (f, seg) in enumerate(wsrc.segments)
            tier, rule = quad_segment_rule(T, seg, tiers, opts)
            ca, cb = segment_block(T, seg, rule, opts.test)
            ia, ib = wsrc.filaments[f]
            for a in 1:4
                buf[ns+ia, a] += ca[a]
                buf[ns+ib, a] += cb[a]
            end
            pairs[tier] += 1
            points[tier] += length(rule)
        end
    end
    for (p, wq) in zip(g3.p, g3.w) # freestream (exact: a polynomial) and external velocity
        nJ = area_normal(T, p)
        u = isnothing(external) ? V : V + external(point(T, p))
        rows -= (wq * dot(u, nJ)) * square_test_weights(opts.test, p)
    end
    return rows
end

"""The rows of test triangle t into buf (columns 1-3), returning its right hand side"""
function test_tri!(buf, pairs, points, mesh::MixedMesh{F}, t, V, σt, rules, opts, touch, ns, wsrc, external, ext_rule) where {F}
    T, tn = mesh.tris[t], mesh.tri_nodes[t]
    nq = length(mesh.quads)
    rows = zero(SVector{3,F})
    for (s, S) in enumerate(mesh.quads)
        tier, rule = tri_pair_rule(T, S, tn, mesh.quad_nodes[s], false, rules, opts, (nq + t, s) in touch)
        blk, rb = block(T, S, rule, V, opts.test)
        cs = mesh.quad_slots[s]
        for bb in 1:4, a in 1:3
            buf[cs[bb], a] += blk[a, bb]
        end
        rows += rb
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    for (s, S) in enumerate(mesh.tris)
        tier, rule = tri_pair_rule(T, S, tn, mesh.tri_nodes[s], s == t, rules, opts, (nq + t, nq + s) in touch)
        blk, rb = pair_block(T, S, rule, σt[s], opts.test, s == t ? 1 : 0, opts.kernel_far)
        cs = mesh.tri_slots[s]
        for bb in 1:3, a in 1:3
            buf[cs[bb], a] += blk[a, bb]
        end
        rows += rb
        pairs[tier] += 1
        points[tier] += length(rule)
    end
    if !isnothing(wsrc)
        sub = view(buf, :, 1:3)
        add_wake!(sub, pairs, points, T, ns, wsrc, rules, opts)
    end
    rows -= dot(V, T.n) * T.area * test_integrals(opts.test)
    if !isnothing(external)
        for (λ, wq) in zip(ext_rule.λ, ext_rule.w)
            rows -= (wq * T.area * dot(external(point(T, λ)), T.n)) * test_weights(opts.test, λ)
        end
    end
    return rows
end
