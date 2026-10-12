"""Welded node pairs of the input edge tags named in `names`"""
function welded_tag_edges(raw::RawSurface, input_node::AbstractVector{Int}, names)
    s = Set{Tuple{Int,Int}}()
    for name in names
        haskey(raw.edge_tags, name) || continue
        for (a, b) in raw.edge_tags[name]
            wa, wb = input_node[a], input_node[b]
            wa != wb && push!(s, minmax(wa, wb))
        end
    end
    return s
end

"""Node of triangle t that is not on edge (a, b)"""
opposite_node(tri::SVector{3,Int}, a::Integer, b::Integer) = tri[findfirst(v -> v != a && v != b, tri)]

"""
Whether two outward faces with normals n1, n2 meeting along an edge form a sharp, convex wedge (included angle at most
`te_angle`) pointing along `wake_dir`. `opp2` is a point of the second face off the edge and `xa` a point on it.
"""
function is_te_wedge(n1, n2, xa, opp2, cte, wake)
    dot(n1, n2) <= -cte || return false
    dot(opp2 - xa, n1) < 0 || return false # convex: the second face lies behind the first
    return dot(n1 + n2, wake) > 0
end

"""
    classify_kutta_edges(x, tris, he, tagged, input_node, opts, log)

Kutta edges (sorted node pairs) among the edges shared by two triangles: the dihedral test, the input tags, or both
(`opts.kutta_source`), then the include and exclude lists. Returns the set and the number of edges where the dihedral
test and the tags disagree (when there are tags).
"""
function classify_kutta_edges(x::AbstractVector{SVector{3,T}}, tris, he::HalfEdges, tagged::Set{Tuple{Int,Int}},
    input_node::AbstractVector{Int}, opts::GeometryOptions, log::Vector{RepairEvent{T}}) where {T}
    cte = cos(deg2rad(T(opts.te_angle_deg)))
    wake = SVector{3,T}(opts.wake_dir)
    dihedral = Set{Tuple{Int,Int}}()
    for e in eachindex(he.edge_nodes)
        h1, h2 = he.edge_he[e]
        h2 == 0 && continue
        t1, t2 = he_tri(h1), he_tri(h2)
        a, b = he.edge_nodes[e]
        n1, n2 = tri_normal(x, tris[t1]), tri_normal(x, tris[t2])
        is_te_wedge(n1, n2, x[a], x[opposite_node(tris[t2], a, b)], cte, wake) && push!(dihedral, (a, b))
    end
    manifold_tags = Set(p for p in tagged if (e = edge_id(he, p...); e > 0 && he.edge_he[e][2] != 0))
    mismatch = 0
    if !isempty(tagged)
        for p in symdiff(dihedral, manifold_tags)
            mismatch += 1
            push!(log, RepairEvent(:kutta_tag_mismatch, (x[p[1]] + x[p[2]]) / 2, zero(T),
                "edge $p is $(p in dihedral ? "" : "not ")a trailing edge by the dihedral test but $(p in manifold_tags ? "is" : "is not") tagged"))
        end
    end
    kutta = opts.kutta_source == :dihedral ? dihedral : opts.kutta_source == :tags ? manifold_tags : union(dihedral, manifold_tags)
    for (a, b) in opts.kutta_exclude
        delete!(kutta, minmax(input_node[a], input_node[b]))
    end
    for (a, b) in opts.kutta_include
        p = minmax(input_node[a], input_node[b])
        e = edge_id(he, p...)
        if e > 0 && he.edge_he[e][2] != 0
            push!(kutta, p)
        else
            push!(log, RepairEvent(:kutta_include_ignored, x[p[1]], zero(T), "kutta_include edge ($a, $b) is not shared by two triangles"))
        end
    end
    return kutta, mismatch
end

"""
    classify_pair(x, tris, he, c, opts)

`:kutta` when the two sides of a chain pair fold back on each other like a trailing edge wedge (by a majority of the
paired edges), `:gap` otherwise.
"""
function classify_pair(x::AbstractVector{SVector{3,T}}, tris, he::HalfEdges, c::PairCandidate, opts::GeometryOptions) where {T}
    cte = cos(deg2rad(T(opts.te_angle_deg)))
    wake = SVector{3,T}(opts.wake_dir)
    votes = 0
    for (e, f) in c.partners
        te, tf = he_tri(he.edge_he[e][1]), he_tri(he.edge_he[f][1])
        a = he.edge_nodes[e][1]
        r, s = he.edge_nodes[f]
        votes += is_te_wedge(tri_normal(x, tris[te]), tri_normal(x, tris[tf]), x[a], x[opposite_node(tris[tf], r, s)], cte, wake)
    end
    return 2votes > length(c.partners) ? :kutta : :gap
end

"""Chain direction rule: reverse `nodes` unless (wake_dir × τ)·up_ref ≥ 0 for the chord τ of the chain"""
function orient_along_tau(x::AbstractVector{SVector{3,T}}, nodes::Vector{Int}, opts::GeometryOptions) where {T}
    length(nodes) < 2 && return nodes
    τ = x[nodes[end]] - x[nodes[1]]
    norm(τ) > 0 || (τ = x[nodes[2]] - x[nodes[1]]) # closed chain: first segment
    s = dot(cross(SVector{3,T}(opts.wake_dir), τ), SVector{3,T}(opts.up_ref))
    return s < -sqrt(eps(T)) * norm(τ) ? reverse(nodes) : nodes
end

arclength(x, nodes) = [zero(eltype(eltype(x))); cumsum([norm(x[nodes[k+1]] - x[nodes[k]]) for k in 1:length(nodes)-1])]

"""
    kutta_chains(x, tris, he, kutta, opts)

Order the Kutta edges into chains between nodes where they end or branch (closed chains when they loop), directed
along τ. The upper triangle of each segment owns the half-edge that runs along τ. End types are filled in later.
"""
function kutta_chains(x::AbstractVector{SVector{3,T}}, tris, he::HalfEdges, kutta::Set{Tuple{Int,Int}}, opts::GeometryOptions) where {T}
    adj = Dict{Int,Vector{Tuple{Int,Int}}}()
    for p in kutta, v in p
        push!(get!(adj, v, Tuple{Int,Int}[]), p)
    end
    used = Set{Tuple{Int,Int}}()
    raw_chains = Tuple{Vector{Int},Bool}[]
    function walk(v0, p0)
        nodes = [v0]
        v, p = v0, p0
        while true
            push!(used, p)
            w = p[1] == v ? p[2] : p[1]
            push!(nodes, w)
            w == v0 && return nodes, true
            length(adj[w]) == 2 || return nodes, false
            q = adj[w][1] == p ? adj[w][2] : adj[w][1]
            q in used && return nodes, false
            v, p = w, q
        end
    end
    for v in sort!(collect(keys(adj)))
        length(adj[v]) == 2 && continue
        for p in adj[v]
            p in used || push!(raw_chains, walk(v, p))
        end
    end
    for p in sort!(collect(kutta))
        p in used || push!(raw_chains, walk(p[1], p))
    end

    chains = KuttaChain{T}[]
    for (nodes, closed) in raw_chains
        nodes = orient_along_tau(x, nodes, opts)
        edges, tu, tl = Int[], Int[], Int[]
        for k in 1:length(nodes)-1
            a, b = nodes[k], nodes[k+1]
            e = edge_id(he, a, b)
            h1, h2 = he.edge_he[e]
            up, low = he_origin(tris, h1) == a ? (h1, h2) : (h2, h1)
            push!(edges, e)
            push!(tu, he_tri(up))
            push!(tl, he_tri(low))
        end
        push!(chains, KuttaChain(nodes, arclength(x, nodes), edges, tu, tl, closed, (:dangling, :dangling), (0, 0)))
    end
    return chains
end

"""
    kutta_pair_chain(x, tris, he, pair, k, opts)

Upper and lower chains of a non-matching trailing edge (chain pair `k`), directed along τ: the upper side owns boundary
half-edges along τ. The union of their nodes, ordered by arclength along the upper chain, gives the wake's trailing
edge nodes; nodes on both chains appear once.
"""
function kutta_pair_chain(x::AbstractVector{SVector{3,T}}, tris, he::HalfEdges, pair::ChainPair{T}, k::Integer, opts::GeometryOptions) where {T}
    a = orient_along_tau(x, pair.chain_a, opts)
    b = pair.chain_b
    if norm(x[b[1]] - x[a[1]]) + norm(x[b[end]] - x[a[end]]) > norm(x[b[1]] - x[a[end]]) + norm(x[b[end]] - x[a[1]])
        b = reverse(b)
    end
    # boundary half-edge of a chain segment and whether it runs along the chain
    function seg_he(c, i)
        e = edge_id(he, c[i], c[i+1])
        h = he.edge_he[e][1]
        return h, he_origin(tris, h) == c[i]
    end
    along_a = count(i -> seg_he(a, i)[2], 1:length(a)-1)
    upper, lower = 2along_a >= length(a) - 1 ? (a, b) : (b, a)
    tri_upper = [he_tri(seg_he(upper, i)[1]) for i in 1:length(upper)-1]
    tri_lower = [he_tri(seg_he(lower, i)[1]) for i in 1:length(lower)-1]

    s_up = arclength(x, upper)
    entries = Tuple{T,Int,Int8}[] # (s, node, side)
    for (i, v) in enumerate(upper)
        push!(entries, (s_up[i], v, Int8(v in lower ? 0 : 1)))
    end
    for v in lower
        v in upper && continue
        (seg, t), _ = project_on_chain(x, x[v], upper)
        push!(entries, (s_up[seg] + t * (s_up[seg+1] - s_up[seg]), v, Int8(-1)))
    end
    sort!(entries)
    loc(chain, v) = (i = findfirst(==(v), chain); isnothing(i) ? project_on_chain(x, x[v], chain)[1] : i < length(chain) ? (i, zero(T)) : (i - 1, one(T)))
    union_nodes = [e[2] for e in entries]
    return KuttaPairChain(k, upper, lower, tri_upper, tri_lower, union_nodes, [e[3] for e in entries],
        [loc(upper, v) for v in union_nodes], [loc(lower, v) for v in union_nodes], [e[1] for e in entries])
end
