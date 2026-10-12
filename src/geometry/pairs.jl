# Detection of T-junctions and of nearly coincident boundary chains (small gaps and non-matching trailing edges).
# Both work on node positions only, so they run before the orientation is fixed and are recorded once it is.

"""Hanging nodes `hanging` (sorted along a→b, at parameters `t`) on the boundary edge from `a` to `b`"""
struct TJCandidate{T}
    a::Int
    b::Int
    hanging::Vector{Int}
    t::Vector{T}
end

"""Two boundary chains paired edge by edge: node chains, their edges, the paired edges, and the largest distance"""
struct PairCandidate{T}
    chain_a::Vector{Int}
    chain_b::Vector{Int}
    edges_a::Vector{Int}
    edges_b::Vector{Int}
    partners::Vector{Tuple{Int,Int}}
    width::T
end

boundary_edges(he::HalfEdges) = [e for e in eachindex(he.edge_he) if he.edge_he[e][2] == 0]

"""
    find_tjunctions(x, he, opts, log)

Boundary edges (a, b) with boundary nodes within `tj_rel`·|ab| of their interior, such that boundary edges chain a,
h₁, ..., b on the other side.
"""
function find_tjunctions(x::AbstractVector{SVector{3,T}}, he::HalfEdges, opts::GeometryOptions, log::Vector{RepairEvent{T}}) where {T}
    cands = TJCandidate{T}[]
    bedges = boundary_edges(he)
    isempty(bedges) && return cands
    bset = Set{Tuple{Int,Int}}(Tuple(he.edge_nodes[e]) for e in bedges)
    bnodes = unique(reduce(vcat, [collect(he.edge_nodes[e]) for e in bedges]))
    Ls = sort([norm(x[he.edge_nodes[e][2]] - x[he.edge_nodes[e][1]]) for e in bedges])
    grid = HashGrid(x, Ls[cld(length(Ls), 2)]; ids=bnodes)
    tmin = T(opts.tj_tmin)
    for e in bedges
        a, b = he.edge_nodes[e]
        L = norm(x[b] - x[a])
        tol = T(opts.tj_rel) * L
        r = SVector{3,T}(tol, tol, tol)
        hs = Tuple{T,Int}[]
        for v in query(grid, min.(x[a], x[b]) - r, max.(x[a], x[b]) + r)
            (v == a || v == b) && continue
            d, t = point_segment(x[v], x[a], x[b])
            d <= tol && tmin < t < 1 - tmin && push!(hs, (t, v))
        end
        isempty(hs) && continue
        sort!(hs)
        chain = [a; last.(hs); b]
        if all(minmax(chain[i], chain[i+1]) in bset for i in 1:length(chain)-1)
            push!(cands, TJCandidate(a, b, last.(hs), first.(hs)))
        else
            push!(log, RepairEvent(:tjunction_partial, x[hs[1][2]], hs[1][1], "node(s) $(last.(hs)) lie on edge ($a, $b) but do not chain along it"))
        end
    end
    return cands
end

"""Edges (sorted node pairs) taken by T-junctions: each long edge and the short edges along it"""
function tjunction_edges(cands::AbstractVector{<:TJCandidate})
    s = Set{Tuple{Int,Int}}()
    for c in cands
        push!(s, minmax(c.a, c.b))
        chain = [c.a; c.hanging; c.b]
        for i in 1:length(chain)-1
            push!(s, minmax(chain[i], chain[i+1]))
        end
    end
    return s
end

"""
Overlap of segment [r, s] on segment [p, q]: parameters (length units along p→q) of the overlapping part
"""
function overlap_on(xp, xq, xr, xs)
    L = norm(xq - xp)
    u = (xq - xp) / L
    tr, ts = dot(xr - xp, u), dot(xs - xp, u)
    lo, hi = max(zero(L), min(tr, ts)), min(L, max(tr, ts))
    return lo, hi, u
end

"""
Whether boundary edges (p, q) and (r, s) are paired, and the Hausdorff distance between them over their overlap
"""
function edge_pair(x, p, q, r, s, opts::GeometryOptions)
    xp, xq, xr, xs = x[p], x[q], x[r], x[s]
    T = eltype(xp)
    Le, Lf = norm(xq - xp), norm(xs - xr)
    (Le > 0 && Lf > 0) || return false, zero(T)
    abs(dot(xq - xp, xs - xr)) / (Le * Lf) >= cos(deg2rad(T(opts.gap_angle_deg))) || return false, zero(T)
    lo, hi, u = overlap_on(xp, xq, xr, xs)
    hi - lo >= T(opts.gap_overlap) * min(Le, Lf) || return false, zero(T)
    lo2, hi2, u2 = overlap_on(xr, xs, xp, xq)
    hi2 - lo2 >= T(opts.gap_overlap) * min(Le, Lf) || return false, zero(T)
    d = max(point_segment(xp + lo * u, xr, xs)[1], point_segment(xp + hi * u, xr, xs)[1],
        point_segment(xr + lo2 * u2, xp, xq)[1], point_segment(xr + hi2 * u2, xp, xq)[1])
    return d <= T(opts.gap_rel) * min(Le, Lf), d
end

"""Nodes of a set of edges ordered into a path (or a closed loop), or nothing when the edges branch"""
function order_path(edge_nodes::AbstractVector{SVector{2,Int}}, edges::AbstractVector{Int})
    adj = Dict{Int,Vector{Int}}()
    for e in edges, v in edge_nodes[e]
        push!(get!(adj, v, Int[]), e)
    end
    any(length(es) > 2 for es in values(adj)) && return nothing
    ends = sort([v for (v, es) in adj if length(es) == 1])
    start = isempty(ends) ? minimum(keys(adj)) : ends[1]
    nodes = [start]
    used = Set{Int}()
    v = start
    while true
        next_e = 0
        for e in adj[v]
            e in used || (next_e = e; break)
        end
        next_e == 0 && break
        push!(used, next_e)
        a, b = edge_nodes[next_e]
        v = a == v ? b : a
        push!(nodes, v)
    end
    length(used) == length(edges) || return nothing
    return nodes
end

"""Location (segment, parameter) of point p on the polyline through `chain`, and its distance"""
function project_on_chain(x, p::SVector{3,T}, chain::AbstractVector{Int}) where {T}
    best = (T(Inf), 1, zero(T))
    for k in 1:length(chain)-1
        d, t = point_segment(p, x[chain[k]], x[chain[k+1]])
        d < best[1] && (best = (d, k, t))
    end
    return (best[2], best[3]), best[1]
end

"""
    find_chain_pairs(x, he, taken, opts, log)

Boundary edges (not in `taken`) paired with a nearly parallel, overlapping boundary edge within the Hausdorff tolerance,
grouped into two chains: partners lie on opposite chains, and paired edges that share a node without being partners lie
on the same chain.
"""
function find_chain_pairs(x::AbstractVector{SVector{3,T}}, he::HalfEdges, taken::Set{Tuple{Int,Int}}, opts::GeometryOptions,
    log::Vector{RepairEvent{T}}) where {T}
    bedges = [e for e in boundary_edges(he) if !(Tuple(he.edge_nodes[e]) in taken)]
    result = PairCandidate{T}[]
    isempty(bedges) && return result
    L = [norm(x[he.edge_nodes[e][2]] - x[he.edge_nodes[e][1]]) for e in bedges]
    mids = [(x[he.edge_nodes[e][1]] + x[he.edge_nodes[e][2]]) / 2 for e in bedges]
    Lmax = maximum(L)
    grid = HashGrid(mids, Lmax)

    partner = Dict{Int,Vector{Int}}()
    width = Dict{Tuple{Int,Int},T}()
    for (i, e) in enumerate(bedges)
        p, q = he.edge_nodes[e]
        rad = L[i] / 2 + Lmax / 2 + T(opts.gap_rel) * L[i]
        r = SVector{3,T}(rad, rad, rad)
        for j in query(grid, mids[i] - r, mids[i] + r)
            j <= i && continue
            f = bedges[j]
            rr, ss = he.edge_nodes[f]
            ok, d = edge_pair(x, p, q, rr, ss, opts)
            ok || continue
            push!(get!(partner, e, Int[]), f)
            push!(get!(partner, f, Int[]), e)
            width[minmax(e, f)] = d
        end
    end
    isempty(partner) && return result

    # two-colour the paired edges: partners opposite, non-partner neighbours along the boundary alike
    paired = sort!(collect(keys(partner)))
    node_edges = Dict{Int,Vector{Int}}()
    for e in paired, v in he.edge_nodes[e]
        push!(get!(node_edges, v, Int[]), e)
    end
    colour = Dict{Int,Int}()
    group = Dict{Int,Int}()
    ngroup = 0
    conflict = Set{Int}()
    for seed in paired
        haskey(colour, seed) && continue
        ngroup += 1
        colour[seed] = 0
        group[seed] = ngroup
        stack = [seed]
        while !isempty(stack)
            e = pop!(stack)
            nbrs = Tuple{Int,Int}[(f, 1) for f in partner[e]]
            for v in he.edge_nodes[e], f in node_edges[v]
                f == e || f in partner[e] || push!(nbrs, (f, 0))
            end
            for (f, rel) in nbrs
                want = xor(colour[e], rel)
                if !haskey(colour, f)
                    colour[f] = want
                    group[f] = ngroup
                    push!(stack, f)
                elseif colour[f] != want
                    push!(conflict, ngroup)
                end
            end
        end
    end

    for g in 1:ngroup
        es = [e for e in paired if group[e] == g]
        side_a = [e for e in es if colour[e] == 0]
        side_b = [e for e in es if colour[e] == 1]
        chain_a = order_path(he.edge_nodes, side_a)
        chain_b = order_path(he.edge_nodes, side_b)
        if g in conflict || isnothing(chain_a) || isnothing(chain_b)
            e = es[1]
            push!(log, RepairEvent(:pair_irregular, x[he.edge_nodes[e][1]], zero(T),
                "boundary edges near node $(he.edge_nodes[e][1]) look like a gap but do not form two simple chains; left open"))
            continue
        end
        partners = Tuple{Int,Int}[(e, f) for e in side_a for f in partner[e]]
        w = maximum(width[minmax(e, f)] for (e, f) in partners)
        push!(result, PairCandidate(chain_a, chain_b, side_a, side_b, partners, w))
    end
    return result
end

"""Edges (sorted node pairs) of all chain pairs"""
function pair_edges(he::HalfEdges, cands::AbstractVector{<:PairCandidate})
    s = Set{Tuple{Int,Int}}()
    for c in cands, e in [c.edges_a; c.edges_b]
        push!(s, Tuple(he.edge_nodes[e]))
    end
    return s
end

"""
Orientation links between triangles across T-junctions and chain pairs: (t1, t2, need) where `need` says that exactly
one of the two must be flipped for their boundary half-edges to run in opposite directions.
"""
function orientation_links(x, tris, he::HalfEdges, tjc::AbstractVector{<:TJCandidate}, pc::AbstractVector{<:PairCandidate})
    links = Tuple{Int,Int,Bool}[]
    hdir(h) = x[he_dest(tris, h)] - x[he_origin(tris, h)]
    for c in tjc
        hl = he.edge_he[edge_id(he, c.a, c.b)][1]
        chain = [c.a; c.hanging; c.b]
        for i in 1:length(chain)-1
            hs = he.edge_he[edge_id(he, chain[i], chain[i+1])][1]
            push!(links, (he_tri(hl), he_tri(hs), dot(hdir(hl), hdir(hs)) > 0))
        end
    end
    for c in pc, (e, f) in c.partners
        he_, hf = he.edge_he[e][1], he.edge_he[f][1]
        push!(links, (he_tri(he_), he_tri(hf), dot(hdir(he_), hdir(hf)) > 0))
    end
    return links
end
