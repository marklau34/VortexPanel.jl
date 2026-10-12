"""Join lines of wake trailing edge nodes where one ends at the node the next starts from; mark lines that close"""
function join_lines(lines::Vector{Vector{Int}}, closed::Vector{Bool})
    lines, closed = copy(lines), copy(closed)
    merged = true
    while merged
        merged = false
        for i in eachindex(lines), j in eachindex(lines)
            (i == j || closed[i] || closed[j]) && continue
            if lines[i][end] == lines[j][1]
                lines[i] = [lines[i]; lines[j][2:end]]
                deleteat!(lines, j)
                deleteat!(closed, j)
                merged = true
                break
            end
        end
    end
    for i in eachindex(lines)
        if !closed[i] && length(lines[i]) > 2 && lines[i][end] == lines[i][1]
            pop!(lines[i])
            closed[i] = true
        end
    end
    return lines, closed
end

"""
    wake_trailing_edge(g, corner_slot, nslot, log)

The wake's trailing edge nodes and the matrices U and L (nwake_te × nslot) that give their upper and lower doublet
values from the slots, μ_u = U μ and μ_l = L μ, so that the nodal Kutta condition μ_w = μ_u - μ_l reads K = (U - L) C.

- Matching trailing edges (Kutta chains): every chain node is a wake trailing edge node, with μ_u and μ_l its corner
  slots in the upper and lower triangles of an adjacent segment.
- Non-matching trailing edges (Kutta pair chains): the wake trailing edge nodes are the union of both chains' nodes;
  μ_u is interpolated linearly along the upper chain and μ_l along the lower one, so the trailing edge filament cancels
  exactly along every segment (a piecewise linear function stays piecewise linear on any refinement).

This is the single place where the Kutta condition is built (a pressure-type Kutta condition would replace it).
"""
function wake_trailing_edge(g::SurfaceGeometry{T}, corner_slot, nslot::Integer, log::Vector{String}) where {T}
    cslot(t, v) = corner_slot[t][findfirst(==(v), g.tris[t])]
    normal(t) = tri_normal(g.x, g.tris[t])
    nodes = Int[]
    index = Dict{Int,Int}()
    up_terms = Vector{Vector{Tuple{Int,T}}}()
    low_terms = Vector{Vector{Tuple{Int,T}}}()
    wedge = SVector{3,T}[]
    function wake_node!(v, ut, lt, b)
        k = get(index, v, 0)
        if k == 0
            push!(nodes, v)
            push!(up_terms, ut)
            push!(low_terms, lt)
            push!(wedge, b)
            k = index[v] = length(nodes)
        else
            wedge[k] += b # the node joins two trailing edge lines
        end
        return k
    end

    lines = Vector{Vector{Int}}()
    closed = Bool[]
    for ch in g.kutta.chains
        vs = ch.closed && ch.nodes[end] == ch.nodes[1] ? ch.nodes[1:end-1] : ch.nodes
        nseg = length(ch.edges)
        line = Int[]
        for (i, v) in enumerate(vs)
            segs = ch.closed ? unique([mod1(i - 1, nseg), mod1(i, nseg)]) : [k for k in (i - 1, i) if 1 <= k <= nseg]
            k = first(segs)
            b = sum(normal(ch.tri_upper[s]) + normal(ch.tri_lower[s]) for s in segs)
            push!(line, wake_node!(v, [(cslot(ch.tri_upper[k], v), one(T))], [(cslot(ch.tri_lower[k], v), one(T))], b))
        end
        push!(lines, line)
        push!(closed, ch.closed)
    end
    for pch in g.kutta.pair_chains
        line = Int[]
        for (j, v) in enumerate(pch.union_nodes)
            us, ut = pch.union_upper[j]
            ls, lt = pch.union_lower[j]
            tu, tl = pch.tri_upper[us], pch.tri_lower[ls]
            uterms = [(cslot(tu, pch.upper[us]), 1 - ut), (cslot(tu, pch.upper[us+1]), ut)]
            lterms = [(cslot(tl, pch.lower[ls]), 1 - lt), (cslot(tl, pch.lower[ls+1]), lt)]
            push!(line, wake_node!(v, uterms, lterms, normal(tu) + normal(tl)))
        end
        push!(lines, line)
        push!(closed, false)
    end
    lines, closed = join_lines(lines, closed)

    tangent = [zero(SVector{3,T}) for _ in nodes]
    for (line, c) in zip(lines, closed)
        m = length(line)
        for (i, k) in enumerate(line)
            prev = i > 1 ? line[i-1] : (c ? line[end] : k)
            next = i < m ? line[i+1] : (c ? line[1] : k)
            tangent[k] = Geometry.safe_normalize(g.x[nodes[next]] - g.x[nodes[prev]])
        end
    end
    bisector = [Geometry.safe_normalize(b) for b in wedge]
    any(iszero, bisector) && push!(log, "a wake trailing edge node has no defined trailing edge bisector")

    selection(terms) = sparse([k for (k, ts) in enumerate(terms) for _ in ts], [s for ts in terms for (s, _) in ts],
        [w for ts in terms for (_, w) in ts], length(terms), nslot)
    U, L = selection(up_terms), selection(low_terms)
    return WakeTE(nodes, g.x[nodes], tangent, bisector, lines, closed), U, L
end
