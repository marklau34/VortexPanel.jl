"""A boundary loop: nodes in the direction of their boundary half-edges (the body on their left)"""
struct BoundaryLoop
    nodes::Vector{Int}
    closed::Bool # false: an open chain that could not be closed
end

"""
    find_openings(x, tris, he, taken, opts, log)

Trace the boundary half-edges whose edges are not in `taken` (T-junction, gap and non-matching trailing edge edges) into
loops. From a half-edge into node q the next one leaves q in the same fan (found by rotating around q). Chains that end
where another starts, within the gap tolerance, are joined, so an opening bounded partly by a gap still forms one loop.
"""
function find_openings(x::AbstractVector{SVector{3,T}}, tris, he::HalfEdges, taken::Set{Tuple{Int,Int}}, opts::GeometryOptions,
    log::Vector{RepairEvent{T}}) where {T}
    nh = length(he.twin)
    is_open = falses(nh)
    for h in 1:nh
        he.twin[h] == 0 && !(Tuple(he.edge_nodes[he.he_edge[h]]) in taken) && (is_open[h] = true)
    end
    open_he = findall(is_open)
    loops = BoundaryLoop[]
    isempty(open_he) && return loops

    leaving = Dict{Int,Vector{Int}}()
    for h in open_he
        push!(get!(leaving, he_origin(tris, h), Int[]), h)
    end
    next = zeros(Int, nh)
    for h in open_he
        q = he_dest(tris, h)
        cands = get(leaving, q, Int[])
        isempty(cands) && continue
        if length(cands) == 1
            next[h] = cands[1]
            continue
        end
        g = he_next(h) # rotate around q within h's fan
        steps = 0
        while he.twin[g] != 0 && steps < 1000
            g = he_next(he.twin[g])
            steps += 1
        end
        next[h] = is_open[g] ? g : cands[1]
    end
    prev = zeros(Int, nh)
    for h in open_he
        next[h] != 0 && (prev[next[h]] = h)
    end

    visited = falses(nh)
    chains = Vector{Int}[] # open chains, as half-edge lists
    for h0 in open_he
        (visited[h0] || prev[h0] != 0) && continue
        hs = Int[]
        h = h0
        while h != 0 && !visited[h]
            visited[h] = true
            push!(hs, h)
            h = next[h]
        end
        push!(chains, hs)
    end
    for h0 in open_he
        visited[h0] && continue
        nodes = Int[]
        h = h0
        while !visited[h]
            visited[h] = true
            push!(nodes, he_origin(tris, h))
            h = next[h]
        end
        push!(loops, BoundaryLoop(nodes, true))
    end

    # join open chains end to start within the gap tolerance
    isempty(chains) && return loops
    first_node(c) = he_origin(tris, c[1])
    last_node(c) = he_dest(tris, c[end])
    function local_len(c)
        return (norm(x[he_dest(tris, c[1])] - x[he_origin(tris, c[1])]) + norm(x[he_dest(tris, c[end])] - x[he_origin(tris, c[end])])) / 2
    end
    succ = zeros(Int, length(chains))
    taken_start = falses(length(chains))
    for i in eachindex(chains)
        best, bestd = 0, T(Inf)
        for j in eachindex(chains)
            taken_start[j] && continue
            d = norm(x[first_node(chains[j])] - x[last_node(chains[i])])
            d <= T(opts.gap_rel) * min(local_len(chains[i]), local_len(chains[j])) && d < bestd && ((best, bestd) = (j, d))
        end
        if best != 0
            succ[i] = best
            taken_start[best] = true
        end
    end
    done = falses(length(chains))
    for i in eachindex(chains)
        done[i] && continue
        # walk back to the start of this sequence
        start = i
        seen = Set{Int}([i])
        while true
            p = findfirst(==(start), succ)
            (isnothing(p) || p in seen) && break
            push!(seen, p)
            start = p
        end
        nodes = Int[]
        j = start
        lastj = start
        closed = false
        while j != 0 && !done[j]
            done[j] = true
            append!(nodes, [he_origin(tris, h) for h in chains[j]])
            lastj = j
            j = succ[j]
            j == start && (closed = true)
        end
        if !closed
            push!(nodes, last_node(chains[lastj]))
            push!(log, RepairEvent(:opening_unclosed, x[nodes[1]], zero(T), "boundary chain of $(length(nodes)) nodes from node $(nodes[1]) does not close into a loop; left open"))
        end
        push!(loops, BoundaryLoop(nodes, closed))
    end
    return loops
end
