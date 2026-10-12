# Half-edge h = 3(t-1)+k of triangle t runs from tris[t][k] to tris[t][mod1(k+1,3)]; the triangle lies on its left
# when seen from the side its normal points to.

he_index(t::Integer, k::Integer) = 3 * (t - 1) + k
he_tri(h::Integer) = (h - 1) ÷ 3 + 1
he_local(h::Integer) = (h - 1) % 3 + 1
he_next(h::Integer) = 3 * ((h - 1) ÷ 3) + mod1(he_local(h) + 1, 3)
he_prev(h::Integer) = 3 * ((h - 1) ÷ 3) + mod1(he_local(h) + 2, 3)
he_origin(tris, h::Integer) = tris[he_tri(h)][he_local(h)]
he_dest(tris, h::Integer) = tris[he_tri(h)][mod1(he_local(h) + 1, 3)]

"""Half-edge connectivity of a triangulation: twins, edge of every half-edge, and the unique edges"""
struct HalfEdges
    twin::Vector{Int}
    he_edge::Vector{Int}
    edge_nodes::Vector{SVector{2,Int}}
    edge_he::Vector{SVector{2,Int}}
    edge_of::Dict{Tuple{Int,Int},Int}
end

"""
    build_halfedges(x, tris)

Pair the half-edges of `tris` through a Dict on sorted node pairs. An edge with more than two half-edges throws a
`NonManifoldEdgeError`. Twins are paired regardless of direction (an inconsistent orientation is fixed later).
"""
function build_halfedges(x::AbstractVector{SVector{3,T}}, tris::AbstractVector{SVector{3,Int}}) where {T}
    edge_of = Dict{Tuple{Int,Int},Int}()
    lists = Vector{Vector{Int}}()
    edge_nodes = SVector{2,Int}[]
    for t in eachindex(tris), k in 1:3
        a, b = tris[t][k], tris[t][mod1(k + 1, 3)]
        key = minmax(a, b)
        e = get(edge_of, key, 0)
        if e == 0
            push!(edge_nodes, SVector(key))
            push!(lists, Int[])
            e = length(edge_nodes)
            edge_of[key] = e
        end
        push!(lists[e], he_index(t, k))
    end
    nh = 3 * length(tris)
    twin = zeros(Int, nh)
    he_edge = zeros(Int, nh)
    edge_he = Vector{SVector{2,Int}}(undef, length(lists))
    for (e, hs) in enumerate(lists)
        if length(hs) > 2
            a, b = edge_nodes[e]
            throw(NonManifoldEdgeError((a, b), (x[a], x[b]), length(hs)))
        end
        for h in hs
            he_edge[h] = e
        end
        if length(hs) == 2
            twin[hs[1]] = hs[2]
            twin[hs[2]] = hs[1]
            edge_he[e] = SVector(hs[1], hs[2])
        else
            edge_he[e] = SVector(hs[1], 0)
        end
    end
    return HalfEdges(twin, he_edge, edge_nodes, edge_he, edge_of)
end

"""Half-edges with no twin"""
boundary_halfedges(he::HalfEdges) = findall(iszero, he.twin)

"""Edge id of the node pair (a, b) in either order, 0 if there is none"""
edge_id(he::HalfEdges, a::Integer, b::Integer) = get(he.edge_of, minmax(a, b), 0)

"""
    vertex_fans(nnode, tris, twin)

Number of fans (groups of triangles connected through shared edges around the vertex) at every node; more than one
marks a bow-tie vertex.
"""
function vertex_fans(nnode::Integer, tris::AbstractVector{SVector{3,Int}}, twin::AbstractVector{Int})
    node_tris = [Int[] for _ in 1:nnode]
    for t in eachindex(tris), k in 1:3
        push!(node_tris[tris[t][k]], t)
    end
    fans = zeros(Int, nnode)
    for v in 1:nnode
        ts = node_tris[v]
        isempty(ts) && continue
        uf = UnionFind(length(ts))
        local_of = Dict(t => i for (i, t) in enumerate(ts))
        for (i, t) in enumerate(ts), k in 1:3
            h = he_index(t, k)
            (tris[t][k] == v || tris[t][mod1(k + 1, 3)] == v) || continue
            g = twin[h]
            g == 0 && continue
            j = get(local_of, he_tri(g), 0)
            j > 0 && unite!(uf, i, j)
        end
        fans[v] = length(components(uf))
    end
    return fans
end
