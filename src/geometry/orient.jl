"""Triangle flips that orient every body consistently and outward, the body of every triangle, and body volumes"""
struct Orientation{T}
    flip::BitVector
    body::Vector{Int}
    volume::Vector{T}
end

"""Signed volume of a set of triangles closed by a fan of its boundary to `centre` (the fan terms vanish relative to it)"""
function signed_volume(x::AbstractVector{SVector{3,T}}, tris, ts, centre::SVector{3,T}) where {T}
    V = zero(T)
    for t in ts
        a, b, c = tris[t]
        V += dot(x[a] - centre, cross(x[b] - centre, x[c] - centre))
    end
    return V / 6
end

"""
    orient(x, tris, he, links, log)

1. Within each patch (triangles connected through shared edges) flip triangles so that twins run in opposite
   directions; a contradiction means a non-orientable surface and is an error.
2. Flip whole patches linked through T-junctions and chain pairs (`links`) so that the linked boundary half-edges run in
   opposite directions; contradicting links are logged and ignored. Linked patches form a body.
3. Flip each body whose signed volume, closed by a fan of its boundary to its node centroid, is negative.
"""
function orient(x::AbstractVector{SVector{3,T}}, tris::AbstractVector{SVector{3,Int}}, he::HalfEdges,
    links::AbstractVector{Tuple{Int,Int,Bool}}, log::Vector{RepairEvent{T}}) where {T}
    nt = length(tris)
    flip = falses(nt)
    patch = zeros(Int, nt)
    np = 0
    for seed in 1:nt
        patch[seed] != 0 && continue
        np += 1
        patch[seed] = np
        stack = [seed]
        while !isempty(stack)
            t = pop!(stack)
            for k in 1:3
                h = he_index(t, k)
                g = he.twin[h]
                g == 0 && continue
                u = he_tri(g)
                want = xor(flip[t], he_origin(tris, h) == he_origin(tris, g))
                if patch[u] == 0
                    patch[u] = np
                    flip[u] = want
                    push!(stack, u)
                elseif flip[u] != want
                    a, b = he_origin(tris, h), he_dest(tris, h)
                    throw(ArgumentError("Non-orientable surface: cannot orient triangles $t and $u consistently across edge ($a, $b) at $(x[a])"))
                end
            end
        end
    end

    # link patches
    padj = [Tuple{Int,Bool}[] for _ in 1:np]
    for (t1, t2, need) in links
        p1, p2 = patch[t1], patch[t2]
        rel = xor(need, flip[t1], flip[t2])
        if p1 == p2
            rel && push!(log, RepairEvent(:link_conflict, x[tris[t1][1]], zero(T), "triangles $t1 and $t2 are linked across a gap or T-junction with inconsistent orientation"))
            continue
        end
        push!(padj[p1], (p2, rel))
        push!(padj[p2], (p1, rel))
    end
    pflip = falses(np)
    pbody = zeros(Int, np)
    nb = 0
    for seed in 1:np
        pbody[seed] != 0 && continue
        nb += 1
        pbody[seed] = nb
        stack = [seed]
        while !isempty(stack)
            p = pop!(stack)
            for (q, rel) in padj[p]
                want = xor(pflip[p], rel)
                if pbody[q] == 0
                    pbody[q] = nb
                    pflip[q] = want
                    push!(stack, q)
                elseif pflip[q] != want
                    push!(log, RepairEvent(:link_conflict, zero(SVector{3,T}), zero(T), "patches $p and $q are linked with inconsistent orientations"))
                end
            end
        end
    end
    flip .= xor.(flip, pflip[patch])
    body = pbody[patch]

    # outward
    oriented = [flip[t] ? SVector(tr[1], tr[3], tr[2]) : tr for (t, tr) in enumerate(tris)]
    volume = zeros(T, nb)
    for b in 1:nb
        ts = findall(==(b), body)
        nodes = unique(reduce(vcat, [collect(oriented[t]) for t in ts]))
        centre = sum(x[nodes]) / length(nodes)
        V = signed_volume(x, oriented, ts, centre)
        if V < 0
            flip[ts] .= .!flip[ts]
            V = -V
        end
        lo = reduce((p, q) -> min.(p, q), x[nodes])
        hi = reduce((p, q) -> max.(p, q), x[nodes])
        V <= sqrt(eps(T)) * maximum(hi - lo)^3 && push!(log, RepairEvent(:flat_body, centre, V, "body $b encloses (almost) no volume; its outward direction is arbitrary"))
        volume[b] = V
    end
    return Orientation(flip, body, volume)
end
