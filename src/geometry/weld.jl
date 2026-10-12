"""Welded nodes and the faces that survive welding, in welded node ids, with their input element ids"""
struct WeldResult{T}
    x::Vector{SVector{3,T}}
    input_node::Vector{Int}
    faces::Vector{Vector{Int}}
    face_parent::Vector{Int}
end

"""
    weld(raw, rel, log)

Merge nodes closer than `rel` times their shortest incident edge (the smaller tolerance of the pair). Candidate pairs
come from a hash grid and are merged closest first with union-find; a merge that would join two nodes of one face
farther apart than their tolerance is refused, so chains of near nodes cannot collapse an element. Each cluster takes
the lowest input id as representative and the mean position. Welded ids follow the input order of the representatives.

Faces are then renumbered: repeated nodes collapse (a quad with one repeated node becomes a triangle), faces with fewer
than three distinct nodes or repeating a non-adjacent node are dropped, and duplicate faces (same node set) are dropped.
"""
function weld(raw::RawSurface{T}, rel::Real, log::Vector{RepairEvent{T}}) where {T}
    x = raw.x
    n = length(x)

    # shortest incident edge of every node, ignoring degenerate edges (shorter than `rel` times its longest incident
    # edge, which welding removes), and the faces of every node
    incident = [T[] for _ in 1:n]
    node_faces = [Int[] for _ in 1:n]
    for (f, face) in enumerate(raw.faces)
        m = length(face)
        for k in 1:m
            a, b = face[k], face[mod1(k + 1, m)]
            push!(node_faces[a], f)
            a == b && continue
            L = norm(x[a] - x[b])
            push!(incident[a], L)
            push!(incident[b], L)
        end
    end
    ℓ = map(incident) do Ls
        isempty(Ls) && return T(Inf)
        Lmax = maximum(Ls)
        return minimum((L for L in Ls if L > T(rel) * Lmax); init=T(Inf))
    end
    finite = sort!(filter(isfinite, ℓ))
    ℓ_default = isempty(finite) ? one(T) : finite[cld(length(finite), 2)]
    tol = [T(rel) * (isfinite(l) ? l : ℓ_default) for l in ℓ]

    # candidate pairs within tolerance
    grid = HashGrid(x, max(maximum(tol; init=zero(T)), eps(T)))
    pairs = Tuple{T,Int,Int}[]
    for i in 1:n
        r = SVector{3,T}(tol[i], tol[i], tol[i])
        for j in query(grid, x[i] - r, x[i] + r)
            j <= i && continue
            d = norm(x[i] - x[j])
            d <= min(tol[i], tol[j]) && push!(pairs, (d, i, j))
        end
    end
    sort!(pairs)

    uf = UnionFind(n)
    members = [[i] for i in 1:n]
    for (d, i, j) in pairs
        ri, rj = find_root(uf, i), find_root(uf, j)
        ri == rj && continue
        refused = false
        for p in members[ri], q in members[rj]
            if !isdisjoint(node_faces[p], node_faces[q]) && norm(x[p] - x[q]) > min(tol[p], tol[q])
                refused = true
                break
            end
        end
        if refused
            push!(log, RepairEvent(:weld_refused, x[i], d, "nodes $i and $j are within tolerance but would merge nodes of one element"))
            continue
        end
        r = unite!(uf, ri, rj)
        s = r == ri ? rj : ri
        append!(members[r], members[s])
        empty!(members[s])
    end

    # welded nodes: representatives in input order, mean positions
    reps = [i for i in 1:n if find_root(uf, i) == i]
    new_id = zeros(Int, n)
    for (k, r) in enumerate(reps)
        new_id[r] = k
    end
    input_node = [new_id[find_root(uf, i)] for i in 1:n]
    xw = Vector{SVector{3,T}}(undef, length(reps))
    for (k, r) in enumerate(reps)
        mem = members[r]
        xw[k] = sum(x[mem]) / length(mem)
        if length(mem) > 1
            spread = maximum(norm(x[p] - x[q]) for p in mem, q in mem)
            push!(log, RepairEvent(:weld, xw[k], spread, "merged input nodes $(sort(mem))"))
            spread > minimum(tol[mem]) &&
                push!(log, RepairEvent(:weld_wide, xw[k], spread, "merged cluster wider than its tolerance $(minimum(tol[mem]))"))
        end
    end

    # renumber faces, collapse repeated nodes, drop degenerate and duplicate faces
    faces = Vector{Int}[]
    face_parent = Int[]
    seen = Dict{Vector{Int},Int}()
    for (f, face) in enumerate(raw.faces)
        c = Int[]
        for i in face
            v = input_node[i]
            (isempty(c) || c[end] != v) && push!(c, v)
        end
        length(c) > 1 && c[1] == c[end] && pop!(c)
        centre = sum(xw[c]) / max(length(c), 1)
        if length(c) < 3 || !allunique(c)
            push!(log, RepairEvent(:drop_degenerate, centre, zero(T), "input element $f collapses to nodes $c"))
            continue
        end
        length(c) < length(face) && push!(log, RepairEvent(:collapse, centre, zero(T), "input element $f: $(length(face)) → $(length(c)) nodes"))
        key = sort(c)
        if haskey(seen, key)
            push!(log, RepairEvent(:drop_duplicate, centre, zero(T), "input element $f duplicates input element $(seen[key])"))
            continue
        end
        seen[key] = f
        push!(faces, c)
        push!(face_parent, f)
    end
    return WeldResult(xw, input_node, faces, face_parent)
end
