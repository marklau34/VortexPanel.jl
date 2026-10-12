# Seeded, reproducible mesh pathologies on MeshLibrary meshes (validation stage 1). Element ids, boundaries and element
# sets of the input are kept; new nodes and elements are appended (new elements join the sets of their parent).

"""Record of an abuse operation: human readable events and counts of the pathologies created"""
struct AbuseLog
    events::Vector{String}
    counts::Dict{Symbol,Int}
end

AbuseLog() = AbuseLog(String[], Dict{Symbol,Int}())

function record!(log::AbuseLog, key::Symbol, msg::String; n::Integer=1)
    push!(log.events, msg)
    log.counts[key] = get(log.counts, key, 0) + n
    return log
end

Base.merge(a::AbuseLog, b::AbuseLog) = AbuseLog([a.events; b.events], mergewith(+, a.counts, b.counts))

coords(mesh) = [SVector{3}(p.coord) for p in mesh.nodes]
elem_nodes(mesh) = [collect(Int, MeshLibrary.get_nodes(e)) for e in mesh.elements]

"""Mesh with new node coordinates and faces; same element type when possible, else Poly2D{Int,4}"""
function rebuild(mesh::MeshLibrary.Mesh, X::AbstractVector{<:SVector{3}}, faces::AbstractVector{<:AbstractVector{Int}};
    new_parents::AbstractVector{Int}=Int[])
    nodes = MeshLibrary.StructVector([MeshLibrary.Point(p) for p in X])
    n = unique(length.(faces))
    elems = if n == [4]
        [MeshLibrary.Quad(f...) for f in faces]
    elseif n == [3]
        [MeshLibrary.Tri(f...) for f in faces]
    else
        [MeshLibrary.Poly2D{Int,4}(f) for f in faces]
    end
    nold = length(faces) - length(new_parents)
    sets = map(mesh.elemset) do s
        extra = [nold + k for (k, p) in enumerate(new_parents) if p in s.elems]
        MeshLibrary.MeshElementSet(s.name, [collect(s.elems); extra])
    end
    return MeshLibrary.Mesh(nodes, MeshLibrary.StructVector(elems); boundaries=deepcopy(mesh.boundaries),
        element_sets=isempty(sets) ? copy(mesh.elemset) : collect(sets))
end

"""Mesh with node coordinates X and the same elements"""
with_nodes(mesh::MeshLibrary.Mesh, X) = rebuild(mesh, X, elem_nodes(mesh))

"""Undirected edges with the elements that use them"""
function edge_elements(faces)
    d = Dict{Tuple{Int,Int},Vector{Int}}()
    for (e, f) in enumerate(faces), k in eachindex(f)
        push!(get!(d, minmax(f[k], f[mod1(k + 1, length(f))]), Int[]), e)
    end
    return d
end

"""Nodes on open boundary edges (used by a single element)"""
open_nodes(faces) = Set(v for (k, es) in edge_elements(faces) if length(es) == 1 for v in k)

"""Nodes on boundaries named trailing_edge*"""
function te_nodes(mesh::MeshLibrary.Mesh)
    s = Set{Int}()
    for b in mesh.boundaries
        startswith(b.name, "trailing_edge") || continue
        for f in b.facets
            con = MeshLibrary.get_nodes(mesh.elements[f.element])
            push!(s, con[f.side], con[mod1(f.side + 1, length(con))])
        end
    end
    return s
end

"""Nodes that abuse leaves alone by default: trailing edge nodes and open boundary (tip) nodes"""
default_fixed(mesh::MeshLibrary.Mesh) = union(te_nodes(mesh), open_nodes(elem_nodes(mesh)))

"""Elements abuse leaves alone by default: those touching a trailing edge or open boundary node, or on a boundary"""
function default_avoid(mesh::MeshLibrary.Mesh)
    fixed = default_fixed(mesh)
    faces = elem_nodes(mesh)
    avoid = Set(e for (e, f) in enumerate(faces) if any(in(fixed), f))
    for b in mesh.boundaries, f in b.facets
        push!(avoid, f.element)
    end
    return avoid
end

"""Unit normal at every node (area-weighted face normals) and the shortest incident edge"""
function node_frames(X, faces)
    n = [zero(SVector{3,eltype(eltype(X))}) for _ in X]
    ℓ = fill(eltype(eltype(X))(Inf), length(X))
    for f in faces
        m = length(f)
        a = cross(X[f[2]] - X[f[1]], X[f[3]] - X[f[1]])
        m == 4 && (a += cross(X[f[3]] - X[f[1]], X[f[4]] - X[f[1]]))
        for k in 1:m
            v, w = f[k], f[mod1(k + 1, m)]
            n[v] += a
            L = norm(X[w] - X[v])
            ℓ[v] = min(ℓ[v], L)
            ℓ[w] = min(ℓ[w], L)
        end
    end
    return [norm(v) > 0 ? normalize(v) : v for v in n], ℓ
end

"""Closest point to p on the triangles (quads split 1-3) of `faces` among `cand`"""
function project_on_faces(p, X, faces, cand)
    best, bestd = p, Inf
    for e in cand
        f = faces[e]
        for tri in (length(f) == 4 ? ((f[1], f[2], f[3]), (f[1], f[3], f[4])) : ((f[1], f[2], f[3]),))
            q = closest_on_triangle(p, X[tri[1]], X[tri[2]], X[tri[3]])
            d = norm(q - p)
            d < bestd && ((best, bestd) = (q, d))
        end
    end
    return best
end

"""Closest point to p on triangle abc (Ericson, Real-Time Collision Detection, 5.1.5)"""
function closest_on_triangle(p, a, b, c)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    (d1 <= 0 && d2 <= 0) && return a
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    (d3 >= 0 && d4 <= d3) && return b
    vc = d1 * d4 - d3 * d2
    (vc <= 0 && d1 >= 0 && d3 <= 0) && return a + d1 / (d1 - d3) * ab
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    (d6 >= 0 && d5 <= d6) && return c
    vb = d5 * d2 - d1 * d6
    (vb <= 0 && d2 >= 0 && d6 <= 0) && return a + d2 / (d2 - d6) * ac
    va = d3 * d6 - d5 * d4
    (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) && return b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)
    denom = 1 / (va + vb + vc)
    return a + ab * (vb * denom) + ac * (vc * denom)
end

"""Elements of every node"""
function node_elements(nnode, faces)
    ne = [Int[] for _ in 1:nnode]
    for (e, f) in enumerate(faces), v in f
        push!(ne[v], e)
    end
    return ne
end

"""
    jitter_nodes(mesh, δ; rng, project=nothing, fixed=default_fixed(mesh))

Move every free node tangentially in a random direction by up to δ times its shortest incident edge, then back onto the
surface: onto `project(p)` when given (e.g. `naca_projector` for an exact wing), else onto the original faceted
surface around the node.
"""
function jitter_nodes(mesh::MeshLibrary.Mesh, δ::Real; rng::AbstractRNG=StableRNG(1), project=nothing, fixed=default_fixed(mesh))
    X = coords(mesh)
    faces = elem_nodes(mesh)
    normals, ℓ = node_frames(X, faces)
    nel = node_elements(length(X), faces)
    Y = copy(X)
    log = AbuseLog()
    for v in eachindex(X)
        (v in fixed || isempty(nel[v])) && continue
        n = normals[v]
        t1 = normalize(abs(n[1]) < 0.9 ? cross(n, SVector(1.0, 0, 0)) : cross(n, SVector(0, 1.0, 0)))
        t2 = cross(n, t1)
        θ = 2π * rand(rng)
        p = X[v] + δ * ℓ[v] * rand(rng) * (cos(θ) * t1 + sin(θ) * t2)
        Y[v] = isnothing(project) ? project_on_faces(p, X, faces, nel[v]) : project(p)
    end
    record!(log, :jittered_nodes, "jittered $(length(X) - length(fixed)) nodes by up to $δ of their shortest edge"; n=count(v -> !(v in fixed), eachindex(X)))
    return with_nodes(mesh, Y), log
end

"""Pick up to n elements from `cand` in random order, no two sharing a node"""
function pick_separated(rng, faces, cand, n)
    picked = Int[]
    used = Set{Int}()
    for e in shuffle(rng, sort!(collect(cand)))
        length(picked) == n && break
        any(in(used), faces[e]) && continue
        push!(picked, e)
        union!(used, faces[e])
    end
    return picked
end

"""
    insert_tjunctions(mesh, n; rng, avoid=default_avoid(mesh))

Split n elements (no two sharing a node) without splitting their neighbours: a quad into two quads across the middle
of a random pair of opposite sides, a triangle into two triangles from the middle of a random side. The new nodes are
the exact midpoints of the straight sides, so each hangs on the neighbour's long edge.
"""
function insert_tjunctions(mesh::MeshLibrary.Mesh, n::Integer; rng::AbstractRNG=StableRNG(1), avoid=default_avoid(mesh))
    X = coords(mesh)
    faces = elem_nodes(mesh)
    ee = edge_elements(faces)
    interior(a, b) = length(get(ee, minmax(a, b), Int[])) == 2
    cand = [e for e in eachindex(faces) if !(e in avoid) && all(k -> interior(faces[e][k], faces[e][mod1(k + 1, length(faces[e]))]), eachindex(faces[e]))]
    log = AbuseLog()
    parents = Int[]
    for e in pick_separated(rng, faces, cand, n)
        f = faces[e]
        if length(f) == 4
            s = rand(rng, 1:2) # split sides s and s+2
            a, b, c, d = f[s], f[s+1], f[mod1(s + 2, 4)], f[mod1(s + 3, 4)]
            push!(X, (X[a] + X[b]) / 2, (X[c] + X[d]) / 2)
            mab, mcd = length(X) - 1, length(X)
            faces[e] = [a, mab, mcd, d]
            push!(faces, [mab, b, c, mcd])
            record!(log, :tjunctions, "element $e split across sides $s and $(s + 2): hanging nodes $mab, $mcd"; n=2)
        else
            s = rand(rng, 1:3)
            a, b, c = f[s], f[mod1(s + 1, 3)], f[mod1(s + 2, 3)]
            push!(X, (X[a] + X[b]) / 2)
            m = length(X)
            faces[e] = [a, m, c]
            push!(faces, [m, b, c])
            record!(log, :tjunctions, "element $e split from side $s: hanging node $m")
        end
        push!(parents, e)
    end
    length(parents) < n && push!(log.events, "only $(length(parents)) of $n elements could be split")
    return rebuild(mesh, X, faces; new_parents=parents), log
end

"""The straightest path of k edges from node v0 through nodes not in `blocked` (nothing if there is none)"""
function straight_path(X, faces, nel, v0, k, blocked)
    nbrs(v) = unique(w for e in nel[v] for (i, u) in enumerate(faces[e]) if u == v for w in (faces[e][mod1(i + 1, length(faces[e]))], faces[e][mod1(i - 1, length(faces[e]))]))
    path = [v0]
    for _ in 1:k
        v = path[end]
        cands = [w for w in nbrs(v) if !(w in path) && !(w in blocked)]
        isempty(cands) && return nothing
        if length(path) == 1
            push!(path, cands[1])
        else
            dir = normalize(X[v] - X[path[end-1]])
            push!(path, cands[argmax([dot(normalize(X[w] - X[v]), dir) for w in cands])])
        end
    end
    return path
end

"""
    left_fan(faces, nel, p_prev, p, p_next)

Elements around node p on the left of the path p_prev → p → p_next, for consistently oriented elements: from the
element holding the directed edge p_prev → p, across the elements' edges at p, to the element holding p → p_next.
Nothing when the fan cannot be traced (a boundary or inconsistent orientation).
"""
function left_fan(faces, nel, p_prev::Integer, p::Integer, p_next::Integer)
    after(f, v) = (i = findfirst(==(v), f); f[mod1(i + 1, length(f))])
    before(f, v) = (i = findfirst(==(v), f); f[mod1(i - 1, length(f))])
    start = findfirst(e -> before(faces[e], p) == p_prev, nel[p])
    isnothing(start) && return nothing
    fan = [nel[p][start]]
    for _ in 1:64
        f = fan[end]
        w = after(faces[f], p)
        w == p_next && return fan
        g = findfirst(e -> e != f && before(faces[e], p) == w, nel[p])
        isnothing(g) && return nothing
        push!(fan, nel[p][g])
    end
    return nothing
end

"""
    open_gaps(mesh, n, ε; rng, k=3, avoid=default_avoid(mesh), project=nothing)

Open n slits: along a straight path of k edges, every interior path node is duplicated for the elements on the left of
the path (traced through the element orientation), and the copy moves tangentially by ε times the local edge length
into those elements. The path ends stay shared, so each slit is a gap between two nearly coincident chains (choose ε
between the weld and the gap tolerances). The elements must be consistently oriented.
"""
function open_gaps(mesh::MeshLibrary.Mesh, n::Integer, ε::Real; rng::AbstractRNG=StableRNG(1), k::Integer=3,
    avoid=default_avoid(mesh), project=nothing)
    X = coords(mesh)
    faces = elem_nodes(mesh)
    normals, ℓ = node_frames(X, faces)
    nel = node_elements(length(X), faces)
    blocked = Set(v for e in avoid for v in faces[e])
    union!(blocked, default_fixed(mesh))
    log = AbuseLog()
    starts = shuffle(rng, [v for v in eachindex(X) if !(v in blocked)])
    for v0 in starts
        get(log.counts, :gaps, 0) == n && break
        path = straight_path(X, faces, nel, v0, k, blocked)
        isnothing(path) && continue
        fans = [left_fan(faces, nel, path[i-1], path[i], path[i+1]) for i in 2:k] # before any element changes
        any(isnothing, fans) && continue
        for (i, fan) in zip(2:k, fans)
            p = path[i]
            d = sum(sum(X[faces[e]]) / length(faces[e]) for e in fan) / length(fan) - X[p]
            d -= dot(d, normals[p]) * normals[p]
            q = X[p] + ε * ℓ[p] * normalize(d)
            push!(X, isnothing(project) ? project_on_faces(q, X, faces, nel[p]) : project(q))
            for e in fan
                faces[e] = [u == p ? length(X) : u for u in faces[e]]
            end
        end
        # keep later slits away from this one
        for v in path, e in nel[v]
            union!(blocked, faces[e])
        end
        record!(log, :gaps, "slit along nodes $path, opened by $ε of the edge length")
    end
    get(log.counts, :gaps, 0) < n && push!(log.events, "only $(get(log.counts, :gaps, 0)) of $n gaps could be opened")
    return rebuild(mesh, X, faces), log
end

"""
    insert_slivers(mesh, n, ratio; rng, avoid=default_avoid(mesh), project=nothing)

Move n free nodes (no two in one element) along one of their edges until that edge is `ratio` of its length, so the
elements around it become slivers. The panel count does not change.
"""
function insert_slivers(mesh::MeshLibrary.Mesh, n::Integer, ratio::Real; rng::AbstractRNG=StableRNG(1), avoid=default_avoid(mesh),
    project=nothing)
    X = coords(mesh)
    faces = elem_nodes(mesh)
    nel = node_elements(length(X), faces)
    blocked = union(Set(v for e in avoid for v in faces[e]), default_fixed(mesh))
    log = AbuseLog()
    for v in shuffle(rng, [v for v in eachindex(X) if !(v in blocked)])
        get(log.counts, :slivers, 0) == n && break
        v in blocked && continue
        nbrs = unique(faces[e][mod1(i + s, length(faces[e]))] for e in nel[v] for (i, u) in enumerate(faces[e]) if u == v for s in (-1, 1))
        isempty(nbrs) && continue
        w = nbrs[rand(rng, 1:length(nbrs))]
        q = X[w] + ratio * (X[v] - X[w])
        X[v] = isnothing(project) ? project_on_faces(q, X, faces, nel[v]) : project(q)
        for e in nel[v]
            union!(blocked, faces[e])
        end
        record!(log, :slivers, "node $v moved to $ratio of its edge to node $w")
    end
    return with_nodes(mesh, X), log
end

"""
    abuse(mesh, spec; seed=1)

Apply the pathologies in `spec` (a NamedTuple; any of `jitter = δ`, `slivers = (n, ratio)`, `gaps = (n, ε)`,
`tjunctions = n`, plus `project` and `avoid`) in that order with one seeded random stream. Returns the mesh and an
`AbuseLog`.
"""
function abuse(mesh::MeshLibrary.Mesh, spec::NamedTuple; seed::Integer=1)
    rng = StableRNG(seed)
    project = get(spec, :project, nothing)
    log = AbuseLog()
    if haskey(spec, :jitter)
        mesh, l = jitter_nodes(mesh, spec.jitter; rng=rng, project=project)
        log = merge(log, l)
    end
    if haskey(spec, :slivers)
        mesh, l = insert_slivers(mesh, spec.slivers...; rng=rng, project=project, avoid=get(spec, :avoid, default_avoid(mesh)))
        log = merge(log, l)
    end
    if haskey(spec, :gaps)
        mesh, l = open_gaps(mesh, spec.gaps...; rng=rng, project=project, avoid=get(spec, :avoid, default_avoid(mesh)))
        log = merge(log, l)
    end
    if haskey(spec, :tjunctions)
        mesh, l = insert_tjunctions(mesh, spec.tjunctions; rng=rng, avoid=get(spec, :avoid, default_avoid(mesh)))
        log = merge(log, l)
    end
    return mesh, log
end
