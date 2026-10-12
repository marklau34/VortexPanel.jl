"""Triangles from splitting the welded faces, with their parent faces and the parent quads"""
struct SplitResult
    tris::Vector{SVector{3,Int}}
    tri_parent::Vector{Int}
    tri_slot::Vector{Int8}
    quads::Vector{SVector{4,Int}}
    quad_parent::Vector{Int}
    quad_diag::Vector{Int8}
    quad_reason::Vector{Symbol}
end

"""Angle between the normals of triangles (a, b, c) and (a, c, d), the fold of quad a-b-c-d along diagonal a-c"""
function fold_angle(a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}, d::SVector{3,T}) where {T}
    n1 = cross(b - a, c - a)
    n2 = cross(c - a, d - a)
    m = norm(n1) * norm(n2)
    m > 0 || return zero(T)
    return acos(clamp(dot(n1, n2) / m, -one(T), one(T)))
end

"""
    choose_diagonal(x, q, parent, opts)

Diagonal of quad `q` = (a, b, c, d) to split along: 1 for a-c, 2 for b-d. See `GeometryOptions` for the rules. The tie
break compares, for each diagonal v, v[m]·v[p]·(c[m] - m₀), with m the mirror axis, m₀ the mirror plane, c the quad
centroid and p the other axis along which the quad extends most. This depends only on the geometry of the quad, not on
its node order or orientation, and is the same for a quad and its mirror image.
"""
function choose_diagonal(x::AbstractVector{SVector{3,T}}, q::AbstractVector{Int}, parent::Integer, opts::GeometryOptions) where {T}
    opts.split == :d13 && return (Int8(1), :fixed)
    opts.split == :d24 && return (Int8(2), :fixed)
    if opts.split == :pattern
        d = opts.split_pattern[parent]
        d in (1, 2) || throw(ArgumentError("split_pattern[$parent] = $d is not 1 or 2"))
        return (Int8(d), :fixed)
    end
    a, b, c, d = x[q[1]], x[q[2]], x[q[3]], x[q[4]]
    L13, L24 = norm(c - a), norm(d - b)
    tie = T(opts.split_tie_rel)
    if opts.split == :min_fold || opts.split == :auto
        f13 = fold_angle(a, b, c, d)
        f24 = fold_angle(b, c, d, a)
        if opts.split == :min_fold || max(f13, f24) > deg2rad(T(opts.split_fold_deg))
            abs(f13 - f24) > tie && return ((f13 < f24 ? Int8(1) : Int8(2)), :fold)
        end
    end
    abs(L13 - L24) > tie * max(L13, L24) && return ((L13 < L24 ? Int8(1) : Int8(2)), :shortest)

    # tied: geometric, mirror symmetric rule
    if !isnothing(opts.split_mirror)
        m, m0 = opts.split_mirror
        centre = (a + b + c + d) / 4
        size = max(L13, L24)
        offset = centre[m] - T(m0)
        if abs(offset) > tie * size
            lo = min.(a, b, c, d)
            hi = max.(a, b, c, d)
            others = [k for k in 1:3 if k != m]
            p = hi[others[1]] - lo[others[1]] >= hi[others[2]] - lo[others[2]] ? others[1] : others[2]
            v13, v24 = c - a, d - b
            k13 = v13[m] * v13[p] * offset
            k24 = v24[m] * v24[p] * offset
            abs(k13 - k24) > tie * size^2 * abs(offset) && return ((k13 > k24 ? Int8(1) : Int8(2)), :tie_mirror)
        end
    end
    return ((min(q[1], q[3]) < min(q[2], q[4]) ? Int8(1) : Int8(2)), :tie_id)
end

"""
    split_faces(x, weld, opts)

Triangulate the welded faces: triangles are kept, quads are split by `choose_diagonal` (both halves keep the quad's
orientation), and larger polygons are fanned from their first node. Records the parent face and the half of the quad.
"""
function split_faces(x::AbstractVector{SVector{3,T}}, w::WeldResult{T}, opts::GeometryOptions, log::Vector{RepairEvent{T}}) where {T}
    tris = SVector{3,Int}[]
    tri_parent = Int[]
    tri_slot = Int8[]
    quads = SVector{4,Int}[]
    quad_parent = Int[]
    quad_diag = Int8[]
    quad_reason = Symbol[]
    for (f, face) in enumerate(w.faces)
        parent = w.face_parent[f]
        n = length(face)
        if n == 3
            push!(tris, SVector{3,Int}(face[1], face[2], face[3]))
            push!(tri_parent, parent)
            push!(tri_slot, Int8(0))
        elseif n == 4
            a, b, c, d = face
            diag, reason = choose_diagonal(x, face, parent, opts)
            halves = diag == 1 ? (SVector(a, b, c), SVector(a, c, d)) : (SVector(a, b, d), SVector(b, c, d))
            for (s, t) in enumerate(halves)
                push!(tris, t)
                push!(tri_parent, parent)
                push!(tri_slot, Int8(s))
            end
            push!(quads, SVector(a, b, c, d))
            push!(quad_parent, parent)
            push!(quad_diag, diag)
            push!(quad_reason, reason)
        else
            push!(log, RepairEvent(:fan, sum(x[face]) / n, zero(T), "input element $parent: $n-gon fanned into $(n - 2) triangles"))
            for k in 2:n-1
                push!(tris, SVector(face[1], face[k], face[k+1]))
                push!(tri_parent, parent)
                push!(tri_slot, Int8(0))
            end
        end
    end
    return SplitResult(tris, tri_parent, tri_slot, quads, quad_parent, quad_diag, quad_reason)
end
