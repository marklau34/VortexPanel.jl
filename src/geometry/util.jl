"""Unit vector along v, or zero for a zero vector"""
function safe_normalize(v::SVector{3,T}) where {T}
    n = norm(v)
    return n > 0 ? v / n : zero(v)
end

"""Area vector (x₂-x₁)×(x₃-x₁)/2 of a triangle"""
tri_area_vector(x::AbstractVector, t::SVector{3,Int}) = cross(x[t[2]] - x[t[1]], x[t[3]] - x[t[1]]) / 2

"""Outward unit normal of a triangle (zero when degenerate)"""
tri_normal(x::AbstractVector, t::SVector{3,Int}) = safe_normalize(tri_area_vector(x, t))

"""Distance from point p to the segment from a to b, and the parameter of the closest point"""
function point_segment(p::SVector{3,T}, a::SVector{3,T}, b::SVector{3,T}) where {T}
    ab = b - a
    L2 = dot(ab, ab)
    t = L2 > 0 ? clamp(dot(p - a, ab) / L2, zero(T), one(T)) : zero(T)
    return norm(p - (a + t * ab)), t
end

"""Parameter of the projection of p on the line through a and b (not clamped)"""
function line_param(p::SVector{3,T}, a::SVector{3,T}, b::SVector{3,T}) where {T}
    ab = b - a
    return dot(p - a, ab) / dot(ab, ab)
end

"""
Uniform hash grid over points, for neighbour queries within boxes. Cell size `h`.
"""
struct HashGrid{T}
    h::T
    cells::Dict{NTuple{3,Int},Vector{Int}}
end

cell_key(h, p) = (floor(Int, p[1] / h), floor(Int, p[2] / h), floor(Int, p[3] / h))

function HashGrid(points::AbstractVector{SVector{3,T}}, h::Real; ids=eachindex(points)) where {T}
    h = T(h)
    h > 0 || (h = one(T))
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    for i in ids
        push!(get!(cells, cell_key(h, points[i]), Int[]), i)
    end
    return HashGrid(h, cells)
end

"""Ids stored in the cells overlapping the box [lo, hi] (a superset of the points inside it)"""
function query(grid::HashGrid, lo::SVector{3}, hi::SVector{3})
    klo, khi = cell_key(grid.h, lo), cell_key(grid.h, hi)
    out = Int[]
    ncell = prod(khi .- klo .+ 1)
    if ncell > 4 * length(grid.cells) # box larger than the occupied cells: scan them instead
        for (k, ids) in grid.cells
            all(klo .<= k .<= khi) && append!(out, ids)
        end
        return out
    end
    for i in klo[1]:khi[1], j in klo[2]:khi[2], k in klo[3]:khi[3]
        ids = get(grid.cells, (i, j, k), nothing)
        isnothing(ids) || append!(out, ids)
    end
    return out
end

"""Simple union-find over 1:n with the smallest id as representative"""
struct UnionFind
    parent::Vector{Int}
end

UnionFind(n::Integer) = UnionFind(collect(1:n))

function find_root(uf::UnionFind, i::Integer)
    p = uf.parent
    while p[i] != i
        p[i] = p[p[i]]
        i = p[i]
    end
    return i
end

function unite!(uf::UnionFind, i::Integer, j::Integer)
    ri, rj = find_root(uf, i), find_root(uf, j)
    ri == rj && return ri
    r, s = minmax(ri, rj)
    uf.parent[s] = r
    return r
end

"""Connected components (lists of ids) of a union-find, ordered by their smallest id"""
function components(uf::UnionFind)
    groups = Dict{Int,Vector{Int}}()
    for i in eachindex(uf.parent)
        push!(get!(groups, find_root(uf, i), Int[]), i)
    end
    return [groups[r] for r in sort!(collect(keys(groups)))]
end
