# Small meshes and comparisons shared by the tests

"""Triangles in a canonical form (rotated to start at their smallest node, sorted), to compare triangulations"""
canonical(tris) = sort([(i = argmin(t); SVector(t[i], t[mod1(i + 1, 3)], t[mod1(i + 2, 3)])) for t in tris])

"""Unit cube of 8 nodes and 6 outward quads, with the quad on face x = 1 optionally split into two quads at z = 1/2"""
function quad_cube(; split_face::Bool=false)
    X = [SVector{3,Float64}(i, j, k) for k in 0:1 for j in 0:1 for i in 0:1]
    id(i, j, k) = 1 + i + 2j + 4k
    faces = [[id(0, 0, 0), id(0, 1, 0), id(1, 1, 0), id(1, 0, 0)], # z = 0
        [id(0, 0, 1), id(1, 0, 1), id(1, 1, 1), id(0, 1, 1)],      # z = 1
        [id(0, 0, 0), id(1, 0, 0), id(1, 0, 1), id(0, 0, 1)],      # y = 0
        [id(0, 1, 0), id(0, 1, 1), id(1, 1, 1), id(1, 1, 0)],      # y = 1
        [id(0, 0, 0), id(0, 0, 1), id(0, 1, 1), id(0, 1, 0)]]      # x = 0
    if split_face
        push!(X, SVector(1.0, 0.0, 0.5), SVector(1.0, 1.0, 0.5))
        m0, m1 = length(X) - 1, length(X)
        push!(faces, [id(1, 0, 0), id(1, 1, 0), m1, m0], [m0, m1, id(1, 1, 1), id(1, 0, 1)])
    else
        push!(faces, [id(1, 0, 0), id(1, 1, 0), id(1, 1, 1), id(1, 0, 1)])
    end
    return X, faces
end

"""
Unit cube whose faces are strips of quads along x: the top (z = 1) and back (y = 1) faces in halves, the front
(y = 0) and bottom (z = 0) faces in thirds, the ends single quads. The top-front and bottom-back seams have
non-matching nodes on coincident lines (zero-width gaps); every face is built on its own nodes and welded.
"""
function seam_cube()
    X = SVector{3,Float64}[]
    faces = Vector{Int}[]
    function strip!(n, p) # quads along x of the face p(x, s), s ∈ [0, 1] across it
        for i in 0:n-1
            x0, x1 = i / n, (i + 1) / n
            base = length(X)
            append!(X, [p(x0, 0.0), p(x1, 0.0), p(x1, 1.0), p(x0, 1.0)])
            push!(faces, collect(base+1:base+4))
        end
    end
    strip!(2, (x, s) -> SVector(x, s, 1.0))       # top, halves
    strip!(3, (x, s) -> SVector(x, 0.0, s))       # front, thirds
    strip!(3, (x, s) -> SVector(x, s, 0.0))       # bottom, thirds
    strip!(2, (x, s) -> SVector(x, 1.0, s))       # back, halves
    for x0 in (0.0, 1.0)                          # ends
        base = length(X)
        append!(X, [SVector(x0, 0.0, 0.0), SVector(x0, 1.0, 0.0), SVector(x0, 1.0, 1.0), SVector(x0, 0.0, 1.0)])
        push!(faces, collect(base+1:base+4))
    end
    return X, faces
end

"""Prism of rhombus section (sharp leading and trailing edges) spanning y = ±1, as raw arrays, tips open"""
function diamond_wing(; n::Integer=4, t::Real=0.1)
    X = SVector{3,Float64}[]
    for j in 0:n
        y = -1 + 2j / n
        append!(X, [SVector(0.0, y, 0.0), SVector(0.5, y, t / 2), SVector(1.0, y, 0.0), SVector(0.5, y, -t / 2)])
    end
    id(j, k) = 4j + k
    faces = Vector{Int}[]
    for j in 0:n-1, k in 1:4
        k2 = mod1(k + 1, 4)
        push!(faces, [id(j, k), id(j + 1, k), id(j + 1, k2), id(j, k2)])
    end
    return X, faces
end
