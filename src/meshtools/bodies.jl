# Closed test bodies, returned as MeshLibrary meshes like the meshes of the geometry packages.

"""Orient triangle or quad node list `f` so its normal points along `outward`"""
function orient_outward(X, f::AbstractVector{Int}, outward)
    n = cross(X[f[2]] - X[f[1]], X[f[3]] - X[f[1]])
    length(f) == 4 && (n += cross(X[f[3]] - X[f[1]], X[f[4]] - X[f[1]]))
    return dot(n, outward) < 0 ? reverse(f) : f
end

points(X) = [MeshLibrary.Point(p) for p in X]

"""
    icosphere(level; r=1.0, center=(0,0,0))

Watertight, outward triangulated sphere: an icosahedron subdivided `level` times with its new nodes pushed onto the
sphere (`10·4^level + 2` nodes). Port of `tm_icosphere` from CutCellMethods.jl's tests.
"""
function icosphere(level::Integer; r::Real=1.0, center=SA[0.0, 0.0, 0.0])
    φ = (1 + sqrt(5)) / 2
    P = [SVector(-1.0, φ, 0), SVector(1.0, φ, 0), SVector(-1.0, -φ, 0), SVector(1.0, -φ, 0),
        SVector(0.0, -1, φ), SVector(0.0, 1, φ), SVector(0.0, -1, -φ), SVector(0.0, 1, -φ),
        SVector(φ, 0, -1), SVector(φ, 0, 1), SVector(-φ, 0, -1), SVector(-φ, 0, 1)]
    P = [normalize(p) for p in P]
    F = [(1, 12, 6), (1, 6, 2), (1, 2, 8), (1, 8, 11), (1, 11, 12), (2, 6, 10), (6, 12, 5),
        (12, 11, 3), (11, 8, 7), (8, 2, 9), (4, 10, 5), (4, 5, 3), (4, 3, 7), (4, 7, 9),
        (4, 9, 10), (5, 10, 6), (3, 5, 12), (7, 3, 11), (9, 7, 8), (10, 9, 2)]
    for _ in 1:level
        mid = Dict{Tuple{Int,Int},Int}()
        m(a, b) = get!(mid, minmax(a, b)) do
            push!(P, normalize(P[a] + P[b]))
            length(P)
        end
        F = [f for (a, b, c) in F for f in (let ab = m(a, b), bc = m(b, c), ca = m(c, a)
            ((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca))
        end)]
    end
    c = SVector{3,Float64}(center)
    X = [r * p + c for p in P]
    tris = [MeshLibrary.Tri(orient_outward(X, collect(f), sum(X[collect(f)]) / 3 - c)...) for f in F]
    return MeshLibrary.Mesh(points(X), tris; element_sets=[MeshLibrary.MeshElementSet("sphere", collect(1:length(tris)))])
end

"""
    cube_sphere(n; r=1.0)

Sphere of quads: each face of a cube is an n×n equiangular grid projected onto the sphere. The faces are built
independently, so the nodes on the cube edges and corners are duplicated (12n + 4 extra nodes) and need welding; the
welded mesh has 6n² + 2 nodes. Quads are oriented outward.
"""
function cube_sphere(n::Integer; r::Real=1.0)
    X = SVector{3,Float64}[]
    quads = MeshLibrary.Quad{Int}[]
    face_sets = MeshLibrary.MeshElementSet[]
    for axis in 1:3, s in (-1.0, 1.0)
        a, b = [k for k in 1:3 if k != axis]
        base = length(X)
        ξ = [tan(π / 4 * (2i / n - 1)) for i in 0:n]
        for j in 0:n, i in 0:n
            p = zeros(3)
            p[axis] = s
            p[a] = ξ[i+1]
            p[b] = ξ[j+1]
            push!(X, r * normalize(SVector{3}(p)))
        end
        id(i, j) = base + j * (n + 1) + i + 1
        first = length(quads) + 1
        for j in 0:n-1, i in 0:n-1
            f = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]
            push!(quads, MeshLibrary.Quad(orient_outward(X, f, sum(X[f]) / 4)...))
        end
        push!(face_sets, MeshLibrary.MeshElementSet("face_$(axis)$(s > 0 ? "+" : "-")", collect(first:length(quads))))
    end
    return MeshLibrary.Mesh(points(X), quads; element_sets=face_sets)
end

"""
    uv_sphere(nθ, nφ; r=1.0)

Latitude-longitude sphere of quads with every node repeated where it coincides: the poles appear nφ+1 times and the
φ = 0 and 2π meridians twice, so the pole quads collapse to triangles once welded. Quads are oriented outward.
"""
function uv_sphere(nθ::Integer, nφ::Integer; r::Real=1.0)
    X = SVector{3,Float64}[]
    for i in 0:nθ, j in 0:nφ
        θ, φ = π * i / nθ, 2π * j / nφ
        push!(X, r * SVector(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)))
    end
    id(i, j) = i * (nφ + 1) + j + 1
    quads = MeshLibrary.Quad{Int}[]
    for i in 0:nθ-1, j in 0:nφ-1
        f = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]
        # the centroid direction is outward even for the pole quads, two of whose nodes coincide
        push!(quads, MeshLibrary.Quad(orient_outward(X, f, sum(X[f]) / 4)...))
    end
    return MeshLibrary.Mesh(points(X), quads)
end

"""
    scale_nodes(mesh, s)

Copy of `mesh` with node coordinates scaled by the 3-vector `s` (e.g. a sphere into an ellipsoid) or a scalar.
"""
function scale_nodes(mesh::MeshLibrary.Mesh, s)
    sv = s isa Real ? SVector(s, s, s) : SVector{3}(s)
    nodes = [MeshLibrary.Point(sv .* p.coord) for p in mesh.nodes]
    return MeshLibrary.Mesh(MeshLibrary.StructVector(nodes), deepcopy(mesh.elements); boundaries=deepcopy(mesh.boundaries),
        element_sets=deepcopy(mesh.elemset))
end

"""Copy of `mesh` with coordinates converted to `T` (e.g. Float32)"""
function convert_nodes(T::Type{<:Real}, mesh::MeshLibrary.Mesh)
    nodes = [MeshLibrary.Point(SVector{3,T}(p.coord)) for p in mesh.nodes]
    return MeshLibrary.Mesh(MeshLibrary.StructVector(nodes), deepcopy(mesh.elements); boundaries=deepcopy(mesh.boundaries),
        element_sets=deepcopy(mesh.elemset))
end

"""Copy of `mesh` with every element's node order reversed (inward normals); boundaries are remapped to the reversed sides"""
function reverse_elements(mesh::MeshLibrary.Mesh)
    elems = [reverse(e) for e in mesh.elements]
    bounds = isempty(mesh.boundaries) ? copy(mesh.boundaries) : map(mesh.boundaries) do b
        facets = map(b.facets) do f
            n = length(MeshLibrary.get_nodes(mesh.elements[f.element]))
            # side s ran from node s to s+1; reversed, that edge is side n - s (from node n-s to n-s+1)
            MeshLibrary.BoundaryFacet(f.element, mod1(n - f.side, n))
        end
        MeshLibrary.MeshBoundary(b.name, collect(facets))
    end
    return MeshLibrary.Mesh(MeshLibrary.StructVector(deepcopy(mesh.nodes)), MeshLibrary.StructVector(elems);
        boundaries=bounds, element_sets=deepcopy(mesh.elemset))
end

"""Mirror image of `mesh` about the plane x[axis] = 0, with elements reversed so normals stay outward"""
function mirror_nodes(mesh::MeshLibrary.Mesh; axis::Integer=2)
    m = reverse_elements(mesh)
    flipv = SVector{3,Float64}([k == axis ? -1.0 : 1.0 for k in 1:3])
    nodes = [MeshLibrary.Point(flipv .* p.coord) for p in m.nodes]
    return MeshLibrary.Mesh(MeshLibrary.StructVector(nodes), m.elements; boundaries=m.boundaries, element_sets=m.elemset)
end
