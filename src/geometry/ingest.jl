"""
Polygonal surface as given, before any repair: node coordinates, faces (node lists of 3 or more nodes), named edge
tags (node pairs, from mesh boundaries) and named face tags (element sets).
"""
struct RawSurface{T}
    x::Vector{SVector{3,T}}
    faces::Vector{Vector{Int}}
    edge_tags::Dict{String,Vector{Tuple{Int,Int}}}
    face_tags::Dict{String,Vector{Int}}
end

"""
    ingest(mesh::MeshLibrary.Mesh)

Read the nodes, the Tri, Quad and Poly2D elements, the boundaries (as edge tags) and the element sets (as face tags)
of a MeshLibrary mesh.
"""
function ingest(mesh::MeshLibrary.Mesh)
    x = [SVector{3}(p.coord) for p in mesh.nodes]
    faces = [Vector{Int}(collect(MeshLibrary.get_nodes(e))) for e in mesh.elements]
    edge_tags = Dict{String,Vector{Tuple{Int,Int}}}()
    for b in mesh.boundaries
        pairs = Tuple{Int,Int}[]
        for f in b.facets
            con = faces[f.element]
            push!(pairs, (con[f.side], con[mod1(f.side + 1, length(con))]))
        end
        edge_tags[b.name] = pairs
    end
    face_tags = Dict{String,Vector{Int}}(s.name => collect(Int, s.elems) for s in mesh.elemset)
    return RawSurface(x, faces, edge_tags, face_tags)
end

"""
    ingest(nodes, faces; edge_tags, face_tags)

Raw surface from node coordinates (any 3-vectors) and faces (node index lists of 3 or more nodes).
"""
function ingest(nodes::AbstractVector, faces::AbstractVector; edge_tags=Dict{String,Vector{Tuple{Int,Int}}}(),
    face_tags=Dict{String,Vector{Int}}())
    x = [SVector{3}(p[1], p[2], p[3]) for p in nodes]
    return RawSurface(x, [Vector{Int}(collect(f)) for f in faces], Dict{String,Vector{Tuple{Int,Int}}}(edge_tags),
        Dict{String,Vector{Int}}(face_tags))
end
