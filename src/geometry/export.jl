"""
    MeshLibrary.write_vtk(g::SurfaceGeometry, prefix; cell_data=Dict(), point_data=Dict())

Write `prefix_tris.vtu` (triangles with their kind, body, parent element, quality and trailing edge side, and node
flags, plus any extra `cell_data` per triangle and `point_data` per node, e.g. a solution) and, when there are any,
`prefix_edges.vtu` (every edge that is not interior, as a line with its class). Returns the file names.
"""
function MeshLibrary.write_vtk(g::SurfaceGeometry, prefix::AbstractString; cell_data=Dict(), point_data=Dict())
    pts = reduce(hcat, g.x)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, collect(t)) for t in g.tris]
    te_side = [maximum(UInt8.(cs)) for cs in g.kutta.corner_side]
    node_te = zeros(Int, length(g.x))
    node_te[g.kutta.te_nodes] .= 1
    node_hanging = zeros(Int, length(g.x))
    for tj in g.tjunctions
        node_hanging[tj.hanging] = 1
    end
    multiplicity = zeros(Int, length(g.x))
    for v in g.input_node
        multiplicity[v] += 1
    end
    f1 = vtk_grid(prefix * "_tris", pts, cells) do vtk
        vtk["kind", VTKCellData()] = Int.(g.tri_kind)
        vtk["body", VTKCellData()] = g.tri_body
        vtk["parent", VTKCellData()] = g.tri_parent
        vtk["min_angle", VTKCellData()] = g.quality.min_angle
        vtk["aspect", VTKCellData()] = g.quality.aspect
        vtk["skew", VTKCellData()] = g.quality.skew
        vtk["sliver", VTKCellData()] = Int.(g.quality.sliver)
        vtk["te_side", VTKCellData()] = Int.(te_side)
        vtk["normal", VTKCellData()] = reduce(hcat, [tri_normal(g.x, t) for t in g.tris])
        vtk["trailing_edge", VTKPointData()] = node_te
        vtk["hanging", VTKPointData()] = node_hanging
        vtk["weld_multiplicity", VTKPointData()] = multiplicity
        vtkdata(v) = eltype(v) <: AbstractVector ? reduce(hcat, v) : v
        for (name, v) in cell_data
            vtk[string(name), VTKCellData()] = vtkdata(v)
        end
        for (name, v) in point_data
            vtk[string(name), VTKPointData()] = vtkdata(v)
        end
    end
    special = findall(!=(Interior), g.edges.class)
    isempty(special) && return f1
    lines = [MeshCell(VTKCellTypes.VTK_LINE, collect(g.edges.nodes[e])) for e in special]
    f2 = vtk_grid(prefix * "_edges", pts, lines) do vtk
        vtk["class", VTKCellData()] = Int.(g.edges.class[special])
    end
    return [f1; f2]
end

"""
    add_caps(mesh, g)

The input MeshLibrary mesh with the cap triangles of `g` appended as elements, as a mixed `Poly2D{Int,4}` mesh. Existing
element ids, boundaries and element sets are unchanged, and an element set `tip_caps` holds the caps. Cap nodes map
to the lowest input node of each welded node (the upper copy of a split trailing edge node), and the caps follow the
orientation of the input mesh. This feeds the same closed geometry to solvers that read MeshLibrary meshes.
"""
function add_caps(mesh::MeshLibrary.Mesh, g::SurfaceGeometry)
    PT = MeshLibrary.Poly2D{Int,4}
    rep = zeros(Int, length(g.x))
    for (i, v) in enumerate(g.input_node)
        rep[v] == 0 && (rep[v] = i)
    end
    # does g keep the input orientation? compare a body triangle with its parent element
    reversed = false
    for t in eachindex(g.tris)
        g.tri_kind[t] == BodyTri || continue
        parent = collect(MeshLibrary.get_nodes(mesh.elements[g.tri_parent[t]]))
        welded = [g.input_node[i] for i in parent]
        a, b = g.tris[t][1], g.tris[t][2]
        n = length(welded)
        reversed = !any(welded[k] == a && welded[mod1(k + 1, n)] == b for k in 1:n)
        break
    end
    elems = [PT(e) for e in mesh.elements]
    caps = findall(==(CapTri), g.tri_kind)
    for t in caps
        a, b, c = rep[g.tris[t][1]], rep[g.tris[t][2]], rep[g.tris[t][3]]
        push!(elems, reversed ? PT(a, c, b) : PT(a, b, c))
    end
    sets = deepcopy(mesh.elemset)
    push!(sets, MeshLibrary.MeshElementSet("tip_caps", collect(length(mesh.elements)+1:length(elems))))
    return MeshLibrary.Mesh(copy(mesh.nodes), MeshLibrary.StructVector(elems); boundaries=deepcopy(mesh.boundaries), element_sets=sets)
end
