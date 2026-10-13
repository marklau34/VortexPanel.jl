"""
    MeshLibrary.write_vtk(w::RigidWake, prefix; point_data=Dict())

Write `prefix_wake.vtu` (the wake triangles with each node's row and trailing edge node, plus any extra `point_data`
per wake node, e.g. the doublet strength) and `prefix_filaments.vtu` (the tip filament segments). Returns the file
names.
"""
function MeshLibrary.write_vtk(w::RigidWake, prefix::AbstractString; point_data=Dict())
    pts = reduce(hcat, w.x)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, collect(t)) for t in w.tris]
    f1 = vtk_grid(prefix * "_wake", pts, cells) do vtk
        vtk["row", VTKPointData()] = w.row
        vtk["te_node", VTKPointData()] = w.te
        for (name, v) in point_data
            vtk[string(name), VTKPointData()] = v
        end
    end
    isempty(w.filaments) && return [f1]
    lines = [MeshCell(VTKCellTypes.VTK_LINE, [a, b]) for (a, b) in w.filaments]
    f2 = vtk_grid(prefix * "_filaments", pts, lines) do vtk
        vtk["te_node", VTKCellData()] = [w.te[a] for (a, _) in w.filaments]
    end
    return [f1; f2]
end
