const CONSTRAINT_CODE = Dict(:free => 0, :tjunction => 1, :gap => 2, :zip => 3, :tip => 4)

"""
    MeshLibrary.write_vtk(dm::DOFMapping, g::SurfaceGeometry, prefix)

Write `prefix_slots.vtu` (a vertex per slot with its constraint kind: 0 free, 1 T-junction, 2 gap, 3 zip, 4 tip; its
trailing edge side: 0 none, 1 upper, 2 lower; and ∫φ dA) and, when there is a trailing edge, `prefix_wake_te.vtu` (the
wake trailing edge lines with their tangent and bisector). Returns the file names.
"""
function MeshLibrary.write_vtk(dm::DOFMapping, g::SurfaceGeometry, prefix::AbstractString)
    pts = reduce(hcat, g.x[dm.slot_node])
    verts = [MeshCell(VTKCellTypes.VTK_VERTEX, [s]) for s in 1:dm.nslot]
    f1 = vtk_grid(prefix * "_slots", pts, verts) do vtk
        vtk["constraint", VTKCellData()] = [CONSTRAINT_CODE[c] for c in dm.constraint]
        vtk["te_side", VTKCellData()] = Int.(UInt8.(dm.slot_side))
        vtk["area", VTKCellData()] = dm.slot_area
    end
    w = dm.wake_te
    nwake_te(w) == 0 && return [f1]
    lines = [MeshCell(VTKCellTypes.VTK_POLY_LINE, c ? [line; line[1]] : line) for (line, c) in zip(w.lines, w.closed)]
    f2 = vtk_grid(prefix * "_wake_te", reduce(hcat, w.x), lines) do vtk
        vtk["tangent", VTKPointData()] = reduce(hcat, w.tangent)
        vtk["bisector", VTKPointData()] = reduce(hcat, w.bisector)
    end
    return [f1; f2]
end
