"""
    write_free_wake_vtk(sim, prefix)

Write `prefix_free.vtu`, the free panel rows (their triangles, with each node's frozen doublet and row, 0 the newest),
and `prefix_particles.vtu`, the particles (their strength vectors and cores). Returns the file names.
"""
function write_free_wake_vtk(sim::FreeWakeSim, prefix::AbstractString)
    st, topo = sim.state, sim.topo
    files = String[]
    nr = length(st.x)
    if nr >= 2
        pts = reduce(hcat, reduce(vcat, st.x))
        id(i, k) = (i - 1) * topo.nte + k
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, c) for i in 1:nr-1 for (A, B) in topo.segments
                 for c in ([id(i, A), id(i + 1, A), id(i, B)], [id(i + 1, A), id(i + 1, B), id(i, B)])]
        append!(files, vtk_grid(prefix * "_free", pts, cells) do vtk
            vtk["mu", VTKPointData()] = reduce(vcat, st.μ)
            vtk["row", VTKPointData()] = repeat(0:nr-1, inner=topo.nte)
        end)
    end
    if !isempty(st.particles)
        pts = reduce(hcat, [p.x for p in st.particles])
        cells = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in eachindex(st.particles)]
        append!(files, vtk_grid(prefix * "_particles", pts, cells) do vtk
            vtk["alpha", VTKPointData()] = reduce(hcat, [p.α for p in st.particles])
            vtk["core", VTKPointData()] = [p.δ for p in st.particles]
        end)
    end
    return files
end
