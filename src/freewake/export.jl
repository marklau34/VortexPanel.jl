"""
    write_free_wake_vtk(sim, prefix)

Write `prefix_free.vtu`, the free panel rows (their triangles, with each node's frozen doublet and row, 0 the newest),
and `prefix_particles.vtu`, the particles (their strength vectors, magnitudes and cores). Returns the file names.
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
            vtk["strength", VTKPointData()] = [norm(p.α) for p in st.particles]
            vtk["core", VTKPointData()] = [p.δ for p in st.particles]
            vtk["pending", VTKPointData()] = Int.(st.tag .> 0) # 1: a fine particle waiting for its strip's batch to merge
        end)
    end
    return files
end

"""
    write_body_vtk(sim, prefix)

Write `prefix_body.vtu`, the body of the last solved step: the pressure coefficient on each triangle ("Cp", unsteady
Bernoulli) and area-averaged to the nodes ("Cp_nodes"), and the doublet at the nodes ("mu", the upper copy at the
trailing edge); and `prefix_te_panel.vtu`, the trailing edge panel with its doublets ("mu", μ_w). Returns the file names.
"""
function write_body_vtk(sim::FreeWakeSim, prefix::AbstractString)
    g = sim.g
    nn = length(g.x)
    body = vtk_grid(prefix * "_body", reduce(hcat, g.x), [MeshCell(VTKCellTypes.VTK_TRIANGLE, collect(t)) for t in g.tris]) do vtk
        vtk["Cp", VTKCellData()] = sim.cp
        vtk["Cp_nodes", VTKPointData()] = nodal_average(g, sim.els, sim.cp)
        vtk["mu", VTKPointData()] = sim.μ[1:nn]
    end
    te = vtk_grid(prefix * "_te_panel", reduce(hcat, sim.te.x), [MeshCell(VTKCellTypes.VTK_TRIANGLE, collect(t)) for t in sim.te.tris]) do vtk
        vtk["mu", VTKPointData()] = sim.μ_w[sim.te.te]
    end
    return [body; te]
end

"""
    FreeWakeSeries(dir; every=1)

A ParaView time series of a free wake run, written to `dir` by `run!(sim, n; series)`: every `every` steps (and the
last), the body with its pressure (`write_body_vtk`), the trailing edge panel, the free rows and the particles
(`write_free_wake_vtk`), as the solve of that step saw them, and the collections body.pvd, te_panel.pvd, free.pvd and
particles.pvd that tie the frames to their times.
"""
struct FreeWakeSeries
    dir::String
    every::Int
    frames::Dict{String,Vector{Tuple{Float64,String}}}
end

FreeWakeSeries(dir::AbstractString; every::Integer=1) = (mkpath(dir); FreeWakeSeries(String(dir), every, Dict{String,Vector{Tuple{Float64,String}}}()))

"""Add files to the series at time t, each to the collection its suffix names (body, te_panel, free, particles)"""
function add_frame!(s::FreeWakeSeries, t::Real, files)
    for f in files
        name = something(findfirst(n -> endswith(f, "_$n.vtu"), ("te_panel", "body", "free", "particles")), 0)
        key = name == 0 ? splitext(basename(f))[1] : ("te_panel", "body", "free", "particles")[name]
        push!(get!(s.frames, key, Tuple{Float64,String}[]), (Float64(t), f))
    end
    return s
end

"""Write the series' .pvd collections (ParaView reads the frames' times from them); returns their file names"""
function write_collections(s::FreeWakeSeries)
    return map(collect(s.frames)) do (name, entries)
        file = joinpath(s.dir, name * ".pvd")
        open(file, "w") do io
            println(io, "<?xml version=\"1.0\"?>")
            println(io, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
            println(io, "  <Collection>")
            for (t, f) in entries
                println(io, "    <DataSet timestep=\"$t\" group=\"\" part=\"0\" file=\"$(basename(f))\"/>")
            end
            println(io, "  </Collection>")
            println(io, "</VTKFile>")
        end
        file
    end
end
