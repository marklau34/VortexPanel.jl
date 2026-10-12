# Geometry layer and DOF map report on the reference wing family and abused variants: quality, repair and constraint
# tables, and VTK files of the triangulation (cell data: kind, body, parent, quality, trailing edge side), of the
# non-interior edges (class), of the slots (constraint kind, side, area) and of the wake trailing edge lines.
#
#   julia --project=validation validation/geometry/report_geometry.jl

using MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs
using VortexPanel
using VortexPanel.MeshTools
using VortexPanel.DOFMap

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "geometry")

function main()
    mkpath(OUTDIR)
    nspan, nchord = reference_wing_size(1)
    wing = reference_wing(1; split_te=false)
    project = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS)
    cases = Pair{String,Any}[]
    for level in 1:3
        push!(cases, "L$(level)_clean" => (reference_wing(level; split_te=true), GeometryOptions()))
    end
    push!(cases, "L1_open_tips" => (reference_wing(1), GeometryOptions(openings=:open)))
    push!(cases, "L1_zip_tips" => (reference_wing(1), GeometryOptions(openings=:zip)))
    push!(cases, "L1_abused" => (first(abuse(wing, (jitter=0.2, slivers=(2, 0.02), gaps=(2, 0.01), tjunctions=3, project=project); seed=1)), GeometryOptions()))
    push!(cases, "L1_lean_opposite_te" => (first(skew_wing(wing, :lean_opposite, 0.3; nspan=nspan, nchord=nchord, te_rows=true)), GeometryOptions()))
    push!(cases, "cube_sphere_8" => (cube_sphere(8), GeometryOptions()))

    rows = Any[]
    names = nothing
    for (name, (mesh, opts)) in cases
        g = build_geometry(mesh; opts=opts)
        dm = build_dofmap(g)
        println("\n== $name: ", g)
        print_report(g)
        print_dofmap(dm, g)
        MeshLibrary.write_vtk(g, joinpath(OUTDIR, name))
        MeshLibrary.write_vtk(dm, g, joinpath(OUTDIR, name))
        r = merge(report_row(g), (slots=dm.nslot, unknowns=dm.nq, tjunction_slots=count(==(:tjunction), dm.constraint),
            gap_slots=count(==(:gap), dm.constraint), zip_slots=count(==(:zip), dm.constraint), wake_te_nodes=nwake_te(dm.wake_te),
            stitch_offset=stitch_offset(dm, g)))
        names = collect(String.(keys(r)))
        push!(rows, permutedims(Any[name; collect(values(r))]))
    end
    writedlm(joinpath(OUTDIR, "geometry_report.csv"), vcat(["case" permutedims(names)], reduce(vcat, rows)), ',')
    println("\nWrote ", normpath(OUTDIR))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
