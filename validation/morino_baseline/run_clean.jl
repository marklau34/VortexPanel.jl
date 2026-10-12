# MESA.jl rigid wake baseline on the clean reference wing family (MESA.jl/examples/wing/run_wing_surfpan.jl, in free air).
#
#   julia -t auto --project=validation validation/morino_baseline/run_clean.jl [levels=1,2,3] [tips=open,capped]
#
# Primary configuration: wake_len = 20 spans, trailing edge panels along the freestream (the Galerkin solver's default
# wake). At levels 1-2 also MESA's example wake length (8 m) and MESA's default trailing edge bisector direction.
# Capped meshes close the open tips with VortexPanel's no-new-node caps (`add_caps`), the same caps the Galerkin
# solver uses. Writes validation/output/morino/.

using MESA, MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles
using VortexPanel

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "common_mesa.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "morino")
const WAKE_LEN = 20 * WING_SPAN
const CP_STATIONS = (0.0, 0.5, 0.8, 0.95)

function parse_args(args)
    levels, tips = [1, 2, 3], [:open, :capped]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "tips" && (tips = Symbol.(split(val, ",")))
    end
    return levels, tips
end

"""Short git hash of a checkout, with +dirty when it has uncommitted changes"""
function git_sha(path)
    try
        sha = readchomp(pipeline(`git -C $path rev-parse --short HEAD`, stderr=devnull))
        dirty = !isempty(readchomp(pipeline(`git -C $path status --porcelain --untracked-files=no`, stderr=devnull)))
        return sha * (dirty ? "+dirty" : "")
    catch
        return "none"
    end
end

function provenance()
    gh = joinpath(@__DIR__, "..", "..", "..")
    return ["MESA.jl" => git_sha(joinpath(gh, "MESA.jl")), "MeshLibrary.jl" => git_sha(joinpath(gh, "MeshLibrary.jl")),
        "_tools_release (PAGE)" => git_sha(joinpath(gh, "_tools_release")), "VortexPanel.jl" => git_sha(joinpath(@__DIR__, "..", "..")),
        "julia" => string(VERSION), "threads" => string(Threads.nthreads())]
end

"""Mesh of one case: the PAGE wing (split trailing edge, as MESA needs), with the tips open or capped"""
function case_mesh(level::Integer, tips::Symbol)
    mesh = reference_wing(level; split_te=true)
    tips == :open && return mesh
    tips == :capped && return VortexPanel.Geometry.add_caps(mesh, VortexPanel.Geometry.build_geometry(mesh))
    error("Unknown tips $tips")
end

function write_case_files(name, r)
    writedlm(joinpath(OUTDIR, "spanwise_$name.csv"), vcat(["eta" "y" "dy" "cl" "cd" "mu_w" "Gamma"], r.spanwise), ',')
    te = hcat(r.y_edges[1:end-1], r.y_edges[2:end], r.μ_w, -r.μ_w)
    writedlm(joinpath(OUTDIR, "te_mu_$name.csv"), vcat(["y1" "y2" "mu_w" "Gamma"], te), ',')
    for η in CP_STATIONS
        up = interp_stations(r.η, r.cp_up, η)
        low = interp_stations(r.η, r.cp_low, η)
        xc = interp_stations(r.η, r.xc_up, η)
        writedlm(joinpath(OUTDIR, @sprintf("cp_%s_eta%.2f.csv", name, η)), vcat(["x_c" "cp_upper" "cp_lower"], hcat(xc, up, low)), ',')
    end
end

function main(args=ARGS)
    levels, tips_modes = parse_args(args)
    mkpath(OUTDIR)

    cases = []
    for tips in tips_modes, level in levels
        push!(cases, (level=level, tips=tips, wake_len=WAKE_LEN, te_dir=:freestream))
        if level <= 2
            push!(cases, (level=level, tips=tips, wake_len=8.0, te_dir=:freestream))
            push!(cases, (level=level, tips=tips, wake_len=WAKE_LEN, te_dir=:bisector))
        end
    end

    # warm up so the timings exclude compilation
    mesa_case(case_mesh(1, :open); wake_len=WAKE_LEN, nspan=17, nchord=19)

    header = ["case" "level" "tips" "wake_len" "te_dir" "npan" "CL" "CL_KJ" "CDi_p" "CDi_T" "CY" "Cm_LE" "Cm_c4" "alpha" "cond1_raw" "cond1_jacobi" "t_assembly" "t_lu" "t_solve" "check_mu"]
    rows = Any[]
    results = Dict()
    @printf("\n%-34s %6s %8s %8s %9s %9s %8s %10s %7s\n", "case", "npan", "CL", "CL_KJ", "CDi_p", "CDi_T", "Cm_c4", "cond1", "t [s]")
    for c in cases
        nspan, nchord = reference_wing_size(c.level)
        name = @sprintf("L%d_%s_w%.0f_%s", c.level, c.tips, c.wake_len, c.te_dir)
        vtk_dir = joinpath(OUTDIR, "vtk", name)
        mkpath(vtk_dir)
        r = mesa_case(case_mesh(c.level, c.tips); wake_len=c.wake_len, te_dir=c.te_dir, nspan=nspan, nchord=nchord, vtk_dir=vtk_dir)
        results[c] = r
        write_case_files(name, r)
        push!(rows, [name c.level string(c.tips) c.wake_len string(c.te_dir) r.npan r.CL r.CL_KJ r.CDi r.CDi_T r.CY r.Cm_LE r.Cm_c4 r.alpha r.cond_raw r.cond_jacobi r.timings.assembly r.timings.lu r.timings.solve r.check])
        @printf("%-34s %6d %8.5f %8.5f %9.6f %9.6f %8.5f %10.3e %7.2f\n", name, r.npan, r.CL, r.CL_KJ, r.CDi, r.CDi_T, r.Cm_c4, r.cond_raw,
            r.timings.assembly + r.timings.solve)
        r.check > 1e-10 && @warn "Rebuilt system does not reproduce MESA's solution" name r.check
    end

    open(joinpath(OUTDIR, "summary.csv"), "w") do io
        for (k, v) in provenance()
            println(io, "# ", k, ": ", v)
        end
        writedlm(io, vcat(header, reduce(vcat, rows)), ',')
    end

    # convergence of the primary configuration
    for tips in tips_modes
        prim = [results[(level=l, tips=tips, wake_len=WAKE_LEN, te_dir=:freestream)] for l in levels if haskey(results, (level=l, tips=tips, wake_len=WAKE_LEN, te_dir=:freestream))]
        length(prim) < 3 && continue
        println("\nRichardson extrapolation, $tips tips, wake 20b along the freestream (levels $(join(levels, ", ")))")
        for q in (:CL, :CL_KJ, :CDi, :CDi_T, :Cm_c4)
            vals = [getfield(r, q) for r in prim]
            p, ext = richardson(vals)
            @printf("  %-6s %s   p = %5.2f, extrapolated %.6f\n", q, join([@sprintf("%.6f", v) for v in vals], "  "), p, ext)
        end
    end
    println("\nWrote ", normpath(OUTDIR))
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
