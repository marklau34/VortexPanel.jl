# Validation stage 5: the rectangular reference wing (NACA 0015, aspect ratio 3.88) with a rigid wake, at α = 5° as the
# MESA baseline (validation/morino_baseline). Lift, spanwise loading and induced drag by both force routes (pressure
# integration; Kutta–Joukowski and the Trefftz plane), under refinement, on clean and abused meshes (doc §11).
#
#   julia -t auto --project=validation validation/wing/run_wing.jl [levels=1,2,3]
#
# 1. Convergence of the capped wing (the default: caps, wake along the freestream to 20 spans) with Richardson
#    extrapolation, the tip loading, cond₁ and timings.
# 2. Tip treatments: caps (their corners at the trailing edge's end node on the upper or the lower copy), zipped and
#    open tips. The solver models no filaments on the body (a filament on a test triangle's edge makes the Galerkin
#    integral diverge), so zipped and open tips leave their boundary vorticity out.
# 3. Wake options: bisector first row, wake length, semi-infinite strips, first row, growth.
# 4. Quadrature sensitivity, and zero lift at α = 0 (the tips' up-down asymmetry shows there).
# 5. Abused meshes: the MESA smoke set (jitter, slivers, gaps, T-junctions, spanwise skew, misaligned and non-matching
#    trailing edges) and all at once.
# 6. MESA's baseline next to ours (validation/output/morino/summary.csv, when it exists).
# Writes validation/output/wing/ (summary CSV, spanwise loading, chordwise Cp at four stations, VTK of the body and wake).

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs, MeshLibrary, PAGE
using VortexPanel, VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "wing")
const ALPHA = 5.0
const REFERENCE = AssemblyOptions(far_ratio=8.0, near_ratio=4.0, near_depth=5, touch_levels=10, touch_order=4, kernel_far=Inf)
const CP_STATIONS = (0.0, 0.5, 0.8, 0.95)
const COLUMNS = [:CL, :CL_KJ, :CDi_p, :CDi_T, :Cm_c4]

function parse_args(args)
    levels = [1, 2, 3]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
    end
    return levels
end

wing(level) = reference_wing(level; split_te=false)

const HEADER = ["case" "level" "ntri" "nq" "CL" "CL_KJ" "CDi_p" "CDi_T" "CY" "Cm_LE" "Cm_c4" "tip_loading" "lambda" "colsum" "cond1" "t_assembly" "t_solve"]
row(name, level, r) = [name level r.ntri r.nq r.CL r.CL_KJ r.CDi_p r.CDi_T r.CY r.Cm_LE r.Cm_c4 r.tip r.λ r.columns r.cond r.t_assembly r.t_solve]

const TABLE = "%-38s %6s %9s %9s %10s %10s %9s %7s %9s %9s %9s %7s\n"
print_header() = @printf("%-38s %6s %9s %9s %10s %10s %9s %7s %9s %9s %9s %7s\n", "case", "nq", "CL", "CL_KJ", "CDi_p", "CDi_T", "Cm_c4", "tip",
    "λ", "colsum", "cond₁", "time")
print_row(name, r) = @printf("%-38s %6d %9.5f %9.5f %10.6f %10.6f %9.5f %7.3f %9.1e %9.1e %9.1e %6.1fs\n", name, r.nq, r.CL, r.CL_KJ, r.CDi_p,
    r.CDi_T, r.Cm_c4, r.tip, r.λ, r.columns, r.cond, r.t_assembly + r.t_solve)


"""MESA's baseline rows (validation/output/morino/summary.csv), if it exists"""
function mesa_baseline()
    file = joinpath(@__DIR__, "..", "output", "morino", "summary.csv")
    isfile(file) || return nothing
    d, h = readdlm(file, ',', header=true, comments=true)
    col(name) = findfirst(==(name), vec(h))
    return Dict(string(d[i, col("case")]) => (CL=d[i, col("CL")], CL_KJ=d[i, col("CL_KJ")], CDi_p=d[i, col("CDi_p")], CDi_T=d[i, col("CDi_T")],
        Cm_c4=d[i, col("Cm_c4")], cond=d[i, col("cond1_raw")]) for i in axes(d, 1))
end

function main(args=ARGS)
    levels = parse_args(args)
    mkpath(OUTDIR)
    println("Reference wing at α = $(ALPHA)°, threads $(Threads.nthreads())")
    rows = Any[]

    println("\n1. Convergence: capped tips, wake along the freestream to 20 spans")
    print_header()
    conv = Dict{Int,Any}()
    for l in levels
        r = solve_wing(wing(l); α=ALPHA)
        conv[l] = r
        print_row("L$l cap", r)
        push!(rows, row("cap", l, r))
        load = r.load
        writedlm(joinpath(OUTDIR, "spanwise_L$l.csv"), vcat(["s0" "s1" "eta" "Gamma" "cl_kj" "cl_p" "cd_p"],
            hcat(load.s0, load.s1, (load.s0 + load.s1) / WING_SPAN .- 1, load.Γ, load.lift_kj / WING_CHORD, load.lift_p / WING_CHORD,
                load.drag_p / WING_CHORD)), ',')
        for η in CP_STATIONS
            up, low = chordwise_cp(r, l, η)
            writedlm(joinpath(OUTDIR, @sprintf("cp_L%d_eta%.2f.csv", l, η)), vcat(["x_c" "cp" "surface"], up, low), ',')
        end
    end
    if length(levels) >= 3
        println("Richardson extrapolation (levels $(join(levels[end-2:end], ", ")))")
        for q in COLUMNS
            vals = [getfield(conv[l], q) for l in levels[end-2:end]]
            p, ext = richardson(vals)
            @printf("  %-6s %s   p = %5.2f, extrapolated %.6f\n", q, join([@sprintf("%.6f", v) for v in vals], "  "), p, ext)
        end
    end

    lv = levels[1:min(end, 2)]
    println("\n2. Tip treatments")
    print_header()
    for l in lv, (name, kw) in (("cap, upper copy", (;)), ("cap, lower copy", (dofmap=DOFMapOptions(cap_te_copy=:lower),)),
        ("zip", (geometry=GeometryOptions(openings=:zip),)), ("open", (geometry=GeometryOptions(openings=:open),)))
        r = name == "cap, upper copy" ? conv[l] : solve_wing(wing(l); α=ALPHA, kw...)
        print_row("L$l $name", r)
        name == "cap, upper copy" || push!(rows, row(name, l, r))
    end

    println("\n3. Wake options (L$(lv[end]))")
    print_header()
    for (name, o) in (("bisector first row", WakeOptions(direction=:bisector)), ("50 spans", WakeOptions(length=50.0)),
        ("semi-infinite strips", WakeOptions(ending=:semi_infinite)), ("first row 0.045 (MESA's)", WakeOptions(first_row=0.045)),
        ("growth 1.1", WakeOptions(growth=1.1)), ("growth 1.5", WakeOptions(growth=1.5)))
        r = solve_wing(wing(lv[end]); α=ALPHA, wake=o)
        print_row(name, r)
        push!(rows, row(name, lv[end], r))
    end

    println("\n4. Quadrature (L$(lv[1])) and zero lift at α = 0")
    print_header()
    for (name, o) in (("reference quadrature", REFERENCE), ("far ratio 6", AssemblyOptions(far_ratio=6.0)),
        ("touch levels 8", AssemblyOptions(touch_levels=8)), ("near depth 5", AssemblyOptions(near_depth=5)))
        r = solve_wing(wing(lv[1]); α=ALPHA, assembly=o)
        print_row(name, r)
        push!(rows, row(name, lv[1], r))
    end
    for l in lv
        r = solve_wing(wing(l); α=0.0)
        print_row("L$l α = 0", r)
        @printf("  L%d at α = 0: |μ_w| at most %.2e (at a tip: %.2e)\n", l, maximum(abs, r.μ_w), abs(r.μ_w[first(r.dm.wake_te.lines[1])]))
        push!(rows, row("alpha 0", l, r))
    end

    println("\n5. Abused meshes (seed 1), lift change from the clean mesh")
    @printf("%-42s %6s %9s %9s %8s %10s %7s %9s %7s   %s\n", "mesh", "level", "CL", "CL_KJ", "ΔCL %", "CDi_T", "tip", "cond₁", "time", "geometry")
    for l in lv
        clean = conv[l]
        for (name, mesh) in abused_wings(l)
            r = solve_wing(mesh; α=ALPHA)
            q = r.g.quality.counts
            @printf("%-42s %6d %9.5f %9.5f %8.3f %10.6f %7.3f %9.1e %6.1fs   TJ %d, gaps %d, non-matching TE %d, slivers %d\n", name, l, r.CL, r.CL_KJ,
                100 * (r.CL / clean.CL - 1), r.CDi_T, r.tip, r.cond, r.t_assembly + r.t_solve, q[:tjunctions], q[:gaps], q[:kutta_pairs], q[:slivers])
            push!(rows, row("abused: " * replace(name, "," => ";"), l, r))
        end
    end

    mesa = mesa_baseline()
    if !isnothing(mesa)
        println("\n6. MESA baseline (Morino, capped and open tips, 20-span wake along the freestream) next to ours (capped)")
        @printf("%-6s %-24s %9s %9s %10s %10s %9s %10s\n", "level", "solver", "CL", "CL_KJ", "CDi_p", "CDi_T", "Cm_c4", "cond₁")
        for l in levels, (label, key) in (("MESA, capped", "L$(l)_capped_w54_freestream"), ("MESA, open", "L$(l)_open_w54_freestream"))
            haskey(mesa, key) || continue
            m = mesa[key]
            @printf("%-6s %-24s %9.5f %9.5f %10.6f %10.6f %9.5f %10.1e\n", "L$l", label, m.CL, m.CL_KJ, m.CDi_p, m.CDi_T, m.Cm_c4, m.cond)
        end
        for l in levels
            r = conv[l]
            @printf("%-6s %-24s %9.5f %9.5f %10.6f %10.6f %9.5f %10.1e\n", "L$l", "Galerkin, capped", r.CL, r.CL_KJ, r.CDi_p, r.CDi_T, r.Cm_c4, r.cond)
        end
    end

    writedlm(joinpath(OUTDIR, "summary.csv"), vcat(HEADER, reduce(vcat, rows)), ',')
    r = conv[lv[end]]
    nn = length(r.g.x)
    MeshLibrary.write_vtk(r.g, joinpath(OUTDIR, "wing_L$(lv[end])"); cell_data=Dict("Cp" => r.cp),
        point_data=Dict("mu" => r.sol.μ[1:nn]))
    MeshLibrary.write_vtk(r.wake, joinpath(OUTDIR, "wing_L$(lv[end])"); point_data=Dict("mu" => r.μ_w[r.wake.te]))
    println("\nWrote ", normpath(OUTDIR))
    return conv
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
