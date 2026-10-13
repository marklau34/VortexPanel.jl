# Validation stage 8: the test-weight comparison (doc §6, "test-weight switch"). The Galerkin rows ∫ψᵢ n·u dA = 0 with
# three test weights ψᵢ, on the reference wing (NACA 0015, aspect ratio 3.88, α = 5°, capped tips, rigid wake) with
# clean and abused meshes. Same quadrature points; only the weights change (AssemblyOptions(test=...)):
#   hat        the nodal hat function φᵢ (Galerkin, the default);
#   one-ring   one on every triangle of node i;
#   dual cell  one on its median dual cell (the part of each triangle where corner i has the largest barycentric
#              coordinate).
#
#   julia -t auto --project=validation validation/test_weights/run_test_weights.jl [levels=1,2,3] [section=19,38,76]
#
# 1. Refinement on clean meshes: CL, CDi and Cm by both force routes, the quadrature diagnostics (column defect, λ),
#    the condition number (raw and with the rows scaled by the diagonal) and the time.
# 2. Sectional accuracy: the midspan section of wings of aspect ratio 15.5, 31 and 62 (68 strips), extrapolated in
#    1/AR to the infinite wing, against the exact section (validation/common/panel2d.jl), at each chordwise resolution.
# 3. Chordwise Cp at midspan (L2), with the 2-D section at the same lift.
# 4. Abused meshes (MESA's smoke set and all at once, L1 and L2): the lift change from the clean mesh of the same weight.
# 5. Zero incidence on clean and abused meshes (L2).
# Writes validation/output/test_weights/ (CSV tables, SVG plots).

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs, MeshLibrary, PAGE
using VortexPanel, VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))
include(joinpath(@__DIR__, "..", "common", "panel2d.jl"))
include(joinpath(@__DIR__, "..", "common", "svgplot.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "test_weights")
const ALPHA = 5.0
const WEIGHTS = [(:hat, "hat"), (:onering, "one-ring"), (:dual, "dual cell")]
const QUANTITIES = (:CL, :CL_KJ, :CDi_p, :CDi_T, :Cm_c4)
const AR_FACTORS = (4, 8, 16) # spans of the sectional study, in reference spans

function parse_args(args)
    levels, section = [1, 2, 3], [19, 38, 76]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "section" && (section = parse.(Int, split(val, ",")))
    end
    return levels, section
end

pct(a, b) = 100 * (a / b - 1)

"""The wing `mesh` solved with test weight `test`"""
solve(mesh, test::Symbol; kw...) = solve_wing(mesh; assembly=AssemblyOptions(test=test), kw...)

# ---------------------------------------------------------------------------------------------------------------------
# 1. refinement

function refinement(levels)
    println("\n1. Refinement on clean meshes, (nspan, nchord) = (17, 19)·2^(level-1), α = $(ALPHA)°")
    @printf("%-10s %5s %9s %9s %10s %10s %9s %9s %9s %8s %8s %7s\n", "weight", "level", QUANTITIES..., "col. def.", "λ", "cond₁", "Jacobi", "time")
    res, rows = Dict{Tuple{Symbol,Int},Any}(), Any[]
    for l in levels, (test, label) in WEIGHTS
        r = solve(reference_wing(l; split_te=false), test; α=ALPHA, jacobi=true)
        res[(test, l)] = r
        v = [getproperty(r, q) for q in QUANTITIES]
        @printf("%-10s %5d %9.5f %9.5f %10.6f %10.6f %+9.5f %9.1e %9.1e %8.1e %8.1e %6.1fs\n", label, l, v..., r.columns, r.λ, r.cond, r.cond_jacobi,
            r.t_assembly + r.t_solve)
        push!(rows, [label l r.nq v... r.columns r.λ r.cond r.cond_jacobi r.t_assembly r.t_solve])
        GC.gc()
    end
    if length(levels) >= 3
        ls = levels[end-2:end]
        println("Richardson extrapolation from levels $(join(ls, ", ")) (p the observed order; NaN when not monotone or p < 0.3)")
        for (test, label) in WEIGHTS, q in QUANTITIES
            v = [getproperty(res[(test, l)], q) for l in ls]
            p, e = extrapolate(v)
            @printf("  %-10s %-6s %s   p = %5.2f, extrapolated %.6f\n", label, q, join([@sprintf("%.6f", x) for x in v], "  "), p, e)
        end
    end
    writedlm(joinpath(OUTDIR, "refinement.csv"), vcat(["weight" "level" "unknowns" "CL" "CL_KJ" "CDi_p" "CDi_T" "Cm_c4" "column_defect" "lambda" "cond1" "cond1_jacobi" "t_assembly" "t_solve"],
        reduce(vcat, rows)), ',')
    h = [2.0^-(l - 1) for l in levels]
    for (q, qk, title, ylabel) in ((:CL, :CL_KJ, "Lift under refinement by test weight", "CL"), (:CDi_p, :CDi_T, "Induced drag under refinement by test weight", "CDi"))
        series = Any[]
        for (i, (test, label)) in enumerate(WEIGHTS)
            push!(series, (x=h, y=[getproperty(res[(test, l)], q) for l in levels], label="$label, pressure", color=SVG_COLORS[i]))
            push!(series, (x=h, y=[getproperty(res[(test, l)], qk) for l in levels], label="$label, " * (qk == :CL_KJ ? "Kutta-Joukowski" : "Trefftz"),
                color=SVG_COLORS[i], dash=true))
        end
        svgplot(joinpath(OUTDIR, "refinement_$q.svg"); series=series, title=title, xlabel="mesh spacing relative to L1", ylabel=ylabel)
    end
    return res
end

# ---------------------------------------------------------------------------------------------------------------------
# 2. sectional accuracy

function section_accuracy(nchords)
    ex = panel2d(naca00_section(1600), ALPHA)
    @printf("\n2. Sectional accuracy: midspan of wings of 68 strips at AR 15.5, 31 and 62, extrapolated to AR = ∞, against the exact section\n")
    @printf("   (cl %.5f, cm_c4 %.5f); order from successive chordwise levels\n", ex.cl_gamma, ex.cm)
    println("   A weight whose midspan lift is off by more than half is reported as unstable (its largest cond₁ shows why).")
    @printf("%-10s %7s %10s %8s %10s %8s %10s %8s %9s %8s\n", "weight", "nchord", "cl", "order", "cl (KJ)", "order", "cm", "order", "cond₁", "time")
    rows, err = Any[], Dict{Tuple{Symbol,Int},NTuple{3,Float64}}()
    stable(test) = all(abs(err[(test, n)][1]) < 0.5 for n in nchords)
    for (test, label) in WEIGHTS, (i, nchord) in enumerate(nchords)
        if i > 1 && !(abs(err[(test, nchords[1])][1]) < 0.5) # unstable on the coarsest: the finer meshes add nothing
            err[(test, nchord)] = (NaN, NaN, NaN)
            @printf("%-10s %7d %10s\n", label, nchord, "skipped")
            continue
        end
        t = @elapsed pts = map(AR_FACTORS) do m
            b = WING_SPAN * m
            r = solve(reference_wing(1; split_te=false, nspan=68, nchord=nchord, span=b), test; α=ALPHA, span=b)
            GC.gc()
            merge((AR=b / WING_CHORD, cond=r.cond), midspan_section(r))
        end
        ar = [p.AR for p in pts]
        cl, kj, cm = (extrapolate_ar(ar, [getproperty(p, q) for p in pts]) for q in (:cl_p, :cl_kj, :cm))
        e = (cl / ex.cl_gamma - 1, kj / ex.cl_gamma - 1, cm - ex.cm)
        err[(test, nchord)] = e
        cmax = maximum(p.cond for p in pts)
        if abs(e[1]) >= 0.5
            @printf("%-10s %7d %10s %8s %10s %8s %10s %8s %9.1e %7.0fs\n", label, nchord, "unstable", "", "", "", "", "", cmax, t)
        else
            order = i == 1 || abs(err[(test, nchords[i-1])][1]) >= 0.5 ? fill(NaN, 3) :
                    [log(abs(a / b)) / log(nchord / nchords[i-1]) for (a, b) in zip(err[(test, nchords[i-1])], e)]
            @printf("%-10s %7d %+9.2f%% %8.2f %+9.2f%% %8.2f %+10.5f %8.2f %9.1e %7.0fs\n", label, nchord, 100e[1], order[1], 100e[2], order[2], e[3], order[3],
                cmax, t)
        end
        push!(rows, [label nchord cl kj cm e... cmax])
    end
    writedlm(joinpath(OUTDIR, "section.csv"), vcat(["weight" "nchord" "cl" "cl_kj" "cm" "cl_error" "cl_kj_error" "cm_error" "max_cond1"], reduce(vcat, rows),
        ["exact" 0 ex.cl_gamma ex.cl_gamma ex.cm 0 0 0 NaN]), ',')
    for (k, name, ylabel) in ((1, "cl", "|cl error| [%]"), (3, "cm", "|cm error|"))
        series = [(x=1 ./ nchords, y=[abs(err[(test, n)][k]) * (k == 1 ? 100 : 1) for n in nchords], label=label, color=SVG_COLORS[i])
                  for (i, (test, label)) in enumerate(WEIGHTS) if stable(test)]
        svgplot(joinpath(OUTDIR, "section_$name.svg"); series=series, xlabel="1 / panels per surface", ylabel=ylabel,
            title="Sectional $name error against the exact section, stable test weights")
    end
    return ex, err
end

# ---------------------------------------------------------------------------------------------------------------------
# 3. chordwise pressure

function midspan_pressure(res, ex, level)
    println("\n3. Chordwise Cp at midspan (L$level): suction peak, and the largest difference from the hat weight between 10% and 90% chord")
    interp(x, y, xi) = (k = clamp(searchsortedlast(x, xi), 1, length(x) - 1); w = (xi - x[k]) / (x[k+1] - x[k]); (1 - w) * y[k] + w * y[k+1])
    hu, hl = chordwise_cp(res[(:hat, level)], level, 0.0)
    rows, series = Any[], Any[]
    for (i, (test, label)) in enumerate(WEIGHTS)
        up, low = chordwise_cp(res[(test, level)], level, 0.0)
        du = maximum(abs(interp(up[:, 1], up[:, 2], x) - interp(hu[:, 1], hu[:, 2], x)) for x in 0.1:0.02:0.9)
        dl = maximum(abs(interp(low[:, 1], low[:, 2], x) - interp(hl[:, 1], hl[:, 2], x)) for x in 0.1:0.02:0.9)
        @printf("   %-10s suction peak %.3f at x/c %.3f; |ΔCp| from hat at most %.3f (upper), %.3f (lower); trailing edge panels %.3f / %.3f\n", label,
            minimum(up[:, 2]), up[argmin(up[:, 2]), 1], du, dl, up[end, 2], low[end, 2])
        push!(rows, hcat(fill(label, size(up, 1) + size(low, 1)), vcat(up, low)))
        push!(series, (x=up[:, 1], y=-up[:, 2], label=label, color=SVG_COLORS[i], marker=false))
        push!(series, (x=low[:, 1], y=-low[:, 2], label="", color=SVG_COLORS[i], marker=false))
    end
    cl_mid = midspan_section(res[(:hat, level)]).cl_p
    αe = ALPHA * cl_mid / ex.cl_gamma
    for _ in 1:3
        αe *= cl_mid / panel2d(naca00_section(400), αe).cl_gamma
    end
    sec = panel2d(naca00_section(400), αe)
    for upper in (true, false)
        pts = [(p[1], c) for (p, c) in zip(sec.x, sec.cp) if (p[2] >= 0) == upper]
        push!(series, (x=first.(pts), y=-last.(pts), label=upper ? @sprintf("2-D section at cl %.3f", cl_mid) : "", color="#52514e", dash=true, marker=false))
    end
    writedlm(joinpath(OUTDIR, "cp_L$(level)_midspan.csv"), vcat(["weight" "x_c" "cp" "surface"], reduce(vcat, rows)), ',')
    svgplot(joinpath(OUTDIR, "cp_L$(level)_midspan.svg"); series=series, title="Chordwise pressure at midspan by test weight, L$level",
        xlabel="x/c", ylabel="-Cp", xlims=(0.0, 1.0))
end

# ---------------------------------------------------------------------------------------------------------------------
# 4. abused meshes, 5. zero incidence

function abused(res, levels)
    println("\n4. Abused meshes (MESA's smoke set and all at once): lift change (%) from the clean mesh of the same weight")
    @printf("%-42s %6s %12s %12s %12s\n", "mesh", "level", (label for (_, label) in WEIGHTS)...)
    rows, summary = Any[], Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for l in levels
        for (name, mesh) in abused_wings(l)
            d = map(WEIGHTS) do (test, _)
                r = solve(mesh, test; α=ALPHA)
                push!(get!(summary, (test, l), Float64[]), pct(r.CL, res[(test, l)].CL))
                (pct(r.CL, res[(test, l)].CL), pct(r.CL_KJ, res[(test, l)].CL_KJ), r.columns)
            end
            @printf("%-42s %6d %+12.3f %+12.3f %+12.3f\n", name, l, first.(d)...)
            push!(rows, [replace(name, "," => ";") l reduce(hcat, [collect(x) for x in d])...])
            GC.gc()
        end
    end
    for l in levels
        @printf("   L%d, mean (largest) |ΔCL|: %s\n", l, join([@sprintf("%s %.3f%% (%.3f%%)", label, sum(abs, summary[(test, l)]) / length(summary[(test, l)]),
            maximum(abs, summary[(test, l)])) for (test, label) in WEIGHTS], ", "))
    end
    header = ["mesh" "level" reduce(hcat, [["$(t)_dCL_pct" "$(t)_dCL_KJ_pct" "$(t)_column_defect"] for (t, _) in WEIGHTS])]
    writedlm(joinpath(OUTDIR, "abused.csv"), vcat(header, reduce(vcat, rows)), ',')
    for l in levels
        n = length(summary[(:hat, l)])
        series = [(x=collect(1.0:n), y=summary[(test, l)], label=label, color=SVG_COLORS[i], line=false) for (i, (test, label)) in enumerate(WEIGHTS)]
        svgplot(joinpath(OUTDIR, "abused_L$l.svg"); series=series, hline=0.0, ylims=(-3.0, 3.0), xlabel="abused mesh (the order of abused.csv)",
            ylabel="ΔCL [%]", title="Lift change on abused meshes by test weight, L$l (beyond ±3%: clipped)")
    end
end

function zero_incidence(level)
    println("\n5. Zero incidence (L$level): CL by test weight")
    meshes = ["clean" => reference_wing(level; split_te=false); [p for p in abused_wings(level) if p.first in
        ("jitter independent 0.3", "lean opposite 0.3", "TE misaligned 0.3, non-matching", "all (jitter, slivers, gaps, T-junctions)")]]
    @printf("%-42s %12s %12s %12s\n", "mesh", (label for (_, label) in WEIGHTS)...)
    rows = Any[]
    for (name, mesh) in meshes
        cl = [solve(mesh, test; α=0.0).CL for (test, _) in WEIGHTS]
        @printf("%-42s %12.2e %12.2e %12.2e\n", name, cl...)
        push!(rows, [replace(name, "," => ";") cl...])
    end
    writedlm(joinpath(OUTDIR, "zero_incidence.csv"), vcat(["mesh" [string(t) for (t, _) in WEIGHTS]...], reduce(vcat, rows)), ',')
end

function main(args=ARGS)
    levels, section = parse_args(args)
    mkpath(OUTDIR)
    println("Test weights on the reference wing, threads $(Threads.nthreads())")
    solve(reference_wing(1; split_te=false), :dual) # compile first, so the timings exclude it
    res = refinement(levels)
    ex, _ = section_accuracy(section)
    lv = filter(<=(2), levels)
    isempty(lv) && (lv = levels[1:1])
    midspan_pressure(res, ex, maximum(lv))
    abused(res, lv)
    zero_incidence(maximum(lv))
    println("\nWrote ", normpath(OUTDIR))
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
