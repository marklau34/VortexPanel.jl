# Validation stage 6: the Galerkin solver against the Morino code (MESA.jl) on the reference wing (doc §11,
# "reference-case comparison"), with an exact two-dimensional referee for the section. Both codes solve the same
# meshes: the PAGE reference wing (NACA 0015, aspect ratio 3.88) at α = 5°, the tips closed with the same caps (no new
# nodes), a rigid wake along the freestream to 20 spans.
#
#   julia -t auto --project=validation validation/comparison/run_comparison.jl [levels=1,2,3]
#
# 1. Refinement on clean meshes: CL (pressure and Kutta–Joukowski), CDi (pressure and Trefftz) and Cm about c/4, with
#    Richardson extrapolation, our zero-tip Kutta variant (tip_kutta = :zero), and the condition number (raw, and with
#    the rows scaled by the diagonal) and wall-clock time at matched resolution.
# 2. Refinement by direction from the L2 mesh, chordwise (19-152 panels per surface at 34 strips) and spanwise (17-136
#    strips at 38 panels): where each code's error comes from, and an estimate at infinite resolution.
# 3. The two-dimensional limit: the midspan section's cl and cm at aspect ratios 15.5, 31 and 62 (68 strips, the
#    chordwise panels of each level), extrapolated in 1/AR to the infinite wing, against an exact 2-D linear-vortex
#    panel solution of the section (validation/common/panel2d.jl): each code's sectional error.
# 4. Chordwise Cp at η = 0, 0.5, 0.8 and 0.95 and the spanwise loading on the L2 mesh, with the 2-D section at the
#    midspan's cl.
# 5. Skew sensitivity, the patterns of MESA/dev/skew_sensitivity.jl: spanwise lean of the same or opposite sign on the
#    two surfaces, or jitter, the same or independent, by σ = 0.05-0.3 of the station spacing (jitter: 5 seeds), and
#    the ratio of the codes' lift changes.
# 6. Trailing edge misalignment: the rows aft of 80% chord offset along the span by ±δ/2 with the trailing edge shared,
#    or the upper rows by δ with a non-matching trailing edge (which MESA cannot represent).
# 7. Zero incidence on skewed meshes: the symmetric section's lift must stay near zero.
# 8. Robustness: the abused meshes of MESA's smoke test (validation/morino_baseline/smoke_abuse.jl), capped: which code
#    runs, and the lift change.
# Writes validation/output/comparison/ (CSV tables, SVG plots, provenance).

using MESA, MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs, Logging
using VortexPanel, VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "..", "morino_baseline", "common_mesa.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))
include(joinpath(@__DIR__, "..", "common", "panel2d.jl"))
include(joinpath(@__DIR__, "..", "common", "svgplot.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "comparison")
const AR0 = WING_SPAN / WING_CHORD
const QUANTITIES = (:CL, :CL_KJ, :CDi_p, :CDi_T, :Cm_c4)
const CODES = [(:ours, "Galerkin"), (:zero, "Galerkin, zero tip"), (:mesa, "MESA (Morino)")]
const CODE_LABEL = Dict(CODES)
const CHORD_STUDY = [19, 38, 76, 152] # panels per surface, at 34 strips
const SPAN_STUDY = [17, 34, 68, 136]  # strips, at 38 panels per surface
const AR_FACTORS = (4, 8, 16)         # spans of the two-dimensional limit, in reference spans
const CP_STATIONS = (0.0, 0.5, 0.8, 0.95)
const PATTERNS = [(:lean_same, "lean, same"), (:lean_opposite, "lean, opposite"), (:jitter_same, "jitter, same"),
    (:jitter_independent, "jitter, independent")]
const SIGMAS = [0.05, 0.1, 0.2, 0.3]
const SEEDS = 1:5
const OFFSETS = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]

function parse_args(args)
    levels = [1, 2, 3]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
    end
    return levels
end

pct(a, b) = 100 * (a / b - 1)

"""
Richardson extrapolation of three values on meshes refined by 2 each time: the observed order and the extrapolated value,
NaN unless the order is at least 0.3 (below it the values are not in the asymptotic range, and 1/(2^p - 1) amplifies
their differences more than threefold)
"""
function extrapolate(vals)
    p, e = richardson(vals)
    return p, (isfinite(p) && p >= 0.3 ? e : NaN)
end

"""The exception's type and the MESA, PAGE or VortexPanel function it came from, when found"""
function error_line(e, bt)
    frames = [string(fr.func) for fr in stacktrace(bt)]
    k = findfirst(in(("calc_CHTLS_pan!", "split_trailing_edge!", "create_trailing_edge", "init", "solve!", "build_geometry",
            "build_dofmap", "build_wake", "assemble", "solve_galerkin")), frames)
    return string(nameof(typeof(e))) * (isnothing(k) ? "" : " in " * frames[k])
end

"""
    mesa_wing(mesh; nspan, nchord, alpha=ALPHA_DEG, span=WING_SPAN)

MESA's solution of a wing mesh in the `split_te=false` form (clean, skewed or abused): its trailing edge split as MESA
needs, the tips capped with our caps, solved at `alpha` with a wake of 20 spans. Returns the result and "", or nothing
and why: "n/a" when MESA cannot represent the mesh (the trailing edge split needs matching upper and lower facets),
"fails" when its solve throws or gives NaN.
"""
function mesa_wing(mesh; nspan::Integer, nchord::Integer, alpha::Real=ALPHA_DEG, span::Real=WING_SPAN)
    m = deepcopy(mesh)
    try
        PAGE.split_trailing_edge!(m, "trailing_edge")
    catch e
        return nothing, "n/a: " * error_line(e, catch_backtrace())
    end
    try
        capped = VortexPanel.Geometry.add_caps(m, build_geometry(m))
        ref = WingReference(span=span, chord=WING_CHORD, rho=RHO, vinf=VINF)
        r = with_logger(NullLogger()) do
            mesa_case(capped; wake_len=20span, nspan=nspan, nchord=nchord, ref=ref, alpha=alpha)
        end
        return isfinite(r.CL) ? (r, "") : (nothing, "fails: NaN lift")
    catch e
        return nothing, "fails: " * error_line(e, catch_backtrace())
    end
end

"""Our solution of a wing mesh (`solve_wing`) and "", or nothing and why"""
function ours_wing(mesh; kw...)
    try
        return solve_wing(mesh; kw...), ""
    catch e
        return nothing, "fails: " * error_line(e, catch_backtrace())
    end
end

"""A coefficient of either code's result by our names (MESA's pressure drag is `CDi`)"""
coef(r, q::Symbol) = q == :CDi_p && !hasproperty(r, :CDi_p) ? r.CDi : getproperty(r, q)

"""Relative change (%) of coefficient `q` from the clean result `r0`; NaN when the solve failed"""
change(r, r0, q::Symbol) = isnothing(r) ? NaN : pct(coef(r, q), coef(r0, q))

"""Change of coefficient `q` from the clean result `r0`; NaN when the solve failed"""
delta(r, r0, q::Symbol) = isnothing(r) ? NaN : coef(r, q) - coef(r0, q)

"""
The diagonals the default split chose for the quads of `mesh`: with `GeometryOptions(split=:pattern, split_pattern=...)`
a perturbed copy of the mesh keeps the triangulation, which separates the effect of panel skew from that of diagonal
flips
"""
function split_pattern(mesh)
    g = build_geometry(mesh)
    p = zeros(Int8, MeshLibrary.nelements(mesh))
    p[g.quad_parent] .= g.quad_diag
    return p
end

"""Our solution of `mesh` with the triangulation `pattern` (`split_pattern` of the clean mesh)"""
ours_clean_split(mesh, pattern; kw...) = ours_wing(mesh; geometry=GeometryOptions(split=:pattern, split_pattern=pattern), kw...)

"""
The clean reference wing with `nspan` strips and `nchord` panels per surface solved by each of `codes` (:ours, :zero
for our zero-tip variant, :mesa), as a Dict
"""
function clean_wing(nspan::Integer, nchord::Integer, codes; jacobi::Bool=false)
    w = reference_wing(1; split_te=false, nspan=nspan, nchord=nchord)
    out = Dict{Symbol,Any}()
    for code in codes
        r, err = code == :mesa ? mesa_wing(w; nspan=nspan, nchord=nchord) :
                 ours_wing(w; jacobi=jacobi && code == :ours, dofmap=DOFMapOptions(tip_kutta=code == :zero ? :zero : :free))
        isnothing(r) && error("$(CODE_LABEL[code]) failed on the clean $nspan × $nchord wing: $err")
        out[code] = r
        GC.gc()
    end
    return out
end

"""Versions of the codes compared (git commit, +dirty with uncommitted changes), Julia and the thread counts"""
function provenance()
    function sha(dir)
        try
            id = readchomp(`git -C $dir rev-parse --short HEAD`)
            return isempty(readchomp(`git -C $dir status --porcelain --untracked-files=no -- .`)) ? id : id * "+dirty"
        catch
            return "unknown"
        end
    end
    open(joinpath(OUTDIR, "provenance.txt"), "w") do io
        for (name, mod) in (("VortexPanel", VortexPanel), ("MESA", MESA), ("PAGE", PAGE), ("MeshLibrary", MeshLibrary))
            println(io, name, " ", sha(pkgdir(mod)))
        end
        println(io, "Julia ", VERSION, ", threads ", Threads.nthreads(), ", BLAS threads ", BLAS.get_num_threads())
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# 1. refinement

function refinement(levels)
    println("\n1. Refinement on clean meshes: (nspan, nchord) = (17, 19)·2^(level-1), capped tips, α = $(ALPHA_DEG)°. cond₁ of the")
    println("   solved matrix, raw and with the rows scaled by the diagonal; time of assembly and of factorization and solve")
    @printf("%-26s %8s %9s %9s %10s %10s %9s %9s %9s %9s %8s\n", "case", "unknowns", QUANTITIES..., "cond₁", "Jacobi", "assembly", "solve")
    res, rows = Dict{Tuple{Symbol,Int},Any}(), Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        out = clean_wing(nspan, nchord, first.(CODES); jacobi=true)
        for (code, label) in CODES
            r = res[(code, l)] = out[code]
            mesa = code == :mesa
            n, cr = mesa ? (r.npan, r.cond_raw) : (r.nq, r.cond)
            ta, ts = mesa ? (r.timings.assembly, r.timings.lu) : (r.t_assembly, r.t_solve)
            v = [coef(r, q) for q in QUANTITIES]
            @printf("%-26s %8d %9.5f %9.5f %10.6f %10.6f %9.5f %9.1e %9.1e %8.2fs %7.2fs\n", "L$l $label", n, v..., cr, r.cond_jacobi, ta, ts)
            push!(rows, [label l nspan nchord n v... cr r.cond_jacobi ta ts])
        end
    end
    if length(levels) >= 3
        ls = levels[end-2:end]
        println("Richardson extrapolation from levels $(join(ls, ", ")) (p the observed order; NaN when not monotone or p < 0.3)")
        for (code, label) in CODES, q in QUANTITIES
            v = [coef(res[(code, l)], q) for l in ls]
            p, e = extrapolate(v)
            @printf("  %-20s %-6s %s   p = %5.2f, extrapolated %.6f\n", label, q, join([@sprintf("%.6f", x) for x in v], "  "), p, e)
        end
    end
    writedlm(joinpath(OUTDIR, "refinement.csv"), vcat(["code" "level" "nspan" "nchord" "unknowns" "CL" "CL_KJ" "CDi_p" "CDi_T" "Cm_c4" "cond1" "cond1_jacobi" "t_assembly" "t_solve"],
        reduce(vcat, rows)), ',')
    h = [2.0^-(l - 1) for l in levels]
    for (q, qk, title, ylabel) in ((:CL, :CL_KJ, "Lift under refinement, α = 5°", "CL"), (:CDi_p, :CDi_T, "Induced drag under refinement", "CDi"))
        series = Any[]
        for (i, (code, label)) in enumerate(CODES)
            push!(series, (x=h, y=[coef(res[(code, l)], q) for l in levels], label="$label, pressure", color=SVG_COLORS[i]))
            push!(series, (x=h, y=[coef(res[(code, l)], qk) for l in levels], label="$label, " * (qk == :CL_KJ ? "Kutta-Joukowski" : "Trefftz"),
                color=SVG_COLORS[i], dash=true))
        end
        svgplot(joinpath(OUTDIR, "refinement_$q.svg"); series=series, title=title, xlabel="mesh spacing relative to L1", ylabel=ylabel)
    end
    svgplot(joinpath(OUTDIR, "refinement_Cm_c4.svg"); title="Pitching moment about c/4 under refinement", xlabel="mesh spacing relative to L1",
        ylabel="Cm about c/4", hline=0.0, series=[(x=h, y=[res[(code, l)].Cm_c4 for l in levels], label=label) for (code, label) in CODES])
    return res
end

# ---------------------------------------------------------------------------------------------------------------------
# 2. refinement by direction

"""Pressure minus circulation lift coefficient of our solution `r`, in all and from the strips with |η| ≥ η0"""
function lift_gap(r; η0::Real=0.95)
    ld = r.load
    η = (ld.s0 + ld.s1) / WING_SPAN .- 1
    d = (ld.lift_p - ld.lift_kj) .* (ld.s1 - ld.s0) / (WING_SPAN * WING_CHORD)
    return r.CL - r.CL_KJ, sum(d[abs.(η) .>= η0])
end

function directions(res)
    println("\n2. Refinement by direction from the L2 mesh (34 strips, 38 panels per surface)")
    l2 = Dict(code => res[(code, 2)] for (code, _) in CODES)
    study = Dict{Tuple{Symbol,Int},Any}()
    for n in CHORD_STUDY
        study[(:chord, n)] = n == 38 ? l2 : clean_wing(34, n, first.(CODES))
    end
    for n in SPAN_STUDY
        study[(:span, n)] = n == 34 ? l2 : clean_wing(n, 38, first.(CODES))
    end
    rows, ext = Any[], Dict{Tuple{Symbol,Symbol,Symbol},Float64}()
    for (dir, ns, title) in ((:chord, CHORD_STUDY, "chordwise: n panels per surface, 34 strips"), (:span, SPAN_STUDY, "spanwise: n strips, 38 panels per surface"))
        println("  $title (∞: Richardson extrapolation from the last three, NaN when not monotone or the observed order is below 0.3)")
        @printf("  %-20s %6s %9s %9s %10s %10s %9s\n", "", "n", QUANTITIES...)
        for (code, label) in CODES
            for n in ns
                v = [coef(study[(dir, n)][code], q) for q in QUANTITIES]
                @printf("  %-20s %6d %9.5f %9.5f %10.6f %10.6f %9.5f\n", label, n, v...)
                push!(rows, [string(dir) label n v...])
            end
            pe = [extrapolate([coef(study[(dir, n)][code], q) for n in ns[end-2:end]]) for q in QUANTITIES]
            for (q, (_, e)) in zip(QUANTITIES, pe)
                ext[(dir, code, q)] = e
            end
            @printf("  %-20s %6s %9.5f %9.5f %10.6f %10.6f %9.5f\n", "", "∞", last.(pe)...)
            @printf("  %-20s %6s %9.2f %9.2f %10.2f %10.2f %9.2f\n", "", "order", first.(pe)...)
            push!(rows, [string(dir) label Inf last.(pe)...])
        end
    end
    println("  Our pressure minus circulation lift under spanwise refinement, in all and from the outer 5% of the span (|η| ≥ 0.95)")
    for code in (:ours, :zero), n in SPAN_STUDY
        @printf("  %-20s %6d %+9.5f %+9.5f\n", CODE_LABEL[code], n, lift_gap(study[(:span, n)][code])...)
    end
    println("  Estimate at infinite resolution, Q(∞, ∞) ≈ Q(34, ∞) + Q(∞, 38) - Q(34, 38) (the errors of the two directions add;")
    println("  NaN where either direction does not extrapolate)")
    @printf("  %-20s %9s %9s %10s %10s %9s\n", "", QUANTITIES...)
    est = Dict{Tuple{Symbol,Symbol},Float64}()
    for (code, label) in CODES
        v = [ext[(:chord, code, q)] + ext[(:span, code, q)] - coef(l2[code], q) for q in QUANTITIES]
        for (q, x) in zip(QUANTITIES, v)
            est[(code, q)] = x
        end
        @printf("  %-20s %9.5f %9.5f %10.6f %10.6f %9.5f\n", label, v...)
        push!(rows, ["estimate" label Inf v...])
    end
    writedlm(joinpath(OUTDIR, "directions.csv"), vcat(["direction" "code" "n" "CL" "CL_KJ" "CDi_p" "CDi_T" "Cm_c4"], reduce(vcat, rows)), ',')
    for (dir, ns, n0, what) in ((:chord, CHORD_STUDY, 38, "chordwise spacing (38/n)"), (:span, SPAN_STUDY, 34, "spanwise spacing (34/n)"))
        for (q, ylabel) in ((:CL, "CL"), (:Cm_c4, "Cm about c/4"), (:CDi_T, "CDi, Trefftz"))
            series = [(x=[0.0; n0 ./ ns], y=[ext[(dir, code, q)]; [coef(study[(dir, n)][code], q) for n in ns]], label=label, color=SVG_COLORS[i])
                      for (i, (code, label)) in enumerate(CODES)]
            svgplot(joinpath(OUTDIR, "direction_$(dir)_$q.svg"); series=series, xlabel="$what; 0: extrapolated", ylabel=ylabel,
                title="$ylabel under $(dir == :chord ? "chordwise" : "spanwise") refinement from L2")
        end
    end
    return est
end

# ---------------------------------------------------------------------------------------------------------------------
# 3. the two-dimensional limit

"""The structured strips at the middle of the span: the two either side of the middle station (an even number of strips) or the middle one"""
middle_strips(n::Integer) = iseven(n) ? (n ÷ 2, n ÷ 2 + 1) : (n ÷ 2 + 1,)

"""Midspan section of our solution: cl from pressure and from μ_w at the middle, and cm about c/4"""
function ours_midspan(r)
    ld = r.load
    n = length(ld.s0)
    ks = middle_strips(n)
    line = r.dm.wake_te.lines[1]
    nodes = iseven(n) ? (line[n÷2+1],) : (line[n÷2+1], line[n÷2+2])
    y_first = r.dm.wake_te.x[line[1]][2]
    y0, y1 = y_first + ld.s0[ks[1]], y_first + ld.s1[ks[end]]
    M = 0.0
    for (el, c) in zip(r.els, r.cp)
        y0 <= el.centroid[2] <= y1 || continue
        M += cross(el.centroid - SVector(0.25WING_CHORD, el.centroid[2], 0.0), -c * el.area * el.n)[2]
    end
    return (cl_p=sum(ld.lift_p[j] for j in ks) / length(ks) / WING_CHORD,
        cl_kj=2 * sum(r.μ_w[k] for k in nodes) / length(nodes) / norm(r.V) / WING_CHORD, cm=M / (WING_CHORD^2 * (y1 - y0)))
end

"""Midspan section of MESA's solution"""
function mesa_midspan(r, nspan)
    ks = middle_strips(nspan)
    mean(f) = sum(f(j) for j in ks) / length(ks)
    return (cl_p=mean(j -> r.spanwise[j, 4]), cl_kj=mean(j -> r.spanwise[j, 7]) * 2 / (VINF * WING_CHORD), cm=mean(j -> r.cm_strip[j]))
end

"""Value at 1/AR = 0 of the quadratic in 1/AR through three points"""
extrapolate_ar(ar, v) = (hcat(ones(3), 1 ./ ar, 1 ./ ar .^ 2) \ v)[1]

"""Section of PAGE's chordwise spacing (cosine, clustered at the leading edge only), n panels per surface, for the 2-D method"""
function page_section(n::Integer; t::Real=WING_THICKNESS)
    ξ = [1 - cos(π * k / (2n)) for k in 0:n]
    return [[SVector(ξ[k], -naca00_half_thickness(ξ[k]; t=t)) for k in n+1:-1:1]; [SVector(ξ[k], naca00_half_thickness(ξ[k]; t=t)) for k in 2:n+1]]
end

function two_d_limit(levels)
    ex = panel2d(naca00_section(1600), ALPHA_DEG)
    @printf("\n3. The two-dimensional limit. Exact section (2-D linear-vortex panels, 3200): cl = %.5f, cm_c4 = %.5f\n", ex.cl_gamma, ex.cm)
    nchords = [reference_wing_size(l)[2] for l in levels]
    for n in nchords
        s = panel2d(page_section(n), ALPHA_DEG)
        @printf("   for scale, the 2-D method on PAGE's chordwise spacing, %3d panels per surface: cl %+.2f%%, cm %+.4f\n", n,
            pct(s.cl_gamma, ex.cl_gamma), s.cm - ex.cm)
    end
    println("   Midspan section of wings of 68 strips; G Galerkin, M MESA")
    @printf("%-8s %6s %9s %9s %9s   %9s %9s %9s %8s %8s\n", "nchord", "AR", "G cl_p", "G cl_KJ", "G cm", "M cl_p", "M cl_KJ", "M cm", "G time", "M time")
    rows, data = Any[], Dict{Int,Any}()
    for nchord in nchords
        nspan, pts = 68, Any[]
        for m in AR_FACTORS
            b = WING_SPAN * m
            w = reference_wing(1; split_te=false, nspan=nspan, nchord=nchord, span=b)
            o, erro = ours_wing(w; span=b)
            isnothing(o) && error("Galerkin failed on the AR $(AR0 * m) wing: $erro")
            om = ours_midspan(o)
            to = o.t_assembly + o.t_solve
            o = nothing
            GC.gc()
            mr, errm = mesa_wing(w; nspan=nspan, nchord=nchord, span=b)
            isnothing(mr) && error("MESA failed on the AR $(AR0 * m) wing: $errm")
            mm = mesa_midspan(mr, nspan)
            tm = mr.timings.assembly + mr.timings.lu
            GC.gc()
            @printf("%-8d %6.2f %9.5f %9.5f %9.5f   %9.5f %9.5f %9.5f %7.1fs %7.1fs\n", nchord, AR0 * m, om.cl_p, om.cl_kj, om.cm, mm.cl_p, mm.cl_kj, mm.cm, to, tm)
            push!(pts, (AR=AR0 * m, ours=om, mesa=mm))
            push!(rows, [nchord AR0 * m om.cl_p om.cl_kj om.cm mm.cl_p mm.cl_kj mm.cm])
        end
        ar = [p.AR for p in pts]
        lim(f) = extrapolate_ar(ar, [f(p) for p in pts])
        L = (ours_cl=lim(p -> p.ours.cl_p), ours_kj=lim(p -> p.ours.cl_kj), ours_cm=lim(p -> p.ours.cm),
            mesa_cl=lim(p -> p.mesa.cl_p), mesa_kj=lim(p -> p.mesa.cl_kj), mesa_cm=lim(p -> p.mesa.cm))
        @printf("%-8d %6s %9.5f %9.5f %9.5f   %9.5f %9.5f %9.5f\n", nchord, "∞", L.ours_cl, L.ours_kj, L.ours_cm, L.mesa_cl, L.mesa_kj, L.mesa_cm)
        data[nchord] = (points=pts, limit=L)
        push!(rows, [nchord Inf L.ours_cl L.ours_kj L.ours_cm L.mesa_cl L.mesa_kj L.mesa_cm])
    end
    println("   Sectional error at AR = ∞ against the exact section (order: from successive chordwise levels)")
    @printf("   %-8s %12s %8s %12s %8s   %12s %8s %12s %8s\n", "nchord", "G cl", "order", "G cm", "order", "M cl", "order", "M cm", "order")
    err(n) = (data[n].limit.ours_cl / ex.cl_gamma - 1, data[n].limit.ours_cm - ex.cm, data[n].limit.mesa_cl / ex.cl_gamma - 1, data[n].limit.mesa_cm - ex.cm)
    for (i, n) in enumerate(nchords)
        e = err(n)
        order = i == 1 ? fill(NaN, 4) : [log(abs(a / b)) / log(n / nchords[i-1]) for (a, b) in zip(err(nchords[i-1]), e)]
        @printf("   %-8d %+11.2f%% %8.2f %+12.5f %8.2f   %+11.2f%% %8.2f %+12.5f %8.2f\n", n, 100e[1], order[1], e[2], order[2], 100e[3], order[3], e[4], order[4])
    end
    writedlm(joinpath(OUTDIR, "two_d_limit.csv"), vcat(["nchord" "AR" "ours_cl_p" "ours_cl_kj" "ours_cm" "mesa_cl_p" "mesa_cl_kj" "mesa_cm"],
        reduce(vcat, rows), [0 Inf ex.cl_gamma ex.cl_gamma ex.cm ex.cl_gamma ex.cl_gamma ex.cm]), ',')
    for (q, ylabel, title) in ((:cl, "cl at midspan", "Midspan section lift toward the infinite wing"), (:cm, "cm about c/4 at midspan", "Midspan section moment toward the infinite wing"))
        series = Any[]
        for (i, n) in enumerate(nchords)
            d = data[n]
            x = [0.0; [1 / p.AR for p in d.points]]
            go = q == :cl ? [d.limit.ours_cl; [p.ours.cl_p for p in d.points]] : [d.limit.ours_cm; [p.ours.cm for p in d.points]]
            gm = q == :cl ? [d.limit.mesa_cl; [p.mesa.cl_p for p in d.points]] : [d.limit.mesa_cm; [p.mesa.cm for p in d.points]]
            push!(series, (x=x, y=go, label="Galerkin, $n chordwise", color=SVG_COLORS[i]))
            push!(series, (x=x, y=gm, label="MESA, $n chordwise", color=SVG_COLORS[i], dash=true))
        end
        push!(series, (x=[0.0], y=[q == :cl ? ex.cl_gamma : ex.cm], label="exact 2-D section", color="#1f1e1c", line=false))
        svgplot(joinpath(OUTDIR, "two_d_limit_$q.svg"); series=series, title="$title (1/AR = 0: extrapolated)", xlabel="1/AR", ylabel=ylabel)
    end
    return ex
end

# ---------------------------------------------------------------------------------------------------------------------
# 4. chordwise pressure and spanwise loading

function pressure_and_loading(res, ex, level)
    println("\n4. Chordwise Cp and spanwise loading (L$level)")
    o, z, m = res[(:ours, level)], res[(:zero, level)], res[(:mesa, level)]
    nspan, _ = reference_wing_size(level)
    cl_mid = ours_midspan(o).cl_p
    αe = ALPHA_DEG * cl_mid / ex.cl_gamma # the section's lift is nearly linear in α
    for _ in 1:3
        αe *= cl_mid / panel2d(naca00_section(400), αe).cl_gamma
    end
    sec = panel2d(naca00_section(400), αe)
    @printf("   2-D section at the Galerkin midspan cl %.4f (MESA's %.4f): α = %.3f°\n", cl_mid, mesa_midspan(m, nspan).cl_p, αe)
    up2d = [(p[1], c) for (p, c) in zip(sec.x, sec.cp) if p[2] >= 0]
    low2d = [(p[1], c) for (p, c) in zip(sec.x, sec.cp) if p[2] < 0]
    for η in CP_STATIONS
        up, low = chordwise_cp(o, level, η)
        zu, zl = chordwise_cp(z, level, η)
        mu, ml, mx = interp_stations(m.η, m.cp_up, η), interp_stations(m.η, m.cp_low, η), interp_stations(m.η, m.xc_up, η)
        writedlm(joinpath(OUTDIR, @sprintf("cp_L%d_eta%.2f.csv", level, η)), vcat(["code" "x_c" "cp" "surface"],
            hcat(fill("Galerkin", size(up, 1) + size(low, 1)), vcat(up, low)),
            hcat(fill("Galerkin, zero tip", size(zu, 1) + size(zl, 1)), vcat(zu, zl)),
            hcat(fill("MESA", 2length(mx)), [mx; mx], [mu; ml], [fill(1, length(mx)); fill(-1, length(mx))])), ',')
        series = Any[(x=up[:, 1], y=-up[:, 2], label="Galerkin", color=SVG_COLORS[1], marker=false),
            (x=low[:, 1], y=-low[:, 2], label="", color=SVG_COLORS[1], marker=false),
            (x=mx, y=-mu, label="MESA", color=SVG_COLORS[3], dash=true, marker=false), (x=mx, y=-ml, label="", color=SVG_COLORS[3], dash=true, marker=false)]
        if η >= 0.9
            push!(series, (x=zu[:, 1], y=-zu[:, 2], label="Galerkin, zero tip", color=SVG_COLORS[2], marker=false))
            push!(series, (x=zl[:, 1], y=-zl[:, 2], label="", color=SVG_COLORS[2], marker=false))
        end
        if η == 0.0
            push!(series, (x=first.(up2d), y=-last.(up2d), label=@sprintf("2-D section at cl %.3f", cl_mid), color="#52514e", marker=false))
            push!(series, (x=first.(low2d), y=-last.(low2d), label="", color="#52514e", marker=false))
        end
        svgplot(joinpath(OUTDIR, @sprintf("cp_L%d_eta%.2f.svg", level, η)); series=series, title=@sprintf("Chordwise pressure at η = %.2f, L%d", η, level),
            xlabel="x/c", ylabel="-Cp", xlims=(0.0, 1.0))
    end
    η_of(ld) = (ld.s0 + ld.s1) / WING_SPAN .- 1
    sp = m.spanwise
    writedlm(joinpath(OUTDIR, "spanwise_L$level.csv"), vcat(["code" "eta" "cl_kj" "cl_p"],
        hcat(fill("Galerkin", length(o.load.s0)), η_of(o.load), o.load.lift_kj / WING_CHORD, o.load.lift_p / WING_CHORD),
        hcat(fill("Galerkin, zero tip", length(z.load.s0)), η_of(z.load), z.load.lift_kj / WING_CHORD, z.load.lift_p / WING_CHORD),
        hcat(fill("MESA", size(sp, 1)), sp[:, 1], sp[:, 7] * 2 / (VINF * WING_CHORD), sp[:, 4])), ',')
    svgplot(joinpath(OUTDIR, "spanwise_L$level.svg"); title="Spanwise loading, L$level", xlabel="η = 2y/b", ylabel="sectional cl", ylims=(0.0, 0.5),
        series=[(x=η_of(o.load), y=o.load.lift_kj / WING_CHORD, label="Galerkin, circulation", color=SVG_COLORS[1], marker=false),
            (x=η_of(o.load), y=o.load.lift_p / WING_CHORD, label="Galerkin, pressure", color=SVG_COLORS[1], dash=true, marker=false),
            (x=η_of(z.load), y=z.load.lift_kj / WING_CHORD, label="Galerkin zero tip, circulation", color=SVG_COLORS[2], marker=false),
            (x=sp[:, 1], y=sp[:, 7] * 2 / (VINF * WING_CHORD), label="MESA, circulation", color=SVG_COLORS[3], marker=false),
            (x=sp[:, 1], y=sp[:, 4], label="MESA, pressure", color=SVG_COLORS[3], dash=true, marker=false)])
end

# ---------------------------------------------------------------------------------------------------------------------
# 5-8. skew, trailing edge misalignment, zero incidence, robustness

function skew_sensitivity(res, levels)
    println("\n5. Skew sensitivity: changes from the clean mesh of the same level (lean: signed; jitter: mean |Δ| over $(length(SEEDS)) seeds).")
    println("   G Galerkin (default split; the skewed rows' diagonals may flip), Gc Galerkin with the clean mesh's diagonals,")
    println("   M MESA; ratio |ΔCL| M / G; CDi from pressure")
    @printf("%-22s %5s %6s %9s %9s %9s %7s %10s %10s %10s %10s\n", "pattern", "σ", "level", "ΔCL% G", "ΔCL% Gc", "ΔCL% M", "ratio", "ΔCDi% G",
        "ΔCDi% M", "ΔCm G", "ΔCm M")
    rows, summary = Any[], Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        pat = split_pattern(base)
        o0, m0 = res[(:ours, l)], res[(:mesa, l)]
        for (pattern, label) in PATTERNS, σ in SIGMAS
            jitter = pattern in (:jitter_same, :jitter_independent)
            d = Vector{Float64}[]
            for seed in (jitter ? SEEDS : 1:1)
                mesh = first(skew_wing(base, pattern, σ; nspan=nspan, nchord=nchord, rng=StableRNG(seed)))
                o, _ = ours_wing(mesh)
                c, _ = ours_clean_split(mesh, pat)
                m, _ = mesa_wing(mesh; nspan=nspan, nchord=nchord)
                x = [change(o, o0, :CL), change(c, o0, :CL), change(m, m0, :CL), change(o, o0, :CDi_p), change(m, m0, :CDi_p), delta(o, o0, :Cm_c4),
                    delta(m, m0, :Cm_c4)]
                push!(d, x)
                push!(rows, [string(pattern) σ l seed x...])
            end
            a = jitter ? [sum(abs(x[k]) for x in d) / length(d) for k in eachindex(d[1])] : d[1]
            @printf("%-22s %5.2f %6d %+9.3f %+9.3f %+9.3f %7.1f %+10.3f %+10.3f %+10.5f %+10.5f\n", label, σ, l, a[1], a[2], a[3], abs(a[3] / a[1]), a[4:7]...)
            push!(summary, (l, pattern, σ, a))
        end
    end
    for l in levels
        s = [x[4] for x in summary if x[1] == l]
        mean_abs(k) = sum(abs(a[k]) for a in s) / length(s)
        @printf("   L%d, mean |ΔCL| over the patterns and amplitudes: Galerkin %.3f%%, with the clean split %.3f%%, MESA %.3f%%; largest %.3f%%, %.3f%%, %.3f%%\n",
            l, mean_abs(1), mean_abs(2), mean_abs(3), (maximum(abs(a[k]) for a in s) for k in 1:3)...)
        @printf("       mean |ΔCDi| (pressure): Galerkin %.2f%%, MESA %.2f%%; mean |ΔCm|: Galerkin %.5f, MESA %.5f\n", mean_abs(4), mean_abs(5), mean_abs(6), mean_abs(7))
    end
    writedlm(joinpath(OUTDIR, "skew.csv"), vcat(["pattern" "sigma" "level" "seed" "ours_dCL_pct" "ours_clean_split_dCL_pct" "mesa_dCL_pct" "ours_dCDi_pct" "mesa_dCDi_pct" "ours_dCm" "mesa_dCm"],
        reduce(vcat, rows)), ',')
    for l in levels, (k, which) in ((1, ""), (2, ", clean split"))
        series = Any[]
        for (i, (pattern, label)) in enumerate(PATTERNS)
            a = [x[4] for x in summary if x[1] == l && x[2] == pattern]
            push!(series, (x=[0.0; SIGMAS], y=[0.0; [v[k] for v in a]], label="Galerkin$which, $label", color=SVG_COLORS[i]))
            push!(series, (x=[0.0; SIGMAS], y=[0.0; [v[3] for v in a]], label="MESA, $label", color=SVG_COLORS[i], dash=true))
        end
        svgplot(joinpath(OUTDIR, "skew_L$l$(k == 2 ? "_clean_split" : "").svg"); series=series, hline=0.0, xlabel="skew amplitude σ (fraction of the station spacing)",
            ylabel="ΔCL [%]", title="Lift change from panel skew, L$l$which (jitter: mean |ΔCL| over 5 seeds)")
    end
end

"""
Our clean wing of each level with the diagonals of its trailing edge row flipped (both surfaces): the share of the
triangulation in the skew sensitivity, since the skew patterns move the rows next to the trailing edge
"""
function te_row_triangulation(res, levels)
    println("\n   The triangulation alone: the clean mesh with the diagonals of the trailing edge row flipped on both surfaces")
    rows = Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        pat = split_pattern(base)
        for e in eachindex(pat)
            pat[e] != 0 && structured_chord(e, nspan, nchord)[1] == nchord && (pat[e] = 3 - pat[e])
        end
        r, _ = ours_clean_split(base, pat)
        x = (change(r, res[(:ours, l)], :CL), change(r, res[(:ours, l)], :CDi_p), delta(r, res[(:ours, l)], :Cm_c4))
        @printf("   L%d: ΔCL %+.3f%%, ΔCDi %+.3f%%, ΔCm %+.5f\n", l, x...)
        push!(rows, [l x...])
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "te_row_triangulation.csv"), vcat(["level" "ours_dCL_pct" "ours_dCDi_pct" "ours_dCm"], reduce(vcat, rows)), ',')
end

function te_misalignment(res, levels)
    println("\n6. Trailing edge misalignment: the rows aft of 80% chord offset along the span by ±δ/2 (trailing edge shared) or the")
    println("   upper rows by δ (non-matching trailing edge); ΔCL in %; G, Gc and M as in 5.")
    @printf("%-6s %6s %10s %10s %10s %7s %14s %14s   %s\n", "level", "δ", "G shared", "Gc shared", "M shared", "ratio", "G non-match.", "Gc non-match.", "M non-matching")
    rows = Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        pat = split_pattern(base)
        o0, m0 = res[(:ours, l)], res[(:mesa, l)]
        for δ in OFFSETS
            sh = first(skew_wing(base, :te_misaligned, δ; nspan=nspan, nchord=nchord, shared_te=true))
            nm = first(skew_wing(base, :te_misaligned, δ; nspan=nspan, nchord=nchord, shared_te=false))
            d = (change(first(ours_wing(sh)), o0, :CL), change(first(ours_clean_split(sh, pat)), o0, :CL), change(first(mesa_wing(sh; nspan=nspan, nchord=nchord)), m0, :CL),
                change(first(ours_wing(nm)), o0, :CL), change(first(ours_clean_split(nm, pat)), o0, :CL))
            mn, why = mesa_wing(nm; nspan=nspan, nchord=nchord)
            @printf("%-6d %+6.2f %+10.3f %+10.3f %+10.3f %7.1f %+14.3f %+14.3f   %s\n", l, δ, d[1], d[2], d[3], abs(d[3] / d[1]), d[4], d[5],
                isnothing(mn) ? why : @sprintf("%+.3f", change(mn, m0, :CL)))
            push!(rows, [l δ d... change(mn, m0, :CL) replace(why, "," => ";")])
        end
    end
    writedlm(joinpath(OUTDIR, "te_misalignment.csv"), vcat(["level" "delta" "ours_shared_dCL_pct" "ours_clean_split_shared_dCL_pct" "mesa_shared_dCL_pct" "ours_nonmatching_dCL_pct" "ours_clean_split_nonmatching_dCL_pct" "mesa_nonmatching_dCL_pct" "mesa_nonmatching"],
        reduce(vcat, rows)), ',')
    series = Any[]
    for (i, l) in enumerate(levels)
        sel = [r for r in rows if r[1] == l]
        x = [r[2] for r in sel]
        push!(series, (x=x, y=[r[3] for r in sel], label="Galerkin, shared TE, L$l", color=SVG_COLORS[3i-2]))
        push!(series, (x=x, y=[r[4] for r in sel], label="Galerkin clean split, shared TE, L$l", color=SVG_COLORS[3i-2], dash=true))
        push!(series, (x=x, y=[r[6] for r in sel], label="Galerkin, non-matching TE, L$l", color=SVG_COLORS[3i-1]))
        push!(series, (x=x, y=[r[5] for r in sel], label="MESA, shared TE, L$l", color=SVG_COLORS[3i]))
    end
    svgplot(joinpath(OUTDIR, "te_misalignment.svg"); series=series, hline=0.0, title="Lift change from trailing edge misalignment",
        xlabel="spanwise offset δ of the rows aft of 80% chord (fraction of the station spacing)", ylabel="ΔCL [%]")
end

function zero_incidence(level)
    println("\n7. Zero incidence, symmetric section: the lift must stay near zero (L$level)")
    nspan, nchord = reference_wing_size(level)
    base = reference_wing(level; split_te=false)
    skew(p, σ; kw...) = first(skew_wing(base, p, σ; nspan=nspan, nchord=nchord, kw...))
    cases = ["clean" => base, "TE misaligned 0.3, shared" => skew(:te_misaligned, 0.3; shared_te=true),
        "TE misaligned 0.3, non-matching" => skew(:te_misaligned, 0.3; shared_te=false), "lean opposite 0.3" => skew(:lean_opposite, 0.3),
        "lean opposite 0.3, TE rows" => skew(:lean_opposite, 0.3; te_rows=true), "jitter independent 0.3" => skew(:jitter_independent, 0.3)]
    @printf("%-34s %12s %12s %12s %12s   %s\n", "mesh", "Galerkin CL", "max |μ_w|", "zero tip CL", "max |μ_w|", "MESA CL")
    rows = Any[]
    for (name, mesh) in cases
        o, _ = ours_wing(mesh; α=0.0)
        z, _ = ours_wing(mesh; α=0.0, dofmap=DOFMapOptions(tip_kutta=:zero))
        m, why = mesa_wing(mesh; nspan=nspan, nchord=nchord, alpha=0.0)
        @printf("%-34s %12.2e %12.2e %12.2e %12.2e   %s\n", name, o.CL, maximum(abs, o.μ_w), z.CL, maximum(abs, z.μ_w), isnothing(m) ? why : @sprintf("%.2e", m.CL))
        push!(rows, [name o.CL maximum(abs, o.μ_w) z.CL maximum(abs, z.μ_w) (isnothing(m) ? NaN : m.CL)])
    end
    writedlm(joinpath(OUTDIR, "zero_incidence.csv"), vcat(["mesh" "ours_CL" "ours_max_mu_w" "ours_zero_tip_CL" "ours_zero_tip_max_mu_w" "mesa_CL"], reduce(vcat, rows)), ',')
end

function robustness(res, levels)
    println("\n8. Robustness on abused meshes (MESA's smoke set, capped): ΔCL (%) from the clean mesh, or why a code does not run")
    @printf("%-42s %6s %12s %12s\n", "mesh", "level", "Galerkin", "MESA")
    rows = Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        o0, m0 = res[(:ours, l)], res[(:mesa, l)]
        for (name, mesh) in abused_wings(l)
            o, whyo = ours_wing(mesh)
            m, whym = mesa_wing(mesh; nspan=nspan, nchord=nchord)
            so = isnothing(o) ? whyo : @sprintf("%+.3f", change(o, o0, :CL))
            sm = isnothing(m) ? whym : @sprintf("%+.3f", change(m, m0, :CL))
            @printf("%-42s %6d %12s %12s\n", name, l, so, sm)
            push!(rows, [replace(name, "," => ";") l change(o, o0, :CL) change(m, m0, :CL) replace(whyo, "," => ";") replace(whym, "," => ";")])
        end
    end
    writedlm(joinpath(OUTDIR, "robustness.csv"), vcat(["mesh" "level" "ours_dCL_pct" "mesa_dCL_pct" "ours_status" "mesa_status"], reduce(vcat, rows)), ',')
end

function main(args=ARGS)
    levels = parse_args(args)
    mkpath(OUTDIR)
    provenance()
    println("Galerkin (VortexPanel) against Morino (MESA) on the reference wing; threads $(Threads.nthreads()), BLAS threads $(BLAS.get_num_threads())")
    let w = reference_wing(1; split_te=false) # compile both codes first, so the timings exclude compilation
        ours_wing(w; jacobi=true)
        mesa_wing(w; nspan=17, nchord=19)
    end
    res = refinement(levels)
    2 in levels && directions(res)
    ex = two_d_limit(levels)
    lv = filter(<=(2), levels)
    isempty(lv) && (lv = levels[1:1])
    pressure_and_loading(res, ex, maximum(lv))
    skew_sensitivity(res, lv)
    te_row_triangulation(res, levels)
    te_misalignment(res, lv)
    zero_incidence(maximum(lv))
    robustness(res, lv)
    println("\nWrote ", normpath(OUTDIR))
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
