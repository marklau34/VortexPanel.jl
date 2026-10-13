# Validation stage 7: the element comparison (doc §11). Split triangles (the default: flat triangles, linear μ, the
# analytic source kernel) against isoparametric bilinear quads (`MixedMesh`: the input's quadrilaterals kept whole,
# bilinear μ, inner integrals by numerical quadrature; tip caps and input triangles stay flat triangles), on the
# cube-sphere and on the reference wing (NACA 0015, aspect ratio 3.88, α = 5°, capped tips, rigid wake), clean and
# abused. Both elements have the same nodes and DOF map, so the unknowns are the same.
#
#   julia -t auto --project=validation validation/elements/run_elements.jl [cube=3,6,12,24] [levels=1,2,3] [section=19,38]
#
# 1. Cube-sphere in uniform flow, clean and jittered: μ and Cp (the gradient route) against the exact solution, λ, the
#    column defect, cond₁ and the assembly time.
# 2. Wing refinement on clean meshes: CL, CDi and Cm by both force routes, the column defect, λ, the condition number
#    (raw and with the rows scaled by the diagonal), the quadrature points and the times.
# 3. Sectional accuracy: the midspan section of wings of aspect ratio 15.5, 31 and 62 (68 strips), extrapolated in 1/AR
#    to the infinite wing, against the exact section (validation/common/panel2d.jl), for the triangles of the default
#    split (whose diagonals mirror at midspan), the triangles of one diagonal throughout, and the quads; Γ around the
#    midspan; and chordwise Cp at midspan (L2).
# 4. Skew sensitivity, the patterns of stage 6 (L1 and L2): triangles with the default split (the skewed rows'
#    diagonals may flip), triangles with the clean mesh's diagonals, and quads (nothing to flip); then triangles and
#    quads with free tips, which shows whether the zero tip's extra sensitivity (stage 6) came from the flips.
# 5. Trailing edge misalignment, shared and non-matching.
# 6. Zero incidence on clean and skewed meshes (L2).
# 7. Abused meshes (MESA's smoke set and all at once, L1 and L2): the lift change and the column defect.
# Writes validation/output/elements/ (CSV tables, SVG plots).

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs, MeshLibrary, PAGE
using VortexPanel, VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))
include(joinpath(@__DIR__, "..", "common", "panel2d.jl"))
include(joinpath(@__DIR__, "..", "common", "exact_flows.jl"))
include(joinpath(@__DIR__, "..", "common", "svgplot.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "elements")
const ALPHA = 5.0
const ELEMENTS = [(:triangle, "triangles"), (:quad, "quads")]
const QUANTITIES = (:CL, :CL_KJ, :CDi_p, :CDi_T, :Cm_c4)
const AR_FACTORS = (4, 8, 16) # spans of the sectional study, in reference spans
const SPHERE = EllipsoidFlow((1.0, 1.0, 1.0), normalize(SVector(1.0, 0.3, -0.2))) # the generic V∞ of stage 3
const ONSPHERE = p -> normalize(p) # moved nodes stay on the unit sphere
const PATTERNS = [(:lean_same, "lean, same"), (:lean_opposite, "lean, opposite"), (:jitter_same, "jitter, same"),
    (:jitter_independent, "jitter, independent")]
const SIGMAS = [0.05, 0.1, 0.2, 0.3]
const SEEDS = 1:5
const OFFSETS = [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]

function parse_args(args)
    cubes, levels, section = [3, 6, 12, 24], [1, 2, 3], [19, 38]
    for a in args
        key, val = split(a, "=")
        key == "cube" && (cubes = parse.(Int, split(val, ",")))
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "section" && (section = parse.(Int, split(val, ",")))
    end
    return cubes, levels, section
end

pct(a, b) = 100 * (a / b - 1)
csv_label(s) = replace(string(s), "," => ";") # writedlm does not quote

"""The wing `mesh` solved with `element` (:triangle or :quad)"""
solve(mesh, element::Symbol; kw...) = solve_wing(mesh; element=element, kw...)

"""
The diagonals the default split chose for the quads of `mesh`: with `GeometryOptions(split=:pattern, split_pattern=...)`
a perturbed copy of the mesh keeps the triangulation (as in stage 6)
"""
function split_pattern(mesh)
    g = build_geometry(mesh)
    p = zeros(Int8, MeshLibrary.nelements(mesh))
    p[g.quad_parent] .= g.quad_diag
    return p
end

"""Triangles with the triangulation `pattern` (`split_pattern` of the clean mesh)"""
solve_clean_split(mesh, pattern; kw...) = solve_wing(mesh; geometry=GeometryOptions(split=:pattern, split_pattern=pattern), kw...)

# ---------------------------------------------------------------------------------------------------------------------
# 1. cube-sphere

"""
The sphere meshed by `mesh` solved with `element` against the exact flow: μ errors at the slots (relative to the
largest exact |μ|, without the mean, which is the null space); Cp from the element gradient against the exact Cp at
the element centres (the triangles' centroids, the quads' parameter centres: `cpL2`, `cpmax`) and over the elements
(`cpL2_surface`: the triangles' constant Cp at 7 points of each, the quads' Cp at their 3×3 Gauss points), all
area-weighted; λ, the column defect, cond₁ and the assembly time. A flat element's normal is constant while the
sphere's turns across it, so away from the centre the error of the tangential velocity is first order.
"""
function solve_sphere(mesh, element::Symbol)
    V = SPHERE.V∞
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    if element == :quad
        els = MixedMesh(g, dm)
        A, b, stats = assemble(g, dm, els, V)
    else
        els = build_elements(g)
        A, b, stats = assemble(g, dm, els, source_strengths(els, V), V)
    end
    sol = solve_galerkin(dm, A, b)
    err = sol.μ - [exact_doublet(SPHERE, g.x[dm.slot_node[s]]) for s in 1:dm.nslot]
    err .-= sum(dm.slot_area .* err) / sum(dm.slot_area)
    scale = doublet_scale(SPHERE)
    rms(e, a) = sqrt(sum(a .* e .^ 2) / sum(a))
    if element == :quad
        c, s = surface_pressure(els, sol.μ, V; m=1), surface_pressure(els, sol.μ, V)
        ec = c.cp - [exact_cp(SPHERE, p) for p in c.x]
        cpL2, cpmax = rms(ec, norm.(c.dA)), maximum(abs, ec)
        cpL2_surface = rms(s.cp - [exact_cp(SPHERE, p) for p in s.x], norm.(s.dA))
    else
        cp = pressure_coefficient(gradient_velocity(els, dm, sol.μ, V), V)
        a = [el.area for el in els]
        ec = cp - [exact_cp(SPHERE, el.centroid) for el in els]
        cpL2, cpmax = rms(ec, a), maximum(abs, ec)
        rule = VortexPanel.Elements.triangle_rule(7)
        cpL2_surface = sqrt(sum(el.area * w * (cp[t] - exact_cp(SPHERE, VortexPanel.Elements.point(el, λ)))^2 for (t, el) in enumerate(els)
                                for (λ, w) in zip(rule.λ, rule.w)) / sum(a))
    end
    return (nodes=dm.nslot, elements=element == :quad ? length(els.quads) + length(els.tris) : length(els),
        μL2=sqrt(sum(dm.slot_area .* err .^ 2) / sum(dm.slot_area)) / scale, μmax=maximum(abs, err) / scale,
        cpL2=cpL2, cpmax=cpmax, cpL2_surface=cpL2_surface, λ=sol.λ[1], columns=consistency(dm, A, b).columns, cond=sol.cond,
        t_assembly=stats.time)
end

function sphere(cubes)
    println("\n1. Cube-sphere in uniform flow (exact μ = V∞·x/2), V∞ = $(round.(SPHERE.V∞, digits=4)): μ relative to its largest")
    println("   exact value, Cp (gradient route) absolute at the element centres and (L2) over the elements; n quads per cube")
    println("   edge; jitter 0.3 with seed 1")
    @printf("%-11s %-10s %4s %7s %9s %9s %9s %9s %9s %9s %9s %8s %7s\n", "mesh", "element", "n", "nodes", "μ L2", "μ max", "Cp L2", "Cp max",
        "Cp L2 el.", "λ", "col. def.", "cond₁", "time")
    families = [("clean", n -> cube_sphere(n)), ("jitter 0.3", n -> first(abuse(cube_sphere(n), (jitter=0.3, project=ONSPHERE); seed=1)))]
    rows, res = Any[], Dict{Tuple{String,Symbol,Int},Any}()
    for (family, make) in families, n in cubes
        mesh = make(n)
        for (element, label) in ELEMENTS
            r = solve_sphere(mesh, element)
            res[(family, element, n)] = r
            @printf("%-11s %-10s %4d %7d %9.2e %9.2e %9.2e %9.2e %9.2e %9.1e %9.1e %8.1e %6.1fs\n", family, label, n, r.nodes, r.μL2, r.μmax, r.cpL2,
                r.cpmax, r.cpL2_surface, r.λ, r.columns, r.cond, r.t_assembly)
            push!(rows, [family label n r.nodes r.elements r.μL2 r.μmax r.cpL2 r.cpmax r.cpL2_surface r.λ r.columns r.cond r.t_assembly])
            GC.gc()
        end
    end
    if length(cubes) >= 2
        n1, n2 = cubes[end-1], cubes[end]
        println("   Observed orders between n = $n1 and $n2, and the quads' error over the triangles' at n = $n2:")
        for (family, _) in families, (element, label) in ELEMENTS
            a, b = res[(family, element, n1)], res[(family, element, n2)]
            o(q) = log(getproperty(a, q) / getproperty(b, q)) / log(n2 / n1)
            @printf("   %-11s %-10s μ L2 %5.2f, μ max %5.2f, Cp L2 %5.2f, Cp L2 over the elements %5.2f", family, label, o(:μL2), o(:μmax), o(:cpL2),
                o(:cpL2_surface))
            if element == :quad
                t = res[(family, :triangle, n2)]
                @printf("   ratio to triangles: μ L2 %.2f, μ max %.2f, Cp L2 %.2f, over the elements %.2f", b.μL2 / t.μL2, b.μmax / t.μmax,
                    b.cpL2 / t.cpL2, b.cpL2_surface / t.cpL2_surface)
            end
            println()
        end
    end
    writedlm(joinpath(OUTDIR, "sphere.csv"), vcat(["mesh" "element" "n" "nodes" "elements" "mu_L2" "mu_max" "cp_L2" "cp_max" "cp_L2_elements" "lambda" "column_defect" "cond1" "t_assembly"],
        reduce(vcat, rows)), ',')
    for (q, name) in ((:μL2, "μ L2 error"), (:cpL2, "Cp L2 error"))
        series = Any[]
        for (i, (family, _)) in enumerate(families), (j, (element, label)) in enumerate(ELEMENTS)
            push!(series, (x=log10.(1 ./ cubes), y=[log10(getproperty(res[(family, element, n)], q)) for n in cubes], label="$label, $family",
                color=SVG_COLORS[j], dash=i == 2))
        end
        svgplot(joinpath(OUTDIR, "sphere_$(q == :μL2 ? "mu" : "cp").svg"); series=series, xlabel="log10(1 / n)", ylabel="log10($name)",
            title="Cube-sphere: $name by element")
    end
    return res
end

# ---------------------------------------------------------------------------------------------------------------------
# 2. wing refinement

function refinement(levels)
    println("\n2. Wing refinement on clean meshes, (nspan, nchord) = (17, 19)·2^(level-1), α = $(ALPHA)°")
    @printf("%-10s %5s %9s %9s %10s %10s %9s %9s %8s %8s %8s %8s %8s %7s\n", "element", "level", QUANTITIES..., "col. def.", "λ", "cond₁", "Jacobi",
        "points", "t asm", "t solve")
    res, rows = Dict{Tuple{Symbol,Int},Any}(), Any[]
    for l in levels, (element, label) in ELEMENTS
        r = solve(reference_wing(l; split_te=false), element; α=ALPHA, jacobi=true)
        res[(element, l)] = r
        v = [getproperty(r, q) for q in QUANTITIES]
        npts = sum(r.stats.points)
        @printf("%-10s %5d %9.5f %9.5f %10.6f %10.6f %+9.5f %9.1e %8.1e %8.1e %8.1e %7.1fM %7.1fs %6.2fs\n", label, l, v..., r.columns, r.λ, r.cond,
            r.cond_jacobi, npts / 1e6, r.t_assembly, r.t_solve)
        push!(rows, [label l r.nq v... r.columns r.λ r.cond r.cond_jacobi npts r.t_assembly r.t_solve])
        GC.gc()
    end
    if length(levels) >= 3
        ls = levels[end-2:end]
        println("Richardson extrapolation from levels $(join(ls, ", ")) (p the observed order; NaN when not monotone or p < 0.3)")
        for (element, label) in ELEMENTS, q in QUANTITIES
            v = [getproperty(res[(element, l)], q) for l in ls]
            p, e = extrapolate(v)
            @printf("  %-10s %-6s %s   p = %5.2f, extrapolated %.6f\n", label, q, join([@sprintf("%.6f", x) for x in v], "  "), p, e)
        end
    end
    for l in levels
        t, q = res[(:triangle, l)], res[(:quad, l)]
        @printf("   L%d, quads against triangles: ΔCL %+.3f%%, ΔCL (KJ) %+.3f%%, ΔCDi %+.2f%% (pressure) %+.2f%% (Trefftz), ΔCm %+.5f;", l,
            pct(q.CL, t.CL), pct(q.CL_KJ, t.CL_KJ), pct(q.CDi_p, t.CDi_p), pct(q.CDi_T, t.CDi_T), q.Cm_c4 - t.Cm_c4)
        @printf(" assembly %.1f× as long, %.1f× the points\n", q.t_assembly / t.t_assembly, sum(q.stats.points) / sum(t.stats.points))
    end
    writedlm(joinpath(OUTDIR, "refinement.csv"), vcat(["element" "level" "unknowns" "CL" "CL_KJ" "CDi_p" "CDi_T" "Cm_c4" "column_defect" "lambda" "cond1" "cond1_jacobi" "points" "t_assembly" "t_solve"],
        reduce(vcat, rows)), ',')
    h = [2.0^-(l - 1) for l in levels]
    for (q, qk, title, ylabel) in ((:CL, :CL_KJ, "Lift under refinement by element", "CL"), (:CDi_p, :CDi_T, "Induced drag under refinement by element", "CDi"))
        series = Any[]
        for (i, (element, label)) in enumerate(ELEMENTS)
            push!(series, (x=h, y=[getproperty(res[(element, l)], q) for l in levels], label="$label, pressure", color=SVG_COLORS[i]))
            push!(series, (x=h, y=[getproperty(res[(element, l)], qk) for l in levels], label="$label, " * (qk == :CL_KJ ? "Kutta-Joukowski" : "Trefftz"),
                color=SVG_COLORS[i], dash=true))
        end
        svgplot(joinpath(OUTDIR, "refinement_$q.svg"); series=series, title=title, xlabel="mesh spacing relative to L1", ylabel=ylabel)
    end
    return res
end

# ---------------------------------------------------------------------------------------------------------------------
# 3. sectional accuracy and chordwise pressure

"""
The variants of the sectional study: the triangles with the default split, whose diagonals mirror at y = 0 (the
midspan, where the section is measured), the triangles with one diagonal throughout (`split=:d13`), and the quads
"""
const SECTION_VARIANTS = [(:triangle, "triangles"), (:one_diagonal, "triangles, one diagonal"), (:quad, "quads")]

section_solve(mesh, variant::Symbol; kw...) =
    variant == :one_diagonal ? solve_wing(mesh; geometry=GeometryOptions(split=:d13), kw...) : solve(mesh, variant; kw...)

function section_accuracy(nchords)
    ex = panel2d(naca00_section(1600), ALPHA)
    @printf("\n3. Sectional accuracy: midspan of wings of 68 strips at AR 15.5, 31 and 62, extrapolated to AR = ∞, against the exact section\n")
    @printf("   (cl %.5f, cm_c4 %.5f); order from successive chordwise levels. The default split mirrors its diagonals at y = 0, the\n", ex.cl_gamma, ex.cm)
    println("   midspan; \"one diagonal\" splits every quad along the same diagonal (GeometryOptions(split=:d13)).")
    @printf("%-24s %7s %10s %8s %10s %8s %10s %8s %9s %8s\n", "element", "nchord", "cl", "order", "cl (KJ)", "order", "cm", "order", "cond₁", "time")
    rows, raw, err = Any[], Any[], Dict{Tuple{Symbol,Int},NTuple{3,Float64}}()
    for (element, label) in SECTION_VARIANTS, (i, nchord) in enumerate(nchords)
        t = @elapsed pts = map(AR_FACTORS) do m
            b = WING_SPAN * m
            r = section_solve(reference_wing(1; split_te=false, nspan=68, nchord=nchord, span=b), element; α=ALPHA, span=b)
            GC.gc()
            merge((AR=b / WING_CHORD, cond=r.cond, columns=r.columns), midspan_section(r))
        end
        for p in pts
            push!(raw, [csv_label(label) nchord p.AR p.cl_p p.cl_kj p.cm p.columns p.cond])
        end
        ar = [p.AR for p in pts]
        cl, kj, cm = (extrapolate_ar(ar, [getproperty(p, q) for p in pts]) for q in (:cl_p, :cl_kj, :cm))
        e = (cl / ex.cl_gamma - 1, kj / ex.cl_gamma - 1, cm - ex.cm)
        err[(element, nchord)] = e
        order = i == 1 ? fill(NaN, 3) : [log(abs(a / b)) / log(nchord / nchords[i-1]) for (a, b) in zip(err[(element, nchords[i-1])], e)]
        cmax = maximum(p.cond for p in pts)
        @printf("%-24s %7d %+9.2f%% %8.2f %+9.2f%% %8.2f %+10.5f %8.2f %9.1e %7.0fs\n", label, nchord, 100e[1], order[1], 100e[2], order[2], e[3], order[3],
            cmax, t)
        push!(rows, [csv_label(label) nchord cl kj cm e... cmax t])
    end
    writedlm(joinpath(OUTDIR, "section.csv"), vcat(["element" "nchord" "cl" "cl_kj" "cm" "cl_error" "cl_kj_error" "cm_error" "max_cond1" "time"], reduce(vcat, rows),
        ["exact" 0 ex.cl_gamma ex.cl_gamma ex.cm 0 0 0 NaN NaN]), ',')
    writedlm(joinpath(OUTDIR, "section_by_aspect_ratio.csv"), vcat(["element" "nchord" "AR" "cl" "cl_kj" "cm" "column_defect" "cond1"], reduce(vcat, raw)), ',')
    for (k, name, ylabel) in ((1, "cl", "|cl error| [%]"), (3, "cm", "|cm error|"))
        series = [(x=1 ./ nchords, y=[abs(err[(element, n)][k]) * (k == 1 ? 100 : 1) for n in nchords], label=label, color=SVG_COLORS[i])
                  for (i, (element, label)) in enumerate(SECTION_VARIANTS)]
        svgplot(joinpath(OUTDIR, "section_$name.svg"); series=series, xlabel="1 / panels per surface", ylabel=ylabel,
            title="Sectional $name error against the exact section, by element")
    end
    return ex, err
end

"""
Γ along the trailing edge near midspan (the wing of AR 62, 68 strips, `nchord` panels per surface) for each section
variant: the default split's mirror line at y = 0 shows as a local dip
"""
function midspan_circulation(nchord)
    b = WING_SPAN * AR_FACTORS[end]
    mesh = reference_wing(1; split_te=false, nspan=68, nchord=nchord, span=b)
    rs = [section_solve(mesh, v; α=ALPHA, span=b) for (v, _) in SECTION_VARIANTS]
    line = rs[1].dm.wake_te.lines[1]
    mid = length(line) ÷ 2 + 1
    println("\n   Γ at the trailing edge nodes around midspan (AR 62, $nchord panels per surface), relative to the quads':")
    @printf("   %-24s %s\n", "y / chord", join([@sprintf("%8.2f", rs[1].dm.wake_te.x[line[j]][2] / WING_CHORD) for j in mid-3:mid+3], ""))
    rows = Any[]
    for (r, (_, label)) in zip(rs, SECTION_VARIANTS)
        d = [100 * (r.μ_w[line[j]] / rs[end].μ_w[line[j]] - 1) for j in mid-3:mid+3]
        @printf("   %-24s %s  (%%)\n", label, join([@sprintf("%+8.3f", x) for x in d], ""))
        push!(rows, [csv_label(label) d...])
    end
    writedlm(joinpath(OUTDIR, "midspan_circulation.csv"), vcat(["variant" ["dGamma_pct_node_$k" for k in -3:3]...], reduce(vcat, rows)), ',')
end

function midspan_pressure(res, ex, level)
    println("\n   Chordwise Cp at midspan (L$level): suction peak, and the largest difference from the triangles between 10% and 90% chord")
    interp(x, y, xi) = (k = clamp(searchsortedlast(x, xi), 1, length(x) - 1); w = (xi - x[k]) / (x[k+1] - x[k]); (1 - w) * y[k] + w * y[k+1])
    tu, tl = chordwise_cp(res[(:triangle, level)], level, 0.0)
    rows, series = Any[], Any[]
    for (i, (element, label)) in enumerate(ELEMENTS)
        up, low = chordwise_cp(res[(element, level)], level, 0.0)
        du = maximum(abs(interp(up[:, 1], up[:, 2], x) - interp(tu[:, 1], tu[:, 2], x)) for x in 0.1:0.02:0.9)
        dl = maximum(abs(interp(low[:, 1], low[:, 2], x) - interp(tl[:, 1], tl[:, 2], x)) for x in 0.1:0.02:0.9)
        @printf("   %-10s suction peak %.3f at x/c %.3f; |ΔCp| from triangles at most %.3f (upper), %.3f (lower); trailing edge panels %.3f / %.3f\n",
            label, minimum(up[:, 2]), up[argmin(up[:, 2]), 1], du, dl, up[end, 2], low[end, 2])
        push!(rows, hcat(fill(label, size(up, 1) + size(low, 1)), vcat(up, low)))
        push!(series, (x=up[:, 1], y=-up[:, 2], label=label, color=SVG_COLORS[i], marker=false))
        push!(series, (x=low[:, 1], y=-low[:, 2], label="", color=SVG_COLORS[i], marker=false))
    end
    cl_mid = midspan_section(res[(:triangle, level)]).cl_p
    αe = ALPHA * cl_mid / ex.cl_gamma
    for _ in 1:3
        αe *= cl_mid / panel2d(naca00_section(400), αe).cl_gamma
    end
    sec = panel2d(naca00_section(400), αe)
    for upper in (true, false)
        pts = [(p[1], c) for (p, c) in zip(sec.x, sec.cp) if (p[2] >= 0) == upper]
        push!(series, (x=first.(pts), y=-last.(pts), label=upper ? @sprintf("2-D section at cl %.3f", cl_mid) : "", color="#52514e", dash=true, marker=false))
    end
    writedlm(joinpath(OUTDIR, "cp_L$(level)_midspan.csv"), vcat(["element" "x_c" "cp" "surface"], reduce(vcat, rows)), ',')
    svgplot(joinpath(OUTDIR, "cp_L$(level)_midspan.svg"); series=series, title="Chordwise pressure at midspan by element, L$level",
        xlabel="x/c", ylabel="-Cp", xlims=(0.0, 1.0))
end

# ---------------------------------------------------------------------------------------------------------------------
# 4. skew, 5. trailing edge misalignment

change(r, r0, q::Symbol) = pct(getproperty(r, q), getproperty(r0, q))

function skew_sensitivity(res, levels)
    println("\n4. Skew sensitivity: changes from the clean mesh of the same level and element (lean: signed; jitter: mean |Δ| over")
    println("   $(length(SEEDS)) seeds). T triangles (default split; the skewed rows' diagonals may flip), Tc triangles with the clean")
    println("   mesh's diagonals, Q quads; CDi from pressure")
    @printf("%-22s %5s %6s %9s %9s %9s %10s %10s %10s %10s\n", "pattern", "σ", "level", "ΔCL% T", "ΔCL% Tc", "ΔCL% Q", "ΔCDi% T", "ΔCDi% Q", "ΔCm T", "ΔCm Q")
    rows, summary = Any[], Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        pat = split_pattern(base)
        t0, q0 = res[(:triangle, l)], res[(:quad, l)]
        for (pattern, label) in PATTERNS, σ in SIGMAS
            jitter = pattern in (:jitter_same, :jitter_independent)
            d = Vector{Float64}[]
            for seed in (jitter ? SEEDS : 1:1)
                mesh = first(skew_wing(base, pattern, σ; nspan=nspan, nchord=nchord, rng=StableRNG(seed)))
                t = solve(mesh, :triangle; α=ALPHA)
                c = solve_clean_split(mesh, pat; α=ALPHA)
                q = solve(mesh, :quad; α=ALPHA)
                x = [change(t, t0, :CL), change(c, t0, :CL), change(q, q0, :CL), change(t, t0, :CDi_p), change(q, q0, :CDi_p), t.Cm_c4 - t0.Cm_c4,
                    q.Cm_c4 - q0.Cm_c4]
                push!(d, x)
                push!(rows, [string(pattern) σ l seed x...])
                GC.gc()
            end
            a = jitter ? [sum(abs(x[k]) for x in d) / length(d) for k in eachindex(d[1])] : d[1]
            @printf("%-22s %5.2f %6d %+9.3f %+9.3f %+9.3f %+10.3f %+10.3f %+10.5f %+10.5f\n", label, σ, l, a...)
            push!(summary, (l, pattern, σ, a))
        end
    end
    for l in levels
        s = [x[4] for x in summary if x[1] == l]
        mean_abs(k) = sum(abs(a[k]) for a in s) / length(s)
        @printf("   L%d, mean |ΔCL| over the patterns and amplitudes: triangles %.3f%%, with the clean split %.3f%%, quads %.3f%%; largest %.3f%%, %.3f%%, %.3f%%\n",
            l, mean_abs(1), mean_abs(2), mean_abs(3), (maximum(abs(a[k]) for a in s) for k in 1:3)...)
        @printf("       mean |ΔCDi| (pressure): triangles %.2f%%, quads %.2f%%; mean |ΔCm|: triangles %.5f, quads %.5f\n", mean_abs(4), mean_abs(5), mean_abs(6), mean_abs(7))
    end
    writedlm(joinpath(OUTDIR, "skew.csv"), vcat(["pattern" "sigma" "level" "seed" "tri_dCL_pct" "tri_clean_split_dCL_pct" "quad_dCL_pct" "tri_dCDi_pct" "quad_dCDi_pct" "tri_dCm" "quad_dCm"],
        reduce(vcat, rows)), ',')
    for l in levels
        series = Any[]
        for (i, (pattern, label)) in enumerate(PATTERNS)
            a = [x[4] for x in summary if x[1] == l && x[2] == pattern]
            push!(series, (x=[0.0; SIGMAS], y=[0.0; [v[1] for v in a]], label="triangles, $label", color=SVG_COLORS[i], dash=true))
            push!(series, (x=[0.0; SIGMAS], y=[0.0; [v[3] for v in a]], label="quads, $label", color=SVG_COLORS[i]))
        end
        svgplot(joinpath(OUTDIR, "skew_L$l.svg"); series=series, hline=0.0, xlabel="skew amplitude σ (fraction of the station spacing)",
            ylabel="ΔCL [%]", title="Lift change from panel skew by element, L$l (jitter: mean |ΔCL| over 5 seeds)")
    end
    return summary
end

"""
Skew sensitivity with free tips (`DOFMapOptions(tip_kutta=:free)`), the same patterns as `skew_sensitivity`. Stage 6
found that the zero tip doubles the lift cost of the split's diagonal flips next to the trailing edge; the quads have
nothing to flip, so their sensitivity with either tip shows whether that extra cost came from the flips. With
`zero` (the summary `skew_sensitivity` returns) the zero-tip figures are printed alongside.
"""
function free_tip_skew(levels; zero=nothing)
    println("\n4b. Skew sensitivity with free tips: changes from the clean mesh of the same level, element and tip condition")
    @printf("%-22s %5s %6s %12s %12s\n", "pattern", "σ", "level", "ΔCL% T free", "ΔCL% Q free")
    free = DOFMapOptions(tip_kutta=:free)
    rows, summary = Any[], Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        t0, q0 = solve(base, :triangle; α=ALPHA, dofmap=free), solve(base, :quad; α=ALPHA, dofmap=free)
        for (pattern, label) in PATTERNS, σ in SIGMAS
            jitter = pattern in (:jitter_same, :jitter_independent)
            d = Vector{Float64}[]
            for seed in (jitter ? SEEDS : 1:1)
                mesh = first(skew_wing(base, pattern, σ; nspan=nspan, nchord=nchord, rng=StableRNG(seed)))
                x = [change(solve(mesh, :triangle; α=ALPHA, dofmap=free), t0, :CL), change(solve(mesh, :quad; α=ALPHA, dofmap=free), q0, :CL)]
                push!(d, x)
                push!(rows, [string(pattern) σ l seed x...])
                GC.gc()
            end
            a = jitter ? [sum(abs(x[k]) for x in d) / length(d) for k in 1:2] : d[1]
            @printf("%-22s %5.2f %6d %+12.3f %+12.3f\n", label, σ, l, a...)
            push!(summary, (l, pattern, σ, a))
        end
    end
    println("   Mean (largest) |ΔCL| over the patterns and amplitudes, by element and tip condition:")
    for l in levels
        s = [x[4] for x in summary if x[1] == l]
        stat(v) = @sprintf("%.3f%% (%.3f%%)", sum(abs, v) / length(v), maximum(abs, v))
        line = "   L$l: triangles, free tips $(stat([a[1] for a in s])); quads, free tips $(stat([a[2] for a in s]))"
        if !isnothing(zero)
            z = [x[4] for x in zero if x[1] == l]
            line *= "; zero tips: triangles $(stat([a[1] for a in z])), quads $(stat([a[3] for a in z]))"
        end
        println(line)
    end
    writedlm(joinpath(OUTDIR, "skew_free_tip.csv"), vcat(["pattern" "sigma" "level" "seed" "tri_free_tip_dCL_pct" "quad_free_tip_dCL_pct"], reduce(vcat, rows)), ',')
    return summary
end

function te_misalignment(res, levels)
    println("\n5. Trailing edge misalignment: the rows aft of 80% chord offset along the span by ±δ/2 (trailing edge shared) or the")
    println("   upper rows by δ (non-matching trailing edge); ΔCL in %; T, Tc and Q as in 4.")
    @printf("%-6s %6s %10s %10s %10s %14s %14s %14s\n", "level", "δ", "T shared", "Tc shared", "Q shared", "T non-match.", "Tc non-match.", "Q non-match.")
    rows = Any[]
    for l in levels
        nspan, nchord = reference_wing_size(l)
        base = reference_wing(l; split_te=false)
        pat = split_pattern(base)
        t0, q0 = res[(:triangle, l)], res[(:quad, l)]
        for δ in OFFSETS
            d = Float64[]
            for shared in (true, false)
                mesh = first(skew_wing(base, :te_misaligned, δ; nspan=nspan, nchord=nchord, shared_te=shared))
                push!(d, change(solve(mesh, :triangle; α=ALPHA), t0, :CL), change(solve_clean_split(mesh, pat; α=ALPHA), t0, :CL),
                    change(solve(mesh, :quad; α=ALPHA), q0, :CL))
                GC.gc()
            end
            @printf("%-6d %+6.2f %+10.3f %+10.3f %+10.3f %+14.3f %+14.3f %+14.3f\n", l, δ, d...)
            push!(rows, [l δ d...])
        end
    end
    writedlm(joinpath(OUTDIR, "te_misalignment.csv"), vcat(["level" "delta" "tri_shared_dCL_pct" "tri_clean_split_shared_dCL_pct" "quad_shared_dCL_pct" "tri_nonmatching_dCL_pct" "tri_clean_split_nonmatching_dCL_pct" "quad_nonmatching_dCL_pct"],
        reduce(vcat, rows)), ',')
    series = Any[]
    for (i, l) in enumerate(levels)
        sel = [r for r in rows if r[1] == l]
        x = [r[2] for r in sel]
        push!(series, (x=x, y=[r[3] for r in sel], label="triangles, shared TE, L$l", color=SVG_COLORS[2i-1], dash=true))
        push!(series, (x=x, y=[r[5] for r in sel], label="quads, shared TE, L$l", color=SVG_COLORS[2i-1]))
        push!(series, (x=x, y=[r[6] for r in sel], label="triangles, non-matching TE, L$l", color=SVG_COLORS[2i], dash=true))
        push!(series, (x=x, y=[r[8] for r in sel], label="quads, non-matching TE, L$l", color=SVG_COLORS[2i]))
    end
    svgplot(joinpath(OUTDIR, "te_misalignment.svg"); series=series, hline=0.0, title="Lift change from trailing edge misalignment by element",
        xlabel="spanwise offset δ of the rows aft of 80% chord (fraction of the station spacing)", ylabel="ΔCL [%]")
end

# ---------------------------------------------------------------------------------------------------------------------
# 6. zero incidence, 7. abused meshes

function zero_incidence(level)
    println("\n6. Zero incidence, symmetric section: the lift must stay near zero (L$level)")
    nspan, nchord = reference_wing_size(level)
    base = reference_wing(level; split_te=false)
    skew(p, σ; kw...) = first(skew_wing(base, p, σ; nspan=nspan, nchord=nchord, kw...))
    cases = ["clean" => base, "TE misaligned 0.3, shared" => skew(:te_misaligned, 0.3; shared_te=true),
        "TE misaligned 0.3, non-matching" => skew(:te_misaligned, 0.3; shared_te=false), "lean opposite 0.3" => skew(:lean_opposite, 0.3),
        "lean opposite 0.3, TE rows" => skew(:lean_opposite, 0.3; te_rows=true), "jitter independent 0.3" => skew(:jitter_independent, 0.3)]
    @printf("%-34s %12s %12s\n", "mesh", "triangles CL", "quads CL")
    rows = Any[]
    for (name, mesh) in cases
        cl = [solve(mesh, element; α=0.0).CL for (element, _) in ELEMENTS]
        @printf("%-34s %12.2e %12.2e\n", name, cl...)
        push!(rows, [csv_label(name) cl...])
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "zero_incidence.csv"), vcat(["mesh" "tri_CL" "quad_CL"], reduce(vcat, rows)), ',')
end

function abused(res, levels)
    println("\n7. Abused meshes (MESA's smoke set and all at once): lift change (%) from the clean mesh of the same element, and the column defect")
    @printf("%-42s %6s %12s %12s %12s %12s %10s %10s\n", "mesh", "level", "T ΔCL", "Q ΔCL", "T ΔCL (KJ)", "Q ΔCL (KJ)", "T col.", "Q col.")
    rows, summary = Any[], Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for l in levels
        for (name, mesh) in abused_wings(l)
            d = map(ELEMENTS) do (element, _)
                r = solve(mesh, element; α=ALPHA)
                r0 = res[(element, l)]
                push!(get!(summary, (element, l), Float64[]), change(r, r0, :CL))
                (change(r, r0, :CL), change(r, r0, :CL_KJ), r.columns)
            end
            @printf("%-42s %6d %+12.3f %+12.3f %+12.3f %+12.3f %10.1e %10.1e\n", name, l, d[1][1], d[2][1], d[1][2], d[2][2], d[1][3], d[2][3])
            push!(rows, [csv_label(name) l d[1]... d[2]...])
            GC.gc()
        end
    end
    for l in levels
        @printf("   L%d, mean (largest) |ΔCL|: %s\n", l, join([@sprintf("%s %.3f%% (%.3f%%)", label, sum(abs, summary[(element, l)]) / length(summary[(element, l)]),
            maximum(abs, summary[(element, l)])) for (element, label) in ELEMENTS], ", "))
    end
    writedlm(joinpath(OUTDIR, "abused.csv"), vcat(["mesh" "level" "tri_dCL_pct" "tri_dCL_KJ_pct" "tri_column_defect" "quad_dCL_pct" "quad_dCL_KJ_pct" "quad_column_defect"],
        reduce(vcat, rows)), ',')
    for l in levels
        n = length(summary[(:triangle, l)])
        series = [(x=collect(1.0:n), y=summary[(element, l)], label=label, color=SVG_COLORS[i], line=false) for (i, (element, label)) in enumerate(ELEMENTS)]
        svgplot(joinpath(OUTDIR, "abused_L$l.svg"); series=series, hline=0.0, xlabel="abused mesh (the order of abused.csv)", ylabel="ΔCL [%]",
            title="Lift change on abused meshes by element, L$l")
    end
end

function main(args=ARGS)
    cubes, levels, section = parse_args(args)
    mkpath(OUTDIR)
    println("Split triangles against isoparametric quads; threads $(Threads.nthreads()), BLAS threads $(BLAS.get_num_threads())")
    for (element, _) in ELEMENTS # compile first, so the timings exclude it
        solve(reference_wing(1; split_te=false), element; jacobi=true)
        solve_sphere(cube_sphere(2), element)
    end
    sphere(cubes)
    res = refinement(levels)
    ex, _ = section_accuracy(section)
    midspan_circulation(first(section))
    lv = filter(<=(2), levels)
    isempty(lv) && (lv = levels[1:1])
    midspan_pressure(res, ex, maximum(lv))
    sk = skew_sensitivity(res, lv)
    free_tip_skew(lv; zero=sk)
    te_misalignment(res, lv)
    zero_incidence(maximum(lv))
    abused(res, lv)
    println("\nWrote ", normpath(OUTDIR))
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
