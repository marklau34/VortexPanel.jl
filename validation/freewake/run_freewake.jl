# Free wake with vortex particles (phase 2: the design doc's first deferred item), with MESA's time stepping. An
# implicit trailing edge panel carries the Kutta condition and is factored into the matrix once; free panel rows shed
# every step keep the doublets they were shed with, move with the flow (forward Euler) and are converted into vortex
# particles beyond `max_rows` rows, by the design doc's conversion, local to each wake triangle (one particle each,
# two per quad). The free rows and the particles couple into the body through their velocity (explicit_rhs). The
# reference wing (NACA 0015, aspect ratio 3.88, α = 5°, capped tips), started impulsively.
#
#   julia -t auto --project=validation validation/freewake/run_freewake.jl [levels=1,2] [mesa=true] [abused=1]
#
# A. Consistency: the trailing edge panel with the other rows of a rigid wake as explicit elements, frozen at the
#    implicit solution, gives that solution back (te_len 0.25, 1 and 4 trailing edge panels).
# B. Impulsive start with a prescribed wake (rows moved with V∞), no conversion (L1): the lift by pressure (unsteady
#    Bernoulli), Kutta–Joukowski and impulse, Kelvin's theorem, steady solves on closed rigid wakes of the same rows,
#    the time step and the trailing edge panel's length.
# C. Conversion with the prescribed wake (L1): max_rows 0, 1, 2, 4, 8 and 16 against no conversion.
# D. Free wake to 5 spans (each level): roll-up, lift histories, the total vorticity's drift, cost; sensitivity to the
#    core, the time step and max_rows (L1); VTK of the wake.
# E. MESA's free wake on the same mesh (L1) with the same time step, trailing edge panel, max_rows and particle core,
#    prescribed and free; and both codes' lift against the trailing edge panel's length.
# F. Abused meshes (MESA's smoke set and all at once): the free wake runs, and its lift change from the clean mesh.
# Writes validation/output/freewake/ (CSV tables, SVG plots, VTK).

using MESA, MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs, Logging
using VortexPanel, VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "..", "morino_baseline", "common_mesa.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))
include(joinpath(@__DIR__, "..", "common", "svgplot.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "freewake")
const ALPHA = 5.0
const S_REF = WING_SPAN * WING_CHORD
const VF, LIFT_DIR = wing_freestream(ALPHA)
const SPANS = 5.0 # wake length of the runs, in spans
const NO_CONVERSION = 100_000

function parse_args(args)
    levels, mesa, abused = [1, 2], true, [1]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "mesa" && (mesa = parse(Bool, val))
        key == "abused" && (abused = val == "none" ? Int[] : parse.(Int, split(val, ",")))
    end
    return levels, mesa, abused
end

csv_label(s) = replace(string(s), "," => ";")
te_size(g) = VortexPanel.Wake.te_panel_size(g)

"""Geometry and DOF map of a mesh"""
function setup(mesh; dofmap=DOFMapOptions())
    g = build_geometry(mesh)
    return g, build_dofmap(g; opts=dofmap)
end

"""Lift coefficients of a run's history by the three routes, and the wake length (spans) of each record"""
function routes(sim)
    h = sim.history
    c(f) = [dot(getproperty(r, f), LIFT_DIR) / S_REF for r in h]
    return (p=c(:F_pressure), kj=c(:F_kj), imp=c(:F_impulse), L=[r.t * norm(VF) / sim.te.span for r in h])
end

"""A free wake run of `spans` wake lengths with the options `kw` (FreeWakeOptions)"""
function run_free(g, dm, dt; spans=SPANS, kw...)
    sim = FreeWakeSim(g, dm, VF; opts=FreeWakeOptions(; dt=dt, kw...))
    nsteps = ceil(Int, spans * sim.te.span / (norm(VF) * dt))
    t = @elapsed run!(sim, nsteps)
    return (sim=sim, cl=routes(sim), time=t, nsteps=nsteps)
end

"""Index of the first record whose wake is at least L spans long"""
at_length(cl, L) = something(findfirst(>=(L - 1e-9), cl.L), length(cl.L))

"""
The steady lift of the rigid wake with the trailing edge panel of length `te_len` and `n` rows of length |V∞| dt behind
it, closed by the starting vortex (`ending = :closed`): what the prescribed wake of record n would give if every row
carried the latest μ_w
"""
function closed_cl(g, dm, te_len, dt, n)
    w = build_wake(g, dm, VF; opts=WakeOptions(rows=[0.0; te_len .+ (0:n) .* (norm(VF) * dt)], ending=:closed))
    els = build_elements(g)
    A, b, _ = assemble(g, dm, els, source_strengths(els, VF), VF; wake=w)
    s = solve_galerkin(dm, A, b; W=wake_rows(dm, w.te))
    cp = pressure_coefficient(gradient_velocity(els, dm, s.μ, VF), VF)
    F, _ = pressure_force(els, cp)
    return dot(F, LIFT_DIR) / S_REF, dot(kutta_joukowski_force(dm, dm.K * s.q, VF), LIFT_DIR) / S_REF
end

# ---------------------------------------------------------------------------------------------------------------------
# A. consistency

function consistency_check(levels)
    println("\nA. Consistency: the trailing edge panel (implicit, factored once) with the other rows of the rigid wake as")
    println("   explicit elements frozen at the implicit solution μ_w*, against that solution")
    @printf("%6s %10s %16s %14s\n", "level", "te_len", "max|Δμ_w|/μ_w", "ΔCL (KJ)")
    rows = Any[]
    for l in levels
        g, dm = setup(reference_wing(l; split_te=false))
        els = build_elements(g)
        σ = source_strengths(els, VF)
        for f in (0.25, 1.0, 4.0)
            wf = build_wake(g, dm, VF; opts=WakeOptions(first_row=f * te_size(g)))
            A, b, _ = assemble(g, dm, els, σ, VF; wake=wf)
            ref = solve_galerkin(dm, A, b; W=wake_rows(dm, wf.te))
            μw = dm.K * ref.q
            te = te_panel(g, dm, VF; length=wf.rows[2])
            At, bt, _ = assemble(g, dm, els, σ, VF; wake=te)
            F = factorize_galerkin(dm, At; W=wake_rows(dm, te.te))
            topo = VortexPanel.FreeWake.Topology(dm)
            st = VortexPanel.FreeWake.FreeWakeState{Float64}(topo.nte)
            for j in 1:length(wf.rows)-1
                push!(st.x, wf.x[j*topo.nte+1:(j+1)*topo.nte])
                push!(st.μ, copy(μw))
            end
            st.μ_line .= μw
            src = VortexPanel.FreeWake.ExplicitSources(VortexPanel.FreeWake.vortex_elements(st, topo), VortexPanel.Kernels.Blob{Float64}[])
            s = solve_galerkin(F, bt + explicit_rhs(els, dm, src))
            d = maximum(abs, dm.K * s.q - μw) / maximum(abs, μw)
            dcl = dot(kutta_joukowski_force(dm, dm.K * s.q - μw, VF), LIFT_DIR) / S_REF
            @printf("%6d %10.4f %16.1e %14.1e\n", l, wf.rows[2], d, dcl)
            push!(rows, [l wf.rows[2] d dcl])
        end
    end
    writedlm(joinpath(OUTDIR, "consistency.csv"), vcat(["level" "te_len" "max_rel_dmu_w" "dCL_KJ"], reduce(vcat, rows)), ',')
end

# ---------------------------------------------------------------------------------------------------------------------
# B. impulsive start, prescribed wake, no conversion

const MARKS = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0)

function impulsive_start()
    println("\nB. Impulsive start, prescribed wake (rows moved with V∞), no conversion, L1; dt = the trailing edge panel size,")
    println("   te_len = 0.3 |V∞| dt. Lift by pressure (unsteady Bernoulli), Kutta–Joukowski and impulse, against the steady")
    println("   solve on the closed rigid wake of the same rows (starting vortex at its end) and the open rigid wake")
    g, dm = setup(reference_wing(1; split_te=false))
    dt = te_size(g)
    base = run_free(g, dm, dt; advection=:freestream, max_rows=NO_CONVERSION)
    cl = base.cl
    open = solve_wing(reference_wing(1; split_te=false))
    @printf("%8s %6s %9s %9s %9s %11s %11s %11s %11s\n", "spans", "step", "CL p", "CL KJ", "CL imp", "closed p", "closed KJ", "KJ - closed",
        "closed-open")
    rows = Any[]
    for L in MARKS
        n = at_length(cl, L)
        cp_, ck = closed_cl(g, dm, base.sim.te.rows[2], dt, n - 1)
        @printf("%8.2f %6d %9.5f %9.5f %9.5f %11.5f %11.5f %+11.1e %+11.1e\n", cl.L[n], n, cl.p[n], cl.kj[n], cl.imp[n], cp_, ck, cl.kj[n] - ck,
            ck - open.CL_KJ)
        push!(rows, [cl.L[n] n cl.p[n] cl.kj[n] cl.imp[n] cp_ ck])
    end
    Ωmax = maximum(norm(r.Ω) for r in base.sim.history)
    @printf("   open rigid wake: CL %.5f, CL (KJ) %.5f; largest |Ω_total| over the run %.1e; %d steps in %.1f s\n", open.CL, open.CL_KJ, Ωmax,
        base.nsteps, base.time)
    writedlm(joinpath(OUTDIR, "impulsive.csv"), vcat(["spans" "step" "CL_p" "CL_KJ" "CL_impulse" "closed_CL_p" "closed_CL_KJ"], reduce(vcat, rows)), ',')
    hist = hcat(cl.L, cl.p, cl.kj, cl.imp)
    writedlm(joinpath(OUTDIR, "impulsive_history.csv"), vcat(["spans" "CL_p" "CL_KJ" "CL_impulse"], hist), ',')
    series = [(x=cl.L[2:end], y=cl.p[2:end], label="pressure (unsteady Bernoulli)", color=SVG_COLORS[1], marker=false),
        (x=cl.L[2:end], y=cl.kj[2:end], label="Kutta-Joukowski", color=SVG_COLORS[2], marker=false),
        (x=cl.L[2:end], y=cl.imp[2:end], label="impulse", color=SVG_COLORS[3], marker=false),
        (x=[r[1] for r in rows], y=[r[7] for r in rows], label="closed rigid wake, KJ (steady)", color="#52514e", line=false)]
    svgplot(joinpath(OUTDIR, "impulsive.svg"); series=series, hline=open.CL_KJ, xlabel="wake length [spans]", ylabel="CL",
        title="Impulsive start, prescribed wake (L1); dashed: open rigid wake, KJ")

    println("\n   Time step and trailing edge panel (prescribed wake, no conversion), CL (KJ) at 1 and 2 spans:")
    @printf("   %-28s %10s %10s\n", "", "1 span", "2 spans")
    trows = Any[]
    for (label, kw) in (("dt = 2 × panel", (dt=2dt,)), ("dt = panel", (dt=dt,)), ("dt = panel / 2", (dt=dt / 2,)),
        ("te_len = 0.1 |V∞| dt", (dt=dt, te_len=0.1dt)), ("te_len = |V∞| dt", (dt=dt, te_len=dt)))
        r = run_free(g, dm, kw.dt; spans=2.0, advection=:freestream, max_rows=NO_CONVERSION, (k => v for (k, v) in pairs(kw) if k != :dt)...)
        c1, c2 = r.cl.kj[at_length(r.cl, 1.0)], r.cl.kj[at_length(r.cl, 2.0)]
        @printf("   %-28s %10.5f %10.5f\n", label, c1, c2)
        push!(trows, [csv_label(label) c1 c2])
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "impulsive_dt.csv"), vcat(["case" "CL_KJ_1span" "CL_KJ_2spans"], reduce(vcat, trows)), ',')
    return base, open
end

# ---------------------------------------------------------------------------------------------------------------------
# C. conversion

function conversion(base)
    println("\nC. Conversion, prescribed wake (L1): max_rows against no conversion; ΔCL at 5 spans and the largest |ΔCL (KJ)| over the run")
    g, dm = base.sim.g, base.sim.dm
    dt = base.sim.opts.dt
    @printf("%9s %12s %12s %16s %10s %8s\n", "max_rows", "ΔCL p (%)", "ΔCL KJ (%)", "max|ΔCL KJ|", "particles", "time")
    rows, series = Any[], Any[(x=base.cl.L[2:end], y=base.cl.kj[2:end], label="no conversion", color="#52514e", marker=false)]
    for (i, m) in enumerate((0, 1, 2, 4, 8, 16))
        r = run_free(g, dm, dt; advection=:freestream, max_rows=m)
        n = min(length(r.cl.kj), length(base.cl.kj))
        dp, dk = 100(r.cl.p[n] / base.cl.p[n] - 1), 100(r.cl.kj[n] / base.cl.kj[n] - 1)
        dmax = maximum(abs, r.cl.kj[1:n] - base.cl.kj[1:n])
        @printf("%9d %+12.3f %+12.3f %16.1e %10d %7.1fs\n", m, dp, dk, dmax, length(r.sim.state.particles), r.time)
        push!(rows, [m dp dk dmax length(r.sim.state.particles) r.time])
        push!(series, (x=r.cl.L[2:end], y=r.cl.kj[2:end], label="max_rows $m", color=SVG_COLORS[mod1(i, 8)], marker=false))
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "conversion.csv"), vcat(["max_rows" "dCL_p_pct" "dCL_KJ_pct" "max_abs_dCL_KJ" "particles" "time"], reduce(vcat, rows)), ',')
    svgplot(joinpath(OUTDIR, "conversion.svg"); series=series, xlabel="wake length [spans]", ylabel="CL (KJ)",
        title="Conversion into particles, prescribed wake (L1)")
end

# ---------------------------------------------------------------------------------------------------------------------
# D. free wake

"""Fit CL(L) = a + b/L² to the records beyond L0 spans: the limit a"""
function extrapolate_length(cl, L0)
    k = findall(>=(L0), cl.L)
    M = hcat(ones(length(k)), 1 ./ cl.L[k] .^ 2)
    return (M \ cl.kj[k])[1], (M \ cl.p[k])[1]
end

function free_wake(levels, base)
    println("\nD. Free wake to $SPANS spans: forward Euler on the free nodes and particles, max_rows = 10, core 1.5 × element size")
    rows, hist_rows, series, res = Any[], Any[], Any[], Dict{Int,Any}()
    for (i, l) in enumerate(levels)
        mesh = reference_wing(l; split_te=false)
        g, dm = setup(mesh)
        dt = te_size(g)
        r = run_free(g, dm, dt)
        res[l] = r
        open = solve_wing(mesh)
        h = r.sim.history
        scale = maximum(maximum(abs, rec.μ_w) for rec in h) * r.sim.te.span
        kj_inf, p_inf = extrapolate_length(r.cl, 2.5)
        @printf("   L%d: %d steps in %.0f s (%.2f s per step at the end), %d particles, max |u| %.2f |V∞|, |Ω_total| %.1e of max μ_w × span\n", l,
            r.nsteps, r.time, h[end].t_rhs + h[end].t_vel, length(r.sim.state.particles), maximum(rec.umax for rec in h) / norm(VF),
            norm(h[end].Ω) / scale)
        @printf("       %8s %9s %9s %9s\n", "spans", "CL p", "CL KJ", "CL imp")
        for L in MARKS[3:end]
            n = at_length(r.cl, L)
            @printf("       %8.2f %9.5f %9.5f %9.5f\n", r.cl.L[n], r.cl.p[n], r.cl.kj[n], r.cl.imp[n])
            push!(hist_rows, [l r.cl.L[n] r.cl.p[n] r.cl.kj[n] r.cl.imp[n]])
        end
        @printf("       extrapolated in 1/L² from 2.5 spans on: CL %.5f, CL (KJ) %.5f; the rigid wake: CL %.5f, CL (KJ) %.5f (%+.2f%%, %+.2f%%)\n", p_inf, kj_inf,
            open.CL, open.CL_KJ, 100(p_inf / open.CL - 1), 100(kj_inf / open.CL_KJ - 1))
        push!(rows, [l r.nsteps r.time length(r.sim.state.particles) norm(h[end].Ω) / scale p_inf kj_inf open.CL open.CL_KJ])
        push!(series, (x=r.cl.L[2:end], y=r.cl.kj[2:end], label="free wake, L$l, KJ", color=SVG_COLORS[2i-1], marker=false))
        push!(series, (x=r.cl.L[2:end], y=r.cl.p[2:end], label="free wake, L$l, pressure", color=SVG_COLORS[2i-1], dash=true, marker=false))
        write_free_wake_vtk(r.sim, joinpath(OUTDIR, "free_L$l"))
        MeshLibrary.write_vtk(r.sim.te, joinpath(OUTDIR, "free_L$(l)_te_panel"))
        GC.gc()
    end
    if 1 in levels
        p = run_free(base.sim.g, base.sim.dm, base.sim.opts.dt; advection=:freestream)
        push!(series, (x=p.cl.L[2:end], y=p.cl.kj[2:end], label="prescribed wake, L1, KJ", color="#52514e", marker=false))
        n = at_length(p.cl, SPANS)
        f = res[1].cl
        @printf("   L1 at %.0f spans, free against prescribed (max_rows 10): CL %+.3f%%, CL (KJ) %+.3f%%\n", SPANS, 100(f.p[n] / p.cl.p[n] - 1),
            100(f.kj[n] / p.cl.kj[n] - 1))
    end
    writedlm(joinpath(OUTDIR, "free.csv"), vcat(["level" "steps" "time" "particles" "omega_drift" "CL_p_extrapolated" "CL_KJ_extrapolated" "rigid_CL" "rigid_CL_KJ"],
        reduce(vcat, rows)), ',')
    writedlm(joinpath(OUTDIR, "free_history.csv"), vcat(["level" "spans" "CL_p" "CL_KJ" "CL_impulse"], reduce(vcat, hist_rows)), ',')
    svgplot(joinpath(OUTDIR, "free.svg"); series=series, xlabel="wake length [spans]", ylabel="CL", title="Free wake: lift against wake length")
    return res
end

function sensitivity(res)
    haskey(res, 1) || return
    println("\n   Sensitivity (L1, free wake), CL at $SPANS spans against the defaults:")
    r0 = res[1]
    g, dm, dt = r0.sim.g, r0.sim.dm, r0.sim.opts.dt
    n0 = length(r0.cl.kj)
    rows = Any[]
    for (label, kw) in (("core 1.0 × size", (core_factor=1.0,)), ("core 2.5 × size", (core_factor=2.5,)), ("max_rows 4", (max_rows=4,)),
        ("max_rows 20", (max_rows=20,)), ("dt × 2", (dtf=2.0,)), ("dt / 2", (dtf=0.5,)))
        dtf = get(kw, :dtf, 1.0)
        r = run_free(g, dm, dtf * dt; (k => v for (k, v) in pairs(kw) if k != :dtf)...)
        dp, dk = 100(r.cl.p[end] / r0.cl.p[n0] - 1), 100(r.cl.kj[end] / r0.cl.kj[n0] - 1)
        @printf("   %-18s CL %+.3f%%, CL (KJ) %+.3f%%, %d particles, %.0f s\n", label, dp, dk, length(r.sim.state.particles), r.time)
        push!(rows, [csv_label(label) dp dk length(r.sim.state.particles) r.time])
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "sensitivity.csv"), vcat(["case" "dCL_p_pct" "dCL_KJ_pct" "particles" "time"], reduce(vcat, rows)), ',')
end

# ---------------------------------------------------------------------------------------------------------------------
# E. MESA

"""The L1 wing as MESA needs it: trailing edge split, tips capped with our caps"""
function mesa_mesh(mesh)
    m = deepcopy(mesh)
    PAGE.split_trailing_edge!(m, "trailing_edge")
    return VortexPanel.Geometry.add_caps(m, build_geometry(m))
end

function mesa_comparison()
    println("\nE. MESA's free wake on the same L1 mesh: dt = the trailing edge panel size / |V∞|, te_len = $(TE_LEN) (MESA's), max_rows 10,")
    println("   one Rosenhead–Moore core for every particle (1.5 × the trailing edge panel size) in both codes")
    mesh = reference_wing(1; split_te=false)
    g, dm = setup(mesh)
    ℓ = te_size(g)
    δ = 1.5ℓ
    mm = mesa_mesh(mesh)
    rows, series = Any[], Any[]
    for (i, prescribed) in enumerate((true, false))
        ours = run_free(g, dm, ℓ; te_len=TE_LEN, max_rows=10, core_factor=1e-9, core_floor=δ, advection=prescribed ? :freestream : :free)
        mcl, mt = with_logger(NullLogger()) do
            mesa_freewake(mm; dt=ℓ / VINF, nsteps=ours.nsteps, maxwakelen=10, te_len=TE_LEN, r_vortex=δ, prescribed=prescribed)
        end
        mode = prescribed ? "prescribed" : "free"
        n = min(length(mcl), length(ours.cl.p))
        @printf("   %-10s at %.1f spans: ours CL %.5f, CL (KJ) %.5f (%.0f s); MESA CL %.5f (%.0f s)\n", mode, ours.cl.L[n], ours.cl.p[n], ours.cl.kj[n], ours.time,
            mcl[n], mt)
        for L in MARKS[3:end]
            k = at_length(ours.cl, L)
            k <= length(mcl) && push!(rows, [mode ours.cl.L[k] ours.cl.p[k] ours.cl.kj[k] mcl[k]])
        end
        push!(series, (x=ours.cl.L[2:end], y=ours.cl.p[2:end], label="ours, $mode, pressure", color=SVG_COLORS[2i-1], marker=false))
        push!(series, (x=ours.cl.L[2:n], y=mcl[2:n], label="MESA, $mode", color=SVG_COLORS[2i-1], dash=true, marker=false))
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "mesa.csv"), vcat(["mode" "spans" "ours_CL_p" "ours_CL_KJ" "mesa_CL"], reduce(vcat, rows)), ',')
    svgplot(joinpath(OUTDIR, "mesa.svg"); series=series, xlabel="wake length [spans]", ylabel="CL", title="Free wake: ours against MESA (L1)")

    println("\n   Lift against the trailing edge panel's length (prescribed wake, 2 spans):")
    @printf("   %10s %12s %12s\n", "te_len", "ours CL", "MESA CL")
    trows = Any[]
    for f in (0.25, 0.5, 1.0, 2.0)
        tl = f * ℓ
        ours = run_free(g, dm, ℓ; spans=2.0, te_len=tl, max_rows=10, core_factor=1e-9, core_floor=δ, advection=:freestream)
        mcl, _ = with_logger(NullLogger()) do
            mesa_freewake(mm; dt=ℓ / VINF, nsteps=ours.nsteps, maxwakelen=10, te_len=tl, r_vortex=δ, prescribed=true)
        end
        @printf("   %10.4f %12.5f %12.5f\n", tl, ours.cl.p[end], mcl[end])
        push!(trows, [tl ours.cl.p[end] ours.cl.kj[end] mcl[end]])
        GC.gc()
    end
    writedlm(joinpath(OUTDIR, "mesa_te_len.csv"), vcat(["te_len" "ours_CL_p" "ours_CL_KJ" "mesa_CL"], reduce(vcat, trows)), ',')
end

# ---------------------------------------------------------------------------------------------------------------------
# F. abused meshes

function abused(levels, res)
    println("\nF. Abused meshes, free wake to $SPANS spans: lift change from the clean mesh (%), or why the run failed")
    @printf("%-42s %6s %10s %10s %8s %8s\n", "mesh", "level", "ΔCL p", "ΔCL KJ", "max|u|", "time")
    rows = Any[]
    for l in levels
        r0 = haskey(res, l) ? res[l] : run_free(setup(reference_wing(l; split_te=false))..., te_size(build_geometry(reference_wing(l; split_te=false))))
        for (name, mesh) in abused_wings(l)
            try
                g, dm = setup(mesh)
                r = run_free(g, dm, r0.sim.opts.dt)
                n = min(length(r.cl.p), length(r0.cl.p))
                dp, dk = 100(r.cl.p[n] / r0.cl.p[n] - 1), 100(r.cl.kj[n] / r0.cl.kj[n] - 1)
                um = maximum(rec.umax for rec in r.sim.history) / norm(VF)
                @printf("%-42s %6d %+10.3f %+10.3f %8.2f %7.0fs\n", name, l, dp, dk, um, r.time)
                push!(rows, [csv_label(name) l dp dk um ""])
            catch e
                msg = sprint(showerror, e)[1:min(end, 80)]
                @printf("%-42s %6d   fails: %s\n", name, l, msg)
                push!(rows, [csv_label(name) l NaN NaN NaN csv_label(msg)])
            end
            GC.gc()
        end
    end
    isempty(rows) || writedlm(joinpath(OUTDIR, "abused.csv"), vcat(["mesh" "level" "dCL_p_pct" "dCL_KJ_pct" "max_u" "status"], reduce(vcat, rows)), ',')
end

function main(args=ARGS)
    levels, mesa, abused_levels = parse_args(args)
    mkpath(OUTDIR)
    println("Free wake with vortex particles on the reference wing; threads $(Threads.nthreads())")
    let (g, dm) = setup(reference_wing(1; split_te=false)) # compile first
        run_free(g, dm, te_size(g); spans=0.2)
    end
    consistency_check(levels)
    base, _ = impulsive_start()
    conversion(base)
    res = free_wake(levels, base)
    sensitivity(res)
    mesa && mesa_comparison()
    abused(abused_levels, res)
    println("\nWrote ", normpath(OUTDIR))
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
