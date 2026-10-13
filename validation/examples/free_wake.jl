# Example: the free wake with vortex particles on the reference wing (NACA 0015, aspect ratio 3.88, capped tips),
# started impulsively, as a time series for ParaView: the body coloured by pressure, the trailing edge panel, the free
# panel rows and the particles, frame by frame as the wake rolls up.
#
#   julia -t auto --project=validation validation/examples/free_wake.jl [level=1] [spans=5] [max_rows=auto] [merge_rows=1]
#                                                                       [alpha=5] [every=4] [advection=free] [dt=1]
#
# level: the reference wing's refinement (1: 1,364 triangles, about 30 s; 2: about 8 min). spans: how far the starting
# vortex travels, in spans. dt: the time step in trailing edge panel lengths travelled (1: one step sheds a row as long
# as a trailing edge panel). max_rows: free panel rows kept before the oldest is converted into particles (auto: as many
# as put the new particles three cores behind the trailing edge, whatever dt; 0: every new row at once). merge_rows:
# converted rows merged into one particle pair per strip (auto: about as many as make it as long as the strip is wide,
# which keeps the particle count down at a small dt). every: write a frame every this many steps. advection: free
# (forward Euler with the local velocity) or freestream (a prescribed wake).
#
# Writes to validation/output/examples/free_wake_L<level>_rows<max_rows>[_merge<merge_rows>][_dt<dt>]/:
#   body.pvd        the body: colour by "Cp" (on each triangle, with the unsteady Bernoulli term) or "Cp_nodes"
#                   (smooth); also "mu", the doublet strength at the nodes
#   te_panel.pvd    the implicit trailing edge panel, with "mu" (μ_w, the Kutta condition's jump)
#   free.pvd        the free panel rows, with "mu" (frozen when shed) and "row" (0 the newest)
#   particles.pvd   the vortex particles, with "alpha" (strength vector), "strength" (its magnitude) and "core"
#   lift.csv        CL by pressure, Kutta–Joukowski and impulse against time and the wake's length in spans
# In ParaView: open the four .pvd files (File > Open, select them together), press Apply and play. Colour the body by
# Cp and the free rows by mu; show the particles as "Point Gaussian" (Gaussian Radius about 0.01, or scale by "core")
# coloured by "strength", or apply Filters > Glyph with "alpha" as the orientation array to see the vortex direction.

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, MeshLibrary, PAGE
using VortexPanel

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))

const OUTROOT = joinpath(@__DIR__, "..", "output", "examples")

function parse_args(args)
    opts = Dict("level" => "1", "spans" => "5", "max_rows" => "auto", "merge_rows" => "1", "alpha" => "5", "every" => "4", "advection" => "free",
        "dt" => "1")
    for a in args
        key, val = split(a, "=")
        haskey(opts, key) || error("Unknown argument $key; use level, spans, max_rows, merge_rows, alpha, every, advection or dt")
        opts[key] = val
    end
    advection = Symbol(opts["advection"])
    advection in (:free, :freestream) || error("advection must be free or freestream")
    count_or_auto(v) = v == "auto" ? :auto : parse(Int, v)
    return (level=parse(Int, opts["level"]), spans=parse(Float64, opts["spans"]), max_rows=count_or_auto(opts["max_rows"]),
        merge_rows=count_or_auto(opts["merge_rows"]), alpha=parse(Float64, opts["alpha"]), every=parse(Int, opts["every"]), advection=advection,
        dt=parse(Float64, opts["dt"]))
end

function main(args=ARGS)
    o = parse_args(args)
    V, lift = wing_freestream(o.alpha)
    S = WING_SPAN * WING_CHORD
    mesh = reference_wing(o.level; split_te=false)
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    dt = o.dt * VortexPanel.Wake.te_panel_size(g) / norm(V) # rows as long as o.dt trailing edge panels
    sim = FreeWakeSim(g, dm, V; opts=FreeWakeOptions(dt=dt, max_rows=o.max_rows, merge_rows=o.merge_rows, advection=o.advection))
    nsteps = ceil(Int, o.spans * sim.te.span / (norm(V) * dt))
    dir = joinpath(OUTROOT, "free_wake_L$(o.level)_rows$(o.max_rows)" * (o.merge_rows == 1 ? "" : "_merge$(o.merge_rows)") * (o.dt == 1 ? "" : "_dt$(o.dt)"))
    if isdir(dir) # a fresh series: remove an earlier run's frames
        foreach(f -> rm(joinpath(dir, f)), filter(f -> endswith(f, ".vtu") || endswith(f, ".pvd") || f == "lift.csv", readdir(dir)))
    end
    series = FreeWakeSeries(dir; every=o.every)
    println("Free wake on the reference wing L$(o.level): α = $(o.alpha)°, dt = $(round(dt, sigdigits=3)) (|V∞| = 1), $nsteps steps to $(o.spans) spans, ",
        "$(sim.keep_rows) free rows (max_rows = $(o.max_rows)), merge_rows = $(o.merge_rows), advection = $(o.advection), threads $(Threads.nthreads())")
    cl(r, f) = dot(getproperty(r, f), lift) / S
    report = max(1, nsteps ÷ 10)
    progress = (s, r) -> (r.n % report == 0 || r.n == nsteps - 1) &&
        @printf("  step %4d  wake %.2f spans  CL %.5f  CL (KJ) %.5f  %5d particles  %.2f s/step\n", r.n, r.t * norm(V) / s.te.span, cl(r, :F_pressure),
            cl(r, :F_kj), r.particles, r.t_rhs + r.t_vel)
    t = @elapsed run!(sim, nsteps; progress=progress, series=series)

    h = sim.history
    rows = [[r.t r.t * norm(V) / sim.te.span cl(r, :F_pressure) cl(r, :F_kj) cl(r, :F_impulse) r.rows r.particles] for r in h]
    writedlm(joinpath(dir, "lift.csv"), vcat(["t" "wake_spans" "CL_p" "CL_KJ" "CL_impulse" "free_rows" "particles"], reduce(vcat, rows)), ',')
    rigid = solve_wing(mesh; α=o.alpha)
    @printf("\n%d steps in %.0f s. At %.1f spans: CL %.5f, CL (KJ) %.5f; the rigid wake (steady): CL %.5f, CL (KJ) %.5f\n", nsteps, t,
        h[end].t * norm(V) / sim.te.span, cl(h[end], :F_pressure), cl(h[end], :F_kj), rigid.CL, rigid.CL_KJ)
    println("(The lift still rises slowly at this length: the starting vortex's downwash fades like (span/length)².)")
    nframes = length(get(series.frames, "body", []))
    println("\nWrote $nframes frames to ", normpath(dir), ": body.pvd, te_panel.pvd, free.pvd, particles.pvd and lift.csv")
    return sim
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
