# Example: the Galerkin solver on an abused reference wing (NACA 0015, aspect ratio 3.88, capped tips, rigid wake), with
# the surface coloured by pressure in ParaView. Solves the clean mesh too, for comparison.
#
#   julia -t auto --project=validation validation/examples/abused_wing.jl [level=2] [abuse=all] [alpha=5] [seed=1] [count=4]
#
# abuse: all (jitter, slivers, gaps and T-junctions at once), jitter, slivers, gaps, tjunctions (`count` quads split in
# two, each leaving two hanging nodes on its neighbours' edges), tjunctions_midspan (the same, all within |η| ≤ 0.2 of
# midspan), skew (opposite spanwise lean of the two surfaces), te_misaligned (a non-matching trailing edge, which MESA
# cannot run) or clean.
#
# Writes to validation/output/examples/, for the clean and the abused mesh (prefix wing_L<level>_<abuse>, with the
# count for T-junctions):
#   *_tris.vtu       the body. Colour by "Cp" (constant on each triangle, the solver's pressure) or "Cp_nodes" (averaged
#                    to the nodes for smooth shading; at the trailing edge it mixes the two surfaces). Also "mu" (the
#                    doublet strength at the nodes; at the trailing edge the upper copy's), the triangle "kind" (1 wing,
#                    2 cap), quality ("min_angle", "aspect", "sliver") and node flags ("hanging": T-junctions).
#   *_edges.vtu      the edges that are not plain interior edges, by "class": 2 trailing edge, 3 non-matching trailing
#                    edge, 4 gap, 5/6 T-junction (long edge / short edges), 8 cap seam.
#   *_wake.vtu       the wake strips with "mu" (μ_w); *_filaments.vtu the tip filaments.
# In ParaView: open the *_tris.vtu files, colour by Cp, and show "Surface With Edges" to see the abuse.

using LinearAlgebra, StaticArrays, Printf, StableRNGs, MeshLibrary, PAGE
using VortexPanel, VortexPanel.MeshTools
using VortexPanel.Post: nodal_average

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "wing_case.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "examples")
const ABUSES = (:all, :jitter, :slivers, :gaps, :tjunctions, :tjunctions_midspan, :skew, :te_misaligned, :clean)
const MIDSPAN_BAND = 0.2 # tjunctions_midspan: elements whose centroid has |η| ≤ this

function parse_args(args)
    opts = Dict("level" => "2", "abuse" => "all", "alpha" => "5", "seed" => "1", "count" => "4")
    for a in args
        key, val = split(a, "=")
        haskey(opts, key) || error("Unknown argument $key; use level, abuse, alpha, seed or count")
        opts[key] = val
    end
    abuse_kind = Symbol(opts["abuse"])
    abuse_kind in ABUSES || error("Unknown abuse $abuse_kind; use one of $(join(ABUSES, ", "))")
    return parse(Int, opts["level"]), abuse_kind, parse(Float64, opts["alpha"]), parse(Int, opts["seed"]), parse(Int, opts["count"])
end

"""Elements of `mesh` whose centroid lies farther than `band` from midspan, in η = 2y/b"""
function outside_band(mesh, band)
    η(e) = 2 * sum(p -> p[2], get_coords(mesh, mesh.elements[e])) / (length(get_nodes(mesh.elements[e])) * WING_SPAN)
    return Set(e for e in eachindex(mesh.elements) if abs(η(e)) > band)
end

"""
The reference wing of `level` (in the `split_te=false` form the abuse generator expects) with abuse `kind`; `count` is
the number of quads split into T-junctions
"""
function abused_wing(level::Integer, kind::Symbol; seed::Integer=1, count::Integer=4)
    nspan, nchord = reference_wing_size(level)
    wing = reference_wing(level; split_te=false)
    project = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS) # moved nodes stay on the exact surface
    rng = StableRNG(seed)
    kind == :clean && return wing
    kind == :jitter && return first(jitter_nodes(wing, 0.2; rng=rng, project=project))
    kind == :slivers && return first(insert_slivers(wing, 3, 0.02; rng=rng, project=project))
    kind == :gaps && return first(open_gaps(wing, 3, 0.01; rng=rng, project=project))
    kind == :tjunctions && return first(insert_tjunctions(wing, count; rng=rng))
    if kind == :tjunctions_midspan # by default the generator also avoids elements at the trailing edge
        avoid = union(VortexPanel.MeshTools.default_avoid(wing), outside_band(wing, MIDSPAN_BAND))
        return first(insert_tjunctions(wing, count; rng=rng, avoid=avoid))
    end
    kind == :skew && return first(skew_wing(wing, :lean_opposite, 0.3; nspan=nspan, nchord=nchord, rng=rng))
    kind == :te_misaligned && return first(skew_wing(wing, :te_misaligned, 0.3; nspan=nspan, nchord=nchord, rng=rng, shared_te=false))
    return first(abuse(wing, (jitter=0.2, slivers=(3, 0.02), gaps=(3, 0.01), tjunctions=4, project=project); seed=seed))
end

"""Write the body and the wake of a `solve_wing` result as VTK, coloured by Cp; returns the file names"""
function write_solution(r, prefix)
    nn = length(r.g.x)
    body = MeshLibrary.write_vtk(r.g, prefix; cell_data=Dict("Cp" => r.cp),
        point_data=Dict("Cp_nodes" => nodal_average(r.g, r.els, r.cp), "mu" => r.sol.μ[1:nn]))
    wake = MeshLibrary.write_vtk(r.wake, prefix; point_data=Dict("mu" => r.μ_w[r.wake.te]))
    return [body; wake]
end

function main(args=ARGS)
    level, kind, α, seed, count = parse_args(args)
    tj = kind in (:tjunctions, :tjunctions_midspan)
    name = tj ? "$(kind)_$(count)" : string(kind)
    mkpath(OUTDIR)
    println("Reference wing L$level, abuse: $name (seed $seed), α = $(α)°, threads $(Threads.nthreads())\n")

    clean = solve_wing(reference_wing(level; split_te=false); α=α)
    r = kind == :clean ? clean : solve_wing(abused_wing(level, kind; seed=seed, count=count); α=α)

    print_report(stdout, r.g) # what the geometry layer found and repaired
    if !isempty(r.g.tjunctions)
        η = [2 * r.g.x[t.hanging][2] / WING_SPAN for t in r.g.tjunctions]
        @printf("  %-22s %d, at η from %+.3f to %+.3f\n", "hanging nodes", length(η), extrema(η)...)
    end
    println()
    @printf("%-22s %9s %9s %10s %10s %9s %8s %10s %9s %7s\n", "mesh", "CL", "CL_KJ", "CDi_p", "CDi_T", "Cm_c4", "unknowns", "col. defect", "cond₁", "time")
    for (label, s) in (("clean", clean), (name, r))
        @printf("%-22s %9.5f %9.5f %10.6f %10.6f %9.5f %8d %10.1e %9.1e %6.1fs\n", label, s.CL, s.CL_KJ, s.CDi_p, s.CDi_T, s.Cm_c4, s.nq,
            s.columns, s.cond, s.t_assembly + s.t_solve)
    end
    @printf("lift change from the clean mesh: %+.4f%% (pressure), %+.4f%% (Kutta–Joukowski)\n\n", 100 * (r.CL / clean.CL - 1), 100 * (r.CL_KJ / clean.CL_KJ - 1))

    files = write_solution(clean, joinpath(OUTDIR, "wing_L$(level)_clean"))
    kind == :clean || append!(files, write_solution(r, joinpath(OUTDIR, "wing_L$(level)_$(name)")))
    println("Wrote:")
    foreach(f -> println("  ", normpath(f)), files)
    return r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
