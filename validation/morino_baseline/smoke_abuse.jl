# One MESA solve per mesh abuse operation on the level 1 reference wing (open tips): does MESA run, and how much does
# the lift change from the clean mesh? Also what VortexPanel's geometry layer finds in each mesh. The full skew study
# runs with the Galerkin solver at validation stage 6.
#
#   julia -t auto --project=validation validation/morino_baseline/smoke_abuse.jl

using MESA, MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles, StableRNGs
using VortexPanel
using VortexPanel.MeshTools

include(joinpath(@__DIR__, "..", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "common", "loads.jl"))
include(joinpath(@__DIR__, "common_mesa.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "morino")

"""Abused level 1 wings (split_te=false form, as the generator expects), by name"""
function abused_wings(; seed::Integer=1)
    nspan, nchord = reference_wing_size(1)
    wing = reference_wing(1; split_te=false)
    project = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS)
    rng() = StableRNG(seed)
    skew(pattern, σ; kw...) = first(skew_wing(wing, pattern, σ; nspan=nspan, nchord=nchord, rng=rng(), kw...))
    return [
        "clean" => wing,
        "jitter 0.2" => first(jitter_nodes(wing, 0.2; rng=rng(), project=project)),
        "slivers (3 nodes, 0.02)" => first(insert_slivers(wing, 3, 0.02; rng=rng(), project=project)),
        "sliver stations (2, 0.02)" => first(sliver_stations(wing, 2, 0.02; nspan=nspan, nchord=nchord, rng=rng())),
        "gaps (3, 0.01)" => first(open_gaps(wing, 3, 0.01; rng=rng(), project=project)),
        "T-junctions (4)" => first(insert_tjunctions(wing, 4; rng=rng())),
        "lean opposite 0.3" => skew(:lean_opposite, 0.3),
        "jitter independent 0.3" => skew(:jitter_independent, 0.3),
        "TE misaligned 0.3, shared" => skew(:te_misaligned, 0.3; shared_te=true),
        "TE misaligned 0.3, non-matching" => skew(:te_misaligned, 0.3; shared_te=false),
        "lean opposite 0.3, TE rows" => skew(:lean_opposite, 0.3; te_rows=true),
    ]
end

function main()
    mkpath(OUTDIR)
    nspan, nchord = reference_wing_size(1)
    rows = Any[]
    CL0 = NaN
    @printf("\n%-34s %-6s %9s %9s %8s   %s\n", "mesh", "MESA", "CL", "CL_KJ", "ΔCL %", "VortexPanel geometry")
    for (name, mesh) in abused_wings()
        g = build_geometry(mesh)
        c = g.quality.counts
        found = @sprintf("TJ %d, gaps %d, non-matching TE %d, Kutta edges %d, slivers %d", c[:tjunctions], c[:gaps], c[:kutta_pairs],
            c[:kutta_edges], c[:slivers])
        status, reason, CL, CL_KJ = "ok", "", NaN, NaN
        stage = "trailing edge split"
        try
            m = deepcopy(mesh)
            PAGE.split_trailing_edge!(m, "trailing_edge")
            stage = "MESA"
            r = mesa_case(m; wake_len=20WING_SPAN, nspan=nspan, nchord=nchord)
            CL, CL_KJ = r.CL, r.CL_KJ
            isfinite(CL) || (status = "NaN")
        catch e
            frames = [string(fr.func) for fr in stacktrace(catch_backtrace())]
            where = something(findfirst(f -> f in ("calc_CHTLS_pan!", "split_trailing_edge!", "create_trailing_edge", "init", "solve!"), frames), 0)
            status = stage == "MESA" ? "fails" : "n/a"
            reason = string(nameof(typeof(e))) * (where > 0 ? " in " * frames[where] : "")
        end
        name == "clean" && (CL0 = CL)
        dCL = 100 * (CL / CL0 - 1)
        @printf("%-34s %-6s %9.5f %9.5f %8.3f   %s%s\n", name, status, CL, CL_KJ, dCL, found, isempty(reason) ? "" : "; MESA: $reason")
        push!(rows, [name status reason CL CL_KJ dCL c[:tjunctions] c[:gaps] c[:kutta_pairs] c[:kutta_edges] c[:slivers]])
    end
    header = ["mesh" "mesa" "mesa_error" "CL" "CL_KJ" "dCL_percent" "tjunctions" "gaps" "kutta_pairs" "kutta_edges" "slivers"]
    writedlm(joinpath(OUTDIR, "smoke_abuse.csv"), vcat(header, reduce(vcat, rows)), ',')
    println("\nWrote ", normpath(joinpath(OUTDIR, "smoke_abuse.csv")))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
