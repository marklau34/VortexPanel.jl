# Validation stage 3: a sphere in uniform flow (doc §11). The exact surface doublet is μ = V∞·x/2 (linear in the
# streamwise coordinate; the outer perturbation potential, the inner one being zero for σ = -V∞·n) and
# Cp = 1 - (9/4)sin²θ: the ellipsoid of validation/common/exact_flows.jl with equal axes.
#
#   julia -t auto --project=validation validation/sphere/run_sphere.jl [levels=1,2,3,4] [cube=3,6,12,24]
#
# 1. Convergence under refinement on icospheres and on cube-spheres (quad meshes split into triangles), clean and
#    jittered: μ, Cp from the singularities and from the element gradient (the two routes of doc §9), their agreement,
#    the leak (normal velocity at the centroids), λ, the column-sum defect of A C (the quadrature's consistency),
#    cond₁ and timings.
# 2. Abused meshes (doc §11: every stage on a clean and an abused mesh): jitter, slivers, gaps, T-junctions and all of
#    them at once, against the clean mesh, with the quadrature error of each.
# 3. Sensitivity to the quadrature tiers: each setting against a quadrature-converged solution on the same mesh, on
#    the clean meshes and on an abused one. λ measures the quadrature only on the abused mesh: the clean meshes are
#    centrally symmetric, which makes the right hand side compatible to rounding whatever the quadrature.
# 4. Test weights: hat (Galerkin), one-ring and median dual cell.
# Writes validation/output/sphere/ (CSV tables, VTK of the finest icosphere and of the fully abused icosphere 3).

include(joinpath(@__DIR__, "..", "common", "closed_body.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "sphere")
const SPHERE = EllipsoidFlow((1.0, 1.0, 1.0), normalize(SVector(1.0, 0.3, -0.2))) # a generic V∞: no symmetry hides errors
const ONSPHERE = p -> normalize(p) # the meshes are unit spheres: moved nodes stay on the exact surface
const JITTER = (jitter=0.3, project=ONSPHERE)

function parse_args(args)
    levels, cubes = [1, 2, 3, 4], [3, 6, 12, 24]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "cube" && (cubes = parse.(Int, split(val, ",")))
    end
    return levels, cubes
end

function test_weights(levels)
    println("\nTest weights: μ L2 error (and Cp from the gradient) on icospheres")
    rows = Any[]
    for test in (:hat, :onering, :dual)
        rs = [solve_case(SPHERE, icosphere(l); opts=AssemblyOptions(test=test)) for l in levels]
        @printf("%-8s %s\n", test, join([@sprintf("%5d nodes: μ %.2e Cp %.2e", r.nodes, r.μL2, r.cpgL2) for r in rs], "   "))
        for r in rs
            push!(rows, [string(test) r.nodes r.μL2 r.cpgL2])
        end
    end
    return rows
end

function main(args=ARGS)
    levels, cubes = parse_args(args)
    mkpath(OUTDIR)
    println("Sphere in uniform flow, V∞ = $(round.(SPHERE.V∞, digits=4)), threads $(Threads.nthreads())")
    rows_ico, ico = convergence("icosphere", SPHERE, [icosphere(l) for l in levels])
    rows_cube, _ = convergence("cube-sphere", SPHERE, [cube_sphere(n) for n in cubes])
    rows_ico_j, _ = convergence("icosphere, jitter 0.3", SPHERE, [first(abuse(icosphere(l), JITTER; seed=1)) for l in levels])
    rows_cube_j, _ = convergence("cube-sphere, jitter 0.3", SPHERE, [first(abuse(cube_sphere(n), JITTER; seed=1)) for n in cubes])
    writedlm(joinpath(OUTDIR, "convergence.csv"), vcat(CONVERGENCE_HEADER, reduce(vcat, [rows_ico; rows_cube; rows_ico_j; rows_cube_j])), ',')

    rows_ab_ico, ab_ico = abused_meshes("icosphere 3", SPHERE, icosphere(3), ONSPHERE)
    rows_ab_cube, _ = abused_meshes("cube-sphere 12", SPHERE, cube_sphere(12), ONSPHERE)
    writedlm(joinpath(OUTDIR, "abuse.csv"), vcat(ABUSE_HEADER, reduce(vcat, [rows_ab_ico; rows_ab_cube])), ',')

    abused = first(abuse(icosphere(3), (jitter=0.3, tjunctions=12, project=ONSPHERE); seed=7))
    sens = [quadrature_sensitivity("icosphere 3", SPHERE, icosphere(3)); quadrature_sensitivity("cube-sphere 12", SPHERE, cube_sphere(12));
        quadrature_sensitivity("icosphere 3, jitter + T-junctions", SPHERE, abused)]
    writedlm(joinpath(OUTDIR, "quadrature.csv"), vcat(QUADRATURE_HEADER, reduce(vcat, sens)), ',')

    tw = test_weights(levels[1:min(end, 3)])
    writedlm(joinpath(OUTDIR, "test_weights.csv"), vcat(["test" "nodes" "mu_L2" "cpgrad_L2"], reduce(vcat, tw)), ',')

    write_case(SPHERE, ico[end], OUTDIR, "icosphere_$(levels[end])")
    write_case(SPHERE, last(ab_ico[end]), OUTDIR, "icosphere_3_abused")
    println("\nWrote ", normpath(OUTDIR))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
