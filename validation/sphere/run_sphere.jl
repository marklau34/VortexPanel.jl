# Validation stage 3: a sphere in uniform flow (doc §11). The exact surface doublet is μ = V∞·x/2 (linear in the
# streamwise coordinate; the outer perturbation potential with the interior at rest) and Cp = 1 - (9/4)sin²θ.
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

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, MeshLibrary
using VortexPanel, VortexPanel.MeshTools

const OUTDIR = joinpath(@__DIR__, "..", "output", "sphere")
const VINF = normalize(SVector(1.0, 0.3, -0.2)) # a generic direction: no symmetry hides errors
const REFERENCE = AssemblyOptions(far_ratio=8.0, near_ratio=4.0, near_depth=5, touch_levels=10, touch_order=4, kernel_far=Inf)
const ONSPHERE = p -> normalize(p) # the meshes are unit spheres: moved nodes stay on the exact surface
const JITTER = (jitter=0.3, project=ONSPHERE)
const ABUSES = [("jitter 0.3", JITTER), ("slivers 10×0.05", (slivers=(10, 0.05), project=ONSPHERE)),
    ("gaps 6×0.01", (gaps=(6, 0.01), project=ONSPHERE)), ("T-junctions 12", (tjunctions=12,)),
    ("all", (jitter=0.3, slivers=(10, 0.05), gaps=(6, 0.01), tjunctions=12, project=ONSPHERE))]

function parse_args(args)
    levels, cubes = [1, 2, 3, 4], [3, 6, 12, 24]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
        key == "cube" && (cubes = parse.(Int, split(val, ",")))
    end
    return levels, cubes
end

"""Solve the sphere problem and measure everything; `opts` are the assembly options"""
function sphere_case(mesh; opts=AssemblyOptions(), velocity=true)
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    els = build_elements(g)
    σ = source_strengths(els, VINF)
    A, b, stats = assemble(g, dm, els, σ, VINF; opts=opts)
    tsolve = @elapsed sol = solve_galerkin(dm, A, b)
    c = consistency(dm, A, b)
    μex = [dot(VINF, g.x[dm.slot_node[s]]) / 2 for s in 1:dm.nslot]
    err = sol.μ - μex
    err .-= sum(dm.slot_area .* err) / sum(dm.slot_area)
    res = (g=g, dm=dm, els=els, sol=sol, nodes=dm.nslot, tris=length(els), μerr=err,
        μL2=sqrt(sum(dm.slot_area .* err .^ 2) / sum(dm.slot_area)) / 0.5, μmax=maximum(abs, err) / 0.5,
        λ=sol.λ[1], colsum=c.columns, cond=sol.cond, tA=stats.time, tsolve=tsolve, pairs=stats.pairs, points=stats.points)
    velocity || return res
    tu = @elapsed u = surface_velocity(els, dm, σ, sol.μ, VINF)
    ug = gradient_velocity(els, dm, sol.μ, VINF)
    area = [el.area for el in els]
    rms(v) = sqrt(sum(area .* v .^ 2) / sum(area))
    cpex = [1 - 9 / 4 * (1 - dot(VINF, normalize(el.centroid))^2) for el in els]
    cp, cpg = pressure_coefficient(u, VINF), pressure_coefficient(ug, VINF)
    leak = [dot(v, el.n) for (v, el) in zip(u, els)]
    tangential(v, n) = v - dot(v, n) * n
    return merge(res, (cp=cp, cpg=cpg, cpex=cpex, leak=leak, cpL2=rms(cp - cpex), cpmax=maximum(abs, cp - cpex),
        cpgL2=rms(cpg - cpex), cpgmax=maximum(abs, cpg - cpex),
        agree=maximum(norm(tangential(u[t], els[t].n) - ug[t]) for t in eachindex(els)),
        leak_rms=rms(leak), leak_max=maximum(abs, leak), tu=tu))
end

"""Largest change of μ between two solutions on the same mesh, relative to the peak of μ (V∞/2), without the mean"""
function μchange(r, ref)
    d = r.sol.μ - ref.sol.μ
    d .-= sum(r.dm.slot_area .* d) / sum(r.dm.slot_area)
    return maximum(abs, d) / 0.5
end

order(e1, e2, n1, n2) = log(e1 / e2) / log(sqrt(n2 / n1)) # observed order in h ∝ 1/√nodes

function convergence(family, meshes)
    println("\n$family: convergence (μ relative to its peak V∞/2; Cp absolute; velocities relative to V∞)")
    @printf("%8s %8s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %8s %7s %7s\n", "nodes", "triangles", "μ L2", "μ max", "Cp L2",
        "Cp max", "Cpg L2", "Cpg max", "agree", "leak rms", "λ", "colsum", "cond₁", "t asm", "t vel")
    rows, results = Any[], Any[]
    for mesh in meshes
        r = sphere_case(mesh)
        push!(results, r)
        @printf("%8d %8d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.1e %9.1e %8.1e %6.1fs %6.1fs\n", r.nodes, r.tris, r.μL2, r.μmax,
            r.cpL2, r.cpmax, r.cpgL2, r.cpgmax, r.agree, r.leak_rms, r.λ, r.colsum, r.cond, r.tA, r.tu)
        push!(rows, [family r.nodes r.tris r.μL2 r.μmax r.cpL2 r.cpmax r.cpgL2 r.cpgmax r.agree r.leak_rms r.leak_max r.λ r.colsum r.cond r.tA r.tsolve r.tu])
    end
    if length(results) > 1
        @printf("observed orders in h: ")
        for k in 2:length(results)
            a, b = results[k-1], results[k]
            @printf("[%d→%d nodes] μ L2 %.2f, Cp(grad) L2 %.2f, leak %.2f   ", a.nodes, b.nodes, order(a.μL2, b.μL2, a.nodes, b.nodes),
                order(a.cpgL2, b.cpgL2, a.nodes, b.nodes), order(a.leak_rms, b.leak_rms, a.nodes, b.nodes))
        end
        println()
    end
    return rows, results
end

function abused_meshes(name, base; seed=7)
    println("\nAbused meshes from $name (seed $seed): errors as in the convergence tables; Δμ quad is the change of μ from",
        " a quadrature-converged solution on the same mesh")
    @printf("%-16s %6s %3s %3s %4s %7s %9s %9s %9s %9s %9s %9s %9s %8s %8s %10s\n", "abuse", "nodes", "TJ", "gap", "slv", "contact",
        "μ L2", "μ max", "Cpg L2", "Cp L2", "leak rms", "λ", "colsum", "cond₁", "max tilt", "Δμ quad")
    rows, results = Any[], Any[]
    for (label, spec) in [("clean", nothing); ABUSES]
        mesh = isnothing(spec) ? base : first(abuse(base, spec; seed=seed))
        r = sphere_case(mesh)
        dq = μchange(r, sphere_case(mesh; opts=REFERENCE, velocity=false))
        q = r.g.quality.counts
        tilt = maximum(acosd(clamp(dot(el.n, normalize(el.centroid)), -1, 1)) for el in r.els) # facet normal against the sphere's
        @printf("%-16s %6d %3d %3d %4d %7d %9.2e %9.2e %9.2e %9.2e %9.2e %9.1e %9.1e %8.1e %7.1f° %10.2e\n", label, r.nodes,
            get(q, :tjunctions, 0), get(q, :gaps, 0), get(q, :slivers, 0), r.pairs.contact, r.μL2, r.μmax, r.cpgL2, r.cpL2,
            r.leak_rms, r.λ, r.colsum, r.cond, tilt, dq)
        push!(rows, [name label r.nodes get(q, :tjunctions, 0) get(q, :gaps, 0) get(q, :slivers, 0) r.pairs.contact r.μL2 r.μmax r.cpgL2 r.cpL2 r.leak_rms r.λ r.colsum r.cond tilt dq])
        push!(results, (label, r))
    end
    return rows, results
end

function quadrature_sensitivity(name, mesh)
    ref = sphere_case(mesh; opts=REFERENCE, velocity=false)
    settings = [("default", AssemblyOptions()),
        ("touch levels 3", AssemblyOptions(touch_levels=3)), ("touch levels 4", AssemblyOptions(touch_levels=4)),
        ("touch levels 8", AssemblyOptions(touch_levels=8)), ("touch order 2", AssemblyOptions(touch_order=2)),
        ("touch order 5", AssemblyOptions(touch_order=5)), ("touch ratio 0.25", AssemblyOptions(touch_ratio=0.25)),
        ("contact 0.05", AssemblyOptions(contact_rel=0.05)), ("contact 0.2", AssemblyOptions(contact_rel=0.2)),
        ("far ratio 2", AssemblyOptions(far_ratio=2.0)), ("far ratio 4", AssemblyOptions(far_ratio=4.0)),
        ("far ratio 6", AssemblyOptions(far_ratio=6.0)), ("far points 1", AssemblyOptions(far_points=1)),
        ("near points 6", AssemblyOptions(near_points=6)), ("near depth 0", AssemblyOptions(near_depth=0)),
        ("near ratio 3", AssemblyOptions(near_ratio=3.0)), ("kernel far 5", AssemblyOptions(kernel_far=5.0)),
        ("kernel far Inf", AssemblyOptions(kernel_far=Inf))]
    println("\nQuadrature sensitivity on $name ($(ref.nodes) nodes): change of μ from a quadrature-converged solution, relative",
        " to its peak; discretization error (max) $(@sprintf("%.1e", sphere_case(mesh; velocity=false).μmax)); reference λ ",
        @sprintf("%.1e", ref.λ), ", colsum ", @sprintf("%.1e", ref.colsum))
    @printf("%-18s %12s %10s %10s %8s %10s %8s\n", "setting", "|Δμ|max", "λ", "colsum", "contact", "points", "time")
    rows = Any[]
    for (label, o) in settings
        r = sphere_case(mesh; opts=o, velocity=false)
        npts = sum(r.points)
        @printf("%-18s %12.2e %10.1e %10.1e %8d %10d %7.2fs\n", label, μchange(r, ref), r.λ, r.colsum, r.pairs.contact, npts, r.tA)
        push!(rows, [name label μchange(r, ref) r.λ r.colsum r.pairs.contact npts r.tA])
    end
    return rows
end

function test_weights(levels)
    println("\nTest weights: μ L2 error (and Cp from the gradient) on icospheres")
    rows = Any[]
    for test in (:hat, :onering, :dual)
        rs = [sphere_case(icosphere(l); opts=AssemblyOptions(test=test)) for l in levels]
        @printf("%-8s %s\n", test, join([@sprintf("%5d nodes: μ %.2e Cp %.2e", r.nodes, r.μL2, r.cpgL2) for r in rs], "   "))
        for r in rs
            push!(rows, [string(test) r.nodes r.μL2 r.cpgL2])
        end
    end
    return rows
end

"""VTK of a solved case: Cp (both routes) and its exact value, the leak, μ, its exact value and its error"""
function write_case(r, prefix)
    nn = length(r.g.x)
    MeshLibrary.write_vtk(r.g, joinpath(OUTDIR, prefix); cell_data=Dict("Cp" => r.cp, "Cp_gradient" => r.cpg,
            "Cp_exact" => r.cpex, "leak" => r.leak), point_data=Dict("mu" => r.sol.μ[1:nn],
            "mu_exact" => [dot(VINF, x) / 2 for x in r.g.x], "mu_error" => r.μerr[1:nn]))
end

function main(args=ARGS)
    levels, cubes = parse_args(args)
    mkpath(OUTDIR)
    println("Sphere in uniform flow, V∞ = $(round.(VINF, digits=4)), threads $(Threads.nthreads())")
    rows_ico, ico = convergence("icosphere", [icosphere(l) for l in levels])
    rows_cube, _ = convergence("cube-sphere", [cube_sphere(n) for n in cubes])
    rows_ico_j, _ = convergence("icosphere, jitter 0.3", [first(abuse(icosphere(l), JITTER; seed=1)) for l in levels])
    rows_cube_j, _ = convergence("cube-sphere, jitter 0.3", [first(abuse(cube_sphere(n), JITTER; seed=1)) for n in cubes])
    header = ["family" "nodes" "triangles" "mu_L2" "mu_max" "cp_L2" "cp_max" "cpgrad_L2" "cpgrad_max" "velocity_agreement" "leak_rms" "leak_max" "lambda" "colsum" "cond1" "t_assembly" "t_solve" "t_velocity"]
    writedlm(joinpath(OUTDIR, "convergence.csv"), vcat(header, reduce(vcat, [rows_ico; rows_cube; rows_ico_j; rows_cube_j])), ',')

    rows_ab_ico, ab_ico = abused_meshes("icosphere 3", icosphere(3))
    rows_ab_cube, _ = abused_meshes("cube-sphere 12", cube_sphere(12))
    writedlm(joinpath(OUTDIR, "abuse.csv"), vcat(["mesh" "abuse" "nodes" "tjunctions" "gaps" "slivers" "contact_pairs" "mu_L2" "mu_max" "cpgrad_L2" "cp_L2" "leak_rms" "lambda" "colsum" "cond1" "max_tilt_deg" "dmu_quadrature"],
        reduce(vcat, [rows_ab_ico; rows_ab_cube])), ',')

    abused = first(abuse(icosphere(3), (jitter=0.3, tjunctions=12, project=ONSPHERE); seed=7))
    sens = [quadrature_sensitivity("icosphere 3", icosphere(3)); quadrature_sensitivity("cube-sphere 12", cube_sphere(12));
        quadrature_sensitivity("icosphere 3, jitter + T-junctions", abused)]
    writedlm(joinpath(OUTDIR, "quadrature.csv"), vcat(["mesh" "setting" "dmu_max" "lambda" "colsum" "contact_pairs" "points" "time"], reduce(vcat, sens)), ',')

    tw = test_weights(levels[1:min(end, 3)])
    writedlm(joinpath(OUTDIR, "test_weights.csv"), vcat(["test" "nodes" "mu_L2" "cpgrad_L2"], reduce(vcat, tw)), ',')

    write_case(ico[end], "icosphere_$(levels[end])")
    write_case(last(ab_ico[end]), "icosphere_3_abused")
    println("\nWrote ", normpath(OUTDIR))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
