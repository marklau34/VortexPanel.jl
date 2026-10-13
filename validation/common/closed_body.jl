# Shared machinery of the closed-body validation stages (the sphere, stage 3, and the ellipsoids, stage 4): solve a
# case against the exact EllipsoidFlow, and print and tabulate convergence, abused meshes and quadrature sensitivity.

using LinearAlgebra, StaticArrays, Printf, DelimitedFiles, MeshLibrary
using VortexPanel, VortexPanel.MeshTools
include(joinpath(@__DIR__, "exact_flows.jl"))

"""A quadrature-converged reference for the sensitivity studies"""
const REFERENCE = AssemblyOptions(far_ratio=8.0, near_ratio=4.0, near_depth=5, touch_levels=10, touch_order=4, kernel_far=Inf)

"""The mesh abuses of the validation stages (doc §11); `project` puts a moved node back on the exact surface"""
abuses(project) = [("jitter 0.3", (jitter=0.3, project=project)), ("slivers 10×0.05", (slivers=(10, 0.05), project=project)),
    ("gaps 6×0.01", (gaps=(6, 0.01), project=project)), ("T-junctions 12", (tjunctions=12,)),
    ("all", (jitter=0.3, slivers=(10, 0.05), gaps=(6, 0.01), tjunctions=12, project=project))]

"""
    solve_case(flow, mesh; opts, velocity=true)

Solve the flow about the body meshed by `mesh` and measure it against the exact `flow`: μ errors at the slots
(relative to the largest exact |μ|, without the mean, which is the null space), Cp from the singularity sums and
from the element gradient at the centroids (absolute), their agreement and the leak (relative to |V∞|), λ, the
column-sum defect of A C, cond₁ and timings.
"""
function solve_case(flow::EllipsoidFlow, mesh; opts=AssemblyOptions(), velocity=true)
    V∞ = flow.V∞
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    els = build_elements(g)
    σ = source_strengths(els, V∞)
    A, b, stats = assemble(g, dm, els, σ, V∞; opts=opts)
    tsolve = @elapsed sol = solve_galerkin(dm, A, b)
    c = consistency(dm, A, b)
    err = sol.μ - [exact_doublet(flow, g.x[dm.slot_node[s]]) for s in 1:dm.nslot]
    err .-= sum(dm.slot_area .* err) / sum(dm.slot_area)
    scale = doublet_scale(flow)
    res = (g=g, dm=dm, els=els, sol=sol, nodes=dm.nslot, tris=length(els), μerr=err,
        μL2=sqrt(sum(dm.slot_area .* err .^ 2) / sum(dm.slot_area)) / scale, μmax=maximum(abs, err) / scale,
        λ=sol.λ[1], colsum=c.columns, cond=sol.cond, tA=stats.time, tsolve=tsolve, pairs=stats.pairs, points=stats.points)
    velocity || return res
    tu = @elapsed u = surface_velocity(els, dm, σ, sol.μ, V∞)
    ug = gradient_velocity(els, dm, sol.μ, V∞)
    area = [el.area for el in els]
    rms(v) = sqrt(sum(area .* v .^ 2) / sum(area))
    cpex = [exact_cp(flow, el.centroid) for el in els]
    cp, cpg = pressure_coefficient(u, V∞), pressure_coefficient(ug, V∞)
    leak = [dot(v, el.n) / norm(V∞) for (v, el) in zip(u, els)]
    tangential(v, n) = v - dot(v, n) * n
    return merge(res, (cp=cp, cpg=cpg, cpex=cpex, leak=leak, cpL2=rms(cp - cpex), cpmax=maximum(abs, cp - cpex),
        cpgL2=rms(cpg - cpex), cpgmax=maximum(abs, cpg - cpex),
        agree=maximum(norm(tangential(u[t], els[t].n) - ug[t]) for t in eachindex(els)) / norm(V∞),
        leak_rms=rms(leak), leak_max=maximum(abs, leak), tu=tu))
end

"""Largest change of μ between two solutions on the same mesh, relative to the largest exact |μ|, without the mean"""
function μchange(flow::EllipsoidFlow, r, ref)
    d = r.sol.μ - ref.sol.μ
    d .-= sum(r.dm.slot_area .* d) / sum(r.dm.slot_area)
    return maximum(abs, d) / doublet_scale(flow)
end

"""Largest angle between a facet normal and the exact surface normal where the ray through the centroid meets it"""
max_tilt(flow::EllipsoidFlow, els) = maximum(acosd(clamp(dot(el.n, surface_normal(flow, surface_point(flow, el.centroid))), -1, 1)) for el in els)

order(e1, e2, n1, n2) = log(e1 / e2) / log(sqrt(n2 / n1)) # observed order in h ∝ 1/√nodes

"""A label as a CSV field: writedlm does not quote, so commas become semicolons"""
csv_label(s) = replace(string(s), "," => ";")

const CONVERGENCE_HEADER = ["family" "nodes" "triangles" "mu_L2" "mu_max" "cp_L2" "cp_max" "cpgrad_L2" "cpgrad_max" "velocity_agreement" "leak_rms" "leak_max" "lambda" "colsum" "cond1" "t_assembly" "t_solve" "t_velocity"]

function convergence(family, flow::EllipsoidFlow, meshes)
    println("\n$family: convergence (μ relative to its largest exact value; Cp absolute; velocities relative to V∞)")
    @printf("%8s %8s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %8s %7s %7s\n", "nodes", "triangles", "μ L2", "μ max", "Cp L2",
        "Cp max", "Cpg L2", "Cpg max", "agree", "leak rms", "λ", "colsum", "cond₁", "t asm", "t vel")
    rows, results = Any[], Any[]
    for mesh in meshes
        r = solve_case(flow, mesh)
        push!(results, r)
        @printf("%8d %8d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.1e %9.1e %8.1e %6.1fs %6.1fs\n", r.nodes, r.tris, r.μL2, r.μmax,
            r.cpL2, r.cpmax, r.cpgL2, r.cpgmax, r.agree, r.leak_rms, r.λ, r.colsum, r.cond, r.tA, r.tu)
        push!(rows, [csv_label(family) r.nodes r.tris r.μL2 r.μmax r.cpL2 r.cpmax r.cpgL2 r.cpgmax r.agree r.leak_rms r.leak_max r.λ r.colsum r.cond r.tA r.tsolve r.tu])
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

const ABUSE_HEADER = ["mesh" "abuse" "nodes" "tjunctions" "gaps" "slivers" "contact_pairs" "mu_L2" "mu_max" "cpgrad_L2" "cp_L2" "leak_rms" "lambda" "colsum" "cond1" "max_tilt_deg" "dmu_quadrature"]

function abused_meshes(name, flow::EllipsoidFlow, base, project; seed=7)
    println("\nAbused meshes from $name (seed $seed): errors as in the convergence tables; Δμ quad is the change of μ from",
        " a quadrature-converged solution on the same mesh")
    @printf("%-16s %6s %3s %3s %4s %7s %9s %9s %9s %9s %9s %9s %9s %8s %8s %10s\n", "abuse", "nodes", "TJ", "gap", "slv", "contact",
        "μ L2", "μ max", "Cpg L2", "Cp L2", "leak rms", "λ", "colsum", "cond₁", "max tilt", "Δμ quad")
    rows, results = Any[], Any[]
    for (label, spec) in [("clean", nothing); abuses(project)]
        mesh = isnothing(spec) ? base : first(abuse(base, spec; seed=seed))
        r = solve_case(flow, mesh)
        dq = μchange(flow, r, solve_case(flow, mesh; opts=REFERENCE, velocity=false))
        q = r.g.quality.counts
        tilt = max_tilt(flow, r.els)
        @printf("%-16s %6d %3d %3d %4d %7d %9.2e %9.2e %9.2e %9.2e %9.2e %9.1e %9.1e %8.1e %7.1f° %10.2e\n", label, r.nodes,
            get(q, :tjunctions, 0), get(q, :gaps, 0), get(q, :slivers, 0), r.pairs.contact, r.μL2, r.μmax, r.cpgL2, r.cpL2,
            r.leak_rms, r.λ, r.colsum, r.cond, tilt, dq)
        push!(rows, [csv_label(name) csv_label(label) r.nodes get(q, :tjunctions, 0) get(q, :gaps, 0) get(q, :slivers, 0) r.pairs.contact r.μL2 r.μmax r.cpgL2 r.cpL2 r.leak_rms r.λ r.colsum r.cond tilt dq])
        push!(results, (label, r))
    end
    return rows, results
end

const QUADRATURE_HEADER = ["mesh" "setting" "dmu_max" "lambda" "colsum" "contact_pairs" "points" "time"]

function quadrature_sensitivity(name, flow::EllipsoidFlow, mesh)
    ref = solve_case(flow, mesh; opts=REFERENCE, velocity=false)
    settings = [("default", AssemblyOptions()),
        ("touch levels 3", AssemblyOptions(touch_levels=3)), ("touch levels 4", AssemblyOptions(touch_levels=4)),
        ("touch levels 8", AssemblyOptions(touch_levels=8)), ("touch order 2", AssemblyOptions(touch_order=2)),
        ("touch order 5", AssemblyOptions(touch_order=5)), ("touch ratio 0.25", AssemblyOptions(touch_ratio=0.25)),
        ("corner 180", AssemblyOptions(corner_angle=180.0)), ("corner 45", AssemblyOptions(corner_angle=45.0)),
        ("contact 0.05", AssemblyOptions(contact_rel=0.05)), ("contact 0.2", AssemblyOptions(contact_rel=0.2)),
        ("far ratio 2", AssemblyOptions(far_ratio=2.0)), ("far ratio 4", AssemblyOptions(far_ratio=4.0)),
        ("far ratio 6", AssemblyOptions(far_ratio=6.0)), ("far points 1", AssemblyOptions(far_points=1)),
        ("near points 6", AssemblyOptions(near_points=6)), ("near depth 0", AssemblyOptions(near_depth=0)),
        ("near ratio 3", AssemblyOptions(near_ratio=3.0)), ("kernel far 5", AssemblyOptions(kernel_far=5.0)),
        ("kernel far Inf", AssemblyOptions(kernel_far=Inf))]
    println("\nQuadrature sensitivity on $name ($(ref.nodes) nodes): change of μ from a quadrature-converged solution, relative",
        " to its largest exact value; discretization error (max) $(@sprintf("%.1e", solve_case(flow, mesh; velocity=false).μmax));",
        " reference λ ", @sprintf("%.1e", ref.λ), ", colsum ", @sprintf("%.1e", ref.colsum))
    @printf("%-18s %12s %10s %10s %8s %10s %8s\n", "setting", "|Δμ|max", "λ", "colsum", "contact", "points", "time")
    rows = Any[]
    for (label, o) in settings
        r = solve_case(flow, mesh; opts=o, velocity=false)
        npts = sum(r.points)
        @printf("%-18s %12.2e %10.1e %10.1e %8d %10d %7.2fs\n", label, μchange(flow, r, ref), r.λ, r.colsum, r.pairs.contact, npts, r.tA)
        push!(rows, [csv_label(name) csv_label(label) μchange(flow, r, ref) r.λ r.colsum r.pairs.contact npts r.tA])
    end
    return rows
end

"""VTK of a solved case: Cp (both routes) and its exact value, the leak, μ, its exact value and its error"""
function write_case(flow::EllipsoidFlow, r, outdir, prefix)
    nn = length(r.g.x)
    MeshLibrary.write_vtk(r.g, joinpath(outdir, prefix); cell_data=Dict("Cp" => r.cp, "Cp_gradient" => r.cpg,
            "Cp_exact" => r.cpex, "leak" => r.leak), point_data=Dict("mu" => r.sol.μ[1:nn],
            "mu_exact" => [exact_doublet(flow, x) for x in r.g.x], "mu_error" => r.μerr[1:nn]))
end
