# Validation stage 3 in brief: a sphere in uniform flow. The exact doublet is μ = V∞·x/2 (the outer perturbation
# potential, linear in the streamwise coordinate) and Cp = 1 - (9/4)sin²θ. The full study is validation/sphere.

using VortexPanel.DOFMap: bordered_rows

const VINF = normalize(SVector(1.0, 0.3, -0.2)) # a generic direction, so no symmetry hides errors

"""Solve the sphere problem on `mesh` and measure the errors"""
function sphere_solution(mesh; opts=AssemblyOptions())
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    els = build_elements(g)
    σ = source_strengths(els, VINF)
    A, b, stats = assemble(g, dm, els, σ, VINF; opts=opts)
    sol = solve_galerkin(dm, A, b)
    μex = [dot(VINF, g.x[dm.slot_node[s]]) / 2 for s in 1:dm.nslot]
    err = sol.μ - μex
    err .-= sum(dm.slot_area .* err) / sum(dm.slot_area) # the null space: compare without the mean
    u = surface_velocity(els, dm, σ, sol.μ, VINF)
    ug = gradient_velocity(els, dm, sol.μ, VINF)
    area = [el.area for el in els]
    rms(v) = sqrt(sum(area .* v .^ 2) / sum(area))
    cpex = [1 - 9 / 4 * (1 - dot(VINF, normalize(el.centroid))^2) for el in els]
    tangential(v, n) = v - dot(v, n) * n
    return (g=g, dm=dm, els=els, A=A, b=b, sol=sol, stats=stats, cp=pressure_coefficient(u, VINF),
        μL2=sqrt(sum(dm.slot_area .* err .^ 2) / sum(dm.slot_area)) / 0.5, μmax=maximum(abs, err) / 0.5,
        cpL2=rms(pressure_coefficient(u, VINF) - cpex), cpgL2=rms(pressure_coefficient(ug, VINF) - cpex),
        agree=maximum(norm(tangential(u[t], els[t].n) - ug[t]) for t in eachindex(els)),
        leak=rms([dot(v, el.n) for (v, el) in zip(u, els)]))
end

const S1, S2, S3 = (sphere_solution(icosphere(level)) for level in 1:3)

@testset "icosphere" begin
    r = S2
    @test r.μL2 < 3e-3 && r.μmax < 1e-2
    @test S1.μL2 / S2.μL2 > 3 && S2.μL2 / S3.μL2 > 3.5 # second order
    @test r.cpgL2 < 0.025 && r.cpL2 < 0.08
    @test S2.cpgL2 / S3.cpgL2 > 1.7 # first order
    @test r.agree < 0.05 && r.leak < 0.045 # the two velocity routes agree, the leak is O(h)
    @test S2.leak / S3.leak > 1.7
end

@testset "system diagnostics" begin
    c = consistency(S2.dm, S2.A, S2.b)
    @test c.rows < 1e-12 # a uniform doublet sheds nothing, exactly
    @test c.columns < 1e-5 # no net flux of a hat function's sheet: the quadrature error
    @test c.rhs < 1e-12 # no net source flux
    @test abs(S2.sol.λ[1]) < 1e-14 && S2.sol.residual < 1e-12 # compatible to rounding: the mesh is centrally symmetric
    shifted = solve_galerkin(S2.dm, S2.A, S2.b .+ 1e-3) # λ takes up an incompatible right hand side
    @test shifted.λ[1] ≈ 1e-3
    @test abs(dot(bordered_rows(S2.dm)[1].c, S2.sol.q)) < 1e-12 # zero mean doublet
    @test S2.sol.cond < 1e3
    @test S2.stats.pairs.self == length(S2.els) && sum(S2.stats.pairs) == length(S2.els)^2
end

@testset "quadrature" begin
    fine = sphere_solution(icosphere(2); opts=AssemblyOptions(far_ratio=8.0, near_ratio=4.0, near_depth=5, touch_levels=10, touch_order=4, kernel_far=Inf))
    @test maximum(abs, S2.sol.μ - fine.sol.μ) / maximum(abs, fine.sol.μ) < 5e-5 # far below the discretization error
    serial = sphere_solution(icosphere(2); opts=AssemblyOptions(threaded=false))
    @test serial.A ≈ S2.A && serial.b ≈ S2.b
end

@testset "test weights" begin
    for test in (:onering, :dual)
        r = sphere_solution(icosphere(2); opts=AssemblyOptions(test=test))
        @test r.μL2 < 1.5e-2 && r.cpgL2 < 0.04
    end
    @test_throws ArgumentError assemble(S1.g, S1.dm, S1.els, source_strengths(S1.els, VINF), VINF; opts=AssemblyOptions(test=:box))
end

@testset "cube-sphere (split quads)" begin
    r = sphere_solution(cube_sphere(4))
    @test r.μL2 < 1.2e-2 && r.cpgL2 < 0.08
    @test r.g.quality.counts[:welded_nodes] > 0 && r.g.bodies[1].closed
end

@testset "abused meshes" begin
    onsphere(p) = normalize(p) # moved nodes stay on the unit sphere
    abused(spec; level=2, seed=7) = first(abuse(icosphere(level), spec; seed=seed))
    # T-junctions reproduce a linear μ, and the contact tier keeps the quadrature as consistent as on the clean mesh
    tj = sphere_solution(abused((tjunctions=6,)))
    @test tj.g.quality.counts[:tjunctions] == 6 && tj.stats.pairs.contact > 0
    @test consistency(tj.dm, tj.A, tj.b).columns < 1e-5
    @test tj.μL2 < 1.05 * S2.μL2
    # gaps: the error stays at the clean level; the flux through the slits makes the right hand side incompatible
    gp = sphere_solution(abused((gaps=(3, 0.01), project=onsphere)))
    @test gp.g.quality.counts[:gaps] == 3 && gp.g.bodies[1].closed
    @test gp.μL2 < 1.05 * S2.μL2
    @test abs(gp.sol.λ[1]) > 1e-8
    # slivers: a local loss of accuracy, but the system is as well conditioned and as consistent as the clean one
    sl = sphere_solution(abused((slivers=(5, 0.05), project=onsphere)))
    @test sl.g.quality.counts[:slivers] > 0 && sl.stats.pairs.contact > 0
    @test sl.μL2 < 5 * S2.μL2 && sl.sol.cond < 1.05 * S2.sol.cond
    @test consistency(sl.dm, sl.A, sl.b).columns < 1e-5
    # jitter: still second order; without central symmetry λ measures the quadrature and shrinks as it is refined
    j2, j3 = (sphere_solution(abused((jitter=0.3, project=onsphere); level=l, seed=1)) for l in 2:3)
    @test j2.μL2 / j3.μL2 > 3.2 && j3.μL2 < 4 * S3.μL2
    fine = sphere_solution(abused((jitter=0.3, project=onsphere); seed=1); opts=AssemblyOptions(far_ratio=8.0, near_ratio=4.0,
        near_depth=5, touch_levels=10, touch_order=4, kernel_far=Inf))
    @test abs(j2.sol.λ[1]) > 1e-10 && abs(fine.sol.λ[1]) < abs(j2.sol.λ[1]) / 100
    # everything at once
    al = sphere_solution(abused((jitter=0.3, slivers=(5, 0.05), gaps=(3, 0.01), tjunctions=6, project=onsphere)))
    @test length(al.g.bodies) == 1 && al.g.bodies[1].closed
    @test al.g.quality.counts[:tjunctions] == 6 && al.g.quality.counts[:gaps] == 3
    @test al.μL2 < 0.03 && al.sol.cond < 1.05 * S2.sol.cond
end

@testset "export" begin
    files = MeshLibrary.write_vtk(S2.g, joinpath(mktempdir(), "sphere"); cell_data=Dict("Cp" => S2.cp),
        point_data=Dict("mu" => S2.sol.μ[1:length(S2.g.x)]))
    @test length(files) == 1 && isfile(files[1])
end
