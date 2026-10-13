# Validation stage 4 in brief: ellipsoids in uniform flow, whose exact Cp follows the curvature (the exact solution is
# validation/common/exact_flows.jl). The full study is validation/ellipsoid.

using QuadGK

const VE = normalize(SVector(1.0, 0.3, -0.2))
const TRIAXIAL = EllipsoidFlow((2.0, 1.0, 0.5), VE)
const PROLATE = EllipsoidFlow((3.0, 0.5, 0.5), VE)
const THIN = EllipsoidFlow((1.0, 1.0, 0.1), VE)

ellipsoid_mesh(flow, level) = scale_nodes(icosphere(level), flow.axes)

"""Solve the flow about the body meshed by `mesh` and measure it against the exact `flow`"""
function ellipsoid_solution(flow, mesh; opts=AssemblyOptions())
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    els = build_elements(g)
    σ = source_strengths(els, flow.V∞)
    A, b, stats = assemble(g, dm, els, σ, flow.V∞; opts=opts)
    sol = solve_galerkin(dm, A, b)
    err = sol.μ - [exact_doublet(flow, g.x[dm.slot_node[s]]) for s in 1:dm.nslot]
    err .-= sum(dm.slot_area .* err) / sum(dm.slot_area) # the null space: compare without the mean
    cperr = pressure_coefficient(gradient_velocity(els, dm, sol.μ, flow.V∞), flow.V∞) - [exact_cp(flow, el.centroid) for el in els]
    area = [el.area for el in els]
    return (g=g, sol=sol, stats=stats, c=consistency(dm, A, b),
        μL2=sqrt(sum(dm.slot_area .* err .^ 2) / sum(dm.slot_area)) / doublet_scale(flow), cpgL2=sqrt(sum(area .* cperr .^ 2) / sum(area)))
end

@testset "exact solution" begin
    @test carlson_rd(0, 2, 1) ≈ 1.7972103521034 rtol = 1e-12 # Carlson (1995)'s test values
    @test carlson_rd(2, 3, 4) ≈ 0.16510527294261 rtol = 1e-12
    @test carlson_rd(1, 1, 1) ≈ 1
    @test carlson_rd(0.3, 2.0, 0.7) ≈ 1.5 * quadgk(t -> 1 / ((t + 0.7) * sqrt((t + 0.3) * (t + 2.0) * (t + 0.7))), 0, Inf; rtol=1e-14)[1] rtol = 1e-12
    for flow in (TRIAXIAL, PROLATE, THIN)
        @test sum(flow.factors) ≈ 2
    end
    e = sqrt(1 - (0.5 / 3)^2) # the prolate spheroid's eccentricity
    @test PROLATE.factors[1] ≈ 2 * (1 - e^2) / e^3 * (atanh(e) - e)
    # the sphere: μ = V∞·x/2 and Cp = 1 - (9/4)sin²θ
    sphere = EllipsoidFlow((1.0, 1.0, 1.0), VE)
    @test all(sphere.factors .≈ 2 / 3)
    x = normalize(SVector(0.3, -0.5, 0.8))
    @test exact_doublet(sphere, x) ≈ dot(VE, x) / 2 && exact_cp(sphere, x) ≈ 1 - 9 / 4 * (1 - dot(VE, x)^2)
    # the surface formulas against the exterior potential: equal on the surface, the Neumann condition, harmonic
    φ(y) = perturbation_potential(TRIAXIAL, y)
    for p in (SVector(0.3, 0.5, 0.7), SVector(-0.8, 0.1, 0.2), SVector(1.0, 0.0, 0.0))
        x = surface_point(TRIAXIAL, p)
        n = surface_normal(TRIAXIAL, x)
        @test φ(x) ≈ exact_doublet(TRIAXIAL, x) atol = 1e-14
        h = 1e-4
        @test (-3φ(x) + 4φ(x + h * n) - φ(x + 2h * n)) / (2h) ≈ -dot(VE, n) atol = 1e-6
        y = x + 0.3n
        @test abs(sum((φ(y + 1e-3 * d) - 2φ(y) + φ(y - 1e-3 * d)) / 1e-6 for d in (SVector(1.0, 0, 0), SVector(0, 1.0, 0), SVector(0, 0, 1.0)))) < 1e-5
    end
    far = normalize(SVector(1.0, 2.0, 3.0))
    @test 1e3^2 * φ(1e3 * far) ≈ 1e4^2 * φ(1e4 * far) rtol = 1e-3 # a dipole far away
end

const E2, E3 = (ellipsoid_solution(TRIAXIAL, ellipsoid_mesh(TRIAXIAL, l)) for l in 2:3)

@testset "triaxial ellipsoid" begin
    @test E3.μL2 < 1e-3 && E2.μL2 / E3.μL2 > 3.5 # second order
    @test E3.cpgL2 < 8e-3 && E2.cpgL2 / E3.cpgL2 > 1.6 # first order
    @test E3.c.rows < 1e-12 && E3.c.columns < 1e-5
    @test E3.sol.cond < 1e3
    @test ellipsoid_solution(TRIAXIAL, scale_nodes(cube_sphere(6), TRIAXIAL.axes)).μL2 < 5e-3
end

@testset "slender and thin bodies" begin
    for (flow, gain) in ((PROLATE, 10), (THIN, 3))
        mesh = ellipsoid_mesh(flow, 2)
        r = ellipsoid_solution(flow, mesh)
        @test r.μL2 < 4e-3 && r.c.columns < 1e-5
        # the corners near 150° of the stretched meshes need the vertex rule split; the contact tolerance is free
        @test ellipsoid_solution(flow, mesh; opts=AssemblyOptions(corner_angle=180.0)).c.columns > gain * r.c.columns
        wide = ellipsoid_solution(flow, mesh; opts=AssemblyOptions(contact_rel=0.3))
        @test wide.stats.pairs.contact > 0
        @test wide.c.columns ≈ r.c.columns rtol = 0.2
    end
end

@testset "abused triaxial ellipsoid" begin
    project = p -> surface_point(TRIAXIAL, p) # along the ray from the centre
    spec = (jitter=0.3, slivers=(5, 0.05), gaps=(3, 0.01), tjunctions=6, project=project)
    r = ellipsoid_solution(TRIAXIAL, first(abuse(ellipsoid_mesh(TRIAXIAL, 2), spec; seed=7)))
    @test length(r.g.bodies) == 1 && r.g.bodies[1].closed
    @test r.g.quality.counts[:tjunctions] == 6 && r.g.quality.counts[:gaps] == 3
    @test r.μL2 < 0.03 && r.sol.cond < 1.05 * E2.sol.cond
end
