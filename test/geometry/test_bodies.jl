function check_closed_sphere(g; nodes, tris, volume, rtol)
    @test G.nnodes(g) == nodes
    @test G.ntris(g) == tris
    @test length(g.bodies) == 1
    @test g.bodies[1].closed
    @test g.bodies[1].euler == 2
    @test all(==(G.Interior), g.edges.class)
    @test isempty(g.openings) && isempty(g.tjunctions) && isempty(g.gaps)
    @test g.quality.counts[:kutta_edges] == 0
    @test g.bodies[1].volume ≈ volume rtol = rtol
    # every normal points away from the centre
    @test all(dot(tri_normal(g.x, t), sum(g.x[t]) / 3) > 0 for t in g.tris)
end

@testset "icosphere" begin
    g = build_geometry(icosphere(3))
    check_closed_sphere(g; nodes=642, tris=1280, volume=4π / 3, rtol=0.02)
    @test g.quality.counts[:welded_nodes] == 0
    files = MeshLibrary.write_vtk(g, joinpath(mktempdir(), "sphere")) # no special edges: triangles only
    @test length(files) == 1 && isfile(files[1])
end

@testset "cube-sphere" begin
    n = 4
    g = build_geometry(cube_sphere(n))
    check_closed_sphere(g; nodes=6n^2 + 2, tris=12n^2, volume=4π / 3, rtol=0.1)
    @test g.quality.counts[:welded_nodes] == 12n - 4 # seam clusters: 12(n-1) pairs and 8 corner triples
    @test g.quality.counts[:split_fold] + g.quality.counts[:split_shortest] + g.quality.counts[:split_tie] == 6n^2
end

@testset "uv sphere with collapsed poles" begin
    nθ, nφ = 8, 12
    g = build_geometry(uv_sphere(nθ, nφ))
    check_closed_sphere(g; nodes=(nθ - 1) * nφ + 2, tris=2nφ * (nθ - 2) + 2nφ, volume=4π / 3, rtol=0.1)
    @test g.quality.counts[:collapsed_elements] == 2nφ
end

@testset "ellipsoid" begin
    g = build_geometry(scale_nodes(icosphere(3), (2.0, 1.0, 0.5)))
    check_closed_sphere(g; nodes=642, tris=1280, volume=4π / 3, rtol=0.02)
end

@testset "inward input is flipped outward" begin
    gout = build_geometry(icosphere(2))
    gin = build_geometry(reverse_elements(icosphere(2)))
    @test gin.quality.counts[:flipped_triangles] == G.ntris(gin)
    @test canonical(gin.tris) == canonical(gout.tris)
end

@testset "generic float type" begin
    g32 = build_geometry(convert_nodes(Float32, cube_sphere(3)))
    @test g32 isa SurfaceGeometry{Float32}
    g64 = build_geometry(cube_sphere(3))
    @test canonical(g32.tris) == canonical(g64.tris)
    @test g32.bodies[1].volume ≈ g64.bodies[1].volume rtol = 1e-5
end
