const WING1 = reference_wing(1; split_te=true)
const G1 = build_geometry(WING1)

@testset "counts" begin
    g = G1
    @test G.nnodes(g) == 684 # 702 with the split trailing edge, welded
    @test g.quality.counts[:welded_nodes] == 18
    @test G.ntris(g) == 1292 + 72 # 646 quads and two 36-triangle caps
    @test G.nedges(g) == 2046
    @test length(g.bodies) == 1 && g.bodies[1].closed && g.bodies[1].euler == 2
    @test isempty(g.tjunctions) && isempty(g.gaps)
    @test length(g.openings) == 2 && all(op -> op.mode == :cap && length(op.cap_tris) == 36, g.openings)
    @test g.quality.counts[:flipped_triangles] == 0 # PAGE's normals are outward
    @test g.bodies[1].volume ≈ reference_wing_volume() rtol = 0.005
end

@testset "trailing edge" begin
    g = G1
    te1 = WING1.boundaries[findfirst(b -> b.name == "trailing_edge_1", WING1.boundaries)]
    tagged = Set(minmax(g.input_node[a], g.input_node[b]) for (a, b) in G.ingest(WING1).edge_tags["trailing_edge_1"])
    kutta = Set(Tuple(g.edges.nodes[e]) for e in findall(==(G.Kutta), g.edges.class))
    @test length(te1) == 17
    @test kutta == tagged
    @test g.kutta.tag_mismatch == 0
    @test length(g.kutta.chains) == 1 && length(g.kutta.te_nodes) == 18
    ch = g.kutta.chains[1]
    @test ch.end_type == (:cap, :cap) && all(>(0), ch.end_cap_tri)
    @test (g.x[ch.nodes[end]] - g.x[ch.nodes[1]])[2] > 0 # τ along +y, so "upper" faces +z
    @test all(p -> 1 <= p <= 323, g.tri_parent[ch.tri_upper])
    @test all(p -> 324 <= p <= 646, g.tri_parent[ch.tri_lower])
    # every corner at a trailing edge node is labelled with its own surface; cap corners as caps
    for t in eachindex(g.tris), k in 1:3
        side = g.kutta.corner_side[t][k]
        if g.tris[t][k] in g.kutta.te_nodes
            expected = g.tri_kind[t] == CapTri ? TECap : g.tri_parent[t] <= 323 ? TEUpper : TELower
            @test side == expected
        else
            @test side == NoTE
        end
    end
    # the same edges from the tags alone, or both
    for source in (:tags, :both)
        gs = build_geometry(WING1; opts=GeometryOptions(kutta_source=source))
        @test Set(Tuple(gs.edges.nodes[e]) for e in findall(==(G.Kutta), gs.edges.class)) == kutta
    end
end

@testset "tip caps" begin
    g = G1
    for op in g.openings
        @test length(op.loop) == 38
        @test g.x[op.le][1] ≈ 0 atol = 1e-12 # leading edge
        @test g.x[op.te][1] ≈ WING_CHORD # trailing edge
        @test op.planarity < 1e-12
        sgn = sign(g.x[op.te][2])
        @test all(t -> tri_normal(g.x, g.tris[t]) ≈ SVector(0, sgn, 0), op.cap_tris) # outward, flat
        @test count(e -> G.CapSeam == g.edges.class[e], 1:G.nedges(g)) == 2 * 38
    end
end

@testset "triangulation invariance" begin
    ref = canonical(G1.tris)
    @test canonical(build_geometry(reference_wing(1; split_te=false)).tris) == ref
    @test canonical(build_geometry(reverse_elements(WING1)).tris) == ref
    for s in (1e-3, 1e3)
        gs = build_geometry(scale_nodes(WING1, s))
        @test canonical(gs.tris) == ref
        @test gs.bodies[1].volume ≈ s^3 * G1.bodies[1].volume rtol = 1e-10
    end
    g32 = build_geometry(convert_nodes(Float32, WING1))
    @test g32 isa SurfaceGeometry{Float32}
    @test canonical(g32.tris) == ref
    # mirror symmetric about y = 0, except the strip that straddles it (tie broken by node id)
    key(p) = round.(p, digits=9)
    cent = Set(key(sum(G1.x[t]) / 3) for t in G1.tris)
    straddles(t) = (y = [G1.x[v][2] for v in G1.tris[t]]; minimum(y) < 0 < maximum(y))
    asym = [t for t in eachindex(G1.tris) if !(key((sum(G1.x[G1.tris[t]]) / 3) .* SVector(1, -1, 1)) in cent)]
    @test all(straddles, asym)
    @test length(asym) <= 2 * 19 * 2 # at most the two triangles of each quad of the middle strip
end

@testset "zip and open modes" begin
    go = build_geometry(WING1; opts=GeometryOptions(openings=:open))
    @test all(op -> op.mode == :open && isempty(op.cap_tris), go.openings)
    @test !go.bodies[1].closed
    @test count(==(G.Open), go.edges.class) == 2 * 38
    @test G.ntris(go) == 1292
    gz = build_geometry(WING1; opts=GeometryOptions(openings=:zip))
    zips = filter(p -> p.kind == :zip, gz.gaps)
    @test length(zips) == 2 && all(op -> op.mode == :zip, gz.openings)
    @test all(z -> z.chain_a[1] == z.chain_b[1] && z.chain_a[end] == z.chain_b[end], zips) # leading to trailing edge
    @test !gz.bodies[1].closed
    @test gz.kutta.chains[1].end_type == (:open, :open)
end

@testset "finer levels" begin
    for level in (2, 3)
        nspan, nchord = reference_wing_size(level)
        g = build_geometry(reference_wing(level))
        @test G.ntris(g) == 4nspan * nchord + 2 * (2nchord - 2)
        @test g.bodies[1].closed && g.bodies[1].euler == 2
        @test g.quality.counts[:kutta_edges] == nspan && length(g.kutta.te_nodes) == nspan + 1
        @test g.kutta.tag_mismatch == 0
        @test g.bodies[1].volume ≈ reference_wing_volume() rtol = 0.005
    end
end

@testset "export" begin
    dir = mktempdir()
    files = MeshLibrary.write_vtk(G1, joinpath(dir, "wing"))
    @test all(isfile, files)
    capped = add_caps(WING1, G1)
    @test nelements(capped) == 646 + 72
    @test get_boundary_names(capped) == get_boundary_names(WING1)
    gc = build_geometry(capped)
    @test isempty(gc.openings)
    @test gc.quality.counts[:flipped_triangles] == 0 # caps follow the input orientation
    @test canonical(gc.tris) == canonical(G1.tris)
end
