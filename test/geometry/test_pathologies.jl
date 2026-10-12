@testset "T-junctions" begin
    X, faces = quad_cube(split_face=true)
    g = build_geometry(X, faces)
    @test length(g.tjunctions) == 2
    @test all(tj -> tj.t ≈ 0.5, g.tjunctions)
    @test Set(tj.hanging for tj in g.tjunctions) == Set(g.input_node[[9, 10]])
    @test count(==(G.TJLong), g.edges.class) == 2 && count(==(G.TJShort), g.edges.class) == 4
    @test isempty(g.openings) && isempty(g.gaps)
    @test length(g.bodies) == 1 && g.bodies[1].closed
    @test g.bodies[1].volume ≈ 1
    @test length(g.touch_pairs) >= 2
    @test all(dot(tri_normal(g.x, t), sum(g.x[t]) / 3 .- 0.5) > 0 for t in g.tris)
end

@testset "small gap" begin
    mesh = cube_sphere(4)
    top = mesh.elemset[findfirst(s -> s.name == "face_3+", mesh.elemset)].elems
    topnodes = unique(reduce(vcat, [collect(mesh.elements[e].con) for e in top]))
    ε = 0.005 # beyond the weld tolerance, within the gap tolerance (edges ≈ 0.4)
    for v in topnodes
        mesh.nodes[v] = MeshLibrary.Point((1 + ε) * mesh.nodes[v].coord)
    end
    g = build_geometry(mesh)
    @test length(g.gaps) == 1 && g.gaps[1].kind == :gap
    @test g.gaps[1].width ≈ ε rtol = 0.05
    @test length(g.gaps[1].chain_a) == length(g.gaps[1].chain_b) == 17 # closed loops of 16 edges
    @test isempty(g.openings) && isempty(g.tjunctions)
    @test count(==(G.Gap), g.edges.class) == 32
    @test length(g.bodies) == 1 && g.bodies[1].closed
    @test g.quality.counts[:kutta_edges] == 0
end

@testset "gap beyond tolerance stays an opening" begin
    X, faces = quad_cube()
    # the x = 1 face on its own nodes, 0.2 away: its loop and the hole's are too far apart to be a gap, so both are capped
    lid = faces[6]
    Xs = [X; [X[v] + SVector(0.2, 0.0, 0.0) for v in lid]]
    g = build_geometry(Xs, [faces[1:5]; [collect(length(X)+1:length(X)+4)]])
    @test isempty(g.gaps)
    @test length(g.openings) == 2
    @test length(g.bodies) == 2
end

@testset "non-matching trailing edge" begin
    mesh = reference_wing(1; split_te=true)
    copies = 685:702 # lower surface copies of the trailing edge nodes added by PAGE's split
    nspan = 17
    ys = sort(unique(round.(getindex.(mesh.nodes.coord[copies], 2), digits=12)))
    dy = minimum(diff(ys))
    for v in copies
        p = mesh.nodes[v].coord
        abs(abs(p[2]) - WING_SPAN / 2) < 1e-9 && continue # tips stay shared
        mesh.nodes[v] = MeshLibrary.Point(p + SVector(0.0, 0.3dy, 0.0))
    end
    g = build_geometry(mesh)
    @test g.quality.counts[:kutta_edges] == 0
    @test length(g.kutta.pairs) == 1 && isempty(g.gaps)
    pc = g.kutta.pair_chains[1]
    @test length(pc.upper) == nspan + 1 && length(pc.lower) == nspan + 1
    @test length(pc.union_nodes) == 2nspan # tips shared
    @test count(==(0), pc.union_side) == 2
    @test issorted(pc.s)
    @test all(p -> 1 <= p <= 323, g.tri_parent[pc.tri_upper])
    @test all(p -> 324 <= p <= 646, g.tri_parent[pc.tri_lower])
    @test count(==(G.KuttaPair), g.edges.class) == 2nspan
    @test length(g.openings) == 2 && all(op -> op.mode == :cap, g.openings)
    @test g.bodies[1].closed
    # corners on each side of the trailing edge are labelled with their own surface
    for t in eachindex(g.tris), k in 1:3
        side = g.kutta.corner_side[t][k]
        side == TEUpper && @test g.tri_parent[t] <= 323
        side == TELower && @test g.tri_parent[t] > 323
    end
end

@testset "sharp leading edge is not a trailing edge" begin
    X, faces = diamond_wing(n=4)
    g = build_geometry(X, faces)
    kutta = findall(==(G.Kutta), g.edges.class)
    @test length(kutta) == 4
    @test all(e -> all(v -> g.x[v][1] ≈ 1, g.edges.nodes[e]), kutta)
    @test length(g.openings) == 2 && g.bodies[1].closed
    @test all(op -> g.x[op.te][1] ≈ 1 && g.x[op.le][1] ≈ 0, g.openings)
end
