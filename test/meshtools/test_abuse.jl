const WING0 = reference_wing(1; split_te=false)
const NSPAN, NCHORD = reference_wing_size(1)
const CLEAN = build_geometry(WING0)
const PROJECT = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS)

on_naca(p) = abs(p[3] - PROJECT(p)[3]) < 1e-12

@testset "wing grid" begin
    grid = wing_grid(WING0; nspan=NSPAN, nchord=NCHORD)
    X = [p.coord for p in WING0.nodes]
    @test all(grid.upper[1, :] .== grid.lower[1, :]) # shared leading edge
    @test all(grid.upper[end, :] .== grid.lower[end, :]) # shared trailing edge (split_te=false)
    for j in 1:NSPAN+1
        @test all(k -> X[grid.upper[k, j]][2] ≈ X[grid.upper[1, j]][2], 1:NCHORD+1)
        @test all(k -> X[grid.upper[k, j]][3] >= -1e-12 && X[grid.lower[k, j]][3] <= 1e-12, 1:NCHORD+1)
    end
end

@testset "reproducible" begin
    spec = (jitter=0.2, slivers=(2, 0.02), gaps=(2, 0.01), tjunctions=3, project=PROJECT)
    m1, l1 = abuse(WING0, spec; seed=7)
    m2, l2 = abuse(WING0, spec; seed=7)
    m3, _ = abuse(WING0, spec; seed=8)
    @test [p.coord for p in m1.nodes] == [p.coord for p in m2.nodes]
    @test [e.con for e in m1.elements] == [e.con for e in m2.elements]
    @test l1.counts == l2.counts
    @test [p.coord for p in m1.nodes] != [p.coord for p in m3.nodes]
    # every requested pathology is created and found by the geometry layer
    g = build_geometry(m1)
    @test l1.counts[:tjunctions] == 6 && length(g.tjunctions) == 6
    @test l1.counts[:gaps] == 2 && count(p -> p.kind == :gap, g.gaps) == 2
    @test l1.counts[:slivers] == 2
    @test g.quality.counts[:kutta_edges] == NSPAN # trailing edge untouched
    @test length(g.bodies) == 1 && g.bodies[1].closed
    @test get_boundary_names(m1) == get_boundary_names(WING0)
    @test all(on_naca, [p.coord for p in m1.nodes][1:length(WING0.nodes)])
end

@testset "many seeds" begin
    spec = (jitter=0.2, gaps=(2, 0.01), tjunctions=3, project=PROJECT)
    for seed in 1:25
        m, log = abuse(WING0, spec; seed=seed)
        g = build_geometry(m) # no non-manifold edges, whatever the seed
        @test count(p -> p.kind == :gap, g.gaps) == log.counts[:gaps] == 2
        @test length(g.tjunctions) == log.counts[:tjunctions] == 6
        @test g.bodies[1].closed && g.quality.counts[:kutta_edges] == NSPAN
    end
end

@testset "T-junctions" begin
    m, log = insert_tjunctions(WING0, 4; rng=StableRNG(3))
    g = build_geometry(m)
    @test length(g.tjunctions) == 8 == log.counts[:tjunctions]
    @test all(tj -> tj.t ≈ 0.5, g.tjunctions)
    @test nelements(m) == nelements(WING0) + 4
    # the new elements join their parents' element sets
    @test all(s -> length(s.elems) >= length(WING0.elemset[findfirst(t -> t.name == s.name, WING0.elemset)].elems), m.elemset)
    mi, li = insert_tjunctions(icosphere(2), 5; rng=StableRNG(3))
    gi = build_geometry(mi)
    @test length(gi.tjunctions) == 5 == li.counts[:tjunctions]
    @test gi.bodies[1].closed && isempty(gi.openings)
end

@testset "gaps" begin
    m, log = open_gaps(WING0, 3, 0.01; rng=StableRNG(4), project=PROJECT)
    g = build_geometry(m)
    @test log.counts[:gaps] == 3
    @test count(p -> p.kind == :gap, g.gaps) == 3
    @test length(g.openings) == 2 # only the tips
    @test g.bodies[1].closed
    mi, _ = open_gaps(icosphere(2), 1, 0.01; rng=StableRNG(4))
    gi = build_geometry(mi)
    @test count(p -> p.kind == :gap, gi.gaps) == 1 && isempty(gi.openings)
end

@testset "slivers" begin
    m, log = insert_slivers(WING0, 3, 0.02; rng=StableRNG(5))
    g = build_geometry(m)
    @test log.counts[:slivers] == 3
    @test g.quality.stats[:min_angle] < CLEAN.quality.stats[:min_angle]
    @test g.quality.counts[:welded_nodes] == 0
    ms, ls = sliver_stations(WING0, 2, 0.02; nspan=NSPAN, nchord=NCHORD, rng=StableRNG(5))
    gs = build_geometry(ms)
    @test ls.counts[:slivers] == 4NCHORD
    @test gs.quality.counts[:slivers] >= CLEAN.quality.counts[:slivers] + NCHORD
    @test gs.quality.stats[:max_aspect] > CLEAN.quality.stats[:max_aspect]
    @test gs.bodies[1].volume ≈ CLEAN.bodies[1].volume rtol = 1e-12 # stations move along the span only
end

@testset "spanwise skew" begin
    for pattern in (:lean_same, :lean_opposite, :jitter_same, :jitter_independent)
        m, _ = skew_wing(WING0, pattern, 0.3; nspan=NSPAN, nchord=NCHORD, rng=StableRNG(6))
        g = build_geometry(m)
        @test G.nnodes(g) == G.nnodes(CLEAN) # leading and trailing edges stay shared
        @test g.quality.counts[:kutta_edges] == NSPAN && isempty(g.kutta.pairs)
        @test g.bodies[1].volume ≈ CLEAN.bodies[1].volume rtol = 1e-12 # exact geometry
    end
    # opposite skew of the trailing edge rows: non-matching trailing edge
    m, log = skew_wing(WING0, :lean_opposite, 0.3; nspan=NSPAN, nchord=NCHORD, te_rows=true)
    g = build_geometry(m)
    @test log.counts[:split_nodes] == NSPAN - 1
    @test g.quality.counts[:kutta_edges] == 0 && length(g.kutta.pairs) == 1
    @test length(g.kutta.pair_chains[1].union_nodes) == 2NSPAN
    @test g.bodies[1].closed && length(g.openings) == 2
    # MESA's trailing edge misalignment study
    for (shared, npairs, nkutta) in ((true, 0, NSPAN), (false, 1, 0))
        mt, _ = skew_wing(WING0, :te_misaligned, 0.3; nspan=NSPAN, nchord=NCHORD, shared_te=shared)
        gt = build_geometry(mt)
        @test length(gt.kutta.pairs) == npairs && gt.quality.counts[:kutta_edges] == nkutta
    end
end

@testset "jitter stays on the surface" begin
    m, _ = jitter_nodes(WING0, 0.3; rng=StableRNG(9), project=PROJECT)
    @test all(on_naca, [p.coord for p in m.nodes])
    g = build_geometry(m)
    @test g.quality.counts[:kutta_edges] == NSPAN
    @test g.bodies[1].volume ≈ CLEAN.bodies[1].volume rtol = 0.01
    mf, _ = jitter_nodes(icosphere(3), 0.3; rng=StableRNG(9))
    @test build_geometry(mf).bodies[1].closed
end
