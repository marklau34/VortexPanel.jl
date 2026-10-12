@testset "triangle quality" begin
    a, b, c = SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.5, sqrt(3) / 2, 0.0)
    θmin, θmax, aspect, skew, area = G.triangle_quality(a, b, c)
    @test θmin ≈ 60 && θmax ≈ 60
    @test aspect ≈ 1
    @test abs(skew) < 1e-12
    @test area ≈ sqrt(3) / 4
    # a 1° sliver is flagged
    sliver = SVector(0.5, tand(1) / 2, 0.0)
    qmin, _, asp, _, _, flags = G.quality_arrays([a, b, sliver], [SVector(1, 2, 3)], GeometryOptions())
    @test qmin[1] ≈ 1 atol = 1e-9
    @test flags[1]
end

@testset "quad diagonal" begin
    opts = GeometryOptions()
    # planar trapezoid: the shorter diagonal 1-3
    x = [SVector(0.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0), SVector(2.0, 1.0, 0.0), SVector(0.0, 3.0, 0.0)]
    @test G.choose_diagonal(x, [1, 2, 3, 4], 1, opts) == (1, :shortest)
    @test G.choose_diagonal(x, [2, 3, 4, 1], 1, opts) == (2, :shortest) # same diagonal, nodes rotated
    # warped quad: the diagonal with the smaller fold, whatever their lengths
    w = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.2, 1.0, 0.4), SVector(0.0, 1.0, 0.0)]
    d, reason = G.choose_diagonal(w, [1, 2, 3, 4], 1, opts)
    f13, f24 = G.fold_angle(w[1], w[2], w[3], w[4]), G.fold_angle(w[2], w[3], w[4], w[1])
    @test reason == :fold
    @test d == (f13 < f24 ? 1 : 2)
    # square: tied, the mirror rule picks mirror images on either side of y = 0
    sq(y0) = [SVector(0.0, y0, 0.0), SVector(1.0, y0, 0.0), SVector(1.0, y0 + 1, 0.0), SVector(0.0, y0 + 1, 0.0)]
    dr, rr = G.choose_diagonal(sq(1.0), [1, 2, 3, 4], 1, opts)
    dl, rl = G.choose_diagonal(sq(-2.0), [1, 2, 3, 4], 1, opts)
    @test rr == :tie_mirror && rl == :tie_mirror
    @test dr != dl # the mirror image of diagonal 1-3 at y > 0 is diagonal 2-4 at y < 0
    @test G.choose_diagonal(sq(1.0), [1, 2, 3, 4], 1, GeometryOptions(split=:d24)) == (2, :fixed)
end

@testset "weld tolerance is relative" begin
    # two triangles sharing an edge, with one shared node duplicated at a distance d (edge lengths about 1)
    for (scale, d, merged) in ((1.0, 5e-4, true), (1.0, 2e-3, false), (1e-3, 5e-7, true), (1e3, 2.0, false))
        X = scale .* [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(1.0, 1.0, 0.0),
            SVector(1.0, d / scale, 0.0)]
        log = G.RepairEvent{Float64}[]
        w = G.weld(G.ingest(X, [[1, 2, 3], [5, 4, 3]]), 1e-3, log)
        @test (length(w.x) == 4) == merged
    end
end

@testset "non-manifold edge" begin
    X = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0), SVector(0.0, 0.0, 1.0)]
    @test_throws NonManifoldEdgeError build_geometry(X, [[1, 2, 3], [2, 1, 4], [1, 2, 5]])
end
