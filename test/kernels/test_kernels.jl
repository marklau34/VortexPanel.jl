using ForwardDiff
using VortexPanel.Kernels

include(joinpath(@__DIR__, "reference.jl"))

# test triangles in a general position: canonical shapes, rotated and translated
const ROT = let a = normalize(SVector(0.3, -0.7, 0.65)), θ = 1.1
    K = SMatrix{3,3}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)
    SMatrix{3,3}(I + sin(θ) * K + (1 - cos(θ)) * K * K)
end
const SHIFT = SVector(0.4, -1.2, 2.5)
place(p) = ROT * p + SHIFT
const SHAPES = [
    "equilateral" => (SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(0.5, sqrt(3) / 2, 0)),
    "right" => (SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(0.0, 1, 0)),
    "sliver (aspect 100)" => (SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(0.5, 0.01, 0)),
    "skewed (175° angle)" => (SVector(0.0, 0, 0), SVector(1.0, 0, 0), SVector(1.9, 0.08, 0)),
]

"""Test points around a triangle: (label, point, side) by category"""
function test_points(tri)
    a, b = tri.v[1], tri.v[2]
    n, h = tri.n, tri.diam
    ν1 = tri.ν[1]
    mid = (a + b) / 2
    bis = normalize(tri.ν[1] + tri.ν[2]) # outward at vertex b
    return [
        ("far", tri.centroid + 20h * normalize(ROT * SVector(1.0, 2.0, 0.5)), 0),
        ("near above", tri.centroid + 0.1h * n, 0),
        ("near below", tri.centroid - 0.1h * n + 0.2h * normalize(b - a), 0),
        ("very near", tri.centroid + 1e-6h * n, 0),
        ("in-plane inside, side +1", tri.centroid, 1),
        ("in-plane inside, side -1", tri.centroid, -1),
        ("in-plane inside, principal value", tri.centroid + 0.1h * normalize(b - a), 0),
        ("in-plane outside", tri.centroid + 2h * normalize(ν1 + 0.4 * normalize(b - a)), 0),
        ("in-plane on an edge's extension", b + 0.5 * (b - a), 0),
        ("edge-adjacent above", mid + 1e-3h * n, 0),
        ("edge-adjacent, 1e-8 above", mid + 1e-8h * n, 0),
        ("edge-adjacent in-plane inside", mid - 1e-3h * ν1, 1),
        ("edge-adjacent in-plane outside", mid + 1e-3h * ν1, 0),
        ("vertex-adjacent above", b + 1e-3h * normalize(n + bis), 0),
        ("vertex-adjacent in-plane outside", b + 1e-3h * bis, 0),
    ]
end

@testset "closed form against polar quadrature" begin
    for (name, (a, b, c)) in SHAPES
        tri = SourceTriangle(place(a), place(b), place(c))
        for (label, x, side) in test_points(tri)
            w = source_velocity(tri, x; side=side, far=Inf)
            r = w_polar(tri.v..., x; side=side)
            # 1e-8 diameters from an edge, rounding the placed point (about 1e-16 of its distance from the origin)
            # already moves it by ~1e-8 of its height, which bounds any agreement (the field varies like its log)
            rtol = label == "edge-adjacent, 1e-8 above" ? 1e-7 : 1e-9
            @test isapprox(w, r; rtol=rtol, atol=1e-13) || (println("$name, $label: $w vs $r"); false)
        end
    end
end

@testset "closed form against cubature" begin
    for (name, (a, b, c)) in SHAPES
        tri = SourceTriangle(place(a), place(b), place(c))
        for x in (tri.centroid + 0.5tri.diam * tri.n, tri.centroid - 2tri.diam * normalize(tri.n + tri.ν[2]), tri.v[2] + 0.2tri.diam * normalize(tri.n - tri.ν[3]))
            @test isapprox(source_velocity(tri, x; far=Inf), w_cubature(tri.v..., x); rtol=1e-10)
        end
    end
end

@testset "limits on the triangle's plane" begin
    for (name, (a, b, c)) in SHAPES
        tri = SourceTriangle(place(a), place(b), place(c))
        n = tri.n
        # inside: the normal velocity is ±1/2, so it jumps by 1 across the triangle; the principal value has none
        up, down, pv = (source_velocity(tri, tri.centroid; side=s) for s in (1, -1, 0))
        @test dot(up, n) ≈ 0.5 && dot(down, n) ≈ -0.5 && abs(dot(pv, n)) < 1e-14
        @test up - down ≈ n
        # the side limit agrees with points just off the plane
        @test isapprox(up, source_velocity(tri, tri.centroid + 1e-9tri.diam * n); rtol=1e-6)
        # the one-sided normal limit is (fraction of the angle seeing the triangle)/2: 1/4 on an edge, α/4π at a vertex
        # of interior angle α, 0 outside (the tangential part is singular there, so the fraction is tested directly)
        @test Kernels.inplane_fraction(tri, (tri.v[1] + tri.v[2]) / 2, 1e-12) ≈ 0.5
        α = acos(dot(normalize(tri.v[3] - tri.v[1]), normalize(tri.v[2] - tri.v[1])))
        @test Kernels.inplane_fraction(tri, tri.v[1], 1e-12) ≈ α / 2π
        @test abs(dot(source_velocity(tri, tri.centroid + 3tri.diam * tri.ν[1]; side=1), n)) < 1e-15
    end
end

@testset "far field" begin
    for (name, (a, b, c)) in SHAPES
        tri = SourceTriangle(place(a), place(b), place(c))
        for u in (normalize(SVector(1.0, 0.2, -0.3)), normalize(SVector(-0.4, 1.0, 0.9)), tri.n)
            x = tri.centroid + FAR_DEFAULT * tri.diam * u
            exact = source_velocity_exact(tri, x, 0, 1e-12)
            @test norm(source_velocity_far(tri, x) - exact) / norm(exact) < 5e-6 # continuity at the switch
            x2 = tri.centroid + 2FAR_DEFAULT * tri.diam * u
            ex2 = source_velocity_exact(tri, x2, 0, 1e-12)
            @test norm(source_velocity_far(tri, x2) - ex2) / norm(ex2) < 5e-6 / 10 # fourth order: ÷16 per doubling
            @test source_velocity(tri, x2) == source_velocity_far(tri, x2)
        end
    end
end

@testset "consistency" begin
    a, b, c = (place(p) for p in SHAPES[4][2])
    tri = SourceTriangle(a, b, c)
    flipped = SourceTriangle(a, c, b)
    @test flipped.n ≈ -tri.n
    for (label, x, side) in test_points(tri)
        label == "far" && continue
        rtol = label == "edge-adjacent, 1e-8 above" ? 1e-7 : 1e-12 # conditioned by rounding, see above
        @test isapprox(source_velocity(flipped, x; side=-side, far=Inf), source_velocity(tri, x; side=side, far=Inf); rtol=rtol, atol=1e-14)
    end
    # Gauss: the flux out of a cube around the triangle is its area
    L = 1.5tri.diam
    xg, wg = gauss(24, -L, L)
    flux = 0.0
    for axis in 1:3, s in (-1.0, 1.0), (xi, wi) in zip(xg, wg), (xj, wj) in zip(xg, wg)
        i, j = [k for k in 1:3 if k != axis]
        p = zeros(3)
        p[axis], p[i], p[j] = s * L, xi, xj
        flux += wi * wj * s * source_velocity(tri, tri.centroid + SVector{3}(p))[axis]
    end
    @test flux ≈ tri.area rtol = 1e-10
    @test_throws ArgumentError SourceTriangle(a, b, a + 2(b - a))
end

@testset "generic number types" begin
    a, b, c = (place(p) for p in SHAPES[3][2])
    tri = SourceTriangle(a, b, c)
    # off the triangle w is a gradient of a harmonic function: divergence and curl vanish
    for x in (tri.centroid + 0.05tri.diam * tri.n, b + 0.02tri.diam * normalize(tri.n + tri.ν[1]), tri.centroid + 30tri.diam * normalize(tri.n + tri.ν[3]))
        J = ForwardDiff.jacobian(y -> source_velocity(tri, SVector{3}(y)), x)
        @test abs(tr(J)) < 1e-12 * norm(J) && norm(J - J') < 1e-12 * norm(J)
    end
    # derivatives with respect to the vertices
    x0 = tri.centroid + 0.3tri.diam * tri.n
    f(p) = source_velocity(SourceTriangle(SVector{3}(p), b, c), x0; far=Inf)
    J = ForwardDiff.jacobian(f, Vector(a))
    h = 1e-6
    Jfd = reduce(hcat, [(f(Vector(a) + h * e) - f(Vector(a) - h * e)) / 2h for e in eachcol(Matrix(1.0I, 3, 3))])
    @test isapprox(J, Jfd; rtol=1e-7)
    # single precision stays single precision, and agrees with double precision on the same (rounded) input
    a32, b32, c32, x32 = Float32.(a), Float32.(b), Float32.(c), Float32.(x0)
    t32 = SourceTriangle(a32, b32, c32)
    @test source_velocity(t32, x32) isa SVector{3,Float32}
    @test source_velocity(t32, x32 + 30t32.diam * t32.n) isa SVector{3,Float32} # far field
    @test isapprox(source_velocity(t32, x32), source_velocity(SourceTriangle(Float64.(a32), Float64.(b32), Float64.(c32)), Float64.(x32)); rtol=1e-4)
end

@testset "vortex sheet" begin
    a, b, c = (place(p) for p in SHAPES[1][2])
    tri = SourceTriangle(a, b, c)
    γ = 2.0 * normalize(b - a) - 0.7 * tri.ν[1]
    x = tri.centroid + 0.3tri.diam * tri.n
    @test sheet_velocity(tri, γ, x) ≈ cross(γ, source_velocity(tri, x))
    # the tangential velocity jumps by γ × n across the sheet
    @test sheet_velocity(tri, γ, tri.centroid; side=1) - sheet_velocity(tri, γ, tri.centroid; side=-1) ≈ cross(γ, tri.n)
end
