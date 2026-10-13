# The isoparametric bilinear quad (validation stage 7): shape functions and geometry, the closest point, the inner
# integrals (quad_field) against adaptive cubature, the principal value against the two-sided limits, and the rules on
# the parameter square.

using HCubature
using VortexPanel.Elements
using VortexPanel.Assembly: vertex_square, square_contact_rule, square_tiers

const QEL = VortexPanel.Elements

const QUADS = [
    "square" => IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0), SVector(0.0, 1.0, 0.0)),
    "trapezoid" => IsoQuad(SVector(0.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0), SVector(1.4, 0.8, 0.0), SVector(0.4, 0.8, 0.0)),
    "warped" => IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.1, 0.2), SVector(1.2, 1.1, -0.1), SVector(-0.1, 0.9, 0.15)),
    "sliver" => IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 0.02, 0.0), SVector(0.0, 0.02, 0.0)),
    "sheared sliver" => IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.6, 0.06, 0.01), SVector(0.6, 0.06, 0.0)), # 5.7° corners
]

integrate_square(r, f) = sum(w * f(p) for (p, w) in zip(r.p, r.w))

"""quad_field's sheet and source parts as one vector"""
field5(Q, x, nx, V; kw...) = (f = quad_field(Q, x, nx, V; kw...); SVector(f[1]..., f[2]))

"""quad_field by adaptive cubature: the parameter square split at the closest point P and graded toward it"""
function field_cubature(Q, x, nx, V; rtol=1e-9)
    P, _ = project(Q, x)
    total = zero(SVector{5,Float64})
    for C in QEL.SQUARE
        D = C - P
        (abs(D[1]) < 1e-15 || abs(D[2]) < 1e-15) && continue
        f(u) = begin
            c, s = QEL.kernel_terms(Q, x, nx, V, P + D .* SVector(u[1]^3, u[2]^3))
            abs(D[1] * D[2]) * 9 * u[1]^2 * u[2]^2 * SVector(c..., s)
        end
        total += first(hcubature(f, SA[0.0, 0.0], SA[1.0, 1.0]; rtol=rtol, atol=1e-13, maxevals=5_000_000))
    end
    return total / 4π
end

@testset "shape functions and geometry" begin
    for (ξ, η) in ((0.3, 0.7), (0.0, 1.0), (0.91, 0.12))
        @test sum(shape(ξ, η)) ≈ 1
        Nξ, Nη = shape_derivatives(ξ, η)
        @test abs(sum(Nξ)) < 1e-15 && abs(sum(Nη)) < 1e-15
    end
    for (_, Q) in QUADS
        @test all(k -> point(Q, QEL.SQUARE[k]) ≈ Q.x[k], 1:4)
        J(u) = norm(area_normal(Q, SVector(u[1], u[2])))
        @test Q.area ≈ first(hcubature(J, SA[0.0, 0.0], SA[1.0, 1.0]; rtol=1e-12)) rtol = 1e-6
        @test Q.centroid ≈ point(Q, SVector(0.5, 0.5))
        p = SVector(0.3, 0.6)
        @test norm(sum(sheet_densities(Q, p))) < 1e-14 # a uniform doublet sheds nothing
        # μ = k·x is bilinear in (ξ, η) on the quad's own surface, so its interpolant is exact there: ∇μ is the
        # tangential part of k at every point, and γ = n × ∇μ
        k = SVector(0.7, -1.3, 0.4)
        μ = SVector(ntuple(a -> dot(k, Q.x[a]), 4))
        for p in (SVector(0.3, 0.6), SVector(0.95, 0.05), SVector(0.0, 1.0))
            n = unit_normal(Q, p)
            @test QEL.gradient(Q, μ, p) ≈ k - dot(k, n) * n
            @test sheet_strength(Q, μ, p) ≈ cross(n, k - dot(k, n) * n)
        end
    end
end

@testset "closest point" begin
    rng = StableRNG(3)
    for (name, Q) in QUADS
        # the closest point on a 401 × 401 parameter grid, refined around its best node
        function brute(x)
            best = Inf
            for a in 0:400, b in 0:400
                best = min(best, norm(point(Q, SVector(a / 400, b / 400)) - x))
            end
            return best
        end
        for _ in 1:30
            x = Q.centroid + Q.diam * (SVector(rand(rng, 3)...) .- 0.5) .* SVector(1.6, 1.6, 0.4)
            P, d = project(Q, x)
            @test all(0 .<= P .<= 1)
            @test norm(point(Q, P) - x) ≈ d
            @test d <= brute(x) + 1e-12 # never misses the closest point (the grid's own error is about 1e-6 here)
        end
    end
end

@testset "inner integrals against cubature" begin
    V = SVector(0.8, -0.3, 0.5)
    for (name, Q) in QUADS
        for (p0, offset) in ((SVector(0.37, 0.58), SVector(0.0, 0.0)), (SVector(0.5, 0.0), SVector(0.0, -0.02)), (SVector(1.0, 1.0), SVector(0.01, 0.015)))
            # above an interior point, just outside the middle of side 1, just outside corner 3
            n = unit_normal(Q, p0)
            base = point(Q, p0 + offset)
            nx = normalize(n + 0.4 * normalize(Q.x[2] - Q.x[1]))
            for d in (3.0, 0.3, 1e-2, 1e-4, 1e-6)
                x = base + d * Q.diam * n
                ref = field_cubature(Q, x, nx, V)
                @test norm(field5(Q, x, nx, V) - ref) <= 1e-6 * norm(ref)
            end
        end
    end
end

@testset "principal value and one-sided limits" begin
    V = SVector(0.8, -0.3, 0.5)
    for (name, Q) in QUADS, p0 in (SVector(0.37, 0.58), SVector(0.08, 0.9))
        n = unit_normal(Q, p0)
        x0 = point(Q, p0)
        nx = normalize(n + 0.5 * normalize(Q.x[2] - Q.x[1])) # tilted: the tangential jump of the sheet enters
        h = 1e-7 * Q.diam
        up, down = field5(Q, x0 + h * n, nx, V), field5(Q, x0 - h * n, nx, V)
        pv, plus, minus = field5(Q, x0, nx, V), field5(Q, x0, nx, V; side=1), field5(Q, x0, nx, V; side=-1)
        scale = norm(up) + norm(down)
        @test norm(pv - (up + down) / 2) <= 1e-5 * scale
        @test norm(plus - up) <= 1e-5 * scale && norm(minus - down) <= 1e-5 * scale
        @test plus[5] - minus[5] ≈ -dot(V, n) * dot(nx, n) # the source sheet's jump σ (nx·n)
    end
end

@testset "parameter-square rules" begin
    tiers = SquareTiers(; levels=6)
    mono2(a, b) = 1 / ((a + 1) * (b + 1))
    for m in (2, 3, 5)
        r = square_gauss(m)
        @test all(integrate_square(r, p -> p[1]^a * p[2]^b) ≈ mono2(a, b) for a in 0:2m-1 for b in 0:2m-1)
    end
    for r in (tiers.self, tiers.edge..., tiers.vertex...)
        @test sum(r.w) ≈ 1
        @test all(p -> all(0 .<= p .<= 1), r.p)
        @test integrate_square(r, p -> p[1]^3 * p[2]^2) ≈ mono2(3, 2) rtol = 1e-12
    end
    # log singularities along side k (from corner k to k + 1), at corner k, and along all four sides
    dist_side = (p -> p[2], p -> 1 - p[1], p -> 1 - p[2], p -> p[1])
    for k in 1:4
        @test integrate_square(tiers.edge[k], p -> log(dist_side[k](p))) ≈ -1 rtol = 1e-6
        @test integrate_square(tiers.vertex[k], p -> log(sum(abs2, p - QEL.SQUARE[k]))) ≈ log(2) - 3 + π / 2 rtol = 1e-6
    end
    @test integrate_square(tiers.self, p -> sum(log(f(p)) for f in dist_side)) ≈ -4 rtol = 1e-6
end

@testset "wide sheared corners and contact" begin
    opts = AssemblyOptions()
    tiers = square_tiers(opts)
    mono2(a, b) = 1 / ((a + 1) * (b + 1))
    # a corner sheared to 150° with unequal sides: the stored rule's diagonal leaves it a 135° triangle corner, which its
    # angular coordinate under-samples; the split into physical angles of at most 75° does not
    a4 = 0.4 * SVector(cosd(150), sind(150), 0.0)
    Q = IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.1, 0.0, 0.0) + a4, a4)
    r = vertex_square(Q, 1, tiers, opts)
    @test sum(r.w) ≈ 1 && all(p -> all(-1e-14 .<= p .<= 1 + 1e-14), r.p)
    @test integrate_square(r, p -> p[1]^3 * p[2]^2) ≈ mono2(3, 2) rtol = 1e-12
    f(p) = log(norm(point(Q, p) - Q.x[1])) * norm(area_normal(Q, p)) # ∫ log|x - x₁| dA
    ref = first(hcubature(u -> f(SVector(u[1]^3, u[2]^3)) * 9 * u[1]^2 * u[2]^2, SA[0.0, 0.0], SA[1.0, 1.0]; rtol=1e-12, maxevals=10^7))
    @test abs(integrate_square(r, f) / ref - 1) < 1e-6
    @test abs(integrate_square(tiers.vertex[1], f) / ref - 1) > 1e-5
    @test vertex_square(Q, 2, tiers, opts) === tiers.vertex[2] # an acute corner: the stored rule
    # a T-junction: S covers half of T's side 1, folded by 30°; T is cut at the hanging node, the half along S takes
    # the edge rule and the other half the corner rule toward the hanging node
    T = IsoQuad(SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0), SVector(0.0, 1.0, 0.0))
    S = IsoQuad(SVector(0.5, 0.0, 0.0), SVector(0.0, 0.0, 0.0), SVector(0.0, -0.5cosd(30), 0.5sind(30)), SVector(0.5, -0.5cosd(30), 0.5sind(30)))
    c = square_contact_rule(T, S, (true, false, false, false), tiers, opts)
    @test !isnothing(c) && sum(c.w) ≈ 1
    @test integrate_square(c, p -> p[1]^3 * p[2]^2) ≈ mono2(3, 2) rtol = 1e-12
    @test count(p -> p[1] < 0.5, c.p) == length(tiers.edge[1])
    # neighbours sharing a whole side or just a corner need no contact rule
    @test isnothing(square_contact_rule(T, IsoQuad(T.x[2], T.x[1], SVector(0.0, -1.0, 0.2), SVector(1.0, -1.0, 0.2)), (true, true, false, false), tiers, opts))
end
