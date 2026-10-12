using VortexPanel.Elements
using VortexPanel.Kernels

"""∫_T λ₁^a λ₂^b λ₃^c dA / A"""
mono(a, b, c) = 2 * factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 2)
integrate(rule, f) = sum(w * f(λ) for (λ, w) in zip(rule.λ, rule.w))

@testset "Gauss-Legendre" begin
    for n in (1, 2, 5, 9)
        x, w = gauss_legendre(n)
        @test all(k -> isapprox(sum(w .* x .^ k), 1 / (k + 1); atol=1e-14), 0:2n-1)
    end
    x, w = graded_gauss(6, 0.15, 3)
    @test sum(w) ≈ 1
    @test sum(w .* log.(x)) ≈ -1 rtol = 1e-6 # ∫₀¹ ln t dt
end

@testset "triangle rules" begin
    for (n, degree) in ((1, 1), (3, 2), (6, 4), (7, 5))
        r = triangle_rule(n)
        @test sum(r.w) ≈ 1
        for a in 0:degree, b in 0:degree-a, c in 0:degree-a-b
            @test integrate(r, λ -> λ[1]^a * λ[2]^b * λ[3]^c) ≈ mono(a, b, c) atol = 1e-14
        end
    end
end

@testset "graded rules" begin
    e, v, s = edge_graded_rule(6, 0.15, 3), vertex_graded_rule(6, 0.15, 3), self_rule(6, 0.15, 3)
    for r in (e, v, s)
        @test sum(r.w) ≈ 1
        @test all(l -> all(>(0), l) && sum(l) ≈ 1, r.λ) # interior points only
        @test integrate(r, λ -> λ[1]^2 * λ[2] * λ[3]) ≈ mono(2, 1, 1) rtol = 1e-12
    end
    # log singularities: along the edge between corners 1 and 2 (λ₃ → 0), at corner 1, along all three edges
    @test integrate(e, λ -> log(λ[3])) ≈ -3 / 2 rtol = 1e-6
    @test integrate(v, λ -> log(1 - λ[1])) ≈ -1 / 2 rtol = 1e-6
    @test integrate(s, λ -> log(λ[1]) + log(λ[2]) + log(λ[3])) ≈ -9 / 2 rtol = 1e-6
    # more levels converge
    errs = [abs(integrate(edge_graded_rule(L, 0.15, 3), λ -> log(λ[3])) + 3 / 2) for L in (2, 4, 6)]
    @test errs[3] < errs[2] < errs[1]
    # rotated to edge k (corners k, k+1) the singular coordinate is the third corner's
    for k in 1:3
        @test integrate(rotate_rule(e, k), λ -> log(λ[mod1(k + 2, 3)])) ≈ -3 / 2 rtol = 1e-6
        @test integrate(rotate_rule(v, k), λ -> log(1 - λ[k])) ≈ -1 / 2 rtol = 1e-6
    end
end

@testset "flat triangle sign convention" begin
    a, b, c = SVector(0.3, -0.2, 0.1), SVector(1.4, 0.1, -0.2), SVector(0.5, 0.9, 0.4)
    el = FlatTriangle(a, b, c)
    n = el.n
    @test n ≈ normalize(cross(b - a, c - a))
    # γ_b = n × ∇λ_b = -e_b/(2A), e_b the edge opposite b (counterclockwise)
    @test el.γ[1] ≈ -(c - b) / (2el.area) && el.γ[2] ≈ -(a - c) / (2el.area) && el.γ[3] ≈ -(b - a) / (2el.area)
    @test norm(sum(el.γ)) < 1e-14 && norm(sum(el.gradλ)) < 1e-14 # a uniform doublet sheds nothing
    # a linear doublet μ(x) = k·x: ∇μ is the tangential part of k, γ = n × ∇μ, and γ × n = ∇μ
    k = SVector(0.7, -1.3, 0.4)
    μ = SVector(dot(k, a), dot(k, b), dot(k, c))
    ∇μ = k - dot(k, n) * n
    @test gradient(el, μ) ≈ ∇μ
    @test sheet_strength(el, μ) ≈ cross(n, ∇μ)
    @test cross(sheet_strength(el, μ), n) ≈ ∇μ
    # the sheet's tangential velocity jumps by ∇μ across it, as a doublet layer's potential jumps by μ
    jump = sheet_velocity(el, μ, el.centroid; side=1) - sheet_velocity(el, μ, el.centroid; side=-1)
    @test jump ≈ ∇μ
    @test basis(el, SVector(0.2, 0.3, 0.5)) == SVector(0.2, 0.3, 0.5)
    @test point(el, SVector(0.0, 1.0, 0.0)) == b
    @test distance_to(el, el.centroid + 0.3n) ≈ 0.3
end
