# The straight vortex filament of linearly varying strength, against adaptive quadrature of Biot–Savart

using QuadGK, ForwardDiff
using VortexPanel.Kernels: segment_velocity

"""Biot–Savart of the filament from a to b with end circulations Γa and Γb, by adaptive quadrature"""
function segment_quadgk(a, b, x, Γa, Γb)
    e = b - a
    f(s) = (Γa * (1 - s) + Γb * s) * cross(e, x - (a + s * e)) / norm(x - (a + s * e))^3
    return first(quadgk(f, 0.0, 1.0; rtol=1e-12)) / 4π
end

@testset "linear strength" begin
    a, b = SVector(0.2, -0.3, 0.1), SVector(1.1, 0.4, -0.2)
    off = normalize(cross(b - a, SVector(0.0, 0.0, 1.0)))
    for x in (SVector(0.5, 0.8, 0.3), SVector(-1.0, 0.2, 0.0), (a + b) / 2 + 1e-3 * off, a + 1e-3 * off, SVector(30.0, -20.0, 10.0))
        va, vb = segment_velocity(a, b, x)
        @test va ≈ segment_quadgk(a, b, x, 1.0, 0.0) rtol = 1e-9
        @test vb ≈ segment_quadgk(a, b, x, 0.0, 1.0) rtol = 1e-9
    end
end

@testset "constant strength" begin
    # the sum is the constant-strength segment; by the middle of a long one, the infinite line's Γ/(2πd), right-handed
    L, d = 1e4, 0.3
    va, vb = segment_velocity(SVector(-L, 0.0, 0.0), SVector(L, 0.0, 0.0), SVector(0.0, d, 0.0))
    @test va + vb ≈ SVector(0.0, 0.0, 1 / (2π * d)) rtol = 1e-6
    # solenoidal (the constant-strength field), zero on the extension, zero by convention on the segment
    a, b = SVector(0.2, -0.3, 0.1), SVector(1.1, 0.4, -0.2)
    J = ForwardDiff.jacobian(x -> sum(segment_velocity(a, b, x)), SVector(0.5, 0.8, 0.3))
    @test abs(tr(J)) < 1e-12 * norm(J)
    @test norm(sum(segment_velocity(a, b, a + 2(b - a)))) < 1e-14
    @test all(iszero, segment_velocity(a, b, (a + b) / 2)[1])
    # generic over the number type
    @test eltype(first(segment_velocity(Float32.(a), Float32.(b), SVector(0.5f0, 0.8f0, 0.3f0)))) == Float32
end
