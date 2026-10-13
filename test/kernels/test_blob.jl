# The vortex particle (Rosenhead–Moore blob) of the free wake: the point-vortex limit, solenoidality, and the 2-point
# Gauss particles that stand for a filament.

using ForwardDiff
using VortexPanel.Kernels

@testset "blob kernel" begin
    xp, α, δ = SVector(0.1, -0.2, 0.3), SVector(0.4, 1.0, -0.7), 0.05
    point_vortex(x) = cross(α, x - xp) / (4π * norm(x - xp)^3)
    for r in (20δ, 100δ)
        x = xp + r * normalize(SVector(0.3, -0.5, 0.8))
        # (r² + δ²)^{-3/2} against r⁻³: a relative difference of 1.5 (δ/r)²
        @test norm(blob_velocity(x, xp, α, δ) - point_vortex(x)) <= 2 * (δ / r)^2 * norm(point_vortex(x))
    end
    @test blob_velocity(xp, xp, α, δ) == zero(xp) # regular at the particle
    for x in (xp + SVector(0.01, 0.02, -0.01), xp + SVector(0.3, -0.1, 0.2))
        J = ForwardDiff.jacobian(y -> blob_velocity(y, xp, α, δ), x)
        @test abs(J[1, 1] + J[2, 2] + J[3, 3]) < 1e-12 * opnorm(J) # solenoidal
    end
    @test blob_velocity(Blob(xp, α, δ), xp + SVector(1.0, 0.0, 0.0)) == blob_velocity(xp + SVector(1.0, 0.0, 0.0), xp, α, δ)
end

@testset "filament as particles" begin
    a, b, Γa, Γb = SVector(0.0, 0.0, 0.0), SVector(1.0, 0.2, -0.1), 0.7, 1.3
    pair = segment_blobs(a, b, Γa, Γb, 0.01)
    # total strength ∫Γ dl and first moment ∫y ⊗ Γ dl (quadratic along the segment) exactly
    @test sum(p.α for p in pair) ≈ (Γa + Γb) / 2 * (b - a)
    Γ(s) = (1 - s) * Γa + s * Γb
    moment = sum(w * (a + s * (b - a)) * Γ(s)' for (s, w) in zip((0.0, 0.5, 1.0), (1 / 6, 4 / 6, 1 / 6))) # Simpson: exact for quadratics
    @test sum(p.x * p.α' for p in pair) ≈ moment * (b - a)'
    # a chain of them converges to the filament's velocity away from it (constant strength: segment_velocity's sum)
    x = SVector(0.4, 0.5, 0.3)
    va, vb = VortexPanel.Kernels.segment_velocity(a, b, x)
    exact = va + vb
    err(m) = norm(sum(blob_velocity(p, x) for k in 0:m-1 for p in segment_blobs(a + k / m * (b - a), a + (k + 1) / m * (b - a), 1.0, 1.0, 1e-4)) - exact)
    @test err(16) < err(4) < err(1)
    @test err(16) < 1e-4 * norm(exact)
end
