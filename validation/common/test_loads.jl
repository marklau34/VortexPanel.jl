# Unit tests of the load post-processing shared by the baselines (common/loads.jl).
#
#   julia --project=validation validation/common/test_loads.jl

using Test, LinearAlgebra, StaticArrays

include(joinpath(@__DIR__, "loads.jl"))

@testset "loads" begin
    ref = WingReference(span=2.0, chord=0.5, rho=1.2, vinf=10.0)
    AR = ref.span^2 / ref.area

    @testset "Trefftz drag of an elliptic loading" begin
        ratio(n) = begin
            y = [-cos(π * k / n) for k in 0:n] .* (ref.span / 2)
            yc = (y[1:end-1] .+ y[2:end]) ./ 2
            Γ = 1.3 .* sqrt.(1 .- (2 .* yc ./ ref.span) .^ 2)
            CL = kj_lift(Γ, diff(y), ref) / (ref.qinf * ref.area)
            CD = trefftz_drag(y, Γ, ref) / (ref.qinf * ref.area)
            CD / (CL^2 / (π * AR))
        end
        @test abs(ratio(320) - 1) < 0.005
        @test abs(ratio(320) - 1) < abs(ratio(80) - 1) < abs(ratio(20) - 1) # converges
    end

    @testset "Kutta-Joukowski lift" begin
        @test kj_lift(fill(2.0, 4), fill(0.5, 4), ref) ≈ ref.rho * ref.vinf * 2.0 * 2.0
    end

    @testset "Richardson extrapolation" begin
        p, ext = richardson([3.0 + 1 / 4, 3.0 + 1 / 16, 3.0 + 1 / 64])
        @test p ≈ 2 && ext ≈ 3
        @test all(isnan, richardson([1.0, 2.0, 1.5]))
    end

    @testset "pitch from a fitted rotation" begin
        α = deg2rad(5)
        R = SMatrix{3,3}([cos(α) 0 sin(α); 0 1 0; -sin(α) 0 cos(α)]) # nose up: the trailing edge drops
        xb = [SVector(rand(3)...) for _ in 1:20]
        @test fit_rotation(xb, [R * p for p in xb]) ≈ R
        @test pitch_angle(R) ≈ α
    end

    @testset "moment transfer" begin
        F, M = SA[1.0, 0.0, 10.0], SA[0.0, 2.0, 0.0]
        c = force_coefficients(F, M, ref; r_c4=SA[0.125, 0.0, 0.0])
        qSc = ref.qinf * ref.area * ref.chord
        @test c.Cm_LE ≈ 2 / qSc
        @test c.Cm_c4 ≈ (2 - (0.0 * 1.0 - 0.125 * 10.0)) / qSc # M - r × F, y component: -(z Fx - x Fz)
        @test c.CL ≈ 10 / (ref.qinf * ref.area)
    end

    @testset "structured strips" begin
        @test [structured_strip(e, 3, 2) for e in 1:13] == [1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 0]
        @test structured_chord(8, 3, 2) == (2, -1)
    end
end
