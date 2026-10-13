# Validation stage 5 in brief: the reference wing with a rigid wake, by both force routes. The full study is
# validation/wing (convergence to CL ≈ 0.340, abused meshes, tip treatments, wake options, MESA's baseline).

using VortexPanel.DOFMap: nwake_te

const WING1 = reference_wing(1; split_te=false)
const R1 = solve_wing(WING1)

@testset "lift and induced drag" begin
    @test 0.33 < R1.CL < 0.36 # level 1; converges to 0.340
    @test abs(R1.CL_KJ / R1.CL - 1) < 0.02 # the two lift routes
    @test abs(R1.CDi_T / R1.CDi_p - 1) < 0.2 # the two induced drag routes, converging under refinement
    @test abs(R1.CY) < 1e-5
    @test norm(R1.μ_w - reverse(R1.μ_w)) < 0.01 * norm(R1.μ_w) # symmetric loading
    @test 0 < R1.tip < 0.3 # the tip's share of the loading falls under refinement
end

@testset "Kutta condition and system" begin
    # μ_w = μ_u - μ_l at every trailing edge node (a matching trailing edge: the lower copy is the second slot)
    wte = R1.dm.wake_te
    @test all(R1.μ_w[k] ≈ R1.sol.μ[v] - R1.sol.μ[R1.dm.lower_slot[v]] for (k, v) in enumerate(wte.node))
    @test R1.columns < 1e-4 && abs(R1.λ) < 1e-8
    # zero lift at zero incidence; the caps' corners at the trailing edge's end nodes take the upper copy, which
    # breaks the up-down symmetry at the tips a little
    r0 = solve_wing(WING1; α=0.0)
    @test abs(r0.CL) < 1e-3 && maximum(abs, r0.μ_w) < 5e-3
end

@testset "wake options" begin
    for o in (WakeOptions(ending=:semi_infinite), WakeOptions(direction=:bisector), WakeOptions(length=50.0), WakeOptions(first_row=0.045))
        @test abs(solve_wing(WING1; wake=o).CL / R1.CL - 1) < 1e-3
    end
end

@testset "abused wings" begin
    nspan, nchord = reference_wing_size(1)
    tj = solve_wing(first(insert_tjunctions(WING1, 4; rng=StableRNG(1))))
    @test tj.g.quality.counts[:tjunctions] == 8 && abs(tj.CL / R1.CL - 1) < 1e-3
    # a non-matching trailing edge sheds its wake from the union of both edges' nodes
    nm = solve_wing(first(skew_wing(WING1, :te_misaligned, 0.3; nspan=nspan, nchord=nchord, rng=StableRNG(1), shared_te=false)))
    @test nm.g.quality.counts[:kutta_pairs] == 1 && nwake_te(nm.dm.wake_te) > nwake_te(R1.dm.wake_te)
    @test abs(nm.CL / R1.CL - 1) < 0.005
end
