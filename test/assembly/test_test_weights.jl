# Validation stage 8 in brief: the test weights on the reference wing. The full study is validation/test_weights
# (refinement, sectional accuracy against the exact section, Cp, abused meshes, zero incidence).

solve_test(mesh, test; kw...) = solve_wing(mesh; assembly=AssemblyOptions(test=test), kw...)

const W1 = reference_wing(1; split_te=false)
const HAT1 = solve_test(W1, :hat)

@testset "dual cell" begin
    d = solve_test(W1, :dual)
    @test d.columns < 1e-4 && abs(d.λ) < 1e-8 # the dual cells partition the surface, as the hat functions do
    @test abs(d.CL / HAT1.CL - 1) < 0.01 && abs(d.CL_KJ / HAT1.CL_KJ - 1) < 0.003 # the same circulation, to 0.2%
    @test d.cond < 2HAT1.cond
    @test abs(solve_test(W1, :dual; α=0.0).CL) < 5e-5
end

@testset "one-ring is unstable" begin
    # one on every triangle of a node: neighbouring rows nearly coincide on thin elements, and the system degenerates
    nspan, nchord = reference_wing_size(1)
    slivers = first(sliver_stations(W1, 2, 0.02; nspan=nspan, nchord=nchord, rng=StableRNG(1)))
    h, o = solve_test(slivers, :hat), solve_test(slivers, :onering)
    @test o.cond > 50h.cond && abs(o.CL / h.CL - 1) > 1
    @test o.columns < 1e-2 # the quadrature is not to blame
end
