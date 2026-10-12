# The contact tier: test and source triangles that touch beyond their shared nodes (T-junctions, gaps, slivers),
# against adaptive cubature over the test triangle.

using HCubature
using VortexPanel.Elements
using VortexPanel.Assembly: contact_rule, tier_rules, pair_block, near_rule

"""The pair block ∫ψ_a n_T·(γ_b × w_S) dA by adaptive cubature: (u, v) ↦ λ = (1-u, u(1-v), uv), dA = 2A u du dv"""
function block_cubature(T, S; rtol=1e-10)
    f(uv) = begin
        u, v = uv
        λ = SVector(1 - u, u * (1 - v), u * v)
        wn = cross(source_velocity(S.src, point(T, λ)), T.n)
        vec(λ * SVector(dot(wn, S.γ[1]), dot(wn, S.γ[2]), dot(wn, S.γ[3]))') .* (2 * T.area * u)
    end
    return reshape(first(hcubature(f, (0.0, 0.0), (1.0, 1.0); rtol=rtol, maxevals=10^8)), 3, 3)
end

shared_corners(T, S) = (T.x[1] in S.x, T.x[2] in S.x, T.x[3] in S.x)
relative_error(T, S, rule, ref) = norm(first(pair_block(T, S, rule, 0.0, :hat, 0, Inf)) - ref) / norm(ref)

const OPTS = AssemblyOptions()
const RULES = tier_rules(OPTS)

# a T-junction: T's long edge a-b with the hanging node h at its middle; S = (h, a, c) across it, folded by 20°
const TA, TB, TH = SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.5, 0.0, 0.0)
const TLONG = FlatTriangle(TA, TB, SVector(0.45, -0.8, 0.0))
const TSHORT = FlatTriangle(TH, TA, SVector(0.2, 0.7cosd(20), 0.7sind(20)))

@testset "pieces" begin
    r = contact_rule(TLONG, TSHORT, shared_corners(TLONG, TSHORT), RULES, OPTS)
    # split at h: the piece along a-h carries the edge rule, the piece h-b-c the corner rule toward h
    @test length(r) == length(RULES.edge[1]) + length(RULES.vertex[1])
    @test sum(r.w) ≈ 1
    @test all(l -> all(>=(-1e-15), l) && sum(l) ≈ 1, r.λ)
    for (p, q, s) in ((1, 0, 0), (2, 1, 0), (1, 1, 1), (3, 0, 2))
        @test integrate(r, λ -> λ[1]^p * λ[2]^q * λ[3]^s) ≈ mono(p, q, s) rtol = 1e-12
    end
    # the reverse pair: S's whole edge h-a lies on T's long edge, so S takes the edge rule unsplit
    @test length(contact_rule(TSHORT, TLONG, shared_corners(TSHORT, TLONG), RULES, OPTS)) == length(RULES.edge[1])
    # a point contact: S touches T's long edge at the hanging node only
    point_contact = FlatTriangle(TH, SVector(0.3, 0.5, 0.1), SVector(0.7, 0.5, 0.1))
    @test length(contact_rule(TLONG, point_contact, shared_corners(TLONG, point_contact), RULES, OPTS)) == 2 * length(RULES.vertex[1])
    # regular neighbours (a shared edge, a shared corner) and distant triangles need no contact rule
    for S in (FlatTriangle(TB, TA, SVector(0.4, 0.7, 0.2)), FlatTriangle(TB, SVector(1.6, 0.4, 0.0), SVector(1.2, 0.8, 0.1)),
        FlatTriangle(SVector(0.0, 0.0, 2.0), SVector(1.0, 0.0, 2.0), SVector(0.5, 1.0, 2.0)))
        @test isnothing(contact_rule(TLONG, S, shared_corners(TLONG, S), RULES, OPTS))
    end
end

@testset "accuracy" begin
    for (T, S) in ((TLONG, TSHORT), (TSHORT, TLONG))
        ref = block_cubature(T, S)
        @test relative_error(T, S, contact_rule(T, S, shared_corners(T, S), RULES, OPTS), ref) < 1e-6
        # the rule of the shared corner alone misses the singular segment
        @test relative_error(T, S, RULES.vertex[findfirst(x -> x in S.x, T.x)], ref) > 1e-3
    end
    # a gap: S faces T's edge a-b across a slit of 1% and 0.1% of its length, sharing no node. The graded edge rule
    # resolves the near singularity whatever the width; the distance tiers alone would take the near rule
    for gap in (0.01, 0.001)
        S = FlatTriangle(TB + SVector(0, gap, 0), TA + SVector(0, gap, 0), SVector(0.5, 0.7, 0.3))
        ref = block_cubature(TLONG, S)
        @test relative_error(TLONG, S, contact_rule(TLONG, S, shared_corners(TLONG, S), RULES, OPTS), ref) < 1e-6
        @test relative_error(TLONG, S, near_rule(TLONG, S, RULES.near, OPTS.near_ratio, OPTS.near_depth), ref) > 1e-4
    end
end
