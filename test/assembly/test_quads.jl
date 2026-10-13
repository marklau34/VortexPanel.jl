# Validation stage 7 in brief: the isoparametric quads (`MixedMesh`) against the split triangles. The full study is
# validation/elements (cube-sphere, wing refinement, sectional accuracy, skew, misalignment, abused meshes).

using HCubature
using VortexPanel.Elements
using VortexPanel.Assembly: quad_pair_rule, square_tiers

"""
The pair block of source quad S on test quad T (hat weights; the 4×4 sheet block and the source part of the right hand
side as one vector) by adaptive cubature over T's parameter square, in pieces: each map takes u ∈ [lo, 1]² to a
point graded toward the pair's singular corner or side, with its Jacobian. With `self`, T = S and the one-sided
limits on the side the normal points to.
"""
function quad_block_cubature(T, S, V, maps; self::Bool=false, rtol=1e-8, lo=1e-3)
    total = zero(SVector{20,Float64})
    for m in maps
        f(u) = begin
            p, jac = m(u)
            nJ = area_normal(T, p)
            J = norm(nJ)
            x = point(T, p)
            c, s = self ? VortexPanel.Elements.polar_field(S, x, nJ / J, V, p, 0.0, 1) : quad_field(S, x, nJ / J, V)
            ψ = shape(p[1], p[2]) * (J * jac)
            SVector{20}(vec(ψ * c')..., (-s * ψ)...)
        end
        total += first(hcubature(f, SA[lo, lo], SA[1.0, 1.0]; rtol=rtol, maxevals=3_000_000))
    end
    return total
end

"""The assembly's pair block with the rule, as `quad_block_cubature`'s vector"""
quad_block_vector(T, S, rule, V, side) = (b = VortexPanel.Assembly.block(T, S, rule, V, :hat, side); SVector{20}(vec(b[1])..., b[2]...))

const VQ = SVector(0.8, -0.3, 0.5)
const QOPTS = AssemblyOptions()
const QTIERS = square_tiers(QOPTS)

@testset "touching quad pairs against cubature" begin
    x1, x2 = SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)
    # a shared side (T's side 1, singular at its two ends too), the neighbour folded by 25° and warped
    T = IsoQuad(x1, x2, SVector(0.9, 0.8, 0.0), SVector(0.1, 0.7, 0.0))
    S = IsoQuad(x2, x1, x1 + SVector(0.05, -0.6cosd(25), 0.6sind(25)), x2 + SVector(-0.05, -0.65cosd(25), 0.65sind(25) + 0.05))
    tier, rule = quad_pair_rule(T, S, SVector(1, 2, 3, 4), SVector(2, 1, 5, 6), false, QTIERS, QOPTS, false)
    @test tier == 2
    halves = (u -> (SVector(0.5u[1]^3, u[2]^3), 4.5 * u[1]^2 * u[2]^2), u -> (SVector(1 - 0.5u[1]^3, u[2]^3), 4.5 * u[1]^2 * u[2]^2))
    ref = quad_block_cubature(T, S, VQ, halves)
    @test norm(quad_block_vector(T, S, rule, VQ, 0) - ref) < 1e-5 * norm(ref)
    # a shared corner sheared to 150° with unequal sides: split by physical angle, it is 500× more accurate than the
    # stored rule
    a4 = 0.4 * SVector(cosd(150), sind(150), 0.0)
    T = IsoQuad(x1, x2, SVector(1.1, 0.0, 0.0) + a4, a4)
    S = IsoQuad(x1, SVector(-0.3, -0.7, 0.25), SVector(0.4, -0.9, 0.3), SVector(0.7, -0.4, 0.1))
    tier, rule = quad_pair_rule(T, S, SVector(1, 2, 3, 4), SVector(1, 7, 8, 9), false, QTIERS, QOPTS, false)
    @test tier == 3
    ref = quad_block_cubature(T, S, VQ, (u -> (SVector(u[1]^3, u[2]^3), 9 * u[1]^2 * u[2]^2),))
    @test norm(quad_block_vector(T, S, rule, VQ, 0) - ref) < 1e-6 * norm(ref)
    @test norm(quad_block_vector(T, S, QTIERS.vertex[1], VQ, 0) - ref) > 1e-5 * norm(ref)
    # the self pair of a warped quad (principal values; the source's one-sided limit), graded toward all four sides
    W = IsoQuad(x1, SVector(1.0, 0.1, 0.2), SVector(1.2, 1.1, -0.1), SVector(-0.1, 0.9, 0.15))
    quarters = [u -> (SVector(i == 0 ? 0.5u[1]^3 : 1 - 0.5u[1]^3, j == 0 ? 0.5u[2]^3 : 1 - 0.5u[2]^3), 2.25 * u[1]^2 * u[2]^2) for i in 0:1 for j in 0:1]
    ref = quad_block_cubature(W, W, VQ, quarters; self=true, rtol=1e-7)
    @test norm(quad_block_vector(W, W, QTIERS.self, VQ, 1) - ref) < 2e-6 * norm(ref)
end

@testset "cube-sphere" begin
    flow = EllipsoidFlow((1.0, 1.0, 1.0), normalize(SVector(1.0, 0.3, -0.2)))
    V = flow.V∞
    g = build_geometry(cube_sphere(6))
    dm = build_dofmap(g)
    m = MixedMesh(g, dm)
    @test length(m.quads) == 216 && isempty(m.tris)
    Aq, bq, _ = assemble(g, dm, m, V)
    els = build_elements(g)
    At, bt, _ = assemble(g, dm, els, source_strengths(els, V), V)
    exact = [exact_doublet(flow, g.x[dm.slot_node[s]]) for s in 1:dm.nslot]
    function μ_error(A, b)
        s = solve_galerkin(dm, A, b)
        e = s.μ - exact
        e .-= sum(dm.slot_area .* e) / sum(dm.slot_area)
        return sqrt(sum(dm.slot_area .* e .^ 2) / sum(dm.slot_area)) / doublet_scale(flow), s
    end
    eq, sol = μ_error(Aq, bq)
    et, _ = μ_error(At, bt)
    @test eq < 0.7et # 0.54× the triangles' error
    @test consistency(dm, Aq, bq).columns < 1e-5 && abs(sol.λ[1]) < 1e-12
    # Cp at the quads' centres, where the flat quad's normal is the sphere's to second order (0.025 against the
    # triangles' 0.11 at their centroids)
    pts = surface_pressure(m, sol.μ, V; m=1)
    @test maximum(abs, pts.cp - [exact_cp(flow, x) for x in pts.x]) < 0.03
end

@testset "two-dimensional limit" begin
    # Every quad split along the same diagonal: for a spanwise-uniform μ the triangles' sheets and summed test functions
    # are the bilinear quad's, so at the midspan of a slender wing the two elements agree. The default split mirrors its
    # diagonals at y = 0, which dents the circulation there (0.6% here).
    b = 4WING_SPAN
    mesh = reference_wing(1; split_te=false, nspan=18, nchord=9, span=b)
    q, d, t = (midspan_section(solve_wing(mesh; span=b, kw...)) for kw in ((element=:quad,), (geometry=GeometryOptions(split=:d13),), NamedTuple()))
    @test abs(q.cl_kj / d.cl_kj - 1) < 5e-4 && abs(q.cl_p / d.cl_p - 1) < 5e-4 && abs(q.cm - d.cm) < 1e-4
    @test t.cl_kj / q.cl_kj - 1 < -3e-3
end

@testset "wing" begin
    w = reference_wing(1; split_te=false)
    t, q = solve_wing(w), solve_wing(w; element=:quad)
    @test abs(q.CL / t.CL - 1) < 1e-3 && abs(q.CL_KJ / t.CL_KJ - 1) < 1e-3 && abs(q.CDi_T / t.CDi_T - 1) < 2e-3
    @test q.columns < 5e-5 && q.cond < t.cond
    @test abs(solve_wing(w; α=0.0, element=:quad).CL) < 2e-5
    # nothing to flip: offsetting the rows aft of 80% chord (shared trailing edge) flips 34 diagonals of the split,
    # which costs the triangles 0.59% of their lift, the quads 0.08%
    nspan, nchord = reference_wing_size(1)
    mis = first(skew_wing(w, :te_misaligned, 0.3; nspan=nspan, nchord=nchord, shared_te=true))
    dq, dt = solve_wing(mis; element=:quad).CL / q.CL - 1, solve_wing(mis).CL / t.CL - 1
    @test abs(dq) < 1.5e-3 && abs(dt) > 3abs(dq)
end
