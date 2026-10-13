# The two force routes on prescribed loadings: pressure integration, Kutta–Joukowski and the Trefftz plane

using VortexPanel.DOFMap: WakeTE
using VortexPanel.Elements: FlatTriangle

"""A straight trailing edge along +y at the given stations, as the wake sees it"""
function straight_te(y)
    n = length(y)
    x = [SVector(0.0, yk, 0.0) for yk in y]
    return WakeTE(collect(1:n), x, fill(SVector(0.0, 1.0, 0.0), n), fill(SVector(1.0, 0.0, 0.0), n), [collect(1:n)], [false])
end

@testset "pressure" begin
    # a uniform pressure on a closed surface exerts no force and no moment
    els = build_elements(build_geometry(icosphere(2)))
    F, M = pressure_force(els, fill(0.7, length(els)); origin=SVector(0.3, -0.2, 0.1))
    @test norm(F) < 1e-13 && norm(M) < 1e-13
    # a unit pressure difference on a flat plate: its area times the normal
    plate = [FlatTriangle(SVector(0.0, 0, 0), SVector(2.0, 0, 0), SVector(0.0, 3, 0))]
    @test first(pressure_force(plate, [-1.0])) ≈ SVector(0.0, 0.0, 3.0)
end

@testset "Kutta-Joukowski" begin
    # uniform circulation 0.8 over a span of 2: lift 2·0.8·2/|V∞| per q∞, along V∞ × τ
    wte = straight_te(range(-1.0, 1.0, length=9))
    V = SVector(2.0, 0.0, 0.0)
    @test kutta_joukowski_force(wte, fill(0.8, 9), V) ≈ SVector(0.0, 0.0, 2 * 0.8 * 2.0 / 2.0)
end

@testset "Trefftz plane" begin
    # an elliptic loading gives CDi = CL²/(π AR), approached as the cosine-spaced strips refine
    b, c = 2.0, 0.5
    AR = b / c
    V = SVector(1.0, 0.0, 0.0)
    function ratio(n)
        y = [-b / 2 * cos(π * k / n) for k in 0:n]
        wte = straight_te(y)
        μ = [sqrt(max(0.0, 1 - (2yk / b)^2)) for yk in y]
        CL = kutta_joukowski_force(wte, μ, V)[3] / (b * c)
        CD = trefftz_drag(wte, μ, V) / (b * c)
        return CD / (CL^2 / (π * AR))
    end
    @test abs(ratio(160) - 1) < 0.01
    @test abs(ratio(160) - 1) < abs(ratio(40) - 1) < abs(ratio(10) - 1)
    # rotating the freestream about the span leaves the drag of the projected loading unchanged
    y = [-cos(π * k / 40) for k in 0:40]
    μ = [sqrt(max(0.0, 1 - yk^2)) for yk in y]
    @test trefftz_drag(straight_te(y), μ, SVector(cosd(5), 0.0, sind(5))) ≈ trefftz_drag(straight_te(y), μ, V)
end
