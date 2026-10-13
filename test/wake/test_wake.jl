# Rigid wake strips on the reference wing

using VortexPanel.Wake
using VortexPanel.DOFMap: nwake_te

const WG = build_geometry(reference_wing(1; split_te=false))
const WDM = build_dofmap(WG)
const WV = SVector(cosd(5), 0.0, sind(5))

@testset "strips" begin
    w = build_wake(WG, WDM, WV)
    n = nwake_te(WDM.wake_te)
    nrow = length(w.rows) - 1
    @test nwake(w) == n * (nrow + 1) && length(w.tris) == 2 * (n - 1) * nrow
    @test w.span ≈ WING_SPAN
    @test w.rows[end-1] < 20w.span <= w.rows[end] # 20 spans
    @test w.rows[3] - w.rows[2] ≈ 1.2 * (w.rows[2] - w.rows[1]) # geometric growth
    # the trailing edge row sits on the body's trailing edge nodes, and the strips run along the freestream
    @test all(w.x[i] == WDM.wake_te.x[w.te[i]] for i in eachindex(w.x) if w.row[i] == 0)
    @test all(normalize(w.x[i] - WDM.wake_te.x[w.te[i]]) ≈ WV for i in eachindex(w.x) if w.row[i] > 0)
    # every triangle faces the upper side, wake_dir × τ (+z), so the trailing edge filaments cancel under μ_u - μ_l
    @test all(el -> el.n[3] > 0.99, wake_elements(w))
end

@testset "tip filaments" begin
    w = build_wake(WG, WDM, WV)
    nrow = length(w.rows) - 1
    line = WDM.wake_te.lines[1] # along τ = +y: from the tip at y = -b/2 to the tip at y = +b/2
    @test WDM.wake_te.x[line[1]][2] < 0 < WDM.wake_te.x[line[end]][2]
    @test length(w.filaments) == 2nrow
    for (a, b) in w.filaments # clockwise about the normal: upstream at y = -b/2, downstream at y = +b/2 (right-handed)
        @test w.te[a] == w.te[b] && abs(w.row[a] - w.row[b]) == 1
        @test w.te[a] == (w.row[b] > w.row[a] ? line[end] : line[1])
    end
end

@testset "options" begin
    w = build_wake(WG, WDM, WV)
    n = nwake_te(WDM.wake_te)
    wb = build_wake(WG, WDM, WV; opts=WakeOptions(direction=:bisector))
    @test all(normalize(wb.x[n+k] - wb.x[k]) ≈ WDM.wake_te.bisector[k] for k in 1:n)
    @test all(normalize(wb.x[2n+k] - wb.x[n+k]) ≈ WV for k in 1:n) # then the freestream
    ws = build_wake(WG, WDM, WV; opts=WakeOptions(ending=:semi_infinite))
    @test length(ws.rows) == length(w.rows) + 1 && ws.rows[end] ≈ 1e3 * w.rows[end]
    @test build_wake(WG, WDM, WV; opts=WakeOptions(first_row=0.045)).rows[2] == 0.045
    @test w.rows[2] ≈ VortexPanel.Wake.te_panel_size(WG) # the trailing edge panels' size
    @test_throws ArgumentError build_wake(WG, WDM, WV; opts=WakeOptions(ending=:open))
    @test_throws ArgumentError build_wake(build_geometry(icosphere(1)), build_dofmap(build_geometry(icosphere(1))), WV)
end

@testset "export" begin
    w = build_wake(WG, WDM, WV)
    files = MeshLibrary.write_vtk(w, joinpath(mktempdir(), "wing"); point_data=Dict("mu" => w.row .* 1.0))
    @test length(files) == 2 && all(isfile, files)
end
