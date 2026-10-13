# Free wake with vortex particles: the trailing edge panel with explicit rows reproduces the implicit rigid wake, the
# right hand side of explicit elements equals assembled columns times their doublets, conversion conserves the
# vorticity and its first moment, and short runs keep Kelvin's theorem and stay bounded.

using VortexPanel.FreeWake
using VortexPanel.Elements: FlatTriangle, sheet_strength
const FWK = VortexPanel.FreeWake

const FW_WING = reference_wing(1; split_te=false)
const FW_V, FW_LIFT = wing_freestream(5.0)
const FW_G = build_geometry(FW_WING)
const FW_DM = build_dofmap(FW_G)
const FW_ELS = build_elements(FW_G)
const FW_SIGMA = source_strengths(FW_ELS, FW_V)

@testset "sheet vorticity" begin
    p, q, r, μ = SVector(0.1, 0.0, 0.2), SVector(1.0, 0.3, -0.1), SVector(0.4, 0.9, 0.3), SVector(0.7, -0.2, 1.1)
    el = FlatTriangle(p, q, r)
    @test FWK.sheet_vorticity(p, q, r, μ) ≈ sheet_strength(el, μ) * el.area # ∮μ dl = γA
end

@testset "trailing edge panel with explicit rows" begin
    # the rigid wake solved implicitly ...
    wf = build_wake(FW_G, FW_DM, FW_V)
    A, b, _ = assemble(FW_G, FW_DM, FW_ELS, FW_SIGMA, FW_V; wake=wf)
    ref = solve_galerkin(FW_DM, A, b; W=wake_rows(FW_DM, wf.te))
    μw = FW_DM.K * ref.q
    # ... equals its first row as the implicit trailing edge panel (factored once) plus the other rows as explicit
    # elements frozen at μ_w (a free wake state on the same nodes; μ_line = μ_w leaves no far end filament)
    te = te_panel(FW_G, FW_DM, FW_V; length=wf.rows[2])
    At, bt, _ = assemble(FW_G, FW_DM, FW_ELS, FW_SIGMA, FW_V; wake=te)
    F = factorize_galerkin(FW_DM, At; W=wake_rows(FW_DM, te.te))
    topo = FWK.Topology(FW_DM)
    n = topo.nte
    st = FWK.FreeWakeState{Float64}(n)
    for j in 1:length(wf.rows)-1
        push!(st.x, wf.x[j*n+1:(j+1)*n])
        push!(st.μ, copy(μw))
    end
    st.μ_line .= μw
    we = vortex_elements(st, topo)
    no_blobs = VortexPanel.Kernels.Blob{Float64}[]
    rhs = explicit_rhs(FW_ELS, FW_DM, FWK.ExplicitSources(we, no_blobs))
    s = solve_galerkin(F, bt + rhs)
    @test maximum(abs, FW_DM.K * s.q - μw) < 1e-12 * maximum(abs, μw)
    @test maximum(abs, s.q - ref.q) < 1e-12 * maximum(abs, ref.q)
    # explicit_rhs is the assembled columns times the doublets: the trailing edge panel's own columns
    cols = At[:, FW_DM.nslot+1:end] * μw[te.te]
    @test maximum(abs, explicit_rhs(FW_ELS, FW_DM, FWK.ExplicitSources(vortex_elements(te, μw), no_blobs)) + cols) < 1e-13 * maximum(abs, cols)
end

@testset "conversion conserves the vorticity and its first moment" begin
    for (tip, max_rows) in ((:zero, 3), (:free, 1))
        dm = build_dofmap(FW_G; opts=DOFMapOptions(tip_kutta=tip))
        sim = FreeWakeSim(FW_G, dm, FW_V; opts=FreeWakeOptions(dt=0.06, max_rows=10_000, advection=:free))
        run!(sim, 6)
        st = sim.state
        core = h -> FWK.core_of(sim, h)
        before = [FWK.free_wake_blobs(vortex_elements(st, sim.topo), core); st.particles]
        Ω0, I0 = vorticity_moments(before)
        scale_Ω, scale_I = sum(norm(b.α) for b in before), sum(norm(b.x) * norm(b.α) for b in before)
        for k in 1:3 # the first conversion (with the starting vortex) and later ones
            convert_oldest!(st, sim.topo, core)
            Ω, I = vorticity_moments([FWK.free_wake_blobs(vortex_elements(st, sim.topo), core); st.particles])
            @test norm(Ω - Ω0) < 1e-14 * scale_Ω
            @test norm(I - I0) < 1e-14 * scale_I
        end
        @test st.conversions == 3 && length(st.x) == 4
    end
end

@testset "impulsive start" begin
    dt = VortexPanel.Wake.te_panel_size(FW_G)
    # prescribed wake: the rows sit exactly V∞ dt apart, Kelvin's theorem holds and the force routes agree
    sim = FreeWakeSim(FW_G, FW_DM, FW_V; opts=FreeWakeOptions(dt=dt, max_rows=4, advection=:freestream))
    hist = run!(sim, 20)
    st = sim.state
    @test all(k -> st.x[k+1][3] ≈ sim.edge[3] + k * dt * FW_V, 0:length(st.x)-2)
    @test all(r -> norm(r.Ω) < 1e-12, hist)
    S = WING_SPAN * WING_CHORD
    cl(r, f) = dot(getproperty(r, f), FW_LIFT) / S
    kj = [cl(r, :F_kj) for r in hist]
    @test kj[1] < kj[5] < kj[10] < kj[20] # the circulation builds up
    @test abs(cl(hist[end], :F_pressure) - cl(hist[end], :F_impulse)) < 0.02 * cl(hist[end], :F_pressure)
    @test length(st.particles) > 0 && st.conversions == 20 - 4 # a conversion every step once 4 rows exist
    # free wake: bounded velocities and a lift close to the prescribed wake's
    free = FreeWakeSim(FW_G, FW_DM, FW_V; opts=FreeWakeOptions(dt=dt, max_rows=4))
    fh = run!(free, 20)
    @test all(r -> r.umax < 3 * norm(FW_V), fh)
    @test abs(cl(fh[end], :F_kj) / kj[end] - 1) < 0.03
end
