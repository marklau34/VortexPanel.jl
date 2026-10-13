"""
    FreeWakeSim(g, dm, V∞; opts, assembly=AssemblyOptions())

A free wake run on the body of geometry g (flat triangles) at rest in the freestream V∞, started impulsively. The
trailing edge panel (`te_panel`) is assembled with the body and the bordered matrix factored once
(`factorize_galerkin`); every `step!` then changes only the right hand side. `history` holds one record per step.
"""
mutable struct FreeWakeSim{T}
    g::SurfaceGeometry{T}
    dm::DOFMapping{T}
    els::Vector{FlatTriangle{T}}
    σ::Vector{T}
    V∞::SVector{3,T}
    opts::FreeWakeOptions
    assembly::AssemblyOptions
    te::RigidWake{T}
    edge::Vector{SVector{3,T}}
    topo::Topology
    b0::Vector{T}
    F::GalerkinFactorization
    state::FreeWakeState{T}
    core::CoreRule{T}
    keep_rows::Int
    merge::Union{Nothing,Vector{Int}}
    n::Int
    t::T
    μ::Vector{T}
    μ_w::Vector{T}
    cp::Vector{T}
    I_prev::SVector{3,T}
    history::Vector{Any}
    t_setup::Float64
end

function FreeWakeSim(g::SurfaceGeometry{T}, dm::DOFMapping, V∞::SVector{3}; opts::FreeWakeOptions,
    assembly::AssemblyOptions=AssemblyOptions()) where {T}
    check_options(opts)
    V = SVector{3,T}(V∞)
    local els, σ, te, b0, F
    t_setup = @elapsed begin
        els = build_elements(g)
        σ = source_strengths(els, V)
        te = te_panel(g, dm, V; length=something(opts.te_len, 0.3 * norm(V) * opts.dt), direction=opts.te_direction)
        A, b0, _ = assemble(g, dm, els, σ, V; opts=assembly, wake=te)
        F = factorize_galerkin(dm, A; W=wake_rows(dm, te.te))
    end
    topo = Topology(dm)
    edge = [te.x[topo.nte+k] for k in 1:topo.nte] # the trailing edge panel's downstream edge (its row 1)
    core = CoreRule(opts, V)
    row = norm(V) * opts.dt
    keep = opts.max_rows === :auto ? rows_for_distance(3, core, edge, topo, te.rows[2], row) : opts.max_rows
    opts.max_rows === :auto || check_conversion_distance(keep, core, edge, topo, te.rows[2], row)
    widths = [norm(edge[B] - edge[A]) for (A, B) in topo.segments]
    merge = opts.merge_rows === :auto ? [max(1, round(Int, w / row)) for w in widths] :
            opts.merge_rows == 1 ? nothing : fill(opts.merge_rows, length(widths))
    return FreeWakeSim(g, dm, els, σ, V, opts, assembly, te, edge, topo, b0, F, FreeWakeState{T}(topo.nte, length(topo.segments)), core,
        keep, merge, 0, zero(T), zeros(T, dm.nslot), zeros(T, topo.nte), zeros(T, length(els)), zero(SVector{3,T}), Any[], t_setup)
end

Base.show(io::IO, s::FreeWakeSim) = print(io, "FreeWakeSim: step $(s.n), t = $(s.t), $(max(length(s.state.x) - 1, 0)) free rows, ",
    "$(length(s.state.particles)) particles")

"""The core of the widest strip's first particles (a triangle of its first converted quad, rows `row` long)"""
function first_core(core::CoreRule, edge, topo::Topology, row::Real)
    w = maximum(norm(edge[B] - edge[A]) for (A, B) in topo.segments)
    o = zero(edge[1])
    return core(o, o + SVector(w, 0, 0), o + SVector(w, row, 0))
end

"""
    conversion_distance(max_rows, core, edge, topo, te_len, row)

Where the first particles are created behind the trailing edge (the middle of the oldest row, `max_rows` rows of length
`row` behind the trailing edge panel), in cores of the widest strip's first particles (`first_core`). Below about two
cores the particles' smoothed field misrepresents the sheet next to the wing and biases the lift (validation: 0.25
cores +2.3%, 0.55 +0.9%, 1.5 +0.14%, 3.3 +0.02%). Returns the distance in cores and the core.
"""
function conversion_distance(max_rows::Integer, core::CoreRule, edge, topo::Topology, te_len::Real, row::Real)
    δ = first_core(core, edge, topo, row)
    return (te_len + (max_rows + 0.5) * row) / δ, δ
end

"""The fewest rows that create the first particles at least `cores` cores behind the trailing edge"""
rows_for_distance(cores::Real, core::CoreRule, edge, topo::Topology, te_len::Real, row::Real) =
    max(0, ceil(Int, (cores * first_core(core, edge, topo, row) - te_len) / row - 0.5))

"""Warn when particles would be created within two cores of the trailing edge, with the `max_rows` that avoids it"""
function check_conversion_distance(max_rows::Integer, core::CoreRule, edge, topo::Topology, te_len::Real, row::Real)
    d, δ = conversion_distance(max_rows, core, edge, topo, te_len, row)
    if d < 2
        need = rows_for_distance(2, core, edge, topo, te_len, row)
        @warn "Particles would be created $(round(d, digits=2)) cores behind the trailing edge (core $(round(δ, sigdigits=3))), which biases the lift; max_rows = $need keeps them two cores away"
    end
    return d
end

"""
    step!(sim)

Advance one time step, in MESA's order:
1. the right hand side from the free wake and the particles (`explicit_rhs`), and the solve with the factored matrix,
   giving the doublets and the trailing edge's μ_w;
2. the loads of the state just solved: pressure (unsteady Bernoulli, ∂μ/∂t by backward difference; zero at the start),
   Kutta–Joukowski from μ_w, and the impulse route -2/|V∞|² dI/dt with the total vorticity Ω and impulse I;
3. velocities at the free nodes and particles: V∞, the body's exact kernels, and every wake element as blobs (the
   trailing edge panel, the free rows and their filaments, the particles), or V∞ alone with `advection = :freestream`;
   the first step first sheds R0 so that it is advected (MESA's first step);
4. forward Euler on the free nodes (R0 included) and the particles;
5. a new R0 on the trailing edge panel's downstream edge with this step's μ_w, conversion of the oldest rows while more
   than `max_rows` remain, and removal of particles beyond `delete_beyond` spans.
Returns the step's record (also pushed to `sim.history`).
"""
function step!(sim::FreeWakeSim{T}) where {T}
    st, o, dm = sim.state, sim.opts, sim.dm
    core = sim.core
    we = vortex_elements(st, sim.topo)
    t_rhs = @elapsed rhs = explicit_rhs(sim.els, dm, ExplicitSources(we, st.particles); opts=sim.assembly)
    sol = solve_galerkin(sim.F, sim.b0 + rhs)
    μ_w = dm.K * sol.q
    loads = step_loads(sim, sol, μ_w, we, core)
    isempty(st.x) && (push!(st.x, copy(sim.edge)); push!(st.μ, copy(μ_w)))
    X = [reduce(vcat, st.x); [p.x for p in st.particles]]
    t_vel = @elapsed U = wake_velocities(sim, sol, μ_w, we, X, core)
    k = 0
    for row in st.x, j in eachindex(row)
        k += 1
        row[j] += o.dt * U[k]
    end
    for (j, p) in enumerate(st.particles)
        k += 1
        st.particles[j] = Blob(p.x + o.dt * U[k], p.α, p.δ)
    end
    pushfirst!(st.x, copy(sim.edge))
    pushfirst!(st.μ, copy(μ_w))
    while length(st.x) - 1 > sim.keep_rows
        convert_oldest!(st, sim.topo, core, sim.merge)
    end
    isfinite(o.delete_beyond) && delete_far!(sim)
    rec = merge((n=sim.n, t=sim.t), loads, (rows=length(st.x) - 1, particles=length(st.particles), λ=sol.λ[1],
        umax=maximum(norm, U; init=zero(T)), t_rhs=t_rhs, t_vel=t_vel))
    push!(sim.history, rec)
    sim.μ, sim.μ_w = sol.μ, μ_w
    sim.n += 1
    sim.t += o.dt
    return rec
end

"""Loads of the solved state (see `step!`); keeps the pressure coefficient of the last step in `sim.cp`"""
function step_loads(sim::FreeWakeSim{T}, sol, μ_w, we, core) where {T}
    dm, els, V = sim.dm, sim.els, sim.V∞
    u = gradient_velocity(els, dm, sol.μ, V)
    dμ = sim.n == 0 ? zeros(T, dm.nslot) : (sol.μ - sim.μ) / sim.opts.dt
    cp = pressure_coefficient(u, V, [sum(m) / 3 for m in element_values(dm, dμ)])
    sim.cp = cp
    Fp, Mp = pressure_force(els, cp)
    Fkj = kutta_joukowski_force(dm, μ_w, V)
    body = [Blob(el.centroid, γ * el.area, zero(T)) for (el, γ) in zip(els, sheet_strengths(els, dm, sol.μ))]
    blobs = [body; free_wake_blobs(vortex_elements(sim.te, μ_w), core); free_wake_blobs(we, core); sim.state.particles]
    Ω, I = vorticity_moments(blobs)
    Fimp = sim.n == 0 ? SVector{3,T}(NaN, NaN, NaN) : -(2 / dot(V, V)) * (I - sim.I_prev) / sim.opts.dt
    sim.I_prev = I
    return (μ_w=μ_w, F_pressure=Fp, M_pressure=Mp, F_kj=Fkj, F_impulse=Fimp, Ω=Ω, I=I)
end

"""Velocities at the points X (free nodes, particles): V∞ + body (exact kernels) + every wake element as blobs"""
function wake_velocities(sim::FreeWakeSim{T}, sol, μ_w, we, X, core) where {T}
    V = sim.V∞
    sim.opts.advection == :freestream && return fill(V, length(X))
    ub = body_velocity(sim.els, sim.σ, sheet_strengths(sim.els, sim.dm, sol.μ), X; kernel_far=sim.assembly.kernel_far)
    blobs = [free_wake_blobs(vortex_elements(sim.te, μ_w), core); free_wake_blobs(we, core); sim.state.particles]
    uw = blob_velocities(blobs, X)
    return [V + a + c for (a, c) in zip(ub, uw)]
end

"""Velocity of the blobs at every point of X, direct summation over threads"""
function blob_velocities(blobs::AbstractVector{<:Blob}, X::AbstractVector{SVector{3,T}}) where {T}
    U = similar(X)
    Threads.@threads for i in eachindex(X)
        x = X[i]
        acc = zero(x)
        for b in blobs
            acc += blob_velocity(b, x)
        end
        U[i] = acc
    end
    return U
end

"""Remove the particles further than `delete_beyond` spans downstream of the trailing edge"""
function delete_far!(sim::FreeWakeSim)
    d = normalize(sim.V∞)
    limit = maximum(dot(x, d) for x in sim.edge) + sim.opts.delete_beyond * sim.te.span
    keep_particles!(sim.state, [dot(p.x, d) <= limit for p in sim.state.particles])
    return sim
end

"""
    run!(sim, nsteps; progress=nothing, series=nothing)

Take `nsteps` steps (`step!`), calling `progress(sim, record)` after each if given. With a `FreeWakeSeries`, every
`series.every` steps and the last are written as a ParaView frame: the wake as the step's solve saw it (written before
the step moves it) with the body's pressure and the trailing edge panel from that solve, at the step's time; the .pvd
collections are written at the end. Returns the history.
"""
function run!(sim::FreeWakeSim, nsteps::Integer; progress=nothing, series=nothing)
    for k in 1:nsteps
        save = !isnothing(series) && (sim.n % series.every == 0 || k == nsteps)
        stem = save ? joinpath(series.dir, "frame_" * lpad(sim.n, 5, '0')) : ""
        wake_files = save ? write_free_wake_vtk(sim, stem) : String[]
        rec = step!(sim)
        save && add_frame!(series, rec.t, [write_body_vtk(sim, stem); wake_files])
        isnothing(progress) || progress(sim, rec)
    end
    isnothing(series) || write_collections(series)
    return sim.history
end
