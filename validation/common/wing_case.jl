# Solve a wing (the reference wing of reference_wing.jl, or a mesh of it) with a rigid wake and measure it by both force
# routes of doc §9. Shared by the stage-5 validation script and the tests. Requires VortexPanel, LinearAlgebra and
# StaticArrays, and reference_wing.jl for the default reference lengths. (VortexPanel.trefftz_drag is qualified: the MESA
# baseline's loads.jl defines its own.)

"""Unit freestream at angle of attack α (degrees) in the wing's frame (x chordwise, z up), and the lift direction"""
wing_freestream(α::Real) = SVector(cosd(α), 0.0, sind(α)), SVector(-sind(α), 0.0, cosd(α))

"""
1-norm condition number estimate of the bordered matrix `M` (`bordered_system`, `nq` unknowns) with its rows scaled by
the diagonal (Jacobi), the border rows by the mean diagonal entry
"""
function jacobi_condition(M::AbstractMatrix, nq::Integer)
    d = abs.(diag(M))
    d[nq+1:end] .= sum(d[1:nq]) / nq
    S = Diagonal(1 ./ d) * M
    return 1 / LinearAlgebra.LAPACK.gecon!('1', lu(S).factors, opnorm(S, 1))
end

"""
    solve_wing(mesh; α=5, geometry, dofmap, wake, assembly, span, chord, jacobi=false)

Solve the flow about a wing mesh with a rigid wake at angle of attack `α` (degrees; the wing's frame, freestream of
unit speed) and measure it: pressure force and moments (the gradient route's Cp, constant on each triangle), the
Kutta–Joukowski force and the Trefftz induced drag from the wake's μ_w, as coefficients on the reference area
span·chord (moments about the leading edge's midpoint, `Cm_LE`, and about the quarter chord, `Cm_c4`, on the chord),
with the spanwise loading, the tip loading (μ_w at the first trailing edge node over its largest value) and the
solver's diagnostics; with `jacobi=true` also the condition number with the rows scaled by the diagonal (`cond_jacobi`,
one more factorization; NaN otherwise).
"""
function solve_wing(mesh; α::Real=5.0, geometry=GeometryOptions(), dofmap=DOFMapOptions(), wake=WakeOptions(),
    assembly=AssemblyOptions(), span::Real=WING_SPAN, chord::Real=WING_CHORD, jacobi::Bool=false)
    V, eL = wing_freestream(α)
    g = build_geometry(mesh; opts=geometry)
    dm = build_dofmap(g; opts=dofmap)
    els = build_elements(g)
    w = build_wake(g, dm, V; opts=wake)
    σ = source_strengths(els, V)
    A, b, stats = assemble(g, dm, els, σ, V; opts=assembly, wake=w)
    W = wake_rows(dm, w.te)
    tsolve = @elapsed sol = solve_galerkin(dm, A, b; W=W)
    μ_w = dm.K * sol.q
    cp = pressure_coefficient(gradient_velocity(els, dm, sol.μ, V), V)
    S = span * chord
    Fp, Mp = pressure_force(els, cp)
    Fkj = kutta_joukowski_force(dm, μ_w, V)
    M_c4 = Mp - cross(SVector(0.25chord, 0.0, 0.0), Fp)
    c = consistency(dm, A, b; W=W)
    cond_jacobi = jacobi ? jacobi_condition(first(bordered_system(dm, A, b; W=W)), dm.nq) : NaN
    return (g=g, dm=dm, els=els, wake=w, sol=sol, μ_w=μ_w, cp=cp, stats=stats, V=V,
        load=spanwise_loading(els, cp, dm, μ_w, V), CL=dot(Fp, eL) / S, CL_KJ=dot(Fkj, eL) / S, CDi_p=dot(Fp, V) / S,
        CDi_T=VortexPanel.trefftz_drag(dm, μ_w, V) / S, CY=Fp[2] / S, Cm_LE=Mp[2] / (S * chord), Cm_c4=M_c4[2] / (S * chord),
        tip=μ_w[first(dm.wake_te.lines[1])] / maximum(abs, μ_w), λ=sol.λ[1], cond=sol.cond, cond_jacobi=cond_jacobi,
        columns=c.columns, nq=dm.nq, ntri=length(els), t_assembly=stats.time, t_solve=tsolve)
end

"""
    chordwise_cp(r, level, η; nspan, nchord, span=WING_SPAN)

Chordwise Cp of a `solve_wing` solution on the structured strip nearest the station η (0 root, 1 tip) of a PAGE wing
mesh, per parent quad (the area-weighted mean of its two triangles): rows [x/c Cp surface] for the upper (+1) and the
lower (-1) surface, sorted by x. Needs `structured_strip` and `structured_chord` (validation/common/loads.jl).
"""
function chordwise_cp(r, level::Integer, η::Real; nspan::Integer=reference_wing_size(level)[1], nchord::Integer=reference_wing_size(level)[2],
    span::Real=WING_SPAN)
    ys = [span * (-0.5 + (1 - cos(π * j / nspan)) / 2) for j in 0:nspan] # CosineFull stations
    mids = (ys[1:end-1] + ys[2:end]) / 2
    j = argmin(abs.(mids .- η * span / 2))
    rows = Dict{Tuple{Int,Int},Vector{Float64}}() # (chord index, surface) => [area, Cp area, x area]
    for (t, el) in enumerate(r.els)
        r.g.tri_kind[t] == VortexPanel.Geometry.BodyTri || continue
        e = r.g.tri_parent[t]
        structured_strip(e, nspan, nchord) == j || continue
        i, side = structured_chord(e, nspan, nchord)
        acc = get!(rows, (i, side), zeros(3))
        acc .+= (el.area, el.area * r.cp[t], el.area * el.centroid[1])
    end
    data = [(k[1], k[2], v[3] / v[1] / WING_CHORD, v[2] / v[1]) for (k, v) in rows]
    up = sort([(x, cp) for (i, s, x, cp) in data if s == 1])
    low = sort([(x, cp) for (i, s, x, cp) in data if s == -1])
    return reduce(vcat, [[x cp 1] for (x, cp) in up]), reduce(vcat, [[x cp -1] for (x, cp) in low])
end

"""
    abused_wings(level; seed=1)

Abused level `level` reference wings (`split_te=false` form): the MESA smoke set (validation/morino_baseline/
smoke_abuse.jl), then all of jitter, slivers, gaps and T-junctions at once. Needs VortexPanel.MeshTools and StableRNGs.
"""
function abused_wings(level::Integer; seed::Integer=1)
    nspan, nchord = reference_wing_size(level)
    base = reference_wing(level; split_te=false)
    project = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS)
    rng() = StableRNG(seed)
    skew(pattern, σ; kw...) = first(skew_wing(base, pattern, σ; nspan=nspan, nchord=nchord, rng=rng(), kw...))
    return [
        "jitter 0.2" => first(jitter_nodes(base, 0.2; rng=rng(), project=project)),
        "slivers (3 nodes, 0.02)" => first(insert_slivers(base, 3, 0.02; rng=rng(), project=project)),
        "sliver stations (2, 0.02)" => first(sliver_stations(base, 2, 0.02; nspan=nspan, nchord=nchord, rng=rng())),
        "gaps (3, 0.01)" => first(open_gaps(base, 3, 0.01; rng=rng(), project=project)),
        "T-junctions (4)" => first(insert_tjunctions(base, 4; rng=rng())),
        "lean opposite 0.3" => skew(:lean_opposite, 0.3),
        "jitter independent 0.3" => skew(:jitter_independent, 0.3),
        "TE misaligned 0.3, shared" => skew(:te_misaligned, 0.3; shared_te=true),
        "TE misaligned 0.3, non-matching" => skew(:te_misaligned, 0.3; shared_te=false),
        "lean opposite 0.3, TE rows" => skew(:lean_opposite, 0.3; te_rows=true),
        "all (jitter, slivers, gaps, T-junctions)" =>
            first(abuse(base, (jitter=0.2, slivers=(3, 0.02), gaps=(3, 0.01), tjunctions=4, project=project); seed=seed)),
    ]
end
