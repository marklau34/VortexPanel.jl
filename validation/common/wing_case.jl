# Solve a wing (the reference wing of reference_wing.jl, or a mesh of it) with a rigid wake and measure it by both force
# routes of doc §9. Shared by the stage-5 validation script and the tests. Requires VortexPanel, LinearAlgebra and
# StaticArrays, and reference_wing.jl for the default reference lengths. (VortexPanel.trefftz_drag is qualified: the MESA
# baseline's loads.jl defines its own.)

"""Unit freestream at angle of attack α (degrees) in the wing's frame (x chordwise, z up), and the lift direction"""
wing_freestream(α::Real) = SVector(cosd(α), 0.0, sind(α)), SVector(-sind(α), 0.0, cosd(α))

"""
    solve_wing(mesh; α=5, geometry, dofmap, wake, assembly, span, chord)

Solve the flow about a wing mesh with a rigid wake at angle of attack `α` (degrees; the wing's frame, freestream of
unit speed) and measure it: pressure force and moments (the gradient route's Cp, constant on each triangle), the
Kutta–Joukowski force and the Trefftz induced drag from the wake's μ_w, as coefficients on the reference area
span·chord (moments about the leading edge's midpoint, `Cm_LE`, and about the quarter chord, `Cm_c4`, on the chord),
with the spanwise loading, the tip loading (μ_w at the first trailing edge node over its largest value) and the
solver's diagnostics.
"""
function solve_wing(mesh; α::Real=5.0, geometry=GeometryOptions(), dofmap=DOFMapOptions(), wake=WakeOptions(),
    assembly=AssemblyOptions(), span::Real=WING_SPAN, chord::Real=WING_CHORD)
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
    return (g=g, dm=dm, els=els, wake=w, sol=sol, μ_w=μ_w, cp=cp, stats=stats, V=V,
        load=spanwise_loading(els, cp, dm, μ_w, V), CL=dot(Fp, eL) / S, CL_KJ=dot(Fkj, eL) / S, CDi_p=dot(Fp, V) / S,
        CDi_T=VortexPanel.trefftz_drag(dm, μ_w, V) / S, CY=Fp[2] / S, Cm_LE=Mp[2] / (S * chord), Cm_c4=M_c4[2] / (S * chord),
        tip=μ_w[first(dm.wake_te.lines[1])] / maximum(abs, μ_w), λ=sol.λ[1], cond=sol.cond, columns=c.columns,
        nq=dm.nq, ntri=length(els), t_assembly=stats.time, t_solve=tsolve)
end
