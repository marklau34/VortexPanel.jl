# Helpers to run MESA.jl (Morino: constant source + doublet, Dirichlet) with a rigid wake on the reference wing and
# extract the quantities the Galerkin solver is compared against. MESA's numerics are not modified: this only calls
# MESA's own functions, in the order of src/solver/solver_rigidwake.jl, to read the system matrix before it is factored.
#
# Requires MESA, MeshLibrary, PAGE, LinearAlgebra, StaticArrays, Printf, DelimitedFiles and the common/ helpers.

struct PitchFrame <: AbstractSimFrames end

const ALPHA_DEG = 5.0 # run_wing_surfpan.jl: RotY(-5°), nose up
const VINF = 30.0
const RHO = 1.225
const TE_LEN = 30 * 0.005 * 0.3 # 0.045 m, as in run_wing_surfpan.jl

const FREESTREAM = ConstantFreestream(SA[VINF, 0.0, 0.0], 0.0, RHO, 1.81e-5, 343.0)
const WING_REF = WingReference(span=WING_SPAN, chord=WING_CHORD, rho=RHO, vinf=VINF)

"""
    setup_mesa(mesh; wake_len, te_dir=:freestream, vtk_dir=nothing)

MESA solver system for a copy of the body frame `mesh` pitched to `ALPHA_DEG`, in free air, on the CPU, with a rigid
wake. The trailing edge panels (length `TE_LEN`) leave along the freestream (`te_dir=:freestream`) or along the
trailing edge bisector (`:bisector`, MESA's default); the rigid wake then runs `wake_len` along the freestream.
Returns the initialized solver system, the component (whose mesh MESA rotates into the global frame) and the loads
monitor.
"""
function setup_mesa(mesh::Mesh; wake_len::Real, te_dir::Symbol=:freestream, vtk_dir=nothing)
    refframes = (MESAFrame(FixedFrame(SA[0.0, 0.0, 0.0], RotY(deg2rad(-ALPHA_DEG), PitchFrame, MESAGlobalFrame))),)
    te_vec = te_dir == :freestream ? SA[1.0, 0.0, 0.0] : te_dir == :bisector ? SA[NaN, NaN, NaN] : error("Unknown te_dir $te_dir")
    comp = SurfacePanelComponent(deepcopy(mesh), PitchFrame; name=:wing, te_len=TE_LEN, te_dir=te_vec)
    loadmon = LoadMon(name=:wing_loads, comp=:wing, frame=MESAGlobalFrame)
    params = MESASettings(ige=false, use_gpu=false, calc_vertsol=true, no_warn=true)
    vis = isnothing(vtk_dir) ? nothing :
        Visualizer(vtks=(VisVTK(:fullmesh, vtk_dir), VisVTK(:fullwake, vtk_dir), VisVTK(:te, vtk_dir)))
    prob = MESAProblem(params, RigidWake(wake_len=wake_len), (comp,), refframes, FREESTREAM; loads_monitor=(loadmon,), visualizer=vis)
    return init(prob), comp, loadmon
end

"""
    rigidwake_system!(sys)

The complete rigid wake system matrix A and right hand side of MESA's solve, built with MESA's own functions in the
order of `solve!(::SolverSys{<:RigidWakeCache})` (src/solver/solver_rigidwake.jl), copied before MESA factors A in
place. Returns (A, RHS, assembly time).
"""
function rigidwake_system!(sys)
    (; simsettings, pansys, tesys, wakesys, aerovars) = sys
    t = @elapsed begin
        MESA.compute_influence!(pansys, sys.freestream, 0, 0.0)
        MESA.reset_linear_solver_AB!(aerovars)
        MESA.compute_influence!(aerovars, pansys, simsettings.reg_params; ige=simsettings.ige)
        MESA.compute_influence!(aerovars, pansys, tesys, simsettings.reg_params; ige=simsettings.ige)
        MESA.update_RHS!(aerovars, pansys)
        MESA.compute_influence!(aerovars, pansys, tesys, wakesys, sys.reg_params; ige=simsettings.ige)
    end
    return copy(aerovars.A), copy(aerovars.RHS), t
end

"""1-norm condition number estimate (LAPACK gecon) of A, raw and with the rows scaled by the diagonal (Jacobi)"""
function condition_numbers(A::AbstractMatrix)
    cond1(M) = 1 / LinearAlgebra.LAPACK.gecon!('1', lu(M).factors, opnorm(M, 1))
    return (raw=cond1(A), jacobi=cond1(Diagonal(1 ./ diag(A)) * A))
end

"""
    mesa_case(mesh; wake_len, te_dir, nspan, nchord, vtk_dir=nothing)

Solve one case and extract every quantity of the baseline. `nspan` and `nchord` give the structured layout of the
PAGE wing (elements appended after it, such as tip caps, belong to no strip).
"""
function mesa_case(mesh::Mesh; wake_len::Real, te_dir::Symbol=:freestream, nspan::Integer, nchord::Integer, vtk_dir=nothing)
    ref = WING_REF
    sys, comp, loadmon = setup_mesa(mesh; wake_len=wake_len, te_dir=te_dir, vtk_dir=vtk_dir)

    # system matrix, condition number, and a check that it is MESA's system
    A, b, t_asm = rigidwake_system!(sys)
    t_lu = @elapsed μ_check = lu(A) \ b
    conds = condition_numbers(A)
    t_solve = @elapsed solve!(sys)
    μ = sys.aerovars.gamSol
    check = maximum(abs, μ_check - μ) / maximum(abs, μ)

    # pitch from the body frame mesh to the mesh MESA rotated into the global frame
    R = fit_rotation(collect.(mesh.nodes), collect.(comp.mesh.nodes))
    r_c4 = R * SA[0.25 * ref.chord, 0.0, 0.0]
    coeffs = force_coefficients(loadmon.force[end], loadmon.moment[end], ref; r_c4=r_c4)

    # shed circulation per trailing edge panel: μ_w = μ_up - μ_low, Γ = -μ_w (MESA's convention, see circulation_lift)
    (; pansys, tesys) = sys
    te_bound = get_boundary(mesh, "trailing_edge_1")
    nte = tesys.nte
    @assert nte == length(te_bound)
    μ_w = tesys.panels.gam[1:nte]
    te_y = [extrema(p -> p[2], get_coords(mesh, get_element_facet(mesh, f))) for f in te_bound.facets]
    order = sortperm(te_y, by=y -> y[1] + y[2])
    y_edges = [te_y[order[1]][1]; [te_y[k][2] for k in order]]
    Γ = -μ_w[order]
    CL_KJ = kj_lift(Γ, diff(y_edges), ref) / (ref.qinf * ref.area)
    CDi_T = trefftz_drag(y_edges, Γ, ref) / (ref.qinf * ref.area)

    # spanwise loading by structured strip (body frame y extents of each strip)
    strip = [structured_strip(e, nspan, nchord) for e in 1:nelements(mesh)]
    strip_y = [extrema(p -> p[2], reduce(vcat, [collect(get_coords(mesh, mesh.elements[e])) for e in findall(==(j), strip)])) for j in 1:nspan]
    dy = [y[2] - y[1] for y in strip_y]
    η = [(y[1] + y[2]) / ref.span for y in strip_y]
    Fz = [sum(e -> pansys.pan_force[e][3], findall(==(j), strip)) for j in 1:nspan]
    Fx = [sum(e -> pansys.pan_force[e][1], findall(==(j), strip)) for j in 1:nspan]
    te_strip = [strip[tesys.te_pan_idx[1, k]] for k in 1:nte]
    μ_w_strip = zeros(nspan)
    μ_w_strip[te_strip] .= μ_w
    spanwise = hcat(η, (η .* ref.span ./ 2), dy, Fz ./ (ref.qinf * ref.chord .* dy), Fx ./ (ref.qinf * ref.chord .* dy), μ_w_strip, -μ_w_strip)

    # chordwise pressure coefficient per strip, upper and lower surface (body frame x/c of the panel centroids)
    cp = pansys.pan_pres ./ ref.qinf
    xc = [sum(p -> p[1], get_coords(mesh, mesh.elements[e])) / (length(get_nodes(mesh.elements[e])) * ref.chord) for e in 1:2nspan*nchord]
    cp_up = [cp[(j - 1) * nchord + i] for j in 1:nspan, i in 1:nchord]
    cp_low = [cp[nspan * nchord + (j - 1) * nchord + i] for j in 1:nspan, i in 1:nchord]
    xc_up = [xc[(j - 1) * nchord + i] for j in 1:nspan, i in 1:nchord]

    timings = (assembly=t_asm, lu=t_lu, solve=t_solve)
    return (; coeffs..., CL_KJ, CDi_T, alpha=rad2deg(pitch_angle(R)), cond_raw=conds.raw, cond_jacobi=conds.jacobi,
        check, timings, npan=pansys.npan, spanwise, μ_w=μ_w[order], y_edges, η, cp_up, cp_low, xc_up)
end
