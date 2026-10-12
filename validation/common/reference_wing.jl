# Reference wing of MESA.jl/examples/wing/run_wing_surfpan.jl (create_simple_wing_surfpan in MESA.jl/examples/util.jl):
# rectangular NACA 0015 wing, untwisted, spanning y = -b/2 .. b/2 with the leading edge on x = 0.
# Requires PAGE, MeshLibrary and StaticArrays to be loaded.

const WING_SPAN = 2.708707
const WING_CHORD = 0.69776
const WING_THICKNESS = 0.15

"""
    reference_wing_size(level)

Number of spanwise and chordwise (per surface) elements of refinement `level`: (17, 19) doubled per level. The cosine
spacings nest exactly, so every node of a level is a node of the next.
"""
reference_wing_size(level::Integer) = (17 * 2^(level - 1), 19 * 2^(level - 1))

"""
    reference_wing(level=1; split_te=true, nspan, nchord)

PAGE mesh of the reference wing with the parameters of `create_simple_wing_surfpan`. All quads, upper surface
elements first, chordwise index fastest from the leading edge, strips in +y. The tips are open. The trailing edge is
tagged with the boundary `trailing_edge` and the element sets `trailing_edge_up`/`trailing_edge_low`; with
`split_te=true` the trailing edge nodes are duplicated and the boundaries `trailing_edge_1`/`trailing_edge_2` added,
as MESA requires.
"""
function reference_wing(level::Integer=1; split_te::Bool=true, nspan::Integer=reference_wing_size(level)[1],
    nchord::Integer=reference_wing_size(level)[2])
    airfoil = AirfoilNACA4(0.0, 0.0, WING_THICKNESS, :airfoil1)
    airfoils = LinearAirfoilSections([0.0, 1.0], [airfoil, airfoil])
    rig = ParametricWingRig(posy=LinearSpanInput(WING_SPAN), chord=WING_CHORD, refloc=0.0, translation=SA[0.0, -WING_SPAN/2, 0.0])
    wing = ParametricWing(rig, airfoils)
    mesh_params = LiftSurfMeshParams(span_space=CosineFull(), chord_space=CosineHalf(option=:start), nspan=nspan, nchord=nchord)
    te_params = TrailingEdgeParams(te_name="trailing_edge", te_angle=30)
    return PAGE.generate_mesh(wing, mesh_params, trailing_edge=te_params, split_te=split_te)
end

"""Half thickness of the NACA 00xx section (sharp trailing edge, as PAGE) at x/c = ξ"""
function naca00_half_thickness(ξ::Real; t::Real=WING_THICKNESS)
    ξ = clamp(ξ, zero(ξ), one(ξ))
    return 5t * (0.2969sqrt(ξ) - 0.1260ξ - 0.3516ξ^2 + 0.2843ξ^3 - 0.1036ξ^4)
end

"""Volume of the exact reference wing: section area ∫2y_t dξ c² times the span"""
function reference_wing_volume(; t::Real=WING_THICKNESS)
    area = 10t * (0.2969 * 2/3 - 0.1260/2 - 0.3516/3 + 0.2843/4 - 0.1036/5)
    return area * WING_CHORD^2 * WING_SPAN
end
