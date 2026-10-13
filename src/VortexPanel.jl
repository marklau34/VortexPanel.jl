"""
Galerkin linear-vortex panel solver: nodal doublet strengths on flat triangles, constant vortex sheets γ = n×∇μ,
prescribed sources, a nodal Kutta condition, and a rigid wake or a free wake shed into vortex particles. See the design
document in the repository root.
"""
module VortexPanel

include("geometry/Geometry.jl")
include("dofmap/DOFMap.jl")
include("kernels/Kernels.jl")
include("elements/Elements.jl")
include("wake/Wake.jl")
include("assembly/Assembly.jl")
include("post/Post.jl")
include("freewake/FreeWake.jl")
include("meshtools/MeshTools.jl")

using .Geometry
using .DOFMap
using .Elements
using .Wake
using .Assembly
using .Post
using .FreeWake
export Geometry, DOFMap, Kernels, Elements, Wake, Assembly, Post, FreeWake, MeshTools
export build_geometry, GeometryOptions, SurfaceGeometry, print_report, report_row, add_caps
export build_dofmap, DOFMapOptions, DOFMapping, print_dofmap, wake_rows
export build_wake, WakeOptions, te_panel
export build_elements, AssemblyOptions, assemble, source_strengths, bordered_system, solve_galerkin, consistency, MixedMesh
export factorize_galerkin, explicit_rhs, ExplicitSources
export surface_velocity, gradient_velocity, pressure_coefficient, body_velocity
export pressure_force, kutta_joukowski_force, trefftz_drag, spanwise_loading, surface_pressure
export FreeWakeOptions, FreeWakeSim, step!, run!, write_free_wake_vtk, write_body_vtk, FreeWakeSeries

end # module VortexPanel
