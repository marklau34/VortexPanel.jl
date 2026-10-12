"""
Galerkin linear-vortex panel solver: nodal doublet strengths on flat triangles, constant vortex sheets γ = n×∇μ,
prescribed sources, a nodal Kutta condition and a rigid wake. See the design document in the repository root.
"""
module VortexPanel

include("geometry/Geometry.jl")
include("dofmap/DOFMap.jl")
include("kernels/Kernels.jl")
include("meshtools/MeshTools.jl")

using .Geometry
using .DOFMap
export Geometry, DOFMap, Kernels, MeshTools
export build_geometry, GeometryOptions, SurfaceGeometry, print_report, report_row, add_caps
export build_dofmap, DOFMapOptions, DOFMapping, print_dofmap

end # module VortexPanel
