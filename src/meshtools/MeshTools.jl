"""
Mesh generators and mesh abuse for testing and validation: closed bodies (sphere, cube-sphere, ellipsoid) and
seeded, reproducible pathologies (node jitter and skew, T-junctions, small gaps, slivers, non-matching trailing edges).
All operate on and return MeshLibrary meshes, so the same meshes feed this package and MESA.jl.
"""
module MeshTools

using LinearAlgebra
using StaticArrays
using Random
using StableRNGs
import MeshLibrary

export icosphere, cube_sphere, uv_sphere, scale_nodes, convert_nodes, reverse_elements, mirror_nodes
export AbuseLog, abuse, jitter_nodes, insert_tjunctions, open_gaps, insert_slivers, default_avoid, default_fixed
export WingGrid, wing_grid, skew_shifts, skew_wing, sliver_stations, naca_projector

include("bodies.jl")
include("abuse.jl")
include("wing_skew.jl")

end # module MeshTools
