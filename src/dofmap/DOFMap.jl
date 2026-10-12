"""
DOF map and constraint matrix C (doc §4): nodal doublet values on the body ("slots", with trailing edge nodes
duplicated into upper and lower copies) in terms of the reduced unknowns, μ = C q, with T-junction hanging nodes,
stitched gaps and zipped openings as constraint rows, and the nodal Kutta condition as the rows K of the wake's
trailing edge nodes, μ_w = K q. The single place where topology decisions about the unknowns live.
"""
module DOFMap

using LinearAlgebra
using SparseArrays
using StaticArrays
using Printf
using WriteVTK
import MeshLibrary
using ..Geometry
using ..Geometry: SurfaceGeometry, TESide, NoTE, TEUpper, TELower, TECap, he_tri, tri_normal

export DOFMapOptions, DOFMapping, WakeTE, build_dofmap, nwake_te
export slot_values, kutta_values, restrict, wake_rows, reduce_system, bordered_rows, reproduction_error, stitch_offset, print_dofmap

include("types.jl")
include("build.jl")
include("kutta.jl")
include("system.jl")
include("export.jl")

end # module DOFMap
