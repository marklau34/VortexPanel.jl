"""
Rigid wake (doc §7): fixed streamwise strips of linear-μ triangles shed from the wake's trailing edge nodes, every node
of a strip carrying the doublet strength of its trailing edge node (the Kutta rows of the DOF map), and the tip
filaments along the strips' outer sides.
"""
module Wake

using LinearAlgebra
using StaticArrays
using WriteVTK
import MeshLibrary
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping, nwake_te
using ..Elements: FlatTriangle

export WakeOptions, RigidWake, build_wake, te_panel, te_panel_size, wake_elements, wake_values, nwake

include("rigid.jl")
include("export.jl")

end # module Wake
