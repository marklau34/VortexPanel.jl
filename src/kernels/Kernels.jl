"""
The one analytic kernel of the method (doc §5): the velocity of a unit-strength source triangle,
w(x) = (1/4π)∫_T (x - y)/|x - y|³ dA. Sources induce σw and constant vortex sheets γ × w.
"""
module Kernels

using LinearAlgebra
using StaticArrays

export SourceTriangle, source_velocity, sheet_velocity, source_velocity_exact, source_velocity_far, FAR_DEFAULT

"""Distance from the centroid, in triangle diameters, beyond which `source_velocity` uses the multipole expansion"""
const FAR_DEFAULT = 10.0

include("source_triangle.jl")

end # module Kernels
