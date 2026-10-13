"""
The one analytic kernel of the method (doc §5), the velocity of a unit-strength source triangle,
w(x) = (1/4π)∫_T (x - y)/|x - y|³ dA: sources induce σw and constant vortex sheets γ × w. And the straight vortex
filament of linearly varying strength, for the edges where the doublet layer ends (the wake's tips), and the
regularized vortex particle of the free wake (`Blob`, Rosenhead–Moore).
"""
module Kernels

using LinearAlgebra
using StaticArrays

export SourceTriangle, source_velocity, sheet_velocity, source_velocity_exact, source_velocity_far, segment_velocity, FAR_DEFAULT
export Blob, blob_velocity, segment_blobs

"""Distance from the centroid, in triangle diameters, beyond which `source_velocity` uses the multipole expansion"""
const FAR_DEFAULT = 10.0

include("source_triangle.jl")
include("segment.jl")
include("blob.jl")

end # module Kernels
