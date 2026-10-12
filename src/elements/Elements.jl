"""
Elements (doc §6): the element interface (basis, sheet velocity, quadrature tiers) and the flat triangle with nodal
hat functions; an isoparametric bilinear quad can be added behind the same interface.
"""
module Elements

using LinearAlgebra
using StaticArrays
using ..Geometry: SurfaceGeometry
using ..Kernels
using ..Kernels: SourceTriangle, source_velocity, sheet_velocity

export AbstractElement, FlatTriangle, build_elements, basis, point, sheet_strength, gradient, sheet_velocity, distance_to
export BaryRule, TierRules, triangle_rule, edge_graded_rule, vertex_graded_rule, self_rule, rotate_rule
export gauss_legendre, graded_gauss, graded_gauss_both

include("rules.jl")
include("flat_triangle.jl")

end # module Elements
