"""
Elements (doc §6): the element interface (basis, sheet velocity, quadrature tiers), the flat triangle with nodal hat
functions, and the isoparametric bilinear quad with its numerical inner integrals and parameter-square rules.
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
export IsoQuad, shape, shape_derivatives, tangents, area_normal, unit_normal, sheet_densities, project, reach, quad_field
export SquareRule, SquareTiers, square_gauss, square_edge_rule, square_vertex_rule, square_self_rule, square_near_rule, turn

include("rules.jl")
include("flat_triangle.jl")
include("iso_quad.jl")
include("quad_field.jl")
include("quad_rules.jl")

end # module Elements
