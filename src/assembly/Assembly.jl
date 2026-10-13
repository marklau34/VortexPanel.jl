"""
Galerkin assembly and solution (doc §6, §8): local 3×3 blocks of every test and source triangle pair with quadrature
tiers chosen by their distance (graded rules for touching pairs), the right hand side from the freestream (exact) and
the prescribed sources, the test-weight switch, the external-velocity hook, the bordered dense solve (factored once for
many right hand sides), and the right hand side of explicit vortex elements (a free wake and its particles).
"""
module Assembly

using LinearAlgebra
using StaticArrays
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping, reduce_system, bordered_rows
using ..Kernels: source_velocity, segment_velocity, FAR_DEFAULT, Blob, blob_velocity
using ..Elements
using ..Wake: RigidWake, nwake

export AssemblyOptions, assemble, source_strengths, bordered_system, solve_galerkin, consistency, near_rule, contact_rule, pair_block
export factorize_galerkin, GalerkinFactorization, bordered_matrix, bordered_rhs, ExplicitSources, explicit_rhs, Segment
export MixedMesh, square_contact_rule

include("assemble.jl")
include("mixed.jl")
include("explicit.jl")
include("solve.jl")

end # module Assembly
