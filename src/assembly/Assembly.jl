"""
Galerkin assembly and solution (doc §6, §8): local 3×3 blocks of every test and source triangle pair with quadrature
tiers chosen by their distance (graded rules for touching pairs), the right hand side from the freestream (exact) and
the prescribed sources, the test-weight switch, the external-velocity hook, and the bordered dense solve.
"""
module Assembly

using LinearAlgebra
using StaticArrays
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping, reduce_system, bordered_rows
using ..Kernels: source_velocity, FAR_DEFAULT
using ..Elements

export AssemblyOptions, assemble, source_strengths, solve_galerkin, consistency, near_rule, contact_rule, pair_block

include("assemble.jl")
include("solve.jl")

end # module Assembly
