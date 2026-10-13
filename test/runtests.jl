using Test
using LinearAlgebra
using StaticArrays
using StableRNGs
using MeshLibrary
using PAGE
using VortexPanel
using VortexPanel.Geometry
using VortexPanel.MeshTools

const G = VortexPanel.Geometry

include(joinpath(@__DIR__, "..", "validation", "common", "reference_wing.jl"))
include(joinpath(@__DIR__, "..", "validation", "common", "exact_flows.jl"))
include(joinpath(@__DIR__, "..", "validation", "common", "wing_case.jl"))
include(joinpath(@__DIR__, "..", "validation", "common", "panel2d.jl"))
include(joinpath(@__DIR__, "helpers.jl"))

@testset "VortexPanel" begin
    @testset "Geometry units" include("geometry/test_units.jl")
    @testset "Geometry closed bodies" include("geometry/test_bodies.jl")
    @testset "Geometry reference wing" include("geometry/test_wing.jl")
    @testset "Geometry pathologies" include("geometry/test_pathologies.jl")
    @testset "Mesh abuse" include("meshtools/test_abuse.jl")
    @testset "DOF map" include("dofmap/test_dofmap.jl")
    @testset "Kernels" include("kernels/test_kernels.jl")
    @testset "Segment kernel" include("kernels/test_segment.jl")
    @testset "Vortex particle kernel" include("kernels/test_blob.jl")
    @testset "Elements" include("elements/test_elements.jl")
    @testset "Isoparametric quads" include("elements/test_iso_quad.jl")
    @testset "Touching and contact quadrature" include("assembly/test_contact.jl")
    @testset "Sphere (validation stage 3)" include("assembly/test_sphere.jl")
    @testset "Ellipsoids (validation stage 4)" include("assembly/test_ellipsoid.jl")
    @testset "Wake" include("wake/test_wake.jl")
    @testset "Force routes" include("post/test_forces.jl")
    @testset "Wing (validation stage 5)" include("assembly/test_wing.jl")
    @testset "Two-dimensional section (validation stage 6)" include("assembly/test_section2d.jl")
    @testset "Element comparison (validation stage 7)" include("assembly/test_quads.jl")
    @testset "Test weights (validation stage 8)" include("assembly/test_test_weights.jl")
    @testset "Free wake with vortex particles" include("freewake/test_freewake.jl")
end
