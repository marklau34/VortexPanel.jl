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
include(joinpath(@__DIR__, "helpers.jl"))

@testset "VortexPanel" begin
    @testset "Geometry units" include("geometry/test_units.jl")
    @testset "Geometry closed bodies" include("geometry/test_bodies.jl")
    @testset "Geometry reference wing" include("geometry/test_wing.jl")
    @testset "Geometry pathologies" include("geometry/test_pathologies.jl")
    @testset "Mesh abuse" include("meshtools/test_abuse.jl")
    @testset "DOF map" include("dofmap/test_dofmap.jl")
end
