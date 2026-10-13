"""
Free wake with vortex particles (the design doc's first deferred item), with MESA's time stepping: an implicit trailing
edge panel (`te_panel`) carries the Kutta condition and is factored into the matrix once; behind it, free panel rows
shed every time step keep the doublets they were shed with and move with the flow (forward Euler), and the oldest row
is converted into vortex particles once more than `max_rows` rows exist. The conversion is the design doc's, local to
each wake triangle (`convert_oldest!`): one particle per triangle from its own corners and doublets, none across the
interface edge. The free rows and the particles couple into the body only through their velocity, in the right hand
side (`explicit_rhs`). No vortex stretching.
"""
module FreeWake

using LinearAlgebra
using StaticArrays
using WriteVTK
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping, nwake_te, wake_rows
using ..Kernels: Blob, blob_velocity, segment_blobs
using ..Elements: FlatTriangle, sheet_strength, build_elements
using ..Wake: RigidWake, te_panel
using ..Assembly: AssemblyOptions, assemble, source_strengths, factorize_galerkin, solve_galerkin, GalerkinFactorization,
    explicit_rhs, Segment
import ..Assembly: ExplicitSources
using ..Post: body_velocity, gradient_velocity, pressure_coefficient, pressure_force, kutta_joukowski_force,
    element_values, sheet_strengths

export FreeWakeOptions, FreeWakeSim, FreeWakeState, step!, run!, vortex_elements, convert_oldest!, free_wake_blobs,
    vorticity_moments, write_free_wake_vtk

include("state.jl")
include("step.jl")
include("export.jl")

end # module FreeWake
