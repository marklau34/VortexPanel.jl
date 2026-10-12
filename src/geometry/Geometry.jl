"""
Geometry and topology layer: turns an imperfect surface mesh of triangles and quads into a closed, consistently
oriented triangulation with classified edges, repairing (welding, quad splitting, T-junctions, gaps, tip caps) and
logging rather than rejecting. Its output, `SurfaceGeometry`, is the single input of the DOF map, wake and assembly.
"""
module Geometry

using LinearAlgebra
using StaticArrays
using Printf
using WriteVTK
import MeshLibrary

export build_geometry, GeometryOptions, SurfaceGeometry, NonManifoldEdgeError
export print_report, report_row, add_caps
export TriKind, BodyTri, CapTri, EdgeClass, Interior, Kutta, KuttaPair, Gap, TJLong, TJShort, Open, CapSeam
export TESide, NoTE, TEUpper, TELower, TECap
export he_index, he_tri, he_local, he_next, he_prev, he_origin, he_dest, tri_normal

include("types.jl")
include("options.jl")
include("util.jl")
include("ingest.jl")
include("weld.jl")
include("quadsplit.jl")
include("halfedge.jl")
include("pairs.jl")
include("orient.jl")
include("kutta.jl")
include("openings.jl")
include("cap.jl")
include("quality.jl")
include("pipeline.jl")
include("export.jl")

end # module Geometry
