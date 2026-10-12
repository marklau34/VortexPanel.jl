# Cross-check of VortexPanel's source triangle kernel against the source panel velocity of AeroSingularityElements
# (vel_srcpan, DUST-derived, used by MESA.jl): the doc asks to reuse that routine if it fits. Both should give the same
# field at regular points. Where they differ, VortexPanel's kernel is the one validated against independent quadrature
# (test/kernels): AeroSingularityElements zeroes an edge's log term within an absolute 1e-8 of the edge, has no
# one-sided limits on the panel's plane (it returns the principal value), and has no far-field switch.
#
#   julia --project=validation validation/kernels/compare_ase.jl

using MESA, LinearAlgebra, StaticArrays, Printf, StableRNGs
using VortexPanel.Kernels

const ASE = MESA.AeroSingularityElements

ase_velocity(a, b, c, x) = ASE.calc_vel_source(x, 1.0, ASE.DoubletSourcePanel([a, b, c]), ASE.RegularizationParams())

function main(; ntri=200, seed=1)
    rng = StableRNG(seed)
    categories = Dict{String,Vector{Float64}}()
    record!(k, e) = push!(get!(categories, k, Float64[]), e)
    for _ in 1:ntri
        a, b, c = (SVector{3}(randn(rng, 3)) for _ in 1:3)
        tri = SourceTriangle(a, b, c)
        h, n = tri.diam, tri.n
        u = normalize(SVector{3}(randn(rng, 3)))
        pts = [("far (5 diameters)", tri.centroid + 5h * u), ("near (0.2 diameters)", tri.centroid + 0.2h * u),
            ("just above (1e-4 diameters)", tri.centroid + 1e-4h * n),
            ("edge-adjacent (1e-4 diameters above)", (a + b) / 2 + 1e-4h * n),
            ("edge-adjacent (1e-9 diameters above)", (a + b) / 2 + 1e-9h * n),
            ("in-plane inside (principal value)", tri.centroid)]
        for (k, x) in pts
            ours = source_velocity(tri, x; far=Inf)
            theirs = ase_velocity(a, b, c, x)
            record!(k, norm(ours - theirs) / norm(ours))
        end
    end
    println("Relative difference |w_VortexPanel - w_ASE| / |w_VortexPanel| over $ntri random triangles")
    @printf("%-40s %10s %10s\n", "points", "median", "max")
    for k in ("far (5 diameters)", "near (0.2 diameters)", "just above (1e-4 diameters)", "edge-adjacent (1e-4 diameters above)",
        "edge-adjacent (1e-9 diameters above)", "in-plane inside (principal value)")
        e = sort(categories[k])
        @printf("%-40s %10.1e %10.1e\n", k, e[cld(length(e), 2)], e[end])
    end
    return categories
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
