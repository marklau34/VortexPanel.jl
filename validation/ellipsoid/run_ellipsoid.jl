# Validation stage 4: ellipsoids in uniform flow (doc §11), for a non-uniform analytic Cp. The exact solution (Lamb §114,
# validation/common/exact_flows.jl) keeps μ linear on the surface, μ = Σ αᵢ/(2 - αᵢ) V∞ᵢ xᵢ (αᵢ the depolarization
# factors), while the surface velocity is the tangential part of the constant vector k, kᵢ = 2V∞ᵢ/(2 - αᵢ), so Cp
# follows the curvature.
#
#   julia -t auto --project=validation validation/ellipsoid/run_ellipsoid.jl [levels=1,2,3,4]
#
# Bodies, at the generic V∞ of the sphere study: a triaxial ellipsoid 2 : 1 : 0.5; a slender prolate spheroid 6 : 1 : 1,
# whose scaled meshes carry panels of aspect ratio up to 6 with corners near 150°; thin oblate spheroids 1 : 1 : 0.2
# and 1 : 1 : 0.1, with opposite faces close together and a tight rim (radius of curvature 0.04 and 0.01).
# 1. The exact solution: the factors, and the surface formulas against the exterior potential.
# 2. Convergence on the triaxial ellipsoid with scaled icospheres (clean and jittered), cube-spheres and
#    latitude-longitude spheres (pole fans of thin triangles), and on the spheroids with icospheres.
# 3. Abused meshes of the triaxial ellipsoid (doc §11): jitter, slivers, gaps, T-junctions and all at once, moved
#    nodes put back on the surface along the ray from the centre.
# 4. Sensitivity to the quadrature tiers on the triaxial, the slender and the thinner oblate body.
# 5. Cp along the cut y = 0 of the triaxial ellipsoid against the exact one (CSV for plots), and VTK of the finest
#    meshes.
# Writes validation/output/ellipsoid/.

include(joinpath(@__DIR__, "..", "common", "closed_body.jl"))

const OUTDIR = joinpath(@__DIR__, "..", "output", "ellipsoid")
const VINF = normalize(SVector(1.0, 0.3, -0.2))
const TRIAXIAL = EllipsoidFlow((2.0, 1.0, 0.5), VINF)
const PROLATE = EllipsoidFlow((3.0, 0.5, 0.5), VINF)
const OBLATE = EllipsoidFlow((1.0, 1.0, 0.2), VINF)
const THIN = EllipsoidFlow((1.0, 1.0, 0.1), VINF)
const BODIES = [("triaxial 2:1:0.5", TRIAXIAL), ("prolate 6:1:1", PROLATE), ("oblate 1:1:0.2", OBLATE), ("oblate 1:1:0.1", THIN)]

function parse_args(args)
    levels = [1, 2, 3, 4]
    for a in args
        key, val = split(a, "=")
        key == "levels" && (levels = parse.(Int, split(val, ",")))
    end
    return levels
end

"""The exact solution of every body, and its surface formulas checked against the exterior potential"""
function exact_checks()
    println("\nExact solutions: depolarization factors α, the vectors k_μ (μ = k_μ·x) and k (velocity = tangential part of k),",
        " and the surface formulas against the exterior potential φ at four surface points (finite differences)")
    rows = Any[]
    e(i) = SVector(ntuple(j -> Float64(j == i), 3))
    for (name, flow) in BODIES
        φ(y) = perturbation_potential(flow, y)
        worst = zeros(3)
        for p in (SVector(0.3, 0.5, 0.7), SVector(-0.8, 0.1, 0.2), SVector(0.0, -1.0, 0.1), SVector(0.6, -0.4, -0.5))
            x = surface_point(flow, p)
            n = surface_normal(flow, x)
            h = 1e-4 * minimum(flow.axes)
            neumann = (-3φ(x) + 4φ(x + h * n) - φ(x + 2h * n)) / (2h) + dot(flow.V∞, n)
            y, hh = x + 0.3 * minimum(flow.axes) * n, 1e-3 * minimum(flow.axes)
            laplacian = sum((φ(y + hh * e(i)) - 2φ(y) + φ(y - hh * e(i))) / hh^2 for i in 1:3)
            worst .= max.(worst, abs.((φ(x) - exact_doublet(flow, x), neumann, laplacian)))
        end
        @printf("%-17s α = (%.6f, %.6f, %.6f), Σα - 2 = %8.1e, k_μ = %s, k = %s;  |φ - μ| %.0e, |∂φ/∂n + V∞·n| %.0e, |Δφ| %.0e\n",
            name, flow.factors..., sum(flow.factors) - 2, round.(flow.kμ, digits=4), round.(flow.k, digits=4), worst...)
        push!(rows, [name flow.axes... flow.factors... flow.kμ... flow.k... worst...])
    end
    return rows
end

"""
Cp along the cut y = 0: for every triangle crossing it, at the middle of the crossing (the angle about the y axis, x
and z), the exact value there (the smooth curve) and at the centroid, where the computed values (both routes) apply
"""
function cut_y0(flow::EllipsoidFlow, r)
    rows = Any[]
    for (t, tri) in enumerate(r.g.tris)
        p = [r.g.x[v] for v in tri]
        crossings = SVector{3,Float64}[]
        for k in 1:3
            a, b = p[k], p[mod1(k + 1, 3)]
            (a[2] < 0) != (b[2] < 0) && push!(crossings, a + (b - a) * (a[2] / (a[2] - b[2])))
        end
        length(crossings) == 2 || continue
        m = sum(crossings) / 2
        s = surface_point(flow, m)
        push!(rows, [atand(s[3], s[1]) s[1] s[3] exact_cp(flow, m) r.cpex[t] r.cpg[t] r.cp[t]])
    end
    return sortslices(reduce(vcat, rows), dims=1)
end

function main(args=ARGS)
    levels = parse_args(args)
    mkpath(OUTDIR)
    println("Ellipsoids in uniform flow, V∞ = $(round.(VINF, digits=4)), threads $(Threads.nthreads())")
    writedlm(joinpath(OUTDIR, "exact.csv"), vcat(["body" "a" "b" "c" "alpha" "beta" "gamma" "kmu_x" "kmu_y" "kmu_z" "k_x" "k_y" "k_z" "phi_minus_mu" "neumann" "laplacian"],
        reduce(vcat, exact_checks())), ',')

    project(flow) = p -> surface_point(flow, p)
    nθ = [6 * 2^(l - 1) for l in levels] # latitude-longitude spheres nθ × 2nθ
    conv = Any[]
    rows, tri_ico = convergence("triaxial, icospheres", TRIAXIAL, [scale_nodes(icosphere(l), TRIAXIAL.axes) for l in levels])
    append!(conv, rows)
    rows, _ = convergence("triaxial, icospheres jittered 0.3", TRIAXIAL,
        [first(abuse(scale_nodes(icosphere(l), TRIAXIAL.axes), (jitter=0.3, project=project(TRIAXIAL)); seed=1)) for l in levels])
    append!(conv, rows)
    rows, _ = convergence("triaxial, cube-spheres", TRIAXIAL, [scale_nodes(cube_sphere(3 * 2^(l - 1)), TRIAXIAL.axes) for l in levels])
    append!(conv, rows)
    rows, _ = convergence("triaxial, latitude-longitude", TRIAXIAL, [scale_nodes(uv_sphere(n, 2n), TRIAXIAL.axes) for n in nθ])
    append!(conv, rows)
    finest = Dict{String,Any}()
    for (name, flow) in BODIES[2:end]
        rows, res = convergence("$name, icospheres", flow, [scale_nodes(icosphere(l), flow.axes) for l in levels])
        append!(conv, rows)
        finest[name] = res[end]
    end
    writedlm(joinpath(OUTDIR, "convergence.csv"), vcat(CONVERGENCE_HEADER, reduce(vcat, conv)), ',')

    rows_ab, ab = abused_meshes("triaxial icosphere 3", TRIAXIAL, scale_nodes(icosphere(3), TRIAXIAL.axes), project(TRIAXIAL))
    rows_ab_cube, _ = abused_meshes("triaxial cube-sphere 12", TRIAXIAL, scale_nodes(cube_sphere(12), TRIAXIAL.axes), project(TRIAXIAL))
    writedlm(joinpath(OUTDIR, "abuse.csv"), vcat(ABUSE_HEADER, reduce(vcat, [rows_ab; rows_ab_cube])), ',')

    sens = Any[]
    for (name, flow) in (BODIES[1], BODIES[2], BODIES[4])
        append!(sens, quadrature_sensitivity("$name icosphere 3", flow, scale_nodes(icosphere(3), flow.axes)))
    end
    writedlm(joinpath(OUTDIR, "quadrature.csv"), vcat(QUADRATURE_HEADER, reduce(vcat, sens)), ',')

    cut = cut_y0(TRIAXIAL, tri_ico[end])
    writedlm(joinpath(OUTDIR, "cp_cut_y0.csv"), vcat(["angle_deg" "x" "z" "cp_exact" "cp_exact_centroid" "cp_gradient" "cp_sums"], cut), ',')
    @printf("\nCp on the cut y = 0 of the finest triaxial mesh (%d crossing triangles): largest error at the centroids %.2e (gradient), %.2e (sums)\n",
        size(cut, 1), maximum(abs, cut[:, 6] - cut[:, 5]), maximum(abs, cut[:, 7] - cut[:, 5]))

    write_case(TRIAXIAL, tri_ico[end], OUTDIR, "triaxial_icosphere_$(levels[end])")
    write_case(TRIAXIAL, last(ab[end]), OUTDIR, "triaxial_icosphere_3_abused")
    write_case(PROLATE, finest["prolate 6:1:1"], OUTDIR, "prolate_icosphere_$(levels[end])")
    write_case(THIN, finest["oblate 1:1:0.1"], OUTDIR, "thin_oblate_icosphere_$(levels[end])")
    println("\nWrote ", normpath(OUTDIR))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
