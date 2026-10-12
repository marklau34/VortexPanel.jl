"""
Angles (degrees), aspect ratio R/(2r) (1 for an equilateral triangle) and equiangle skew
max((θmax - 60)/120, (60 - θmin)/60) of a triangle, and its area
"""
function triangle_quality(a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}) where {T}
    la, lb, lc = norm(c - b), norm(a - c), norm(b - a) # side opposite each vertex
    area = norm(cross(b - a, c - a)) / 2
    angle(u, v) = rad2deg(acos(clamp(dot(u, v) / (norm(u) * norm(v)), -one(T), one(T))))
    θ = (angle(b - a, c - a), angle(a - b, c - b), angle(a - c, b - c))
    θmin, θmax = minimum(θ), maximum(θ)
    s = (la + lb + lc) / 2
    aspect = area > 0 ? la * lb * lc * s / (8 * area^2) : T(Inf)
    skew = max((θmax - 60) / 120, (60 - θmin) / 60)
    return θmin, θmax, aspect, skew, area
end

"""Per triangle quality and the slivers (smallest angle or aspect ratio beyond the thresholds)"""
function quality_arrays(x::AbstractVector{SVector{3,T}}, tris, opts::GeometryOptions) where {T}
    n = length(tris)
    qmin, qmax, asp, skew, area = zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n)
    for (k, t) in enumerate(tris)
        qmin[k], qmax[k], asp[k], skew[k], area[k] = triangle_quality(x[t[1]], x[t[2]], x[t[3]])
    end
    sliver = (qmin .< T(opts.sliver_angle_deg)) .| (asp .> T(opts.sliver_aspect))
    return qmin, qmax, asp, skew, area, BitVector(sliver)
end

const REPORT_COUNTS = (:nodes, :triangles, :cap_triangles, :edges, :bodies, :closed_bodies, :input_nodes, :input_elements,
    :welded_nodes, :dropped_elements, :collapsed_elements, :split_shortest, :split_fold, :split_tie, :split_fixed,
    :flipped_triangles, :bowtie_vertices, :tjunctions, :gaps, :kutta_pairs, :openings_capped, :openings_zipped,
    :openings_open, :kutta_edges, :kutta_chains, :te_nodes, :kutta_tag_mismatch, :slivers)

"""
    print_report(io=stdout, g)

Mesh quality and repair report of a `SurfaceGeometry`.
"""
function print_report(io::IO, g::SurfaceGeometry)
    c, s = g.quality.counts, g.quality.stats
    println(io, "Geometry report")
    @printf(io, "  %-22s %d → %d (%d welded), %d elements → %d triangles (%d caps)\n", "nodes / elements",
        c[:input_nodes], c[:nodes], c[:welded_nodes], c[:input_elements], c[:triangles], c[:cap_triangles])
    @printf(io, "  %-22s %d bodies (%d closed), Euler characteristic %s\n", "topology", c[:bodies], c[:closed_bodies],
        join([string(b.euler) for b in g.bodies], ", "))
    @printf(io, "  %-22s shortest %d, fold %d, tie %d, fixed %d\n", "quad split", c[:split_shortest], c[:split_fold], c[:split_tie], c[:split_fixed])
    @printf(io, "  %-22s dropped %d, collapsed %d, flipped %d, bow-ties %d, max weld displacement %.3g\n", "repairs",
        c[:dropped_elements], c[:collapsed_elements], c[:flipped_triangles], c[:bowtie_vertices], s[:max_weld])
    @printf(io, "  %-22s T-junctions %d, gaps %d (max width %.3g), non-matching trailing edges %d\n", "pathologies",
        c[:tjunctions], c[:gaps], s[:max_gap], c[:kutta_pairs])
    @printf(io, "  %-22s capped %d, zipped %d, open %d, max non-planarity %.3g\n", "openings",
        c[:openings_capped], c[:openings_zipped], c[:openings_open], s[:max_cap_planarity])
    @printf(io, "  %-22s %d edges in %d chains, %d nodes, %d tag mismatches\n", "trailing edge",
        c[:kutta_edges], c[:kutta_chains], c[:te_nodes], c[:kutta_tag_mismatch])
    @printf(io, "  %-22s min angle %.2f°, max angle %.2f°, max aspect %.3g, max skew %.3f, %d slivers\n", "quality",
        s[:min_angle], s[:max_angle], s[:max_aspect], s[:max_skew], c[:slivers])
    @printf(io, "  %-22s %s\n", "volume", join([@sprintf("%.6g", b.volume) for b in g.bodies], ", "))
    return nothing
end

print_report(g::SurfaceGeometry) = print_report(stdout, g)

"""Report counts and magnitudes as a flat NamedTuple, for tables"""
function report_row(g::SurfaceGeometry)
    c, s = g.quality.counts, g.quality.stats
    counts = NamedTuple{REPORT_COUNTS}(Tuple(c[k] for k in REPORT_COUNTS))
    stats = (min_angle=s[:min_angle], max_angle=s[:max_angle], max_aspect=s[:max_aspect], max_skew=s[:max_skew],
        max_weld=s[:max_weld], max_gap=s[:max_gap], max_cap_planarity=s[:max_cap_planarity],
        volume=sum(b.volume for b in g.bodies; init=zero(s[:min_angle])))
    return merge(counts, stats)
end
