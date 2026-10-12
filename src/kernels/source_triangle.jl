"""
    SourceTriangle(a, b, c)

A flat triangle with vertices a, b, c (counterclockwise about its unit normal n ∝ (b-a)×(c-a)) and the data its
source kernel needs: centroid, area, diameter (longest edge), and the in-plane outward unit normal and length of
each edge (edge k runs from vertex k to vertex k+1).
"""
struct SourceTriangle{T<:Real}
    v::SVector{3,SVector{3,T}}
    n::SVector{3,T}
    centroid::SVector{3,T}
    area::T
    diam::T
    ν::SVector{3,SVector{3,T}}
    len::SVector{3,T}
end

function SourceTriangle(a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}) where {T<:Real}
    nA = cross(b - a, c - a)
    A2 = norm(nA)
    A2 > 0 || throw(ArgumentError("degenerate triangle $a, $b, $c"))
    n = nA / A2
    len = SVector(norm(b - a), norm(c - b), norm(a - c))
    ν = SVector(cross(b - a, n) / len[1], cross(c - b, n) / len[2], cross(a - c, n) / len[3])
    return SourceTriangle(SVector(a, b, c), n, (a + b + c) / 3, A2 / 2, maximum(len), ν, len)
end

SourceTriangle(a::AbstractVector, b::AbstractVector, c::AbstractVector) =
    SourceTriangle(promote(SVector{3}(a...), SVector{3}(b...), SVector{3}(c...))...)

"""
    edge_log(ri, rj, d)

∫ ds/|x - y| along the straight edge from vertex i to vertex j (length d), from ri = x - vᵢ and rj = x - vⱼ:
ln((rᵢ + rⱼ + d)/(rᵢ + rⱼ - d)). The denominator is formed without cancellation, (rᵢ + rⱼ)² - d² = 2(rᵢrⱼ + ri·rj) and,
when ri·rj < 0, rᵢrⱼ + ri·rj = |ri × rj|²/(rᵢrⱼ - ri·rj), so points next to the edge keep full precision. Infinite on
the edge itself (the log singularity of the tangential field).
"""
@inline function edge_log(ri::SVector{3}, rj::SVector{3}, d)
    a, b = norm(ri), norm(rj)
    p = dot(ri, rj)
    s2d2 = p < 0 ? 2 * sum(abs2, cross(ri, rj)) / (a * b - p) : 2 * (a * b + p)
    return log((a + b + d)^2 / s2d2)
end

"""
    inplane_fraction(tri, x, tol)

For a point in the triangle's plane: the fraction of the full angle (2π) under which it sees the triangle, which sets
the one-sided limit of the normal velocity: 1 inside, 1/2 on an edge, α/2π at a vertex of interior angle α, 0 outside.
"""
function inplane_fraction(tri::SourceTriangle, x::SVector{3}, tol)
    a, b, c = tri.v
    A2 = 2 * tri.area
    λ = (dot(cross(b - x, c - x), tri.n) / A2, dot(cross(c - x, a - x), tri.n) / A2, dot(cross(a - x, b - x), tri.n) / A2)
    any(l -> l < -tol, λ) && return zero(eltype(x))
    nzero = count(l -> abs(l) <= tol, λ)
    nzero == 0 && return one(eltype(x))
    nzero == 1 && return one(eltype(x)) / 2
    k = findfirst(l -> abs(l) > tol, λ) # the vertex the point sits on
    e1, e2 = tri.v[mod1(k + 1, 3)] - tri.v[k], tri.v[mod1(k + 2, 3)] - tri.v[k]
    cosα = dot(e1, e2) / (norm(e1) * norm(e2))
    return acos(clamp(cosα, -one(cosα), one(cosα))) / (2 * oftype(cosα, π))
end

"""
    source_velocity_exact(tri, x, side, plane_tol)

Closed form of w(x) = (1/4π)∫_T (x - y)/|x - y|³ dA (a unit source of the triangle): tangential part
(1/4π)Σₑ νₑ ln((rᵢ + rⱼ + dₑ)/(rᵢ + rⱼ - dₑ)) over the edges (outward in-plane normals νₑ), normal part Ω/(4π) with Ω the
solid angle under which x sees the triangle, positive on the side n points to (Van Oosterom and Strackee). In the
triangle's plane (|z| ≤ plane_tol·diam) the normal part is the one-sided limit `side`·Ω₀/(4π) (`side` = +1 on the n
side, -1 on the other, 0 for the principal value).
"""
function source_velocity_exact(tri::SourceTriangle, x::SVector{3}, side::Integer, plane_tol::Real)
    ra, rb, rc = x - tri.v[1], x - tri.v[2], x - tri.v[3]
    wt = tri.ν[1] * edge_log(ra, rb, tri.len[1]) + tri.ν[2] * edge_log(rb, rc, tri.len[2]) + tri.ν[3] * edge_log(rc, ra, tri.len[3])
    z = dot(x - tri.centroid, tri.n)
    π_ = oftype(z, π) # π in the working number type (Float32, Dual, ...)
    Ω = if abs(z) <= plane_tol * tri.diam
        side * 2 * π_ * inplane_fraction(tri, x, plane_tol)
    else
        la, lb, lc = norm(ra), norm(rb), norm(rc)
        N = dot(ra, cross(rb, rc))
        D = la * lb * lc + dot(ra, rb) * lc + dot(ra, rc) * lb + dot(rb, rc) * la
        2 * atan(N, D)
    end
    return (wt + Ω * tri.n) / (4 * π_)
end

"""
    source_velocity_far(tri, x)

Multipole expansion of the unit source about the centroid up to the octupole (the dipole vanishes). With K = r/r³,
r = x - centroid and dₐ the vertex offsets from the centroid, the exact moments of a triangle,
∫s⊗s dA = (A/12)Σₐdₐ⊗dₐ and ∫s⊗s⊗s dA = (A/30)Σₐdₐ⊗dₐ⊗dₐ, give
w ≈ A/(4π)[K + (1/24)Σₐ(dₐ·∇)²K - (1/180)Σₐ(dₐ·∇)³K], with (d·∇)ⁿK = -∇(d·∇)ⁿ(1/r) from the Legendre form of
(d·∇)ⁿ(1/r). Every term is harmonic. The truncation error is O((diam/r)⁴) relative to the leading term.
"""
function source_velocity_far(tri::SourceTriangle, x::SVector{3})
    r = x - tri.centroid
    R2 = dot(r, r)
    R = sqrt(R2)
    iR2 = 1 / R2
    iR3 = iR2 / R
    iR5 = iR3 * iR2
    iR7 = iR5 * iR2
    iR9 = iR7 * iR2
    w = r * iR3
    for va in tri.v
        d = va - tri.centroid
        dr = dot(d, r)
        dd = dot(d, d)
        quad = 15 * dr^2 * iR7 * r - 6 * dr * iR5 * d - 3 * dd * iR5 * r
        oct = 45 * dr^2 * iR7 * d - 105 * dr^3 * iR9 * r - 9 * dd * iR5 * d + 45 * dd * dr * iR7 * r
        w += quad / 24 - oct / 180
    end
    return tri.area / (4 * oftype(R, π)) * w
end

"""
    source_velocity(tri, x; side=0, far=FAR_DEFAULT, plane_tol=1e-12)

Velocity w(x) = (1/4π)∫_T (x - y)/|x - y|³ dA induced at x by a unit-strength source spread uniformly over the triangle
(a source σ induces σw). Generic over the number type of x and of the triangle (ForwardDiff works).

- Off the triangle's plane: the exact closed form (`source_velocity_exact`), or beyond `far` diameters from the
  centroid the multipole expansion (`source_velocity_far`, `far = Inf` disables it).
- In the triangle's plane (|z| ≤ `plane_tol`·diam): the tangential part is exact (log-singular on the edges and at the
  vertices, which must not be evaluated); the normal part, which jumps by 1 across the triangle, is the limit from the
  side `side` = +1 (where n points), -1, or 0 (principal value): ±1/2 inside, ±1/4 on an edge, ±α/4π at a vertex of
  interior angle α, 0 outside.
"""
function source_velocity(tri::SourceTriangle, x::SVector{3}; side::Integer=0, far::Real=FAR_DEFAULT, plane_tol::Real=1e-12)
    r = x - tri.centroid
    dot(r, r) > (far * tri.diam)^2 && return source_velocity_far(tri, x)
    return source_velocity_exact(tri, x, side, plane_tol)
end

"""
    sheet_velocity(tri, γ, x; side=0, far, plane_tol)

Velocity induced at x by a constant vortex sheet γ (tangent to the triangle) on the triangle: (1/4π)∫ γ × (x - y)/|x - y|³ dA
= γ × w(x). Across the sheet the tangential velocity jumps by γ × n.
"""
sheet_velocity(tri::SourceTriangle, γ::SVector{3}, x::SVector{3}; kwargs...) = cross(γ, source_velocity(tri, x; kwargs...))
