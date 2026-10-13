"""
    IsoQuad(x₁, x₂, x₃, x₄)

Isoparametric bilinear quadrilateral (doc §6, the second element type) with corners x₁…x₄ counterclockwise about the
outward normal: the surface x(ξ, η) = Σ N_a x_a over the parameter square (ξ, η) ∈ [0, 1]², with the bilinear shape
functions N₁ = (1-ξ)(1-η), N₂ = ξ(1-η), N₃ = ξη, N₄ = (1-ξ)η, and the doublet μ = Σ N_a μ_a on the same functions. The
surface is a hyperbolic paraboloid with straight edges (a plane when the corners are coplanar), so neighbouring
elements meet along their edges exactly. With the tangents a₁ = ∂x/∂ξ and a₂ = ∂x/∂η, dA = |a₁ × a₂| dξ dη, and the
vortex sheet γ = n × ∇μ satisfies γ dA = (μ_ξ a₂ - μ_η a₁) dξ dη, so the sheet's density per unit parameter area is a
polynomial: `sheet_densities`. `warp` is ∂²x/∂ξ∂η = x₁ - x₂ + x₃ - x₄ (zero for a parallelogram), `centroid` the point
x(½, ½), `n` the unit normal there, `diam` the longest edge or diagonal, `reach` the largest distance from the centroid
to a corner.
"""
struct IsoQuad{T} <: AbstractElement
    x::SVector{4,SVector{3,T}}
    warp::SVector{3,T}
    centroid::SVector{3,T}
    n::SVector{3,T}
    area::T
    diam::T
    reach::T
end

function IsoQuad(a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}, d::SVector{3,T}) where {T}
    x = SVector(a, b, c, d)
    nc = cross(c - a, d - b) # 2 a₁ × a₂ at the centre
    norm(nc) > 0 || throw(ArgumentError("degenerate quadrilateral $a, $b, $c, $d"))
    ξs, ws = gauss_legendre(4)
    area = zero(T)
    for (ξ, wξ) in zip(ξs, ws), (η, wη) in zip(ξs, ws)
        a1 = (1 - η) * (b - a) + η * (c - d)
        a2 = (1 - ξ) * (d - a) + ξ * (c - b)
        area += wξ * wη * norm(cross(a1, a2))
    end
    diam = max(norm(b - a), norm(c - b), norm(d - c), norm(a - d), norm(c - a), norm(d - b))
    centroid = (a + b + c + d) / 4
    return IsoQuad(x, a - b + c - d, centroid, nc / norm(nc), area, diam, maximum(norm(p - centroid) for p in x))
end

"""Bilinear shape functions N₁…N₄ at (ξ, η)"""
@inline shape(ξ, η) = SVector((1 - ξ) * (1 - η), ξ * (1 - η), ξ * η, (1 - ξ) * η)

"""∂N_a/∂ξ and ∂N_a/∂η at (ξ, η)"""
@inline shape_derivatives(ξ, η) = (SVector(-(1 - η), 1 - η, η, -η), SVector(-(1 - ξ), -ξ, ξ, 1 - ξ))

"""Hat (shape) function values at the parameters p = (ξ, η)"""
basis(::IsoQuad, p::SVector{2}) = shape(p[1], p[2])

"""Point at the parameters p = (ξ, η)"""
@inline function point(Q::IsoQuad, p::SVector{2})
    N = shape(p[1], p[2])
    return N[1] * Q.x[1] + N[2] * Q.x[2] + N[3] * Q.x[3] + N[4] * Q.x[4]
end

"""Tangents a₁ = ∂x/∂ξ and a₂ = ∂x/∂η at the parameters p"""
@inline function tangents(Q::IsoQuad, p::SVector{2})
    ξ, η = p
    a1 = (1 - η) * (Q.x[2] - Q.x[1]) + η * (Q.x[3] - Q.x[4])
    a2 = (1 - ξ) * (Q.x[4] - Q.x[1]) + ξ * (Q.x[3] - Q.x[2])
    return a1, a2
end

"""a₁ × a₂ at the parameters p: the outward normal times the area element (linear in ξ and η)"""
@inline area_normal(Q::IsoQuad, p::SVector{2}) = cross(tangents(Q, p)...)

"""Unit normal at the parameters p"""
@inline unit_normal(Q::IsoQuad, p::SVector{2}) = normalize(area_normal(Q, p))

"""Vortex sheet densities per unit parameter area of the four unit nodal doublet values at p: N_b,ξ a₂ - N_b,η a₁"""
@inline function sheet_densities(Q::IsoQuad, p::SVector{2})
    a1, a2 = tangents(Q, p)
    Nξ, Nη = shape_derivatives(p[1], p[2])
    return SVector(ntuple(b -> Nξ[b] * a2 - Nη[b] * a1, 4))
end

"""Vortex sheet γ = n × ∇μ at p for the nodal doublet values μ"""
function sheet_strength(Q::IsoQuad, μ::SVector{4}, p::SVector{2})
    a1, a2 = tangents(Q, p)
    Nξ, Nη = shape_derivatives(p[1], p[2])
    return (dot(Nξ, μ) * a2 - dot(Nη, μ) * a1) / norm(cross(a1, a2))
end

"""Surface gradient ∇μ at p of the nodal doublet values μ: μ_ξ a¹ + μ_η a² with a¹ = a₂ × n/J, a² = n × a₁/J"""
function gradient(Q::IsoQuad, μ::SVector{4}, p::SVector{2})
    a1, a2 = tangents(Q, p)
    nJ = cross(a1, a2)
    J2 = dot(nJ, nJ)
    Nξ, Nη = shape_derivatives(p[1], p[2])
    return (dot(Nξ, μ) * cross(a2, nJ) + dot(Nη, μ) * cross(nJ, a1)) / J2
end

"""Closest point of the segment from a to b to p: its parameter t ∈ [0, 1] and the distance"""
@inline function segment_closest(a::SVector{3}, b::SVector{3}, p::SVector{3})
    e = b - a
    t = clamp(dot(p - a, e) / dot(e, e), zero(eltype(e)), one(eltype(e)))
    return t, norm(p - (a + t * e))
end

"""
    project(Q, p)

The point of Q closest to p: its parameters (ξ, η) and the distance. Newton's method on |x(ξ, η) - p|², from the
nearest point of a 5×5 parameter grid, with each step cut short to stay within [-½, 3/2]² and then halved until the
distance falls, gives an interior candidate (on a sheared sliver the first step can overshoot far); the four straight
edges give the boundary ones; the nearest wins.
"""
function project(Q::IsoQuad{T}, p::SVector{3}) where {T}
    best_d, best = T(Inf), SVector(T(0.5), T(0.5))
    q = best
    f0 = T(Inf)
    for a in 0:4, b in 0:4
        g = SVector(T(a / 4), T(b / 4))
        r = point(Q, g) - p
        dot(r, r) < f0 && ((f0, q) = (dot(r, r), g))
    end
    for _ in 1:60
        r = point(Q, q) - p
        a1, a2 = tangents(Q, q)
        g1, g2 = dot(r, a1), dot(r, a2)
        H11, H22, H12 = dot(a1, a1), dot(a2, a2), dot(a1, a2) + dot(r, Q.warp)
        det = H11 * H22 - H12^2
        step = det > 0 ? SVector((H22 * g1 - H12 * g2) / det, (H11 * g2 - H12 * g1) / det) : SVector(g1, g2) / (H11 + H22)
        t = one(T)
        for i in 1:2 # stay within [-½, 3/2]²
            step[i] > 0 && (t = min(t, (q[i] + T(0.5)) / step[i]))
            step[i] < 0 && (t = min(t, (q[i] - T(1.5)) / step[i]))
        end
        f = dot(r, r)
        while t > 1e-10
            rn = point(Q, q - t * step) - p
            dot(rn, rn) <= f && break
            t /= 2
        end
        q -= t * step
        t * (abs(step[1]) + abs(step[2])) < 1e-14 && break
    end
    if 0 <= q[1] <= 1 && 0 <= q[2] <= 1
        best_d, best = norm(point(Q, q) - p), q
    end
    for k in 1:4
        t, d = segment_closest(Q.x[k], Q.x[mod1(k + 1, 4)], p)
        if d < best_d
            best_d = d
            best = k == 1 ? SVector(t, zero(T)) : k == 2 ? SVector(one(T), t) : k == 3 ? SVector(1 - t, one(T)) : SVector(zero(T), 1 - t)
        end
    end
    return best, best_d
end

distance_to(Q::IsoQuad, p::SVector{3}) = last(project(Q, p))

"""Largest distance from the centroid to a corner"""
reach(Q::IsoQuad) = Q.reach
