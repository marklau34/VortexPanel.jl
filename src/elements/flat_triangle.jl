"""Element interface: `basis`, `point`, `sheet_strength`, `gradient`, `sheet_velocity` (an isoparametric quad later)"""
abstract type AbstractElement end

"""
    FlatTriangle(x₁, x₂, x₃)

Flat triangle element with nodal hat functions φ_a = λ_a (barycentric coordinates). A unit nodal value at corner b
carries the constant vortex sheet γ_b = n × ∇λ_b = -e_b/(2A), with e_b = x_{b+2} - x_{b+1} the edge opposite b,
counterclockwise about the outward normal n (∇λ_b = n × e_b/(2A) points from that edge toward b). `src` is the
triangle's source kernel data.
"""
struct FlatTriangle{T} <: AbstractElement
    x::SVector{3,SVector{3,T}}
    n::SVector{3,T}
    area::T
    centroid::SVector{3,T}
    diam::T
    gradλ::SVector{3,SVector{3,T}}
    γ::SVector{3,SVector{3,T}}
    src::SourceTriangle{T}
end

function FlatTriangle(a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}) where {T}
    src = SourceTriangle(a, b, c)
    n, A2 = src.n, 2 * src.area
    e = (c - b, a - c, b - a)
    gradλ = SVector(cross(n, e[1]) / A2, cross(n, e[2]) / A2, cross(n, e[3]) / A2)
    γ = SVector(-e[1] / A2, -e[2] / A2, -e[3] / A2)
    return FlatTriangle(SVector(a, b, c), n, src.area, src.centroid, src.diam, gradλ, γ, src)
end

"""Hat function values at barycentric coordinates λ (the coordinates themselves)"""
basis(::FlatTriangle, λ::SVector{3}) = λ

"""Point at barycentric coordinates λ"""
point(el::FlatTriangle, λ::SVector{3}) = λ[1] * el.x[1] + λ[2] * el.x[2] + λ[3] * el.x[3]

"""Vortex sheet γ = n × ∇μ of the nodal doublet values μ (one per corner)"""
sheet_strength(el::FlatTriangle, μ::SVector{3}) = μ[1] * el.γ[1] + μ[2] * el.γ[2] + μ[3] * el.γ[3]

"""Surface gradient ∇μ of the nodal doublet values μ (= γ × n)"""
gradient(el::FlatTriangle, μ::SVector{3}) = μ[1] * el.gradλ[1] + μ[2] * el.gradλ[2] + μ[3] * el.gradλ[3]

"""
    sheet_velocity(el::FlatTriangle, μ, x; side, far)

Velocity induced at x by the element's vortex sheet γ = n × ∇μ for the nodal doublet values μ (one per corner):
γ × w(x). (A method of the kernel's `sheet_velocity`, which takes the sheet γ of a source triangle directly.)
"""
Kernels.sheet_velocity(el::FlatTriangle, μ::SVector{3}, x::SVector{3}; kwargs...) = cross(sheet_strength(el, μ), source_velocity(el.src, x; kwargs...))

"""Flat triangle elements of a geometry, in triangle order"""
build_elements(g::SurfaceGeometry) = [FlatTriangle(g.x[t[1]], g.x[t[2]], g.x[t[3]]) for t in g.tris]

"""Distance from point p to the triangle (closest point after Ericson, Real-Time Collision Detection, 5.1.5)"""
function distance_to(el::FlatTriangle, p::SVector{3})
    a, b, c = el.x
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    (d1 <= 0 && d2 <= 0) && return norm(ap)
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    (d3 >= 0 && d4 <= d3) && return norm(bp)
    vc = d1 * d4 - d3 * d2
    (vc <= 0 && d1 >= 0 && d3 <= 0) && return norm(p - (a + d1 / (d1 - d3) * ab))
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    (d6 >= 0 && d5 <= d6) && return norm(cp)
    vb = d5 * d2 - d1 * d6
    (vb <= 0 && d2 >= 0 && d6 <= 0) && return norm(p - (a + d2 / (d2 - d6) * ac))
    va = d3 * d6 - d5 * d4
    (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) && return norm(p - (b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)))
    denom = 1 / (va + vb + vc)
    return norm(p - (a + ab * (vb * denom) + ac * (vc * denom)))
end
