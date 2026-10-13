# Exact potential flow about a fixed ellipsoid in uniform flow (Lamb, Hydrodynamics, §114), the sphere being the special
# case a = b = c. Shared by the validation scripts and the tests.

using LinearAlgebra, StaticArrays

"""
    carlson_rd(x, y, z)

Carlson's symmetric elliptic integral R_D(x, y, z) = (3/2)∫₀^∞ dt / ((t+z)√((t+x)(t+y)(t+z))), for x, y ≥ 0 (at most
one of them zero) and z > 0, by the duplication algorithm of Carlson (1995, Numer. Algorithms 10), to rounding.
"""
function carlson_rd(x::Real, y::Real, z::Real)
    x0, y0, z0 = float(x), float(y), float(z)
    A0 = (x0 + y0 + 3z0) / 5
    Q = (eps(A0) / 4)^(-1 / 6) * max(abs(A0 - x0), abs(A0 - y0), abs(A0 - z0))
    x, y, z, A = x0, y0, z0, A0
    s, f = zero(A0), one(A0) # the sum and 4^-m
    while f * Q >= abs(A)
        λ = sqrt(x) * sqrt(y) + sqrt(x) * sqrt(z) + sqrt(y) * sqrt(z)
        s += f / (sqrt(z) * (z + λ))
        f /= 4
        x, y, z, A = (x + λ) / 4, (y + λ) / 4, (z + λ) / 4, (A + λ) / 4
    end
    X, Y = f * (A0 - x0) / A, f * (A0 - y0) / A
    Z = -(X + Y) / 3
    # (no juxtaposed coefficients below: 3E2 would parse as the literal 300.0)
    e2, e3 = X * Y - 6 * Z^2, (3 * X * Y - 8 * Z^2) * Z
    e4, e5 = 3 * (X * Y - Z^2) * Z^2, X * Y * Z^3
    series = 1 - 3 * e2 / 14 + e3 / 6 + 9 * e2^2 / 88 - 3 * e4 / 22 - 9 * e2 * e3 / 52 + 3 * e5 / 26
    return f * A^(-3 / 2) * series + 3 * s
end

"""
    ellipsoid_factors(a, b, c)

Depolarization factors (α₀, β₀, γ₀) of the ellipsoid with semi-axes a, b, c: α₀ = abc∫₀^∞ dλ/((a²+λ)Δ(λ)) with
Δ = √((a²+λ)(b²+λ)(c²+λ)), and likewise for b and c. They sum to 2 (2/3 each for a sphere).
"""
function ellipsoid_factors(a::Real, b::Real, c::Real)
    k = 2 * a * b * c / 3
    return SVector(k * carlson_rd(b^2, c^2, a^2), k * carlson_rd(c^2, a^2, b^2), k * carlson_rd(a^2, b^2, c^2))
end

"""
    EllipsoidFlow(axes, V∞)

Uniform flow V∞ about the fixed ellipsoid Σ xᵢ²/aᵢ² = 1. On the surface the perturbation potential is linear,
φ = Σ αᵢ/(2 - αᵢ) V∞ᵢ xᵢ = k_μ·x (αᵢ the depolarization factors). That is the doublet μ of the Galerkin formulation,
whose inner perturbation potential is zero (σ = -V∞·n). The surface velocity is the tangential part of the constant
vector k = V∞ + k_μ, kᵢ = 2V∞ᵢ/(2 - αᵢ), so Cp = 1 - |k - (k·n)n|²/|V∞|² follows the curvature. A sphere has
αᵢ = 2/3: μ = V∞·x/2 and k = 3V∞/2.
"""
struct EllipsoidFlow
    axes::SVector{3,Float64}
    V∞::SVector{3,Float64}
    factors::SVector{3,Float64}
    kμ::SVector{3,Float64}
    k::SVector{3,Float64}
end

function EllipsoidFlow(axes, V∞)
    ax, V = SVector{3,Float64}(axes), SVector{3,Float64}(V∞)
    α = ellipsoid_factors(ax...)
    return EllipsoidFlow(ax, V, α, α ./ (2 .- α) .* V, 2 .* V ./ (2 .- α))
end

"""Point where the ray from the centre through p meets the surface"""
surface_point(e::EllipsoidFlow, p::SVector{3}) = p / sqrt(sum(abs2, p ./ e.axes))

"""Outward unit normal at the surface point x"""
surface_normal(e::EllipsoidFlow, x::SVector{3}) = normalize(x ./ e.axes .^ 2)

"""Exact doublet μ (the outer perturbation potential) at a surface point; linear, so also at the nodes of constraints"""
exact_doublet(e::EllipsoidFlow, x::SVector{3}) = dot(e.kμ, x)

"""Largest |μ| on the surface, the scale of the doublet errors"""
doublet_scale(e::EllipsoidFlow) = norm(e.kμ .* e.axes)

"""Exact surface velocity where the ray through p meets the surface"""
function exact_velocity(e::EllipsoidFlow, p::SVector{3})
    n = surface_normal(e, surface_point(e, p))
    return e.k - dot(e.k, n) * n
end

"""Exact pressure coefficient where the ray through p meets the surface"""
exact_cp(e::EllipsoidFlow, p::SVector{3}) = 1 - sum(abs2, exact_velocity(e, p)) / sum(abs2, e.V∞)

"""The ellipsoidal coordinate of x: the largest root λ of Σ xᵢ²/(aᵢ² + λ) = 1, zero on and inside the surface"""
function ellipsoidal_coordinate(e::EllipsoidFlow, x::SVector{3})
    a2 = e.axes .^ 2
    F(λ) = sum(x .^ 2 ./ (a2 .+ λ)) - 1
    F(0.0) <= 0 && return 0.0
    lo, hi = 0.0, max(1.0, sum(abs2, x)) # F(hi) < 0
    for _ in 1:200
        mid = (lo + hi) / 2
        F(mid) > 0 ? (lo = mid) : (hi = mid)
        hi - lo <= 2eps(hi) && break
    end
    return (lo + hi) / 2
end

"""
    perturbation_potential(e, x)

Perturbation potential at a point x on or outside the ellipsoid, Σ Aᵢ xᵢ Iᵢ(λ) with Iᵢ(λ) = ∫_λ^∞ dλ'/((aᵢ²+λ')Δ(λ'))
(Carlson's R_D), λ the ellipsoidal coordinate of x and Aᵢ = abc V∞ᵢ/(2 - αᵢ). On the surface it is k_μ·x; it is
harmonic outside, decays like r⁻², and meets ∂φ/∂n = -V∞·n on the surface (the checks of the surface formulas).
"""
function perturbation_potential(e::EllipsoidFlow, x::SVector{3})
    λ = ellipsoidal_coordinate(e, x)
    a2 = e.axes .^ 2 .+ λ
    I = 2 / 3 .* SVector(carlson_rd(a2[2], a2[3], a2[1]), carlson_rd(a2[3], a2[1], a2[2]), carlson_rd(a2[1], a2[2], a2[3]))
    return sum(prod(e.axes) .* e.V∞ ./ (2 .- e.factors) .* x .* I)
end
