# Independent references for the source triangle kernel w(x) = (1/4π)∫_T (x - y)/|x - y|³ dA.
#
# 1. Polar quadrature about the projection p of x on the triangle's plane: the radial integrals are done analytically
#    and the angle adaptively (QuadGK), split at the vertex directions where the integrand has kinks. Works at any height
#    z, including the plane (principal value for p inside the triangle, one-sided normal limit chosen by `side`).
# 2. Plain adaptive cubature over the triangle (HCubature on the Duffy square), for points off the plane.

using QuadGK
using HCubature

"""[ρ_in, ρ_out] of the ray ρ ↦ ρe (ρ ≥ 0, from the origin) inside the 2D triangle q (counterclockwise), or nothing"""
function ray_interval(q, e::SVector{2,Float64})
    lo, hi = 0.0, Inf
    for k in 1:3
        a, b = q[k], q[mod1(k + 1, 3)]
        t = b - a
        ν = SVector(t[2], -t[1]) # outward
        c0, c1 = -dot(ν, a), dot(ν, e) # inside: c0 + ρ c1 ≤ 0
        if abs(c1) < 1e-300
            c0 > 0 && return nothing
        elseif c1 > 0
            hi = min(hi, -c0 / c1)
        else
            lo = max(lo, -c0 / c1)
        end
    end
    return lo < hi ? (lo, hi) : nothing
end

"""Polar-quadrature reference of w(x) for the triangle (a, b, c) (counterclockwise about n), see the file header"""
function w_polar(a::SVector{3,Float64}, b::SVector{3,Float64}, c::SVector{3,Float64}, x::SVector{3,Float64}; side::Integer=0, rtol=1e-13)
    n = normalize(cross(b - a, c - a))
    e1 = normalize(b - a)
    e2 = cross(n, e1)
    z = dot(x - a, n)
    size = max(norm(b - a), norm(c - b), norm(a - c))
    abs(z) <= 1e-12 * size && (z = 0.0) # in the plane up to rounding: the side is chosen by `side`
    p = x - z * n
    loc(y) = SVector(dot(y - p, e1), dot(y - p, e2))
    q = (loc(a), loc(b), loc(c))
    Z = abs(z)
    function integrand(θ)
        e = SVector(cos(θ), sin(θ))
        iv = ray_interval(q, e)
        isnothing(iv) && return zero(SVector{3,Float64})
        ρ1, ρ2 = iv
        if Z > 0
            F(ρ) = asinh(ρ / Z) - ρ / sqrt(ρ^2 + Z^2)
            tang = -(F(ρ2) - F(ρ1))
            normal = z * (1 / sqrt(ρ1^2 + Z^2) - 1 / sqrt(ρ2^2 + Z^2))
        else
            # in the plane: from inside, the principal value drops the ln(ρ1 → 0) part, whose angular integral vanishes
            tang = ρ1 == 0 ? -log(ρ2) : -log(ρ2 / ρ1)
            normal = ρ1 == 0 ? Float64(side) : 0.0
        end
        return SVector(tang * e[1], tang * e[2], normal)
    end
    breaks = sort(unique([0.0; [mod(atan(v[2], v[1]), 2π) for v in q]; 2π]))
    I = sum(quadgk(integrand, breaks[k], breaks[k+1]; rtol=rtol, maxevals=10^6)[1] for k in 1:length(breaks)-1)
    return (I[1] * e1 + I[2] * e2 + I[3] * n) / 4π
end

"""Adaptive cubature reference of w(x) over the Duffy square, for x off the triangle"""
function w_cubature(a::SVector{3,Float64}, b::SVector{3,Float64}, c::SVector{3,Float64}, x::SVector{3,Float64}; rtol=1e-11)
    A2 = norm(cross(b - a, c - a))
    f(uv) = begin
        u, v = uv
        y = a + u * (b - a) + u * v * (c - b)
        r = x - y
        (A2 * u) * r / norm(r)^3
    end
    I, _ = hcubature(f, SVector(0.0, 0.0), SVector(1.0, 1.0); rtol=rtol, maxevals=5 * 10^6)
    return I / 4π
end
