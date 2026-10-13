# Velocities induced by an isoparametric quad's vortex and source sheets: the inner integrals of the Galerkin form,
# numerical since a curved sheet of varying strength has no closed form.

const INV4PI = 1 / (4π)

"""Gauss–Legendre rules on [0, 1] with 1…48 points, tabulated once"""
const GAUSS = [gauss_legendre(n) for n in 1:48]

"""Corners of the parameter square, counterclockwise; side k runs from corner k to corner k + 1"""
const SQUARE = (SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0))

"""
The integrands per unit parameter area at the parameters p, without the factor 1/4π: g_b·(r × nx)/|r|³ for the four
unit nodal sheets (= nx·(g_b × r)/|r|³) and σJ nx·r/|r|³ for the source sheet σ = -V∞·n, where r = x - y(p) and
σJ = -V∞·(a₁ × a₂)
"""
@inline function kernel_terms(Q::IsoQuad, x::SVector{3}, nx::SVector{3}, V∞::SVector{3}, p::SVector{2})
    y = point(Q, p)
    a1, a2 = tangents(Q, p)
    r = x - y
    r2 = dot(r, r)
    ir3 = 1 / (r2 * sqrt(r2))
    t = cross(r, nx) * ir3
    Nξ, Nη = shape_derivatives(p[1], p[2])
    return Nξ * dot(t, a2) - Nη * dot(t, a1), -dot(V∞, cross(a1, a2)) * dot(nx, r) * ir3
end

"""Order of the tensor Gauss rule at `r` diameters from the quad"""
gauss_order(r::Real) = r >= 8 ? 3 : r >= 4 ? 4 : r >= 2 ? 6 : 8

"""
    quad_field(Q, x, nx, V∞; side=0, polar=1)

Normal components along the unit vector `nx` of the velocities the quad Q induces at x: of its vortex sheet per unit
nodal doublet value at each corner, c_b = nx·v_b(x) with v_b = (1/4π)∫ γ_b × (x - y)/|x - y|³ dA, and of its source
sheet σ = -V∞·n, s = nx·(1/4π)∫ σ (x - y)/|x - y|³ dA (zero when V∞ is). Returns (c, s), c::SVector{4}.

At least `polar` diameters from Q a tensor Gauss rule whose order grows as the distance shrinks. Closer, polar
coordinates in the parameter square about the point P of Q closest to x (`polar_field`). For x on Q (inside it), the
principal values with `side=0`, the one-sided limits with `side=±1` (+1 the side the normal points to).
"""
function quad_field(Q::IsoQuad{T}, x::SVector{3,T}, nx::SVector{3,T}, V∞::SVector{3,T}; side::Integer=0, polar::Real=1.0) where {T}
    dc = norm(x - Q.centroid) - Q.reach # a lower bound on the distance
    dc >= polar * Q.diam && return gauss_field(Q, x, nx, V∞, gauss_order(dc / Q.diam))
    P, d = project(Q, x)
    d >= polar * Q.diam && return gauss_field(Q, x, nx, V∞, gauss_order(d / Q.diam))
    return polar_field(Q, x, nx, V∞, P, d, side)
end

"""quad_field by the m×m tensor Gauss rule"""
function gauss_field(Q::IsoQuad{T}, x::SVector{3,T}, nx::SVector{3,T}, V∞::SVector{3,T}, m::Integer) where {T}
    ξs, ws = GAUSS[m]
    c, s = zero(SVector{4,T}), zero(T)
    for i in 1:m, j in 1:m
        dc, ds = kernel_terms(Q, x, nx, V∞, SVector(ξs[i], ξs[j]))
        w = ws[i] * ws[j]
        c += w * dc
        s += w * ds
    end
    return INV4PI * c, INV4PI * s
end

"""
The pieces of [a, b] (split at `brk` when it lies inside) with the number of Gauss points of each, `base` plus
`slope` per unit length, at most 48: in the sinh variables the integrands are smooth and level off along the long
stretches, so the points need to grow only slowly with the length
"""
function pieces(a::Real, b::Real, brk::Real, base::Integer, slope::Real)
    ab = a < brk < b ? ((a, brk), (brk, b)) : ((a, b),)
    return map(((lo, hi),) -> (lo, hi, min(48, base + ceil(Int, slope * (hi - lo)))), ab)
end

"""
    polar_field(Q, x, nx, V∞, P, d, side; angle_points=5, angle_slope=1.5, ray_points=6, ray_slope=1, pv_points=12)

`quad_field` in polar coordinates about the point P of Q closest to x, at distance d.

The parameters are first mapped by z = L(p - P), with L the Cholesky factor of the metric at P (LᵀL = [aᵢ·aⱼ]), so
that |z| is the physical distance to first order and the polar angle is physical too: without it a sliver's sectors
would crowd its whole angular variation into a few degrees. dξdη = dz/J with J = det L = |a₁ × a₂| at P. The square
maps to a parallelogram, which gives one sector per side (a side through P has none), integrated with sinh
substitutions and Gauss rules whose points grow slowly with the substituted length (`pieces`):
- across each sector u = h sinh v along the side (h the distance from P to the side's line, u measured from the foot
  of the perpendicular), so a sector whose apex lies next to its side stays smooth in v: dθ = dv/cosh v; split at the
  foot, `angle_points` + `angle_slope` per unit v;
- along each ray ρ = ρ₀ sinh s, with ρ₀ = d/|A(θ)| the near-singular scale (A(θ) the physical tangent of the ray at
  P, of unit length to first order), `ray_points` + `ray_slope` per unit s; plain Gauss (`pv_points`) when the ray
  is shorter than ρ₀.
For x on Q (the self pair) the sheet terms are principal values: each ray's 1/ρ term is subtracted and added back
analytically (Guiggiani), PV = Σ∫dθ [∫₀ᴿ (F - F₋₁/ρ) dρ + F₋₁ ln(R|A|)] with F₋₁(θ) = -nx·(g₀ × A)/(J|A|³) and g₀ the
sheet density at P; the logarithms of the excluded radius cancel over the full circle since F₋₁ is odd in θ. The
source's own field is regular when nx is the normal at P, and its F₋₁ is subtracted the same way otherwise. `side`
≠ 0 then gives the one-sided limit (+1: the side the normal points to), adding half the jumps across the layer:
±(γ × n)/2 for the sheets (normal to n, so nothing when nx = n) and ±σn/2 for the source.
"""
function polar_field(Q::IsoQuad{T}, x::SVector{3,T}, nx::SVector{3,T}, V∞::SVector{3,T}, P::SVector{2,T}, d::T, side::Integer;
    angle_points::Integer=5, angle_slope::Real=1.5, ray_points::Integer=6, ray_slope::Real=1.0, pv_points::Integer=12) where {T}
    onq = d <= 1e-12 * Q.diam
    onq && minimum((P[1], P[2], 1 - P[1], 1 - P[2])) <= 1e-12 &&
        throw(ArgumentError("quad_field: x lies on the quad's boundary, where the sheet's field is singular"))
    a1s, a2s = tangents(Q, P)
    g11, g12, g22 = dot(a1s, a1s), dot(a1s, a2s), dot(a2s, a2s)
    l11 = sqrt(g11)
    l12, l22 = g12 / l11, sqrt(g22 - g12^2 / g11)
    J = l11 * l22
    toz(p) = SVector(l11 * p[1] + l12 * p[2], l22 * p[2]) # z = L p
    top(z) = SVector((z[1] - l12 * z[2] / l22) / l11, z[2] / l22) # p = L⁻¹ z
    g0 = sheet_densities(Q, P)
    σJ0 = -dot(V∞, cross(a1s, a2s))
    c, s = zero(SVector{4,T}), zero(T)
    for k in 1:4
        A, B = toz(SQUARE[k] - P), toz(SQUARE[mod1(k + 1, 4)] - P)
        L = norm(B - A)
        et = (B - A) / L
        nk = SVector(et[2], -et[1]) # outward normal of side k (L keeps the orientation)
        h = dot(A, nk)
        h <= 1e-14 * L && continue
        for (va, vb, nv) in pieces(asinh(dot(A, et) / h), asinh(dot(B, et) / h), 0.0, angle_points, angle_slope)
            vs, wvs = GAUSS[nv]
            for (vq, wq) in zip(vs, wvs)
                v = va + (vb - va) * vq
                ch = cosh(v)
                R = h * ch
                e = top((nk + sinh(v) * et) / ch) # parameter step per unit z along the ray
                Aθ = e[1] * a1s + e[2] * a2s
                nA = norm(Aθ)
                W = wq * (vb - va) / (ch * J) # dθ / J
                if onq
                    f = SVector(ntuple(b -> -dot(nx, cross(g0[b], Aθ)) / nA^3, 4))
                    fs = -σJ0 * dot(nx, Aθ) / nA^3 # zero when nx is the normal at P
                    dc, ds = ray_pv(Q, x, nx, V∞, P, e, R, f, fs, pv_points)
                    dc += f * log(R * nA)
                    ds += fs * log(R * nA)
                else
                    dc, ds = ray_near(Q, x, nx, V∞, P, e, R, d / nA, ray_points, ray_slope, pv_points)
                end
                c += W * dc
                s += W * ds
            end
        end
    end
    c, s = INV4PI * c, INV4PI * s
    if onq && side != 0 # half the jumps toward the side n points to: γ × n for the sheets, σ n for the source
        nP = unit_normal(Q, P)
        c += side / 2 * SVector(ntuple(b -> dot(nx, cross(g0[b], nP)) / J, 4))
        s += side * (-dot(V∞, nP)) / 2 * dot(nx, nP)
    end
    return c, s
end

"""
∫₀ᴿ (F - F₋₁/ρ) dρ along a ray from P (on Q), F = ρ × integrand, for the sheets (F₋₁ = f) and the source (F₋₁ = fs,
which is -σJ nx·A/|A|³ and vanishes when nx is the normal at P: the source's own field is then regular)
"""
function ray_pv(Q::IsoQuad{T}, x, nx, V∞, P, e, R, f, fs, n::Integer) where {T}
    ρs, wρ = GAUSS[n]
    c, s = zero(SVector{4,T}), zero(T)
    for (ρq, wq) in zip(ρs, wρ)
        ρ = R * ρq
        dc, ds = kernel_terms(Q, x, nx, V∞, P + ρ * e)
        c += (wq * R) * (ρ * dc - f / ρ)
        s += (wq * R) * (ρ * ds - fs / ρ)
    end
    return c, s
end

"""∫₀ᴿ F dρ along a ray from P (x off Q), with ρ = ρ₀ sinh s when the near-singular scale ρ₀ is below R"""
function ray_near(Q::IsoQuad{T}, x, nx, V∞, P, e, R, ρ0, npts, slope, nplain) where {T}
    c, s = zero(SVector{4,T}), zero(T)
    if ρ0 < R
        for (sa, sb, n) in pieces(0.0, asinh(R / ρ0), -1.0, npts, slope)
            ss, ws = GAUSS[n]
            for (sq, wq) in zip(ss, ws)
                σv = sa + (sb - sa) * sq
                ρ = ρ0 * sinh(σv)
                jac = wq * (sb - sa) * ρ0 * cosh(σv) * ρ
                dc, ds = kernel_terms(Q, x, nx, V∞, P + ρ * e)
                c += jac * dc
                s += jac * ds
            end
        end
    else
        ρs, wρ = GAUSS[nplain]
        for (ρq, wq) in zip(ρs, wρ)
            ρ = R * ρq
            dc, ds = kernel_terms(Q, x, nx, V∞, P + ρ * e)
            c += (wq * R * ρ) * dc
            s += (wq * R * ρ) * ds
        end
    end
    return c, s
end
