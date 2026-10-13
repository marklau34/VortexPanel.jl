# Load post-processing shared by the MESA baseline and, later, the VortexPanel solver. Moves into the Post module once
# that exists. Requires LinearAlgebra and StaticArrays.

"""Reference quantities for force and moment coefficients"""
Base.@kwdef struct WingReference{T<:Real}
    span::T
    chord::T
    area::T = span*chord
    rho::T
    vinf::T
    qinf::T = rho*vinf^2/2
end

"""
    fit_rotation(xb, xg)

Rotation R with xg ≈ R*xb (Kabsch, frames sharing their origin), e.g. from body frame mesh nodes to the same nodes after
a solver rotated them into the global frame.
"""
function fit_rotation(xb::AbstractVector{<:AbstractVector}, xg::AbstractVector{<:AbstractVector})
    H = sum(SVector{3}(b) * SVector{3}(g)' for (b, g) in zip(xb, xg))
    F = svd(Matrix(H))
    d = sign(det(F.V * F.U'))
    return SMatrix{3,3}(F.V * Diagonal([1, 1, d]) * F.U')
end

"""Angle of attack of a body pitched by R about y (body x from leading to trailing edge, freestream along global +x)"""
pitch_angle(R::AbstractMatrix) = atan(-R[3, 1], R[1, 1])

"""
    force_coefficients(F, M, ref; r_c4)

Coefficients from the total force `F` and the moment `M` about the global origin, both in the wind frame (freestream
along +x, lift along +z). `r_c4` is the moment reference point for `Cm_c4`, in the same frame.
"""
function force_coefficients(F::AbstractVector, M::AbstractVector, ref::WingReference; r_c4::AbstractVector=SA[0.0, 0.0, 0.0])
    qS = ref.qinf * ref.area
    M_c4 = SVector{3}(M) - cross(SVector{3}(r_c4), SVector{3}(F))
    return (CL=F[3]/qS, CDi=F[1]/qS, CY=F[2]/qS, Cm_LE=M[2]/(qS*ref.chord), Cm_c4=M_c4[2]/(qS*ref.chord))
end

"""Lift from the shed circulation, L = ρ V∞ Σ Γ_k Δy_k, with Γ positive for positive lift"""
kj_lift(Γ::AbstractVector, dy::AbstractVector, ref::WingReference) = ref.rho * ref.vinf * sum(Γ .* dy)

"""
    trefftz_drag(y_edges, Γ, ref)

Induced drag of a planar wake made of strips [y_edges[k], y_edges[k+1]] with constant circulation Γ[k] (positive for
positive lift). The trailing vortices sit on the strip edges with strengths ΔΓ = Γ[k] - Γ[k-1]; the downwash
w(y) = 1/(4π) Σ ΔΓ_m/(y - y_m) at the lifting line (half the Trefftz plane value) is taken at the strip midpoints, and
D = ρ Σ Γ_k w_k Δy_k. An elliptic loading gives CDi = CL²/(π AR) in the limit of many strips.
"""
function trefftz_drag(y_edges::AbstractVector, Γ::AbstractVector, ref::WingReference)
    n = length(Γ)
    @assert length(y_edges) == n + 1
    ΔΓ = [k == 1 ? Γ[1] : k == n + 1 ? -Γ[n] : Γ[k] - Γ[k-1] for k in 1:n+1]
    D = zero(promote_type(eltype(Γ), eltype(y_edges)))
    for k in 1:n
        yc = (y_edges[k] + y_edges[k+1]) / 2
        w = sum(ΔΓ[m] / (yc - y_edges[m]) for m in 1:n+1) / 4π
        D += Γ[k] * w * (y_edges[k+1] - y_edges[k])
    end
    return ref.rho * D
end

"""
    richardson(vals)

Observed order p and extrapolated value from three values on meshes refined by 2 each time; (NaN, NaN) when the
differences change sign.
"""
function richardson(vals::AbstractVector)
    ratio = (vals[1] - vals[2]) / (vals[2] - vals[3])
    ratio <= 0 && return (NaN, NaN)
    p = log2(ratio)
    return (p, vals[3] + (vals[3] - vals[2]) / (2^p - 1))
end

"""
    extrapolate(vals; min_order=0.3)

`richardson`, with the extrapolated value NaN unless the observed order is at least `min_order`: below it the values
are not in the asymptotic range, and 1/(2^p - 1) amplifies their differences more than threefold.
"""
function extrapolate(vals::AbstractVector; min_order::Real=0.3)
    p, e = richardson(vals)
    return p, (isfinite(p) && p >= min_order ? e : NaN)
end

"""
    structured_strip(e, nspan, nchord)

Spanwise strip of element `e` of a PAGE wing mesh (upper surface elements first, chordwise index fastest, strips in
+y), or 0 for elements appended after the wing (tip caps).
"""
function structured_strip(e::Integer, nspan::Integer, nchord::Integer)
    nsurf = nspan * nchord
    e > 2nsurf && return 0
    return mod(e - 1, nsurf) ÷ nchord + 1
end

"""Chordwise index (1 at the leading edge) and surface (+1 upper, -1 lower) of element `e` of a PAGE wing mesh"""
function structured_chord(e::Integer, nspan::Integer, nchord::Integer)
    nsurf = nspan * nchord
    return mod(e - 1, nchord) + 1, e <= nsurf ? 1 : -1
end

"""
    interp_stations(η_strips, values, η)

Linear interpolation of strip values (rows of `values` per strip) to the spanwise station η, clamped at the ends.
"""
function interp_stations(η_strips::AbstractVector, values::AbstractMatrix, η::Real)
    k = searchsortedlast(η_strips, η)
    k <= 0 && return values[1, :]
    k >= length(η_strips) && return values[end, :]
    w = (η - η_strips[k]) / (η_strips[k+1] - η_strips[k])
    return (1 - w) .* values[k, :] .+ w .* values[k+1, :]
end
