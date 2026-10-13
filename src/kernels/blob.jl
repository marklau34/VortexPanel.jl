"""
    Blob(x, α, δ)

A vortex particle at x with strength α (the vector circulation ∫ω dV it carries) and core radius δ: the free wake's
particles, and the blob representation of wake sheets and filaments where velocities are needed at wake points.
"""
struct Blob{T}
    x::SVector{3,T}
    α::SVector{3,T}
    δ::T
end

"""
    blob_velocity(x, xp, α, δ)
    blob_velocity(b::Blob, x)

Velocity at x of a vortex particle at xp of strength α and core δ, with the Rosenhead–Moore kernel (as MESA):
(1/4π) α × r / (|r|² + δ²)^{3/2}, r = x - xp. It is the Biot–Savart law of a point vortex where |r| ≫ δ, vanishes at
the particle and is solenoidal everywhere.
"""
@inline function blob_velocity(x::SVector{3}, xp::SVector{3}, α::SVector{3}, δ::Real)
    r = x - xp
    d2 = dot(r, r) + δ^2
    return cross(α, r) / (4 * oftype(d2, π) * d2 * sqrt(d2))
end

@inline blob_velocity(b::Blob, x::SVector{3}) = blob_velocity(x, b.x, b.α, b.δ)

"""
Two vortex particles standing for the straight filament from a to b with circulation Γa at a and Γb at b, at its
2-point Gauss positions: they carry its total strength ∫Γ dl and its first moment ∫y ⊗ Γ dl exactly
"""
function segment_blobs(a::SVector{3,T}, b::SVector{3,T}, Γa::Real, Γb::Real, δ::Real) where {T}
    s = T(1 / (2 * sqrt(3)))
    return map((T(0.5) - s, T(0.5) + s)) do ξ
        Blob(a + ξ * (b - a), (((1 - ξ) * Γa + ξ * Γb) / 2) * (b - a), T(δ))
    end
end
