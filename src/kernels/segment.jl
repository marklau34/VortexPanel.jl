"""
    segment_velocity(a, b, x)

Velocity at x of a straight vortex filament from a to b whose circulation falls linearly from 1 at a to 0 at b, and of
the one rising from 0 at a to 1 at b: the pair (va, vb), so a filament with end circulations Γa and Γb (positive by
the right-hand rule about a → b) induces Γa va + Γb vb. From Biot–Savart, u = (1/4π)∫ Γ(s) dl × (x - y)/|x - y|³:

    (va, vb) = (ra × rb) / (4π (|ra||rb| + ra·rb)) · (1/|ra|, 1/|rb|),   ra = x - a, rb = x - b,

whose sum is the constant-strength segment (Katz & Plotkin). The denominator is formed without cancellation next to
the segment, |ra||rb| + ra·rb = |ra × rb|²/(|ra||rb| - ra·rb) when ra·rb < 0. A filament of varying strength is not
divergence free on its own: it is completed by the vortex sheets its varying strength feeds (an edge of a linear
doublet panel). Singular on the segment (zero is returned there and at its ends); zero on its extension.
"""
@inline function segment_velocity(a::SVector{3}, b::SVector{3}, x::SVector{3})
    ra, rb = x - a, x - b
    Ra, Rb = norm(ra), norm(rb)
    c = cross(ra, rb)
    p = dot(ra, rb)
    den = p >= 0 ? Ra * Rb + p : sum(abs2, c) / (Ra * Rb - p)
    z = zero(c)
    (iszero(den) || iszero(Ra) || iszero(Rb)) && return z, z
    f = c / (4 * oftype(den, π) * den)
    return f / Ra, f / Rb
end
