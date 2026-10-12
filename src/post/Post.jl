"""
Post-processing (doc §9): surface velocity two ways, pressure coefficient, the leak residual and nodal averages.
"""
module Post

using LinearAlgebra
using StaticArrays
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping
using ..Kernels: source_velocity, FAR_DEFAULT
using ..Elements: FlatTriangle, sheet_strength, gradient

export element_values, sheet_strengths, surface_velocity, gradient_velocity, pressure_coefficient, nodal_average

"""Nodal doublet values of every element's corners from the slot values"""
element_values(dm::DOFMapping, μ::AbstractVector) = [SVector(μ[cs[1]], μ[cs[2]], μ[cs[3]]) for cs in dm.corner_slot]

"""Vortex sheet γ = n × ∇μ of every element"""
sheet_strengths(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, μ::AbstractVector) =
    [sheet_strength(el, m) for (el, m) in zip(els, element_values(dm, μ))]

"""
    surface_velocity(els, dm, σ, μ, V∞; kernel_far=FAR_DEFAULT, threaded=true)

Velocity on the fluid side at every element centroid from the singularities: V∞ + Σ σ_s w_s + Σ γ_s × w_s, the
element's own sheet taken on the side its normal points to (so it includes ½σn and ½γ × n).
"""
function surface_velocity(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, σ::AbstractVector, μ::AbstractVector, V∞::SVector{3};
    kernel_far::Real=FAR_DEFAULT, threaded::Bool=true)
    γ = sheet_strengths(els, dm, μ)
    u = Vector{typeof(V∞ + els[1].n)}(undef, length(els))
    function at(t)
        x = els[t].centroid
        acc = V∞ + zero(els[t].n)
        for (s, el) in enumerate(els)
            w = source_velocity(el.src, x; side=s == t ? 1 : 0, far=kernel_far)
            acc += σ[s] * w + cross(γ[s], w)
        end
        u[t] = acc
    end
    if threaded
        Threads.@threads for t in eachindex(els)
            at(t)
        end
    else
        foreach(at, eachindex(els))
    end
    return u
end

"""
    gradient_velocity(els, dm, μ, V∞)

Tangential velocity on the fluid side from the element gradient of the doublet, V∞ - (V∞·n)n + ∇μ: μ is the outer
perturbation potential (the interior perturbation is uniform), so this needs no influence sums and no stencils.
"""
gradient_velocity(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, μ::AbstractVector, V∞::SVector{3}) =
    [V∞ - dot(V∞, el.n) * el.n + gradient(el, m) for (el, m) in zip(els, element_values(dm, μ))]

"""Pressure coefficient 1 - |u|²/|V∞|²"""
pressure_coefficient(u::AbstractVector, V∞::SVector{3}) = [1 - dot(v, v) / dot(V∞, V∞) for v in u]

"""Area-weighted average of element values at the nodes, for plots"""
function nodal_average(g::SurfaceGeometry, els::AbstractVector{<:FlatTriangle}, values::AbstractVector)
    acc = zeros(eltype(values), length(g.x))
    wsum = zeros(length(g.x))
    for (t, tri) in enumerate(g.tris), v in tri
        acc[v] += els[t].area * values[t]
        wsum[v] += els[t].area
    end
    return acc ./ max.(wsum, eps())
end

end # module Post
