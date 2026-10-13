"""
Post-processing (doc §9): surface velocity two ways, pressure coefficient, the leak residual and nodal averages.
"""
module Post

using LinearAlgebra
using StaticArrays
using ..Geometry: SurfaceGeometry
using ..DOFMap: DOFMapping, WakeTE
using ..Kernels: source_velocity, segment_velocity, FAR_DEFAULT
using ..Elements: FlatTriangle, IsoQuad, sheet_strength, gradient, gauss_legendre, area_normal, point
using ..Assembly: MixedMesh
using ..Wake: RigidWake, wake_elements

export element_values, sheet_strengths, surface_velocity, gradient_velocity, pressure_coefficient, nodal_average, body_velocity
export pressure_force, kutta_joukowski_force, trefftz_drag, spanwise_loading, surface_pressure

"""Nodal doublet values of every element's corners from the slot values"""
element_values(dm::DOFMapping, μ::AbstractVector) = [SVector(μ[cs[1]], μ[cs[2]], μ[cs[3]]) for cs in dm.corner_slot]

"""Vortex sheet γ = n × ∇μ of every element"""
sheet_strengths(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, μ::AbstractVector) =
    [sheet_strength(el, m) for (el, m) in zip(els, element_values(dm, μ))]

"""
    surface_velocity(els, dm, σ, μ, V∞; wake=nothing, μ_wake=nothing, kernel_far=FAR_DEFAULT, threaded=true)

Velocity on the fluid side at every element centroid from the singularities: V∞ + Σ σ_s w_s + Σ γ_s × w_s, the
element's own sheet taken on the side its normal points to (so it includes ½σn and ½γ × n), plus, with a `wake` and
the doublet strength of its nodes `μ_wake` (`wake_values`), the wake's sheets and tip filaments.
"""
function surface_velocity(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, σ::AbstractVector, μ::AbstractVector, V∞::SVector{3};
    wake::Union{Nothing,RigidWake}=nothing, μ_wake=nothing, kernel_far::Real=FAR_DEFAULT, threaded::Bool=true)
    γ = sheet_strengths(els, dm, μ)
    wels = isnothing(wake) ? eltype(els)[] : wake_elements(wake)
    wγ = isnothing(wake) ? typeof(γ[1])[] : [sheet_strength(el, SVector(μ_wake[t[1]], μ_wake[t[2]], μ_wake[t[3]])) for (el, t) in zip(wels, wake.tris)]
    fil = isnothing(wake) ? Tuple{typeof(els[1].n),typeof(els[1].n),eltype(μ),eltype(μ)}[] :
          [(wake.x[a], wake.x[b], μ_wake[a], μ_wake[b]) for (a, b) in wake.filaments]
    u = Vector{typeof(V∞ + els[1].n)}(undef, length(els))
    function at(t)
        x = els[t].centroid
        acc = V∞ + zero(els[t].n)
        for (s, el) in enumerate(els)
            w = source_velocity(el.src, x; side=s == t ? 1 : 0, far=kernel_far)
            acc += σ[s] * w + cross(γ[s], w)
        end
        for (el, g) in zip(wels, wγ)
            acc += cross(g, source_velocity(el.src, x; far=kernel_far))
        end
        for (a, b, μa, μb) in fil
            va, vb = segment_velocity(a, b, x)
            acc += μa * va + μb * vb
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

"""
    pressure_coefficient(u, V∞)
    pressure_coefficient(u, V∞, dμdt)

Pressure coefficient 1 - |u|²/|V∞|², and with the rate of change of the doublet (the outer perturbation potential) at
the same points, the unsteady Bernoulli equation's 1 - |u|²/|V∞|² - 2 ∂μ/∂t / |V∞|²
"""
pressure_coefficient(u::AbstractVector, V∞::SVector{3}) = [1 - dot(v, v) / dot(V∞, V∞) for v in u]
pressure_coefficient(u::AbstractVector, V∞::SVector{3}, dμdt::AbstractVector) =
    [1 - (dot(v, v) + 2 * d) / dot(V∞, V∞) for (v, d) in zip(u, dμdt)]

"""
    body_velocity(els, σ, γ, xs; kernel_far=FAR_DEFAULT, threaded=true)

Velocity the body's sources σ and vortex sheets γ induce at the points xs, which must lie off its triangles (a free
wake's nodes and particles): Σ σ_s w_s(x) + γ_s × w_s(x), without the freestream
"""
function body_velocity(els::AbstractVector{<:FlatTriangle}, σ::AbstractVector, γ::AbstractVector, xs::AbstractVector{<:SVector{3}};
    kernel_far::Real=FAR_DEFAULT, threaded::Bool=true)
    u = similar(xs)
    function at(i)
        x = xs[i]
        acc = zero(x)
        for (s, el) in enumerate(els)
            w = source_velocity(el.src, x; far=kernel_far)
            acc += σ[s] * w + cross(γ[s], w)
        end
        u[i] = acc
    end
    if threaded
        Threads.@threads for i in eachindex(xs)
            at(i)
        end
    else
        foreach(at, eachindex(xs))
    end
    return u
end

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

include("forces.jl")
include("mixed.jl")

end # module Post
