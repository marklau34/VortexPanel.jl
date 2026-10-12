# Quadrature rules on triangles, in barycentric coordinates with weights as fractions of the triangle's area (sum 1)

"""A quadrature rule on a triangle: barycentric points `λ` and weights `w` that sum to one (fractions of the area)"""
struct BaryRule
    λ::Vector{SVector{3,Float64}}
    w::Vector{Float64}
end

Base.length(r::BaryRule) = length(r.w)

"""Gauss–Legendre nodes and weights on [0, 1]"""
function gauss_legendre(n::Integer)
    x = zeros(n)
    w = zeros(n)
    for i in 1:cld(n, 2)
        z = cos(π * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            dp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / dp
            z -= dz
            abs(dz) < 1e-16 && break
        end
        p1, p2 = 1.0, 0.0
        for j in 1:n
            p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
        end
        dp = n * (z * p1 - p2) / (z^2 - 1)
        x[i], x[n+1-i] = (1 - z) / 2, (1 + z) / 2
        w[i] = w[n+1-i] = 1 / ((1 - z^2) * dp^2)
    end
    return x, w
end

"""
    graded_gauss(levels, ratio, order)

Rule on [0, 1] for integrands with a (log) singularity at 0 (hp-graded): the intervals [0, ratio^levels],
[ratio^levels, ratio^(levels-1)], ..., [ratio, 1] carry Gauss–Legendre rules whose order grows by one per interval away
from the singularity, from `order` on the innermost one. Every outer interval is a scaled copy of [ratio, 1], so a
fixed order would leave the same relative error on each; the growing order makes the error fall exponentially with
the levels. On the innermost interval [0, a] the substitution t = a u² turns a log singularity into u ln u.
"""
function graded_gauss(levels::Integer, ratio::Real, order::Integer)
    breaks = [0.0; [ratio^k for k in levels:-1:1]; 1.0]
    x, w = Float64[], Float64[]
    for k in 1:length(breaks)-1
        xg, wg = gauss_legendre(order + k - 1)
        a, b = breaks[k], breaks[k+1]
        if k == 1
            append!(x, b .* xg .^ 2)
            append!(w, 2b .* xg .* wg)
        else
            append!(x, a .+ (b - a) .* xg)
            append!(w, (b - a) .* wg)
        end
    end
    return x, w
end

"""Rule on [0, 1] graded toward both ends (singularities at 0 and 1)"""
function graded_gauss_both(levels::Integer, ratio::Real, order::Integer)
    x, w = graded_gauss(levels, ratio, order)
    return [x ./ 2; 1 .- reverse(x) ./ 2], [w ./ 2; reverse(w) ./ 2]
end

"""
    triangle_rule(n)

Symmetric rules with `n` points: 1 (degree 1), 3 (degree 2), 6 (degree 4, Strang–Fix) and 7 (degree 5, Radon).
"""
function triangle_rule(n::Integer)
    perm(a, b) = [SVector(a, b, b), SVector(b, a, b), SVector(b, b, a)]
    if n == 1
        return BaryRule([SVector(1 / 3, 1 / 3, 1 / 3)], [1.0])
    elseif n == 3
        return BaryRule(perm(2 / 3, 1 / 6), fill(1 / 3, 3))
    elseif n == 6
        a, b = 0.445948490915965, 0.091576213509771
        return BaryRule([perm(1 - 2a, a); perm(1 - 2b, b)], [fill(0.223381589678011, 3); fill(0.109951743655322, 3)])
    elseif n == 7
        s = sqrt(15)
        a, b = (6 - s) / 21, (6 + s) / 21
        return BaryRule([[SVector(1 / 3, 1 / 3, 1 / 3)]; perm(1 - 2a, a); perm(1 - 2b, b)],
            [9 / 40; fill((155 - s) / 1200, 3); fill((155 + s) / 1200, 3)])
    end
    throw(ArgumentError("no $n-point triangle rule (1, 3, 6 or 7)"))
end

"""Rule `r` with its barycentric coordinates rotated so that base corner 1 becomes corner k (2 → k+1, 3 → k+2)"""
function rotate_rule(r::BaryRule, k::Integer)
    k == 1 && return r
    idx = (mod1(2 - k, 3), mod1(3 - k, 3), mod1(4 - k, 3)) # λ_T[j] = λ_base[idx[j]]
    return BaryRule([SVector(l[idx[1]], l[idx[2]], l[idx[3]]) for l in r.λ], copy(r.w))
end

"""
    edge_graded_rule(levels, ratio, order)

Rule for integrands log-singular along the edge between corners 1 and 2 (and at its ends): collapsed coordinates
x = (1 - t)[(1 - s)x₁ + s x₂] + t x₃ (dA/A = 2(1 - t) ds dt), graded toward t = 0 (the edge) with `levels` levels and
toward both ends in s with half as many (the ends carry only point singularities).
"""
function edge_graded_rule(levels::Integer, ratio::Real, order::Integer)
    ts, wt = graded_gauss(levels, ratio, order)
    ss, ws = graded_gauss_both(max(1, levels ÷ 2), ratio, order)
    λ, w = SVector{3,Float64}[], Float64[]
    for (t, a) in zip(ts, wt), (s, b) in zip(ss, ws)
        push!(λ, SVector((1 - t) * (1 - s), (1 - t) * s, t))
        push!(w, 2 * (1 - t) * a * b)
    end
    return BaryRule(λ, w)
end

"""
    vertex_graded_rule(levels, ratio, order)

Rule for integrands log-singular at corner 1: x = x₁ + ρ[(x₂ - x₁) + s(x₃ - x₂)] (dA/A = 2ρ dρ ds), graded toward
ρ = 0.
"""
function vertex_graded_rule(levels::Integer, ratio::Real, order::Integer)
    ρs, wρ = graded_gauss(levels, ratio, order)
    ss, ws = gauss_legendre(order + levels)
    λ, w = SVector{3,Float64}[], Float64[]
    for (ρ, a) in zip(ρs, wρ), (s, b) in zip(ss, ws)
        push!(λ, SVector(1 - ρ, ρ * (1 - s), ρ * s))
        push!(w, 2 * ρ * a * b)
    end
    return BaryRule(λ, w)
end

"""
    self_rule(levels, ratio, order)

Rule for integrands log-singular along all three edges (the triangle's own field): the triangle is split at its
centroid into three triangles, each with an `edge_graded_rule` toward its outer edge.
"""
function self_rule(levels::Integer, ratio::Real, order::Integer)
    base = edge_graded_rule(levels, ratio, order)
    c = SVector(1 / 3, 1 / 3, 1 / 3)
    λ, w = SVector{3,Float64}[], Float64[]
    for k in 1:3
        ek = SVector(ntuple(j -> Float64(j == k), 3))
        ek1 = SVector(ntuple(j -> Float64(j == mod1(k + 1, 3)), 3))
        for (μ, wb) in zip(base.λ, base.w)
            push!(λ, μ[1] * ek + μ[2] * ek1 + μ[3] * c)
            push!(w, wb / 3)
        end
    end
    return BaryRule(λ, w)
end

"""
Rules of the quadrature tiers on a test triangle: `far` and `near` (symmetric), `self`, and the rules graded toward
each edge (`edge[k]`, the edge from corner k to k+1) and each corner (`vertex[k]`).
"""
struct TierRules
    far::BaryRule
    near::BaryRule
    self::BaryRule
    edge::NTuple{3,BaryRule}
    vertex::NTuple{3,BaryRule}
end

function TierRules(; far_points::Integer=3, near_points::Integer=7, levels::Integer=6, ratio::Real=0.15, order::Integer=3)
    e = edge_graded_rule(levels, ratio, order)
    v = vertex_graded_rule(levels, ratio, order)
    return TierRules(triangle_rule(far_points), triangle_rule(near_points), self_rule(levels, ratio, order),
        (e, rotate_rule(e, 2), rotate_rule(e, 3)), (v, rotate_rule(v, 2), rotate_rule(v, 3)))
end
