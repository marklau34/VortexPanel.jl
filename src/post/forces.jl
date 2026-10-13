# Forces and moments by the two independent routes of doc §9: pressure integration over the body, and Kutta–Joukowski
# with Trefftz-plane induced drag from the circulation the wake carries away (μ_w at its trailing edge nodes).
# All forces are in units of the dynamic pressure q∞ = ρ|V∞|²/2 (divide by the reference area for coefficients).

"""
    pressure_force(els, cp; origin=0)

Force and moment about `origin`, in units of q∞, of a pressure coefficient constant on each triangle (the gradient
route's, `gradient_velocity`): F = -Σ Cp n A and M = Σ (c - origin) × (-Cp n A), c the centroid; exact for such a Cp.
"""
function pressure_force(els::AbstractVector{<:FlatTriangle}, cp::AbstractVector; origin=zero(els[1].centroid))
    F, M = zero(els[1].n), zero(els[1].n)
    for (el, c) in zip(els, cp)
        f = -c * el.area * el.n
        F += f
        M += cross(el.centroid - origin, f)
    end
    return F, M
end

"""Segments (a, b) of the wake's trailing edge lines along τ, as indices of wake trailing edge nodes"""
te_segments(wte) = [(line[i], line[mod1(i + 1, length(line))]) for (line, closed) in zip(wte.lines, wte.closed)
                    for i in 1:(closed ? length(line) : length(line) - 1)]

"""
    kutta_joukowski_force(dm, μ_w, V∞)
    kutta_joukowski_force(wake_te, μ_w, V∞)

Force on the bound vorticity, in units of q∞, from the circulation the wake carries away: F = (2/|V∞|²) V∞ × ∫μ_w dl
along the trailing edge lines, oriented along τ. μ_w is linear between the wake's trailing edge nodes, so the
trapezoidal sum is exact; with μ_w = μ_u - μ_l (the potential jump across the wake) it is positive for a force along
V∞ × τ (lift, for a wing with τ along its span).
"""
function kutta_joukowski_force(wte::WakeTE, μ_w::AbstractVector, V∞::SVector{3})
    G = sum((μ_w[a] + μ_w[b]) / 2 * (wte.x[b] - wte.x[a]) for (a, b) in te_segments(wte); init=zero(V∞))
    return 2 / dot(V∞, V∞) * cross(V∞, G)
end

kutta_joukowski_force(dm::DOFMapping, μ_w::AbstractVector, V∞::SVector{3}) = kutta_joukowski_force(dm.wake_te, μ_w, V∞)

"""
    trefftz_drag(dm, μ_w, V∞)
    trefftz_drag(wake_te, μ_w, V∞)

Induced drag, in units of q∞, from the circulation shed along the trailing edge lines, in the Trefftz plane (normal to
V∞, onto which everything is projected). Each strip between neighbouring wake trailing edge nodes carries the mean of
their μ_w, Γ_k; trailing vortices ΔΓ = Γ_k - Γ_(k-1) leave the nodes (the end strips' ±Γ at the ends of an open line);
their normalwash at the middle of each strip on the lifting line (half its Trefftz plane value), w_k, gives
D = ρ Σ Γ_k w_k ℓ_k. The vector form of `trefftz_drag` in validation/common/loads.jl: an elliptic loading gives
CDi = CL²/(π AR) as the strips refine. A concentrated tip vortex (μ_w ≠ 0 at a line's end) would have infinite induced
drag in the exact Trefftz integral; the strips' midpoint downwash keeps it finite, so watch μ_w at the tips.
"""
function trefftz_drag(wte::WakeTE, μ_w::AbstractVector, V∞::SVector{3})
    d = normalize(V∞)
    P(x) = x - dot(x, d) * d
    vortices = Tuple{typeof(V∞),eltype(μ_w)}[]
    strips = Tuple{typeof(V∞),typeof(V∞),eltype(μ_w)}[]
    for (line, closed) in zip(wte.lines, wte.closed)
        m = length(line)
        nseg = closed ? m : m - 1
        Γ = [(μ_w[line[i]] + μ_w[line[mod1(i + 1, m)]]) / 2 for i in 1:nseg]
        for i in 1:m
            before = i > 1 ? Γ[i-1] : closed ? Γ[nseg] : zero(eltype(Γ))
            after = i <= nseg ? Γ[i] : zero(eltype(Γ))
            push!(vortices, (P(wte.x[line[i]]), after - before))
        end
        for i in 1:nseg
            push!(strips, (P(wte.x[line[i]]), P(wte.x[line[mod1(i + 1, m)]]), Γ[i]))
        end
    end
    D = zero(eltype(μ_w))
    for (a, b, Γk) in strips
        c, ℓ = (a + b) / 2, norm(b - a)
        n = cross(d, (b - a) / ℓ)
        v = sum(ΔΓ / 4π * cross(d, c - p) / sum(abs2, c - p) for (p, ΔΓ) in vortices)
        D += Γk * dot(n, v) * ℓ
    end
    return 2D / dot(V∞, V∞)
end

trefftz_drag(dm::DOFMapping, μ_w::AbstractVector, V∞::SVector{3}) = trefftz_drag(dm.wake_te, μ_w, V∞)

"""
    spanwise_loading(els, cp, dm, μ_w, V∞; line=1)

Loading along trailing edge line `line` (open), per strip between neighbouring wake trailing edge nodes, with the
coordinate s along the line (from its first node, along its end-to-end direction t): the strip's ends `s0`, `s1`, its
circulation Γ (the mean of μ_w), the force per unit span normal to V∞ and t (lift, for a wing) from Kutta–Joukowski,
2Γ/|V∞|, and from pressure, the force of the triangles whose centroid falls in the strip (`lift_p`, and `drag_p` along
V∞), all in units of q∞. Divide by the local chord for sectional coefficients.
"""
function spanwise_loading(els::AbstractVector{<:FlatTriangle}, cp::AbstractVector, dm::DOFMapping, μ_w::AbstractVector, V∞::SVector{3}; line::Integer=1)
    wte = dm.wake_te
    nodes = wte.lines[line]
    x0 = wte.x[nodes[1]]
    t = normalize(wte.x[nodes[end]] - x0)
    d = normalize(V∞)
    eL = normalize(cross(d, t))
    s = [dot(wte.x[k] - x0, t) for k in nodes]
    n = length(nodes) - 1
    lift, drag = zeros(eltype(cp), n), zeros(eltype(cp), n)
    for (el, c) in zip(els, cp)
        k = searchsortedlast(s, dot(el.centroid - x0, t))
        1 <= k <= n || continue
        f = -c * el.area * el.n
        lift[k] += dot(f, eL)
        drag[k] += dot(f, d)
    end
    Γ = [(μ_w[nodes[k]] + μ_w[nodes[k+1]]) / 2 for k in 1:n]
    ds = diff(s)
    return (s0=s[1:n], s1=s[2:end], Γ=Γ, lift_kj=2Γ / norm(V∞), lift_p=lift ./ ds, drag_p=drag ./ ds)
end
