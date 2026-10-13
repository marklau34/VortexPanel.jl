# Post-processing of mixed meshes (isoparametric quads and flat triangles, `MixedMesh`): the pressure of the gradient
# route at Gauss points, and the force routes and spanwise loading built on them.

"""
    surface_pressure(mesh::MixedMesh, μ, V∞; m=3)

Points of a mixed mesh with the pressure coefficient of the gradient route (the tangential velocity
V∞ - (V∞·n)n + ∇μ, as `gradient_velocity`) and their outward area vectors: the m×m Gauss points of each quad, whose
doublet gradient varies across it, then the centroid of each triangle. `element` gives each point's element (quads
first). Returns (x, dA, cp, element).

The m×m points integrate the forces. To compare Cp with a smooth surface pointwise, `m=1` gives each quad's
parameter centre: on a curved surface a nearly flat quad's normal there matches the surface's to second order, while
off the centre it is off by O(h), as a flat triangle's is off its centroid (validation stage 7).
"""
function surface_pressure(mesh::MixedMesh{T}, μ::AbstractVector, V∞::SVector{3}; m::Integer=3) where {T}
    ξs, ws = gauss_legendre(m)
    x, dA, cp, element = SVector{3,T}[], SVector{3,T}[], T[], Int[]
    V2 = dot(V∞, V∞)
    for (e, (Q, sl)) in enumerate(zip(mesh.quads, mesh.quad_slots))
        μq = SVector(μ[sl[1]], μ[sl[2]], μ[sl[3]], μ[sl[4]])
        for i in 1:m, j in 1:m
            p = SVector(ξs[i], ξs[j])
            nJ = area_normal(Q, p)
            n = nJ / norm(nJ)
            u = V∞ - dot(V∞, n) * n + gradient(Q, μq, p)
            push!(x, point(Q, p))
            push!(dA, ws[i] * ws[j] * nJ)
            push!(cp, 1 - dot(u, u) / V2)
            push!(element, e)
        end
    end
    nq = length(mesh.quads)
    for (t, (el, sl)) in enumerate(zip(mesh.tris, mesh.tri_slots))
        u = V∞ - dot(V∞, el.n) * el.n + gradient(el, SVector(μ[sl[1]], μ[sl[2]], μ[sl[3]]))
        push!(x, el.centroid)
        push!(dA, el.area * el.n)
        push!(cp, 1 - dot(u, u) / V2)
        push!(element, nq + t)
    end
    return (x=x, dA=dA, cp=cp, element=element)
end

"""
    pressure_force(pts; origin=0)

Force and moment about `origin`, in units of q∞, of the pressure at the points of `surface_pressure`: F = -Σ Cp dA,
M = Σ (x - origin) × (-Cp dA)
"""
function pressure_force(pts::NamedTuple; origin=zero(pts.x[1]))
    F, M = zero(pts.x[1]), zero(pts.x[1])
    for (x, a, c) in zip(pts.x, pts.dA, pts.cp)
        F -= c * a
        M += cross(x - origin, -c * a)
    end
    return F, M
end

"""
    spanwise_loading(pts, dm, μ_w, V∞; line=1)

`spanwise_loading` for the points of `surface_pressure` (a mixed mesh): each point's pressure force goes to the strip
its position falls in
"""
function spanwise_loading(pts::NamedTuple, dm::DOFMapping, μ_w::AbstractVector, V∞::SVector{3}; line::Integer=1)
    wte = dm.wake_te
    nodes = wte.lines[line]
    x0 = wte.x[nodes[1]]
    t = normalize(wte.x[nodes[end]] - x0)
    d = normalize(V∞)
    eL = normalize(cross(d, t))
    s = [dot(wte.x[k] - x0, t) for k in nodes]
    n = length(nodes) - 1
    lift, drag = zeros(eltype(μ_w), n), zeros(eltype(μ_w), n)
    for (x, a, c) in zip(pts.x, pts.dA, pts.cp)
        k = searchsortedlast(s, dot(x - x0, t))
        1 <= k <= n || continue
        lift[k] -= c * dot(a, eL)
        drag[k] -= c * dot(a, d)
    end
    Γ = [(μ_w[nodes[k]] + μ_w[nodes[k+1]]) / 2 for k in 1:n]
    ds = diff(s)
    return (s0=s[1:n], s1=s[2:end], Γ=Γ, lift_kj=2Γ / norm(V∞), lift_p=lift ./ ds, drag_p=drag ./ ds)
end
