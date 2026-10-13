# Two-dimensional reference: inviscid flow about an airfoil section by a linear-strength vortex panel method (Katz &
# Plotkin, Low-Speed Aerodynamics, §11.4), converged by refinement, to judge the sectional lift and moment of the 3-D
# solvers in the limit of a long wing. Requires LinearAlgebra and StaticArrays.

"""
    naca00_section(n; t=0.15, chord=1)

Nodes of a NACA 00xx section with a sharp trailing edge (the closed-TE coefficient -0.1036, as PAGE), clockwise from
the trailing edge along the lower surface to the leading edge and back along the upper surface: 2n panels with cosine
spacing, node 1 and node 2n+1 both at the trailing edge.
"""
function naca00_section(n::Integer; t::Real=0.15, chord::Real=1.0)
    yt(ξ) = 5t * (0.2969sqrt(ξ) - 0.1260ξ - 0.3516ξ^2 + 0.2843ξ^3 - 0.1036ξ^4)
    ξ = [(1 - cos(π * k / n)) / 2 for k in 0:n] # 0 at the leading edge
    lower = [SVector(chord * ξ[k], -chord * yt(ξ[k])) for k in n+1:-1:1]
    upper = [SVector(chord * ξ[k], chord * yt(ξ[k])) for k in 2:n+1]
    return [lower; upper]
end

"""
    joukowski_section(n; a=1, m=0.1)

Nodes (clockwise from the trailing edge, as `naca00_section`) of the symmetric Joukowski airfoil from the circle of
radius a + m centred at (-m, 0) under z = ζ + a²/ζ, and its exact lift slope: cl = 8π(a + m) sin α / chord.
"""
function joukowski_section(n::Integer; a::Real=1.0, m::Real=0.1)
    R = a + m
    θ = [π * k / n for k in 0:2n] # θ = 0 at the trailing edge ζ = a; clockwise: decreasing angle on the circle
    z(t) = (ζ = complex(-m, 0) + R * cis(-t); w = ζ + a^2 / ζ; SVector(real(w), imag(w)))
    pts = z.(θ)
    pts[end] = pts[1]
    chord = 2a - (-(2m + a) - a^2 / (2m + a))
    return pts, 8π * R / chord
end

"""
    linear_vortex_panel(p1, p2, x)

Velocity at x induced by a vortex panel from p1 to p2 whose strength (clockwise positive) varies linearly from 1 at p1
to 0 at p2, and the one varying from 0 to 1: (v1, v2), from the closed-form integrals in the panel's frame (x along the
panel, z to its left).
"""
function linear_vortex_panel(p1::SVector{2}, p2::SVector{2}, x::SVector{2})
    ℓ = norm(p2 - p1)
    t = (p2 - p1) / ℓ
    n = SVector(-t[2], t[1])
    xl, zl = dot(x - p1, t), dot(x - p1, n)
    θ1, θ2 = atan(zl, xl), atan(zl, xl - ℓ)
    Δθ = θ2 - θ1
    lr = log(hypot(xl, zl) / hypot(xl - ℓ, zl))
    u1 = (Δθ * (1 - xl / ℓ) + zl / ℓ * lr) / 2π
    u2 = (Δθ * xl / ℓ - zl / ℓ * lr) / 2π
    w1 = -((1 - xl / ℓ) * lr + 1 - zl / ℓ * Δθ) / 2π
    w2 = -(xl / ℓ * lr - 1 + zl / ℓ * Δθ) / 2π
    return u1 * t + w1 * n, u2 * t + w2 * n
end

"""
    panel2d(nodes, α; x_ref=0.25)

Linear-strength vortex panel solution about the closed airfoil `nodes` (clockwise, the trailing edge first and last)
at angle of attack α (degrees), unit freestream: flow tangency at the panel midpoints and the Kutta condition
γ₁ + γ_end = 0. Returns the lift coefficient from the circulation (`cl_gamma`, 2Γ/(U c)) and from the pressure
(`cl`), the drag coefficient from the pressure (`cd`, zero up to the discretization), the moment coefficient about
x_ref·chord (`cm`, nose up positive), and the midpoints and Cp of the panels (the outer tangential velocity is γ).
"""
function panel2d(nodes::AbstractVector{<:SVector{2}}, α::Real; x_ref::Real=0.25)
    np = length(nodes) - 1
    V = SVector(cosd(α), sind(α))
    mids = [(nodes[j] + nodes[j+1]) / 2 for j in 1:np]
    normals = [(d = nodes[j+1] - nodes[j]; SVector(-d[2], d[1]) / norm(d)) for j in 1:np]
    A = zeros(np + 1, np + 1)
    b = zeros(np + 1)
    for i in 1:np
        for j in 1:np
            v1, v2 = linear_vortex_panel(nodes[j], nodes[j+1], mids[i])
            A[i, j] += dot(v1, normals[i])
            A[i, j+1] += dot(v2, normals[i])
        end
        b[i] = -dot(V, normals[i])
    end
    A[np+1, 1] = A[np+1, np+1] = 1 # Kutta
    γ = A \ b
    ℓ = [norm(nodes[j+1] - nodes[j]) for j in 1:np]
    chord = maximum(p -> p[1], nodes) - minimum(p -> p[1], nodes)
    Γ = sum((γ[j] + γ[j+1]) / 2 * ℓ[j] for j in 1:np)
    cp = [1 - ((γ[j] + γ[j+1]) / 2)^2 for j in 1:np]
    F = -sum(cp[j] * normals[j] * ℓ[j] for j in 1:np) / chord
    eL = SVector(-sind(α), cosd(α))
    xr = SVector(minimum(p -> p[1], nodes) + x_ref * chord, 0.0)
    M = sum((r = mids[j] - xr; f = -cp[j] * normals[j] * ℓ[j]; r[2] * f[1] - r[1] * f[2]) for j in 1:np) / chord^2
    return (cl_gamma=2Γ / chord, cl=dot(F, eL), cd=dot(F, V), cm=M, x=mids, cp=cp, γ=γ)
end
