# Quadrature rules on the parameter square [0, 1]² of an isoparametric quad (the test side of the Galerkin integrals),
# built from the graded one-dimensional and triangle rules; weights are fractions of the parameter area.

"""A quadrature rule on the parameter square: points (ξ, η) and weights that sum to one"""
struct SquareRule
    p::Vector{SVector{2,Float64}}
    w::Vector{Float64}
end

Base.length(r::SquareRule) = length(r.w)

"""The m×m tensor Gauss–Legendre rule"""
function square_gauss(m::Integer)
    x, w = gauss_legendre(m)
    return SquareRule([SVector(a, b) for a in x for b in x], [wa * wb for wa in w for wb in w])
end

"""Rule r turned by k - 1 quarter turns about the centre, (ξ, η) → (1 - η, ξ), which take side j and corner j to j + k - 1"""
function turn(r::SquareRule, k::Integer)
    p = r.p
    for _ in 2:k
        p = [SVector(1 - q[2], q[1]) for q in p]
    end
    return SquareRule(p, copy(r.w))
end

"""Append the triangle rule r mapped onto the parameter triangle with corners a, b, c (area fraction frac)"""
function append_triangle!(ps, ws, r::BaryRule, a::SVector{2}, b::SVector{2}, c::SVector{2}, frac::Real)
    for (λ, w) in zip(r.λ, r.w)
        push!(ps, λ[1] * a + λ[2] * b + λ[3] * c)
        push!(ws, w * frac)
    end
    return nothing
end

"""
    square_edge_rule(levels, ratio, order)

Rule for integrands log-singular along side 1 (η = 0, corner 1 to corner 2): graded toward η = 0 with `levels` levels
and toward both ends of the side with half as many (`graded_gauss`, `graded_gauss_both`).
"""
function square_edge_rule(levels::Integer, ratio::Real, order::Integer)
    ηs, wη = graded_gauss(levels, ratio, order)
    ξs, wξ = graded_gauss_both(max(1, levels ÷ 2), ratio, order)
    return SquareRule([SVector(ξ, η) for η in ηs for ξ in ξs], [a * b for a in wη for b in wξ])
end

"""Rule for integrands log-singular at corner 1 (0, 0): the two triangles of the diagonal from it, each with the vertex rule"""
function square_vertex_rule(levels::Integer, ratio::Real, order::Integer)
    base = vertex_graded_rule(levels, ratio, order)
    ps, ws = SVector{2,Float64}[], Float64[]
    append_triangle!(ps, ws, base, SQUARE[1], SQUARE[2], SQUARE[3], 0.5)
    append_triangle!(ps, ws, base, SQUARE[1], SQUARE[3], SQUARE[4], 0.5)
    return SquareRule(ps, ws)
end

"""Rule for integrands log-singular along all four sides (the element's own field): four triangles from the centre, each graded toward its side"""
function square_self_rule(levels::Integer, ratio::Real, order::Integer)
    base = edge_graded_rule(levels, ratio, order)
    ps, ws = SVector{2,Float64}[], Float64[]
    centre = SVector(0.5, 0.5)
    for k in 1:4
        append_triangle!(ps, ws, base, SQUARE[k], SQUARE[mod1(k + 1, 4)], centre, 0.25)
    end
    return SquareRule(ps, ws)
end

"""
Rules of the quadrature tiers on a test quad's parameter square: `far` and `near` (tensor Gauss), `self`, the rules
graded toward each side (`edge[k]`, side k from corner k to k + 1) and each corner (`vertex[k]`), and the triangle
vertex rule they are built from (`vertex_base`), for corners a shear makes wide (the assembly's `vertex_rect!`)
"""
struct SquareTiers
    far::SquareRule
    near::SquareRule
    self::SquareRule
    edge::NTuple{4,SquareRule}
    vertex::NTuple{4,SquareRule}
    vertex_base::BaryRule
end

function SquareTiers(; far_points::Integer=2, near_points::Integer=4, levels::Integer=6, ratio::Real=0.15, order::Integer=3)
    e = square_edge_rule(levels, ratio, order)
    v = square_vertex_rule(levels, ratio, order)
    return SquareTiers(square_gauss(far_points), square_gauss(near_points), square_self_rule(levels, ratio, order),
        ntuple(k -> turn(e, k), 4), ntuple(k -> turn(v, k), 4), vertex_graded_rule(levels, ratio, order))
end

"""
    square_near_rule(Q, S, base, ratio, depth)

Rule on test quad Q for a near source S (anything with a `distance_to`): the parameter square is split in four while a
piece's centre is closer to S than `ratio` times the piece's size, up to `depth` levels; the leaves carry `base`.
"""
function square_near_rule(Q::IsoQuad, S, base::SquareRule, ratio::Real, depth::Integer)
    ps, ws = SVector{2,Float64}[], Float64[]
    near_squares!(ps, ws, Q, S, base, ratio, depth, SVector(0.0, 0.0), 1.0, 0)
    return SquareRule(ps, ws)
end

"""Append `base` on the piece of side `size` at `lo`, split while close to S"""
function near_squares!(ps, ws, Q::IsoQuad, S, base::SquareRule, ratio::Real, depth::Integer, lo::SVector{2}, size::Real, level::Integer)
    corners = (point(Q, lo), point(Q, lo + SVector(size, 0.0)), point(Q, lo + SVector(size, size)), point(Q, lo + SVector(0.0, size)))
    extent = max(norm(corners[3] - corners[1]), norm(corners[4] - corners[2]))
    if level < depth && distance_to(S, point(Q, lo .+ size / 2)) < ratio * extent
        h = size / 2
        for o in (SVector(0.0, 0.0), SVector(h, 0.0), SVector(0.0, h), SVector(h, h))
            near_squares!(ps, ws, Q, S, base, ratio, depth, lo + o, h, level + 1)
        end
    else
        for (q, w) in zip(base.p, base.w)
            push!(ps, lo + size * q)
            push!(ws, w * size^2)
        end
    end
    return nothing
end
