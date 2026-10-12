"""Interior angle of the closed polyline `nodes` at position i"""
function loop_angle(x::AbstractVector{SVector{3,T}}, nodes::AbstractVector{Int}, i::Integer) where {T}
    m = length(nodes)
    a = x[nodes[mod1(i - 1, m)]] - x[nodes[i]]
    b = x[nodes[mod1(i + 1, m)]] - x[nodes[i]]
    return acos(clamp(dot(a, b) / (norm(a) * norm(b)), -one(T), one(T)))
end

"""
    loop_corners(x, nodes, te_candidates)

Trailing edge corner (a loop node on a trailing edge, the sharpest if several; otherwise the sharpest loop node),
leading edge corner (the loop node farthest from it), and the two chains from the leading to the trailing edge corner.
"""
function loop_corners(x::AbstractVector{SVector{3,T}}, nodes::AbstractVector{Int}, te_candidates::Set{Int}) where {T}
    m = length(nodes)
    angles = [loop_angle(x, nodes, i) for i in 1:m]
    cand = [i for i in 1:m if nodes[i] in te_candidates]
    ite = isempty(cand) ? argmin(angles) : cand[argmin(angles[cand])]
    ile = argmax([norm(x[nodes[i]] - x[nodes[ite]]) for i in 1:m])
    fwd = [nodes[mod1(ile + k, m)] for k in 0:mod(ite - ile, m)]
    bwd = [nodes[mod1(ile - k, m)] for k in 0:mod(ile - ite, m)]
    return nodes[ile], nodes[ite], fwd, bwd
end

"""
    ladder(x, A, B, tie)

Triangles between two chains that share their first (leading edge) and last (trailing edge) nodes, using only their
nodes. From the leading edge the ladder advances along the chain whose new diagonal is shorter; ties go to the chain
behind in normalised arclength, then alternate.
"""
function ladder(x::AbstractVector{SVector{3,T}}, A::AbstractVector{Int}, B::AbstractVector{Int}, tie::T) where {T}
    p, q = length(A) - 1, length(B) - 1
    tris = SVector{3,Int}[]
    if p == 1 && q == 1
        return tris
    elseif p == 1
        for j in 1:q-1
            push!(tris, SVector(B[j], B[j+1], B[end]))
        end
        return tris
    elseif q == 1
        for i in 1:p-1
            push!(tris, SVector(A[i], A[i+1], A[end]))
        end
        return tris
    end
    sA = arclength(x, A)
    sB = arclength(x, B)
    sA ./= sA[end]
    sB ./= sB[end]
    push!(tris, SVector(A[1], A[2], B[2]))
    i, j = 2, 2
    toggle = false
    while !(i == p && j == q)
        advance_a = if i == p
            false
        elseif j == q
            true
        else
            dA = norm(x[A[i+1]] - x[B[j]])
            dB = norm(x[A[i]] - x[B[j+1]])
            if abs(dA - dB) > tie * max(dA, dB)
                dA < dB
            elseif abs(sA[i+1] - sB[j+1]) > tie
                sA[i+1] < sB[j+1]
            else
                toggle = !toggle
            end
        end
        if advance_a
            push!(tris, SVector(A[i], A[i+1], B[j]))
            i += 1
        else
            push!(tris, SVector(A[i], B[j+1], B[j]))
            j += 1
        end
    end
    push!(tris, SVector(A[p], A[p+1], B[q]))
    return tris
end

"""Orient cap triangles opposite to the loop's boundary half-edges (loop nodes in half-edge order)"""
function orient_cap(tris::AbstractVector{SVector{3,Int}}, loop::AbstractVector{Int})
    m = length(loop)
    along = Set{Tuple{Int,Int}}((loop[k], loop[mod1(k + 1, m)]) for k in 1:m)
    return map(tris) do t
        directed = ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
        any(in(along), directed) ? SVector(t[1], t[3], t[2]) : t
    end
end

"""Largest distance of the points from their best fit plane"""
function plane_deviation(pts::AbstractVector{SVector{3,T}}) where {T}
    c = sum(pts) / length(pts)
    C = sum((p - c) * (p - c)' for p in pts)
    n = eigen(Symmetric(Matrix(C))).vectors[:, 1]
    return maximum(abs(dot(p - c, n)) for p in pts)
end

"""Projection records of `chain` onto `other` (segment, parameter) and the largest distance"""
function chain_projections(x, chain::AbstractVector{Int}, other::AbstractVector{Int})
    locs = [project_on_chain(x, x[v], other) for v in chain]
    return first.(locs), maximum(last.(locs))
end

"""
Projection records of `chain` onto `other` (segment, parameter) at the same normalised arclength, for two chains with
common end points (the two sides of an opening), and the largest distance between paired points
"""
function arclength_projections(x::AbstractVector{SVector{3,T}}, chain::AbstractVector{Int}, other::AbstractVector{Int}) where {T}
    s = arclength(x, chain)
    so = arclength(x, other)
    s ./= s[end]
    so ./= so[end]
    width = zero(T)
    locs = Vector{Tuple{Int,T}}(undef, length(chain))
    for (i, v) in enumerate(chain)
        k = clamp(searchsortedlast(so, s[i]), 1, length(other) - 1)
        t = clamp((s[i] - so[k]) / (so[k+1] - so[k]), zero(T), one(T))
        width = max(width, norm(x[v] - ((1 - t) * x[other[k]] + t * x[other[k+1]])))
        locs[i] = (k, t)
    end
    return locs, width
end

"""
    close_opening(x, loop, te_candidates, mode, opts)

Corners and chains of an opening and, for `mode = :cap`, its cap triangles; for `:zip`, the chain pair that ties its
two sides. Returns the `Opening` (with an empty triangle range) and the triangles or chain pair.
"""
function close_opening(x::AbstractVector{SVector{3,T}}, loop::BoundaryLoop, te_candidates::Set{Int}, mode::Symbol, opts::GeometryOptions) where {T}
    nodes = loop.nodes
    pts = x[nodes]
    diam = maximum(norm(p - q) for p in pts, q in pts)
    planarity = length(nodes) >= 3 ? plane_deviation(pts) / diam : zero(T)
    le, te, A, B = loop_corners(x, nodes, te_candidates)
    up = SVector{3,T}(opts.up_ref)
    meanup(c) = sum(v -> dot(x[v], up), c) / length(c)
    upper, lower = meanup(A) >= meanup(B) ? (A, B) : (B, A)
    tris = SVector{3,Int}[]
    pair = nothing
    if mode == :cap && loop.closed
        # the ladder starts along the upper chain on ties, so mirror image openings get mirror image caps
        tris = orient_cap(ladder(x, upper, lower, T(opts.split_tie_rel)), nodes)
    elseif mode == :zip && loop.closed
        # the two sides are paired by normalised arclength from the leading edge, not by nearest point
        pa, wa = arclength_projections(x, upper, lower)
        pb, wb = arclength_projections(x, lower, upper)
        pair = ChainPair(:zip, upper, lower, pa, pb, max(wa, wb))
    end
    actual = loop.closed ? mode : :open
    return Opening(nodes, le, te, upper, lower, actual, 1:0, planarity), tris, pair
end
