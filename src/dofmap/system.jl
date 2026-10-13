"""Doublet values at every slot, μ = C q"""
slot_values(dm::DOFMapping, q::AbstractVector) = dm.C * q

"""Doublet strengths of the wake's trailing edge nodes, μ_w = K q"""
kutta_values(dm::DOFMapping, q::AbstractVector) = dm.K * q

"""Unknowns from slot values (the values at the free slots)"""
restrict(dm::DOFMapping, μ::AbstractVector) = μ[dm.free_slot]

"""
    wake_rows(dm, te_of_wake_node)

Rows of the constraint matrix for the wake nodes: a wake node carries the doublet strength of the trailing edge node
its strip starts from (`te_of_wake_node[i]`), so its row is that row of K. The full constraint matrix of a system
whose columns are the body slots followed by the wake nodes is `[dm.C; wake_rows(dm, te_of_wake_node)]`.
"""
wake_rows(dm::DOFMapping, te_of_wake_node::AbstractVector{<:Integer}) = dm.K[te_of_wake_node, :]

"""
    reduce_system(dm, A, b; W=nothing)

Galerkin reduction to the unknowns: (Cᵀ A C_full) q = Cᵀ b, where A has a row per slot (the test functions) and a
column per slot followed by a column per wake node, and C_full = [C; W] with W = `wake_rows(...)` (C when A has no wake
columns). Applying Cᵀ on the test side keeps the method a Galerkin method.
"""
function reduce_system(dm::DOFMapping, A::AbstractMatrix, b::AbstractVector; W=nothing)
    Cf = isnothing(W) ? dm.C : vcat(dm.C, W)
    size(A) == (dm.nslot, size(Cf, 1)) || throw(DimensionMismatch("A is $(size(A)), expected ($(dm.nslot), $(size(Cf, 1)))"))
    return dm.C' * (A * Cf), dm.C' * b
end

"""
    bordered_rows(dm)

For every body b, the border vectors of its constraint row (doc §8): `e` is one on the body's unknowns and
`c = Cᵀ(∫φ dA)` restricted to the body, so that cᵀq is the integral of μ over the body; the system becomes
[CᵀAC e; cᵀ 0][q; λ] = [Cᵀb; 0], pinning the area-weighted mean of μ to zero. A uniform μ on a body sheds no vortex
sheet and no wake (its Kutta rows vanish), and the solver models no filaments on the body, so every body needs the
row, closed or not (an open or zipped body would leave uncancelled filaments only in the doublet-layer picture).
"""
function bordered_rows(dm::DOFMapping{T}) where {T}
    rows = NamedTuple{(:body, :e, :c),Tuple{Int,Vector{T},Vector{T}}}[]
    for b in eachindex(dm.body_closed)
        e = zeros(T, dm.nq)
        e[dm.body_q[b]] .= 1
        c = Vector(dm.C' * (dm.slot_area .* (dm.slot_body .== b)))
        push!(rows, (body=b, e=e, c=c))
    end
    return rows
end

"""Position of the node of every slot"""
slot_points(dm::DOFMapping, g::SurfaceGeometry) = g.x[dm.slot_node]

"""
    reproduction_error(dm, g, f)

Largest |(C q)[s] - f(x_s)| over the slots when q holds f at the free slots, with x_s the slot's node. Zero up to
rounding for a linear f through T-junctions and stitches of coincident chains; for a gap of width w it is at most
|∇f| w (the stitched node takes the value at its projection on the other side).
"""
function reproduction_error(dm::DOFMapping, g::SurfaceGeometry, f)
    x = slot_points(dm, g)
    q = [f(x[s]) for s in dm.free_slot]
    return maximum(abs.(dm.C * q .- f.(x)))
end

"""
    stitch_offset(dm, g)

Largest distance between a constrained slot's node and the point its constraint interpolates at (zero for T-junctions
and coincident stitches; the gap width otherwise; the section thickness for zipped tips).
"""
function stitch_offset(dm::DOFMapping{T}, g::SurfaceGeometry{T}) where {T}
    x = slot_points(dm, g)
    xf = x[dm.free_slot]
    anchor = dm.C * reduce(hcat, xf)'
    d = [norm(SVector{3,T}(anchor[s, :]) - x[s]) for s in 1:dm.nslot if dm.constraint[s] != :free]
    return maximum(d; init=zero(T))
end

"""
    print_dofmap(io=stdout, dm[, g])

Slots, unknowns, constraints by kind, wake trailing edge nodes and bordered rows of a `DOFMapping` (and the largest
stitch offset when the geometry is given).
"""
function print_dofmap(io::IO, dm::DOFMapping, g::Union{Nothing,SurfaceGeometry}=nothing)
    println(io, "DOF map")
    @printf(io, "  %-22s %d (%d nodes + %d lower trailing edge copies)\n", "slots", dm.nslot, count(s -> s == dm.slot_node[s], 1:dm.nslot),
        count(>(0), dm.lower_slot))
    @printf(io, "  %-22s %d free; constrained: T-junction %d, gap %d, zip %d\n", "unknowns", dm.nq,
        count(==(:tjunction), dm.constraint), count(==(:gap), dm.constraint), count(==(:zip), dm.constraint))
    isnothing(g) || @printf(io, "  %-22s %.3g\n", "max stitch offset", stitch_offset(dm, g))
    if nwake_te(dm.wake_te) == 0
        @printf(io, "  %-22s none\n", "wake trailing edge")
    else
        @printf(io, "  %-22s %d nodes in %d line(s) of %s nodes\n", "wake trailing edge", nwake_te(dm.wake_te), length(dm.wake_te.lines),
            join(length.(dm.wake_te.lines), ", "))
    end
    @printf(io, "  %-22s %d of %d closed (each has a bordered row); cap corners at trailing edge nodes use the %s copy\n", "bodies",
        count(dm.body_closed), length(dm.body_closed), dm.options.cap_te_copy)
    for m in dm.log
        println(io, "  note: ", m)
    end
    return nothing
end

print_dofmap(dm::DOFMapping, g=nothing) = print_dofmap(stdout, dm, g)
