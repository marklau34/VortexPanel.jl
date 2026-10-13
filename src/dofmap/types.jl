"""
    DOFMapOptions(; cap_te_copy=:upper)

- `cap_te_copy`: which copy of a duplicated trailing edge node the cap triangles touching it use, `:upper` or `:lower`.
  The cap meets both surfaces at the tip's trailing edge node, so it can only be continuous with one of them; the
  jump μ_u - μ_l to the other runs as a filament along the cap's seam with that surface (doc §4 and §7, "wake tips").
"""
Base.@kwdef struct DOFMapOptions
    cap_te_copy::Symbol = :upper
end

function check_options(opts::DOFMapOptions)
    opts.cap_te_copy in (:upper, :lower) || throw(ArgumentError("cap_te_copy must be :upper or :lower, not $(opts.cap_te_copy)"))
    return opts
end

"""
Trailing edge nodes of the wake, the rows of the Kutta condition. Wake trailing edge node k sits on body node
`node[k]` at `x[k]`; `tangent` runs along τ, and `bisector` is the unit outward bisector of the trailing edge wedge
(the direction of a first wake row along the bisector). `lines` lists the nodes of each trailing edge line along τ;
their end nodes are where the wake's tip edges start (`closed` lines have none).
"""
struct WakeTE{T}
    node::Vector{Int}
    x::Vector{SVector{3,T}}
    tangent::Vector{SVector{3,T}}
    bisector::Vector{SVector{3,T}}
    lines::Vector{Vector{Int}}
    closed::Vector{Bool}
end

nwake_te(w::WakeTE) = length(w.node)

"""
    DOFMapping

Nodal doublet values on the body ("slots") in terms of the reduced unknowns q, μ_slots = C q, and the doublet strengths
of the wake's trailing edge nodes, μ_w = K q, from the nodal Kutta condition. This is where every topology decision
lives (doc §4).

Slots and corners
- Slot s ≤ nnodes is body node s; trailing edge nodes that carry both an upper and a lower surface get a second slot
  for the lower copy (`lower_slot[node]`, 0 for other nodes). `slot_node` and `slot_side` give the node and side of
  every slot.
- `corner_slot[t][k]` is the slot of corner k of triangle t: the trial and test functions of the Galerkin system are
  the hat functions of the slots.

Constraints and unknowns
- `constraint[s]` is `:free` for an unknown, or how the slot is constrained: `:tjunction` (a hanging node, linear along
  the long edge), `:gap` (a node of the finer side of a stitched gap, linear along the coarser side) or `:zip` (a node
  of one side of a zipped opening, linear along the other side at the same normalised arclength).
- `free_slot[q]` is the slot of unknown q and `q_of_slot[s]` the unknown of a free slot (0 if constrained).
- `C` (nslot × nq) maps the unknowns to all slots; its rows sum to one.

Wake
- `wake_te`: the wake's trailing edge nodes; `K` (nwake_te × nq) gives their doublet strength μ_w = μ_u - μ_l, with the
  upper and lower values interpolated along their own chains at non-matching trailing edges. Its rows sum to zero, so
  a uniform μ sheds nothing. `implicit[k]` flags rows that are unknowns of the implicit system (all, for a rigid wake).

Bodies
- `slot_area[s]` = ∫φ_s dA, `body_q[b]` the unknowns of body b and `body_closed[b]` whether it is closed (left open
  or zipped otherwise). Every body has a bordered constraint row (`bordered_rows`): a uniform μ on it induces nothing.
"""
struct DOFMapping{T}
    nslot::Int
    nq::Int
    slot_node::Vector{Int}
    slot_side::Vector{TESide}
    lower_slot::Vector{Int}
    corner_slot::Vector{SVector{3,Int}}
    constraint::Vector{Symbol}
    free_slot::Vector{Int}
    q_of_slot::Vector{Int}
    C::SparseMatrixCSC{T,Int}
    wake_te::WakeTE{T}
    K::SparseMatrixCSC{T,Int}
    implicit::BitVector
    slot_area::Vector{T}
    slot_body::Vector{Int}
    body_q::Vector{Vector{Int}}
    body_closed::Vector{Bool}
    log::Vector{String}
    options::DOFMapOptions
end

function Base.show(io::IO, dm::DOFMapping{T}) where {T}
    ncon = count(!=(:free), dm.constraint)
    print(io, "DOFMapping{$T}: $(dm.nslot) slots ($(count(>(0), dm.lower_slot)) lower trailing edge copies), $(dm.nq) unknowns, ",
        "$ncon constrained, $(nwake_te(dm.wake_te)) wake trailing edge nodes in $(length(dm.wake_te.lines)) line(s)")
end
