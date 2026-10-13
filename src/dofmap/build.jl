"""Triangle of every boundary edge (one triangle), keyed by its sorted node pair"""
function boundary_triangles(g::SurfaceGeometry)
    d = Dict{Tuple{Int,Int},Int}()
    for (e, nodes) in enumerate(g.edges.nodes)
        h1, h2 = g.edges.he[e]
        h2 == 0 && (d[Tuple(nodes)] = he_tri(h1))
    end
    return d
end

"""
    assign_slots(g, opts, log)

Duplicate every trailing edge node whose corners lie on both an upper and a lower surface (the nodes of matching
trailing edges, and nodes shared by both chains of a non-matching one): its lower copy is a new slot. Then give every
triangle corner its slot: lower surface corners at a duplicated node take the lower copy, cap corners the copy chosen
by `opts.cap_te_copy`, all others the node's own slot.
"""
function assign_slots(g::SurfaceGeometry, opts::DOFMapOptions, log::Vector{String})
    nn = length(g.x)
    has_up, has_low = falses(nn), falses(nn)
    for (t, cs) in enumerate(g.kutta.corner_side), k in 1:3
        cs[k] == TEUpper && (has_up[g.tris[t][k]] = true)
        cs[k] == TELower && (has_low[g.tris[t][k]] = true)
    end
    slot_node = collect(1:nn)
    slot_side = [has_up[v] ? TEUpper : has_low[v] ? TELower : NoTE for v in 1:nn]
    lower_slot = zeros(Int, nn)
    for v in 1:nn
        if has_up[v] && has_low[v]
            push!(slot_node, v)
            push!(slot_side, TELower)
            lower_slot[v] = length(slot_node)
        end
    end
    corner_slot = map(eachindex(g.tris)) do t
        SVector{3,Int}(ntuple(3) do k
            v = g.tris[t][k]
            lower_slot[v] == 0 && return v
            side = g.kutta.corner_side[t][k]
            side == TELower && return lower_slot[v]
            side == TECap && return opts.cap_te_copy == :lower ? lower_slot[v] : v
            side == NoTE && push!(log, "corner $k of triangle $t at trailing edge node $v has no side label; it takes the upper copy")
            return v
        end)
    end
    return slot_node, slot_side, lower_slot, corner_slot
end

"""
    collect_constraints(g, corner_slot, lower_slot, log)

Constraints slot => [(slot, weight)] (weights summing to one):
- T-junctions: the hanging node is linear along the long edge, between the long edge triangle's corner slots.
- Stitched gaps and zipped openings: the nodes of the side with more nodes are linear along the other side (between
  the corner slots of the other side's boundary triangle), at their projection; nodes on both sides already share a
  slot, and trailing edge nodes keep their two copies (so zipping a wing tip does not tie the Kutta jump).
A slot constrained twice keeps its first constraint.
"""
function collect_constraints(g::SurfaceGeometry{T}, corner_slot, lower_slot, log::Vector{String}) where {T}
    btri = boundary_triangles(g)
    cslot(t, v) = corner_slot[t][findfirst(==(v), g.tris[t])]
    cons = Dict{Int,Vector{Tuple{Int,T}}}()
    kind = Dict{Int,Symbol}()
    function add!(s, terms, k)
        if haskey(cons, s)
            push!(log, "slot $s is constrained twice ($(kind[s]) and $k); keeping the $(kind[s]) constraint")
            return
        end
        cons[s] = terms
        kind[s] = k
    end
    tie = sqrt(eps(T))
    for tj in g.tjunctions
        t = get(btri, minmax(tj.a, tj.b), 0)
        if t == 0
            push!(log, "T-junction long edge ($(tj.a), $(tj.b)) is not a boundary edge; hanging node $(tj.hanging) left free")
            continue
        end
        lower_slot[tj.hanging] > 0 && push!(log, "hanging node $(tj.hanging) is a trailing edge node; only its upper slot is constrained")
        add!(tj.hanging, [(cslot(t, tj.a), 1 - tj.t), (cslot(t, tj.b), tj.t)], :tjunction)
    end
    for cp in g.gaps
        a_master = length(unique(cp.chain_a)) <= length(unique(cp.chain_b))
        M, S, proj = a_master ? (cp.chain_a, cp.chain_b, cp.proj_b) : (cp.chain_b, cp.chain_a, cp.proj_a)
        masters = Set(M)
        seen = Set{Int}()
        for (i, v) in enumerate(S)
            (v in masters || v in seen || lower_slot[v] > 0) && continue
            push!(seen, v)
            seg, t = proj[i]
            tm = get(btri, minmax(M[seg], M[seg+1]), 0)
            if tm == 0
                push!(log, "$(cp.kind) segment ($(M[seg]), $(M[seg+1])) is not a boundary edge; node $v left free")
                continue
            end
            s1, s2 = cslot(tm, M[seg]), cslot(tm, M[seg+1])
            terms = t <= tie ? [(s1, one(T))] : t >= 1 - tie ? [(s2, one(T))] : [(s1, 1 - t), (s2, t)]
            add!(v, terms, cp.kind == :zip ? :zip : :gap)
        end
    end
    return cons, kind
end

"""
    tip_constraints!(cons, kind, wake_te, lower_slot, log)

For `tip_kutta = :zero`: tie the lower copy of the trailing edge node at each end of an open trailing edge line (a tip)
to its upper copy, so μ_u = μ_l there and the Kutta rows give μ_w = 0: the Kutta condition at the tip corner, with no
tip filament in the wake and no jump at the cap seam.
"""
function tip_constraints!(cons::Dict{Int,Vector{Tuple{Int,T}}}, kind, wake_te::WakeTE, lower_slot, log) where {T}
    for (line, closed) in zip(wake_te.lines, wake_te.closed)
        closed && continue
        for k in unique((first(line), last(line)))
            v = wake_te.node[k]
            if lower_slot[v] == 0
                push!(log, "tip trailing edge node $v has no lower copy; its μ_w is left free")
                continue
            end
            haskey(cons, lower_slot[v]) && continue
            cons[lower_slot[v]] = [(v, one(T))]
            kind[lower_slot[v]] = :tip
        end
    end
    return cons, kind
end

struct ConstraintCycle <: Exception
    slot::Int
end

"""Expansion of slot s over the free slots, through chains of constraints"""
function expand!(memo::Dict{Int,Vector{Tuple{Int,T}}}, cons, s::Int, visiting::Set{Int}) where {T}
    haskey(memo, s) && return memo[s]
    haskey(cons, s) || return (memo[s] = [(s, one(T))])
    s in visiting && throw(ConstraintCycle(s))
    push!(visiting, s)
    acc = Dict{Int,T}()
    for (r, w) in cons[s], (f, wf) in expand!(memo, cons, r, visiting)
        acc[f] = get(acc, f, zero(T)) + w * wf
    end
    delete!(visiting, s)
    return (memo[s] = sort!([(f, w) for (f, w) in acc]))
end

"""
    resolve_constraints(nslot, cons, log)

Expand every slot over the free slots. A constraint that closes a cycle is dropped (its slot becomes free) and logged.
"""
function resolve_constraints(::Type{T}, nslot::Integer, cons::Dict{Int,Vector{Tuple{Int,T}}}, kind::Dict{Int,Symbol}, log) where {T}
    while true
        memo = Dict{Int,Vector{Tuple{Int,T}}}()
        try
            for s in 1:nslot
                expand!(memo, cons, s, Set{Int}())
            end
            return memo
        catch e
            e isa ConstraintCycle || rethrow()
            push!(log, "constraints of slot $(e.slot) form a cycle; slot $(e.slot) left free")
            delete!(cons, e.slot)
            delete!(kind, e.slot)
        end
    end
end

"""
    build_dofmap(g::SurfaceGeometry; opts=DOFMapOptions())

Slots, constraints, the constraint matrix C and the Kutta rows K of a geometry (see `DOFMapping`).
"""
function build_dofmap(g::SurfaceGeometry{T}; opts::DOFMapOptions=DOFMapOptions()) where {T}
    check_options(opts)
    log = String[]
    slot_node, slot_side, lower_slot, corner_slot = assign_slots(g, opts, log)
    nslot = length(slot_node)

    cons, kind = collect_constraints(g, corner_slot, lower_slot, log)
    wake_te, U, L = wake_trailing_edge(g, corner_slot, nslot, log)
    opts.tip_kutta == :zero && tip_constraints!(cons, kind, wake_te, lower_slot, log)
    memo = resolve_constraints(T, nslot, cons, kind, log)
    constraint = [get(kind, s, :free) for s in 1:nslot]
    free_slot = findall(==(:free), constraint)
    q_of_slot = zeros(Int, nslot)
    q_of_slot[free_slot] .= eachindex(free_slot)
    nq = length(free_slot)
    I, J, V = Int[], Int[], T[]
    for s in 1:nslot, (f, w) in memo[s]
        push!(I, s)
        push!(J, q_of_slot[f])
        push!(V, w)
    end
    C = sparse(I, J, V, nslot, nq)

    K = (U - L) * C
    dropzeros!(K)

    slot_area = zeros(T, nslot)
    slot_body = zeros(Int, nslot)
    for (t, tri) in enumerate(g.tris)
        a = norm(cross(g.x[tri[2]] - g.x[tri[1]], g.x[tri[3]] - g.x[tri[1]])) / 6
        for s in corner_slot[t]
            slot_area[s] += a
            slot_body[s] = g.tri_body[t]
        end
    end
    body_q = [Int[] for _ in g.bodies]
    for (q, s) in enumerate(free_slot)
        slot_body[s] > 0 && push!(body_q[slot_body[s]], q)
    end
    body_closed = [b.closed for b in g.bodies]

    return DOFMapping(nslot, nq, slot_node, slot_side, lower_slot, corner_slot, constraint, free_slot, q_of_slot, C, wake_te,
        K, trues(nwake_te(wake_te)), slot_area, slot_body, body_q, body_closed, log, opts)
end
