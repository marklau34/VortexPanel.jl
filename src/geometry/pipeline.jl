"""Chain pair with its edges and partners renumbered from one half-edge structure to another (by node pairs)"""
function remap_pair(c::PairCandidate{T}, from::HalfEdges, to::HalfEdges) where {T}
    m(e) = edge_id(to, from.edge_nodes[e]...)
    return PairCandidate(c.chain_a, c.chain_b, m.(c.edges_a), m.(c.edges_b), [(m(e), m(f)) for (e, f) in c.partners], c.width)
end

"""Sorted node pairs of consecutive nodes of a chain"""
chain_edge_keys(chain) = [minmax(chain[k], chain[k+1]) for k in 1:length(chain)-1]

"""
    label_fan!(corner_side, tris, he, tri_kind, stop_edge, t0, v, came, side)

Label the corners at node v of the triangles reached from t0 by rotating around v away from edge `came`, until an edge
where `stop_edge` is true, a boundary, a cap triangle or t0 again.
"""
function label_fan!(corner_side, tris, he::HalfEdges, tri_kind, stop_edge, t0::Integer, v::Integer, came::Integer, side::TESide)
    t = t0
    for _ in 1:256
        kv = findfirst(==(v), tris[t])
        corner_side[t] = setindex(corner_side[t], side, kv)
        hv = he_index(t, kv)
        hp = he_prev(hv)
        hc = he.he_edge[hv] == came ? hp : hv
        ec = he.he_edge[hc]
        stop_edge(ec) && return
        g = he.twin[hc]
        g == 0 && return
        u = he_tri(g)
        (u == t0 || tri_kind[u] == CapTri) && return
        t, came = u, ec
    end
    return
end

"""
    build_geometry(mesh::MeshLibrary.Mesh; opts=GeometryOptions())
    build_geometry(nodes, faces; edge_tags, face_tags, opts=GeometryOptions())

Turn a surface mesh of triangles and quads into a closed, consistently and outwardly oriented triangulation with
classified edges (`SurfaceGeometry`). The mesh is repaired, never rejected, except for edges shared by more than two
faces (`NonManifoldEdgeError`) and non-orientable surfaces. Every repair is logged in `g.log`; see `print_report`.

Stages: ingest; weld; split quads; half-edges; T-junction and chain pair candidates (positions only); orient
(consistently through twins and candidate links, then outward by signed volume) and freeze; record T-junctions and chain
pairs; Kutta edges and non-matching trailing edges; openings; cap, zip or leave them open; final connectivity, trailing
edge fans, bodies and quality.
"""
build_geometry(mesh::MeshLibrary.Mesh; opts::GeometryOptions=GeometryOptions()) = build_geometry(ingest(mesh); opts=opts)

function build_geometry(nodes::AbstractVector, faces::AbstractVector; edge_tags=Dict{String,Vector{Tuple{Int,Int}}}(),
    face_tags=Dict{String,Vector{Int}}(), opts::GeometryOptions=GeometryOptions())
    return build_geometry(ingest(nodes, faces; edge_tags=edge_tags, face_tags=face_tags); opts=opts)
end

function build_geometry(raw::RawSurface{T}; opts::GeometryOptions=GeometryOptions()) where {T}
    check_options(opts)
    log = RepairEvent{T}[]

    # weld and split
    w = weld(raw, opts.weld_rel, log)
    x = w.x
    sp = split_faces(x, w, opts, log)
    tris = copy(sp.tris)

    # twins, bow-ties
    he1 = build_halfedges(x, tris)
    fans = vertex_fans(length(x), tris, he1.twin)
    for v in findall(>(1), fans)
        push!(log, RepairEvent(:bowtie, x[v], T(fans[v]), "node $v joins $(fans[v]) fans"))
    end

    # T-junction and chain pair candidates, then orientation
    tjc = find_tjunctions(x, he1, opts, log)
    tj_taken = tjunction_edges(tjc)
    pc1 = find_chain_pairs(x, he1, tj_taken, opts, log)
    ori = orient(x, tris, he1, orientation_links(x, tris, he1, tjc, pc1), log)
    for t in findall(ori.flip)
        tris[t] = SVector(tris[t][1], tris[t][3], tris[t][2])
    end
    he2 = build_halfedges(x, tris)
    pc = [remap_pair(c, he1, he2) for c in pc1]

    # record T-junctions
    tjunctions = [TJunction(h, c.a, c.b, c.t[k]) for c in tjc for (k, h) in enumerate(c.hanging)]
    tj_long = Set(minmax(c.a, c.b) for c in tjc)
    for tj in tjunctions
        push!(log, RepairEvent(:tjunction, x[tj.hanging], tj.t, "node $(tj.hanging) hangs on edge ($(tj.a), $(tj.b))"))
    end

    # Kutta edges, gaps and non-matching trailing edges
    tagged = welded_tag_edges(raw, w.input_node, opts.kutta_tags)
    kutta, mismatch = classify_kutta_edges(x, tris, he2, tagged, w.input_node, opts, log)
    gaps = ChainPair{T}[]
    kpairs = ChainPair{T}[]
    for c in pc
        kind = classify_pair(x, tris, he2, c, opts)
        pa, _ = chain_projections(x, c.chain_a, c.chain_b)
        pb, _ = chain_projections(x, c.chain_b, c.chain_a)
        cp = ChainPair(kind, c.chain_a, c.chain_b, pa, pb, c.width)
        push!(kind == :kutta ? kpairs : gaps, cp)
        push!(log, RepairEvent(kind == :kutta ? :kutta_pair : :gap, x[c.chain_a[1]], c.width,
            "$(kind == :kutta ? "non-matching trailing edge" : "gap") between chains of $(length(c.chain_a)) and $(length(c.chain_b)) nodes"))
    end
    gap_keys = Set(k for g in gaps for c in (g.chain_a, g.chain_b) for k in chain_edge_keys(c))
    kpair_keys = Set(k for g in kpairs for c in (g.chain_a, g.chain_b) for k in chain_edge_keys(c))

    # openings: cap, zip or leave open
    taken = union(tj_taken, gap_keys, kpair_keys)
    loops = find_openings(x, tris, he2, taken, opts, log)
    te_candidates = Set{Int}(v for p in kutta for v in p)
    for g in kpairs
        union!(te_candidates, g.chain_a, g.chain_b)
    end
    tri_kind = fill(BodyTri, length(tris))
    tri_parent = copy(sp.tri_parent)
    tri_slot = copy(sp.tri_slot)
    openings = Opening{T}[]
    for loop in loops
        op, captris, zp = close_opening(x, loop, te_candidates, opts.openings, opts)
        r = length(tris)+1:length(tris)+length(captris)
        append!(tris, captris)
        append!(tri_kind, fill(CapTri, length(captris)))
        append!(tri_parent, zeros(Int, length(captris)))
        append!(tri_slot, zeros(Int8, length(captris)))
        push!(openings, Opening(op.loop, op.le, op.te, op.chain_upper, op.chain_lower, op.mode, r, op.planarity))
        isnothing(zp) || push!(gaps, zp)
        push!(log, RepairEvent(Symbol(:opening_, op.mode), x[op.te], op.planarity,
            "opening of $(length(op.loop)) nodes ($(op.mode)$(op.mode == :cap ? ", $(length(captris)) triangles" : ""))"))
    end

    # final connectivity and edge classes
    he = build_halfedges(x, tris)
    ne = length(he.edge_nodes)
    class = Vector{EdgeClass}(undef, ne)
    tj_short = setdiff(tj_taken, tj_long)
    for e in 1:ne
        h1, h2 = he.edge_he[e]
        key = Tuple(he.edge_nodes[e])
        class[e] = if h2 == 0
            key in tj_long ? TJLong : key in tj_short ? TJShort : key in kpair_keys ? KuttaPair : key in gap_keys ? Gap : Open
        elseif key in kutta
            Kutta
        elseif tri_kind[he_tri(h1)] != tri_kind[he_tri(h2)]
            CapSeam
        else
            Interior
        end
    end

    # trailing edge chains and fans
    chains = kutta_chains(x, tris, he, kutta, opts)
    pair_chains = [kutta_pair_chain(x, tris, he, kp, k, opts) for (k, kp) in enumerate(kpairs)]
    corner_side = fill(SVector{3,TESide}(NoTE, NoTE, NoTE), length(tris))
    stop(e) = class[e] == Kutta || class[e] == KuttaPair
    for ch in chains, k in eachindex(ch.edges)
        for v in (ch.nodes[k], ch.nodes[k+1])
            label_fan!(corner_side, tris, he, tri_kind, stop, ch.tri_upper[k], v, ch.edges[k], TEUpper)
            label_fan!(corner_side, tris, he, tri_kind, stop, ch.tri_lower[k], v, ch.edges[k], TELower)
        end
    end
    for pch in pair_chains
        for (chain, ts, side) in ((pch.upper, pch.tri_upper, TEUpper), (pch.lower, pch.tri_lower, TELower)), k in eachindex(ts)
            e = edge_id(he, chain[k], chain[k+1])
            label_fan!(corner_side, tris, he, tri_kind, stop, ts[k], chain[k], e, side)
            label_fan!(corner_side, tris, he, tri_kind, stop, ts[k], chain[k+1], e, side)
        end
    end
    te_nodes = sort!(unique!(vcat(Int[], (c.nodes for c in chains)..., (p.union_nodes for p in pair_chains)...)))
    te_set = Set(te_nodes)
    cap_of = Dict{Int,Int}()
    for t in eachindex(tris)
        tri_kind[t] == CapTri || continue
        for k in 1:3
            v = tris[t][k]
            haskey(cap_of, v) || (cap_of[v] = t)
            v in te_set && (corner_side[t] = setindex(corner_side[t], TECap, k))
        end
    end
    open_nodes = Set(v for e in 1:ne if class[e] == Open for v in he.edge_nodes[e])
    pair_nodes = Set(v for p in pair_chains for v in p.union_nodes)
    kdeg = Dict{Int,Int}()
    for p in kutta, v in p
        kdeg[v] = get(kdeg, v, 0) + 1
    end
    function end_info(v)
        haskey(cap_of, v) && return :cap, cap_of[v]
        v in open_nodes && return :open, 0
        v in pair_nodes && return :pair, 0
        get(kdeg, v, 0) > 2 && return :branch, 0
        return :dangling, 0
    end
    chains = map(chains) do ch
        ch.closed && return ch
        (t1, c1), (t2, c2) = end_info(ch.nodes[1]), end_info(ch.nodes[end])
        KuttaChain(ch.nodes, ch.s, ch.edges, ch.tri_upper, ch.tri_lower, ch.closed, (t1, t2), (c1, c2))
    end

    # bodies: triangles connected through edges, T-junctions and chain pairs
    uf = UnionFind(length(tris))
    for h in eachindex(he.twin)
        he.twin[h] > h && unite!(uf, he_tri(h), he_tri(he.twin[h]))
    end
    touch = Set{Tuple{Int,Int}}()
    tri_of_edge(a, b) = (e = edge_id(he, a, b); e == 0 ? 0 : he_tri(he.edge_he[e][1]))
    for c in tjc
        tl = tri_of_edge(c.a, c.b)
        chain = [c.a; c.hanging; c.b]
        for i in 1:length(chain)-1
            ts = tri_of_edge(chain[i], chain[i+1])
            (tl == 0 || ts == 0) && continue
            unite!(uf, tl, ts)
            push!(touch, minmax(tl, ts))
        end
    end
    for c in pc, (e, f) in c.partners
        ta = tri_of_edge(he2.edge_nodes[e]...)
        tb = tri_of_edge(he2.edge_nodes[f]...)
        (ta == 0 || tb == 0) && continue
        unite!(uf, ta, tb)
        push!(touch, minmax(ta, tb))
    end
    open_or_zip = Set(v for op in openings if op.mode != :cap for v in op.loop)
    bodies = Body{T}[]
    tri_body = zeros(Int, length(tris))
    for (b, ts) in enumerate(components(uf))
        tri_body[ts] .= b
        nodes = sort!(unique!(reduce(vcat, [collect(tris[t]) for t in ts])))
        bedges = unique!(reduce(vcat, [[he.he_edge[he_index(t, k)] for k in 1:3] for t in ts]))
        closed = !any(class[e] == Open for e in bedges) && isdisjoint(nodes, open_or_zip)
        centre = sum(x[nodes]) / length(nodes)
        push!(bodies, Body(ts, nodes, closed, length(nodes) - length(bedges) + length(ts), signed_volume(x, tris, ts, centre)))
    end

    # parent quads in the final orientation
    first_half = Dict{Int,Int}()
    for t in eachindex(tri_parent)
        tri_slot[t] == 1 && (first_half[tri_parent[t]] = t)
    end
    quads = map(eachindex(sp.quads)) do q
        a, b, c, d = sp.quads[q]
        half = sp.quad_diag[q] == 1 ? (a, b, c) : (a, b, d)
        t = tris[first_half[sp.quad_parent[q]]]
        same = any(Tuple(t) == (half[mod1(k, 3)], half[mod1(k + 1, 3)], half[mod1(k + 2, 3)]) for k in 1:3)
        same ? sp.quads[q] : SVector(a, d, c, b)
    end

    # quality and counts
    qmin, qmax, asp, skew, area = quality_arrays(x, tris, opts)
    sliver = BitVector((qmin .< T(opts.sliver_angle_deg)) .| (asp .> T(opts.sliver_aspect)))
    nevent(k) = count(ev -> ev.kind == k, log)
    counts = Dict{Symbol,Int}(
        :nodes => length(x), :triangles => length(tris), :cap_triangles => count(==(CapTri), tri_kind), :edges => ne,
        :bodies => length(bodies), :closed_bodies => count(b -> b.closed, bodies), :input_nodes => length(raw.x),
        :input_elements => length(raw.faces), :welded_nodes => nevent(:weld),
        :dropped_elements => nevent(:drop_degenerate) + nevent(:drop_duplicate), :collapsed_elements => nevent(:collapse),
        :split_shortest => count(==(:shortest), sp.quad_reason), :split_fold => count(==(:fold), sp.quad_reason),
        :split_tie => count(r -> r in (:tie_mirror, :tie_id), sp.quad_reason), :split_fixed => count(==(:fixed), sp.quad_reason),
        :flipped_triangles => count(ori.flip), :bowtie_vertices => nevent(:bowtie), :tjunctions => length(tjunctions),
        :gaps => count(g -> g.kind == :gap, gaps), :kutta_pairs => length(kpairs),
        :openings_capped => count(op -> op.mode == :cap, openings), :openings_zipped => count(op -> op.mode == :zip, openings),
        :openings_open => count(op -> op.mode == :open, openings), :kutta_edges => length(kutta), :kutta_chains => length(chains),
        :te_nodes => length(te_nodes), :kutta_tag_mismatch => mismatch, :slivers => count(sliver))
    wmax(k) = maximum((ev.magnitude for ev in log if ev.kind == k); init=zero(T))
    stats = Dict{Symbol,T}(:min_angle => minimum(qmin; init=T(Inf)), :max_angle => maximum(qmax; init=zero(T)),
        :max_aspect => maximum(asp; init=zero(T)), :max_skew => maximum(skew; init=zero(T)), :max_weld => wmax(:weld),
        :max_gap => maximum((g.width for g in gaps if g.kind == :gap); init=zero(T)),
        :max_cap_planarity => maximum((op.planarity for op in openings); init=zero(T)))
    quality = QualityReport(qmin, qmax, asp, skew, area, sliver, counts, stats)

    edges = EdgeTable(he.edge_nodes, he.edge_he, class)
    kdata = KuttaData(chains, kpairs, pair_chains, te_nodes, corner_side, mismatch)
    return SurfaceGeometry(x, tris, tri_kind, tri_parent, tri_slot, tri_body, quads, sp.quad_parent, sp.quad_diag, he.twin,
        he.he_edge, edges, w.input_node, bodies, tjunctions, gaps, openings, kdata, sort!(collect(touch)), quality, log,
        raw.edge_tags, raw.face_tags, opts)
end
