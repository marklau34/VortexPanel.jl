# Structured perturbations of a PAGE wing mesh (`PAGE.generate_mesh(::ParametricWing, ::LiftSurfMeshParams)`,
# split_te=false): upper surface elements first, chordwise index fastest from the leading edge, strips in +y. Nodes move
# only along the span, which keeps the geometry of a prismatic (untwisted, untapered, unswept) wing exact.

"""Node of chordwise row k (1 at the leading edge) and span station j on the upper and lower surfaces"""
struct WingGrid
    nspan::Int
    nchord::Int
    upper::Matrix{Int}
    lower::Matrix{Int}
end

"""
    wing_grid(mesh; nspan, nchord)

Structured node map of a PAGE wing mesh: upper element (j-1)·nchord + i has nodes (i,j), (i+1,j), (i+1,j+1), (i,j+1)
and the lower element nspan·nchord + (j-1)·nchord + i has them reversed.
"""
function wing_grid(mesh::MeshLibrary.Mesh; nspan::Integer, nchord::Integer)
    nsurf = nspan * nchord
    MeshLibrary.nelements(mesh) >= 2nsurf || throw(ArgumentError("mesh has fewer than 2·nspan·nchord elements"))
    U = zeros(Int, nchord + 1, nspan + 1)
    L = zeros(Int, nchord + 1, nspan + 1)
    for j in 1:nspan, i in 1:nchord
        c = MeshLibrary.get_nodes(mesh.elements[(j-1)*nchord+i])
        U[i, j], U[i+1, j], U[i+1, j+1], U[i, j+1] = c[1], c[2], c[3], c[4]
        c = MeshLibrary.get_nodes(mesh.elements[nsurf+(j-1)*nchord+i])
        L[i, j+1], L[i+1, j+1], L[i+1, j], L[i, j] = c[1], c[2], c[3], c[4]
    end
    return WingGrid(nspan, nchord, U, L)
end

"""
    skew_shifts(pattern, σ, nk, nj; rng, rows, te_rows=false, xc=nothing, band_start=0.8, shared_te=true)

Spanwise shifts (as fractions of the smaller adjacent station spacing) of every upper and lower node (row k, station
j), with the semantics of `skew_pattern` in MESA.jl/dev/skew_sensitivity.jl:
- `:lean_same`: σ on both surfaces; `:lean_opposite`: σ on the upper and -σ on the lower surface
- `:jitter_same`: uniform in [-σ, σ], the same on both surfaces; `:jitter_independent`: drawn separately per surface
- `:te_misaligned` (`te_misaligned_mesh`): rows aft of chord fraction `band_start` (`xc` of each row) move by +σ/2 and
  -σ/2 on the two surfaces with the trailing edge shared (`shared_te=true`), or the upper rows by σ including the
  trailing edge row (`shared_te=false`, a non-matching trailing edge).
Only rows with `rows[k]` move: by default all but the leading and trailing edge rows, plus the trailing edge row with
`te_rows=true`. End stations never move.
"""
function skew_shifts(pattern::Symbol, σ::Real, nk::Integer, nj::Integer; rng::AbstractRNG=StableRNG(1),
    rows::Union{Nothing,AbstractVector{Bool}}=nothing, te_rows::Bool=false, xc=nothing, band_start::Real=0.8, shared_te::Bool=true)
    mask = isnothing(rows) ? [1 < k < nk || (te_rows && k == nk) for k in 1:nk] : collect(rows)
    up = zeros(nk, nj)
    low = zeros(nk, nj)
    if pattern == :lean_same
        up[mask, :] .= σ
        low[mask, :] .= σ
    elseif pattern == :lean_opposite
        up[mask, :] .= σ
        low[mask, :] .= -σ
    elseif pattern in (:jitter_same, :jitter_independent)
        r_up = σ .* (2 .* rand(rng, nk, nj) .- 1)
        r_low = pattern == :jitter_same ? r_up : σ .* (2 .* rand(rng, nk, nj) .- 1)
        up[mask, :] .= r_up[mask, :]
        low[mask, :] .= r_low[mask, :]
    elseif pattern == :te_misaligned
        isnothing(xc) && throw(ArgumentError(":te_misaligned needs the chord fraction xc of every row"))
        band = [k > 1 && xc[k] >= band_start - 1e-12 && !(shared_te && k == nk) for k in 1:nk]
        up[band, :] .= shared_te ? σ / 2 : σ
        shared_te && (low[band, :] .= -σ / 2)
    else
        throw(ArgumentError("Unknown skew pattern $pattern"))
    end
    up[:, [1, nj]] .= 0
    low[:, [1, nj]] .= 0
    return up, low
end

"""
    skew_wing(mesh, pattern, σ; nspan, nchord, rng=StableRNG(1), kwargs...)

Wing mesh with its nodes moved along the span by `skew_shifts` (see there for the patterns and keywords). The leading
and trailing edge nodes are shared by both surfaces in the input; where the two surfaces shift a shared node
differently, the lower surface gets its own copy (a non-matching trailing edge with `te_rows=true`, which MESA's
classical trailing edge cannot take).
"""
function skew_wing(mesh::MeshLibrary.Mesh, pattern::Symbol, σ::Real; nspan::Integer, nchord::Integer, rng::AbstractRNG=StableRNG(1), kwargs...)
    grid = wing_grid(mesh; nspan=nspan, nchord=nchord)
    X = coords(mesh)
    faces = elem_nodes(mesh)
    nk, nj = nchord + 1, nspan + 1
    ys = [X[grid.upper[1, j]][2] for j in 1:nj]
    xc = let x = [X[grid.upper[k, 1]][1] for k in 1:nk]
        (x .- x[1]) ./ (x[end] - x[1])
    end
    up, low = skew_shifts(pattern, σ, nk, nj; rng=rng, xc=xc, kwargs...)
    spacing(j) = min(ys[j] - ys[j-1], ys[j+1] - ys[j])
    log = AbuseLog()
    nsurf = nspan * nchord
    Y = copy(X)
    # lower surface copies of shared nodes that the two surfaces move differently
    for j in 2:nj-1, k in 1:nk
        u, l = grid.upper[k, j], grid.lower[k, j]
        if u == l && up[k, j] != low[k, j]
            push!(Y, X[u])
            c = length(Y)
            for e in nsurf+1:2nsurf
                faces[e] = [v == u ? c : v for v in faces[e]]
            end
            grid.lower[k, j] = c
            record!(log, :split_nodes, "row $k station $j: the lower surface gets its own copy of node $u")
        end
    end
    P = copy(Y) # positions before the shifts, including the new copies
    for j in 2:nj-1, k in 1:nk
        u, l = grid.upper[k, j], grid.lower[k, j]
        Y[u] = P[u] + SVector(0.0, up[k, j] * spacing(j), 0.0)
        l != u && (Y[l] = P[l] + SVector(0.0, low[k, j] * spacing(j), 0.0))
    end
    record!(log, :skewed_nodes, "$pattern skew of amplitude $σ"; n=count(!iszero, up) + count(!iszero, low))
    return rebuild(mesh, Y, faces), log
end

"""
    sliver_stations(mesh, n, ratio; nspan, nchord, rng)

Move n interior span stations (not next to a tip, no two adjacent) of a wing toward their outboard neighbour until
their spacing is `ratio` of what it was: conforming strips of sliver panels, with the panel count unchanged.
"""
function sliver_stations(mesh::MeshLibrary.Mesh, n::Integer, ratio::Real; nspan::Integer, nchord::Integer, rng::AbstractRNG=StableRNG(1))
    grid = wing_grid(mesh; nspan=nspan, nchord=nchord)
    X = coords(mesh)
    nj = nspan + 1
    picked = Int[]
    for j in shuffle(rng, collect(3:nj-2))
        length(picked) == n && break
        any(p -> abs(p - j) <= 1, picked) && continue
        push!(picked, j)
    end
    log = AbuseLog()
    Y = copy(X)
    for j in picked
        dy = (1 - ratio) * (X[grid.upper[1, j+1]][2] - X[grid.upper[1, j]][2])
        for v in unique([grid.upper[:, j]; grid.lower[:, j]])
            Y[v] = X[v] + SVector(0.0, dy, 0.0)
        end
        record!(log, :slivers, "station $j moved to $ratio of its spacing to station $(j + 1)"; n=2nchord)
    end
    return with_nodes(mesh, Y), log
end

"""
    naca_projector(; chord, thickness, x_le=0.0)

Projection onto an untwisted NACA 00xx wing with its chord along x from `x_le` and the section symmetric about z = 0
(PAGE's sharp trailing edge form): x and y are kept and z is set to the surface on the side of the point.
"""
function naca_projector(; chord::Real, thickness::Real, x_le::Real=0.0)
    function project(p)
        ξ = clamp((p[1] - x_le) / chord, 0.0, 1.0)
        yt = 5thickness * (0.2969sqrt(ξ) - 0.1260ξ - 0.3516ξ^2 + 0.2843ξ^3 - 0.1036ξ^4)
        return SVector(p[1], p[2], (p[3] >= 0 ? 1 : -1) * chord * yt)
    end
    return project
end
