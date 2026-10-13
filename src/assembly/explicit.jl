# Explicit (known) vortex elements in the right hand side: the free wake's panels and filaments and its particles, which
# couple into the body through their velocity while the matrix (body and trailing edge panel) stays factored.

"""
    ExplicitSources(sheets, sheet_μ, segments, segment_Γ, blobs)

Vortex elements of known strength for `explicit_rhs`: flat triangles carrying the vortex sheet of their corner doublets
`sheet_μ`, straight filaments with circulations `segment_Γ` at their two ends (`Segment`), and vortex particles (`Blob`)
"""
struct ExplicitSources{F}
    sheets::Vector{FlatTriangle{F}}
    sheet_μ::Vector{SVector{3,F}}
    segments::Vector{Segment{F}}
    segment_Γ::Vector{SVector{2,F}}
    blobs::Vector{Blob{F}}
end

ExplicitSources{F}() where {F} = ExplicitSources(FlatTriangle{F}[], SVector{3,F}[], Segment{F}[], SVector{2,F}[], Blob{F}[])

Elements.distance_to(b::Blob, p::SVector{3}) = norm(p - b.x)

"""
    explicit_rhs(els, dm, src; opts=AssemblyOptions())

The right hand side of explicit vortex elements: -∫ψ_i n·u_src dA for every slot i, to be added to `assemble`'s b. The
sheets and filaments take the assembly's own tiers (`wake_pair_rule`, `segment_rule`: contact and near rules where they
come close), so a known wake enters exactly as its assembled columns times its doublets would; the particles take the
far rule when they are more than `far_ratio` times the larger of the test triangle's diameter and their core away, and
the near rule otherwise. Nothing is stored and the matrix is untouched.
"""
function explicit_rhs(els::AbstractVector{<:FlatTriangle}, dm::DOFMapping, src::ExplicitSources; opts::AssemblyOptions=AssemblyOptions())
    check_options(opts)
    nt = length(els)
    rules = tier_rules(opts)
    F = eltype(els[1].n)
    b = zeros(F, dm.nslot)
    lk = ReentrantLock()
    nchunk = opts.threaded ? min(nt, 8 * Threads.nthreads()) : 1
    chunks = [collect(c) for c in Iterators.partition(1:nt, cld(nt, nchunk))]
    tasks = map(chunks) do chunk
        task = () -> begin
            ufar = zeros(SVector{3,F}, length(rules.far.w))
            for t in chunk
                T = els[t]
                rows = zero(SVector{3,F})
                for (S, μS) in zip(src.sheets, src.sheet_μ)
                    _, rule = wake_pair_rule(T, S, rules, opts)
                    blk, _ = pair_block(T, S, rule, zero(F), opts.test, 0, opts.kernel_far)
                    rows -= blk * μS
                end
                for (f, Γ) in zip(src.segments, src.segment_Γ)
                    _, rule = segment_rule(T, f, rules, opts)
                    ca, cb = segment_block(T, f, rule, opts.test)
                    rows -= ca * Γ[1] + cb * Γ[2]
                end
                fill!(ufar, zero(SVector{3,F}))
                for bl in src.blobs
                    if norm(bl.x - T.centroid) > opts.far_ratio * max(T.diam, bl.δ)
                        for (i, λ) in enumerate(rules.far.λ)
                            ufar[i] += blob_velocity(bl, point(T, λ))
                        end
                    else
                        rule = near_rule(T, bl, rules.near, opts.near_ratio, opts.near_depth)
                        for (λ, wq) in zip(rule.λ, rule.w)
                            rows -= (wq * T.area * dot(blob_velocity(bl, point(T, λ)), T.n)) * test_weights(opts.test, λ)
                        end
                    end
                end
                if !isempty(src.blobs)
                    for (i, (λ, wq)) in enumerate(zip(rules.far.λ, rules.far.w))
                        rows -= (wq * T.area * dot(ufar[i], T.n)) * test_weights(opts.test, λ)
                    end
                end
                ct = dm.corner_slot[t]
                lock(lk) do
                    for a in 1:3
                        b[ct[a]] += rows[a]
                    end
                end
            end
        end
        opts.threaded ? Threads.@spawn(task()) : (task(); nothing)
    end
    opts.threaded && foreach(wait, tasks)
    return b
end
