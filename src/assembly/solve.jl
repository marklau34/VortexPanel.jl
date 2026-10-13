"""
    bordered_system(dm, A, b; W=nothing)

The system `solve_galerkin` factors: the reduction to the unknowns, (Cᵀ A [C; W]) q = Cᵀ b, bordered with one row per
body (doc §8, `bordered_rows`),

    [CᵀA C  e] [q]   [Cᵀb]
    [cᵀ     0] [λ] = [ 0 ],

c = Cᵀ∫φ dA pinning the body's area-weighted mean doublet to zero and e one on its unknowns. Both border vectors are
scaled to a 2-norm equal to the mean diagonal entry of the matrix, so the bordered matrix is as well conditioned as the
operator away from its null space (which leaves q unchanged). Returns the matrix, the right hand side and the scales
of the border columns (λ is the scale times the solution's border entry).
"""
function bordered_system(dm::DOFMapping, A::AbstractMatrix, b::AbstractVector; W=nothing)
    M, se = bordered_matrix(dm, A; W=W)
    T = promote_type(eltype(M), eltype(b), eltype(dm.C))
    T === eltype(M) || (M = Matrix{T}(M); se = Vector{T}(se))
    return M, bordered_rhs(dm, b, T), se
end

"""The matrix of `bordered_system` and the scales of its border columns"""
function bordered_matrix(dm::DOFMapping, A::AbstractMatrix; W=nothing)
    Ar, _ = reduce_system(dm, A, zeros(eltype(A), dm.nslot); W=W)
    rows = bordered_rows(dm)
    nq, nb = dm.nq, length(rows)
    T = eltype(Ar)
    M = zeros(T, nq + nb, nq + nb)
    M[1:nq, 1:nq] .= Ar
    scale = sum(abs, diag(Ar)) / nq
    se = zeros(T, nb)
    for (k, r) in enumerate(rows)
        se[k] = scale / norm(r.e)
        M[1:nq, nq+k] .= se[k] .* r.e
        M[nq+k, 1:nq] .= r.c .* (scale / norm(r.c))
    end
    return M, se
end

"""The right hand side of `bordered_system`, [Cᵀb; 0]"""
bordered_rhs(dm::DOFMapping, b::AbstractVector, ::Type{T}=promote_type(eltype(dm.C), eltype(b))) where {T} =
    [Vector{T}(dm.C' * b); zeros(T, length(dm.body_closed))]

"""
    solve_galerkin(dm, A, b; W=nothing)

Solve the bordered system of `bordered_system` by dense LU. Returns the unknowns `q`, the slot values `μ = C q`, the
multipliers `λ` (the per-row incompatibility of the right hand side), the 1-norm condition number estimate `cond` of
the bordered matrix, and the relative residual.
"""
function solve_galerkin(dm::DOFMapping, A::AbstractMatrix, b::AbstractVector; W=nothing)
    M, rhs, se = bordered_system(dm, A, b; W=W)
    T, nq = eltype(M), dm.nq
    F = lu(M)
    rcond = T <: LinearAlgebra.BlasFloat ? LAPACK.gecon!('1', copy(F.factors), opnorm(M, 1)) : inv(cond(M, 1))
    x = F \ rhs
    q = x[1:nq]
    return (q=q, μ=dm.C * q, λ=se .* x[nq+1:end], cond=1 / rcond, residual=norm(M * x - rhs) / norm(rhs))
end

"""
The LU factorization of a bordered system (`factorize_galerkin`), with what `solve_galerkin(F, b)` needs to map a
right hand side to the unknowns
"""
struct GalerkinFactorization{T,L,M}
    lu::L
    se::Vector{T}
    C::M
    nq::Int
    nb::Int
    cond::T
end

"""
    factorize_galerkin(dm, A; W=nothing)

Factor the bordered matrix of `bordered_system` once, for many right hand sides: a body at rest keeps its matrix
while an explicit wake and particles change the right hand side every step (`explicit_rhs`). Solve with
`solve_galerkin(F, b)`.
"""
function factorize_galerkin(dm::DOFMapping, A::AbstractMatrix; W=nothing)
    M, se = bordered_matrix(dm, A; W=W)
    T = eltype(M)
    F = lu(M)
    rcond = T <: LinearAlgebra.BlasFloat ? LAPACK.gecon!('1', copy(F.factors), opnorm(M, 1)) : inv(cond(M, 1))
    return GalerkinFactorization(F, se, dm.C, dm.nq, length(se), 1 / rcond)
end

"""
    solve_galerkin(F::GalerkinFactorization, b)

Solve the factored bordered system for the right hand side b (over the slots, as `assemble` returns it). Returns `q`,
`μ = C q`, the multipliers `λ` and the condition number estimate `cond` of the factored matrix.
"""
function solve_galerkin(F::GalerkinFactorization{T}, b::AbstractVector) where {T}
    x = F.lu \ [Vector{T}(F.C' * b); zeros(T, F.nb)]
    q = x[1:F.nq]
    return (q=q, μ=F.C * q, λ=F.se .* x[F.nq+1:end], cond=F.cond)
end

"""
    consistency(dm, A, b; W=nothing)

Quadrature diagnostics of the assembled system on closed bodies, in terms of A C_full with C_full = [C; W] (W the wake
rows when A has wake columns). With exact integration every row sums to zero (a uniform doublet sheds no vortex sheet,
and its Kutta rows vanish; exact by construction), every column sums to zero (the field of an unknown's vortex sheets
and filaments, continuous through T-junctions and stitches, is solenoidal and has no net flux through the closed
surface), and the entries of b sum to zero (no net source flux). Returns the relative row and column sum defects and
the relative sum of b. A gap of relative width ε adds O(ε log ε) to the columns at its stitch: the flux through the slit.
"""
function consistency(dm::DOFMapping, A::AbstractMatrix, b::AbstractVector; W=nothing)
    ACf = A * (isnothing(W) ? dm.C : vcat(dm.C, W))
    nA = norm(ACf, Inf)
    return (rows=maximum(abs, sum(ACf, dims=2)) / nA, columns=maximum(abs, sum(ACf, dims=1)) / nA, rhs=abs(sum(b)) / sum(abs, b))
end
