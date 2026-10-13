"""
    solve_galerkin(dm, A, b; W=nothing)

Reduce the system to the unknowns, (Cᵀ A [C; W]) q = Cᵀ b, border it with one row per closed body (doc §8),

    [CᵀA C  e] [q]   [Cᵀb]
    [cᵀ     0] [λ] = [ 0 ],

c = Cᵀ∫φ dA pinning the body's area-weighted mean doublet to zero and e one on its unknowns, and solve it by dense LU.
Both border vectors are scaled to a 2-norm equal to the mean diagonal entry of the matrix, so the bordered matrix is as
well conditioned as the operator away from its null space (which leaves q unchanged; λ is reported unscaled). Returns
the unknowns `q`, the slot values `μ = C q`, the multipliers `λ` (the per-row incompatibility of the right hand side),
the 1-norm condition number estimate `cond` of the bordered matrix, and the relative residual.
"""
function solve_galerkin(dm::DOFMapping, A::AbstractMatrix, b::AbstractVector; W=nothing)
    Ar, br = reduce_system(dm, A, b; W=W)
    rows = bordered_rows(dm)
    nq, nb = dm.nq, length(rows)
    T = promote_type(eltype(Ar), eltype(br))
    M = zeros(T, nq + nb, nq + nb)
    M[1:nq, 1:nq] .= Ar
    scale = sum(abs, diag(Ar)) / nq
    se = zeros(T, nb)
    for (k, r) in enumerate(rows)
        se[k] = scale / norm(r.e)
        M[1:nq, nq+k] .= se[k] .* r.e
        M[nq+k, 1:nq] .= r.c .* (scale / norm(r.c))
    end
    rhs = [Vector{T}(br); zeros(T, nb)]
    F = lu(M)
    rcond = T <: LinearAlgebra.BlasFloat ? LAPACK.gecon!('1', copy(F.factors), opnorm(M, 1)) : inv(cond(M, 1))
    x = F \ rhs
    q = x[1:nq]
    return (q=q, μ=dm.C * q, λ=se .* x[nq+1:end], cond=1 / rcond, residual=norm(M * x - rhs) / norm(rhs))
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
