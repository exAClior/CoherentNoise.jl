using Random, LinearAlgebra
using Yao

struct QuantumChannel{T<:AbstractFloat, MT<:AbstractMatrix{Complex{T}}}
    n::Int
    kraus_ops::Vector{MT}
end

function damping_channel(p::T) where {T<:AbstractFloat}
    return QuantumChannel(
        1, [Complex{T}[1.0 0.0; 0.0 sqrt(1.0 - p)], Complex{T}[0.0 sqrt(p); 0.0 0.0]]
    )
end

function unitary_channel(gate::Matrix{Complex{T}}) where {T<:AbstractFloat}
    isunitary(gate) || throw(ArgumentError("Input gate is not unitary"))
    n = Yao.log2i(size(gate, 1))
    return QuantumChannel(n, [gate])
end

function pauli_channel(probs::Vector{T}) where {T<:AbstractFloat}
    iseven(Yao.log2i(length(probs))) ||
        throw(ArgumentError("Input vector length is not a power of 4"))
    sum(probs) ≈ 1.0 || throw(ArgumentError("Input vector is not normalized"))

    n = Yao.log2i(length(probs)) ÷ 2

    return QuantumChannel(
        n,
        [
            sqrt(p) * mat(kron(pauli...)) for
            (p, pauli) in zip(probs, Iterators.product(repeat([[I2, X, Y, Z]], n)...))
        ],
    )
end

function apply(
    qc::QuantumChannel{T}, ρ::MT
) where {T<:AbstractFloat,MT<:AbstractMatrix{Complex{T}}}
    if size(ρ) != (2^qc.n, 2^qc.n)
        throw(DimensionMismatch("Dimension mismatch"))
    end
    ρ_new = mapreduce(x -> x * ρ * x', .+, qc.kraus_ops)
    return ρ_new
end

struct PauliTransferMatrix{T<:AbstractFloat}
    n::Int
    mtx::Matrix{Complex{T}}
end

function PauliTransferMatrix(qc::QuantumChannel{T}) where {T}
    n = qc.n
    paulis = [I2, X, Y, Z]
    mtx = zeros(Complex{T}, 4^n, 4^n)
    for ii in 0:(4^n - 1), jj in 0:(4^n - 1)
        pauli_ii = kron([paulis[kk + 1] for kk in digits(ii; base=4, pad=n)]...)
        pauli_jj = kron([paulis[kk + 1] for kk in digits(jj; base=4, pad=n)]...)
        mtx[ii + 1, jj + 1] = tr(mat(pauli_ii) * apply(qc, mat(pauli_jj))) / 2^n
    end
    return PauliTransferMatrix(n, mtx)
end

nonunital(ptm::PauliTransferMatrix{T}) where {T} = ptm.mtx[2:end, 1]

unital(ptm::PauliTransferMatrix{T}) where {T} = ptm.mtx[2:end, 2:end]

preserve_hermitian(ptm::PauliTransferMatrix) = all(isreal.(ptm.mtx))

function is_tp(ptm::PauliTransferMatrix{T}; atol=1e-10) where {T}
    return (
        isapprox(ptm.mtx[1, 1], one(Complex{T}); atol=atol) &&
        all(isapprox.(ptm.mtx[1, 2:end], zero(Complex{T}), atol=atol))
    )
end

function is_unital(ptm::PauliTransferMatrix{T}; atol=1e-10) where {T}
    return (
        isapprox(ptm.mtx[1, 1], one(Complex{T}); atol=atol) &&
        all(isapprox.(ptm.mtx[2:end, 1], zero(Complex{T}), atol=atol))
    )
end

function is_unitary(ptm::PauliTransferMatrix{T}; atol=1e-10) where {T}
    return (
        is_unital(ptm; atol=atol) && isapprox(
            unital(ptm) * transpose(unital(ptm)), I(4^ptm.n - 1); atol=atol
        )
    )
end

function is_cvx_sum_pauli(ptm::PauliTransferMatrix{T}; atol=1e-10) where {T}
    return all(isapprox.(ptm.mtx .- Diagonal(ptm.mtx), zero(Complex{T}), atol=atol))
end
