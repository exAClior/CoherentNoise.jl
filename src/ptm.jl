using Random, LinearAlgebra 
using Yao

struct QuantumChannel{T<:AbstractFloat}
    n::Int
    kraus_ops::Vector{Matrix{Complex{T}}}
end

damping_channel(p::T) where {T<:AbstractFloat}  = QuantumChannel(1, [Complex{T}[1.0 0.0; 0.0 sqrt(1.0-p)], Complex{T}[0.0 sqrt(p); 0.0 0.0]])

function apply(qc::QuantumChannel{T}, ρ::MT) where {T<:AbstractFloat, MT <: AbstractMatrix{Complex{T}}}
    if size(ρ) != (2^qc.n, 2^qc.n)
        throw(DimensionMismatch("Dimension mismatch"))
    end
    ρ_new = mapreduce(x -> x * ρ * x', .+,  qc.kraus_ops) 
    return ρ_new
end

dc = damping_channel(0.001)

ρ = state(density_matrix(ghz_state(1)))

struct PauliTransferMatrix{T<:AbstractFloat}
    n::Int
    mtx::Matrix{Complex{T}}
end

function PauliTransferMatrix(qc::QuantumChannel{T}) where {T}
    n = qc.n
    paulis = [I2, X, Y, Z]
    mtx = zeros(Complex{T}, 4^n, 4^n)
    for ii in 0:4^n-1, jj in 0:4^n-1
        pauli_ii = kron([paulis[kk+1] for kk in digits(ii, base=4, pad=n)]...) 
        pauli_jj = kron([paulis[kk+1] for kk in digits(jj, base=4, pad=n)]...) 
        mtx[ii+1,jj+1] = tr(mat(pauli_ii) * apply(qc, mat(pauli_jj))) / 2^n
    end
    return PauliTransferMatrix(n, mtx)
end

ptm = PauliTransferMatrix(dc)

display(ptm.mtx)