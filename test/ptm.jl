using Test, CoherentNoise, Yao
using CoherentNoise: damping_channel, pauli_channel, unitary_channel
using CoherentNoise: PauliTransferMatrix, is_tp, preserve_hermitian, is_unital, is_unitary, is_cvx_sum_pauli


@testset "Basic Properties" begin
dc = damping_channel(0.001)

ρ = state(density_matrix(ghz_state(1)))

ptm = PauliTransferMatrix(dc)

display(ptm.mtx)

@test is_tp(ptm) 
@test !is_unital(ptm)
@test !is_unitary(ptm)

uc = unitary_channel(mat(Rx(0.1)))
ptm = PauliTransferMatrix(uc)

@test is_tp(ptm)
@test is_unital(ptm)
@test is_unitary(ptm)

pc= pauli_channel([1/4 for _ in 1:4])

ptm = PauliTransferMatrix(pc)
display(ptm.mtx)
@test is_tp(ptm)
@test is_cvx_sum_pauli(ptm)
@test is_unital(ptm)
@test !is_unitary(ptm)


end

