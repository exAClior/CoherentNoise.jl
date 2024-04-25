using Test, CoherentNoise, Yao
using CoherentNoise: damping_channel, pauli_channel, unitary_channel, depolarizing_channel
using CoherentNoise: QuantumChannel, PauliTransferMatrix, is_tp, preserve_hermitian, is_unital, is_unitary, is_cvx_sum_pauli, average_infidelity
using Plots

@testset "Rotation/dephasing channel" begin
    θ = 0.001*π 
    qD = 0.00001
    qR = 0.5
    ϵ = 2 * qD - qR * (1 - cos(θ)) 
    δ = qR * sin(θ)
    # when ϵ is on the order of δ^2 we will see quadratic accumulation of error
    @show ϵ, δ^2
    m = 100
    mixed_channel = QuantumChannel(1,sqrt.([1.0-qD-qR,qD,qR]).* mat.([I2,X,Rx(θ)]))
    ptm = PauliTransferMatrix(mixed_channel)

    infids = [average_infidelity(ptm^ii) for ii in 1:m]
    # for this insane values, we see quadratic accumulation of error
    # θ = 0.001*π 
    # qD = 0.00001
    # qR = 0.5
    plot(1:m,infids)


end

@testset "Accumulation speed" begin
    m = 100
    uc = unitary_channel(mat(Rx(0.0001)))
    ptm = PauliTransferMatrix(uc)
    infids = [average_infidelity(ptm^ii) for ii in 1:m]
    # it looks quadratic
    plot(1:m,infids)

    dpc = depolarizing_channel(0.001)
    ptm = PauliTransferMatrix(dpc)
    infids = [average_infidelity(ptm^ii) for ii in 1:m]
    # linear 
    plot(1:m, infids)
end


@testset "Basic Properties" begin
dc = damping_channel(0.001)

ρ = state(density_matrix(ghz_state(1)))

ptm = PauliTransferMatrix(dc)

display(ptm.mtx)

@test is_tp(ptm) 
@test !is_unital(ptm)
@test !is_unitary(ptm)

r = average_infidelity(ptm)

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

dpc = depolarizing_channel(0.1)


end

