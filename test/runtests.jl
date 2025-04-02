using ThermoDiff
using Test

@testset "ThermoDiff.jl" begin

@testset "eos.jl" begin
    # BENCHMARK - KYANITE
    ΔfH°        :: AbstractFloat = -2594220.46
    S°          :: AbstractFloat = 82.4300
    V°          :: AbstractFloat = 4.412
    CPbb84_k1   :: AbstractFloat = 262.68478
    CPbb84_k4   :: AbstractFloat = -2001.407
    CPbb84_k3   :: AbstractFloat = -1999740.000
    CPbb84_k8   :: AbstractFloat = -63181880.
    Vtwq_v1     :: AbstractFloat = 2.39725295
    Vtwq_v2     :: AbstractFloat = 0.00000000
    Vtwq_v3     :: AbstractFloat = -0.06459655
    Vtwq_v4     :: AbstractFloat = 0.00000000

    # benchmark against theriak at two temperatures and pressures
    P1         :: AbstractFloat = 12000.0
    T1         :: AbstractFloat = 573.15
    G = apparent_G(P1, T1, ΔfH°, S°, V°, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)
    @test G ≈ -2602600.11 atol=1e-2
    P2          :: AbstractFloat = 2000.0
    T2          :: AbstractFloat = 873.15
    G = apparent_G(P2, T2, ΔfH°, S°, V°, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)
    @test G ≈ -2713087.07 atol=1e-2

    # test with arrays
    P = repeat([P1, P2], 1, 2)
    T = repeat([T1, T2]', 2, 1)

    G = apparent_G(P, T, ΔfH°, S°, V°, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)
    @test G ≈ [-2602600.11 -2.668558377671368e6; -2.6468114979949193e6 -2713087.07] atol=1e-2
end

end
