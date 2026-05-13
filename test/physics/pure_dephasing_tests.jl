import MPSDynamics: chaincoeffs_ohmic, puredephasingmpo
import QuadGK: quadgk

@testset "pure dephasing tests" begin
    d = 10 # number of Fock states of the chain modes
    N = 30 # length of the chain
    α = 0.01 # coupling strength
    ω0 = 0.008 # TLS gap
    s = 1 # ohmicity
    ωc = 0.035 # Cut-off of the spectral density J(ω)

    dt = 1.0 # time step
    tfinal = 10.0 # simulation time
    method = :TDVP1 # time-evolution method

    D = 6 # MPS bond dimension

    # Initial electronic system in a superposition of 1 and 2
    ψ = zeros(2)
    ψ[1] = 1/sqrt(2)
    ψ[2] = 1/sqrt(2)

    Johmic(ω,s) = (2*α*ω^s)/(ωc^(s-1))
    time_analytical = LinRange(0.0, tfinal, Int(tfinal))

    ob1 = OneSiteObservable("sz", sz, 1)

    @testset "0 temperature" begin
        β = Inf
        cpars = chaincoeffs_ohmic(N, α, s; ωc=ωc)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
        H = puredephasingmpo(ω0, d, N, cpars)
        A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
        A, dat = runsim(dt, tfinal, A, H;
                        name = "pure dephasing model with temperature",
                        method = method,
                        obs = [ob1],
                        convobs = [ob1],
                        params = @LogParams(ω0, N, d, α, s),
                        convparams = D,
                        reduceddensity=true,
                        verbose = false,
                        save = false,
                        plot = false,
                        );
        Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))/x^2, 0, ωc)[1]
        Decoherence_ohmic(t) = 0.5*exp(Γohmic(t))
        t, ρ12 = dat["data/times"], abs.(dat["data/Reduced ρ"][1,2,:])
        @test isapprox(ρ12, Decoherence_ohmic.(t); rtol=1e-6)
    end

    @testset "finite temperature" begin
        β = 100
        cpars = chaincoeffs_finiteT(N, β; α=α, s=s, J=nothing, ωc=ωc, mc=4, mp=0, AB=nothing, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false, smooth=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
        H = puredephasingmpo(ω0, d, N, cpars)
        A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>
        A, dat = runsim(dt, tfinal, A, H;
                        name = "pure dephasing model with temperature",
                        method = method,
                        obs = [ob1],
                        convobs = [ob1],
                        params = @LogParams(ω0, N, d, α, s),
                        convparams = D,
                        reduceddensity=true,
                        verbose = false,
                        save = false,
                        plot = false,
                        );
        Γohmic(t) = - quadgk(x -> Johmic(x,s)*(1 - cos(x*t))*coth(β*x/2)/x^2, 0, ωc)[1]
        Decoherence_ohmic(t) = 0.5*exp(Γohmic(t))
        t, ρ12 = dat["data/times"], abs.(dat["data/Reduced ρ"][1,2,:])
        @test isapprox(ρ12, Decoherence_ohmic.(t); rtol=1e-6)
    end
end
