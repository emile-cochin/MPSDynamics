using MPSDynamics
using Test

coverage_test = ["runsim_tests.jl"]
#physics_test = ["chain_coeffs_tests.jl",
physics_test = ["pure_dephasing_tests.jl"]

@testset "MPSDynamics tests" begin
    
    @testset "coverage tests" begin
        for test in coverage_test
            include("coverage/" * test)
        end
    end

    @testset "physics tests" begin
        for test in physics_test
            include("physics/" * test)
        end
    end
end
