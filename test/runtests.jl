using Radiant
using Test
import Radiant

include("tests.jl")

@testset "Radiant.jl" begin
    @testset "Mathematical tools" begin
        @testset "Quadrature" begin
            @testset "Gauss-Legendre quadrature" begin test_gauss_legendre() end
            @testset "Gauss-Lobatto quadrature" begin test_gauss_lobatto() end
            @testset "Carlson quadrature" begin test_carlson() end
            @testset "Lebedev" begin #= to do =# end
            @testset "Gauss-Chebychev" begin #= to do =# end
        end
        @testset "Polynomials" begin
            @testset "Legendre polynomials" begin test_legendre_polynomials() end
            @testset "Real spherical harmonics" begin test_real_spherical_harmonics() end
        end
        @testset "Integrals" begin #= to do =# end
        @testset "Iterative method" begin #= to do =# end
    end

    @testset "Cross sections" begin
        @testset "Energy group structures" begin test_energy_group_structure() end
    end
    
    @testset "Particle transport" begin
        
    end
end
nothing


