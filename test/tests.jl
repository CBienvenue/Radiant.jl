function test_legendre_polynomials()
    n = 4; x = 1/2;
    P = Radiant.legendre_polynomials(n,x)
    @test begin
        P[1] ≈ 1                   &&
        P[2] ≈ x                   &&
        P[3] ≈ (3*x^2-1)/2         &&
        P[4] ≈ (5*x^3-3*x)/2       &&
        P[5] ≈ (35*x^4-30*x^2+3)/8
    end
end

function test_real_spherical_harmonics()

    μ = 1/2; η = -1/3; ξ = sqrt(1-μ^2-η^2)
    ϕ = sign(ξ) * (π - atan(abs(ξ/η)))

    ℓ = 0; m = 0;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,1,atol=1e-6)

    ℓ = 1; m = -1;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,ξ,atol=1e-6)

    ℓ = 1; m = 0;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,μ,atol=1e-6)

    ℓ = 2; m = -2;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,sqrt(3)*η*ξ,atol=1e-6)

    ℓ = 1; m = 1;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,η,atol=1e-6)

    ℓ = 2; m = -1;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,sqrt(3)*μ*ξ,atol=1e-6)

    ℓ = 2; m = 0;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,(3*μ^2-1)/2,atol=1e-6)

    ℓ = 2; m = 1;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,sqrt(3)*μ*η,atol=1e-6)

    ℓ = 2; m = 2;
    Rℓm = Radiant.real_spherical_harmonics(ℓ,m,μ,ϕ)
    @test isapprox(Rℓm,sqrt(3)/2*(η^2-ξ^2),atol=1e-6)

    return true
end

function test_gauss_legendre()

    μ,w = Radiant.gauss_legendre(1)
    @test abs(μ[1]) < eps() && w[1] ≈ 2
    
    μ,w = Radiant.gauss_legendre(2)
    @test begin 
        μ[1] ≈ 1/sqrt(3) && w[1] ≈ 1    &&
        μ[2] ≈ -μ[1]     && w[2] ≈ w[1]
    end

    μ,w = Radiant.gauss_legendre(3)
    @test begin 
        μ[1] ≈ sqrt(3/5)  && w[1] ≈ 5/9  &&
        abs(μ[2]) < eps() && w[2] ≈ 8/9  &&
        μ[3] ≈ -μ[1]      && w[3] ≈ w[1]
    end

    μ,w = Radiant.gauss_legendre(4)
    @test begin 
        μ[1] ≈ sqrt(3/7+2/7*sqrt(6/5)) && w[1] ≈ (18-sqrt(30))/36 &&
        μ[2] ≈ sqrt(3/7-2/7*sqrt(6/5)) && w[2] ≈ (18+sqrt(30))/36 &&
        μ[3] ≈ -μ[2]                   && w[3] ≈ w[2]             &&
        μ[4] ≈ -μ[1]                   && w[4] ≈ w[1]
    end

    μ,w = Radiant.gauss_legendre(64)
    @test begin
        μ[1] ≈ 0.9993050417357722 && w[1] ≈ 0.0017832807216964 &&
        μ[2] ≈ 0.9963401167719553 && w[2] ≈ 0.0041470332605625
    end

end

function test_gauss_lobatto()

    μ,w = Radiant.gauss_lobatto(3)
    @test begin
        μ[1] ≈ 1          && w[1] ≈ 1/3  &&
        abs(μ[2]) < eps() && w[2] ≈ 4/3  &&
        μ[3] ≈ -μ[1]      && w[3] ≈ w[1]
    end

    μ,w = Radiant.gauss_lobatto(6)
    @test begin
        μ[1] ≈ 1 && w[1] ≈ 1/15 &&
        μ[2] ≈ sqrt(1/3+2*sqrt(7)/21) && w[2] ≈ (14-sqrt(7))/30 &&
        μ[3] ≈ sqrt(1/3-2*sqrt(7)/21) && w[3] ≈ (14+sqrt(7))/30 &&
        μ[4] ≈ -μ[3] && w[4] ≈ w[3] &&
        μ[5] ≈ -μ[2] && w[5] ≈ w[2] &&
        μ[6] ≈ -μ[1] && w[6] ≈ w[1]
    end
    
end

function test_carlson()

    Ω,w = Radiant.carlson(2)
    @test begin
        length(w) == 8 &&
        isapprox(Ω[1][1],0.5773503,atol=1e-6)
    end

    Ω,w = Radiant.carlson(10)
    @test begin
        length(w) == 120 &&
        isapprox(Ω[1][1],0.1666667,atol=1e-6)  &&
        isapprox(Ω[1][6],0.4699819,atol=1e-6)  &&
        isapprox(Ω[1][10],0.7059675,atol=1e-6) &&
        isapprox(Ω[1][13],0.8746233,atol=1e-6) &&
        isapprox(Ω[1][15],0.9759494,atol=1e-6)
    end

end

function test_energy_group_structure()

    E = 10; Ec = 0.001; Ng = 10;

    Eᵇ = Radiant.energy_group_structure(Ng,E,Ec,"log")
    @test begin
        isapprox((Eᵇ[1] + Eᵇ[2])/2,E,atol=1e-6)      &&
        isapprox(Eᵇ[end],Ec,atol=1e-6)               &&
        isapprox(Eᵇ[1],14.453910961767651,atol=1e-6)
    end

    Eᵇ = Radiant.energy_group_structure(Ng,E,Ec,"linear")
    @test begin
        isapprox((Eᵇ[1] + Eᵇ[2])/2,E,atol=1e-6)      &&
        isapprox(Eᵇ[end],Ec,atol=1e-6)               &&
        isapprox(Eᵇ[1],10.526263157894737,atol=1e-6)
    end

end