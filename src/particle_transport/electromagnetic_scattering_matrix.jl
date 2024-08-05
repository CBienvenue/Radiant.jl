function electromagnetic_scattering_matrix(𝓔::Vector{Float64},𝓑::Vector{Float64},q::Real,Ω,w,Ndims::Int64,Mn::Array{Float64},Dn::Array{Float64},pℓ::Vector{Int64},pm::Vector{Int64},P::Int64,Ng::Int64,Eb::Vector{Float64},ΔE::Vector{Float64},Qdims::Int64)

c = 1 # Speed of light (to set)
# case $μ = ±1$ to deal with
if Qdims != 3 error("Quadrature on the unit sphere is required for electromagnetic transport.") end

# Derivatives
μ = Ω[1]; η = Ω[2]; ξ = Ω[3];
Nd = length(w)
M = zeros(Nd,P)
@inbounds for n in range(1,Nd)

    # Compute ϕ
    if ξ[n] != 0
        if η[n] > 0
            ϕ = sign(ξ[n]) * atan(abs(ξ[n]/η[n]))
        elseif η[n] == 0
            ϕ = sign(ξ[n]) * π/2
        elseif η[n] < 0
            ϕ = sign(ξ[n]) * (π - atan(abs(ξ[n]/η[n])))
        end
    else
        if η[n] >= 0
            ϕ = 0.0
        elseif η[n] < 0
            ϕ = 1.0 * π
        end
    end

    # Derivative in μ
    for p in range(1,P)
        ℓ = pℓ[p]; m = pm[p]
        Cℓm = sqrt((2-(m == 0)) * exp( sum(log.(1:ℓ-abs(m))) - sum(log.(1:ℓ+abs(m))) ))
        if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
        if abs(μ[n]) != 1
            if ℓ > 0
                if ℓ > abs(m)
                    Pℓm = ferrer_associated_legendre(ℓ,abs(m),μ[n])
                    Pℓ⁻m = ferrer_associated_legendre(ℓ-1,abs(m),μ[n])
                    M[n,p] += ((ℓ+abs(m))*Pℓ⁻m-ℓ*μ[n]*Pℓm)/(1-μ[n]^2) * Cℓm * 𝓣m * (2*ℓ+1)/(4*π) * q*c* (𝓑[3]*η[n]-𝓑[2]*ξ[n])
                elseif ℓ == m
                    Pℓℓ = ferrer_associated_legendre(ℓ,ℓ,μ[n])
                    M[n,p] += -ℓ*μ[n]*Pℓℓ/(1-μ[n]^2) * Cℓm * 𝓣m * (2*ℓ+1)/(4*π) * q*c* (𝓑[3]*η[n]-𝓑[2]*ξ[n])
                end
            end
        else
            if abs(m) == 1
                M[n,p] += -ℓ*(ℓ+1)/2*(μ[n])^ℓ * Cℓm * 𝓣m * (2*ℓ+1)/(4*π) * q*c* (𝓑[3]*cos(ϕ)-𝓑[2]*sin(ϕ))
            end
        end
    end

    # Derivative in ϕ
    for p in range(1,P)
        ℓ = pℓ[p]; m = pm[p]
        Cℓm = sqrt((2-(m == 0)) * exp( sum(log.(1:ℓ-abs(m))) - sum(log.(1:ℓ+abs(m))) ))
        Pℓm = ferrer_associated_legendre(ℓ,abs(m),μ[n])
        if abs(μ[n]) != 1
            if m ≥ 0
                M[n,p] += Pℓm * Cℓm * -m*sin(m*ϕ) * (2*ℓ+1)/(4*π) * q*c/(1-μ[n]^2) * (𝓑[2]*μ[n]*η[n]+𝓑[3]*μ[n]*ξ[n]-𝓑[1]*(1-μ[n]^2))
            else
                M[n,p] += Pℓm * Cℓm * abs(m)*cos(abs(m)*ϕ) * (2*ℓ+1)/(4*π) * q*c/(1-μ[n]^2) * (𝓑[2]*μ[n]*η[n]+𝓑[3]*μ[n]*ξ[n]-𝓑[1]*(1-μ[n]^2))
            end
        else
            if m == 1
                M[n,p] += ℓ*(ℓ+1)/2*(μ[n])^(ℓ+1) * Cℓm * -m*sin(m*ϕ) * (2*ℓ+1)/(4*π) * q*c * (𝓑[2]*cos(ϕ)+𝓑[3]*sin(ϕ))
            elseif m == -1
                M[n,p] += ℓ*(ℓ+1)/2*(μ[n])^(ℓ+1) * Cℓm * abs(m)*cos(abs(m)*ϕ) * (2*ℓ+1)/(4*π) * q*c * (𝓑[2]*cos(ϕ)+𝓑[3]*sin(ϕ))
            end
        end
    end

end

# Energy-dependant coefficients
ℳ_EM = zeros(Ng,P,P)
mₑc² = 0.510999
for ig in range(1,Ng)
    𝒯_EM = 2/ΔE[ig] * log((sqrt(Eb[ig]+2*mₑc²)+sqrt(Eb[ig]))/(sqrt(Eb[ig+1]+2*mₑc²)+sqrt(Eb[ig+1])))
    for n in range(1,Nd), p in range(1,P), q in range(1,P)
        ℳ_EM[ig,q,p] += Dn[p,n] * M[n,q] * 𝒯_EM
    end
end

return ℳ_EM

end