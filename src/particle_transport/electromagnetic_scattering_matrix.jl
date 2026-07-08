"""
    electromagnetic_scattering_matrix(𝓔::Vector{Float64},𝓑::Vector{Float64},charge::Real,
    Ω::Vector{Vector{Float64}},w::Vector{Float64},Ndims::Int64,Mn::Array{Float64},
    Dn::Array{Float64},pl::Vector{Int64},pm::Vector{Int64},P::Int64,Ng::Int64,
    Eb::Vector{Float64},ΔE::Vector{Float64},Qdims::Int64)

Compute the scattering matrix corresponding to the Lorentz force by external
electromagnetic fields.

# Input Argument(s)
- `𝓔::Vector{Float64}` : electric field along x-, y- and z-axis.
- `𝓑::Vector{Float64}` : magnetic field along x-, y- and z-axis.
- `charge::Real` : particle charge.
- `Ω::Vector{Vector{Float64}}` : angular quadrature points.
- `w::Vector{Float64}` : angulare quadrature weights.
- `Ndims::Int64` : dimension of the geometry.
- `Mn::Array{Float64}` : moment-to-discrete matrix.
- `Dn::Array{Float64}` : discrete-to-moment matrix.
- `pl::Vector{Int64}` : legendre order associated with each of the spherical harmonics in
  the interpolation basis.
- `pm::Vector{Int64}` : spherical harmonics order associated with each of the spherical
  harmonics in the interpolation basis.
- `P::Int64` : number of spherical harmonics in the interpolation basis.
- `Ng::Int64` : number of groups.
- `Eb::Vector{Float64}` : energy boundaries.
- `ΔE::Vector{Float64}` : energy widths.
- `Qdims::Int64` : dimension of the quadrature.

# Output Argument(s)
- `ℳ_EM::Array{Float64}` : scattering matrix (per group, moment-to-moment).
- `λ₀_EM::Vector{Float64}` : per-group extended-transport diagonal shift to add to the total
  cross section (zero for the `"spherical-harmonics"` method).

# Reference(s)
- Fan et al. (2013), Modeling Electron Transport in the Presence of Electric and Magnetic
  Fields.
- St-Aubin et al. (2015), A deterministic solution to the first order linear Boltzmann
  transport equation in the presence of external magnetic fields.

"""
function electromagnetic_scattering_matrix(𝓔::Vector{Float64},𝓑::Vector{Float64},charge::Real,Ω::Vector{Vector{Float64}},w::Vector{Float64},Ndims::Int64,Mn::Array{Float64},Dn::Array{Float64},pl::Vector{Int64},pm::Vector{Int64},P::Int64,Ng::Int64,Eb::Vector{Float64},ΔE::Vector{Float64},Qdims::Int64)

if Qdims != 3 error("Quadrature on the unit sphere is required for electromagnetic transport.") end

Nd = length(w)
ℳ_EM = zeros(Ng,P,P)
λ₀_EM = zeros(Ng)
c = 2.99792458
mₑc² = 0.510999
𝒯_EM = zeros(Ng)
for ig in range(1,Ng)
    𝒯_EM[ig] = 2/ΔE[ig] * log((sqrt(Eb[ig]+2*mₑc²)+sqrt(Eb[ig]))/(sqrt(Eb[ig+1]+2*mₑc²)+sqrt(Eb[ig+1])))
end

# Derivatives
μ = Ω[1]; η = Ω[2]; ξ = Ω[3];
M = zeros(Nd,P)
for n in range(1,Nd)

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
        l = pl[p]; m = pm[p]
        Clm = sqrt((2-(m == 0)) * exp( sum(log.(1:l-abs(m))) - sum(log.(1:l+abs(m))) ))
        if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
        if abs(μ[n]) != 1
            if l > 0
                if l > abs(m)
                    Plm = ferrer_associated_legendre(l,abs(m),μ[n])
                    Pl⁻m = ferrer_associated_legendre(l-1,abs(m),μ[n])
                    M[n,p] += ((l+abs(m))*Pl⁻m-l*μ[n]*Plm)/(1-μ[n]^2) * Clm * 𝓣m * (2*l+1)/(4*π) * charge*c* (𝓑[3]*η[n]-𝓑[2]*ξ[n])
                elseif l == m
                    Pll = ferrer_associated_legendre(l,l,μ[n])
                    M[n,p] += -l*μ[n]*Pll/(1-μ[n]^2) * Clm * 𝓣m * (2*l+1)/(4*π) * charge*c* (𝓑[3]*η[n]-𝓑[2]*ξ[n])
                end
            end
        else
            if abs(m) == 1
                M[n,p] += -l*(l+1)/2*(μ[n])^l * Clm * 𝓣m * (2*l+1)/(4*π) * charge*c* (𝓑[3]*cos(ϕ)-𝓑[2]*sin(ϕ))
            end
        end
    end

    # Derivative in ϕ
    for p in range(1,P)
        l = pl[p]; m = pm[p]
        Clm = sqrt((2-(m == 0)) * exp( sum(log.(1:l-abs(m))) - sum(log.(1:l+abs(m))) ))
        Plm = ferrer_associated_legendre(l,abs(m),μ[n])
        if abs(μ[n]) != 1
            if m ≥ 0
                M[n,p] += Plm * Clm * -m*sin(m*ϕ) * (2*l+1)/(4*π) * charge*c/(1-μ[n]^2) * (𝓑[2]*μ[n]*η[n]+𝓑[3]*μ[n]*ξ[n]-𝓑[1]*(1-μ[n]^2))
            else
                M[n,p] += Plm * Clm * abs(m)*cos(abs(m)*ϕ) * (2*l+1)/(4*π) * charge*c/(1-μ[n]^2) * (𝓑[2]*μ[n]*η[n]+𝓑[3]*μ[n]*ξ[n]-𝓑[1]*(1-μ[n]^2))
            end
        else
            if m == 1
                M[n,p] += l*(l+1)/2*(μ[n])^(l+1) * Clm * -m*sin(m*ϕ) * (2*l+1)/(4*π) * charge*c * (𝓑[2]*cos(ϕ)+𝓑[3]*sin(ϕ))
            elseif m == -1
                M[n,p] += l*(l+1)/2*(μ[n])^(l+1) * Clm * abs(m)*cos(abs(m)*ϕ) * (2*l+1)/(4*π) * charge*c * (𝓑[2]*cos(ϕ)+𝓑[3]*sin(ϕ))
            end
        end
    end

end

G = -(Dn * M)
for ig in range(1,Ng)
    ℳ_EM[ig,:,:] .= 𝒯_EM[ig] .* G
end

return ℳ_EM, λ₀_EM
end