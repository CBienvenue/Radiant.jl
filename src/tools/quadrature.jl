"""
    quadrature(N::Int64)

Computing the weights and evaluation points of requested quadrature.

See also [`compute_flux`](@ref Radiant.compute_flux), [`gauss_legendre`](@ref),
[`gauss_lobatto`](@ref),[`gauss_legendre_chebychev`](@ref), [`lebedev`](@ref),
[`carlson`](@ref).

# Input Argument(s)
- 'N::Int64': quadrature order.
- 'type::String': type of quadrature.
- 'Ndims::Int64=1': geometry dimension.

# Output Argument(s)
- 'Ω::Vector{Vector{Float64}}': evaluation points.
- 'w::Vector{Float64}': weights.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function quadrature(N::Int64,type::String,Ndims::Int64=1)

# Quadrature over the [-1,1] domain of integration
if Ndims == 1

    if type == "gauss-legendre"
        Ω,w = gauss_legendre(N)
    elseif type == "gauss-lobatto"
        Ω,w = gauss_lobatto(N)
    else
        error("Error in quadrature.jl: Unknown ",type," quadrature in 1D geometry.")
    end

# Quadrature over the unit sphere
else

    if type == "gauss-legendre-chebychev"
        Ω,w = gauss_legendre_chebychev(N)
    elseif type == "lebedev"
        Ω,w = lebedev(N)
    elseif type == "carlson"
        Ω,w = carlson(N)
    else
        error("Error in quadrature.jl: Unknown ",type," quadrature in ",Ndims,"D geometry.")
    end

    # Adapt the 3D unit sphere quadrature for 2D calculations
    if Ndims == 2

        # Exclude the values with ξ < 0
        Ndir = length(w)
        μ₀ = zeros(0); η₀ = zeros(0); ξ₀ = zeros(0); w₀ = zeros(0)
        μ₀_left = zeros(0); η₀_left = zeros(0); ξ₀_left = zeros(0); w₀_left = zeros(0)
        for n in range(1,Ndir)
            if Ω[3][n] ≥ 0
                push!(μ₀,Ω[1][n]); push!(η₀,Ω[2][n]); push!(ξ₀,Ω[3][n]); push!(w₀,w[n])
            else
                push!(μ₀_left,Ω[1][n]); push!(η₀_left,Ω[2][n]); push!(ξ₀_left,Ω[3][n]); push!(w₀_left,w[n])
            end
        end

        # Adjust the weights
        for n in range(1,length(ξ₀_left))
            index = intersect( findall(x -> x ≈ μ₀_left[n],μ₀), findall(x -> x ≈ η₀_left[n],η₀), findall(x -> x ≈ -ξ₀_left[n],ξ₀) )
            w₀[index[1]] += w₀_left[n]
        end
        
        # Reactualize the points and weights
        Ω = [μ₀,η₀,ξ₀]
        w = w₀

    end

end

# Output values
return Ω,w

end