"""
    quadrature(N::Int64)

Computing the weights and evaluation points of requested quadrature.

# Input Argument(s)
- 'N::Int64': quadrature order.
- 'type::String': type of quadrature.
- 'Ndims::Int64=1': geometry dimension.

# Output Argument(s)
- 'Ω::Vector{Vector{Float64}}': evaluation points.
- 'w::Vector{Float64}': weights.

# Reference(s)
N/A

"""
function quadrature(N::Int64,type::String,Ndims::Int64=1,Qdims::Int64=Ndims)

if type == "gauss-legendre" && Ndims == 1
    Ω,w = gauss_legendre(N)
elseif type == "gauss-lobatto" && Ndims == 1
    Ω,w = gauss_lobatto(N)
elseif type == "gauss-legendre-chebychev"
    Ω,w = gauss_legendre_chebychev(N)
elseif type == "lebedev"
    Ω,w = lebedev(N)
elseif type == "carlson"
    Ω,w = carlson(N)
else
    error("Error in quadrature.jl: Unknown $type quadrature in $(Ndims)D geometry.")
end

# Adapt the 3D unit sphere quadrature for 2D calculations (filter out the values with ξ < 0)
if Qdims == 2
    condition = Ω[3] .≥ 0
    μ₀ = Ω[1][condition]
    η₀ = Ω[2][condition]
    ξ₀ = Ω[3][condition]
    w₀ = copy(w)
    Nd = length(w)
    for i in range(1,Nd)
        if ~condition[i]
            ϵ_ξ = zeros(Nd)
            for j in range(1,Nd)
                ϵ_ξ[j] = abs(Ω[1][i] - Ω[1][j])^2 + abs(Ω[2][i] - Ω[2][j])^2 + abs(Ω[3][i] + Ω[3][j])^2
            end
            w₀[argmin(ϵ_ξ)] += w₀[i]
        end
    end
    w = w₀[condition]
    Ω = [μ₀,η₀,ξ₀]
end

# Output values
return Ω,w

end