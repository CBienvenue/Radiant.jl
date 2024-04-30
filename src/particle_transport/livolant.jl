"""
    livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

Estimate the solution using estimate from Livolant acceleration method. 

# Input Argument(s)
- '𝚽₂::Array{Float64}': flux at iteration (i)
- '𝚽₁::Array{Float64}': flux at iteration (i-1)
- '𝚽₀::Array{Float64}': flux at iteration (i-2)

# Output Argument(s)
- '𝚽::Array{Float64}': flux estimated by Livolant acceleration method for iteration (i+1).

# Reference(s)
- Hébert (2016) : Applied Reactor Physics (Sect. C.1.3 - Iterative approach).

"""
function livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

# Errors vectors
e₀ = vec(𝚽₁) - vec(𝚽₀)
e₁ = vec(𝚽₂) -  vec(𝚽₁)
Δe = e₁-e₀

# Acceleration factor
μj = -sum(e₀.*Δe)/(sum(Δe.^2))
if μj ≤ 0 μj = 1 end

# Compute the next iterations
𝚽 = reshape(μj * 𝚽₂ + (1-μj) * 𝚽₁,size(𝚽₂))

return 𝚽

end