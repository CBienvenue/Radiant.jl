"""
    livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

Estimate the solution using estimate from Livolant acceleration method. 

# Input Argument(s)
- `𝚽₂::Array{Float64}`: flux at iteration (i)
- `𝚽₁::Array{Float64}`: flux at iteration (i-1)
- `𝚽₀::Array{Float64}`: flux at iteration (i-2)

# Output Argument(s)
- `𝚽::Array{Float64}`: flux estimated by Livolant acceleration method for iteration (i+1).

# Reference(s)
- Hébert (2016),  Applied Reactor Physics (Sect. C.1.3 - Iterative approach).

"""
function livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

# Compute acceleration factor
sum_eΔe = 0.0
sum_Δe2 = 0.0
for i in eachindex(𝚽₂, 𝚽₁, 𝚽₀)
    e₀ = 𝚽₁[i] - 𝚽₀[i]
    e₁ = 𝚽₂[i] - 𝚽₁[i]
    Δe = e₁ - e₀
    sum_eΔe += e₀ * Δe
    sum_Δe2 += Δe * Δe
end
μj = -sum_eΔe / sum_Δe2
if μj ≤ 0 μj = 1 end

# Compute the next iteration
𝚽 = similar(𝚽₂)
for i in eachindex(𝚽₂, 𝚽₁)
    𝚽[i] = μj * 𝚽₂[i] + (1-μj) * 𝚽₁[i]
end

return 𝚽

end