"""
    map_moments(𝒪i::Vector{Int64},𝒪f::Vector{Int64})

Produce a map between different space/energy moments configuration of the angular flux.

# Input Argument(s)
- `𝒪i::Vector{Int64}`: incoming moment orders in space and energy.
- `𝒪f::Vector{Int64}`: outgoing moment orders in space and energy.

# Output Argument(s)
- `map::Vector{Int64}`: vector of correspondance between the incoming and outgoing moments.

# Reference(s)
N/A

"""
function map_moments(𝒪i::Vector{Int64},𝒪f::Vector{Int64})

Nmi = prod(𝒪i)
Nmf = prod(𝒪f)

M = zeros(max(𝒪i[1],𝒪f[1]),max(𝒪i[2],𝒪f[2]),max(𝒪i[3],𝒪f[3]),max(𝒪i[4],𝒪f[4]))
map = zeros(Int64,Nmf)

for ix in range(1,𝒪i[1]), iy in range(1,𝒪i[2]), iz in range(1,𝒪i[3]), iE in range(1,𝒪i[4])
    i = 𝒪i[2]*𝒪i[1]*𝒪i[4]*(iz-1)+𝒪i[1]*𝒪i[4]*(iy-1)+𝒪i[4]*(ix-1)+iE
    M[ix,iy,iz,iE] = i
end

for ix in range(1,𝒪f[1]), iy in range(1,𝒪f[2]), iz in range(1,𝒪f[3]), iE in range(1,𝒪f[4][1])
    i = 𝒪f[2]*𝒪f[1]*𝒪f[4]*(iz-1)+𝒪f[1]*𝒪f[4]*(iy-1)+𝒪f[4]*(ix-1)+iE
    if ix <= 𝒪i[1] && iy <= 𝒪i[2] && iz <= 𝒪i[3] && iE <= 𝒪i[4]
        map[i] = M[ix,iy,iz,iE]
    end
end

return map

end