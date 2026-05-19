"""
    sommerfield(Ei::Float64,L::Int64)

Gives the Legendre moments of the Sommerfield angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wl::Vector{Wl}` : Sommerfield angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon

"""
function sommerfield(Ei::Float64,L::Int64)
    β = sqrt(Ei*(Ei+2)/(Ei+1)^2)
    Wl = zeros(L+1)
    for l in range(0,L)
        if l == 0
            Wl[l+1] = 1
        elseif l == 1
            Wl[l+1] = (2*β + (1-β^2)*log((1-β)/(1+β)))/(2*β^2)
        else
            Wl[l+1] = ((2*l-1)*Wl[l] - l*β*Wl[l-1])/((l-1)*β)
        end
    end
    return Wl
end