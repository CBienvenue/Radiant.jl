function pn_weights_legendre(L::Int64)
    Np = L+1
    pl = collect(0:L)
    𝒩⁻ = zeros(Np)
    𝒩 = zeros(Np)
    𝒩⁺ = zeros(Np)
    for p in range(1,Np)
        l = pl[p]
        if (l > 0) 𝒩⁻[p] = 0.5 * l/(sqrt(2l-1)*sqrt(2l+1)) end
        𝒩[p] = 0.5
        if (l < L) 𝒩⁺[p] = 0.5 * (l+1)/(sqrt(2l+1)*sqrt(2l+3)) end
    end
    return 𝒩⁻,𝒩,𝒩⁺
end

function pn_weights_spherical_harmonics(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩⁻ = zeros(Np)
    𝒩 = zeros(Np)
    𝒩⁺ = zeros(Np)
    for p in range(1,Np)
        l = pl[p]
        m = pm[p]
        if (l > 0 && abs(m) <= l-1) 𝒩⁻[p] = 0.5 * sqrt((l-m)*(l+m))/(sqrt(2l-1)*sqrt(2l+1)) end
        𝒩[p] = 0.5
        if (l < L) 𝒩⁺[p] = 0.5 * sqrt((l-m+1)*(l+m+1))/(sqrt(2l+1)*sqrt(2l+3)) end
    end
    return 𝒩⁻,𝒩,𝒩⁺
end

function pn_weights_cartesian_harmonics(L::Int64)
    Np = cartesian_harmonics_number_basis(L)
    pl,pa,pb,pc = cartesian_harmonics_indices(L)
    𝒩⁻ = zeros(Np)
    𝒩 = zeros(Np)
    𝒩⁺ = zeros(Np)
    for p in range(1,Np)
        l = pl[p]
        a = pa[p]
        b = pb[p]
        c = pc[p]
        if (l > 0 && a > 0) 𝒩⁻[p] = 0.5 * a/(2*l+1) * cartesian_harmonics_normalization(l,a,b,c)/cartesian_harmonics_normalization(l-1,a-1,b,c) end
        𝒩[p] = 0.5
        if (l < L) 𝒩⁺[p] = 0.5 * cartesian_harmonics_normalization(l,a,b,c)/cartesian_harmonics_normalization(l+1,a+1,b,c) end
    end
    return 𝒩⁻,𝒩,𝒩⁺
end