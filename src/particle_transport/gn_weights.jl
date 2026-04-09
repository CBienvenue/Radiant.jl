function gn_weights_spherical_harmonics(L::Int64,Nv::Int64,Ndims::Int64)
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if Nv <= 0 error("Number of direction cosine patches should be greater than zero.") end
    if ~(1 ≤ Ndims ≤ 3) error("Number of dimensions should be between 1 and 3.") end
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩 = zeros(Np,Np,Ndims,8,Nv,Nv)
    N = 32
    x,weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)

    sx = [1,1,1,1,-1,-1,-1,-1]
    sy = [1,1,-1,-1,1,1,-1,-1]
    sz = [1,-1,1,-1,1,-1,1,-1]

    # Offsets azimutaux par octant (ne dépend que de u)
    ϕ_offset_u = Vector{Float64}(undef, 8)
    for u in 1:8
        ϕ_offset_u[u] = (π / 2) * (2 + (sy[u] + 1) / 2 - (sz[u] + 1) / 2 - (sy[u] + 1) * (sz[u] + 1) / 2)
    end

    # Pré-calcul des valeurs 𝒯m(m,4*φp) avec 4*φp = π*(x+1)
    # et tables azimutales sur les m uniques
    m_unique = sort!(collect(Set(pm)))
    m_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(m_unique))
    m_idx_p = [m_to_idx[pm[p]] for p in 1:Np]

    Tm = Matrix{Float64}(undef, length(m_unique), N)
    for (im, m) in enumerate(m_unique)
        if m ≥ 0
            for n in 1:N
                Tm[im, n] = cos(m * (π * (x[n] + 1.0)))
            end
        else
            for n in 1:N
                Tm[im, n] = sin((-m) * (π * (x[n] + 1.0)))
            end
        end
    end
    TmW = Matrix{Float64}(undef, length(m_unique), N)
    for n in 1:N
        wn = weight[n]
        for im in 1:length(m_unique)
            TmW[im, n] = Tm[im, n] * wn
        end
    end
    az0_table = Matrix{Float64}(undef, length(m_unique), length(m_unique))
    LinearAlgebra.mul!(az0_table, TmW, LinearAlgebra.transpose(Tm))
    az0_pq = Matrix{Float64}(undef, Np, Np)
    for p in 1:Np
        imp = m_idx_p[p]
        for q in 1:Np
            imq = m_idx_p[q]
            az0_pq[p, q] = az0_table[imp, imq]
        end
    end

    # Pré-calcul des polynômes associés Pnm(l,|m|,±x)
    P_pos = Matrix{Float64}(undef, Np, N)
    P_neg = Matrix{Float64}(undef, Np, N)
    for p in 1:Np
        lp = pl[p]
        mp = pm[p]
        amp = abs(mp)
        for n in 1:N
            P_pos[p, n] = Pnm(lp, amp, x[n])
            P_neg[p, n] = Pnm(lp, amp, -x[n])
        end
    end

    # Factorisation des constantes Cp*Cq (en conservant la dépendance originale en mq==0)
    cp_base = Vector{Float64}(undef, Np)
    for p in 1:Np
        lp = pl[p]
        mp = pm[p]
        cp_base[p] = sqrt((2 * lp + 1) * factorial_factor([lp - abs(mp)], [lp + abs(mp)]))
    end
    cq_base = Vector{Float64}(undef, Np)
    mq_factor = Vector{Float64}(undef, Np)
    for q in 1:Np
        lq = pl[q]
        mq = pm[q]
        cq_base[q] = sqrt((2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
        mq_factor[q] = 2 * (2 - (mq == 0)) / π
    end
    sq = Vector{Float64}(undef, Np)
    for q in 1:Np
        sq[q] = mq_factor[q] * cq_base[q]
    end

    # Buffers réutilisés
    Pw = Matrix{Float64}(undef, Np, N)
    cosine_x = Matrix{Float64}(undef, Np, Np)
    cosine_yz = Matrix{Float64}(undef, Np, Np)

    cw = Vector{Float64}(undef, N)
    sw = Vector{Float64}(undef, N)
    cos_dt = Vector{Float64}(undef, N)
    sin_dt = Vector{Float64}(undef, N)

    TmC = Matrix{Float64}(undef, length(m_unique), N)
    az_tmp = Matrix{Float64}(undef, length(m_unique), length(m_unique))

    # Helper μ(v) : évite les closures μ_uv/ϕ_uw
    denom = Float64(Nv * (Nv + 1))
    mu_of_v(s::Int, v::Int) = (s == 1) ? (1.0 - ((Nv + 1 - v) * (Nv + 2 - v)) / denom) : (-1.0 + ((v - 1) * v) / denom)

    # Simplification: (Δμ*Δϕ)/4 * (π/(2ΔμΔϕ)) = π/8
    pref = π / 8

    for u in 1:8
        su = sx[u]
        Psign = (su == 1) ? P_pos : P_neg
        ϕ_offset = ϕ_offset_u[u]

        for v in 1:Nv
            Nw = (su == 1) ? (Nv + 1 - v) : v
            Δϕ = (π / 2) / Nw

            # Cos(Δϕ * t[n]) et Sin(Δϕ * t[n]) (utiles pour Ndims>=2)
            if Ndims ≥ 2
                for n in 1:N
                    θ = Δϕ * t[n]
                    cos_dt[n] = cos(θ)
                    sin_dt[n] = sin(θ)
                end
            end

            # Poids pour les intégrales cosinus (dépend de v, pas de w)
            μ0 = mu_of_v(su, v)
            μ1 = mu_of_v(su, v + 1)
            Δμ = μ1 - μ0

            # cosine_integral_x : poids = weight[n] * μi
            for p in 1:Np
                for n in 1:N
                    μi = μ0 + Δμ * t[n]
                    Pw[p, n] = Psign[p, n] * (weight[n] * μi)
                end
            end
            LinearAlgebra.mul!(cosine_x, Pw, LinearAlgebra.transpose(Psign))

            if Ndims ≥ 2
                # cosine_integral_yz : poids = weight[n] * sqrt(1-μi^2)
                for p in 1:Np
                    for n in 1:N
                        μi = μ0 + Δμ * t[n]
                        Pw[p, n] = Psign[p, n] * (weight[n] * sqrt(1 - μi^2))
                    end
                end
                LinearAlgebra.mul!(cosine_yz, Pw, LinearAlgebra.transpose(Psign))
            end

            # Ndims == 1 : tout est indépendant de w (par (u,v)) sauf la sortie elle-même
            if Ndims == 1
                for w in 1:Nw
                    for p in 1:Np
                        sp = pref * cp_base[p]
                        for q in 1:Np
                            𝒩[p, q, 1, u, v, w] += (sp * sq[q]) * az0_pq[p, q] * cosine_x[p, q]
                        end
                    end
                end
                continue
            end

            # Ndims >= 2 : partie azimutale dépend de w via cos(φi)/sin(φi)
            cosΔ = cos(Δϕ)
            sinΔ = sin(Δϕ)
            cosϕ0 = cos(ϕ_offset)
            sinϕ0 = sin(ϕ_offset)

            for w in 1:Nw
                # cw/sw = weight[n] * cos(φi_n) / sin(φi_n) avec φi_n = ϕ0 + Δϕ*t[n]
                for n in 1:N
                    c = cosϕ0 * cos_dt[n] - sinϕ0 * sin_dt[n]
                    s = sinϕ0 * cos_dt[n] + cosϕ0 * sin_dt[n]
                    cw[n] = weight[n] * c
                    sw[n] = weight[n] * s
                end

                # azimutal_integral_x est constant (az0_pq), azimutal_integral_y/z via cw/sw
                for n in 1:N
                    cwn = cw[n]
                    swn = sw[n]
                    for im in 1:length(m_unique)
                        TmC[im, n] = Tm[im, n] * cwn
                    end
                end
                LinearAlgebra.mul!(az_tmp, TmC, LinearAlgebra.transpose(Tm))

                # Remplissage dim 1 & 2
                for p in 1:Np
                    sp = pref * cp_base[p]
                    imp = m_idx_p[p]
                    for q in 1:Np
                        imq = m_idx_p[q]
                        common = (sp * sq[q])
                        𝒩[p, q, 1, u, v, w] += common * az0_pq[p, q] * cosine_x[p, q]
                        𝒩[p, q, 2, u, v, w] += common * az_tmp[imp, imq] * cosine_yz[p, q]
                    end
                end

                if Ndims == 3
                    for n in 1:N
                        swn = sw[n]
                        for im in 1:length(m_unique)
                            TmC[im, n] = Tm[im, n] * swn
                        end
                    end
                    LinearAlgebra.mul!(az_tmp, TmC, LinearAlgebra.transpose(Tm))

                    for p in 1:Np
                        sp = pref * cp_base[p]
                        imp = m_idx_p[p]
                        for q in 1:Np
                            imq = m_idx_p[q]
                            𝒩[p, q, 3, u, v, w] += (sp * sq[q]) * az_tmp[imp, imq] * cosine_yz[p, q]
                        end
                    end
                end

                # ϕ0 += Δϕ (update trig via récurrence)
                cosϕ1 = cosϕ0 * cosΔ - sinϕ0 * sinΔ
                sinϕ1 = sinϕ0 * cosΔ + cosϕ0 * sinΔ
                cosϕ0 = cosϕ1
                sinϕ0 = sinϕ1
            end
        end
    end
    return 𝒩
end

# function gn_weights_spherical_harmonics(L::Int64,Nv::Int64,Ndims::Int64)
#     if L < 0 error("Legendre order is greater or equal to zero.") end
#     if Nv <= 0 error("Number of direction cosine patches should be greater than zero.") end
#     if ~(1 ≤ Ndims ≤ 3) error("Number of dimensions should be between 1 and 3.") end
#     Np = spherical_harmonics_number_basis(L)
#     pl,pm = spherical_harmonics_indices(L)
#     𝒩 = zeros(Np,Np,Ndims,8,Nv,Nv)
#     N = 32
#     x,weight = gauss_legendre(N)
#     sx = [1,1,1,1,-1,-1,-1,-1]
#     sy = [1,1,-1,-1,1,1,-1,-1]
#     sz = [1,-1,1,-1,1,-1,1,-1]
#     μ_uv(u,v) = sx[u]*(1-(1-v+(sx[u]+1)/2*Nv)*(-v+(sx[u]+1)/2*(Nv+2))/(Nv*(Nv+1)))
#     ϕ_uw(u,v,w) = (π/2)*(w-1)/(-sx[u]*v + (sx[u]+1)/2*(Nv+1)) + π/2 * (2 + (sy[u]+1)/2 - (sz[u]+1)/2 - (sy[u]+1)*(sz[u]+1)/2)
#     for p in range(1,Np), q in range(1,Np)
#         lp = pl[p]
#         lq = pl[q]
#         mp = pm[p]
#         mq = pm[q]
#         Cp = sqrt(2*(2-(mq==0))/π * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
#         Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
#         for u in range(1,8)
#             for v in range(1,Nv)
#                 Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
#                 for w in range(1,Nw)
#                     Δμ_uv = μ_uv(u,v+1) - μ_uv(u,v)
#                     Δϕ_uw = ϕ_uw(u,v,w+1) - ϕ_uw(u,v,w)

#                     # Azimuthal integral
#                     azimutal_integral_x = 0.0
#                     azimutal_integral_y = 0.0
#                     azimutal_integral_z = 0.0
#                     for n in range(1,N)
#                         φi = ϕ_uw(u,v,w) + 0.5 * Δϕ_uw * (x[n]+1)
#                         φp = π/2 * (φi-ϕ_uw(u,v,w))/Δϕ_uw
#                         azimutal_integral_x += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp)
#                         if (Ndims ≥ 2) azimutal_integral_y += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp) * cos(φi) end
#                         if (Ndims == 3) azimutal_integral_z += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp) * sin(φi) end
#                     end

#                     # Cosine integral
#                     cosine_integral_x = 0.0
#                     cosine_integral_yz = 0.0
#                     for n in range(1,N)
#                         μi = μ_uv(u,v) + 0.5 * Δμ_uv * (x[n]+1)
#                         μp = sx[u] * ((sx[u]-1)/2+(μi-μ_uv(u,v))/Δμ_uv)
#                         cosine_integral_x += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μp-1) * μi
#                         if (Ndims ≥ 2) cosine_integral_yz += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μp-1) * sqrt(1-μi^2) end
#                     end

#                     # Factor
#                     Cuvw = π/(2*Δμ_uv*Δϕ_uw)

#                     # Accumulate weights
#                     𝒩[p,q,1,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_x * cosine_integral_x
#                     if (Ndims ≥ 2) 𝒩[p,q,2,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_y * cosine_integral_yz end
#                     if (Ndims == 3) 𝒩[p,q,3,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_z * cosine_integral_yz end

#                 end
#             end
#         end
#     end
#     return 𝒩
# end

