function half_to_full_range_matrix_legendre(L::Int64)
    Np = L+1
    pl = collect(0:L)
    Mll = zeros(Np,Np)
    for ip in range(1,Np), jp in range(1,Np)
        il = pl[ip]
        jl = pl[jp]
        for ik in range(0,div(il,2)), jk in range(0,div(jl,2))
            for j in range(0,jl-2*jk)
                Mll[ip,jp] += sqrt(2*jl+1)/2^(il+jl) * (-1)^(ik+jk) * binomial(il,ik) * binomial(jl,jk) * binomial(2*il-2*ik,il) * binomial(2*jl-2*jk,jl) * binomial(jl-2*jk,j) * (-1)^(jl-2*jk-j) * 2^j / (il-2*ik+j+1)
            end
        end
    end
    return Np,Mll
end

function half_to_full_range_matrix_spherical_harmonics(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    Mll = zeros(Np,Np)
    for ip in range(1,Np), jp in range(1,Np)
        il = pl[ip]
        jl = pl[jp]
        im = pm[ip]
        jm = pm[jp]
        if im == jm
            C = π * (1+(im==0)) * sqrt((2-(im==0))*factorial_factor([il-abs(im)],[il+abs(im)])) * sqrt((2-(jm==0))/(2*π) * (2*jl+1)*factorial_factor([jl-abs(jm)],[jl+abs(jm)]))
            C2 = 1/(2^(abs(im))) * factorial_factor([1],[il-abs(im),jl-abs(im)])
            for ik in range(0,il-abs(im)), jk in range(0,jl-abs(im))
                C3 = (-1)^(ik+jk) * binomial(il-abs(im),ik) * binomial(jl-abs(im),jk) * factorial_factor([il+abs(im)+ik,jl+abs(im)+jk],[abs(im)+ik,abs(im)+jk])/2^(ik)
                if iseven(abs(im))
                    for i in range(0,abs(im)+ik+jk), j in range(0,div(abs(im),2))
                        C4 = (-1)^i * binomial(abs(im)+ik+jk,i) * binomial(div(abs(im),2),j)
                        Mll[ip,jp] += C * C2 * C3 * C4 / (i + j + div(abs(im),2) + 1)
                    end
                else
                    for i in range(0,abs(im)+ik+jk), j in range(0,div(abs(im)-1,2))
                        C4 = (-1)^i * binomial(abs(im)+ik+jk,i) * binomial(div(abs(im)-1,2),j)
                        Mll[ip,jp] += C * C2 * C3 * C4 * (𝒢₈(i + j + div(abs(im)-1,2),0,1,1,1) - 𝒢₈(i + j + div(abs(im)-1,2),0,1,1,0))
                    end
                end
            end
        end
    end
    return Np,Mll
end

function quarter_to_full_range_matrix_spherical_harmonics(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    Mll = zeros(Np,Np)
    N = 32
    x,w = gauss_legendre(N)
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        lq = pl[q]
        mp = pm[p]
        mq = pm[q]
        Cp = sqrt((2-(mp==0)) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
        Cq = (-1)^abs(mq) * sqrt((2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
        azimutal_integral = 0.0
        for n in range(1,N)
            if mp ≥ 0
                τp = cos(0.5*π*x[n]*mp)
            else
                τp = sin(-0.5*π*x[n]*mp)
            end
            if mq ≥ 0
                τq = cos(π*x[n]*mq)
            else
                τq = sin(-π*x[n]*mq)
            end
            azimutal_integral += 0.5 * π * w[n] * τp * τq
        end
        cosine_integral = 0.0
        for n in range(1,N)
            fp = 0.0
            for kp in range(0,lp-abs(mp))
                fp += binomial(lp-abs(mp),kp) * factorial_factor([lp+abs(mp)+kp],[abs(mp)+kp]) * (((x[n]+1)/2-1)/2)^(kp)
            end
            fp *= factorial_factor([lp,lp+abs(mp)],[lp-abs(mp),lp+abs(mp),lp],[(2,-abs(mp))]) * (1-((x[n]+1)/2)^2)^(abs(mp)/2)
            fq = 0.0
            for kq in range(0,lq-abs(mq))
                fq += binomial(lq-abs(mq),kq) * factorial_factor([lq+abs(mq)+kq],[abs(mq)+kq]) * ((x[n]-1)/2)^(kq)
            end
            fq *= factorial_factor([lq,lq+abs(mq)],[lq-abs(mq),lq+abs(mq),lq],[(2,-abs(mq))]) * (1-x[n]^2)^(abs(mq)/2)
            cosine_integral += w[n] * 0.5 * fp * fq
        end
        Mll[p,q] = Cp * Cq * azimutal_integral * cosine_integral
    end
    return Np,Mll
end

function octant_to_full_range_matrix_spherical_harmonics(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    Mll = zeros(Np,Np,2)
    N = 32
    x,w = gauss_legendre(N)
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        lq = pl[q]
        mp = pm[p]
        mq = pm[q]
        Cp = sqrt((2-(mp==0)) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
        Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
        azimutal_integral = zeros(2)
        for n in range(1,N)
            if mp ≥ 0
                τp = [cos(0.25*π*(x[n]+1)*mp),cos(0.25*π*(x[n]+3)*mp)]
            else
                τp = [sin(-0.25*π*(x[n]+1)*mp),sin(-0.25*π*(x[n]+3)*mp)]
            end
            if mq ≥ 0
                τq = cos(π*(x[n]+1)*mq)
            else
                τq = sin(-π*(x[n]+1)*mq)
            end
            azimutal_integral[1] += 0.25 * π * w[n] * τp[1] * τq
            azimutal_integral[2] += 0.25 * π * w[n] * τp[2] * τq
        end
        cosine_integral = 0.0
        for n in range(1,N)
            fp = 0.0
            for kp in range(0,lp-abs(mp))
                fp += binomial(lp-abs(mp),kp) * factorial_factor([lp+abs(mp)+kp],[abs(mp)+kp]) * (((x[n]+1)/2-1)/2)^(kp)
            end
            fp *= factorial_factor([lp,lp+abs(mp)],[lp-abs(mp),lp+abs(mp),lp],[(2,-abs(mp))]) * (1-((x[n]+1)/2)^2)^(abs(mp)/2)
            fq = 0.0
            for kq in range(0,lq-abs(mq))
                fq += binomial(lq-abs(mq),kq) * factorial_factor([lq+abs(mq)+kq],[abs(mq)+kq]) * ((x[n]-1)/2)^(kq)
            end
            fq *= factorial_factor([lq,lq+abs(mq)],[lq-abs(mq),lq+abs(mq),lq],[(2,-abs(mq))]) * (1-x[n]^2)^(abs(mq)/2)
            cosine_integral += w[n] * 0.5 * fp * fq
        end
        Mll[p,q,1] = Cp * Cq * azimutal_integral[1] * cosine_integral
        Mll[p,q,2] = Cp * Cq * azimutal_integral[2] * cosine_integral
    end
    return Np,Mll
end

function patch_to_full_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nq = spherical_harmonics_number_basis(Lq)
    pl_p,pm_p = spherical_harmonics_indices(L)
    pl_q,pm_q = spherical_harmonics_indices(Lq)
    Mll = zeros(Np,Nq,8,Nv,Nv)
    N = 32
    x,weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)

    sx = [1,1,1,1,-1,-1,-1,-1]
    sy = [1,1,-1,-1,1,1,-1,-1]
    sz = [1,-1,1,-1,1,-1,1,-1]

    # Pré-calculs sur les harmoniques (constants pour toute la fonction)
    Cp = Vector{Float64}(undef, Np)
    for p in 1:Np
        lp = pl_p[p]
        mp = pm_p[p]
        Cp[p] = sqrt((2 - (mp == 0)) * factorial_factor([lp - abs(mp)], [lp + abs(mp)]))
    end
    Cq = Vector{Float64}(undef, Nq)
    for q in 1:Nq
        lq = pl_q[q]
        mq = pm_q[q]
        Cq[q] = sqrt(2 * (2 - (mq == 0)) / π * (2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
    end

    # Indices uniques en m pour factoriser l'intégrale azimutale
    mp_unique = sort!(collect(Set(pm_p)))
    mq_unique = sort!(collect(Set(pm_q)))
    mp_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(mp_unique))
    mq_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(mq_unique))
    mp_idx_p = [mp_to_idx[pm_p[p]] for p in 1:Np]
    mq_idx_q = [mq_to_idx[pm_q[q]] for q in 1:Nq]

    # Offsets azimutaux par octant (ne dépend que de u)
    ϕ_offset_u = Vector{Float64}(undef, 8)
    for u in 1:8
        ϕ_offset_u[u] = (π / 2) * (2 + (sy[u] + 1) / 2 - (sz[u] + 1) / 2 - (sy[u] + 1) * (sz[u] + 1) / 2)
    end

    # Pré-calcul de 𝒯m(mq, 4*φq) avec 4*φq = π*(x+1)
    Tq = Matrix{Float64}(undef, length(mq_unique), N)
    for (imq, mq) in enumerate(mq_unique)
        for n in 1:N
            ϕ = π * (x[n] + 1.0)
            Tq[imq, n] = (mq ≥ 0) ? cos(mq * ϕ) : sin(-mq * ϕ)
        end
    end

    # Pré-calcul de Pnm(lq,|mq|, ±x) : pour sx=1 => x, pour sx=-1 => -x
    Pq_pos = Matrix{Float64}(undef, Nq, N)
    Pq_neg = Matrix{Float64}(undef, Nq, N)
    for q in 1:Nq
        lq = pl_q[q]
        mq = pm_q[q]
        amq = abs(mq)
        for n in 1:N
            Pq_pos[q, n] = Pnm(lq, amq, x[n])
            Pq_neg[q, n] = Pnm(lq, amq, -x[n])
        end
    end

    # Buffers (réutilisés pour limiter les allocations)
    Pp = Matrix{Float64}(undef, Np, N)
    PpW = Matrix{Float64}(undef, Np, N)
    cosine_mat = Matrix{Float64}(undef, Np, Nq)

    Tmp = Matrix{Float64}(undef, length(mp_unique), N)
    TmpW = Matrix{Float64}(undef, length(mp_unique), N)
    az_table = Matrix{Float64}(undef, length(mp_unique), length(mq_unique))

    # Helper local pour μ(v) : évite la closure μ_uv et les recalculs inutiles
    denom = Float64(Nv * (Nv + 1))
    mu_of_v(s::Int, v::Int) = (s == 1) ? (1.0 - ((Nv + 1 - v) * (Nv + 2 - v)) / denom) : (-1.0 + ((v - 1) * v) / denom)

    for u in 1:8
        su = sx[u]
        ϕ_offset = ϕ_offset_u[u]
        Pq_sign = (su == 1) ? Pq_pos : Pq_neg

        for v in 1:Nv
            # Δμ et Δϕ sont constants pour un couple (u,v)
            μ0 = mu_of_v(su, v)
            μ1 = mu_of_v(su, v + 1)
            Δμ = μ1 - μ0

            Nw = (su == 1) ? (Nv + 1 - v) : v
            Δϕ = (π / 2) / Nw

            # Pré-calcul de l'intégrale en cosinus pour tous (p,q) (indépendant de w)
            for p in 1:Np
                lp = pl_p[p]
                mp = pm_p[p]
                amp = abs(mp)
                for n in 1:N
                    μp = μ0 + Δμ * t[n]
                    Pp[p, n] = Pnm(lp, amp, μp)
                end
            end
            for n in 1:N
                wn = weight[n]
                for p in 1:Np
                    PpW[p, n] = Pp[p, n] * wn
                end
            end
            LinearAlgebra.mul!(cosine_mat, PpW, LinearAlgebra.transpose(Pq_sign))

            # Facteur (Δμ*Δϕ)/4 * sqrt(π/(2ΔμΔϕ)) = (1/4)*sqrt(π*Δμ*Δϕ/2)
            scale_uv = 0.25 * sqrt((π / 2) * (Δμ * Δϕ))

            ϕ0 = ϕ_offset
            for w in 1:Nw
                # Matrice Tmp (mp_unique × N) : 𝒯m(mp, ϕp) avec ϕp = ϕ0 + Δϕ*t
                for (imp, mp) in enumerate(mp_unique)
                    if mp ≥ 0
                        for n in 1:N
                            Tmp[imp, n] = cos(mp * (ϕ0 + Δϕ * t[n]))
                        end
                    else
                        for n in 1:N
                            Tmp[imp, n] = sin((-mp) * (ϕ0 + Δϕ * t[n]))
                        end
                    end
                end
                for n in 1:N
                    wn = weight[n]
                    for imp in 1:length(mp_unique)
                        TmpW[imp, n] = Tmp[imp, n] * wn
                    end
                end

                # az_table[imp, imq] = Σ_n weight[n] * 𝒯m(mp,ϕp_n) * 𝒯m(mq,π*(x_n+1))
                LinearAlgebra.mul!(az_table, TmpW, LinearAlgebra.transpose(Tq))

                # Remplissage de la tranche Mll[:,:,u,v,w]
                for p in 1:Np
                    scale_p = scale_uv * Cp[p]
                    imp = mp_idx_p[p]
                    for q in 1:Nq
                        imq = mq_idx_q[q]
                        Mll[p, q, u, v, w] += scale_p * Cq[q] * az_table[imp, imq] * cosine_mat[p, q]
                    end
                end

                ϕ0 += Δϕ
            end
        end
    end
    return Np,Nq,Mll
end

function 𝒯m(m::Int64,φ::Real)
    if m ≥ 0
        return cos(m*φ)
    else
        return sin(-m*φ)
    end
end

function Pnm(n::Int64,m::Int64,x::Real)
    return Pnmαβ(n,m,0,0,x)
end

function Pnmαβ(n::Int64,m::Int64,α::Int64,β::Int64,x::Real)
    return (1-x)^((m+α)/2) * (1+x)^((m+β)/2) * factorial_factor([α+β+n+m],[α+β+n],[(2,-m)]) * Pnαβ(n-m,α+m,β+m,x)
end

function Pnαβ(n::Int64,α::Int64,β::Int64,x::Real)
    Pnαβ = 0
    for k in range(0,n)
       Pnαβ += binomial(n,k) * factorial_factor([α+β+n+k],[α+k]) * ((x-1)/2)^k
    end
    Pnαβ *= factorial_factor([α+n],[n,α+β+n])
    return Pnαβ
end

function patch_to_half_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nq = spherical_harmonics_number_basis(Lq)
    pl_p,pm_p = spherical_harmonics_indices(L)
    pl_q,pm_q = spherical_harmonics_indices(Lq)
    Mll = zeros(Np,Nq,8,Nv,Nv,2*Ndims,2)
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

    # Pré-calcul des constantes harmoniques
    Cp = Vector{Float64}(undef, Np)
    for p in 1:Np
        lp = pl_p[p]
        mp = pm_p[p]
        Cp[p] = sqrt((2 - (mp == 0)) / (2 * π) * (2 * lp + 1) * factorial_factor([lp - abs(mp)], [lp + abs(mp)]))
    end
    Cq = Vector{Float64}(undef, Nq)
    for q in 1:Nq
        lq = pl_q[q]
        mq = pm_q[q]
        Cq[q] = sqrt(2 * (2 - (mq == 0)) / π * (2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
    end

    # Indices uniques en m pour factoriser l'azimutal
    mp_unique = sort!(collect(Set(pm_p)))
    mq_unique = sort!(collect(Set(pm_q)))
    mp_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(mp_unique))
    mq_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(mq_unique))
    mp_idx_p = [mp_to_idx[pm_p[p]] for p in 1:Np]
    mq_idx_q = [mq_to_idx[pm_q[q]] for q in 1:Nq]

    # 𝒯m(mq,4*ϕq) avec 4*ϕq = π*(x+1) est indépendant de (u,v,w,b,is)
    Tq = Matrix{Float64}(undef, length(mq_unique), N)
    for (imq, mq) in enumerate(mq_unique)
        if mq ≥ 0
            for n in 1:N
                Tq[imq, n] = cos(mq * (π * (x[n] + 1.0)))
            end
        else
            for n in 1:N
                Tq[imq, n] = sin((-mq) * (π * (x[n] + 1.0)))
            end
        end
    end

    # Pré-calcul de Pnm(lq,|mq|, ±x) car 2*μq-1 = sx[u]*x[n]
    Pq_pos = Matrix{Float64}(undef, Nq, N)
    Pq_neg = Matrix{Float64}(undef, Nq, N)
    for q in 1:Nq
        lq = pl_q[q]
        mq = pm_q[q]
        amq = abs(mq)
        for n in 1:N
            Pq_pos[q, n] = Pnm(lq, amq, x[n])
            Pq_neg[q, n] = Pnm(lq, amq, -x[n])
        end
    end

    # Buffers réutilisés
    μi_vals = Vector{Float64}(undef, N)
    Pp = Matrix{Float64}(undef, Np, N)
    PpW = Matrix{Float64}(undef, Np, N)
    cosine_mat = Matrix{Float64}(undef, Np, Nq)

    Tmp = Matrix{Float64}(undef, length(mp_unique), N)
    TmpW = Matrix{Float64}(undef, length(mp_unique), N)
    az_table = Matrix{Float64}(undef, length(mp_unique), length(mq_unique))

    # Helper μ(v) (identique à celui utilisé dans patch_to_full_range)
    denom = Float64(Nv * (Nv + 1))
    mu_of_v(s::Int, v::Int) = (s == 1) ? (1.0 - ((Nv + 1 - v) * (Nv + 2 - v)) / denom) : (-1.0 + ((v - 1) * v) / denom)

    for is in 1:2
        s = (is == 1) ? -1 : 1
        for b in 1:(2 * Ndims)
            μ⁻, μ⁺, ϕ⁻, ϕ⁺, sb = cartesian_surface_source(b, s)
            Δμ = μ⁺ - μ⁻
            Δϕ = ϕ⁺ - ϕ⁻
            boundary_scale = sqrt(2 * π / (Δμ * Δϕ))
            wrap = (b == 3 && s == -1) || (b == 4 && s == 1)

            for u in 1:8
                # Filtre de face boundary
                if !((b ∈ [1, 2] && sx[u] == sb * s) || (b ∈ [3, 4] && sy[u] == sb * s) || (b ∈ [5, 6] && sz[u] == sb * s))
                    continue
                end

                su = sx[u]
                Pq_sign = (su == 1) ? Pq_pos : Pq_neg
                ϕ_offset = ϕ_offset_u[u]

                for v in 1:Nv
                    # Patch sizes (u,v)
                    Nw = (su == 1) ? (Nv + 1 - v) : v
                    Δϕ_uw = (π / 2) / Nw

                    μ0 = mu_of_v(su, v)
                    μ1 = mu_of_v(su, v + 1)
                    Δμ_uv = μ1 - μ0

                    # Pré-calcul μi aux points de quadrature
                    for n in 1:N
                        μi_vals[n] = μ0 + Δμ_uv * t[n]
                    end

                    # Pré-calcul de l'intégrale cosinus pour tous (p,q) (indépendant de w)
                    for p in 1:Np
                        lp = pl_p[p]
                        mp = pm_p[p]
                        amp = abs(mp)
                        for n in 1:N
                            μi = μi_vals[n]
                            μp = -s * sb * (((-s * sb - 1) / 2) + (μi - μ⁻) / Δμ)
                            Pp[p, n] = Pnm(lp, amp, 2 * μp - 1)
                        end
                    end
                    for n in 1:N
                        wn = weight[n]
                        for p in 1:Np
                            PpW[p, n] = Pp[p, n] * wn
                        end
                    end
                    LinearAlgebra.mul!(cosine_mat, PpW, LinearAlgebra.transpose(Pq_sign))

                    # Facteur combiné : boundary_scale * (Δμ_uv*Δϕ_uw)/4 * sqrt(π/(2Δμ_uvΔϕ_uw))
                    # = boundary_scale * 0.25 * sqrt((π/2) * Δμ_uv * Δϕ_uw)
                    scale_uv = boundary_scale * 0.25 * sqrt((π / 2) * (Δμ_uv * Δϕ_uw))

                    # Boucle w : seulement l'azimutal dépend de w
                    ϕ_base = ϕ_offset
                    for w in 1:Nw
                        for n in 1:N
                            ϕi = ϕ_base + Δϕ_uw * t[n]
                            if wrap && (3 * π / 2 ≤ ϕi ≤ 2 * π)
                                ϕi -= 2 * π
                            end
                            ϕp = 2 * π * (ϕi - ϕ⁻) / Δϕ
                            for (imp, mp) in enumerate(mp_unique)
                                if mp ≥ 0
                                    Tmp[imp, n] = cos(mp * ϕp)
                                else
                                    Tmp[imp, n] = sin((-mp) * ϕp)
                                end
                            end
                        end

                        for n in 1:N
                            wn = weight[n]
                            for imp in 1:length(mp_unique)
                                TmpW[imp, n] = Tmp[imp, n] * wn
                            end
                        end
                        LinearAlgebra.mul!(az_table, TmpW, LinearAlgebra.transpose(Tq))

                        for p in 1:Np
                            scale_p = scale_uv * Cp[p]
                            imp = mp_idx_p[p]
                            for q in 1:Nq
                                imq = mq_idx_q[q]
                                Mll[p, q, u, v, w, b, is] += scale_p * Cq[q] * az_table[imp, imq] * cosine_mat[p, q]
                            end
                        end

                        ϕ_base += Δϕ_uw
                    end
                end
            end
        end
    end
    return Mll
end

function reflection_matrix(L::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    Mll = zeros(Np,Np,2*Ndims)
    N = 32
    x,weight = gauss_legendre(N)
    sp = -1
    sq = 1
    for b in range(1,2*Ndims)
        μ⁻_p,μ⁺_p,ϕ⁻_p,ϕ⁺_p,sb_p = cartesian_surface_source(b,sp)
        μ⁻_q,μ⁺_q,ϕ⁻_q,ϕ⁺_q,sb_q = cartesian_surface_source(b,sq)
        if sb_p != sb_q error("Error in boundary identification") end
        sb = sb_p
        Δμ_p = μ⁺_p-μ⁻_p
        Δμ_q = μ⁺_q-μ⁻_q
        Δϕ_p = ϕ⁺_p-ϕ⁻_p
        Δϕ_q = ϕ⁺_q-ϕ⁻_q
        boundary_scale_p = sqrt(2*π/(Δμ_p*Δϕ_p))
        boundary_scale_q = sqrt(2*π/(Δμ_q*Δϕ_q))
        for p in range(1,Np), q in range(1,Np)
            lp = pl[p]
            lq = pl[q]
            mp = pm[p]
            mq = pm[q]
            Cp = sqrt((2-(mp==0)) / (2*π) * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
            Cq = sqrt((2-(mq==0)) / (2*π) * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))

            # Azimuthal integral
            azimutal_integral = 0
            for n in range(1,N)
                ϕi = ϕ⁻_p + 0.5 * Δϕ_p * (x[n]+1)
                if (b == 3 && (3*π/2 ≤ ϕi ≤ 2*π) && sp == -1) || (b == 4 && (3*π/2 ≤ ϕi ≤ 2*π) && sp == 1) ϕi -= 2*π end
                ϕp = ϕi
                ϕq = ϕi
                if (b == 3 || b == 4) 
                    ϕq = π - ϕq 
                elseif (b == 5 || b == 6) 
                    ϕq = 2*π - ϕq
                end
                ϕp = 2*π * (ϕp-ϕ⁻_p)/Δϕ_p
                ϕq = 2*π * (ϕq-ϕ⁻_q)/Δϕ_q
                azimutal_integral += weight[n] * 𝒯m(mp,ϕp) * 𝒯m(mq,ϕq)
            end

            # Cosine integral
            cosine_integral = 0.0
            for n in range(1,N)
                μi = μ⁻_p + 0.5 * Δμ_p * (x[n]+1)
                μp = μi
                μq = μi
                if (b == 1 || b == 2) μq = -μi end
                μp = sp*sb * ((sp*sb-1)/2+(μp-μ⁻_p)/Δμ_p)
                μq = sq*sb * ((sq*sb-1)/2+(μq-μ⁻_q)/Δμ_q)
                cosine_integral += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μq-1)
            end

            # Update matrix
            Mll[p,q,b] += (Δμ_p*Δϕ_p)/4 * boundary_scale_p * boundary_scale_q * Cp * Cq * azimutal_integral * cosine_integral

        end
    end
    return Mll
end

# function patch_to_full_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64)
#     Np = spherical_harmonics_number_basis(L)
#     Nq = spherical_harmonics_number_basis(Lq)
#     pl_p,pm_p = spherical_harmonics_indices(L)
#     pl_q,pm_q = spherical_harmonics_indices(Lq)
#     Mll = zeros(Np,Nq,8,Nv,Nv)
#     N = 32
#     x,weight = gauss_legendre(N)
#     sx = [1,1,1,1,-1,-1,-1,-1]
#     sy = [1,1,-1,-1,1,1,-1,-1]
#     sz = [1,-1,1,-1,1,-1,1,-1]
#     μ_uv(u,v) = sx[u]*(1-(1-v+(sx[u]+1)/2*Nv)*(-v+(sx[u]+1)/2*(Nv+2))/(Nv*(Nv+1)))
#     ϕ_uw(u,v,w) = (π/2)*(w-1)/(-sx[u]*v + (sx[u]+1)/2*(Nv+1)) + π/2 * (2 + (sy[u]+1)/2 - (sz[u]+1)/2 - (sy[u]+1)*(sz[u]+1)/2)
#     for p in range(1,Np), q in range(1,Nq)
#         lp = pl_p[p]
#         lq = pl_q[q]
#         mp = pm_p[p]
#         mq = pm_q[q]
#         Cp = sqrt((2-(mp==0)) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
#         Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
#         for u in range(1,8)
#             for v in range(1,Nv)
#                 Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
#                 for w in range(1,Nw)
#                     Δμ_uv = μ_uv(u,v+1) - μ_uv(u,v)
#                     Δϕ_uw = ϕ_uw(u,v,w+1) - ϕ_uw(u,v,w)

#                     # Azimuthal integral
#                     azimutal_integral = 0
#                     for n in range(1,N)
#                         φp = ϕ_uw(u,v,w) + 0.5 * Δϕ_uw * (x[n]+1)
#                         φq = π/2 * (φp-ϕ_uw(u,v,w))/Δϕ_uw
#                         azimutal_integral += weight[n] * 𝒯m(mp,φp) * 𝒯m(mq,4*φq)
#                     end

#                     # Cosine integral
#                     cosine_integral = 0.0
#                     for n in range(1,N)
#                         μp = μ_uv(u,v) + 0.5 * Δμ_uv * (x[n]+1)
#                         μq = sx[u] * ((sx[u]-1)/2+(μp-μ_uv(u,v))/Δμ_uv)
#                         cosine_integral += weight[n] * Pnm(lp,abs(mp),μp) * Pnm(lq,abs(mq),2*μq-1)
#                     end

#                     # Factor
#                     Cuvw = sqrt(π/(2*Δμ_uv*Δϕ_uw))

#                     Mll[p,q,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral * cosine_integral

#                 end
#             end
#         end
#     end
#     return Np,Nq,Mll
# end

# function patch_to_half_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64)
#     Np = spherical_harmonics_number_basis(L)
#     Nq = spherical_harmonics_number_basis(Lq)
#     pl_p,pm_p = spherical_harmonics_indices(L)
#     pl_q,pm_q = spherical_harmonics_indices(Lq)
#     Mll = zeros(Np,Nq,8,Nv,Nv,2*Ndims,2)
#     N = 32
#     x,weight = gauss_legendre(N)
#     sx = [1,1,1,1,-1,-1,-1,-1]
#     sy = [1,1,-1,-1,1,1,-1,-1]
#     sz = [1,-1,1,-1,1,-1,1,-1]
#     μ_uv(u,v) = sx[u]*(1-(1-v+(sx[u]+1)/2*Nv)*(-v+(sx[u]+1)/2*(Nv+2))/(Nv*(Nv+1)))
#     ϕ_uw(u,v,w) = (π/2)*(w-1)/(-sx[u]*v + (sx[u]+1)/2*(Nv+1)) + π/2 * (2 + (sy[u]+1)/2 - (sz[u]+1)/2 - (sy[u]+1)*(sz[u]+1)/2)
#     for is in range(1,2)
#         if (is == 1) sd = -1 else sd = 1 end
#         for b in range(1,2*Ndims)
#             μ⁻,μ⁺,ϕ⁻,ϕ⁺,sb = cartesian_surface_source(b,sd)
#             Δμ = μ⁺-μ⁻
#             Δϕ = ϕ⁺-ϕ⁻
#             for p in range(1,Np), q in range(1,Nq)
#                 lp = pl_p[p]
#                 lq = pl_q[q]
#                 mp = pm_p[p]
#                 mq = pm_q[q]
#                 Cp = sqrt((2-(mp==0))/(2*π) * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
#                 Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
#                 for u in range(1,8)
#                     if ~((b ∈ [1,2] && sx[u] == sb*sd) || (b ∈ [3,4] && sy[u] == sb*sd) || (b ∈ [5,6] && sz[u] == sb*sd)) continue end
#                     for v in range(1,Nv)
#                         Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
#                         for w in range(1,Nw)
#                             Δμ_uv = μ_uv(u,v+1) - μ_uv(u,v)
#                             Δϕ_uw = ϕ_uw(u,v,w+1) - ϕ_uw(u,v,w)

#                             # Azimuthal integral
#                             azimutal_integral = 0.0
#                             for n in range(1,N)
#                                 ϕi = ϕ_uw(u,v,w) + 0.5 * Δϕ_uw * (x[n]+1)
#                                 μi = μ_uv(u,v) + 0.5 * Δμ_uv * (x[n]+1)
#                                 if b == 1 || b == 2
#                                     ϕp = mod(ϕi, 2π)
#                                 elseif b == 3 || b == 4
#                                     ξi = sqrt(1-μi^2) * cos(ϕi)
#                                     ϕp = mod(atan(ξi, μi), 2π)
#                                 else
#                                     ηi = sqrt(1-μi^2) * sin(ϕi)
#                                     ϕp = mod(atan(ηi, μi), 2π)
#                                 end
#                                 ϕq = (π/2) * (ϕi-ϕ_uw(u,v,w))/Δϕ_uw
#                                 azimutal_integral += weight[n] * 𝒯m(mp,ϕp) * 𝒯m(mq,4*ϕq)
#                             end

#                             # Cosine integral
#                             cosine_integral = 0.0
#                             for n in range(1,N)
#                                 ϕi = ϕ_uw(u,v,w) + 0.5 * Δϕ_uw * (x[n]+1)
#                                 μi = μ_uv(u,v) + 0.5 * Δμ_uv * (x[n]+1)
#                                 if b == 1 || b == 2
#                                     μp = sd*sb*μi
#                                 elseif b == 3 || b == 4
#                                     μp = sd*sb*sqrt(1-μi^2)*cos(ϕi)
#                                 else
#                                     μp = sd*sb*sqrt(1-μi^2)*sin(ϕi)
#                                 end 
#                                 μq = sx[u] * ((sx[u]-1)/2+(μi-μ_uv(u,v))/Δμ_uv)
#                                 cosine_integral += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μq-1)
#                             end

#                             # Factor
#                             Cuvw = sqrt(π/(2*Δμ_uv*Δϕ_uw))

#                             # Update matrix
#                             Mll[p,q,u,v,w,b,is] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral * cosine_integral

#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return Mll
# end

function pos_to_neg_half_range_matrix_spherical_harmonics(L::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    Rpq = zeros(Np,Np,2*Ndims)
    for b in range(1,2*Ndims)
        for p in range(1,Np)
            if b ∈ [1,2]
                Rpq[p,p,b] = 1.0
            else
                # Exact parity for real spherical-harmonics reflection on Y/Z faces.
                parity = pl[p] + abs(pm[p]) + (pm[p] < 0)
                Rpq[p,p,b] = iseven(parity) ? 1.0 : -1.0
            end
        end
    end
    return Rpq
end