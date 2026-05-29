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

function patch_to_full_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64;tiling::String="polar-anchored")
    if tiling == "symmetric"
        return patch_to_full_range_matrix_spherical_harmonics_symmetric(L,Lq,Nv)
    end
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

"""
    patch_to_full_range_matrix_legendre(L::Int64, L_elem::Int64, Nv::Int64)

Build the patch-to-full-range moment-transfer matrix `Mll[p, q, u, v, 1]` for the
1D Legendre (azimuthally symmetric) GN basis, the Legendre counterpart of
`patch_to_full_range_matrix_spherical_harmonics`. Here
`Mll[p,q,u,v,1] = ∫_band P_{l_p}(μ) ψ_q^loc(μ) dμ`, with the patch-orthonormal
Legendre basis `ψ_q^loc(μ) = sqrt((2 q_l + 1)/Δ) · P_{q_l}(ξ)`. Only octants
`u ∈ {1, 5}` carry patches. With `Nv = 1` this reproduces
`half_to_full_range_matrix_legendre` (per hemisphere).

# Output Argument(s)
- `Np::Int64` : number of full-range moments (`L + 1`).
- `Nq::Int64` : number of per-patch moments (`L_elem + 1`).
- `Mll::Array{Float64,5}` : moment-transfer matrix, shape `(Np, Nq, 8, Nv, 1)`.
"""
function patch_to_full_range_matrix_legendre(L::Int64, L_elem::Int64, Nv::Int64)
    if L < 0 || L_elem < 0 error("Legendre orders must be ≥ 0.") end
    if Nv <= 0 error("Number of μ-bands per hemisphere must be > 0.") end
    Np = L + 1
    Nq = L_elem + 1
    Mll = zeros(Np, Nq, 8, Nv, 1)
    N = 32
    x, weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)
    Ploc = [legendre_polynomials_up_to_L(L_elem, x[n]) for n in 1:N]
    for u in (1, 5), v in 1:Nv
        μ0, μ1 = _gn_legendre_band(u, v, Nv)
        Δ = μ1 - μ0
        for n in 1:N
            μ = μ0 + Δ * t[n]
            Pfull = legendre_polynomials_up_to_L(L, μ)
            wn = weight[n] * (Δ / 2)
            for q in 1:Nq
                ψq = sqrt((2 * (q - 1) + 1) / Δ) * Ploc[n][q]
                for p in 1:Np
                    Mll[p, q, u, v, 1] += wn * Pfull[p] * ψq
                end
            end
        end
    end
    return Np, Nq, Mll
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

function patch_to_half_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64;tiling::String="polar-anchored")
    if tiling == "symmetric"
        return patch_to_half_range_matrix_spherical_harmonics_symmetric(L,Lq,Nv,Ndims)
    end
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

"""
    patch_to_half_range_matrix_legendre(L_surf::Int64, L_elem::Int64, Nv::Int64, Ndims::Int64)

Build the patch-to-half-range surface moment-transfer matrix
`Mll[p, q, u, v, 1, b, is]` for the 1D Legendre GN basis (`Ndims == 1`), the
Legendre counterpart of `patch_to_half_range_matrix_spherical_harmonics`. Here
`Mll[p,q,u,v,1,b,is] = ∫_band R̄_p(sx[u]·μ) ψ_q^loc(μ) dμ`, where the half-range
surface basis `R̄_p(μ̂) = sqrt(2 p_l + 1) · P_{p_l}(2 μ̂ - 1)` over `μ̂ = sx[u]·μ ∈ [0,1]`
matches `half_range_legendre_polynomials_up_to_L`. `is = 1` (incoming, s = -1)
and `is = 2` (outgoing, s = +1) select the hemisphere through the octant filter
`sx[u] == sb·s` (`sb = -1` for face `b = 1`, `sb = +1` for face `b = 2`).

# Output Argument(s)
- `Mll::Array{Float64,7}` : shape `(L_surf+1, L_elem+1, 8, Nv, 1, 2*Ndims, 2)`.
"""
function patch_to_half_range_matrix_legendre(L_surf::Int64, L_elem::Int64, Nv::Int64, Ndims::Int64)
    if Ndims != 1 error("Legendre half-range surface matrix is only available in 1D.") end
    Np = L_surf + 1
    Nq = L_elem + 1
    Mll = zeros(Np, Nq, 8, Nv, 1, 2 * Ndims, 2)
    N = 32
    x, weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)
    sb_of = (-1, 1)  # cartesian_surface_source sb for faces b = 1, 2
    norm_half = sqrt.(2 .* (0:L_surf) .+ 1)
    Ploc = [legendre_polynomials_up_to_L(L_elem, x[n]) for n in 1:N]
    for is in 1:2
        s = (is == 1) ? -1 : 1
        for b in 1:(2 * Ndims)
            sb = sb_of[b]
            for u in (1, 5)
                σ = _GN_SX[u]
                if σ != sb * s continue end
                for v in 1:Nv
                    μ0, μ1 = _gn_legendre_band(u, v, Nv)
                    Δ = μ1 - μ0
                    for n in 1:N
                        μ = μ0 + Δ * t[n]
                        Rhalf = norm_half .* legendre_polynomials_up_to_L(L_surf, 2 * σ * μ - 1)
                        wn = weight[n] * (Δ / 2)
                        for q in 1:Nq
                            ψq = sqrt((2 * (q - 1) + 1) / Δ) * Ploc[n][q]
                            for p in 1:Np
                                Mll[p, q, u, v, 1, b, is] += wn * Rhalf[p] * ψq
                            end
                        end
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

"""
    pos_to_neg_half_range_matrix_legendre(L_surf::Int64, Ndims::Int64)

Reflection matrix `Rpq[p, q, b]` mapping outgoing to incoming half-range surface
moments for the 1D Legendre GN basis, the Legendre counterpart of
`pos_to_neg_half_range_matrix_spherical_harmonics`. In the per-hemisphere
half-range Legendre basis, the μ → -μ reflection on an x-face is the identity
(matching the spherical-harmonics result on the x-faces).

# Output Argument(s)
- `Rpq::Array{Float64,3}` : shape `(L_surf+1, L_surf+1, 2*Ndims)`.
"""
function pos_to_neg_half_range_matrix_legendre(L_surf::Int64, Ndims::Int64)
    Np = L_surf + 1
    Rpq = zeros(Np, Np, 2 * Ndims)
    for b in 1:(2 * Ndims), p in 1:Np
        Rpq[p, p, b] = 1.0
    end
    return Rpq
end

"""
    patch_to_full_range_matrix_spherical_harmonics_symmetric(L, Lq, Nv)

Build the patch-to-full-range moment-transfer matrix `Mll[p, q, u, i, j]` for
the symmetric tiling. The full-range basis is indexed by `p`, the patch (octant
sub-triangle) basis by `q`. Storage shape: `(Np, Nq, 8, Nv, 2Nv-1)`.

Integration is by 32-point Duffy-chart Gauss-Legendre per sub-triangle.
"""
function patch_to_full_range_matrix_spherical_harmonics_symmetric(L::Int64,Lq::Int64,Nv::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nq = spherical_harmonics_number_basis(Lq)
    Nw_max = 2*Nv - 1
    Mll = zeros(Np,Nq,8,Nv,Nw_max)
    N = 32
    x,weight = gauss_legendre(N)
    Nquad = N * N

    Yfull = Matrix{Float64}(undef, Np, Nquad)

    for u in 1:8, i in 1:Nv
        for j in 1:(2i - 1)
            Ψ_orth, Px, Py, Pz, JW = symmetric_patch_orthonormal_basis(Lq, Nv, u, i, j, x, weight)

            for k in 1:Nquad
                μ_g = Px[k]
                ϕ_g = atan(Pz[k], Py[k])
                if ϕ_g < 0.0 ϕ_g += 2π end
                Yf = real_spherical_harmonics_up_to_L(L, μ_g, ϕ_g)
                wfac = JW[k]
                p = 0
                for l in 0:L, m in -l:l
                    p += 1
                    Yfull[p, k] = Yf[l+1][l+m+1] * wfac
                end
            end

            LinearAlgebra.mul!(view(Mll, :, :, u, i, j), Yfull, LinearAlgebra.transpose(Ψ_orth))
        end
    end

    return Np, Nq, Mll
end

"""
    patch_to_half_range_matrix_spherical_harmonics_symmetric(L, Lq, Nv, Ndims)

Build the patch-to-half-range moment-transfer matrix `Mll[p, q, u, i, j, b, is]`
for the symmetric tiling. Shape: `(Np, Nq, 8, Nv, 2Nv-1, 2*Ndims, 2)`.

For each boundary face `b` and orientation `is`, only octants whose sign on the
boundary axis matches contribute. Within an octant, every sub-triangle is fully
on the correct side, so the same octant filter is sufficient.
"""
function patch_to_half_range_matrix_spherical_harmonics_symmetric(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nq = spherical_harmonics_number_basis(Lq)
    Nw_max = 2*Nv - 1
    Mll = zeros(Np,Nq,8,Nv,Nw_max,2*Ndims,2)
    N = 32
    x,weight = gauss_legendre(N)
    Nquad = N * N

    sx = [1,1,1,1,-1,-1,-1,-1]
    sy = [1,1,-1,-1,1,1,-1,-1]
    sz = [1,-1,1,-1,1,-1,1,-1]

    Yhalf = Vector{Float64}(undef, Np)
    Plm_buf = Matrix{Float64}(undef, L+1, L+1)
    Pl_buf  = Vector{Float64}(undef, L+1)

    for u in 1:8, i in 1:Nv
        for j in 1:(2i - 1)
            Ψ_orth, Px, Py, Pz, JW = symmetric_patch_orthonormal_basis(Lq, Nv, u, i, j, x, weight)

            for is in 1:2
                s_bnd = (is == 1) ? -1 : 1
                for b in 1:(2 * Ndims)
                    μ⁻, μ⁺, ϕ⁻, ϕ⁺, sb = cartesian_surface_source(b, s_bnd)

                    if !((b ∈ (1,2) && sx[u] == sb * s_bnd) || (b ∈ (3,4) && sy[u] == sb * s_bnd) || (b ∈ (5,6) && sz[u] == sb * s_bnd))
                        continue
                    end

                    Δμ_b = μ⁺ - μ⁻
                    Δϕ_b = ϕ⁺ - ϕ⁻
                    boundary_scale = sqrt(2 * π / (Δμ_b * Δϕ_b))
                    wrap = (b == 3 && s_bnd == -1) || (b == 4 && s_bnd == 1)

                    for k in 1:Nquad
                        μ_g = Px[k]
                        ϕ_g = atan(Pz[k], Py[k])
                        if ϕ_g < 0.0 ϕ_g += 2π end
                        if wrap && (3π/2 ≤ ϕ_g ≤ 2π) ϕ_g -= 2π end

                        μp_b = -s_bnd * sb * (((-s_bnd * sb - 1) / 2) + (μ_g - μ⁻) / Δμ_b)
                        μp_b = clamp(μp_b, 0.0, 1.0)
                        ϕp_b = 2 * π * (ϕ_g - ϕ⁻) / Δϕ_b

                        real_half_range_spherical_harmonics_up_to_L!(Yhalf, Plm_buf, Pl_buf, L, μp_b, ϕp_b)

                        wfac = boundary_scale * JW[k]
                        for q in 1:Nq
                            ψq_k = Ψ_orth[q, k]
                            for p in 1:Np
                                Mll[p, q, u, i, j, b, is] += wfac * Yhalf[p] * ψq_k
                            end
                        end
                    end
                end
            end
        end
    end

    return Mll
end