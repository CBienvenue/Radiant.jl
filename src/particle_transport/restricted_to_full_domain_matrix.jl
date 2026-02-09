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
    t01 = 0.5 .* (x .+ 1.0)

    # --- Helpers (allocation-free tables) ---
    function _fill_Plm_full!(Plm::Matrix{Float64}, L::Int64, μ::Float64)
        fill!(Plm, 0.0)
        if μ == -1 || μ == 1
            Pl = jacobi_polynomials_up_to_L(L, 0, 0, μ)
            for l in 0:L
                Plm[l+1, 1] = Pl[l+1]
            end
            return Plm
        end
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m, m, m, μ)
            for l in m:L
                Plm[l+1, m+1] = factorial_factor([l+m], [l], [(2,-m),(1-μ^2,m/2)]) * Pl[l-m+1]
            end
        end
        return Plm
    end

    function _fill_Plm_octant!(Plm::Matrix{Float64}, L::Int64, μ01::Float64)
        fill!(Plm, 0.0)
        t = 2*μ01 - 1
        if t == -1 || t == 1
            Pl = jacobi_polynomials_up_to_L(L, 0, 0, t)
            for l in 0:L
                Plm[l+1, 1] = Pl[l+1]
            end
            return Plm
        end
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m, m, m, t)
            for l in m:L
                Plm[l+1, m+1] = factorial_factor([l+m], [l], [(2,-m),(1-t^2,m/2)]) * Pl[l-m+1]
            end
        end
        return Plm
    end

    # Precompute constants per basis index
    abs_m_p = zeros(Int64, Np)
    C_full = zeros(Np)
    for p in 1:Np
        l = pl_p[p]
        m = pm_p[p]
        am = abs(m)
        abs_m_p[p] = am
        C_full[p] = sqrt((2-(m==0)) * factorial_factor([l-am], [l+am]))
    end

    abs_m_q = zeros(Int64, Nq)
    C_oct = zeros(Nq)
    for q in 1:Nq
        l = pl_q[q]
        m = pm_q[q]
        am = abs(m)
        abs_m_q[q] = am
        C_oct[q] = sqrt(2*(2-(m==0))/π * (2*l+1) * factorial_factor([l-am], [l+am]))
    end

    # Patch-range basis evaluation on normalized coordinates is independent of (u,v,w)
    Plm_oct = zeros(Lq+1, Lq+1)
    A_oct = zeros(Nq, N)
    for n in 1:N
        _fill_Plm_octant!(Plm_oct, Lq, t01[n])
        for q in 1:Nq
            lq = pl_q[q]
            A_oct[q,n] = C_oct[q] * Plm_oct[lq+1, abs_m_q[q]+1]
        end
    end

    T_oct = zeros(Nq, N)
    for m in 1:N
        ϕ = (π/4) * (x[m] + 1)
        for q in 1:Nq
            mq = pm_q[q]
            T_oct[q,m] = (mq ≥ 0) ? cos(4*mq*ϕ) : sin(4*abs(mq)*ϕ)
        end
    end

    # Work buffers (reused)
    Plm_full = zeros(L+1, L+1)
    A_full = zeros(Np, N)
    T_full = zeros(Np, N)
    Rvec = zeros(Np)
    ψvec = zeros(Nq)

    sx_u = (1,1,1,1,-1,-1,-1,-1)
    sy_u = (1,1,-1,-1,1,1,-1,-1)
    sz_u = (1,-1,1,-1,1,-1,1,-1)

    for u in 1:8
        sx = sx_u[u]
        sy = sy_u[u]
        sz = sz_u[u]

        sx1 = (sx + 1) ÷ 2
        sy1 = (sy + 1) ÷ 2
        sz1 = (sz + 1) ÷ 2
        ϕ_offset = (π/2) * (2 + sy1 - sz1 - 2*sy1*sz1)

        for v in 1:Nv
            # μ_v(v_) from real_patch_range_spherical_harmonics_up_to_L
            term1 = 1 - v + sx1*Nv
            term2 = -v + sx1*(Nv+2)
            μ0 = sx * (1 - (term1*term2) / (Nv*(Nv+1)))

            vp1 = v + 1
            term1p1 = 1 - vp1 + sx1*Nv
            term2p1 = -vp1 + sx1*(Nv+2)
            μ1 = sx * (1 - (term1p1*term2p1) / (Nv*(Nv+1)))
            Δμ = μ1 - μ0

            Nw = Int(-sx*v + sx1*(Nv+1))
            Δϕ = (π/2) / Nw

            for w in 1:Nw
                ϕ0 = (π/2) * ((w-1) / Nw) + ϕ_offset

                # Jacobian (ΔμΔϕ/4) combined with patch normalization sqrt(π/(2ΔμΔϕ))
                α_patch = sqrt(π) * sqrt(Δμ*Δϕ) / (4*sqrt(2))

                for n in 1:N
                    μn = μ0 + t01[n]*Δμ
                    _fill_Plm_full!(Plm_full, L, μn)
                    for p in 1:Np
                        lp = pl_p[p]
                        A_full[p,n] = C_full[p] * Plm_full[lp+1, abs_m_p[p]+1]
                    end
                end

                for m in 1:N
                    φm = ϕ0 + t01[m]*Δϕ
                    for p in 1:Np
                        mp = pm_p[p]
                        T_full[p,m] = (mp ≥ 0) ? cos(mp*φm) : sin(abs(mp)*φm)
                    end
                end

                for n in 1:N
                    wn = weight[n]
                    for m in 1:N
                        wnm = wn * weight[m] * α_patch

                        for p in 1:Np
                            Rvec[p] = A_full[p,n] * T_full[p,m]
                        end
                        for q in 1:Nq
                            ψvec[q] = A_oct[q,n] * T_oct[q,m]
                        end

                        for p in 1:Np
                            rp = Rvec[p]
                            for q in 1:Nq
                                Mll[p,q,u,v,w] += wnm * rp * ψvec[q]
                            end
                        end
                    end
                end
            end
        end
    end
    return Np,Nq,Mll
end

function patch_to_half_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nq = spherical_harmonics_number_basis(Lq)
    pl_p,pm_p = spherical_harmonics_indices(L)
    pl_q,pm_q = spherical_harmonics_indices(Lq)

    Mll = zeros(Np, Nq, 8, Nv, Nv, 2*Ndims)

    # Quadrature on [0,1]
    N = 32
    x, weight = gauss_legendre(N)
    t01 = 0.5 .* (x .+ 1.0)
    K = N * N

    # --- Allocation-free associated Legendre table on half-range coordinate μ01 ∈ [0,1] ---
    function _fill_Plm_01!(Plm::Matrix{Float64}, Lloc::Int64, μ01::Float64)
        fill!(Plm, 0.0)
        t = 2*μ01 - 1
        if t == -1 || t == 1
            Pl = jacobi_polynomials_up_to_L(Lloc, 0, 0, t)
            for l in 0:Lloc
                Plm[l+1, 1] = Pl[l+1]
            end
            return Plm
        end
        for m in 0:Lloc
            Pl = jacobi_polynomials_up_to_L(Lloc-m, m, m, t)
            for l in m:Lloc
                Plm[l+1, m+1] = factorial_factor([l+m], [l], [(2,-m),(1-t^2,m/2)]) * Pl[l-m+1]
            end
        end
        return Plm
    end

    # Basis constants for half-range (p) and octant-range (q)
    abs_m_p = zeros(Int64, Np)
    C_half = zeros(Np)
    for p in 1:Np
        l = pl_p[p]
        m = pm_p[p]
        am = abs(m)
        abs_m_p[p] = am
        C_half[p] = sqrt((2-(m == 0)) / (2*π) * (2*l+1) * factorial_factor([l-am], [l+am]))
    end

    abs_m_q = zeros(Int64, Nq)
    C_oct = zeros(Nq)
    for q in 1:Nq
        l = pl_q[q]
        m = pm_q[q]
        am = abs(m)
        abs_m_q[q] = am
        C_oct[q] = sqrt(2*(2-(m == 0))/π * (2*l+1) * factorial_factor([l-am], [l+am]))
    end

    # Precompute octant-range basis table with sqrt(weight) included:
    # base_oct[q,k] = sqrt(w_n*w_m) * octant_SH_q(μ01_n, ϕ01_m)
    Plm_oct = zeros(Lq+1, Lq+1)
    A_oct = zeros(Nq, N)
    for n in 1:N
        _fill_Plm_01!(Plm_oct, Lq, t01[n])
        for q in 1:Nq
            lq = pl_q[q]
            A_oct[q, n] = C_oct[q] * Plm_oct[lq+1, abs_m_q[q]+1]
        end
    end
    T_oct = zeros(Nq, N)
    for m in 1:N
        ϕ01 = (π/4) * (x[m] + 1)
        for q in 1:Nq
            mq = pm_q[q]
            T_oct[q, m] = (mq ≥ 0) ? cos(4*mq*ϕ01) : sin(4*abs(mq)*ϕ01)
        end
    end

    base_oct = zeros(Nq, K)
    k = 1
    for n in 1:N
        wn_sqrt = sqrt(weight[n])
        for m in 1:N
            wnm_sqrt = wn_sqrt * sqrt(weight[m])
            for q in 1:Nq
                base_oct[q, k] = wnm_sqrt * A_oct[q, n] * T_oct[q, m]
            end
            k += 1
        end
    end
    base_oct_T = zeros(K, Nq)
    for q in 1:Nq
        for kk in 1:K
            base_oct_T[kk, q] = base_oct[q, kk]
        end
    end

    # Geometry: octant signs per u
    sx_u = (1,1,1,1,-1,-1,-1,-1)
    sy_u = (1,1,-1,-1,1,1,-1,-1)
    sz_u = (1,-1,1,-1,1,-1,1,-1)

    # Buffers
    Plm_half = zeros(L+1, L+1)
    Hw = zeros(Np, K)
    Mtmp = zeros(Np, Nq)

    for nb in 1:(2*Ndims)

        for u in 1:8
            sx = sx_u[u]
            sy = sy_u[u]
            sz = sz_u[u]

            # Keep only directions outgoing from the considered boundary
            if nb == 1 && sx == -1
                continue
            elseif nb == 2 && sx == 1
                continue
            elseif nb == 3 && sy == -1
                continue
            elseif nb == 4 && sy == 1
                continue
            elseif nb == 5 && sz == -1
                continue
            elseif nb == 6 && sz == 1
                continue
            end

            sx1 = (sx + 1) ÷ 2
            sy1 = (sy + 1) ÷ 2
            sz1 = (sz + 1) ÷ 2
            ϕ_offset = (π/2) * (2 + sy1 - sz1 - 2*sy1*sz1)

            for v in 1:Nv
                term1 = 1 - v + sx1*Nv
                term2 = -v + sx1*(Nv+2)
                μ0 = sx * (1 - (term1*term2) / (Nv*(Nv+1)))

                vp1 = v + 1
                term1p1 = 1 - vp1 + sx1*Nv
                term2p1 = -vp1 + sx1*(Nv+2)
                μ1 = sx * (1 - (term1p1*term2p1) / (Nv*(Nv+1)))
                Δμ_uv = μ1 - μ0

                Nw = Int(-sx*v + sx1*(Nv+1))
                Δϕ_uv = (π/2) / Nw

                for w in 1:Nw
                    ϕ0 = (π/2) * ((w-1) / Nw) + ϕ_offset

                    # jacobian (ΔμΔϕ/4) times patch normalization sqrt(π/(2ΔμΔϕ))
                    α_patch = sqrt(π) * sqrt(Δμ_uv*Δϕ_uv) / (4*sqrt(2))

                    # Build Hw[p,k] = sqrt(w_n*w_m) * half_range_SH_p(μ⁺, ϕ⁺)
                    k = 1
                    for n in 1:N
                        μx = μ0 + t01[n]*Δμ_uv
                        s = 1 - μx*μx
                        if s < 0
                            s = 0.0
                        end
                        sinθ = sqrt(s)
                        wn_sqrt = sqrt(weight[n])
                        for m in 1:N
                            ϕ = ϕ0 + t01[m]*Δϕ_uv
                            η = sinθ * cos(ϕ)
                            ξ = sinθ * sin(ϕ)

                            μplus = 0.0
                            ϕplus = 0.0
                            if nb == 1
                                μplus = μx
                                ϕplus = atan(ξ, η)
                            elseif nb == 2
                                μplus = -μx
                                ϕplus = atan(ξ, η)
                            elseif nb == 3
                                μplus = η
                                ϕplus = atan(μx, ξ)
                            elseif nb == 4
                                μplus = -η
                                ϕplus = atan(μx, ξ)
                            elseif nb == 5
                                μplus = ξ
                                ϕplus = atan(η, μx)
                            elseif nb == 6
                                μplus = -ξ
                                ϕplus = atan(η, μx)
                            else
                                error("Invalid boundary index nb = $nb")
                            end

                            if μplus < 0
                                μplus = 0.0
                            elseif μplus > 1
                                μplus = 1.0
                            end

                            _fill_Plm_01!(Plm_half, L, μplus)
                            wnm_sqrt = wn_sqrt * sqrt(weight[m])
                            for p in 1:Np
                                lp = pl_p[p]
                                mp = pm_p[p]
                                A = C_half[p] * Plm_half[lp+1, abs_m_p[p]+1]
                                T = (mp ≥ 0) ? cos(mp*ϕplus) : sin(abs(mp)*ϕplus)
                                Hw[p, k] = wnm_sqrt * A * T
                            end
                            k += 1
                        end
                    end

                    # Mtmp = Hw * base_oct'  then scale by α_patch
                    mul!(Mtmp, Hw, base_oct_T)
                    for p in 1:Np
                        for q in 1:Nq
                            Mll[p, q, u, v, w, nb] = α_patch * Mtmp[p, q]
                        end
                    end
                end
            end
        end
    end

    return Mll
end

# function patch_to_half_range_matrix_spherical_harmonics(L::Int64,Lq::Int64,Nv::Int64,Ndims::Int64)
#     Np = spherical_harmonics_number_basis(L)
#     Nq = spherical_harmonics_number_basis(Lq)
#     pl_p,pm_p = spherical_harmonics_indices(L)
#     pl_q,pm_q = spherical_harmonics_indices(Lq)

#     Mll = zeros(Np,Nq,8,Nv,Nv,2*Ndims)

#     N = 32
#     x,weight = gauss_legendre(N)
#     sx = [1,1,1,1,-1,-1,-1,-1]
#     sy = [1,1,-1,-1,1,1,-1,-1]
#     sz = [1,-1,1,-1,1,-1,1,-1]

#     μ⁻ = [0,-1,-1,-1,-1,-1]
#     μ⁺ = [1,0,1,1,1,1]
#     φ⁻ = [0,0,-π/2,π/2,0,π]
#     φ⁺ = [2*π,2*π,π/2,3*π/2,π,2*π]

#     μ_uv(u,v) = sx[u]*(1-(1-v+(sx[u]+1)/2*Nv)*(-v+(sx[u]+1)/2*(Nv+2))/(Nv*(Nv+1)))
#     ϕ_uw(u,v,w) = (π/2)*(w-1)/(-sx[u]*v + (sx[u]+1)/2*(Nv+1)) + π/2 * (2 + (sy[u]+1)/2 - (sz[u]+1)/2 - (sy[u]+1)*(sz[u]+1)/2)

#     for nb in range(1,2*Ndims)

#         Δμ = μ⁺[nb] - μ⁻[nb]
#         Δϕ = φ⁺[nb] - φ⁻[nb]

#         for p in range(1,Np), q in range(1,Np)
#             lp = pl_p[p]
#             lq = pl_q[q]
#             mp = pm_p[p]
#             mq = pm_q[q]
#             for n in range(1,N), m in range(1,N)
#                 for u in range(1,8)
#                     if nb == 1 && sx[u] == -1 continue end
#                     if nb == 2 && sx[u] == 1 continue end
#                     if nb == 3 && sy[u] == -1 continue end
#                     if nb == 4 && sy[u] == 1 continue end
#                     if nb == 5 && sz[u] == -1 continue end
#                     if nb == 6 && sz[u] == 1 continue end
#                     for v in range(1,Nv)
#                         Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
#                         for w in range(1,Nw)
                        
#                             Δμ_uv = μ_uv(u,v+1) - μ_uv(u,v)
#                             Δϕ_uw = ϕ_uw(u,v,w+1) - ϕ_uw(u,v,w)

#                             Mll[p,q,u,v,w,nb] += weight[n] * weight[m] * (Δμ_uv*Δϕ_uw)/4 * real_half_range_spherical_harmonics_up_to_L(L,((μ_uv(u,v)+Δμ_uv/2*(x[n]+1))-μ⁻[nb])/Δμ,(2*π*(ϕ_uw(u,v,w)+Δϕ_uw/2*(x[m]+1))-φ⁻[nb])/Δϕ)[lp+1][lp+mp+1] * real_patch_range_spherical_harmonics_up_to_L(L,(μ_uv(u,v)+Δμ_uv/2*(x[n]+1)),(ϕ_uw(u,v,w)+Δϕ_uw/2*(x[m]+1)),Nv,u,v,w)[lq+1][lq+mq+1]
#                         end 
#                     end
#                 end
#             end
#         end
#     end
#     return Mll
# end