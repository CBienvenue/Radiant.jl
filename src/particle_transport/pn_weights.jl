function pn_weights_legendre_1D(L::Int64)
    Np = L+1
    pl = collect(0:L)
    𝒩 = zeros(Np,Np)
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        lq = pl[q]
        if lp == lq
            𝒩[p,q] = 0.5
        elseif lq == lp-1 && lp > 0
            𝒩[p,q] = 0.5 * lp/(sqrt(2*lp-1)*sqrt(2*lp+1))
        elseif lq == lp+1 && lp < L
            𝒩[p,q] = 0.5 * (lp+1)/(sqrt(2*lp+1)*sqrt(2*lp+3))
        end
    end
    return [𝒩]
end

function pn_weights_spherical_harmonics_1D(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩 = zeros(Np,Np)
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        mp = pm[p]
        lq = pl[q]
        mq = pm[q]
        if (mp == mq)
            if (lq == lp-1) && (lp > 0 && abs(mp) <= lp-1) 𝒩[p,q] = 0.5 * sqrt((lp-mp)*(lp+mp))/(sqrt(2*lp-1)*sqrt(2*lp+1)) end
            if (lq == lp) 𝒩[p,q] = 0.5 end
            if (lq == lp+1) && (lp < L) 𝒩[p,q] = 0.5 * sqrt((lp-mp+1)*(lp+mp+1))/(sqrt(2*lp+1)*sqrt(2*lp+3)) end
        end
    end
    return [𝒩]
end

function pn_weights_spherical_harmonics_2D(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩x = zeros(Np,Np)
    𝒩y = zeros(Np,Np)

    N = 32
    x,w = gauss_legendre(N)

    # Matrix for x-streaming term
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        mp = pm[p]
        lq = pl[q]
        mq = pm[q]
        if (mp == mq)
            if (lq == lp-1) && (lp > 0 && abs(mp) <= lp-1) 𝒩x[p,q] = 0.5 * sqrt((lp-mp)*(lp+mp))/(sqrt(2*lp-1)*sqrt(2*lp+1)) end
            if (lq == lp) 𝒩x[p,q] = 0.5 end
            if (lq == lp+1) && (lp < L) 𝒩x[p,q] = 0.5 * sqrt((lp-mp+1)*(lp+mp+1))/(sqrt(2*lp+1)*sqrt(2*lp+3)) end
        end
    end

    # Matrix for y-streaming term
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        mp = pm[p]
        lq = pl[q]
        mq = pm[q]
        Cp = (-1)^abs(mp) * sqrt((2-(mp==0))/π * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
        Cq = (-1)^abs(mq) * sqrt((2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
        azimutal_integral = 0.0
        for n in range(1,N)
            if mp ≥ 0
                τp = cos(π*x[n]*mp)
            else
                τp = sin(-π*x[n]*mp)
            end
            if mq ≥ 0
                τq = cos(π*x[n]*mq)
            else
                τq = sin(-π*x[n]*mq)
            end
            azimutal_integral += 0.5 * π * w[n] * cos(0.5*π*x[n]) * τp * τq
        end
        cosine_integral = 0.0
        for n in range(1,N)
            fp = 0.0
            for kp in range(0,lp-abs(mp))
                fp += binomial(lp-abs(mp),kp) * factorial_factor([lp+abs(mp)+kp],[abs(mp)+kp]) * ((x[n]-1)/2)^(kp)
            end
            fp *= factorial_factor([lp,lp+abs(mp)],[lp-abs(mp),lp+abs(mp),lp],[(2,-abs(mp))]) * (1-x[n]^2)^(abs(mp)/2)
            fq = 0.0
            for kq in range(0,lq-abs(mq))
                fq += binomial(lq-abs(mq),kq) * factorial_factor([lq+abs(mq)+kq],[abs(mq)+kq]) * ((x[n]-1)/2)^(kq)
            end
            fq *= factorial_factor([lq,lq+abs(mq)],[lq-abs(mq),lq+abs(mq),lq],[(2,-abs(mq))]) * (1-x[n]^2)^(abs(mq)/2)
            cosine_integral += w[n] * 0.25 * sqrt((1-x[n])*(3+x[n])) * fp * fq
        end
        𝒩y[p,q] = Cp * Cq * azimutal_integral * cosine_integral
    end
    return [𝒩x,𝒩y]
end

function pn_weights_spherical_harmonics_3D(L::Int64)
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩x = zeros(Np,Np)
    𝒩y = zeros(Np,Np)
    𝒩z = zeros(Np,Np)
    N = 32
    x,w = gauss_legendre(N)

    # Matrix for x-streaming term
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        mp = pm[p]
        lq = pl[q]
        mq = pm[q]
        if (mp == mq)
            if (lq == lp-1) && (lp > 0 && abs(mp) <= lp-1) 𝒩x[p,q] = 0.5 * sqrt((lp-mp)*(lp+mp))/(sqrt(2*lp-1)*sqrt(2*lp+1)) end
            if (lq == lp) 𝒩x[p,q] = 0.5 end
            if (lq == lp+1) && (lp < L) 𝒩x[p,q] = 0.5 * sqrt((lp-mp+1)*(lp+mp+1))/(sqrt(2*lp+1)*sqrt(2*lp+3)) end
        end
    end

    # Matrix for y- and z-streaming terms
    for p in range(1,Np), q in range(1,Np)
        lp = pl[p]
        mp = pm[p]
        lq = pl[q]
        mq = pm[q]
        Cp = sqrt(2*(2-(mp==0))/π * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
        Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
        azimutal_integral_cos = 0.0
        azimutal_integral_sin = 0.0
        for n in range(1,N)
            if mp ≥ 0
                τp = cos(π*(x[n]+1)*mp)
            else
                τp = sin(-π*(x[n]+1)*mp)
            end
            if mq ≥ 0
                τq = cos(π*(x[n]+1)*mq)
            else
                τq = sin(-π*(x[n]+1)*mq)
            end
            azimutal_integral_cos += 0.25 * π * w[n] * cos(0.25*π*(x[n]+1)) * τp * τq
            azimutal_integral_sin += 0.25 * π * w[n] * sin(0.25*π*(x[n]+1)) * τp * τq
        end
        cosine_integral = 0.0
        for n in range(1,N)
            fp = 0.0
            for kp in range(0,lp-abs(mp))
                fp += binomial(lp-abs(mp),kp) * factorial_factor([lp+abs(mp)+kp],[abs(mp)+kp]) * ((x[n]-1)/2)^(kp)
            end
            fp *= factorial_factor([lp,lp+abs(mp)],[lp-abs(mp),lp+abs(mp),lp],[(2,-abs(mp))]) * (1-x[n]^2)^(abs(mp)/2)
            fq = 0.0
            for kq in range(0,lq-abs(mq))
                fq += binomial(lq-abs(mq),kq) * factorial_factor([lq+abs(mq)+kq],[abs(mq)+kq]) * ((x[n]-1)/2)^(kq)
            end
            fq *= factorial_factor([lq,lq+abs(mq)],[lq-abs(mq),lq+abs(mq),lq],[(2,-abs(mq))]) * (1-x[n]^2)^(abs(mq)/2)
            cosine_integral += w[n] * 0.25 * sqrt((1-x[n])*(3+x[n])) * fp * fq
        end
        𝒩y[p,q] = Cp * Cq * azimutal_integral_cos * cosine_integral
        𝒩z[p,q] = Cp * Cq * azimutal_integral_sin * cosine_integral
    end
    return [𝒩x,𝒩y,𝒩z]
end

function gn_weights_spherical_harmonics(L::Int64,Nv::Int64,Ndims::Int64)
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if Nv <= 0 error("Number of direction cosine patches should be greater than zero.") end
    if ~(1 ≤ Ndims ≤ 3) error("Number of dimensions should be between 1 and 3.") end
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩 = zeros(Np,Np,Ndims,8,Nv,Nv)
    N = 32
    x,weight = gauss_legendre(N)

    # Patch-range basis evaluation on normalized coordinates is independent of (u,v,w):
    # μ' = 0.5*(x[n]+1) ∈ [0,1], ϕ' = (π/4)*(x[m]+1) ∈ [0,π/2]
    # Scaling * Jacobian collapses to a constant factor: (ΔμΔϕ)/4 * (π/(2ΔμΔϕ)) = π/8
    wnm0 = (π/8)

    μ̂ = [0.5*(x[n]+1) for n in 1:N]
    ϕ̂ = [(π/4)*(x[m]+1) for m in 1:N]

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

    # Precompute constants and octant basis tables
    abs_m = zeros(Int64, Np)
    C_oct = zeros(Np)
    for p in 1:Np
        l = pl[p]
        m = pm[p]
        am = abs(m)
        abs_m[p] = am
        C_oct[p] = sqrt(2*(2-(m==0))/π * (2*l+1) * factorial_factor([l-am], [l+am]))
    end

    Plm_oct = zeros(L+1, L+1)
    A_oct = zeros(Np, N)
    for n in 1:N
        _fill_Plm_octant!(Plm_oct, L, μ̂[n])
        for p in 1:Np
            lp = pl[p]
            A_oct[p,n] = C_oct[p] * Plm_oct[lp+1, abs_m[p]+1]
        end
    end
    T_oct = zeros(Np, N)
    for m in 1:N
        ϕ = ϕ̂[m]
        for p in 1:Np
            mp = pm[p]
            T_oct[p,m] = (mp ≥ 0) ? cos(4*mp*ϕ) : sin(4*abs(mp)*ϕ)
        end
    end

    sx_u = (1,1,1,1,-1,-1,-1,-1)
    sy_u = (1,1,-1,-1,1,1,-1,-1)
    sz_u = (1,-1,1,-1,1,-1,1,-1)

    ψvec = zeros(Np)
    μn_vec = zeros(N)
    sμ_vec = zeros(N)
    cosφ = zeros(N)
    sinφ = zeros(N)

    for u in 1:8
        sx = sx_u[u]; sy = sy_u[u]; sz = sz_u[u]
        for v in 1:Nv
            μ_v(v_) = sx*(1-(1-v_+(sx+1)/2*Nv)*(-v_+(sx+1)/2*(Nv+2))/(Nv*(Nv+1)))
            Nw = Int(-sx*v + (sx+1)/2*(Nv+1))
            for w in 1:Nw
                ϕ_w(w_) = (π/2)*(w_-1)/Nw + π/2 * (2 + (sy+1)/2 - (sz+1)/2 - (sy+1)*(sz+1)/2)
                Δμ = μ_v(v+1) - μ_v(v)
                Δϕ = ϕ_w(w+1) - ϕ_w(w)

                μ0 = μ_v(v)
                ϕ0 = ϕ_w(w)
                for n in 1:N
                    μn = μ0 + μ̂[n]*Δμ
                    μn_vec[n] = μn
                    sμ_vec[n] = sqrt(max(0.0, 1-μn^2))
                end
                for m in 1:N
                    φm = ϕ0 + μ̂[m]*Δϕ
                    cosφ[m] = cos(φm)
                    sinφ[m] = sin(φm)
                end

                for n in 1:N
                    wn = weight[n]
                    μn = μn_vec[n]
                    sμ = sμ_vec[n]
                    for m in 1:N
                        wnm = wn * weight[m] * wnm0

                        for p in 1:Np
                            ψvec[p] = A_oct[p,n] * T_oct[p,m]
                        end

                        c1 = wnm * μn
                        if Ndims ≥ 2
                            c2 = wnm * sμ * cosφ[m]
                        end
                        if Ndims == 3
                            c3 = wnm * sμ * sinφ[m]
                        end

                        for p in 1:Np
                            ψp = ψvec[p]
                            for q in 1:Np
                                ψpq = ψp * ψvec[q]
                                𝒩[p,q,1,u,v,w] += c1 * ψpq
                                if Ndims ≥ 2
                                    𝒩[p,q,2,u,v,w] += c2 * ψpq
                                end
                                if Ndims == 3
                                    𝒩[p,q,3,u,v,w] += c3 * ψpq
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return 𝒩
end