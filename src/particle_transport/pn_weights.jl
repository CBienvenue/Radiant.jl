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
