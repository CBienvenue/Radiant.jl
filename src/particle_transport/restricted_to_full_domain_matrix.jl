function half_to_full_range_matrix_legendre(L)
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

function half_to_full_range_matrix_spherical_harmonics(L)
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

function quarter_to_full_range_matrix_spherical_harmonics(L)
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

function octant_to_full_range_matrix_spherical_harmonics(L)
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