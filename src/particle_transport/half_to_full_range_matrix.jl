function half_to_full_range_matrix_legendre(L)
    Np = L+1
    Nq = Np
    pl = collect(0:L)
    Mll = zeros(Np,Nq)
    for ip in range(1,Np), jp in range(1,Nq)
        il = pl[ip]
        jl = pl[jp]
        for ik in range(0,div(il,2)), jk in range(0,div(jl,2))
            for j in range(0,jl-2*jk)
                Mll[ip,jp] += sqrt(2*jl+1)/2^(il+jl) * (-1)^(ik+jk) * binomial(il,ik) * binomial(jl,jk) * binomial(2*il-2*ik,il) * binomial(2*jl-2*jk,jl) * binomial(jl-2*jk,j) * (-1)^(jl-2*jk-j) * 2^j / (il-2*ik+j+1)
            end
        end
    end
    return Np,Nq,Mll
end

function half_to_full_range_matrix_spherical_harmonics(L)
    Np = spherical_harmonics_number_basis(L)
    Nq = Np
    pl,pm = spherical_harmonics_indices(L)
    Mll = zeros(Np,Nq)
    for ip in range(1,Np), jp in range(1,Nq)
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
    return Np,Nq,Mll
end

function half_to_full_range_matrix_cartesian_harmonics(L)
    Np = spherical_harmonics_number_basis(L)
    Nq = cartesian_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    ql,qa,qb,qc = cartesian_harmonics_indices(L)
    Mll = zeros(Np,Nq)
    for ip in range(1,Np), jp in range(1,Nq)
        il = pl[ip]
        im = pm[ip]
        jl = ql[jp]
        ja = qa[jp]
        jb = qb[jp]
        jc = qc[jp]
        for i in range(0,div(ja,2)), j in range(0,div(jb,2)), k in range(0,div(jc,2)), r in range(0,div(il-abs(im),2)), t in range(0,abs(im))
            Clm = sqrt((2-(im==0))*factorial_factor([il-abs(im)],[il+abs(im)]))
            α = (-1)^r * factorial_factor([2*il-2*r],[r,il-r,il-abs(im)-2*r])/2^il * binomial(abs(im),t)
            if im ≥ 0
                α *= cos(0.5*π*t)
            else
                α *= sin(0.5*π*t)
            end
            Clabc = cartesian_harmonics_normalization(jl,ja,jb,jc)
            β = (-1)^(i+j+k) * double_factorial(2*jl-2*(i+j+k)-1)/double_factorial(2*jl-1)/2^(i+j+k) * factorial_factor([ja,jb,jc],[i,j,k,ja-2*i,jb-2*j,jc-2*k])
            Mll[ip,jp] += Clm * Clabc * α * β * cartesian_harmonics_I(ja-2*i,il-abs(im)-2*r,jb-2*j+abs(im)-t,jc-2*k+t)
        end
    end
    return Np,Nq,Mll
end