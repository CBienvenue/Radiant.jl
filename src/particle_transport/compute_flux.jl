"""
    compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::Solver,
    source::Source)

Solve the transport equation for a given particle.  

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross section informations.
- `geometry::Geometry` : geometry informations.
- `solver::Solver` : solver informations.
- `source::Source` : source informations.

# Output Argument(s)
- `flux::Flux_Per_Particle`: flux informations.

# Reference(s)
N/A

"""
function compute_flux(cross_sections::Cross_Sections,geometry::Geometry,solver::Solver,source::Source)

#----
# Geometry data
#----
Ndims = geometry.get_dimension()
geo_type = geometry.get_type()
if geo_type != "cartesian" error("Transport of particles in",geo_type," is unavailable.") end
Ns = geometry.get_number_of_voxels()
Δs = geometry.get_voxels_width()
mat = geometry.get_material_per_voxel()

#----
# Preparation of angular discretisation
#----
if solver isa Discrete_Ordinates

    is_SN = true
    is_PN = false
    is_SPH = false

    L = solver.get_legendre_order()
    N = solver.get_quadrature_order()
    quadrature_type = solver.get_quadrature_type()
    SN_type = solver.get_angular_boltzmann()
    Qdims = solver.get_quadrature_dimension(Ndims)

    Ω,w = quadrature(N,quadrature_type,Ndims,Qdims)
    if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
    Nd = length(w)
    Np,Mn,Dn,pl,pm = angular_polynomial_basis(Ω,w,L,SN_type,Qdims)
    Np_surf,Mn_surf,Dn_surf,n⁺_to_n,n_to_n⁺,pl_surf,pm_surf = surface_angular_polynomial_basis(Ω,w,L,SN_type,Qdims,Ndims,geo_type)
    Mll = zeros(0,0)

elseif solver isa Spherical_Harmonics

    is_SN = false
    is_PN = true

    Nd = 0
    Ω = Vector{Vector{Float64}}()
    Mn = zeros(Float64,0,0)
    Dn = zeros(Float64,0,0)
    Mn_surf = Vector{Array{Float64}}()
    Dn_surf = Vector{Array{Float64}}()
    n_to_n⁺ = Vector{Vector{Int64}}()
    L = solver.get_legendre_order()
    polynomial_basis = solver.get_polynomial_basis(Ndims)
    if Ndims == 1
        if polynomial_basis == "legendre"
            is_SPH = false
            Np = L+1
            pl = collect(0:L)
            pm = zeros(Int64,Np)
            Mll = zeros(Np,Np)
            for ip in range(1,Np), jp in range(1,Np)
                il = ip - 1
                jl = jp - 1
                for ik in range(0,div(il,2)), jk in range(0,div(jl,2))
                    for j in range(0,jl-2*jk)
                        Mll[ip,jp] += sqrt(2*jl+1)/2^(il+jl) * (-1)^(ik+jk) * binomial(il,ik) * binomial(jl,jk) * binomial(2*il-2*ik,il) * binomial(2*jl-2*jk,jl) * binomial(jl-2*jk,j) * (-1)^(jl-2*jk-j) * 2^j / (il-2*ik+j+1)
                    end
                end
            end
        elseif polynomial_basis == "spherical-harmonics"
            is_SPH = true
            Np = (L+1)^2
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            p = 1
            for l in range(0,L), m in range(-l,l)
                pl[p] = l
                pm[p] = m
                p += 1
            end
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
        else
            error("Unknown polynomial basis.")
        end
    else
        error("Spherical Harmonics method is only available in 1D.")
    end
    Np_surf = Inf
else
    error("No methods available for $(get_type(solver.particle)) particle.")
end

#----
# Preparation of cross sections
#----

part = solver.get_particle()
solver_type,is_CSD = solver.get_solver_type()
Nmat = cross_sections.get_number_of_materials()
Ng = cross_sections.get_number_of_groups(part)
if is_CSD
    ΔE = cross_sections.get_energy_width(part)
    E = cross_sections.get_energies(part)
    Eb = cross_sections.get_energy_boundaries(part)
end

isFC = solver.get_is_full_coupling()
schemes,𝒪,Nm = solver.get_schemes(geometry,isFC)
ω,𝒞,is_adaptive,𝒲 = scheme_weights(𝒪,schemes,Ndims,is_CSD)

println(">>>Particle: $(get_type(part)) <<<")

# Total cross sections
Σtot = zeros(Ng,Nmat)
if solver_type ∈ [4,5] 
    Σtot = cross_sections.get_absorption(part)
else
    Σtot = cross_sections.get_total(part)
end

# Scattering cross sections
Σs = zeros(Nmat,Ng,Ng,L+1)
if solver_type ∉ [4,5]
    Σs = cross_sections.get_scattering(part,part,L)
end

# Stopping powers
if is_CSD
    S⁻ = zeros(Ng,Nmat); S⁺ = zeros(Ng,Nmat)
    Sb = cross_sections.get_boundary_stopping_powers(part)
    for n in range(1,Nmat)
        S⁻[:,n] = Sb[1:Ng,n] ; S⁺[:,n] = Sb[2:Ng+1,n]
    end
    S = zeros(Ng,Nmat,𝒪[4])
    for n in range(1,Nmat), ig in range(1,Ng)
        S[ig,n,1] = (S⁻[ig,n]+S⁺[ig,n])/2
        if (𝒪[4] > 1) S[ig,n,2] = (S⁻[ig,n]-S⁺[ig,n])/(2*sqrt(3)) end
    end
end

# Momentum transfer
if solver_type ∈ [2,4]
    T = zeros(Ng,Nmat)
    T = cross_sections.get_momentum_transfer(part)
    fokker_planck_type = solver.get_angular_fokker_planck()
    ℳ,λ₀ = fokker_planck_scattering_matrix(N,Nd,quadrature_type,Ndims,fokker_planck_type,Mn,Dn,pl,Np,Qdims)
    Σtot .+= T .* λ₀
end

# Elastic-free approximation
if solver_type == 6
    for n in range(1,Nmat), ig in range(1,Ng)
        Σtot[ig,n] -= Σs[n,ig,ig,1]
    end
end

# External electro-magnetic fields
𝓔 = [0.0,0.0,0.0]; 𝓑 = [0.0,0.0,0.0]
is_EM = false
if is_EM
    q = part.get_charge()
    ℳ_EM = electromagnetic_scattering_matrix(𝓔,𝓑,q,Ω,w,Ndims,Mn,Dn,pl,pm,Np,Ng,Eb,ΔE,Qdims)
else
    ℳ_EM = zeros(Ng,Np,Np);
end

#----
# Acceleration solver
#----

𝒜 = solver.get_acceleration()

#----
# Fixed sources
#----

surface_sources = source.get_surface_sources()
volume_sources = source.get_volume_sources()
Np_surf = Int64(min(Np_surf,length(surface_sources[1,:,1])))

#----
# Flux calculations
#----

ϵ_max = solver.get_convergence_criterion()
I_max = solver.get_maximum_iteration()

# Initialization flux
𝚽l = zeros(Ng,Np,Nm[5],Ns[1],Ns[2],Ns[3])
if is_CSD 𝚽cutoff = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3]) end

# All-group iteration
i_out = 1
is_outer_convergence = false
ϵ_out = Inf
is_outer_iteration = false
Ntot = 0
if is_outer_iteration 𝚽l⁻ = zeros(Ng,Ns[1],Ns[2],Ns[3]) end

while ~(is_outer_convergence)

    ρ_in = -ones(Ng) # In-group spectral radius
    is_CSD ? 𝚽E12 = zeros(Nd,Nm[4],Ns[1],Ns[2],Ns[3]) : 𝚽E12 = Array{Float64}(undef)

    # Loop over energy group
    for ig in range(1,Ng)

        # Calculation of the Legendre components of the source (out-scattering)
        Qlout = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
        if solver_type ∉ [4,5] Qlout = scattering_source(Qlout,𝚽l,Σs[:,:,ig,:],mat,Np,pl,Nm[5],Ns,Ng,ig) end

        # Fixed volumic sources
        Qlout .+= volume_sources[ig,:,:,:,:,:]

        # Calculation of the group flux
        if is_CSD
            if (ig != 1) 𝚽E12 = 𝚽E12 .* ΔE[ig]/ΔE[ig-1] end
            Eg = E[ig]
            ΔEg = ΔE[ig]
            Sg⁻ = S⁻[ig,:]/ΔEg
            Sg⁺ = S⁺[ig,:]/ΔEg
            Sg = S[ig,:,:]/ΔEg
            if solver_type ∈ [2,4]
                Tg = T[ig,:]
            else
                Tg = Vector{Float64}()
                ℳ = Array{Float64}(undef)
            end
        else
            Eg = 0.0
            ΔEg = 0.0
            Sg⁻ = Vector{Float64}()
            Sg⁺ = Vector{Float64}()
            Sg = Vector{Float64}()
            Tg = Vector{Float64}()
            ℳ = Array{Float64}(undef)
        end
        𝚽l[ig,:,:,:,:,:],𝚽E12,ρ_in[ig],Ntot = compute_one_speed(𝚽l[ig,:,:,:,:,:],Qlout,Σtot[ig,:],Σs[:,ig,ig,:],mat,Ndims,Nd,ig,Ns,Δs,Ω,Mn,Dn,Np,pl,pm,Mn_surf,Dn_surf,Np_surf,n_to_n⁺,𝒪,Nm,isFC,𝒞,ω,I_max,ϵ_max,surface_sources[ig,:,:],is_adaptive,is_CSD,solver_type,Eg,ΔEg,𝚽E12,Sg⁻,Sg⁺,Sg,Tg,ℳ,𝒜,Ntot,is_EM,ℳ_EM[ig,:,:],𝒲,Mll,is_SN,is_PN,is_SPH)
    end

    # Verification of convergence in all energy groups
    if is_outer_iteration
        ϵ_out = maximum(vec(abs.(𝚽l[:,1,1,:,:,:] .- 𝚽l⁻)))/maximum(vec(abs.(𝚽l[:,1,1,:,:,:])))
        𝚽l⁻ = 𝚽l[:,1,1,:,:,:]
    end
    if (ϵ_out < ϵ_max || i_out >= I_max) || ~is_outer_iteration
        is_outer_convergence = true
        # Calculate the flux at the cutoff energy
        if is_CSD
            for n in range(1,Nd), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[4]), p in range(1,Np)
                𝚽cutoff[p,is,ix,iy,iz] += Dn[p,n] * 𝚽E12[n,is,ix,iy,iz]
            end
        end
    else
        i_out += 1
    end
    
end

# Save flux
flux = Flux_Per_Particle(part)
flux.add_flux(𝚽l)
if is_CSD flux.add_flux_cutoff(𝚽cutoff) end

return flux

end