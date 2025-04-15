"""
    charge_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,
    sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

Calculate and extract the charge deposition. 

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross section informations.
- `geometry::Geometry` : geometry informations.
- `solvers::Solvers` : solvers informations.
- `sources::Fixed_Sources` : sources informations.
- `flux::Flux` : flux informations.
- `particles::Vector{Particle}` : list of particles

# Output Argument(s)
- `Ctot::Array{Float64}` : charge deposition per voxel [MeV/g × cmⁿ per particle], where
   n is the geometry dimension. 

# Reference(s)
- Morel et al. (1996), A Hybrid Multigroup/Continuous-Energy Monte Carlo Method for Solving
  the Boltzmann-Fokker-Planck Equation.
- Liscum-Powell (2000), Finite Element Numerical Solution of a Self-Adjoint Transport
  Equation for Coupled Electron-Photon Problems

"""
function charge_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,sources::Fixed_Sources,flux::Flux,particles::Vector{<:Particle})

#----
# Extract geometry data
#----
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

#----
# Calculate the total energy deposition
#----
Ctot = zeros(Ns[1],Ns[2],Ns[3])
for part in particles

    # Extract discrete_ordinates data
    discrete_ordinates = solvers.get_method(part)
    _,isCSD = discrete_ordinates.get_solver_type()
    charge = get_charge(part)
    norm = sources.get_normalization_factor()
    C = zeros(Ns[1],Ns[2],Ns[3])

    #----
    # Extract flux data
    #----
    𝚽l = flux.get_flux(part)
    if isCSD 𝚽cutoff = flux.get_flux_cutoff(part) end

    #----
    # Extract cross sections data
    #----
    Nmat = cross_sections.get_number_of_materials()
    Ng = cross_sections.get_number_of_groups(part)
    ρ = cross_sections.get_densities()

    # Extract charge deposition cross sections
    Σc = cross_sections.get_charge_deposition(part)

    #----
    # Charge deposition calculations
    #----
    for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])

        # In-group charge deposition
        for ig in range(1,Ng)
            C[ix,iy,iz] += Σc[ig,mat[ix,iy,iz]] * 𝚽l[ig,1,1,ix,iy,iz]
        end
        if isCSD C[ix,iy,iz] += Σc[Ng+1,mat[ix,iy,iz]] * 𝚽cutoff[1,1,ix,iy,iz] end

        # Normalization
        C[ix,iy,iz] /= ρ[mat[ix,iy,iz]] * norm

    end
    Ctot += C
end

if Ndims == 1
    return Ctot[:,1,1]
elseif Ndims == 2
    return Ctot[:,:,1]
elseif Ndims == 3
    return Ctot
end

end