"""
    energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,
    solvers::Solvers,sources::Fixed_Sources,flux::Flux,type::String)

Calculate the energy deposition and print it in the output file. 

See also [`charge_deposition`](@ref).

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross section informations.
- 'geometry::Geometry': geometry informations.
- 'solvers::Solvers': solvers informations.
- 'sources::Sources': sources informations.
- 'flux::Flux': flux informations.
- 'type::String': type of output.

# Output Argument(s)
- 'Dtot::Array{Float64}': energy deposition per voxel [MeV/g × cmⁿ per particle], where
   n is the geometry dimension. 

# Author(s)
Charles Bienvenue

# Reference(s)
- Morel (1996) : A Hybrid Multigroup/Continuous-Energy Monte Carlo Method for Solving the
  Boltzmann-Fokker-Planck Equation

"""
function energy_deposition(cross_sections::Cross_Sections,geometry::Geometry,solvers::Solvers,sources::Fixed_Sources,flux::Flux,type::String)

#----
# Extract geometry data
#----
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

#----
# Calculate the total energy deposition
#----
Dtot = zeros(Ns[1],Ns[2],Ns[3])
particles = flux.get_particles()
if type != "total"
    if type ∉ particles error("Unknown particle type.") end
    particles = [type]
end
for part in particles

    # Extract discrete_ordinates data
    discrete_ordinates = solvers.get_method(part)
    _,isCSD = discrete_ordinates.get_solver_type()

    norm = sources.get_normalization_factor()
    D = zeros(Ns[1],Ns[2],Ns[3])

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

    # Extract energy deposition cross sections
    Σe = zeros(Ng,Nmat)
    Σe = cross_sections.get_energy_deposition(part)

    # Extract stopping power at cutoff energy
    if isCSD
        βcutoff = cross_sections.get_stopping_powers(part)[end,:]
        Ecutoff = cross_sections.get_energy_boundaries(part)[end]
    end

    #----
    # Energy deposition calculations
    #----
    @inbounds for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])

        # In-group energy deposition
        for ig in range(1,Ng)
            D[ix,iy,iz] += Σe[ig,mat[ix,iy,iz]] * 𝚽l[ig,1,1,ix,iy,iz]
        end

        # Under-the-cutoff energy deposition
        if isCSD
            D[ix,iy,iz] += Ecutoff * βcutoff[mat[ix,iy,iz]] * 𝚽cutoff[1,1,ix,iy,iz]
        end

        # Normalization
        D[ix,iy,iz] /= ρ[mat[ix,iy,iz]] * norm

    end
    Dtot += D
end

if Ndims == 1
    return Dtot[:,1,1]
elseif Ndims == 2
    return Dtot[:,:,1]
elseif Ndims == 3
    return Dtot
end

end