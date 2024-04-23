"""
    flux(cross_sections::Cross_Sections,geometry::Geometry,flux::Flux,type::String)

Calculate the flux solution and print it in the output file.

See also [`transport`](@ref).

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross section informations.
- 'geometry::Geometry': geometry informations.
- 'flux::Flux': flux informations.
- 'type::String': type of output.

# Output Argument(s)
- 'F::Array{Float64}': integrated flux [particle per cm² per s] per group and per voxel.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function flux(cross_sections::Cross_Sections,geometry::Geometry,flux::Flux,type::String)

#----
# Extract geometry data
#----
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()

#----
# Print the flux solution
#----
particles = flux.get_particles()
if type ∉ particles error("Unknown particle type.") end
particles = [type]
Ng = cross_sections.get_number_of_groups(type)
F = zeros(Ng,Ns[1],Ns[2],Ns[3])
for part in particles

    𝚽l = flux.get_flux(part)

    @inbounds for ig in range(1,Ng) ,ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
        F[ig,ix,iy,iz] = 𝚽l[ig,1,1,ix,iy,iz]
    end

end

if Ndims == 1
    return F[:,:,1,1]
elseif Ndims == 2
    return F[:,:,:,1]
elseif Ndims == 3
    return F
end

end