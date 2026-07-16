"""
    particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,
    geometry::Geometry,solver_in::Solver,
    solver_out::Solver)

Compute the source of particle produced by interaction of another type of particle with
matter, for any pair of solvers.

# Input Argument(s)
- `flux::Flux_Per_Particle`: flux informations of the incoming particle.
- `cross_sections::Cross_Sections`: cross section informations.
- `geometry::Geometry`: geometry informations.
- `solver_in::Solver`: method informations of the incoming particle.
- `solver_out::Solver`: method informations of the outgoing particle.

# Output Argument(s)
- `ps::Source`: source information of the outgoing particle.

# Reference(s)
N/A

"""
function particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,geometry::Geometry,solver_in::Solver,solver_out::Solver)

# Geometry data
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

# Solver data
particle_in = solver_in.get_particle()
particle_out = solver_out.get_particle()
isFC_in = solver_in.get_is_full_coupling()
isFC_out = solver_out.get_is_full_coupling()
_,𝒪_in,Nm_in = solver_in.get_schemes(geometry,isFC_in)
_,𝒪_out,Nm_out = solver_out.get_schemes(geometry,isFC_out)
Nm_in = Nm_in[5]
Nm_out = Nm_out[5]

# Transfer between the two angular discretizations
basis_in = angular_basis(solver_in,geometry)
basis_out = angular_basis(solver_out,geometry)
T = angular_transfer_matrix(basis_in,basis_out)
P_in = basis_in.P
P_out = basis_out.P

# Cross-sections data
Ng_in = cross_sections.get_number_of_groups(particle_in)
Ng_out = cross_sections.get_number_of_groups(particle_out)
Σs = cross_sections.get_scattering(particle_in,particle_out,maximum(basis_in.pl))

# Flux data
𝚽l = flux.get_flux()

# Compute the scattered particle source
Ql_in = zeros(Ng_out,P_in,Nm_in,Ns[1],Ns[2],Ns[3])
Ql_out = zeros(Ng_out,P_out,Nm_out,Ns[1],Ns[2],Ns[3])
particle_sources(Ql_in,𝚽l,Σs,mat,P_in,basis_in.pl,Nm_in,Ns,Ng_in,Ng_out)

# Adapt the source to the new particle flux expansions
map = map_moments(𝒪_in,𝒪_out,isFC_in,isFC_out)
for i in range(1,length(map))
    m = map[i]
    if m != 0
        for p_in in range(1,P_in), p_out in range(1,P_out)
            if T[p_out,p_in] != 0.0
                Ql_out[:,p_out,i,:,:,:] .+= T[p_out,p_in] .* Ql_in[:,p_in,m,:,:,:]
            end
        end
    end
end

# Save source informations
ps = Source(particle_out,cross_sections,geometry,solver_out)
ps.add_volume_source(Ql_out)
return ps

end

"""
    Angular_Basis

Solver-agnostic description of the angular representation in which a solver stores its flux
and expects its source.

# Field(s)
- `is_discrete::Bool` : true if the flux is collocated on a quadrature (`SN`).
- `P::Int64` : number of angular basis functions.
- `pl::Vector{Int64}` : Legendre order associated with each basis function.
- `pm::Vector{Int64}` : spherical-harmonics order associated with each basis function.
- `Qdims::Int64` : dimension of the angular domain.
- `Ω::Vector{Vector{Float64}}` : director cosines (collocated bases only).
- `w::Vector{Float64}` : quadrature weights (collocated bases only).
- `Dn::Array{Float64}` : discrete-to-moment matrix (collocated bases only).
- `type::String` : type of scattering source treatment (collocated bases only).

"""
struct Angular_Basis
    is_discrete ::Bool
    P           ::Int64
    pl          ::Vector{Int64}
    pm          ::Vector{Int64}
    Qdims       ::Int64
    Ω           ::Vector{Vector{Float64}}
    w           ::Vector{Float64}
    Dn          ::Array{Float64}
    type        ::String
end

"""
    angular_basis(solver::Solver,geometry::Geometry)

Describe the angular representation of any solver in the uniform form of an
[`Angular_Basis`](@ref), so that the coupling between two particles is built from the pair of
descriptions.

# Input Argument(s)
- `solver::Solver`: method informations of the particle.
- `geometry::Geometry`: geometry informations.

# Output Argument(s)
- `basis::Angular_Basis`: description of the angular representation.

# Reference(s)
N/A

"""
function angular_basis(solver::Solver,geometry::Geometry)

Ndims = geometry.get_dimension()
L = solver.get_legendre_order()
no_Ω = Vector{Vector{Float64}}()
no_w = Vector{Float64}()
no_Dn = Array{Float64}(undef,0,0)

if solver isa SN

    Qdims = solver.get_quadrature_dimension(Ndims)
    Ω,w = quadrature(solver.get_quadrature_order(),solver.get_quadrature_type(),Ndims,Qdims)
    if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
    type = solver.get_angular_boltzmann()
    P,_,Dn,pl,pm = angular_polynomial_basis(Ω,w,L,type,Qdims)
    return Angular_Basis(true,P,pl,pm,Qdims,Ω,w,Dn,type)

elseif solver isa GN

    polynomial_basis = solver.get_polynomial_basis(Ndims)
    if polynomial_basis == "legendre"
        if Ndims != 1 error("The GN Legendre basis is only available in 1D.") end
        P = L+1
        pl = collect(0:L)
        pm = zeros(Int64,P)
        Qdims = 1
    elseif polynomial_basis == "spherical-harmonics"
        pl,pm = spherical_harmonics_indices(L)
        P = length(pl)
        Qdims = (Ndims == 1) ? 1 : 3
    else
        error("Unknown polynomial basis for the GN solver: " * String(polynomial_basis))
    end
    return Angular_Basis(false,P,pl,pm,Qdims,no_Ω,no_w,no_Dn,"")

elseif solver isa CP

    if Ndims != 1 error("The CP solver is only available in 1D.") end
    P = L+1
    pl = collect(0:L)
    pm = zeros(Int64,P)
    return Angular_Basis(false,P,pl,pm,1,no_Ω,no_w,no_Dn,"")

else
    error("Unknown angular discretization method.")
end

end

"""
    moments_to_directions(pl::Vector{Int64},pm::Vector{Int64},
    Ω::Vector{Vector{Float64}},Qdims::Int64)

Produce the matrix reconstructing the angular flux at the directions `Ω` from its
spherical-harmonics moments, in the convention set by `Qdims`.

# Input Argument(s)
- `pl::Vector{Int64}`: Legendre order associated with each moment.
- `pm::Vector{Int64}`: spherical-harmonics order associated with each moment.
- `Ω::Vector{Vector{Float64}}`: director cosines at which the flux is reconstructed.
- `Qdims::Int64`: dimension of the angular domain of the reconstructed flux.

# Output Argument(s)
- `M::Array{Float64}`: moment-to-discrete matrix.

# Reference(s)
- Lewis (1984), Computational Methods of Neutron Transport.

"""
function moments_to_directions(pl::Vector{Int64},pm::Vector{Int64},Ω::Vector{Vector{Float64}},Qdims::Int64)

P = length(pl)
Nd = length(Ω[1])
Lmax = maximum(pl)
M = zeros(Nd,P)

if Qdims == 1
    for n in range(1,Nd)
        Pl = legendre_polynomials_up_to_L(Lmax,Ω[1][n])
        for p in range(1,P)
            if pm[p] == 0 M[n,p] = (2*pl[p]+1)/2 * Pl[pl[p]+1] end
        end
    end
else
    for n in range(1,Nd)
        Rlm = real_spherical_harmonics_up_to_L(Lmax,Ω[1][n],atan(Ω[3][n],Ω[2][n]))
        for p in range(1,P)
            M[n,p] = (2*pl[p]+1)/(4*π) * Rlm[pl[p]+1][pl[p]+pm[p]+1]
        end
    end
end
return M

end

"""
    angular_transfer_matrix(basis_in::Angular_Basis,basis_out::Angular_Basis)

Produce the matrix transferring angular moments from the basis of an incoming particle to the
basis of an outgoing one, for any pair of solvers.

# Input Argument(s)
- `basis_in::Angular_Basis`: angular basis of the incoming particle.
- `basis_out::Angular_Basis`: angular basis of the outgoing particle.

# Output Argument(s)
- `T::Array{Float64}`: transfer matrix.

# Reference(s)
N/A

"""
function angular_transfer_matrix(basis_in::Angular_Basis,basis_out::Angular_Basis)

is_same_quadrature = basis_in.is_discrete && basis_out.is_discrete && basis_in.type == basis_out.type && basis_in.Qdims == basis_out.Qdims && length(basis_in.w) == length(basis_out.w) && basis_in.w == basis_out.w && basis_in.Ω == basis_out.Ω

if basis_out.is_discrete && ~is_same_quadrature
    return basis_out.Dn * moments_to_directions(basis_in.pl,basis_in.pm,basis_out.Ω,basis_out.Qdims)
else
    T = zeros(basis_out.P,basis_in.P)
    index_in = Dict{Tuple{Int64,Int64},Int64}()
    for p in range(1,basis_in.P)
        index_in[(basis_in.pl[p],basis_in.pm[p])] = p
    end
    for p_out in range(1,basis_out.P)
        p_in = get(index_in,(basis_out.pl[p_out],basis_out.pm[p_out]),0)
        if p_in != 0 T[p_out,p_in] = 1.0 end
    end
    return T
end

end