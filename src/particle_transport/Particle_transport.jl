# Transport calculations

# Geometry
include("geometry.jl")

# Sources"
include("volume_source.jl")
include("surface_source.jl")
include("scattering_source.jl")
include("particle_source.jl")
include("fokker_planck_source.jl")
include("fokker_planck_scattering_matrix.jl")
include("fokker_planck_finite_difference.jl")
include("fokker_planck_differential_quadrature.jl")
include("fokker_planck_galerkin.jl")

# Transport
include("transport.jl")
include("compute_flux.jl")
include("compute_one_speed.jl")
include("compute_sweep_1D.jl")
include("compute_sweep_2D.jl")
include("compute_sweep_3D.jl")
include("flux_1D_BTE.jl")
include("flux_1D_BFP.jl")
include("flux_2D_BTE.jl")
include("flux_2D_BFP.jl")
include("flux_3D_BTE.jl")
include("flux_3D_BFP.jl")
include("adaptive_1D.jl")
include("adaptive_2D.jl")
include("scheme_weights.jl")
include("map_moments.jl")
include("constant_linear.jl")
include("angular_polynomial_basis.jl")
include("livolant.jl")

# CUDA transport
#include("compute_sweep_2D_cuda.jl")
#include("compute_sweep_3D_cuda.jl")
#include("flux_2D_BTE_cuda.jl")
#include("flux_2D_BFP_cuda.jl")
#include("flux_3D_BFP_cuda.jl")

# Print
include("print_dose.jl")
include("print_charge.jl")
include("print_flux.jl")