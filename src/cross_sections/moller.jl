"""
    moller(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0,Ti::Float64=0.0,
    is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)

Gives the Møller differential cross-sections.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `W::Float64` : energy lost by the incoming electron.
- `Ui::Float64` : binding energy of the electron(s).
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.
- `is_focusing_møller::Bool` : boolean to enable or not the focusing term.
- `is_hydrogenic_distribution_term::Bool`: boolean to enable or not the hydrogenic
  distribution term.

# Output Argument(s)
- `σs::Float64` : Møller differential cross-sections.

# Reference(s)
- Møller (1932), Zur theorie des durchgangs schneller elektronen durch materie.
- Perkins (1989), The Livermore electron impact ionization data base.
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function moller(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0,Ti::Float64=0.0,is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    if (is_focusing_møller) P = moller_focusing_term(Ei,Ui,Ti) else P = 1 end
    if (is_hydrogenic_distribution_term) G = moller_hydrogenic_distribution_term(Ei,W,Ui,Ti) else G = 0 end
    
    # Møller formula
    F = 1/(W+Ui)^2 + 1/(Ei-W)^2 + 1/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * 1/((Ei-W)*(W+Ui)) + G
    σs = 2*π*rₑ^2/β² * Zi * P * F

    return σs
end


"""
    integrate_moller(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=1e-5,Wmax::Float64=Ei,
    is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)

Gives the integration of the Møller differential cross-sections multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Wmin::Float64` : minimum energy lost by the incoming electron.
- `Wmax::Float64` : maximum energy lost by the incoming electron.
- `is_focusing_møller::Bool` : boolean to enable or not the focusing term.
- `is_hydrogenic_distribution_term::Bool`: boolean to enable or not the hydrogenic
  distribution term.

# Output Argument(s)
- `σn::Float64` : integrated Møller cross-sections.

# Reference(s)
- Møller (1932), Zur theorie des durchgangs schneller elektronen durch materie.
- Perkins (1989), The Livermore electron impact ionization data base.
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function integrate_moller(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=0.0,Wmax::Float64=Ei/2,is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)
    Nshells,Zi,Ui,Ti,_,_ = electron_subshells(Z)
    σn = 0
    for δi in range(1,Nshells)
        σn += Zi[δi] * integrate_moller_per_subshell(Ei,n,Ui[δi],Ti[δi],Wmin,Wmax,is_focusing_møller,is_hydrogenic_distribution_term)
    end
    return σn
end

"""
    integrate_moller_per_subshell(Ei::Float64,n::Int64,Ui::Float64=0.0,Ti::Float64=0.0,
    Wmin::Float64=Ui,Wmax::Float64=Ei/2,is_focusing_møller::Bool=false)

Gives the integration of the Møller differential cross-sections for an electron with a given
binding energy, multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Ui::Float64` : binding energy of the subshell.
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.
- `Wmin::Float64` : minimum energy lost by the incoming electron.
- `Wmax::Float64` : maximum energy lost by the incoming electron.
- `is_focusing_møller::Bool` : boolean to enable or not the focusing term.
- `is_hydrogenic_distribution_term::Bool`: boolean to enable or not the hydrogenic
  distribution term.

# Output Argument(s)
- `σni::Float64` : integrated Møller cross-sections in a given subshell.

# Reference(s)
- Møller (1932), Zur theorie des durchgangs schneller elektronen durch materie.
- Perkins (1989), The Livermore electron impact ionization data base.
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function integrate_moller_per_subshell(Ei::Float64,n::Int64,Ui::Float64=0.0,Ti::Float64=0.0,Wmin::Float64=Ui,Wmax::Float64=Ei/2,is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    if (is_focusing_møller) P = moller_focusing_term(Ei,Ui,Ti) else P = 1 end

    # Compute the Møller cross-section per subshell
    σni = 0
    W⁺ = min((Ei-Ui)/2,Wmax)
    W⁻ = max(Ui,Wmin)
    if W⁺ > W⁻
        if (is_hydrogenic_distribution_term) G = integrate_moller_hydrogenic_distribution_term(n,Ei,Ui,Ti,W⁻,W⁺) else G = 0 end
        if n == 0
            J₀(W) = -1/(W+Ui) + 1/(Ei-W) + W/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(W+Ui)-log(Ei-W))/(Ei+Ui)
            σni += (J₀(W⁺) - J₀(W⁻)) + G
        elseif n == 1
            J₁(W) = log(W+Ui) + log(Ei-W) + (Ei+Ui)/(Ei-W) + W*(W+2*Ui)/(2*(Ei+1)^2) + (2*Ei+1)/(Ei+1)^2 * log(Ei-W)
            σni += (J₁(W⁺) - J₁(W⁻)) + G
        else
            error("Integral is given only for n=0 or n=1.")
        end
    end
    σni *= 2*π*rₑ^2/β² * P
    return σni
end

"""
    angular_moller(Ei::Float64,Ef::Float64,L::Int64)

Gives the Legendre moments of the Møller angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing particle energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `σn::Float64` : Legendre moments of the Møller angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.

"""
function angular_moller(Ei::Float64,Ef::Float64,L::Int64)

    # Initialization
    Wℓ = zeros(L+1)

    # Compute the direction cosine
    μ = sqrt((Ef*(Ei+2))/(Ei*(Ef+2)))

    # Compute the Legendre moments angular distribution
    Pℓμ = legendre_polynomials(L,μ)
    for ℓ in range(0,L) Wℓ[ℓ+1] += Pℓμ[ℓ+1] end
    return Wℓ
end

"""
    moller_focusing_term(Ei::Float64,Ui::Float64,Ti::Float64)

Gives the Møller focusing term.

# Input Argument(s)
- `Ei::Float64` : incoming electron energy.
- `Ui::Vector{Float64}`: binding energy per subshell.
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.

# Output Argument(s)
- `P::Float64` : Møller focusing term.

# Reference(s)
- Seltzer (1988), Cross Sections for Bremsstrahlung Production and Electron-Impact
  Ionization.
- Perkins (1989), The Livermore electron impact ionization data base.

"""
function moller_focusing_term(Ei::Float64,Ui::Float64,Ti::Float64)
    return Ei/(Ei+Ui+Ti)
end

"""
    moller_hydrogenic_distribution_term(Ei::Float64,W::Float64,Ui::Float64,Ti::Float64)

Gives the Møller term accounting for the isotropic hydrogenic distribution of orbital
electron velocities.

# Input Argument(s)
- `Ei::Float64` : incoming electron energy.
- `W::Float64` : energy lost by the incoming electron.
- `Ui::Vector{Float64}`: binding energy per subshell.
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.

# Output Argument(s)
- `Gi::Float64` : Møller term accounting for the isotropic hydrogenic distribution of orbital
  electron velocities.

# Reference(s)
- Seltzer (1988), Cross Sections for Bremsstrahlung Production and Electron-Impact
  Ionization.
- Perkins (1989), The Livermore electron impact ionization data base.

"""
function moller_hydrogenic_distribution_term(Ei::Float64,W::Float64,Ui::Float64,Ti::Float64)
	y = W/Ti
	Gi = 8*Ti/(3*π) * (1/(W+Ui)^3 + 1/(Ei-W)^3) * (atan(sqrt(y)) + sqrt(y)*(y-1)/(y+1)^2)
	return Gi
end


"""
    integrate_moller_hydrogenic_distribution_term(n::Int64,Ei::Float64,Ui::Float64,Ti::Float64,Wmin::Float64=0.0,Wmax::Float64=(Ei-Ui)/2)

Gives the integrated Møller term accounting for the isotropic hydrogenic distribution of orbital
electron velocities, multiplied by Wⁿ.

# Input Argument(s)
- `Ei::Float64` : incoming electron energy.
- `W::Float64` : energy lost by the incoming electron.
- `Ui::Vector{Float64}`: binding energy per subshell.
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.

# Output Argument(s)
- `Gi::Float64` : integrated Møller term accounting for the isotropic hydrogenic distribution of orbital
  electron velocities.

# Reference(s)
- Seltzer (1988), Cross Sections for Bremsstrahlung Production and Electron-Impact
  Ionization.
- Perkins (1989), The Livermore electron impact ionization data base.

"""
function integrate_moller_hydrogenic_distribution_term(n::Int64,Ei::Float64,Ui::Float64,Ti::Float64,Wmin::Float64=0.0,Wmax::Float64=(Ei-Ui)/2)
	Gi = 0
	W⁺ = min((Ei-Ui)/2,Wmax)
    W⁻ = max(Ui,Wmin)
    if W⁺ > W⁻
		if n == 0
            if abs(Ui-Ti) ≤ 1e-7
                G₀_hydrogen(W) = (((-5*Ui^2+6*W*Ui+3*W^2)*atan(sqrt(W/Ui))+5*sqrt(W/Ui)*(Ui+3/5*W)*Ui)/(W+Ui)^2/Ui^2/16)+(-(Ei*Ui)^(-1/2)*(3*(Ei+Ui/3)*Ui*(Ei-W)^2*atanh(sqrt(W/Ui)*Ui*(Ei*Ui)^(-1/2))+(-4*Ei*(Ui/2+Ei-W/2)*(W+Ui)*atan(sqrt(W/Ui))+Ui*sqrt(W/Ui)*(Ei-W)*(Ui+Ei))*sqrt(Ei*Ui))/(Ei-W)^2/(Ui+Ei)^2/Ei/4)+(-(Ei*Ui)^(-1/2)*(3*(Ei+Ui/3)*Ui*(Ei-W)^2*atanh(sqrt(W/Ui)*Ui*(Ei*Ui)^(-1/2))+(-4*Ei*(Ui/2+Ei-W/2)*(W+Ui)*atan(sqrt(W/Ui))+Ui*sqrt(W/Ui)*(Ei-W)*(Ui+Ei))*sqrt(Ei*Ui))/(Ei-W)^2/(Ui+Ei)^2/Ei/4)+(3/4*((W+Ui)*(Ei-W)^2*(Ei^3-11*Ui*Ei^2+13/3*Ei*Ui^2+Ui^3/3)*atanh(sqrt(W/Ui)*Ui*(Ui*Ei)^(-1/2))+13/3*sqrt(Ui*Ei)*(-16/13*Ei*(W+Ui)*(Ei-W)^2*(-2*Ui+Ei)*atan(sqrt(W/Ui))+(Ui+Ei)*((Ui+5/13*W)*Ei^3+(-12/13*Ui^2-31/13*W*Ui-3/13*W^2)*Ei^2-Ui*(Ui^2-11*W*Ui-20*W^2)*Ei/13-Ui^2*W*(W+Ui)/13)*sqrt(W/Ui)))*(Ui*Ei)^(-1/2)*Ui/(Ui+Ei)^4/(W+Ui)/(Ei-W)^2/Ei)
                Gi = G₀_hydrogen(Wmax) - G₀_hydrogen(Wmin)
            else
			    G₀(W) = ((Ti*(W+Ui)^2*(-3*Ui+Ti)*atan(sqrt(W/Ti)*Ti*(Ti*Ui)^(-1/2))+sqrt(Ti*Ui)*(-2*Ui*(W+Ti)*(Ti-2*Ui-W)*atan(sqrt(W/Ti))+Ti*sqrt(W/Ti)*(W+Ui)*(Ti-Ui)))*(Ti*Ui)^(-1/2)/(W+Ui)^2/(Ti-Ui)^2/Ui/4)+(-(Ei*Ti)^(-1/2)*(3*Ti*(Ei+Ti/3)*(Ei-W)^2*atanh(sqrt(W/Ti)*Ti*(Ei*Ti)^(-1/2))+(-4*(Ti/2+Ei-W/2)*(W+Ti)*Ei*atan(sqrt(W/Ti))+Ti*sqrt(W/Ti)*(Ei-W)*(Ti+Ei))*sqrt(Ei*Ti))/(Ei-W)^2/(Ti+Ei)^2/Ei/4)+(-Ti*(Ti*Ui)^(-1/2)*((W+Ui)^2*(W+Ti)*(Ti^3-13*Ui*Ti^2-33*Ti*Ui^2-3*Ui^3)*atan(sqrt(W/Ti)*Ti*(Ti*Ui)^(-1/2))-(-32*Ui*(Ti+Ui/2)*(W+Ui)^2*(W+Ti)*atan(sqrt(W/Ti))+((-13*Ti-5*W)*Ui^3+(-12*Ti^2-31*W*Ti-3*W^2)*Ui^2+Ti*(Ti^2-11*W*Ti-20*W^2)*Ui-Ti^2*W*(W+Ti))*sqrt(W/Ti)*(Ti-Ui))*sqrt(Ti*Ui))/(Ti-Ui)^4/(W+Ui)^2/Ui/(W+Ti)/4)+(3/4*Ti*((Ei^3-11*Ti*Ei^2+13/3*Ei*Ti^2+Ti^3/3)*(W+Ti)*(Ei-W)^2*atanh(sqrt(W/Ti)*Ti*(Ti*Ei)^(-1/2))+13/3*(-16/13*Ei*(W+Ti)*(Ei-W)^2*(-2*Ti+Ei)*atan(sqrt(W/Ti))+sqrt(W/Ti)*(Ti+Ei)*((Ti+5/13*W)*Ei^3+(-12/13*Ti^2-31/13*W*Ti-3/13*W^2)*Ei^2-Ti*(Ti^2-11*W*Ti-20*W^2)*Ei/13-Ti^2*W*(W+Ti)/13))*sqrt(Ti*Ei))*(Ti*Ei)^(-1/2)/(Ti+Ei)^4/(Ei-W)^2/Ei/(W+Ti))
			    Gi = G₀(Wmax) - G₀(Wmin)
            end
		elseif n == 1
            if abs(Ui-Ti) ≤ 1e-7
                G₁_hydrogen(W) = (((-3*Ui^2-6*W*Ui+5*W^2)*atan(sqrt(W/Ui))+3*(Ui+5/3*W)*sqrt(W/Ui)*Ui)/(W+Ui)^2/Ui/16)+(-(Ei*Ui)^(-1/2)*(-Ui*(Ei-W)^2*(3*Ui+Ei)*atanh(sqrt(W/Ui)*Ui*(Ei*Ui)^(-1/2))+sqrt(Ei*Ui)*(2*((Ei-2*W)*Ui-Ei*W)*(W+Ui)*atan(sqrt(W/Ui))+Ui*sqrt(W/Ui)*(Ei-W)*(Ui+Ei)))/(Ei-W)^2/(Ui+Ei)^2/4)+((3*(W+Ui)^4*atan(sqrt(W/Ui))-3*sqrt(W/Ui)*(Ui^3+11/3*Ui^2*W+53/3*W^2*Ui-W^3)*Ui)/(W+Ui)^4/Ui/96)+(-((W+Ui)*(Ei-W)^2*(Ei^3+13*Ei^2*Ui-33*Ei*Ui^2+3*Ui^3)*atanh(sqrt(W/Ui)*Ui*(Ui*Ei)^(-1/2))-sqrt(Ui*Ei)*(32*(W+Ui)*(Ei-Ui/2)*Ui*(Ei-W)^2*atan(sqrt(W/Ui))+((3*Ei-5*W)*Ui^3+(-20*Ei^2+31*Ei*W-13*W^2)*Ui^2+Ei*(Ei^2-11*Ei*W+12*W^2)*Ui+Ei^2*W*(Ei+W))*sqrt(W/Ui)*(Ui+Ei)))*(Ui*Ei)^(-1/2)*Ui/(Ui+Ei)^4/(W+Ui)/(Ei-W)^2/4)
                Gi = G₁_hydrogen(Wmax) - G₁_hydrogen(Wmin)
            else
                G₁(W) = (-(-3*Ti*(Ti-Ui/3)*(W+Ui)^2*atan(sqrt(W/Ti)*Ti*(Ui*Ti)^(-1/2))+sqrt(Ui*Ti)*(2*((Ui+2*W)*Ti-W*Ui)*(W+Ti)*atan(sqrt(W/Ti))+Ti*sqrt(W/Ti)*(W+Ui)*(Ti-Ui)))*(Ui*Ti)^(-1/2)/(W+Ui)^2/(Ti-Ui)^2/4)+(-(-Ti*(Ei-W)^2*(3*Ti+Ei)*atanh(sqrt(W/Ti)*Ti*(Ei*Ti)^(-1/2))+(2*(W+Ti)*((Ei-2*W)*Ti-Ei*W)*atan(sqrt(W/Ti))+Ti*sqrt(W/Ti)*(Ei-W)*(Ti+Ei))*sqrt(Ei*Ti))*(Ei*Ti)^(-1/2)/(Ei-W)^2/(Ti+Ei)^2/4)+(-3/4*Ti*((Ti^3+11*Ui*Ti^2+13/3*Ui^2*Ti-Ui^3/3)*(W+Ui)^2*(W+Ti)*atan(sqrt(W/Ti)*Ti*(Ui*Ti)^(-1/2))-(16/3*Ti*(W+Ui)^2*(W+Ti)*(Ti+2*Ui)*atan(sqrt(W/Ti))+sqrt(W/Ti)*((Ui+5/3*W)*Ti^3+(20/3*Ui^2+31/3*W*Ui+13/3*W^2)*Ti^2+Ui*(Ui^2+11*W*Ui+12*W^2)*Ti/3+Ui^2*W*(-W+Ui)/3)*(Ti-Ui))*sqrt(Ui*Ti))*(Ui*Ti)^(-1/2)/(Ti-Ui)^4/(W+Ti)/(W+Ui)^2
                )+(-(Ei*Ti)^(-1/2)*Ti*((W+Ti)*(Ei-W)^2*(Ei^3+13*Ei^2*Ti-33*Ei*Ti^2+3*Ti^3)*atanh(sqrt(W/Ti)*Ti*(Ei*Ti)^(-1/2))-sqrt(Ei*Ti)*(32*(W+Ti)*(Ei-W)^2*Ti*(Ei-Ti/2)*atan(sqrt(W/Ti))+sqrt(W/Ti)*((3*Ei-5*W)*Ti^3+(-20*Ei^2+31*Ei*W-13*W^2)*Ti^2+Ei*(Ei^2-11*Ei*W+12*W^2)*Ti+Ei^2*W*(Ei+W))*(Ti+Ei)))/(Ti+Ei)^4/(W+Ti)/(Ei-W)^2/4)
                Gi = G₁(Wmax) - G₁(Wmin)
            end
		else
			error("Integral is given only for n=0 or n=1.")
		end
		Gi *= 8*Ti/(3*π)
	end
	if isinf(Gi) || isnan(Gi) error([Gi,n,Ti,Ui]) end
	return Gi
end
