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
    F = 1/W^2 + 1/(Ei-W+Ui)^2 + 1/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * 1/(W*(Ei-W+Ui)) + G
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
- `E⁻min::Float64` : minimum energy of the knock-on electron.
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
function integrate_moller(Z::Int64,Ei::Float64,n::Int64,E⁻min::Float64=0.0,is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)
    Nshells,Zi,Ui,Ti,_,_ = electron_subshells(Z)
    σn = 0
    for δi in range(1,Nshells)
        σn += Zi[δi] * integrate_moller_per_subshell(Ei,n,Ui[δi],Ti[δi],E⁻min,is_focusing_møller,is_hydrogenic_distribution_term)
    end
    return σn
end

"""
    integrate_moller_per_subshell(Ei::Float64,n::Int64,Ui::Float64=0.0,Ti::Float64=0.0,
    Wmin::Float64=Ui,Wmax::Float64=Ei/2,is_focusing_møller::Bool=false)

Gives the integration of the Møller differential cross-sections for an electron with a given
binding energy, multiplied by (W+Ui)ⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Ui::Float64` : binding energy of the subshell.
- `Ti::Vector{Float64}`: mean kinetic energy per subshell.
- `E⁻min::Float64` : minimum energy of the knock-on electron.
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
function integrate_moller_per_subshell(Ei::Float64,n::Int64,Ui::Float64=0.0,Ti::Float64=0.0,E⁻min::Float64=0.0,is_focusing_møller::Bool=false,is_hydrogenic_distribution_term::Bool=false)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    if (is_focusing_møller) P = moller_focusing_term(Ei,Ui,Ti) else P = 1 end

    # Compute the Møller cross-section per subshell
    σni = 0
    W⁺ = (Ei+Ui)/2
    W⁻ = E⁻min+Ui
    if W⁺ > W⁻
        if (is_hydrogenic_distribution_term) G = integrate_moller_hydrogenic_distribution_term(n,Ei,Ui,Ti,E⁻min) else G = 0 end
        if n == 0
            J₀(W) = -1/W + 1/(Ei-W+Ui) +  W/(Ei+1)^2 + (2*Ei+1)/((Ei+1)^2*(Ei+Ui))*(log(Ei-W+Ui)-log(W))
            σni += (J₀(W⁺) - J₀(W⁻)) + G
        elseif n == 1
            J₁(W) = log(W) + log(Ei-W+Ui) + (Ei+Ui)/(Ei-W+Ui) + W^2/(2*(Ei+1)^2) + (2*Ei+1)/(Ei+1)^2 * log(Ei-W+Ui)
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
	y = (W-Ui)/Ti
	Gi = 8*Ti/(3*π) * (1/W^3 + 1/(Ei+Ui-W)^3) * (atan(sqrt(y)) + sqrt(y)*(y-1)/(y+1)^2)
	return Gi
end


"""
    integrate_moller_hydrogenic_distribution_term(n::Int64,Ei::Float64,Ui::Float64,
    Ti::Float64,E⁻min::Float64=0.0)

Gives the integrated Møller term accounting for the isotropic hydrogenic distribution of 
orbital electron velocities, multiplied by Wⁿ.

# Input Argument(s)
- `Ei::Float64` : incoming electron energy.
- `W::Float64` : energy lost by the incoming electron.
- `Ui::Float64`: binding energy per subshell.
- `Ti::Float64`: mean kinetic energy per subshell.
- `E⁻min::Float64` : minimum energy of the knock-on electron.

# Output Argument(s)
- `Gi::Float64` : integrated Møller term accounting for the isotropic hydrogenic
  distribution of orbital electron velocities.

# Reference(s)
- Seltzer (1988), Cross Sections for Bremsstrahlung Production and Electron-Impact
  Ionization.
- Perkins (1989), The Livermore electron impact ionization data base.

"""
function integrate_moller_hydrogenic_distribution_term(n::Int64,Ei::Float64,Ui::Float64,Ti::Float64,E⁻min::Float64=0.0)
	Gi = 0
	W⁺ = (Ei+Ui)/2
    W⁻ = E⁻min+Ui
    if W⁺ > W⁻
		if n == 0
            if abs(Ui-Ti) ≤ 1e-7
                G₀_hydrogen(W) = ((3*sqrt((W-Ui)/Ui)*W*Ui+2*sqrt((W-Ui)/Ui)*Ui^2+3*atan(sqrt((W-Ui)/Ui))*W^2-8*atan(sqrt((W-Ui)/Ui))*Ui^2)/W^2/Ui^2/16)+((-3/4*Ui*(Ei+Ui/3)*(Ei-W+Ui)^2*atanh(sqrt((W-Ui)/Ui)*Ui*(Ei*Ui)^(-1/2))+(Ei*W*(Ei-W/2+Ui)*atan(sqrt((W-Ui)/Ui))-sqrt((W-Ui)/Ui)*Ui*(Ui+Ei)*(Ei-W+Ui)/4)*sqrt(Ei*Ui))*(Ei*Ui)^(-1/2)/(Ei-W+Ui)^2/(Ui+Ei)^2/Ei)+((-3*atan(sqrt((W-Ui)/Ui))*W^4-3*Ui*sqrt((W-Ui)/Ui)*(W^3+2/3*W^2*Ui+40/3*W*Ui^2-16*Ui^3))/W^4/Ui^2/96)+(3/4*Ui*((Ei-W+Ui)^2*(Ei^3-11*Ui*Ei^2+13/3*Ui^2*Ei+Ui^3/3)*W*atanh(sqrt((W-Ui)/Ui)*Ui*(Ui*Ei)^(-1/2))+5/3*sqrt(Ui*Ei)*(-16/5*Ei*W*(-2*Ui+Ei)*(Ei-W+Ui)^2*atan(sqrt((W-Ui)/Ui))+sqrt((W-Ui)/Ui)*(Ei+Ui)*((W+8/5*Ui)*Ei^3+(-3/5*W^2-5*Ui*W+16/5*Ui^2)*Ei^2+(4*W^2*Ui-29/5*W*Ui^2+8/5*Ui^3)*Ei-W*Ui^2*(W-Ui)/5)))*(Ui*Ei)^(-1/2)/(Ei+Ui)^4/(Ei-W+Ui)^2/Ei/W)
                Gi = G₀_hydrogen(W⁺) - G₀_hydrogen(W⁻)
            else
                G₀(W) = (-(Ti*Ui)^(-1/2)*(-Ti*W^2*(-3*Ui+Ti)*atan(sqrt((W-Ui)/Ti)*Ti*(Ti*Ui)^(-1/2))/2+(Ui*(W-Ui+Ti)*(Ti-Ui-W)*atan(sqrt((W-Ui)/Ti))-Ti*sqrt((W-Ui)/Ti)*W*(Ti-Ui)/2)*sqrt(Ti*Ui))/W^2/(Ti-Ui)^2/Ui/2)+(-(3*(Ei-W+Ui)^2*Ti*(Ei+Ti/3)*atanh(sqrt((W-Ui)/Ti)*Ti*(Ei*Ti)^(-1/2))+(-4*(Ti/2+Ei-W/2+Ui/2)*(W-Ui+Ti)*Ei*atan(sqrt((W-Ui)/Ti))+Ti*sqrt((W-Ui)/Ti)*(Ei-W+Ui)*(Ti+Ei))*sqrt(Ei*Ti))*(Ei*Ti)^(-1/2)/(Ei-W+Ui)^2/(Ti+Ei)^2/Ei/4)+((Ui*Ti)^(-1/2)*(-W^2*(Ti^3-13*Ti^2*Ui-33*Ti*Ui^2-3*Ui^3)*(W-Ui+Ti)*atan(sqrt(-(-W+Ui)/Ti)*Ti*(Ui*Ti)^(-1/2))/2+(-16*W^2*(Ti+Ui/2)*(W-Ui+Ti)*Ui*atan(sqrt(-(-W+Ui)/Ti))+(-Ui+Ti)*(Ui^4+(-Ti+W/2)*Ui^3+(-Ti^2+9/2*W*Ti-3/2*W^2)*Ui^2+Ti*(Ti^2-9/2*W*Ti-10*W^2)*Ui-Ti^2*W*(W+Ti)/2)*sqrt(-(-W+Ui)/Ti))*sqrt(Ui*Ti))*Ti/(-Ui+Ti)^4/W^2/Ui/(W-Ui+Ti)/2)+(3/4*((Ei-W+Ui)^2*(W-Ui+Ti)*(Ei^3-11*Ti*Ei^2+13/3*Ei*Ti^2+Ti^3/3)*atanh(sqrt(-(-W+Ui)/Ti)*Ti*(Ti*Ei)^(-1/2))+13/3*(-16/13*Ei*(W-Ui+Ti)*(Ei-W+Ui)^2*(-2*Ti+Ei)*atan(sqrt(-(-W+Ui)/Ti))+((Ti-5/13*Ui+5/13*W)*Ei^3+(-12/13*Ti^2+(31/13*Ui-31/13*W)*Ti-3/13*(-W+Ui)^2)*Ei^2-(Ti^2+(11*Ui-11*W)*Ti-20*(-W+Ui)^2)*Ti*Ei/13+Ti^2*(-W+Ui)*(W-Ui+Ti)/13)*(Ti+Ei)*sqrt(-(-W+Ui)/Ti))*sqrt(Ti*Ei))*Ti*(Ti*Ei)^(-1/2)/(Ti+Ei)^4/(Ei-W+Ui)^2/Ei/(W-Ui+Ti))
                Gi = G₀(W⁺) - G₀(W⁻)
            end
		elseif n == 1
            if abs(Ui-Ti) ≤ 1e-7
                G₁_hydrogen(W) = (((W-2*Ui)*atan(sqrt((W-Ui)/Ui))+sqrt((W-Ui)/Ui)*Ui)/W/Ui/2)+(-(-Ui*(Ei-Ui)*(Ei-W+Ui)^2*atanh(sqrt((W-Ui)/Ui)*Ui*(Ei*Ui)^(-1/2))+(-2*Ei*atan(sqrt((W-Ui)/Ui))*W^2+Ui*sqrt((W-Ui)/Ui)*(Ei+Ui)*(Ei-W+Ui))*sqrt(Ei*Ui))*(Ei*Ui)^(-1/2)/(Ei-W+Ui)^2/(Ei+Ui)/Ei/4)+(-2/3*Ui*sqrt((W-Ui)/Ui)*(W-Ui)/W^3)+(-Ui*((Ei-Ui)*(Ei^2+10*Ui*Ei+Ui^2)*(Ei-W+Ui)^2*atanh(sqrt((W-Ui)/Ui)*Ui*(Ui*Ei)^(-1/2))-sqrt(Ui*Ei)*(16*Ui*Ei*(Ei-W+Ui)^2*atan(sqrt((W-Ui)/Ui))+(Ei^3+(W-9*Ui)*Ei^2+(8*W*Ui-9*Ui^2)*Ei-Ui^2*(W-Ui))*sqrt((W-Ui)/Ui)*(Ei+Ui)))*(Ui*Ei)^(-1/2)/(Ei+Ui)^3/(Ei-W+Ui)^2/Ei/4)
                Gi = G₁_hydrogen(W⁺) - G₁_hydrogen(W⁻)
            else
                G₁(W) = (-(Ui*Ti)^(-1/2)*(-Ti*atan(sqrt((W-Ui)/Ti)*Ti*(Ui*Ti)^(-1/2))*W+atan(sqrt((W-Ui)/Ti))*(W-Ui+Ti)*sqrt(Ui*Ti))/W/(Ti-Ui))+(-(-Ti*(Ei^2+(3*Ti-3*Ui)*Ei-Ui*Ti)*(Ei-W+Ui)^2*atanh(sqrt((W-Ui)/Ti)*Ti*(Ei*Ti)^(-1/2))+sqrt(Ei*Ti)*(2*Ei*(W-Ui+Ti)*((Ti-Ui-W)*Ei+(Ui-2*W)*Ti-Ui*(-W+Ui))*atan(sqrt((W-Ui)/Ti))+Ti*sqrt((W-Ui)/Ti)*(Ei+Ui)*(Ei-W+Ui)*(Ti+Ei)))*(Ei*Ti)^(-1/2)/(Ei-W+Ui)^2/(Ti+Ei)^2/Ei/4)+(Ti*(Ui*Ti)^(-1/2)*(-W*(Ti^2+6*Ui*Ti+Ui^2)*(W-Ui+Ti)*atan(sqrt((W-Ui)/Ti)*Ti*(Ui*Ti)^(-1/2))+(4*W*(Ti+Ui)*(W-Ui+Ti)*atan(sqrt((W-Ui)/Ti))+sqrt((W-Ui)/Ti)*(Ti^2+3*Ti*W-(-W+Ui)*Ui)*(-Ui+Ti))*sqrt(Ui*Ti))/(-Ui+Ti)^3/W/(W-Ui+Ti))+(-((W-Ui+Ti)*(Ei-W+Ui)^2*(Ei^4+(13*Ti-3*Ui)*Ei^3-33*Ti*(-Ui+Ti)*Ei^2+(3*Ti^3-13*Ti^2*Ui)*Ei-Ti^3*Ui)*atanh(sqrt(-(-W+Ui)/Ti)*Ti*(Ti*Ei)^(-1/2))-sqrt(Ti*Ei)*(32*(W-Ui+Ti)*(Ei-W+Ui)^2*Ei*((Ti-Ui/2)*Ei-Ti*(Ti-2*Ui)/2)*atan(sqrt(-(-W+Ui)/Ti))+sqrt(-(-W+Ui)/Ti)*(Ti+Ei)*((W-Ui+Ti)*Ei^4+(-20*Ti^2+(24*Ui-11*W)*Ti-4*Ui^2+3*Ui*W+W^2)*Ei^3+(3*Ti^3+(-43*Ui+31*W)*Ti^2+(43*Ui^2-55*Ui*W+12*W^2)*Ti-3*Ui*(-W+Ui)^2)*Ei^2+4*((Ui-5/4*W)*Ti^2-6*(-W+Ui)*(Ui-13/24*W)*Ti+5*Ui*(-W+Ui)^2)*Ti*Ei+Ti^2*Ui*(-W+Ui)*(W-Ui+Ti))))*(Ti*Ei)^(-1/2)*Ti/(Ti+Ei)^4/(Ei-W+Ui)^2/Ei/(W-Ui+Ti)/4)
                Gi = G₁(W⁺) - G₁(W⁻)
            end
		else
			error("Integral is given only for n=0 or n=1.")
		end
		Gi *= 8*Ti/(3*π)
	end
	if isinf(Gi) || isnan(Gi) error("Infinite or NaN values calculating Gi.") end
	return Gi
end
