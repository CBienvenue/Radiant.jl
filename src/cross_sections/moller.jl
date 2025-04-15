"""
    moller(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0)

Gives the Møller differential cross-sections.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `W::Float64` : energy lost by the incoming electron.
- `Ui::Float64` : binding energy of the electron(s).

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
function moller(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    
    # Møller formula
    F = 1/(W+Ui)^2 + 1/(Ei-W)^2 + 1/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * 1/((Ei-W)*(W+Ui))
    σs = 2*π*rₑ^2/β² * Zi * F

    return σs
end


"""
    integrate_moller(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=1e-5,Wmax::Float64=Ei)

Gives the integration of the Møller differential cross-sections multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Wmin::Float64` : minimum energy lost by the incoming electron.
- `Wmax::Float64` : maximum energy lost by the incoming electron.

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
function integrate_moller(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=0.0,Wmax::Float64=Ei/2)
    Nshells,Zi,Ui,_,_,_ = electron_subshells(Z)
    σn = 0
    for δi in range(1,Nshells)
        σn += Zi[δi] * integrate_moller_per_subshell(Z,Ei,n,Ui[δi],Wmin,Wmax)
    end
    return σn
end

"""
    integrate_moller_per_subshell(Z::Int64,Ei::Float64,n::Int64,Ui::Float64,
    Wmin::Float64=1e-5,Wmax::Float64=Ei)

Gives the integration of the Møller differential cross-sections for an electron with a given
binding energy, multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Ui::Float64` : binding energy of the subshell.
- `Wmin::Float64` : minimum energy lost by the incoming electron.
- `Wmax::Float64` : maximum energy lost by the incoming electron.

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
function integrate_moller_per_subshell(Z::Int64,Ei::Float64,n::Int64,Ui::Float64,Wmin::Float64=Ui,Wmax::Float64=Ei/2)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2

    # Compute the Møller cross-section per subshell
    σni = 0
    W⁺ = min((Ei-Ui)/2,Wmax)
    W⁻ = max(Ui,Wmin)
    if W⁺ > W⁻
        if n == 0
            J₀(W) = -1/(W+Ui) + 1/(Ei-W) + W/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(W+Ui)-log(Ei-W))/(Ei+Ui)
            σni += (J₀(W⁺) - J₀(W⁻))
        elseif n == 1
            J₁(W) = log(W+Ui) + log(Ei-W) + (Ei+Ui)/(Ei-W) + W*(W+2*Ui)/(2*(Ei+1)^2) + (2*Ei+1)/(Ei+1)^2 * log(Ei-W)
            σni += (J₁(W⁺) - J₁(W⁻))
        else
            error("Integral is given only for n=0 or n=1.")
        end
    end
    σni *= 2*π*rₑ^2/β²
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