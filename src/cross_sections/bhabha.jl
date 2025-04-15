"""
    bhabha(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0)

Gives the Bhabha differential cross-sections.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `W::Float64` : energy lost by the positron.
- `Ui::Float64` : binding energy of the electron(s).

# Output Argument(s)
- `σs::Float64` : Bhabha differential cross-sections.

# Reference(s)
- Bhabha (1936), The scattering of positrons by electrons with exchange on Dirac’s theory
  of the positron.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function bhabha(Zi::Real,Ei::Float64,W::Float64,Ui::Float64=0.0)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    γ = Ei+1
    
    # Bhabha formula
    b = ((γ-1)/γ)^2
    b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
    b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
    b3 = b * (2*(γ-1)*γ)/(γ+1)^2
    b4 = b * (γ-1)^2/(γ+1)^2
    F = (1 - b1*((W+Ui)/Ei) + b2*((W+Ui)/Ei)^2 - b3*((W+Ui)/Ei)^3 + b4*((W+Ui)/Ei)^4)/(W+Ui)^2
    σs = 2*π*rₑ^2/β² * Zi * F

    return σs
end

"""
    integrate_bhabha(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=1e-5,Wmax::Float64=Ei)

Gives the integration of the Bhabha differential cross-sections multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Wmin::Float64` : minimum energy lost by the positron.
- `Wmax::Float64` : maximum energy lost by the positron.

# Output Argument(s)
- `σn::Float64` : integrated Bhabha cross-sections.

# Reference(s)
- Bhabha (1936), The scattering of positrons by electrons with exchange on Dirac’s theory
  of the positron.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function integrate_bhabha(Z::Int64,Ei::Float64,n::Int64,Wmin::Float64=1e-5,Wmax::Float64=Ei)
  Nshells,Zi,Ui,_,_,_ = electron_subshells(Z)
  σn = 0
  for δi in range(1,Nshells)
      σn += Zi[δi] * integrate_moller_per_subshell(Z,Ei,n,Ui[δi],Wmin,Wmax)
  end
  return σn
end

"""
    integrate_bhabha_per_subshell(Z::Int64,Ei::Float64,n::Int64,Ui::Float64,
    Wmin::Float64=1e-5,Wmax::Float64=Ei)

Gives the integration of the Bhabha differential cross-sections for an electron with a given
binding energy, multiplied by Wⁿ.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `n::Int64` : order of the energy-loss factor.
- `Ui::Float64` : binding energy of the subshell.
- `Wmin::Float64` : minimum energy lost by the incoming electron.
- `Wmax::Float64` : maximum energy lost by the incoming electron.

# Output Argument(s)
- `σni::Float64` : integrated Bhabha cross-sections in a given subshell.

# Reference(s)
- Møller (1932), Zur theorie des durchgangs schneller elektronen durch materie.
- Perkins (1989), The Livermore electron impact ionization data base.
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.
- Salvat and Fernández-Varea (2009), Overview of physical interaction models for photon and
  electron transport used in Monte Carlo codes.

"""
function integrate_bhabha_per_subshell(Z::Int64,Ei::Float64,n::Int64,Ui::Float64,Wmin::Float64=0.0,Wmax::Float64=Ei)

    # Varibles
    rₑ = 2.81794092E-13       # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    γ = Ei+1
    b = ((γ-1)/γ)^2
    b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
    b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
    b3 = b * (2*(γ-1)*γ)/(γ+1)^2
    b4 = b * (γ-1)^2/(γ+1)^2

    # Compute the Bhabha cross-section per subshell
    σni = 0
    W⁺ = min(Ei-Ui,Wmax)
    W⁻ = max(0.0,Wmin)
    if W⁺ > W⁻
        if n == 0
          J₀(W) = -1/(W+Ui) - b1*log(W+Ui)/Ei + b2*W/Ei^2 - b3*(W^2/2+Ui*W)/Ei^3 + b4*(W^3/3+Ui*W^2+Ui^2*W)/Ei^4
          σni = J₀(W⁺) - J₀(W⁻)
        elseif n == 1
          J₁(W) = log(W+Ui) - b1*W/Ei + b2*(W^2/2+Ui*W)/Ei^2 - b3*(W^3/3+Ui*W^2+Ui^2*W)/Ei^3 + b4*(W^4/4+Ui*W^3+3*Ui^2*W^2/2+Ui^3*W)/Ei^4
          σni = J₁(W⁺) - J₁(W⁻)
        else
            error("Integral is given only for n=0 or n=1.")
        end
    end
    σni *= 2*π*rₑ^2/β²
    return σni
end