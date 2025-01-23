# 2.2.4 Annihilation

The annihilation interaction consists of an incoming positron that annihilates with an atomic electron to produce two photons. The combined photon's energy is equivalent to the sum of kinetic and mass energies of the annihilated atomic electron and incoming positron. As the kinetic energy of an incoming positron approaches zero, this interaction becomes inevitable. The model of Nelson et al. [nelson1985egs4](@cite) is employed. The presented annihilation model is an improvement over CEPXS since it is explicitly defined as a positron cross-section.

The following cross-sections describe the annihilation of an incoming positron ($p'=\texttt{e+}$) with an atomic electron producing two photons ($p=\gamma$), assuming that the electrons are free and at rest. The differential cross-sections in the energy of the lowest energy photon [salvat2009overview,heitler1984quantum](@cite)

$$\sigma_{s}(E \rightarrow E_{\gamma_{-}}) = \frac{\pi r_{e}^2}{(\gamma+1)^2(\gamma^2-1)} \left[ S(\zeta) + S(1-\zeta) \vphantom{\frac{1}{2}} \right]$$

with

$$S(\zeta) = -(\gamma+1)^2 + (\gamma^2+4\gamma+1)\frac{1}{\zeta} - \frac{1}{\zeta^2} \,,$$

where $E$ is the incoming positron energy, $E_{\gamma_{-}}$ is the lowest photon energy, $E_{\gamma_{+}} = E + 2 - E_{\gamma_{-}}$ is the highest photon energy and $\zeta = E_{\gamma_{-}}/(E+2)$ is the ratio of the lowest energy photon to the total (kinetic + mass) energy. The value of the lowest photon energy is bounded by

$$E^{\texttt{max}} = \frac{\gamma+1}{2} \quad \text{and} \quad E^{\texttt{min}} = \frac{\gamma+1}{\gamma+1+\sqrt{\gamma^2-1}} \,.$$

The scattering angles of the lowest and highest energy photons are respectively

$$\mu_{-} = \frac{1}{\sqrt{\gamma^2-1}} \left[ \gamma+1-\frac{1}{\zeta} \right] \quad \text{and} \quad \mu_{+} = \frac{1}{\sqrt{\gamma^2-1}} \left[ \gamma+1-\frac{1}{1-\zeta} \right] \,,$$

and the double differential cross-sections for the lowest and highest energy photons are given by

$$\sigma_{s}(E \rightarrow E_{\gamma_{-}},\mu) = \frac{1}{2\pi}\sigma_{s}(E \rightarrow E_{\gamma_{-}})\delta(\mu-\mu_{-})$$

and

$$\sigma_{s}(E \rightarrow E_{\gamma_{+}},\mu) = \frac{1}{2\pi}\sigma_{s}(E \rightarrow E_{\gamma_{+}})\delta(\mu-\mu_{+}) \,.$$

## 2.2.4.1 Scattering cross-sections for the lowest energy photons

The annihilation Legendre moments of the differential scattering cross-sections for the lowest energy photons are given by

$$\sigma_{s,\ell}(E\rightarrow E_{\gamma_{-}}) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}(E\rightarrow E_{\gamma_{-}},\mu) = P_{\ell}(\mu_{-})\sigma_{s}(E\rightarrow E_{\gamma_{-}}) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e+}\rightarrow\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} Z_{i} \displaystyle\int_{\max\left\{E^{\gamma}_{g+1/2},E^{\texttt{min}}\right\}}^{\min\left\{E^{\gamma}_{g-1/2},E^{\texttt{max}}\right\}}dE_{\gamma_{-}}\sigma_{s,\ell}(E\rightarrow E_{\gamma_{-}}) \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

## 2.2.4.2 Scattering cross-sections for the highest energy photons

The annihilation Legendre moments of the differential scattering cross-sections for the highest energy photons are given by

$$\sigma_{s,\ell}(E\rightarrow E_{\gamma_{+}}) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}(E\rightarrow E_{\gamma_{+}},\mu) = P_{\ell}(\mu_{+})\sigma_{s}(E\rightarrow E_{\gamma_{+}}) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e+}\rightarrow\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} Z_{i} \displaystyle\int_{\max\left\{E^{\gamma}_{g+1/2},E^{\texttt{max}}\right\}}^{\min\left\{E^{\gamma}_{g-1/2},(\gamma+1)-E^{\texttt{min}}\right\}}dE_{\gamma_{+}}\sigma_{s,\ell}(E\rightarrow E_{\gamma_{+}}) \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

## 2.2.4.3 Total cross-sections

The annihilation total cross-sections are defined by

$$\Sigma_{t}^{\texttt{e+}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} Z_{i} \int_{E^{\texttt{min}}}^{E^\texttt{max}} dE_{\gamma_{-}} \sigma_{s,0}(E \rightarrow E_{\gamma_{-}})$$

and are evaluated analytically.

## 2.2.4.4 Annihilation when positrons scatter under the cutoff energy

Positrons ($p'=\texttt{e+}$) will inevitably annihilate with atomic electrons and produce two photons ($p=\gamma$). Therefore, all positrons scattered under the cutoff energy $E_{G+1/2}$ should annihilate. The positron energy under the cutoff is small enough that the annihilation photon can be assumed to have isotropic scattering, as done in GEANT4 [collaboration2016physics](@cite). It is then assumed that two 511 keV photons are produced. The positrons are scattered under the cutoff following either an inelastic interaction or a bremsstrahlung interaction, and positrons can also be produced under the cutoff following pair production interaction. The $\ell = 0$ Legendre moments of the scattering cross-sections are given by

$$\Sigma_{s,0,g' \rightarrow g}^{\texttt{e+}\rightarrow\gamma} =  2\left[\Sigma_{a,g'}^{\texttt{inel}}+\Sigma_{a,g'}^{\texttt{brem}}\vphantom{\dfrac{1}{1}}\right] \times \begin{cases} 1 & 1 \in [E'_{g+1/2},E'_{g-1/2}] \\ 0 & \text{otherwise} \end{cases}$$

and

$$\Sigma_{s,0,g' \rightarrow g}^{\gamma\rightarrow\gamma} =  2\Sigma_{a,g'}^{\texttt{pp}} \times \begin{cases} 1 & 1 \in [E'_{g+1/2},E'_{g-1/2}] \\ 0 & \text{otherwise} \end{cases} \,\,\,,$$

where the $\ell \ge 1$ moments are equal to zero since scattering is isotropic. While these scattering cross-sections account for catastrophic impact ionization interactions, they do not include positrons reaching the cutoff energy by soft impact ionization interactions. In order to accurately estimate the production of annihilation photons following all impact ionization interactions, the lower bounds of integrals over energy loss $W$ should be $E-E_{G+1/2}$ rather than $E-\min\left\{E_{G+1/2},E_{c}\right\}$, because the latter excluding soft interactions. However, the Bhabha-based impact ionization model is inaccurate for small energy losses, i.e. for soft interaction, and it can lead to substantial error in the solution when used to describe such interactions. Since these soft events require an improved impact ionization model for positron, which would required enormous effort, the annihilation photons produced by absorption of positrons at the cutoff energy are neglected. This will result in an underestimation of the production of annihilation photon, which have an impact on the accuracy of transport calculations.

