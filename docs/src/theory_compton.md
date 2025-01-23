# 2.1.2 Compton Scattering

The Compton scattering, often called incoherent scattering, consists of the ionization of an atomic electron by an incoming photon, which loses energy in the process. This interaction is dominant for intermediate energies ranging from the hundreds of keV to a few MeV. Penelope and EGSnrc models are based on the subshell-dependent relativistic impulse approximation from Brusa et al. [brusa1996fast](@cite). While using this model is desirable for multigroup deterministic transport, it is challenging to implement to have fast and accurate results. The production of Legendre moments of these cross-sections is an outstanding challenge since two very challenging integrations, in angle and energy of the scattered photon, must be performed, one following the other, in an accurate and fast enough way. Therefore, until such multigroup cross-sections can be developed, the proposed model uses the Waller-Hartree incoherent function, which is applied to the Klein-Nishina cross-section [heitler1984quantum](@cite), in order to encompass most electron binding effects [brusa1996fast](@cite).

The Compton cross-section describes the interaction of an incoming photon ($p'=\gamma$) with atomic electrons, resulting in a scattered photon ($p=\gamma$) and a produced electron ($p=\texttt{e+}$). The Klein-Nishina differential cross-section in the energy of the scattered photon for a single interaction with an assumed free atomic electron is given by [lorence1989physics,heitler1984quantum](@cite)

$$\sigma_{s}(E\rightarrow E') = \frac{\pi r_{e}^{2}}{E^{2}}\left[\frac{E}{E'}+\frac{E'}{E}-2\left(\frac{1}{E'}-\frac{1}{E}\right)+\left(\frac{1}{E'}-\frac{1}{E}\right)^{2}\right] \,,$$

where $E$ is the incoming photon energy, $E'$ is the scattered photon energy and $W = E - E'$ is the produced electron energy. The scattering angles for the scattered photon and the produced electron are respectively

$$\mu_{\gamma} = 1 + \frac{1}{E} - \frac{1}{E'} \quad \text{and} \quad \mu_{e} = \frac{1+E}{E}\left[1+\frac{2}{W}\right]^{-\frac{1}{2}} \,.$$

The double differential cross-sections for the scattered photon and the produced electron, with an incoherent scattering factor $S_{i}(E,\mu)$ extracted from JENDL-5 library [iwamoto2023japanese](@cite) taking into account some binding effects [brusa1996fast](@cite), are given by

$$\sigma_{s}^{i}(E\rightarrow E',\mu) = \dfrac{1}{2\pi}S_{i}(E,\mu)\sigma_{s}(E\rightarrow E')\delta(\mu-\mu_{\gamma})$$

and

$$\sigma_{s}^{i}(E\rightarrow W,\mu) = \dfrac{1}{2\pi}S_{i}(E,\mu)\sigma_{s}(E\rightarrow W)\delta(\mu-\mu_{e}) \,.$$

## 2.1.2.1 Scattering cross-sections of the incoming photon

The Compton Legendre moments of the differential scattering cross-sections of the incoming photon are given by

$$\sigma_{s,\ell}^{i}(E\rightarrow E') = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i}(E\rightarrow E',\mu) = P_{\ell}(\mu_{\gamma})S_{i}(E,\mu_{\gamma})\sigma_{s}(E\rightarrow E') \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\gamma\rightarrow\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{\max\left\{E'_{g+1/2},\frac{E}{1+2E}\right\}}^{\min\left\{E'_{g-1/2},E\right\}} dE' \sigma_{s,\ell}^{i}(E\rightarrow E') \,,$$

which are integratedby numerical quadrature.

## 2.1.2.2 Scattering cross-sections of the produced electron

The Compton Legendre moments of the differential scattering cross-sections of the produced electron are given by

$$\sigma_{s,\ell}^{i}(E\rightarrow W) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i}(E\rightarrow W,\mu) = P_{\ell}(\mu_{e})S_{i}(E,\mu_{e})\sigma_{s}(E\rightarrow W) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\gamma\rightarrow\texttt{e-}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{E'_{g+1/2}}^{\min\left\{E'_{g-1/2},\frac{2E^2}{1+2E}\right\}} dW \sigma_{\ell}^{i}(E\rightarrow W) \,,$$

which are integrated by numerical quadrature.

## 2.1.2.3 Total cross-sections

The Compton total cross-sections are defined by

$$\Sigma_{t}^{\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{\frac{E}{1+2E}}^{E}dE'S_{i}(E,\mu_{\gamma})\sigma_{s,0}(E\rightarrow E') \,,$$

which are integrated by numerical quadrature.