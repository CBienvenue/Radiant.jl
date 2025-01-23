# 2.1.1 Rayleigh scattering

Rayleigh scattering, also often called coherent scattering, consists of the elastic change of direction of an incoming photon by atomic electrons. The Rayleigh scattering model from PENELOPE [salvat2006penelope](@cite), which is based on atomic form factors [hubbell1975atomic](@cite) and anomalous scattering factors [cromer1965anomalous](@cite), is used.

The Rayleigh cross-sections, which described the elastic scattering of photons ($p=p'=\gamma$) are given by [cromer1965anomalous,hubbell1975atomic,kissel1995validity](@cite)

$$\sigma_{s}^{i}(E,\mu) = \pi r_{e}^2 \left(1+\mu^2\right) \left[\left(F_{i}(E,\mu)+f_{i}'(E)\vphantom{\frac{1}{2}}\right)^{2} + \left(f_{i}''(E)\vphantom{\frac{1}{2}}\right)^{2}\right] \,,$$

where $F_{i}(E,\mu)$ is the atomic form factor for the $i^{\text{th}}$-element, where the factors $f_{i}'(E)$ and $f_{i}''(E)$ respectively are the real and imaginary parts of the anomalous scattering factors for the $i^{\text{th}}$-element, which are all tabulated by the Japanese evaluated nuclear data library (JENDL-5) [iwamoto2023japanese](@cite), which are based on the EPDL library [cullen1989tables](@cite). The double differential cross-sections are given by

$$\sigma_{s}^{i}(E\rightarrow E',\mu) = \sigma_{s}^{i}(E,\mu) \delta(E'-E) \,.$$

## 2.1.1.1 Scattering cross-sections of the incoming photon

The Rayleigh Legendre moments of the differential scattering cross-sections of the incoming photon are simply given by

$$\Sigma_{s,\ell}^{\gamma\rightarrow\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{-1}^{1} P_{\ell}(\mu)\sigma_{s}^{i}(E,\mu) \,,$$

which are integrated using numerical quadrature.

## 2.1.1.2 Total cross-sections

The Rayleigh total cross-sections are given by

$$\Sigma_{s,\ell}^{\gamma}(E) = \Sigma_{s,0}^{\gamma\rightarrow\gamma}(E) \,.$$