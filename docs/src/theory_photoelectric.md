# 2.1.3 Photoelectric Effect

The photoelectric effect consists of the emission of electrons by absorption of photons. Tables from the subshell-dependent Evauated Photon Data Library (EPDL) [cullen1997epdl97](@cite), and the Sauter distribution is used.

The photoelectric cross-section describes the absorption of an incoming photon ($p'=\gamma$) and the emission of an atomic electron ($p=\texttt{e-}$). The microscopic absorption cross-sections are given by

$$\sigma_{a}^{i}(E) = \sum_{k=1}^{N_{\texttt{shells}}} \sigma_{a}^{i,k}(E) \,,$$

where $\sigma_{a}^{i,k}(E)$ is given by linear interpolation of the data given by the absorption cross-sections per subshells from the JENDL-5 library [iwamoto2023japanese](@cite), which are available for $Z_{i} \le 100$ and for energies up to 100 GeV. 

## 2.1.3.1 Scattering cross-sections of the produced electron

The photoelectric scattering cross-sections are given by

$$\sigma_{s}^{i,k}(E \rightarrow E',\mu) = \sigma_{a}^{i,k}(E) \delta(E'-E+U_{i,k}) \Theta(E,\mu) \,,$$

where $E'$ is the energy of the photo-electron, $U_{i,k}$ is the binding energy of the k$^{\text{th}}$ shell and the Sauter cross-section for the K-shell, normalized over the angular domain, is given by [sauter1931atomaren](@cite)

$$\Theta(E,\mu) = \Gamma(E) \frac{1-\mu^2}{(1-\beta\mu)^4}\left[1 + \frac{\gamma(\gamma-1)(\gamma-2)}{2}(1-\beta\mu)\right] \,,$$

where the normalization factor is

$$\Gamma(E) = \left\{\frac{4}{3(1-\beta^2)^2} + \frac{\gamma(\gamma-1)(\gamma-2)}{2\beta^3}\left[\frac{2\beta}{1-\beta^2}-\ln\left(\frac{1+\beta}{1-\beta}\right)\right]\right\}^{-1} \,.$$

The Legendre moments of the normalized Sauter cross-section are given by

$$\Theta_{\ell}(E) = \int_{-1}^{1}d\mu P_{\ell}(\mu) \Theta(E,\mu) \,.$$

The Legendre moments are computed analytically by expanding the Legendre polynomials in power of $\mu$ such as

$$\Theta_{\ell}(E) = \frac{\Gamma(E)}{2^{\ell}} \sum_{k=0}^{\lfloor\ell/2\rfloor} C_{\ell,k} \sum_{j=0}^{1} (-1)^j \left[ I^{\ell-2k+2j,4}(\beta) + \frac{\gamma(\gamma-1)(\gamma-2)}{2} I^{\ell-2k+2j,3}(\beta) \right]$$

The scattering cross-sections are given by

$$\Sigma_{s,\ell,g'}^{\gamma\rightarrow\texttt{e-}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{E'_{g+1/2}}^{E'_{g-1/2}} dE' \sum_{k=1}^{N_{\texttt{shells}}} \sigma_{s,\ell}^{i,k}(E \rightarrow E') \,,$$

which can be rewritten as

$$\Sigma_{s,\ell,g'}^{\gamma\rightarrow\texttt{e-}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} \sigma_{a}^{i,k}(E) \Theta_{\ell}(E) \times \begin{cases} 1 & E-U_{i,k} \in [E'_{g+1/2},E'_{g-1/2}] \\ 0 & \text{otherwise} \end{cases} \,\,\,.$$

## 2.1.3.2 Total cross-sections

The photoelectric total cross-sections are simply

$$\Sigma_{t}^{\gamma}(E) = \sum_{i=1}^{N_{e}}  \mathcal{N}_{n,i} f_{i} \sigma_{a}^{i}(E) \,.$$