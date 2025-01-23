# 2.2.3 Bremsstrahlung

The bremsstrahlung interaction occurs when an incoming electron or positron interacts with the field of the atomic nucleus and its electrons. This interaction, vital at high energies, slows the incoming particle by producing a bremsstrahlung photon. The tables of Seltzer and Berger [seltzer1986bremsstrahlung](@cite) are used for both differential cross-sections and stopping powers, which are considered the most accurate and comprehensive energy-dependent bremsstrahlung data available [salvat2022bethe](@cite), and a modified dipole distribution of emitted photons distribution is used [acosta2002monte,kissel1983shape](@cite).

## 2.2.3.1 Catastrophic Bremsstrahlung Cross-Section

The bremsstrahlung cross-section describes the interaction of an incoming electron ($p'=\texttt{e-}$) or positron ($p'=\texttt{e+}$) with the field of the atomic nucleus and its electrons. The incoming particle scatter ($p=\texttt{e-}$ or $p=\texttt{e+}$), while photon are produced ($p=\gamma$). Seltzer and Berger proposed tables to compute the differential cross-section as a function of the energy of the produced photon, which includes both electron-nucleus and electron-electron interactions, and for monoelemental material $i$, it is given by [seltzer1986bremsstrahlung](@cite)

$$\sigma_{s}^{i}(E \rightarrow E_{\gamma}) = \begin{cases} \dfrac{F_{p,i}^{\pm}Z_{i}^2\tilde{\sigma}_{s}(Z_{i},E,E_{\gamma})}{E_{\gamma}\beta^2} & E_{\gamma} \le E \\ 0 & \text{otherwise} \end{cases} \,\,\,,$$

where $E$ is incoming electron or positron energy, $E_{\gamma}$ is the produced photon energy and $\tilde{\sigma}_{s}(Z_{i},E,E_{\gamma})$ is the total scaled bremsstrahlung energy-weighted cross-section, given in $\text{cm}^{2}$, from interpolation of the Seltzer and Berger tables [seltzer1986bremsstrahlung](@cite). These cross-sections are defined for $Z_{i} \le 100$ and energies between 1 keV and 10 GeV. The factor $F_{p,i}^{\pm}$ is given by [salvat2023sbethe](@cite)

$$\begin{split} F_{p,i}^{+} &= 1 - \exp\left\{ -1.2359\times10^{-1}t_{i} + 6.1274\times10^{-2}t_{i}^2 - 3.1516\times10^{-2}t_{i}^3 + 7.7446\times10^{-3}t_{i}^4 \right. \\ &\left. - 1.0595\times10^{-3}t_{i}^5 + 7.0568\times10^{-5}t_{i}^6 - 1.8080\times 10^{-6}t_{i}^7 \right\} \,, \end{split}$$

with

$$t_{i} = \ln\left(1+\frac{10^{6}E}{Z_{i}^2}\right)$$

for incoming positrons, which is based on the tabulated positron-to-electron ratios from Kim et al. [kim1986ratio](@cite), and $F_{p,i}^{-} = 1$ for incoming electrons. The double differential cross-sections are given by [lorence1989physics](@cite)

$$\sigma_{s}^{i}(E\rightarrow E_{e},\mu) = \dfrac{1}{2\pi}\sigma_{s}^{i}(E\rightarrow E_{e})\delta(\mu-1) \quad \text{and} \quad \sigma_{s}^{i}(E\rightarrow E_{\gamma},\mu) = \dfrac{1}{2\pi}\sigma_{s}^{i}(E\rightarrow E_{\gamma})\Theta(E,\mu) \,,$$

where it is assumed, as in CEPXS [lorence1989physics](@cite), that the incoming electron does not deflect from its trajectory after interaction ($\mu_{e} = 1$). The bremsstrahlung photon angular distribution, greatly inspired by the modified dipole distribution from Acosta et al. [acosta2002monte,salvat2006monte](@cite), is given by

$$\Theta(E,E_{\gamma},\mu) = \frac{3(1-C^2)}{4(2A+B)(1-\mu C)^2}\left[ (A+B) + (A-B)\left(\frac{\mu-C}{1-\mu C}\right)^2 \right] \,,$$

where the parameters $A=A(E,E_{\gamma})$, $B=B(E,E_{\gamma})$ and $C=C(E,E_{\gamma})$, which have values strictly between 0 and 1, are adjusted by least-square method to fit the shape function from Po$\check{\text{s}}$kus [povskus2019shape](@cite) for $E \le 3$ MeV, and are set to give the dipole distribution ($A = 1$, $B = 0$ and $C = \beta$) for energy $E > 3$ MeV.

### 2.2.3.1.1 Scattering cross-sections for deflected electron or positron

The bremsstrahlung Legendre moments of the differential scattering cross-sections for the incoming lepton are given by

$$\sigma_{s,\ell}^{i}(E\rightarrow E_{e}) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i}(E\rightarrow E_{e},\mu) = \sigma_{s}^{i}(E\rightarrow E_{e}) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e±}\rightarrow\texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \displaystyle\int_{E'_{g+1/2}}^{\min\left\{E'_{g-1/2},E_{c}(E)\right\}}dE_{e}\sigma_{s,\ell}^{i}(E\rightarrow E_{e}) \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

### 2.2.3.1.2 Scattering cross-sections for produced photon

The Legendre moments of the differential scattering cross-sections for the produced bremsstrahlung photon are given by

$$\sigma_{s,\ell}^{i}(E\rightarrow E_{\gamma}) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i}(E\rightarrow E_{\gamma},\mu) = \sigma_{s}^{i}(E\rightarrow E_{\gamma}) \Theta_{\ell}(E) \,,$$

where the moment of the angular distributionis given by

$$\Theta_{\ell}(E) = \int_{-1}^{1}d\mu\,P_{\ell}(\mu)\Theta(E,\mu) \,.$$

These Legendre moments of the angular distribution can be computed analytically by expanding the Legendre polynomials in power of $\mu$ such as

$$\begin{split} \Theta_{\ell}(E,E_{\gamma}) &= \frac{3(1-C^2)}{4(2A+B)} \frac{1}{2^{\ell}} \sum_{k=0}^{\lfloor\ell/2\rfloor} C_{\ell,k} \left[(A+B)I^{\ell-2k,2}(C) + (A-B) \sum_{j=0}^{2} \alpha_{j} I^{\ell-2k+j,4}(C) \right] \end{split}$$

where $\alpha_{0} = C^2$, $\alpha_{1} = -2C$ and $\alpha_{2} = 1$. The remaining integrals are given by

$$I^{n,m}(a) = \int_{-1}^{1}d\mu \frac{\mu^{n}}{(1-a\mu)^m} = \mathcal{G}_{3}^{n,m}(a,1)-\mathcal{G}_{3}^{n,m}(a,-1)$$

and is solved using integral expressions from Gradshteyn et al. (Sect. 2.153) [gradshteyn2014table](@cite)

$$\mathcal{G}_{3}^{n,m}(a,x) = \begin{cases} -\dfrac{1}{a}\log\left|1-ax\right| & n = 0 \text{ and } m \neq 1 \\ \dfrac{(1-ax)^{1-m}}{1-m} & n = 0 \text{ and } m = 1 \\ -\dfrac{1}{a(n-m+1)}\left[\vphantom{\dfrac{1}{2}} (1-ax)^{1-m}x^{n}-n\mathcal{G}_{3}^{n-1,m}(a,x)\right] & n > 0 \text{ and } n-m \neq -1 \\ -\dfrac{1}{(1-m)}\left[\vphantom{\dfrac{1}{2}} (1-ax)^{1-m}x^{n+1}-(n-m+2)\mathcal{G}_{3}^{n,m+1}(a,x)\right] & n > 0 \text{ and } n-m = -1 \end{cases} \,\,\,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e±}\rightarrow\gamma}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \displaystyle\int_{E'_{g+1/2}}^{\min\left\{E'_{g-1/2},E\right\}}dE_{\gamma}\sigma_{s,\ell}^{i}(E\rightarrow E_{\gamma}) \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

### 2.2.3.1.3 Catastrophic total cross-sections

The catastrophic bremsstrahlung total cross-sections are defined by

$$\Sigma_{t}^{\texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{0}^{E_{c}(E)}dE_{e}\sigma_{s,0}^{i}(E\rightarrow E_{e})$$

and are evaluated using numerical quadrature.

### 2.2.3.1.4 Absorption cross-sections for incoming positrons

The bremsstrahlung catastrophic absorption cross-sections, for annihilation calculations, are given by

$$\Sigma_{a}^{\texttt{brem}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{0}^{E_{G+1/2}}dE_{e}\sigma_{s,0}^{i}(E\rightarrow E_{e})$$

and are evaluated using numerical quadrature.

## 2.2.3.2 Soft stopping power

### 2.2.3.2.1 Total stopping power

The radiative stopping powers of electron are given by

$$S_{t}(E) = \alpha r_{e}^2 (E+1) \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} F_{p,i}^{\pm} Z_{i}^2 \phi_{\texttt{rad}}(Z_{i},E) \,,$$

where $\phi_{\texttt{rad}}(Z_{i},E)$ are also given by the tables of Selzer and Berger [seltzer1986bremsstrahlung](@cite).

### 2.2.3.2.2 Soft stopping power

The soft radiative stopping powers, $S^{\texttt{e±}}(E)$, are given by removing catastrophic bremsstrahlung stopping powers from total stopping powers. They are given by

$$S^{\texttt{e±}}(E) = S_{t}(E) - \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{0}^{E_{c}(E)}dE_{e}(E-E_{e})\sigma_{s,0}^{i}(E\rightarrow E_{e}) \,,$$

where the integral is evaluated with numerical quadrature.