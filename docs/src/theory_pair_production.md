# 2.1.4 Pair Production

An electron-positron pair is produced from the incoming photon energy in the neighbours of the atom. The threshold for pair production is approximately 1.022 MeV, the required energy to constitute at least the electron and positron mass energy, the extra energy going into their kinetic energies. This phenomenon is the dominant photon process at high energy, about 10 MeV and more. A variation of the model of Baró et al. [baro1994analytical](@cite), based on the screening function of Tsai [tsai1974pair](@cite), is employed.

The pair production cross-section describes the absorption of a photon ($p' = \gamma$) and the production of an electron ($p=\texttt{e-}$) and a positron ($p=\texttt{e+}$). A variation of the semi-empirical model of Baró [baro1994analytical](@cite) is used and, in a monoelemental material $i$, the differential cross-section is given by

$$\sigma_{s}^{i}(E \rightarrow E') = \begin{cases} A(Z_{i},E)\left[ 2\left(\dfrac{1}{2}-\dfrac{E'+1}{E}\right)^2 \phi_{i,1}(E') + \phi_{i,2}(E') \right] & E > 2 \\ 0 & \text{otherwise} \end{cases} \,\,\,,$$

where $E$ is the incoming photon energy and $E'$ is the outgoing electron or positron energy. The screening function, derived from the one of Tsai [tsai1974pair](@cite), are given by

$$\phi_{i,1}(E') = \max\left\{g_{i,1}(E') + g_{i,0},0\right\} \quad \text{and} \quad \phi_{i,2}(E') = \max\left\{g_{i,2}(E') + g_{i,0},0\right\} \,,$$

where

$$\begin{split} g_{i,1}(E') &= \frac{7}{3} - 2\ln(1+b_{i}^2) -6b_{i}\arctan\left(\frac{1}{b_{i}}\right) - b_{i}^2\left[ 4 - 4b_{i}\arctan\left(\frac{1}{b_{i}}\right) - 3\ln\left(1+\frac{1}{b_{i}^2}\right) \right] \,, \\ g_{i,2}(E') &= \frac{11}{6} - 2\ln(1+b_{i}^2) -3b_{i}\arctan\left(\frac{1}{b_{i}}\right) - \frac{b_{i}^2}{2}\left[ 4 - 4b_{i}\arctan\left(\frac{1}{b_{i}}\right) - 3\ln\left(1+\frac{1}{b_{i}^2}\right) \right] \,, \end{split}$$

$$g_{i,0} = 4\ln(r_{s,i}) - 4f_{C,i} \,,$$

with

$$b_{i} = \frac{r_{s,i}}{2}\frac{E}{(E'+1)(E-E'-1)} \,.$$

The variable $r_{s,i}$ corresponds to the reduced screening radius, tabulated in Baró [baro1994analytical](@cite) for $Z_{i} \le 92$ and extended up to $Z_{i} \le 99$ in PENELOPE [salvat2006penelope](@cite). The high-energy Coulomb correction $f_{C,i}$ is given by [davies1954theory](@cite)

$$f_{C,i} = \alpha^2 Z_{i}^2 \displaystyle\sum_{k=1}^{\infty}\frac{1}{k(k^2+\alpha^2 Z_{i}^2)} \,.$$

As it is done in EGSnrc [kawrakow2000egsnrc](@cite), a normalization factor $A(Z_{i},E)$ is added to the pair production cross-section, which is defined as the ratio between the total cross-sections obtained with $A(Z_{i},E)=1$ and the tabulated values from JENDL-5 for the pair production in both the nuclear and electron field [iwamoto2023japanese](@cite)

$$A(Z_{i},E) = \frac{\sigma_{t}^{\texttt{JENDL-5}}(Z_{i},E)}{\left.\sigma_{t}(Z_{i},E)\right|_{A=1}} \,.$$

Since pair production and bremsstrahlung are closely related through a substitution rule, the same angular distribution can be used for the electron and positron emission [lorence1989physics,tsai1974pair](@cite). The double differential cross-section is given by

$$\sigma_{s}^{i}(E\rightarrow E',\mu) = \dfrac{1}{2\pi}\sigma_{s}^{i}(E\rightarrow E')\Theta(E,\mu) \,.$$

This distribution is the Sommerfield distribution.

## 2.1.4.1 Scattering cross-sections for the produced electron and positron

The pair production Legendre moments of the differential scattering cross-sections for the produced leptons are given by

$$\sigma_{s,\ell}^{i}(E\rightarrow E') = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i}(E\rightarrow E',\mu) = \sigma_{s}^{i}(E\rightarrow E') \Theta_{\ell}(E) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\gamma\rightarrow\texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \displaystyle\int_{E'_{g+1/2}}^{\min\left\{E'_{g-1/2},E - 2\right\}}dE_{e}\sigma_{s,\ell}^{i}(E\rightarrow E') \mathcal{H}_{b} \,,$$

which are integrated by numerical quadrature. 

## 2.1.4.2 Total cross-sections

The pair production total cross-sections are defined by

$$\Sigma_{t}^{\gamma}(E) =  \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{0}^{E-2}dE'\sigma_{s,0}^{i}(E\rightarrow E') \mathcal{H}_{b}$$

and are integrated using numerical quadrature.

## 2.1.4.3 Absorption cross-sections

The pair production absorption cross-sections, which are required for annihilation calculations, are given by

$$\Sigma_{a}^{\texttt{pp}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i}f_{i}\int_{0}^{\min\left\{E-2,E_{G+1/2}\right\}}dE'\sigma_{s,0}^{i}(E\rightarrow E') \mathcal{H}_{b}$$

and are integrated using numerical quadrature.