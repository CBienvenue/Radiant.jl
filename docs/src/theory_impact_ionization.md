# 2.2.1 Impact Ionization

When travelling through the medium, an electron or a positron ionizes atomic electrons in impact ionization. When an inner-shell electron is ionized, a relaxation cascade is initiated. This interaction is the dominant energy-loss phenomenon up to a few MeV for electrons and positrons [fernandez2005relativistic](@cite). The BFP formalism is also used to treat impact ionization. The catastrophic cross-sections are defined by subshell-dependant Møller [moller1932theorie](@cite) and Bhabha [bhabha1936scattering](@cite) models, for electrons and positrons respectively, based on models developed by Seltzer [seltzer1988cross,perkins1989livermore](@cite) and Salvat et al. [salvat1992semiempirical,fernandez2005relativistic](@cite). The soft stopping powers are developed based on the state-of-the-art stopping power models of Salvat and Andreo [salvat2022bethe](@cite).

## 2.2.1.1 Catastrophic Inelastic Cross-Section

he total inelastic cross-section describes the interaction of an incoming electron ($p'=\texttt{e-}$) or positron ($p'=\texttt{e+}$) with an atomic electron. Such interaction results in the scattering of the incoming particle ($p'=\texttt{e-}$ or $p'=\texttt{e+}$) and the production of knock-on electrons ($p'=\texttt{e-}$). The impact ionization is the sum of three components: the close collisions, the distance collisions and a density effect term [perkins1989livermore](@cite). Distance collisions consist of the interaction of the incident electron or positron perturbing field, which can ionize atomic electrons and be interpreted as virtual photons producing photoionization [seltzer1988cross](@cite). Since catastrophic collisions do not deal with small energy losses and distant collisions result in such energy transfer, they are neglected in the catastrophic cross-section model [salvat2022inelastic](@cite). The density effect component accounts for the material dielectric polarization and can also be neglected in catastrophic cross-sections since by inspection, the corresponding Scofield cross-section [scofield1978k,perkins1989livermore](@cite) is important only when small energy losses occur.

The collisional inelastic differential cross-section in energy of the knock-on electron in a monoelemental material $i$ is given by the sum of the differential cross-section per electron subshell $k$ [salvat2006penelope,lorence1989physics,salvat2009overview](@cite)

$$\sigma_{s}^{i}(E\rightarrow W) = \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \sigma_{s}^{i,k}(E\rightarrow W) \,,$$

where $Z_{i,k}$ is the mean number of electrons in subshells $k$ for element $i$, $E$ is the incoming electron or positron energy, $E'$ is the scattered electron or positron energy and $W = E - U_{i,k} - E'$ is the knock-on electron energy, where $U_{i,k}$ is the binding energy of the $k^{\text{th}}$ subshell. The values of $Z_{i,k}$ and $U_{i,k}$ are extracted from the Evaluated Atomic Data Library (EADL) [perkins1991tables](@cite). The close interaction cross-section is given by

$$\sigma_{s}^{i,k}(E\rightarrow W) = \dfrac{2\pi r_{e}^{2}}{\beta^{2}} F^{\pm}_{i,k}(E,W) \,,$$

where the subshell-dependant Møller factor, for electrons, is given by [moller1932theorie,seltzer1988cross](@cite)

$$F^{-}_{i,k}(E,W) = \frac{1}{(W+U_{i,k})^2} + \frac{1}{(E-W)^2} + \frac{1}{(E+1)^2} - \frac{(2E+1)}{(E+1)^2(E-W)(W+U_{i,k})}$$

and the subshell-dependant Bhabha factor, for positrons, is given by [bhabha1936scattering](@cite)

$$\begin{split}
    F^{+}_{i,k}(E,W) &=  \dfrac{1}{(W+U_{i,k})^2}\left[ 1 - b_{1}\left(\frac{W+U_{i,k}}{E}\right) + b_{2}\left(\frac{W+U_{i,k}}{E}\right)^2 \right. \\ &\left.- b_{3}\left(\frac{W+U_{i,k}}{E}\right)^3 + b_{4}\left(\frac{W+U_{i,k}}{E}\right)^4 \right] \,,
  \end{split}$$

with

$$\begin{split}
        &b_{1} = \left(\frac{\gamma-1}{\gamma}\right)^2 \frac{2(\gamma+1)^2-1}{\gamma^2-1} \, , \, b_{2} = \left(\frac{\gamma-1}{\gamma}\right)^2 \frac{3(\gamma+1)^2+1}{(\gamma+1)^2} \,, \\ &b_{3} = \left(\frac{\gamma-1}{\gamma}\right)^2 \frac{2\gamma(\gamma-1)}{(\gamma+1)^2} \, , \,\,\,\,\,\,\,\,\,\, b_{4} = \left(\frac{\gamma-1}{\gamma}\right)^2 \frac{(\gamma-1)^2}{(\gamma+1)^2} \,.
    \end{split}$$

The scattering angles for the incoming particle and knock-on electron are respectively [lorence1989physics](@cite)

$$\mu_{p} = \sqrt{\frac{E'(E+2)}{E(E'+2)}} \quad \text{and} \quad \mu_{s} = \sqrt{\frac{W(E+2)}{E(W+2)}}  \,,$$

and their double differential cross-sections are respectively given by

$$\sigma_{s}(E\rightarrow E',\mu) = \dfrac{1}{2\pi}\sigma_{s}(E\rightarrow E')\delta(\mu-\mu_{p})$$

and

$$\sigma_{s}(E\rightarrow W,\mu) = \dfrac{1}{2\pi}\sigma_{s}(E\rightarrow W)\delta(\mu-\mu_{s}) \,.$$

The maximum energy of knock-on for incoming electron and positron, which results from classical kinematics and the indistinguishability of electrons, are given respectively by

$$W_{\texttt{max}}^{-} = \frac{E-U_{i,k}}{2} \quad \text{and} \quad  W_{\texttt{max}}^{+} = E-U_{i,k} \,.$$

### 2.2.1.1.1 Scattering cross-sections for incoming electrons or positrons

The Legendre moments of the differential scattering cross-sections are given by

$$\sigma_{s,\ell}^{i,k}(E\rightarrow E') = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i,k}(E\rightarrow E',\mu) = P_{\ell}(\mu_{p})\sigma_{s}^{i,k}(E\rightarrow E') \,.$$

The macroscopic Legendre moments of the scattering cross-sections for particle scattering from group $g'$ are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e±} \rightarrow \texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \displaystyle\int_{\max\left\{E'_{g+1/2},W_{\texttt{max}}^{\pm}\right\}}^{\min\left\{E'_{g-1/2},E_{c}(E),E-U_{i,k}\right\}}dE' \sigma_{s,\ell}^{i,k}(E\rightarrow E') \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

### 2.2.1.1.2 Scattering cross-sections for knock-on electrons

The Legendre moments of the differential scattering cross-sections are given by

$$\sigma_{s,\ell}^{i,k}(E\rightarrow W) = 2\pi\int_{-1}^{1}d\mu\,P_{\ell}(\mu)\sigma_{s}^{i,k}(E\rightarrow W,\mu) = P_{\ell}(\mu_{s})\sigma_{s}^{i,k}(E\rightarrow W) \,.$$

The multigroup Legendre moments of the scattering cross-sections are therefore given by

$$\Sigma_{s,\ell,g'}^{\texttt{e±} \rightarrow \texttt{e-}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \displaystyle\int_{E'_{g+1/2}}^{\min\left\{E'_{g-1/2},W_{\texttt{max}}^{\pm}\right\}}dW\sigma_{s,\ell}^{i,k}(E\rightarrow W)  \mathcal{H}_{b} \,.$$

This equation is evaluated using numerical quadrature.

### 2.2.1.1.3 Catastrophic total cross-sections

The catastrophic inelastic total cross-sections are defined by

$$\Sigma_{t}^{\texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \int_{E-E_{c}(E)}^{W_{\texttt{max}}^{\pm}} dW  \sigma_{s,0}^{i,k}(E\rightarrow W) \mathcal{H}_{b} \,,$$

which is evaluated analytically.

\subsubsection{Absorption cross-sections for incoming positrons}
The inelastic catastrophic absorption cross-sections for inelastic interaction with incoming positrons, which include both soft and catastrophic contributions and are required for annihilation calculations, are given by

$$\Sigma_{a}^{\texttt{inel}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \int_{E-\min\left\{E_{G+1/2},E_{c}\right\}}^{W_{\texttt{max}}^{+}} dW  \sigma_{s,0}^{i,k}(E\rightarrow W) \mathcal{H}_{b} \,,$$

which is evaluated analytically.

## 2.2.1.2 Soft stopping powers

### 2.2.1.2.1 Total stopping power

The total collisional stopping powers of electron and positron for any compound, in $m_{e}c^{2} \times \text{cm}^{-1}$, are given by the density- and shell-corrected Bethe formula [salvat2023sbethe](@cite)

$$S_{t}(E) = \frac{2\pi r_{e}^2}{\beta^2} \mathcal{N}_{e}^{\texttt{eff}} \left[ \ln\left\{\frac{E+2}{2}\left(\frac{E}{I^{\texttt{eff}}}\right)^2\right\}  + f^{(\pm)} - \delta_{F} - 2\mathcal{C}(E) \right] \,,$$

where $E$ is the energy of the incoming electron or positron, $\mathcal{N}_{e}^{\texttt{eff}}$ is the effective electron density in the medium given by

$$\mathcal{N}_{e}^{\texttt{eff}} = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i}Z_{i} \,,$$

with $f_{i}$ is the weight percent of the $i^{\text{th}}$-element of the compound, $Z_{i}$ is the atomic number of the $i^{\text{th}}$-element of the compound. $\delta_{F}$ is the Fermi density effect, defined in the next subsection, and $\mathcal{C}(E)$ is the shell correction, extracted from SBETHE program files [salvat2023sbethe](@cite), a state-of-the-art program to compute stopping power of charged particle. $I^{\texttt{eff}}$ is the effective mean excitation energy of the compound is given by [salvat2006penelope](@cite)

$$I^{\texttt{eff}} = \exp\left\{\frac{1}{Z^{\texttt{eff}}}\displaystyle\sum_{i=1}^{N_{e}} f_{i}Z_{i} \log(I_{i})\right\}  \quad \text{with} \quad Z^{\texttt{eff}} = \sum_{i=1}^{N_{e}} f_{i}Z_{i} \,,$$

where the mean excitation energy of the $i^{\text{th}}$-element of the compound is $I_{i}$ (tabulated by Seltzer and Berger [seltzer1982evaluation](@cite)), unless more accurate value is provided. For water, $I^{\texttt{eff}} = 78$ eV is used, as recommended by the ICRU report 90 [seltzer2016key,salvat2022bethe](@cite). The electron factor $f^{(-)}$ is given by [rohrlich1954positron](@cite)

$$f^{(-)} = 1 - \beta^{2} - \frac{(2E+1)}{(E+1)^{2}}\ln(2) + \frac{1}{8}\left(\frac{E}{E+1}\right)^2$$

and the positron factor $f^{(+)}$ is given by

$$f^{(+)} = 2\ln(2) - \frac{\beta^2}{12}\left[23 + \frac{14}{E+2} + \frac{10}{(E+2)^2} + \frac{4}{(E+2)^3}\right] \,.$$

### 2.2.1.2.2 Fermi density effect

The following calculation of the Fermi density effect is the same as described in SBETHE and the same used in the Monte Carlo transport code PENELOPE [salvat2006penelope](@cite). It is based on the formula from Fano [fano1956atomic,inokuti1982fermi](@cite), which is

$$\delta_{F} = \frac{1}{Z^{\texttt{eff}}}\sum_{i=1}^{N_{e}}\sum_{k=1}^{N_{\texttt{shells}}} f_{i}Z_{i,k}\ln\left(1+\frac{L^2}{W_{i,k}^2}\right) - \frac{L^2}{\Omega_{p}^2}(1-\beta^2) \,,$$

where $L$ is given by solving

$$(1-\beta^2) = \frac{\Omega_{p}^{2}}{Z^{\texttt{eff}}}\sum_{i=1}^{N_{e}}\sum_{k=1}^{N_{\texttt{shells}}}\frac{f_{i}Z_{i,k}}{W_{i,k}^2+L^2} \,.$$

which can be done using the Newton-bisection method. The plasma energy is given by

$$\Omega_{p}^2 = 4\pi\mathcal{N}_{e}^{\texttt{eff}}\hbar^2 c^2r_{e}^{2}$$

and the resonance energy in the $i^{\texttt{th}}$-subshell is [sternheimer1952density](@cite)

$$W_{i,k} = \sqrt{(aU_{i,k})^2+\frac{2}{3}\frac{Z_{i,k}}{Z_{i}}\Omega_{p}^2} \,,$$

where $a$ is given by solving

$$Z_{i}\ln\left(I^{\texttt{eff}}\right) = \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k}\ln\left(\sqrt{(aU_{i,k})^2+\frac{2}{3}\frac{Z_{i,k}}{Z_{i}}\Omega_{p}^2}\right) \,,$$

again using the Newton-bisection method.

### 2.2.1.2.3 Soft stopping powers

The inelastic collisional soft stopping power, $S^{\texttt{e±}}(E)$, is given by removing the catastrophic stopping power from the total stopping power. It is given by

$$S^{\texttt{e±}}(E) = S_{t}(E) - \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \sum_{k=1}^{N_{\texttt{shells}}} Z_{i,k} \displaystyle\int_{E-E_{c}(E)}^{W_{\texttt{max}}^{\pm}}dW (W+U_{i,k}) \sigma_{s,0}^{i,k}(E\rightarrow W) \mathcal{H}_{b} \,,$$

which is evaluated analytically.






