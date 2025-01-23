# 2.2.2 Elastic Scattering

An interaction in which the particle changes direction without losing energy is said to be elastic. Elastic interaction of electron and positron with matter is of utmost importance to accurately model the distribution of this particle in the medium up to hundreds of MeV. A combined Molière screened Mott-Rutherford cross-section is employed, based on an interpolation formula and data proposed by Lijian et al. [lijian1995analytic,boschini2013expression](@cite), jointly with the knock-on dependent correction from Kawrakow [kawrakow1997improved](@cite) and Seltzer's adjustment to Molière screening factor [seltzer1988overview](@cite).

The Mott cross-section describes the elastic interaction of an incoming electron ($p=p'=\texttt{e-}$) or positron ($p=p'=\texttt{e+}$) with the Coulomb field of an atom. The microscopic differential scattering cross-section in a monoelemental material $i$ is given by

$$\sigma_{s}^{i}(E,\mu) = \frac{2\pi Z_{i}(Z_{i}+\xi^{\pm}_{i})r_{e}^2}{\beta^{2}E(E+2)}\frac{1}{(1-\mu+2\eta_{i})^2}\mathcal{R}_{\text{Mott}}^{i} \,,$$

with $\eta$, the Molière screening factor with Seltzer's adjustment, obtained by [lorence1989physics,moliere1947theorie,seltzer1988overview](@cite)

$$\eta_{i} = \frac{\alpha^2Z_{i}^{2/3}\left(1.13 + 3.76\left(\frac{Z_{i}\alpha}{\beta}\right)^2\sqrt{\frac{E}{E+1}}\right)}{4\left(\frac{9\pi^2}{128}\right)^{2/3}E(E+2)} \,,$$

where $\alpha\approx1/137$ is the fine structure constant. Note that when the incoming electron or positron has increasingly more kinetic energy, $\eta_{i}$ becomes tiny, and the scattering cross-section becomes increasingly large as $\mu$ tends to 1, producing highly forward-peaked scattering. Lijian et al. [lijian1995analytic](@cite) proposed the following interpolation formula for $\mathcal{R}_{\text{Mott}}^{i}$, the ratio of the unscreened Mott differential cross-section to Rutherford's cross-section:

$$\mathcal{R}_{\text{Mott}}^{i} = \sum_{j=0}^{4}a_{i,j}(1-\mu)^{j/2} \,, \quad \text{where} \quad a_{i,j} = \sum_{k=1}^{6}b_{k,j}(Z_{i})\left(\beta-\bar{\beta}\right)^{k-1}$$

and $\bar{\beta}=0.7181287$. Boschini et al. [boschini2013expression](@cite) have generated $b_{k,j}(Z_{i})$ parameters for any $Z_{i} \le 118$, for both electron and positron, and these are valid for energies between 1 keV and 900 MeV.

For incoming electrons or positrons, to take into account atomic electron contribution to the multiple scattering of charged particle transport, $Z^{2}$ is often replaced by $Z(Z+1)$ [fano1954inelastic](@cite). However, a double counting issue was observed in Monte Carlo codes due to an overlap with knock-on electrons [li1995electron](@cite), which are already taken into account by the Møller and Bhabha models. Kawrakow proposed a correction for this issue, which is significant for low Z material or when cutoff energy $E_{G+1/2}$ is low [kawrakow1997improved](@cite). This correction is given by

$$\xi^{\pm}_{i} = 1-\frac{g_{\texttt{inel},i}^{\pm}}{g_{\texttt{el},i}} \,,$$

with $g_{\texttt{el},i}$, which is calculated using methodology described by Kawrakow for the Mott cross-sections rather than the screened Rutherford ones, is given by

$$g_{\texttt{el},i} = \frac{\tilde{\sigma}_{s,0}^{i}(E)-\tilde{\sigma}_{s,2}^{i}(E)}{3} \,,$$

where $\tilde{\sigma}_{s,\ell}^{i}(E)$ is given by Eq.~\ref{altsigma}, and with $g_{\texttt{inel},i}$ is given by

$$g_{\texttt{inel},i}^{\pm} = \sum_{k=1}^{N_{\texttt{shells}}} \frac{Z_{i,k}}{Z_{i}} E(E+2) \int_{E_{G+1/2}}^{W_{\texttt{max}}^{\pm}}dW\left(1-\mu_{p}^2\right) F^{\pm}_{i,k}(E,W) \mathcal{H}_{b} \,,$$

where $E_{G+1/2}$ is the minimum energy transfer to knock-on electron. This equation is evaluated analytically.

## 2.2.2.1 Legendre moments of the scattering cross-sections

For deterministic transport calculations, the Legendre moments of the elastic cross-sections are required, namely

$$\Sigma_{s,\ell}^{\texttt{e±}\rightarrow\texttt{e±}}(E) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \int_{-1}^{1}d\mu P_{\ell}(\mu) \sigma_{s}^{i}(E,\mu) = \sum_{i=1}^{N_{e}} \mathcal{N}_{n,i} f_{i} \frac{2\pi Z_{i}(Z_{i}+\xi_{i}^{\pm})r_{e}^2}{\beta^{2}E(E+2)} \tilde{\sigma}_{s,\ell}^{i}(E) \,,$$

with the following definition

$$\tilde{\sigma}_{s,\ell}^{i}(E) = \int_{-1}^{1}d\mu P_{\ell}(\mu) \frac{\mathcal{R}_{\text{Mott}}^{i}}{(1-\mu+2\eta_{i})^2} \,.$$

Calculating these Legendre moments using numerical quadrature is highly inefficient due to the near singularity that often occurs at $\mu = 1$, when $\eta_{i}$ is small. As demonstrated in the following lines, the integration can be fast, done analytically and in a way that ensures numerical stability.

First, the Legendre polynomial can be expressed as a sum of powers of $\mu$ as [koepf1998hypergeometric](@cite)

$$P_{\ell}(\mu) = \frac{1}{2^{\ell}} \sum_{k=0}^{\lfloor\ell/2\rfloor} C_{\ell,k} \mu^{\ell-2k} \,, \quad \text{where} \quad C_{\ell,k} = \frac{(-1)^{k}(2\ell-2k)!}{k!(\ell-k)!(\ell-2k)!} \,,$$

and therefore, the Legendre moments can be expressed as

$$\tilde{\sigma}_{s,\ell}^{i}(E)  = \frac{1}{2^{\ell}}  \sum_{k=0}^{\lfloor\ell/2\rfloor} C_{\ell,k} \int_{-1}^{1} d\mu \frac{\mathcal{R}_{\text{Mott}}^{i}}{(1-\mu+2\eta_{i})^2} \mu^{\ell-2k} = \frac{1}{2^{\ell}}  \sum_{k=0}^{\lfloor\ell/2\rfloor} C_{\ell,k} \sum_{r=1}^{2}I_{r}^{\ell,k} \,,$$

where $I_{1}^{\ell,k}$ is given by

$$\begin{split} I_{1}^{\ell,k} &= \int_{-1}^{1} d\mu \frac{\mu^{\ell-2k}}{(1-\mu+2\eta_{i})^2} \left[ \alpha_{i,0} + \alpha_{i,1}\mu + \alpha_{i,2}\mu^2 \right] \\ &= \sum_{j=0}^{2} \alpha_{i,j} \left( \mathcal{G}_{1}^{\ell - 2k + j}(1 + 2\eta_{i},-1,1) - \mathcal{G}_{1}^{\ell - 2k + j}(1 + 2\eta_{i},-1,-1) \right) \end{split}$$

and $I_{2}^{\ell,k}$ is given by

$$\begin{split} I_{2}^{\ell,k} &= \int_{-1}^{1} d\mu \frac{\mu^{\ell-2k}\sqrt{(1-\mu)}}{(1-\mu+2\eta)^2} \left[ \alpha_{3} + \alpha_{4}\mu \right] \\ &= 2 \sum_{j=0}^{1} \alpha_{i,j+3} \sum_{g=0}^{\ell-2k+j} (-1)^g  \frac{(\ell-2k+j)!}{g!(\ell-2k+j-g)!} \mathcal{G}_{2}^{2+2g}(2\eta_{i},1,\sqrt{2}) \,, \end{split}$$

with $\alpha_{i,0} = a_{i,0}+a_{i,2}+a_{i,4}$, $\alpha_{i,1} = - (a_{i,2} + 2a_{i,4})$, $\alpha_{i,2} = a_{i,4}$, $\alpha_{i,3} = a_{i,1}+a_{i,3}$ and $\alpha_{i,4} = -a_{i,3}$. The first integral is obtained using the following integral from Gradshteyn et al. (Eq.4, Sect. 2.111) [gradshteyn2014table](@cite)

$$\mathcal{G}_{1}^{n}(a,b,x) = \sum_{g=1}^{n-1}(-1)^{g-1}\frac{ga^{g-1}x^{n-g}}{(n-g)b^{g+1}} + (-1)^{n-1}\frac{a^{n}}{b^{n+1}(a+bx)} + (-1)^{n+1}\frac{na^{n-1}}{b^{n+1}}\ln(a+bx) \,,$$

while the second integral is evaluated by applying the change of variable $u=\sqrt{1-\mu}$, using the binomial theorem and then formulae from Gradshteyn et al. (Sect. 2.172, Eq.~1 Sect. 2.173 and Eq.~1 Sect. 2.174)

$$\mathcal{G}_{2}^{n}(a,b,x) = \begin{cases} \dfrac{x}{2aR} + \dfrac{1}{2a\sqrt{ab}}\arctan{\left(\frac{bx}{\sqrt{ab}}\right)} & n = 0 \\ -\dfrac{x^{n-1}}{(3-n)bR} + \dfrac{(n-1)a}{(3-n)b}\mathcal{G}_{2}^{n-2}(a,b,x)  & n \text{ even} > 0 \end{cases} \,\,\,.$$

This fully analytical solution can be numerically unstable, because of catastrophic cancellations for high-order Legendre terms when $\eta_{i}$ is large, which happens when the elastic scattering is becoming more or less isotropic, notably at low energies. Since an isotropic flux is characterized by $\ell \ge 1$ Legendre moments equal to zero, the following correction for $\ell_x \ge 1$ is proposed:

$$\textbf{if} \quad \left|\sigma_{s,\ell_x}^{i}(E)\right| > \left|\sigma_{s,0}^{i}(E)\right| \quad \textbf{then} \quad \sigma_{s,\ell}^{i}(E) = 0 \quad \forall \ \ell \ge \ell_x \,,$$

which are based on the upper bound of high-order Legendre moments.

## 2.2.2.2 Total Cross-Sections

The elastic total cross-sections are simply given by

$$\Sigma_{t}^{\texttt{e±}}(E) = \Sigma_{s,0}^{\texttt{e±}\rightarrow\texttt{e±}}(E) \,.$$

## 2.2.2.3 Transport Correction and Elastic Decomposition in Soft and Catastrophic Components

The elastic scattering cross-sections can be decomposed into soft and catastrophic components [landesman1989angular](@cite). Let $L$ be the order of the cross-sections Legendre expansion and $\Sigma_{s,\ell,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}}$ be the $\ell$-order Legendre moment of the elastic scattering cross-section in group $g$. Let $L_{\texttt{max}} \le L$ be the last non-zero Legendre moments of the scattering in group $g$. The Legendre moments of the soft are assumed to be given by

$$\Sigma_{s,\ell,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}} = \Sigma_{s,0,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}}- T_{g}\ell(\ell+1)$$

for $\ell \in \left\{1,...,L_{\texttt{max}}\right\}$. This expression is obtained by establishing equality between the eigenvalues of the Boltzmann and the AFP operator applied to the Legendre polynomials \cite{morel1981fokker}. This method sets a relation, which depends on two undefined parameters $\Sigma_{s,0,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}}$ and $T_{g}$, such as moments of the Boltzmann operator are preserved by the Fokker-Planck operator. Landesman and Morel proposed to equate $\Sigma_{s,L_{\texttt{max}}-1}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}} = \Sigma_{s,L_{\texttt{max}}-1,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}}$ and $\Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}} = \Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}}$ to define these parameters, which results in

$$T_{g} = \frac{\Sigma_{s,L_{\texttt{max}}-1,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} - \Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}}}{2L_{\texttt{max}}} \,,$$

and

$$\Sigma_{s,0,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±},\texttt{soft}} = \Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} + T_{g}L_{\texttt{max}}(L_{\texttt{max}}+1)\,.$$

The soft component of the elastic cross-sections should then be withdrawn from cross-sections, since the AFP operator, jointly with the momentum transfer, is used to treat that soft scattering. The total elastic cross-sections are redefined as

$$\Sigma_{t,g} \leftarrow \Sigma_{t,g} - \Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} - T_{g}L_{\texttt{max}}(L_{\texttt{max}}+1) \,,$$

and the $\ell$-order Legendre moment of the scattering cross-sections are redefined as

$$\Sigma_{s,\ell,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} \leftarrow \Sigma_{s,\ell,g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} - \Sigma_{s,L_{\texttt{max}},g \rightarrow g}^{\texttt{e±}\rightarrow\texttt{e±}} - T_{g}\left[L_{\texttt{max}}(L_{\texttt{max}}+1)-\ell(\ell+1)\right] \,.$$

The goal of this operation is to transfer the handling of forward-peaked scattering from the Boltzmann operator, which encounters monotonicity issues with such scattering [azmy2010advances](@cite), to the Fokker-Planck operator, which can be tackled by finite-difference discretization schemes.

This method, as presented, includes the extended transport correction [morel1979validity,drumm2007analysis,hebert2009applied](@cite). The transport correction is a method that reduces the amplitude of the elastic cross-sections for the BFP solver while leaving intact the flux solution, which requires the use of Galerkin quadrature [azmy2010advances](@cite). This amplitude reduction aims to reduce the scattering ratio, which is very close to one in charged particle transport. Without this correction, the source iteration process would require an astounding number of iterations, thus a long time, to converge [morel1989hybrid](@cite).