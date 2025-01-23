# 1.1 The Steady-State Transport Equation

The linear Boltzmann transport equation (BTE) describes the statistical behavior of particles, like neutrons, photons, or electrons, as they move through and interact with a medium. It is derived based on the enforcement of particle conservation under a few assumptions: 1) the transported particles interact independently, without collective effects or correlations between them, 2) the medium are isotropic in space, and 3) the material and external sources properties does not evolve over time. The derivation of this transport equation can be found in many nuclear engineering textbooks, such as [lewis1984computational](@citet) or [hebert2009applied](@citet). For many applications, there is no evolution in time of the variables (i.e. sources or material properties) over time. In such case, the steady-state form of the transport equation is employed.

!!! note
    Radiant deals only with steady-state conditions. 

## 1.1.1 The steady-state Linear Boltzmann Transport Equation

Let E be the particle kinetic energy, $\mathbf{r} = (x,y,z)$ be the particle position in $\mathbb{R}^3$, $\mathbf{\Omega} = (\mu,\phi) = (\mu,\eta,\xi)$ be the moving particle direction, where $\mu$ is the principal direction cosine aligned with the x-axis, $\phi$ is the azimuthal angle, while $\eta$ and $\xi$ are the secondary direction cosines under the constraint $\mu^2+\eta^2+\xi^2 = 1$. Let

$$\mathbb{S}^2 = \left\{\mathbf{\Omega} = (\mu,\eta,\xi) \in \mathbb{R}^3 : \left|\left|\mathbf{\Omega}\right|\right| = 1\right\}$$

be the unit sphere. Let $P$ be the set of particles to be considered. The steady-state linear BTE for a particle $p \in P$ is given by

$$\mathbf{\Omega} \cdot \nabla \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) + \Sigma_{t}^{p}(\mathbf{r},E) \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) = Q_{p}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega},E) + Q_{p}^{\texttt{ext}}(\mathbf{r},\mathbf{\Omega},E) \,,$$

which solution is given by $\Phi_{p}(\mathbf{r},\mathbf{\Omega},E)$, the angular flux for particle $p$. It is defined as the product of the population density $n_{p}(\mathbf{r},\mathbf{\Omega},E)$ and the particle speed $v_{p}(E)$, i.e.

$$\Phi_{p}(\mathbf{r},\mathbf{\Omega},E) = v_{p}(E) n_{p}(\mathbf{r},\mathbf{\Omega},E) \,.$$

The variable $\Sigma_{t}^{p}(\mathbf{r},E)$ is the macroscopic total cross-section for particle $p$. There is also two sources terms, 1)
the external fixed volume source, $Q_{p}^{\texttt{ext}}(\mathbf{r},\mathbf{\Omega},E)$, which can represent a radioactive source within the medium for example, and the source associated with the Boltzmann operator, which account for both the particle scattering and production in the medium, that is given by

$$Q_{p}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega},E) = \sum_{p' \in P} \int_{\mathbb{S}^2}d^{2}\Omega'\int_{0}^{\infty}dE' \Sigma_{s}^{p'\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}'\rightarrow\mathbf{\Omega}) \Phi_{p'}(\mathbf{r},\mathbf{\Omega}',E') \,.$$

where $\Sigma_{s}^{p'\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}'\rightarrow\mathbf{\Omega})$ is the macroscopic double differential scattering cross-section describing the probability of a particle $p'$ with energy $E'$ and direction $\mathbf{\Omega}'$ to yield a particle $p$ with energy $E$ and direction $\mathbf{\Omega}$. Due to the assumed isotropic properties of the medium, the scattering cross-sections only dependent on the scattering angle $\mathbf{\Omega}\cdot\mathbf{\Omega}'$ such as the Boltzmann operator can be rewritten as

$$Q_{p}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega},E) = \frac{1}{2\pi}\sum_{p' \in P} \int_{\mathbb{S}^2}d^{2}\Omega'\int_{0}^{\infty}dE' \Sigma_{s}^{p'\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}\cdot\mathbf{\Omega}') \Phi_{p'}(\mathbf{r},\mathbf{\Omega}',E') \,.$$
 
The integrated flux for particle $p$ is given by integrating the angular flux over the unit sphere

$$\Phi_{p}(\mathbf{r},E) = \int_{\mathbb{S}^2}d^{2}\Omega \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) \,.$$

This integrated flux is useful to compute so-called reaction rates, such as energy or charge deposition.

!!! note
    Radiant only deals with medium exhibiting isotropic properties. 

## 1.1.2 The Fokker-Planck Equation

While successful with neutral particles such as neutrons and photons in reactor physics, the steady-state BTE becomes inefficient for charged particle transport. This failure is due to the highly forward-peaked scattering of charged particles, the large amount of scattering events due to interactions with atomic electrons in the medium, and the small amount of energy lost in each of these scattering events. An accurate discretization of the energy domain would require an astronomical number of energy meshes to capture these small energy losses effectively [lorence1989physics](@cite).

Under the assumption of forward-peaked and small energy-loss processes, it is possible to derive an alternate form of the Boltzmann operator. Indeed, using Taylor expansion in the direction cosine and energy variables, the Boltzmann operator reduces to the Fokker-Planck (FP) operator. The steady-state FP equation for a particle $p \in P$ is given by [morel1981fokker,pomraning1992fokker](@cite)

$$\mathbf{\Omega} \cdot \nabla \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) + \Sigma_{a}^{p}(\mathbf{r},E) \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) = Q_{p}^{\texttt{FP}}(\mathbf{r},\mathbf{\Omega},E) + Q_{p}^{\texttt{ext}}(\mathbf{r},\mathbf{\Omega},E) \,,$$

where $\Sigma_{a}^{p}(\mathbf{r},E)$ is the macroscopic absorption cross-section for particle $p$ and where the FP operator is given by

$$Q^{\texttt{FP}}_{p}(\mathbf{r},\mathbf{\Omega},E) = Q^{\texttt{CSD}}_{p}(\mathbf{r},\mathbf{\Omega},E) + Q^{\texttt{AFP}}_{p}(\mathbf{r},\mathbf{\Omega},E) \,,$$

where the continuous slowing-down (CSD) operator is given by

$$Q^{\texttt{CSD}}_{p}(\mathbf{r},\mathbf{\Omega},E) = \frac{\partial}{\partial E}\left[ S^{p}(\mathbf{r},E)\Phi_{p}(\mathbf{r},\mathbf{\Omega},E) \vphantom{\frac{\partial}{\partial E}} \right]$$

and the angular Fokker-Planck (AFP) operator is given by

$$Q^{\texttt{AFP}}_{p}(\mathbf{r},\mathbf{\Omega},E) = T^{p}(\mathbf{r},E) \left[\frac{\partial}{\partial \mu}\left(1-\mu^2\right)\frac{\partial}{\partial \mu} + \frac{1}{1-\mu^2}\frac{\partial^2}{\partial\phi^2}\right]\Phi_{p}(\mathbf{r},\mathbf{\Omega},E) \,.$$

where $S^{p}(\mathbf{r},E)$ is the stopping power and $T^{p}(\mathbf{r},E)$ is the momentum transfer, which are given by [pomraning1992fokker](@cite)

$$S^{p}(\mathbf{r},E) = \int_{\mathbb{S}^2}d^{2}\Omega' \int_{0}^{\infty}dE' (E-E') \Sigma_{s}^{p\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}\cdot\mathbf{\Omega}') \,,$$

$$T^{p}(\mathbf{r},E) = \int_{\mathbb{S}^2}d^{2}\Omega' \int_{0}^{\infty}dE' (1-\mu) \Sigma_{s}^{p\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}\cdot\mathbf{\Omega}') \,.$$

A higher-order FP term from the Taylor expansion is often kept, which is the energy straggling operator, and is important to obtain more accurate results [morel1981fokker,uilkema2012proton](@cite). Generalized FP operators also have been developed [pomraning1996higher,olbrant2010generalized](@cite). While this equation performs well with very fast or heavily charged particles [uilkema2012proton](@cite), electrons and positrons often change direction when interacting. These large deflections often come with large energy transfers, which cannot be ignored without making significant errors.

## 1.1.3 The Boltzmann Fokker-Planck Equation

[przybylski1982numerical](@citet) proposed to divide the scattering events into two domains: soft, which consists of small change of direction and energy-loss, and catastrophic, for large change of direction and energy-loss. The catastrophic interactions are treated with the Boltzmann operator, while the soft interactions are treated with the FP operator. The two domains are defined as a function of the particle energy $E$ by a threshold function $E_{c}(E) \le E$, over which interactions are considered soft and under which they are considered catastrophic. The steady-state Boltzmann Fokker-Planck (BFP) equation is given by

$$\mathbf{\Omega} \cdot \nabla \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) + \Sigma_{t}^{p}(\mathbf{r},E) \Phi_{p}(\mathbf{r},\mathbf{\Omega},E) = Q_{p}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega},E) + Q_{p}^{\texttt{FP}}(\mathbf{r},\mathbf{\Omega},E) + Q_{p}^{\texttt{ext}}(\mathbf{r},\mathbf{\Omega},E) \,.$$

The definition of the cross-sections, stopping powers and momentum transfers change from the BTE and FP equation defined previously. The total and scattering cross-sections are limited to catastrophic events, while the stopping powers and momentum transfers are limited to soft events [azmy2010advances](@cite).

It can be appreciated that the BFP can be seen as a generalization of both BTE and FP equation. The threshold function $E_{c}(E) \le E$ can be choosen, for a given interaction, as 1) $E_{c}(E) = E$, which reduces the BFP to the BTE equation, or 2) $E_{c}(E) = 0$, which reduces the BFP to the FP equation. Therefore, in the following sections, the discretization will be restricted to the BFP equation since it also encompasses the BTE and FP. Otherwise, the choice of $E_{c}(E)$ depends on the needs for a given interaction and will be discussed further in the following sections.

## 1.1.4 The Lorentz Force in Presence of External Electromagnetic Fields

!!! note
    Radiant does not have any capabilities for external electromagnetic fields yet.

## 1.1.5 The Adjoint Linear Transport Equation

!!! note
    Radiant does not containts any adjoint capabilities yet.