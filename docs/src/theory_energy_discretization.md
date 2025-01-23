# 1.2 Energy discretization

The energy domain of the transport equation is discretized using finite-element methods. This discretization defines the structure of the atomic data (cross-sections, stopping powers and momentum transfer) to be prepared for transport calculation.

## 1.2.1 Classical Galerkin Method of Weighted Residuals

```@raw html
<figure>
    <img src="../energy_groups.png" alt="drawing" style="width: min(550px, 100vw);"/>
    <figcaption> Figure 1: Energy Discretization.</figcaption>
</figure>
```

The energy domain is divided $G$ meshes, often refered to as energy group. The convention where $g = 1$ corresponds to the highest energy group is used. The midpoint energy of the group $g$ is $E_{g}$, its width is $\Delta E_{g}$, and its upper and lower boundaries are respectively $E_{g-1/2}$ and $E_{g+1/2}$, as shown in Fig.1. The lower energy boundary of the lowest energy group, $E_{G+1/2}$, is called the cutoff energy. Particles that scatter or are produced under this energy cutoff are considered to be absorbed locally. It also means that charge and energy are deposited locally. This energy should be chosen carefully, being low enough to ensure that the absorbed particle mean free path is lower than the spatial meshes' width to be realistic.

The classical multigroup method can be applied assuming that the angular flux is separable in energy, as describe in [lewis1984computational](@citet). It also can be developed as a special case of the classical Galerkin method of weighted residuals [hesthaven2007nodal](@cite) applied to the energy domain. The change of variable

$$u = \frac{2E-E_{g-1/2}-E_{g+1/2}}{2\Delta E_{g}} \in [-1/2,1/2]$$

can be applied to the energy variable $E$ to simplify the application of the method. The BFP equation then becomes

$$\begin{split} &\underbrace{\mathbf{\Omega} \cdot \nabla \Phi_{p}(\mathbf{r},\mathbf{\Omega},u)}_{\mathsf{A}} + \underbrace{\Sigma_{t}^{p}(\mathbf{r},u) \Phi_{p}(\mathbf{r},\mathbf{\Omega},u)}_{\mathsf{B}} + \underbrace{\frac{-1}{\Delta E_{g}} \frac{\partial}{\partial u}\left[ S^{p}(\mathbf{r},u)\Phi_{p}(\mathbf{r},\mathbf{\Omega},u) \vphantom{\frac{\partial}{\partial u}} \right]}_{\mathsf{C}} \\ &= \underbrace{Q_{p}^{\texttt{AFP}}(\mathbf{r},\mathbf{\Omega},u)}_{\mathsf{D}} + \underbrace{Q_{p}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega},u)}_{\mathsf{E}} + \underbrace{Q_{p}^{\texttt{ext}}(\mathbf{r},\mathbf{\Omega},u)}_{\mathsf{F}} \,,\end{split}$$

where the energy derivative in the CSD operator have also been affected by this change of variable. The angular flux in group $g$ is expanded up to the $M$ order, namely

$$\Phi_{p,g}(\mathbf{r},\mathbf{\Omega},u) = \sum_{m=0}^{M} \tilde{P}_{m}(u) \Phi_{p,g}^{(m)}(\mathbf{r},\mathbf{\Omega}) \,,$$

such as the moment of the angular flux are given by

$$\Phi_{p,g}^{(m)}(\mathbf{r},\mathbf{\Omega}) = \int_{-1/2}^{1/2}du \tilde{P}_{m}(u) \Phi_{p,g}(\mathbf{r},\mathbf{\Omega},u) \,.$$

Nonetheless, the total cross-sections, the scattering cross-sections, the momentum transfers and the stopping powers are also expanded up to the $M$ order,

$$\Sigma_{t,g}^{p}(\mathbf{r},u) = \sum_{m=0}^{M} \tilde{P}_{m}(u) \Sigma_{t,g}^{p,(m)}(\mathbf{r})\,,$$

$$\Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}',u,u') = \frac{1}{\Delta E_{g'}} \sum_{m=0}^{M} \sum_{m'=0}^{M} \tilde{P}_{m}(u) \tilde{P}_{m'}(u') \Sigma_{s,g'\rightarrow g}^{p'\rightarrow p,(m'\rightarrow m)}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') \,,$$

$$T_{g}^{p}(\mathbf{r},u) = \sum_{m=0}^{M} \tilde{P}_{m}(u) T_{g}^{p,(m)}(\mathbf{r})\,,$$

$$S_{g}^{p}(\mathbf{r},u) = \sum_{m=0}^{M} \tilde{P}_{m}(u) S_{g}^{p,(m)}(\mathbf{r})\,,$$

such as their $m$-th moments are given respectively by

$$\Sigma_{t,g}^{p,(m)}(\mathbf{r}) = \int_{-1/2}^{1/2}du \tilde{P}_{m}(u) \Sigma_{t,g}^{p}(\mathbf{r},u)\,,$$

$$\Sigma_{s,g'\rightarrow g}^{p'\rightarrow p,(m'\rightarrow m)}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') = \Delta E_{g'} \int_{-1/2}^{1/2}du' \tilde{P}_{m'}(u')\int_{-1/2}^{1/2}du \tilde{P}_{m}(u) \Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}',u,u') \,,$$

$$T_{g}^{p,(m)}(\mathbf{r}) = \int_{-1/2}^{1/2}du \tilde{P}_{m}(u) T_{g}^{p}(\mathbf{r},u)\,,$$

$$S_{g}^{p,(m)}(\mathbf{r}) = \int_{-1/2}^{1/2}du \tilde{P}_{m}(u) S_{g}^{p}(\mathbf{r},u)\,.$$

The goal is to express the BFP equation as a function of these moments, which can be calculated beforehand. Multiplying the BFP equation by the normalized Legendre polynomials $P_{m}(u)$ and integrating over $u \in [-1/2,1/2]$, the $m$-th moment of the BFP equation, term-per-term, become 

$$\left< \mathsf{A} \right> = \mathbf{\Omega} \cdot \nabla \Phi_{p,g}^{(m)}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\left<  \mathsf{B} \right> = \sqrt{2m+1} \sum_{j_{2}=0}^{M} \sum_{j_{3}=0}^{M} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{m,j_2,j_3} \Sigma_{t,g}^{p,(j_{2})}(\mathbf{r}) \Phi_{p,g}^{(j_3)}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\begin{split} \left<  \mathsf{C} \right> &= -\frac{\sqrt{2m+1}}{\Delta E_{g}} \left\{ S^{p}_{g-1/2}(\mathbf{r})\Phi_{p,g-1/2}(\mathbf{r},\mathbf{\Omega}) - (-1)^{m} S^{p}_{g+1/2}(\mathbf{r})\Phi_{p,g+1/2}(\mathbf{r},\mathbf{\Omega}) \right. \\ &\left.- \sum_{j_{1}=0}^{m-1} \frac{2j_{1}+1}{2} \left[ 1 - (-1)^{m-j_{1}-1} \right] \sum_{j_{2}=0}^{M} \sum_{j_3=0}^{M} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{j_1,j_2,j_3} S_{g}^{p,(j_{2})}(\mathbf{r})\Phi_{p,g}^{(j_3)}(\mathbf{r},\mathbf{\Omega}) \right\} \,,\end{split}$$

$$\left<  \mathsf{D} \right> = \sqrt{2m+1} \sum_{j_{2}=0}^{M} \sum_{j_{3}=0}^{M} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{m,j_2,j_3} T_{g}^{p,(j_{2})}(\mathbf{r}) \nabla^{2}_{\Omega}\Phi_{p,g}^{(j_3)}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\left<  \mathsf{E} \right> = \frac{1}{2\pi} \sum_{p' \in P} \int_{\mathbb{S}^2}d^{2}\Omega' \sum_{g'=1}^{G+1} \sum_{m'=0}^{M} \Sigma_{s,g'\rightarrow g}^{p'\rightarrow p,(m'\rightarrow m)}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') \Phi_{p',g'}^{(m')}(\mathbf{r},\mathbf{\Omega}') \,,$$

and $\left<\mathsf{F}\right>$ is calculated depending on the definition of the fixed external source inside the medium volume. The function $\mathcal{W}_{j_1,j_2,j_3}$ is defined as

$$\mathcal{W}_{j_1,j_2,j_3} = \frac{1}{2}\int_{-1}^{1}dx P_{j_1}(x) P_{j_2}(x) P_{j_3}(x) \,.$$

The variables $S^{p}_{g\pm1/2}(\mathbf{r})$ and $\Phi_{p,g\pm1/2}(\mathbf{r},\mathbf{\Omega})$ are the stopping powers and angular flux at the higher (index g-1/2) and lower (index g+1/2) boundaries. These incoming and outgoing flux at energy group boundaries requires to define a closure relation, which will be defined later.

!!! warning
    RADIANT does not include yet all the capabilities for generalized ($M > 0$) Galerkin discretization in energy. It includes generalization for any $M$ only for the CSD term, while also assuming that the stopping powers moments are given by 
    
    $$S_{g}^{p,(0)}(\mathbf{r}) = \frac{1}{2} \left[ S^{p}_{g-1/2}(\mathbf{r}) + S^{p}_{g+1/2}(\mathbf{r})\right]\,,$$

    $$S_{g}^{p,(1)}(\mathbf{r}) = \frac{1}{2\sqrt{3}} \left[ S^{p}_{g-1/2}(\mathbf{r}) - S^{p}_{g+1/2}(\mathbf{r})\right] \,,$$

    and

    $$S_{g}^{p,(m)}(\mathbf{r}) = 0 \quad \forall \, m \ge 2 \,.$$

    For all other terms, the following multigroup discretisation is used.


## 1.2.2 Multigroup Discretization

The multigroup discetization is the most employed energy discretization method. It is assumed that the particles define in the same energy group share the same energy and averaged energy-dependent variables. The classical multigroup method can be obtained assuming that the angular flux is separable in energy, as describe in [lewis1984computational](@citet). It also correspond to taking the Galerkin method of weighted residuals up to order $M = 0$. The group-averaged angular flux in group $g$ is therefore given by

$$\Phi_{p,g}(\mathbf{r},\mathbf{\Omega},E) = \Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) \,,$$

with

$$\Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) = \frac{1}{\Delta E_{g}} \int_{E_{g+1/2}}^{E_{g-1/2}}dE \Phi_{p,g}(\mathbf{r},\mathbf{\Omega},E) \,,$$

where the order index $(0)$ is and will omitted in the following section. The group-averaged total cross-sections, scattering cross-sections, momentum transfers and stopping powers are then given by

$$\Sigma_{t,g}^{p}(\mathbf{r}) = \frac{1}{\Delta E_{g}}\int_{E_{g+1/2}}^{E_{g-1/2}}dE \Sigma_{t,g}^{p}(\mathbf{r},E) \,,$$

$$\Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') = \frac{1}{\Delta E_{g}} \int_{E_{g+1/2}}^{E_{g-1/2}}dE \int_{E_{g'+1/2}}^{E_{g'-1/2}}dE' \Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},E'\rightarrow E,\mathbf{\Omega}\cdot\mathbf{\Omega}') \,,$$

$$T_{g}^{p}(\mathbf{r}) = \frac{1}{\Delta E_{g}}\int_{E_{g+1/2}}^{E_{g-1/2}}dE  T_{g}^{p}(\mathbf{r},E) \,,$$

$$S_{g}^{p}(\mathbf{r}) = \frac{1}{\Delta E_{g}}\int_{E_{g+1/2}}^{E_{g-1/2}}dE  S_{g}^{p}(\mathbf{r},E) \,.$$

Because of the energy derivative, transport calculations also requires boundary values of the stopping powers, namely

$$S_{g-1/2}^{p}(\mathbf{r}) = S^{p}(\mathbf{r},E_{g-1/2}) \quad\text{and}\quad S_{g+1/2}^{p}(\mathbf{r}) = S^{p}(\mathbf{r},E_{g+1/2}) \,.$$

The zeroth-th moment of the BFP equation, term-per-term, is given by

$$\left< \mathsf{A} \right> = \mathbf{\Omega} \cdot \nabla \Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\left<  \mathsf{B} \right> = \Sigma_{t,g}^{p}(\mathbf{r}) \Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\begin{split} \left<  \mathsf{C} \right> &= -\frac{1}{\Delta E_{g}} \left\{ S^{p}_{g-1/2}(\mathbf{r})\Phi_{p,g-1/2}(\mathbf{r},\mathbf{\Omega}) - S^{p}_{g+1/2}(\mathbf{r})\Phi_{p,g+1/2}(\mathbf{r},\mathbf{\Omega}) \right\} \,,\end{split}$$

$$\left<  \mathsf{D} \right> = T_{g}^{p}(\mathbf{r}) \nabla^{2}_{\Omega}\Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) \,,$$

$$\left<  \mathsf{E} \right> = \frac{1}{2\pi} \sum_{p' \in P} \int_{\mathbb{S}^2}d^{2}\Omega' \sum_{g'=1}^{G+1} \Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') \Phi_{p',g'}(\mathbf{r},\mathbf{\Omega}') \,,$$

and $\left<\mathsf{F}\right>$ is calculated depending on the definition of the fixed external source inside the medium volume.
 
## 1.2.3 Energy Group Structure

The energy domain can be discretized using three pieces of information: the number of groups, $G$, the mid-point energy of the highest energy group, $E_1$, and the cutoff energy $E_{G+1/2}$.

### 1.2.3.1 Linear Structure

The linear structure is characterized by same-size energy group. The mathbftor containing its energy boundaries is given by
   
$$\mathbf{E} = \left\{ E_{G+1/2} + \frac{k(E_{1/2}-E_{G+1/2})}{G} \mid k = 0,1,2,...,G \right\}\,,$$

where the higher bound of the maximum energy group is given by

$$E_{1/2} = E_{1} + \frac{E_{1}-E_{G+1/2}}{2G+1} \,.$$

### 1.2.3.2 Logarithmic Structure

The logarithmic structure is given by

$$\mathbf{E} = \left\{ 10^{\log_{10}{E_{G+1/2}} + \frac{k(\log_{10}{E_{1/2}}-\log_{10}{E_{G+1/2}})}{G}} \mid k = 0,1,2,...,G \right\}\,.$$

The higher bound of the maximum energy group is given by solving

$$2E_{1} - E_{1/2} - 10^{\log_{10}{E_{1/2}} - \left(\log_{10}{E_{1/2}}-\log_{10}{E_{G+1/2}}\right)/G} = 0$$

using a root-finding algorithm.

## 1.2.4 Definition of soft and catastrophic interactions

```@raw html
<figure>
    <img src="../soft_catas_cutoff.jpg" alt="drawing" style="width: min(750px, 100vw);"/>
    <figcaption> Figure 2: Definition of soft and catastrophic interactions.</figcaption>
</figure>
```

The definition of soft and catastrophic domains on the energy spectrum at group boundaries are shown in Fig.2, where two points $(E,E_{c}(E))$ are defined as $(E_{g-1/2},E_{g+1/2})$ and $(E_{g+1/2},E_{g+3/2})$ for any group $g$, with the additional definition $E_{G+3/2} = 0$. A smooth transition between these two coordinates is required and is set to be linear in $E$. The resulting energy dividing the soft and catastrophic domains, for $E \in [E_{g-1/2},E_{g+1/2}]$, is given by

$$E_{c}(E) = \left(\frac{E_{g+1/2}-E_{g+3/2}}{E_{g-1/2}-E_{g+1/2}}\right) E - \frac{E_{g+1/2}^2-E_{g-1/2}E_{g+3/2}}{E_{g-1/2}-E_{g+1/2}}$$

for any $1 \le g \le G$.