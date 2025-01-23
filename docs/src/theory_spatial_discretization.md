# 1.4 Spatial discretization

The spatial domain is discretized using finite-element methods.

## 1.4.1 Cartesian geometry

```@raw html
<figure>
    <img src="../spatial_meshes.png" alt="drawing" style="width: min(550px, 100vw);"/>
    <figcaption> Figure 1: Spatial Discretization.</figcaption>
</figure>
```

The spatial domain is divided $N_{x}$ meshes along the $x$ axis, $N_{y}$ meshes along the $y$ axis and $N_{z}$ meshes along the $z$ axis. The midpoint energy of the mesh $i$ ($j$ or $k$) along the $x$ ($y$ or $z$) axis is $x_{i}$ ($y_{j}$ or $z_{k}$), its width is $\Delta x_{i}$ ($\Delta y_{j}$ or $\Delta z_{k}$), and its upper and lower boundaries are respectively $x_{i+1/2}$ ($y_{j+1/2}$ or $z_{k+1/2}$) and $x_{i-1/2}$ ($y_{j-1/2}$ or $z_{k-1/2}$), as shown in Fig.1. The classical Galerkin method of weighted residuals [hesthaven2007nodal](@cite) is applied to the spatial domain. The changes of variable

$$u_{x} = \frac{2x-x_{i-1/2}-x_{i+1/2}}{2\Delta x_{i}} \in [-1/2,1/2] \,,$$

$$u_{y} = \frac{2y-y_{j-1/2}-y_{j+1/2}}{2\Delta y_{j}} \in [-1/2,1/2] \,,$$

$$u_{z} = \frac{2z-z_{k-1/2}-z_{k+1/2}}{2\Delta z_{k}} \in [-1/2,1/2] \,.$$

is applied. The atomic data are assumed to be constant within each spatial element, such as the resulting discrete ordinates multigroup BFP equation is given by

$$\begin{split} &\frac{\mu_{n}}{\Delta x_{i}} \frac{\partial \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}})}}{\partial u_{x}}(u_{x},u_{y},u_{z}) + \frac{\eta_{n}}{\Delta y_{j}} \frac{\partial \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}})}}{\partial u_{y}}(u_{x},u_{y},u_{z}) + \frac{\xi_{n}}{\Delta z_{k}} \frac{\partial \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}})}}{\partial u_{z}}(u_{x},u_{y},u_{z}) \\ &+ \sqrt{2m_{\texttt{E}}+1} \sum_{j_{2}=0}^{M_{\texttt{E}}} \sum_{j_{3}=0}^{M_{\texttt{E}}} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{m_{\texttt{E}},j_2,j_3} \Sigma_{t,g,i,j,k}^{p,(j_{2})} \Phi_{p,g,n,i,j,k}^{(j_3)}(u_{x},u_{y},u_{z}) \\ &- \frac{\sqrt{2m_{\texttt{E}}+1}}{\Delta E_{g}} \left\{ S^{p}_{g-1/2,i,j,k}\Phi_{p,g-1/2,n,i,j,k}(u_{x},u_{y},u_{z}) - (-1)^{m_{\texttt{E}}} S^{p}_{g+1/2,i,j,k}\Phi_{p,g+1/2,n,i,j,k}(u_{x},u_{y},u_{z}) \right. \\ &\left.- \sum_{j_{1}=0}^{m_{\texttt{E}}-1} \frac{2j_{1}+1}{2} \left[ 1 - (-1)^{m_{E}-j_{1}-1} \right] \sum_{j_{2}=0}^{M_{\texttt{E}}} \sum_{j_3=0}^{M_{\texttt{E}}} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{j_1,j_2,j_3} S_{g,i,j,k}^{p,(j_{2})}\Phi_{p,g,n,i,j,k}^{(j_3)}(u_{x},u_{y},u_{z}) \right\} \\ &= \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}})}(u_{x},u_{y},u_{z}) \end{split}$$

where $\Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}})}(u_{x},u_{y},u_{z})$ includes the Boltzmann, the angular Fokker-Planck and the external sources. The angular flux are expanded up the $(M_{\texttt{E}},M_{x},M_{y},M_{z})$ order such as

$$\Phi_{p,g,n,i,j,k}(u_{E},u_{x},u_{y},u_{z}) = \sum_{m_{\texttt{E}}=0}^{M_{\texttt{E}}} \sum_{m_{x}=0}^{M_{x}} \sum_{m_{y}=0}^{M_{y}} \sum_{m_{z}=0}^{M_{z}} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})} \,,$$

while the group boundary fluxes are given by an expansion up to the $(M_{x},M_{y},M_{z})$ order,

$$\Phi_{p,g\pm1/2,n,i,j,k}(u_{x},u_{y},u_{z}) = \sum_{m_{x}=0}^{M_{x}} \sum_{m_{y}=0}^{M_{y}} \sum_{m_{z}=0}^{M_{z}} \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g\pm1/2,n,i,j,k}^{(m_{x},m_{y},m_{z})} \,,$$

the spatial mesh x-boundary fluxes are given by an expansion up to the $(M_{\texttt{E}},M_{y},M_{z})$ order,

$$\Phi_{p,g,n,i\pm1/2,j,k}(u_{\texttt{E}},u_{y},u_{z}) = \sum_{m_{\texttt{E}}=0}^{M_{\texttt{E}}} \sum_{m_{y}=0}^{M_{y}} \sum_{m_{z}=0}^{M_{z}} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i\pm1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} \,,$$

the spatial mesh y-boundary fluxes are given by an expansion up to the $(M_{\texttt{E}},M_{x},M_{z})$ order,

$$\Phi_{p,g,n,i,j\pm1/2,k}(u_{\texttt{E}},u_{x},u_{z}) = \sum_{m_{\texttt{E}}=0}^{M_{\texttt{E}}} \sum_{m_{x}=0}^{M_{x}} \sum_{m_{z}=0}^{M_{z}} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i,j\pm1/2,k}^{(m_{\texttt{E}},m_{x},m_{z})} \,,$$

the spatial mesh z-boundary fluxes are given by an expansion up to the $(M_{\texttt{E}},M_{x},M_{y})$ order,

$$\Phi_{p,g,n,i,j,k\pm1/2}(u_{\texttt{E}},u_{x},u_{y}) = \sum_{m_{\texttt{E}}=0}^{M_{\texttt{E}}} \sum_{m_{x}=0}^{M_{x}} \sum_{m_{y}=0}^{M_{y}} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \Phi_{p,g,n,i,j,k\pm1/2}^{(m_{\texttt{E}},m_{x},m_{y})} \,.$$

The moment of the angular flux are given by

$$\Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})} = \int_{-1/2}^{1/2}du_{\texttt{E}} \int_{-1/2}^{1/2}du_{x} \int_{-1/2}^{1/2}du_{y} \int_{-1/2}^{1/2}du_{z} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i,j,k}(u_{\texttt{E}},u_{x},u_{y},u_{z}) \,.$$

while the boundary moments are given by

$$\Phi_{p,g\pm1/2,n,i,j,k}^{(m_{x},m_{y},m_{z})} = \int_{-1/2}^{1/2}du_{x} \int_{-1/2}^{1/2}du_{y} \int_{-1/2}^{1/2}du_{z} \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g\pm1/2,n,i,j,k}(u_{x},u_{y},u_{z}) \,,$$

$$\Phi_{p,g,n,i\pm1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} = \int_{-1/2}^{1/2}du_{\texttt{E}} \int_{-1/2}^{1/2}du_{y} \int_{-1/2}^{1/2}du_{z} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{y}}(u_{y}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i\pm1/2,j,k}(u_{\texttt{E}},u_{y},u_{z}) \,,$$

$$\Phi_{p,g,n,i,j\pm1/2,k}^{(m_{\texttt{E}},m_{x},m_{z})} = \int_{-1/2}^{1/2}du_{\texttt{E}} \int_{-1/2}^{1/2}du_{x} \int_{-1/2}^{1/2}du_{z} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{z}}(u_{z}) \Phi_{p,g,n,i,j\pm1/2,k}(u_{\texttt{E}},u_{x},u_{z}) \,,$$

$$\Phi_{p,g,n,i,j,k\pm1/2}^{(m_{\texttt{E}},m_{x},m_{y})} = \int_{-1/2}^{1/2}du_{\texttt{E}} \int_{-1/2}^{1/2}du_{x} \int_{-1/2}^{1/2}du_{y} \tilde{P}_{m_{\texttt{E}}}(u_{\texttt{E}}) \tilde{P}_{m_{x}}(u_{x}) \tilde{P}_{m_{y}}(u_{y}) \Phi_{p,g,n,i,j,k\pm1/2}(u_{\texttt{E}},u_{x},u_{y}) \,.$$

Multiplying the BFP equation by the normalized Legendre polynomials $P_{m_{x}}(u)$ and integrating over $u_{x} \in [-1/2,1/2]$, the $m_{x}$-th moment equation become

$$\begin{split} & \sqrt{2m_{x}+1}\frac{\mu_{n}}{\Delta x_{i}} \left[ \Phi_{p,g,n,i+1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} - (-1)^{m_{x}} \Phi_{p,g,n,i-1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} - \sum_{j_{1}=0}^{m_{x}-1} \sqrt{2j_{1}+1} \left[ 1 - (-1)^{m_{x}-j_{1}-1} \right] \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},j_{1},m_{y},m_{z})} \right] \\ &+ \sqrt{2m_{y}+1}\frac{\eta_{n}}{\Delta y_{j}} \left[ \Phi_{p,g,n,i,j+1/2,k}^{(m_{\texttt{E}},m_{x},m_{z})} - (-1)^{m_{y}} \Phi_{p,g,n,i,j-1/2,k}^{(m_{\texttt{E}},m_{x},m_{z})} - \sum_{j_{1}=0}^{m_{y}-1} \sqrt{2j_{1}+1} \left[ 1 - (-1)^{m_{y}-j_{1}-1} \right] \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},j_{1},m_{z})} \right] \\ &+ \sqrt{2m_{z}+1}\frac{\xi_{n}}{\Delta z_{k}} \left[ \Phi_{p,g,n,i,j,k+1/2}^{(m_{\texttt{E}},m_{x},m_{y})} - (-1)^{m_{z}} \Phi_{p,g,n,i,j,k-1/2}^{(m_{\texttt{E}},m_{x},m_{y})} - \sum_{j_{1}=0}^{m_{z}-1} \sqrt{2j_{1}+1} \left[ 1 - (-1)^{m_{z}-j_{1}-1} \right] \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},j_{1})} \right]\\ &+ \sqrt{2m_{\texttt{E}}+1} \sum_{j_{2}=0}^{M_{\texttt{E}}} \sum_{j_{3}=0}^{M_{\texttt{E}}} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{m_{\texttt{E}},j_2,j_3} \Sigma_{t,g,i,j,k}^{p,(j_{2})} \Phi_{p,g,n,i,j,k}^{(j_3,m_{x},m_{y},m_{z})} \\ &- \frac{\sqrt{2m_{\texttt{E}}+1}}{\Delta E_{g}} \left\{ S^{p}_{g-1/2,i,j,k}\Phi_{p,g-1/2,n,i,j,k}^{(m_{x},m_{y},m_{z})} - (-1)^{m_{\texttt{E}}} S^{p}_{g+1/2,i,j,k}\Phi_{p,g+1/2,n,i,j,k}^{(m_{x},m_{y},m_{z})} \right. \\ &\left.- \sum_{j_{1}=0}^{m_{\texttt{E}}-1} \frac{2j_{1}+1}{2} \left[ 1 - (-1)^{m_{E}-j_{1}-1} \right] \sum_{j_{2}=0}^{M_{\texttt{E}}} \sum_{j_3=0}^{M_{\texttt{E}}} \sqrt{2j_2+1}\sqrt{2j_3+1} \mathcal{W}_{j_1,j_2,j_3} S_{g,i,j,k}^{p,(j_{2})}\Phi_{p,g,n,i,j,k}^{(j_3,m_{x},m_{y},m_{z})} \right\} \\ &= \Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})} \end{split}$$

The incoming and outgoing flux at energy group and spatial mesh boundaries require the definition of closure relations, which is describe in a further section.