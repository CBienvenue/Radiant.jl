# 1.5 Closure Relations

The spatial- and energy-discretized BFP equation (Section 1.4) is written in terms of the volume-averaged moments $\Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})}$ and the boundary moments $\Phi_{p,g\pm 1/2,n,i,j,k}^{(m_{x},m_{y},m_{z})}$, $\Phi_{p,g,n,i\pm 1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})}$ and analogous terms along $y$ and $z$. Closure relations link the *outgoing* boundary moments of a voxel (or energy group) to the volume and *incoming* boundary moments solved in the cell. They are the central ingredient that makes the sweeping process (Section 1.6) solvable on a voxel-by-voxel basis.

Radiant supports several families of closure relations selectable per axis with `sn.set_scheme(axis,scheme,order)`. They all take the generic form, along an axis $x$ for the $\nu$-order moment expansion of order $M_{x} = \mathcal{O}_{x}-1$,

$$\Phi_{p,g,n,i+1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} = (-1)^{M_{x}+1} \omega^{x}_{0}\,\Phi_{p,g,n,i-1/2,j,k}^{(m_{\texttt{E}},m_{y},m_{z})} + \sum_{m_{x}=0}^{M_{x}} \omega^{x}_{m_{x}+1}\,\Phi_{p,g,n,i,j,k}^{(m_{\texttt{E}},m_{x},m_{y},m_{z})} \,,$$

where the weights $\omega^{x}_{0},\ldots,\omega^{x}_{M_{x}+1}$ characterize the scheme. The same structure is reused, with sign conventions appropriate to the upwind direction, along $y$, $z$ and $\texttt{E}$.

## 1.5.1 Diamond Difference (DD)

The diamond-difference closure of order $\mathcal{O}_{x}$ corresponds to

$$\omega^{x}_{0} = (-1)^{M_{x}+1} \,,\qquad \omega^{x}_{m_{x}+1} = 1 - (-1)^{M_{x}+1-m_{x}}\,, \quad m_{x}=0,\ldots,M_{x} \,.$$

It is a generalization of the classical diamond-difference scheme (recovered at order 1) to higher polynomial orders [bienvenue2022high](@cite). The scheme is exact for polynomials up to degree $M_{x}+1$ but can produce negative outgoing moments in optically thick voxels.

## 1.5.2 Discontinuous Galerkin (DG, DG+, DG-)

The DG closure does not enforce continuity of the angular flux between voxels. The outgoing boundary moments are expressed as the trace of the interior polynomial expansion. In Radiant's normalized Legendre basis this is

$$\omega^{x}_{0} = 0 \,,\qquad \omega^{x}_{m_{x}+1} = 1\,, \quad m_{x}=0,\ldots,M_{x} \,.$$

Two upwind-biased variants of DG are available:

- **DG-**: identical to DG except for the highest-order coefficient, which is set to $\omega^{x}_{M_{x}+1} = M_{x}/(2M_{x}+1)$.
- **DG+**: identical to DG except for the highest-order coefficient, which is set to a large value (numerically saturating the slope), effectively forcing the outgoing flux to follow the highest-order interior moment.

DG, DG+ and DG- require an expansion of order $\mathcal{O}_{x} \ge 2$ (i.e. $M_{x} \ge 1$).

## 1.5.3 Adaptive Weighted Difference (AWD)

The adaptive-weighted-difference scheme [voloschenko2011some,bienvenue2022high](@cite) selects $\omega^{x}_{0}$ at each voxel as the value that maximizes the local accuracy while keeping the outgoing moments non-negative. The scheme uses the DD closure as a starting point and falls back to a more diffusive choice (closer to DG) when DD would produce a negative outgoing moment. The adaptive procedure is purely local — it never requires an outer iteration — and is particularly effective when traversing optically thick voxels or sharp material discontinuities.

When any axis of the calculation uses AWD, the adaptive flag is propagated to every other axis: the calculation as a whole becomes adaptive. AWD is currently restricted to 1st-order ($\mathcal{O}=1$) and 2nd-order ($\mathcal{O}=2$) expansions.

## 1.5.4 Multidimensional and Energy Closures

In 2D and 3D the closure relations are applied independently along each axis. Two coupling strategies are available for the high-order moments of the volume flux (Section 6.2.5):

- **Full coupling** (`set_is_full_coupling(true)`, default): the volume expansion contains every monomial $\tilde{P}_{m_{x}}(u_{x})\tilde{P}_{m_{y}}(u_{y})\tilde{P}_{m_{z}}(u_{z})$ with $0\le m_{x}\le M_{x}$, $0\le m_{y}\le M_{y}$, $0\le m_{z}\le M_{z}$. The total number of volume moments is $\mathcal{O}_{x}\mathcal{O}_{y}\mathcal{O}_{z}$.
- **Partial coupling**: only the diagonal and the axis-aligned moments are retained, giving $1+\sum_{i}(\mathcal{O}_{i}-1)$ volume moments. This is a cheaper but less accurate variant.

The same closure relations apply to the continuous slowing-down derivative in energy. The energy axis is treated as an additional direction with its own scheme and order (`set_scheme("E",scheme,order)`), and the same weighting machinery is reused.

The discrete closure relations described here are implemented in `scheme_weights.jl` and consumed by the sweep routines `sn_sweep_1D.jl`, `sn_sweep_2D.jl` and `sn_sweep_3D.jl` (Section 1.6).
