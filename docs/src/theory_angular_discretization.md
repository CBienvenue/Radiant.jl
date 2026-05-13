# 1.3 Angular Discretization

The angular variable $\mathbf{\Omega} = (\mu,\eta,\xi) \in \mathbb{S}^2$ of the BFP equation is discretized in two complementary ways in Radiant:

1. the **Boltzmann operator** is discretized on a quadrature set $\{(\mathbf{\Omega}_n,w_n)\}_{n=1}^{N_d}$ ("discrete ordinates");
2. the **angular Fokker–Planck operator** is discretized either on the same quadrature set using finite differences or differential quadrature, or in a moment basis (Galerkin scheme).

The discrete ordinates and Galerkin treatments share the same machinery — a transformation between *discrete* angular values and *Legendre / spherical-harmonics moments* — and so are presented in a unified way below.

## 1.3.1 The Discrete Ordinates Method

Letting the angular flux be evaluated at the $N_d$ quadrature directions,

$$\Phi_{p,g,n}(\mathbf{r}) \equiv \Phi_{p,g}(\mathbf{r},\mathbf{\Omega}_n) \,,$$

the streaming and total cross-section terms (terms $\mathsf{A}$ and $\mathsf{B}$ of the multigroup equation) reduce to

$$\mathbf{\Omega}_n \cdot \nabla \Phi_{p,g,n}(\mathbf{r}) + \Sigma_{t,g}^{p}(\mathbf{r}) \Phi_{p,g,n}(\mathbf{r}) \,,$$

which is a system of $N_d$ coupled first-order PDEs. The scattering source (term $\mathsf{E}$) requires the angular flux moments and is treated using the polynomial expansion of Section 1.3.2.

The quadrature set must be chosen with care: it must integrate the angular flux moments exactly up to the truncation order, and it must be compatible with the symmetry of the geometry. The quadratures supported by Radiant, their construction and their applicability are described in Section 3.2 of the Theory chapter.

## 1.3.2 Legendre and Spherical Harmonics Expansion of the Scattering Source

Under the isotropic-medium assumption, the differential scattering cross-section depends on the directions only through the cosine $\mathbf{\Omega}\cdot\mathbf{\Omega}'$. The Legendre addition theorem gives

$$\Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mathbf{\Omega}\cdot\mathbf{\Omega}') = \frac{1}{2\pi} \sum_{\ell=0}^{L} \frac{2\ell+1}{2} \Sigma_{s,\ell,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r}) \sum_{m=-\ell}^{\ell} R_{\ell}^{m}(\mathbf{\Omega}) R_{\ell}^{m}(\mathbf{\Omega}') \,,$$

where $R_{\ell}^{m}(\mathbf{\Omega})$ are the real spherical harmonics (Section 3.1.4) and $L$ is the Legendre truncation order. The corresponding Legendre moment of the scattering cross-section is

$$\Sigma_{s,\ell,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r}) = 2\pi \int_{-1}^{1} d\mu_{0}\, P_{\ell}(\mu_{0})\,\Sigma_{s,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r},\mu_{0}) \,.$$

The angular flux moments are

$$\Phi_{p,g}^{\ell,m}(\mathbf{r}) = \int_{\mathbb{S}^{2}} d^{2}\Omega\, R_{\ell}^{m}(\mathbf{\Omega})\, \Phi_{p,g}(\mathbf{r},\mathbf{\Omega}) \,,$$

so that the scattering source becomes

$$Q_{p,g}^{\texttt{B}}(\mathbf{r},\mathbf{\Omega}) = \frac{1}{2\pi}\sum_{p' \in P}\sum_{g'=1}^{G}\sum_{\ell=0}^{L} \frac{2\ell+1}{2} \Sigma_{s,\ell,g'\rightarrow g}^{p'\rightarrow p}(\mathbf{r})\sum_{m=-\ell}^{\ell} R_{\ell}^{m}(\mathbf{\Omega}) \Phi_{p',g'}^{\ell,m}(\mathbf{r}) \,.$$

In 1D Cartesian geometry, only the axial-symmetric moments ($m=0$) are needed and the spherical harmonics reduce to the Legendre polynomials $R_{\ell}^{0}(\mathbf{\Omega}) = P_{\ell}(\mu)$.

## 1.3.3 Moment-to-Discrete and Discrete-to-Moment Matrices

Let $\mathbf{\Phi}_{p,g}^{\bullet}(\mathbf{r}) = \left[\Phi_{p,g,1}(\mathbf{r}),\ldots,\Phi_{p,g,N_{d}}(\mathbf{r})\right]^{T}$ be the vector of discrete angular fluxes, and $\mathbf{\Phi}_{p,g}^{\circ}(\mathbf{r})$ the vector of $N_{m}$ flux moments $\Phi_{p,g}^{\ell,m}(\mathbf{r})$ for $0 \le \ell \le L$, $-\ell \le m \le \ell$. The transformation between the two representations is performed by two rectangular matrices:

- the **discrete-to-moment matrix** $\mathbf{D}\in\mathbb{R}^{N_{m}\times N_{d}}$, with $D_{(\ell,m),n} = w_{n}R_{\ell}^{m}(\mathbf{\Omega}_{n})$, such that
$$\mathbf{\Phi}_{p,g}^{\circ}(\mathbf{r}) = \mathbf{D}\,\mathbf{\Phi}_{p,g}^{\bullet}(\mathbf{r}) \,,$$
- the **moment-to-discrete matrix** $\mathbf{M}\in\mathbb{R}^{N_{d}\times N_{m}}$, with $M_{n,(\ell,m)} = R_{\ell}^{m}(\mathbf{\Omega}_{n})$, such that
$$\Phi_{p,g}(\mathbf{r},\mathbf{\Omega}_{n}) \approx \sum_{\ell,m} \frac{2\ell+1}{4\pi} R_{\ell}^{m}(\mathbf{\Omega}_{n}) \Phi_{p,g}^{\ell,m}(\mathbf{r}) \,.$$

### Standard (SN) Method

In the standard discrete-ordinates method, the moments are computed from the discrete fluxes using $\mathbf{D}$, and the discrete representation of the scattering source is reconstructed using $\mathbf{M}$:

$$\mathbf{Q}_{p,g}^{\texttt{B},\bullet}(\mathbf{r}) = \mathbf{M}\,\mathbf{S}\,\mathbf{D}\,\mathbf{\Phi}_{p,g}^{\bullet}(\mathbf{r}) \,,$$

where $\mathbf{S}$ is the diagonal matrix of the Legendre cross-section moments $\Sigma_{s,\ell,g'\rightarrow g}^{p'\rightarrow p}$. This method does not generally preserve the moments of the angular flux exactly; the product $\mathbf{D}\,\mathbf{M}$ is, in general, not the identity matrix when $N_{m} < N_{d}$.

### Galerkin Method

The Galerkin method enforces the consistency relation $\mathbf{D}\,\mathbf{M} = \mathbf{I}$ between the moment-to-discrete and discrete-to-moment matrices, by choosing the truncation order such that $N_{m} = N_{d}$. There are two equivalent formulations:

- **Galerkin-D**: invert the discrete-to-moment matrix, $\mathbf{M} = \mathbf{D}^{-1}$ (default in Radiant).
- **Galerkin-M**: invert the moment-to-discrete matrix, $\mathbf{D} = \mathbf{M}^{-1}$.

The Galerkin method preserves the Legendre moments of the angular flux exactly up to the truncation order, which is critical for charged-particle transport where high-order moments carry physical information.

The user selects the strategy with `sn.set_angular_boltzmann("standard"|"galerkin-d"|"galerkin-m")`.

## 1.3.4 Discretization of the Angular Fokker–Planck Operator

The angular Fokker–Planck (AFP) operator

$$Q^{\texttt{AFP}}_{p}(\mathbf{r},\mathbf{\Omega},E) = T^{p}(\mathbf{r},E) \left[\frac{\partial}{\partial \mu}\left(1-\mu^2\right)\frac{\partial}{\partial \mu} + \frac{1}{1-\mu^2}\frac{\partial^2}{\partial\phi^2}\right]\Phi_{p}(\mathbf{r},\mathbf{\Omega},E)$$

is an elliptic operator in the angular variable. Radiant provides three discretizations.

### Finite-Difference Discretization

In 1D, where only the principal cosine $\mu$ matters, the AFP operator reduces to

$$Q^{\texttt{AFP}}_{p,g}(\mathbf{r},\mu) = T^{p}_{g}(\mathbf{r}) \frac{\partial}{\partial \mu}\left[ (1-\mu^2) \frac{\partial \Phi_{p,g}(\mathbf{r},\mu)}{\partial \mu} \right]$$

and is discretized at the quadrature nodes using the conservative finite-difference scheme of Morel [morel1985improved](@cite), in which the operator is rewritten under the form

$$Q^{\texttt{AFP}}_{p,g,n}(\mathbf{r}) \approx \frac{T^{p}_{g}(\mathbf{r})}{w_{n}}\left[\beta_{n+1/2}\frac{\Phi_{p,g,n+1}(\mathbf{r}) - \Phi_{p,g,n}(\mathbf{r})}{\mu_{n+1}-\mu_{n}} - \beta_{n-1/2}\frac{\Phi_{p,g,n}(\mathbf{r}) - \Phi_{p,g,n-1}(\mathbf{r})}{\mu_{n}-\mu_{n-1}}\right] \,,$$

with the coefficients $\beta_{n\pm 1/2}$ chosen recursively so that the scheme preserves both particle conservation and the first Legendre moments of an isotropic flux. The recursion is

$$\beta_{n+1/2} = \beta_{n-1/2} - 2\mu_{n}w_{n}\,,\quad \beta_{1/2} = \beta_{N_{d}+1/2} = 0 \,.$$

The same scheme is applied in 2D and 3D using a Carlson-style level-symmetric variant.

### Differential Quadrature

The differential-quadrature scheme (only available in 1D) approximates $\partial/\partial\mu$ as a dense matrix acting on the discrete fluxes:

$$\left[\frac{\partial \Phi}{\partial\mu}\right]_{n} \approx \sum_{n'} \mathbf{A}_{n,n'}\, \Phi_{n'} \,.$$

The matrix $\mathbf{A}$ is built so that the derivatives of the orthonormal Legendre polynomials are reproduced exactly at the quadrature nodes. The full AFP operator is then represented as $\mathbf{A}\,\mathbf{B}\,\mathbf{A}$ with $\mathbf{B} = \mathrm{diag}(1-\mu_{n}^{2})$. This scheme is highly accurate for smooth fluxes but couples every direction to every other direction.

### Galerkin (Moment-Preserving) Discretization

Because the spherical harmonics $R_{\ell}^{m}(\mathbf{\Omega})$ are eigenfunctions of the Laplace–Beltrami operator,

$$\left[\frac{\partial}{\partial \mu}\left(1-\mu^2\right)\frac{\partial}{\partial \mu} + \frac{1}{1-\mu^2}\frac{\partial^2}{\partial\phi^2}\right]R_{\ell}^{m}(\mathbf{\Omega}) = -\ell(\ell+1)R_{\ell}^{m}(\mathbf{\Omega}) \,,$$

the AFP operator is diagonal in the moment representation. Using the discrete-to-moment matrix $\mathbf{D}$ and moment-to-discrete matrix $\mathbf{M}$ to switch between the two representations, the discrete AFP operator becomes

$$\mathbf{Q}^{\texttt{AFP},\bullet}_{p,g}(\mathbf{r}) = -T_{g}^{p}(\mathbf{r})\,\mathbf{M}\,\mathbf{\Lambda}\,\mathbf{D}\,\mathbf{\Phi}_{p,g}^{\bullet}(\mathbf{r}) \,,$$

with $\mathbf{\Lambda} = \mathrm{diag}\!\left[\ell(\ell+1)\right]$. This scheme preserves the Legendre moments of the angular flux exactly up to the truncation order and is the default in `DPN` and `GN` solvers.

The user selects the AFP discretization with
`sn.set_angular_fokker_planck("finite-difference"|"differential-quadrature"|"galerkin")`.

## 1.3.5 Decomposition of Scattering into Boltzmann and Fokker–Planck Parts

For interactions configured with the `"BFP"` scattering model, the differential scattering cross-section is decomposed into a "soft" component, treated by the angular Fokker–Planck operator, and a "catastrophic" component, treated by the Boltzmann operator. The decomposition is parameterized by an interaction-specific cutoff that is computed from the Legendre moments and constrained to keep the cross-sections positive at every order. The decomposition is described per interaction in Sections 2.1 and 2.2.
