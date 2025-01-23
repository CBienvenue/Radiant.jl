# 3.1 Polynomials

## 3.1.1 Legendre polynomials

The Legendre polynomials, noted $P_{\ell}(x)$ for $\ell \ge 0$, are a complete and orthogonal system defined over $x \in [-1,1]$ such as

$$\int_{-1}^{1}dxP_{\ell}(x)P_{\ell'}(x) = \frac{2}{2\ell+1}\delta_{\ell,\ell'}$$

where $\delta_{\ell,\ell'}$ is the Kronecker delta. They also can be defined as the eigenfunctions of the Legendre's differential equation

$$\frac{d}{dx}\left[(1-x^2)\frac{d}{dx}\right]P_{\ell}(x) = -\ell(\ell+1)P_{\ell}(x)$$

The Legendre polynomials can be calculated using the Bonnet's recursion formula with

$$P_{\ell}(x) = \begin{cases} 1 & \ell = 0 \\ x & \ell = 1 \\ \left[\dfrac{2\ell-1}{\ell}\right]xP_{\ell-1}(x) - \left[\dfrac{\ell-1}{\ell}\right]P_{\ell-2}(x) & \ell \ge 2 \end{cases}$$

It can be shown that the Legendre polynomials are bounded by [elbert1994inequality](@cite)

$$\left|P_{\ell}(x)\right| \le 1$$

The $\ell$-order Legendre moments of a function $f(x)$ can be defined as

$$f_{\ell} = \int_{-1}^{1} dx P_{\ell}(x)f(x)$$

and, if the function $f(x)$ is positive, triangle inequality applied to integral and using the Legendre polynomial upper bound, we can show the following inequality for $\ell \ge 1$ Legendre moments

$$\frac{\left|f_{\ell}\right|}{f_{0}} = \frac{\left| \displaystyle\int_{-1}^{1} dx P_{\ell}(x)f(x) \right|}{ f_{0}} \le \frac{\displaystyle\int_{-1}^{1} dx \left|P_{\ell}(x)\right|f(x)}{f_{0}} \le 1$$

## 3.1.2 Normalized Legendre polynomials

The normalized Legendre polynomials, noted $\tilde{P}_{\ell}(x)$ for $\ell \ge 0$, are a complete and orthogonal system defined over $x \in [-1/2,1/2]$ such as

$$\int_{-1/2}^{1/2}dx\tilde{P}_{\ell}(x)\tilde{P}_{\ell'}(x) = \delta_{\ell,\ell'}$$

They can be calculated using the following recursion formula [bienvenue2022high](@cite)

$$\tilde{P}_{\ell}(x) = \begin{cases} 1 & \ell = 0 \\ 2\sqrt{3}x & \ell = 1 \\ 2\sqrt{\dfrac{2\ell+1}{2\ell-1}}\left[\dfrac{2\ell-1}{\ell}\right]x\tilde{P}_{\ell-1}(x) - \sqrt{\dfrac{2\ell+1}{2\ell-3}}\left[\dfrac{\ell-1}{\ell}\right]\tilde{P}_{\ell-2}(x) & \ell \ge 2 \end{cases}$$

An upper bound for the size of the normalized Legendre polynomial can be derived from the one of the classical Legendre polynomial and is given by

$$\left|\tilde{P}_{\ell}(x)\right| \le \sqrt{2\ell+1}$$

The $\ell$-order normalized Legendre moments of a function $f(x)$ can be defined as

$$f_{\ell} = \int_{-1}^{1} dx \tilde{P}_{\ell}(x)f(x)$$

and, if the function $f(x)$ is positive, triangle inequality applied to integral and using the normalized Legendre polynomial upper bound, we can show the following inequality for $\ell \ge 1$ Legendre moments

$$\frac{\left|f_{\ell}\right|}{f_{0}} = \frac{\left| \displaystyle\int_{-1/2}^{1/2} dx \tilde{P}_{\ell}(x)f(x) \right|}{ f_{0}} \le \frac{\displaystyle\int_{-1/2}^{1/2} dx \left|\tilde{P}_{\ell}(x)\right|f(x)}{f_{0}} \le \sqrt{2\ell+1}$$

## 3.1.3 Associated Legendre polynomials (Ferrer definition)

The associated Legendre polynomials, using the Ferrer definition in which the factor $(-1)^m$ is absent, is given by

$$P_{\ell}^{m}(x) = (1-x^2)^{m/2}\frac{d^{m}}{dx^{m}}P_{\ell}(x), \quad \ 0 \le m \le \ell \\$$

To compute this term numerically, one can initially use the following identity (from [arfken2011mathematical](@cite), Eq.~15.93)

$$P_{m}^{m}(x) = (1-x^2)^{m/2} (2m-1)!! = (1-x^2)^{m/2} \prod_{k=1}^{m}(2k+1)$$

If $\ell = m$, then computation is over, otherwise, one should use this second identity (from [arfken2011mathematical](@cite), Eq.~15.94)

$$P^{m}_{m+1}(x) = (2m+1)xP_{m}^{m}(x)$$

If $\ell = m+1$, then computation is over, otherwise, this third identity (from [arfken2011mathematical](@cite), Eq.~15.88) iteratively can, since $m \le \ell$, be use for all other case

$$P^{m}_{\ell}(x) = \frac{(2\ell-1)xP^{m}_{\ell-1}(x)-(\ell+m-1)P^{m}_{\ell-2}(x)}{\ell-m}$$

## 3.1.4 Real spherical harmonics

The real spherical harmonics, noted $R_{\ell}^{m}(\mu,\phi)$ for $\ell \ge 0$ and $|m| < \ell$, are a complete and orthogonal system defined over $\mu \in [-1,1]$ and $\phi \in [0,2\pi]$ such as

$$\int_{4\pi}d^{2}\Omega R_{\ell}^{m}(\vec{\Omega}) R_{\ell'}^{m'}(\vec{\Omega}) = \frac{4\pi}{2\ell+1}\delta_{\ell,\ell'}\delta_{m,m'}$$

They also can be defined as the eigenfunctions of the Laplace equation given by

$$\left[\frac{\partial}{\partial \mu}\left[(1-\mu^2)\frac{\partial}{\partial \mu}\right] + \frac{1}{1-\mu^2}\frac{\partial^2}{\partial\phi^2}\right]R_{\ell}^{m} = -\ell(\ell+1)R_{\ell}^{m}$$

The real spherical harmonics are given by [hebert2009applied](@cite)

$$R_{\ell}^{m}(\vec{\Omega}) = \sqrt{(2-\delta_{m,0})\frac{(\ell-\left|m\right|)!}{(\ell+\left|m\right|)!}} P_{\ell}^{\left|m\right|}(\mu) \mathcal{T}_{m}(\phi)$$

with 

$$\mathcal{T}_{m}(\phi) = \begin{cases} \cos{m\phi}, &\quad m \ge 0 \\ \sin{\left|m\right|\phi}, &\quad m < 0 \end{cases}$$