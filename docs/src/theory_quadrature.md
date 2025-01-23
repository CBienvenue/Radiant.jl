# 3.2 Quadratures

## 3.2.1 Gauss-Legendre quadrature

The $N$-order Gauss-Legendre quadrature is defined over $x \in [-1,1]$. It gives exact result for polynomial of degree $2N-1$ or less. The $i$-th node $x_{i}$ correspond to the $i$-th root of the Legendre polynomial $P_{N}(x)$. It can be found using Newton iteration method, if a good choice of initial node is done. Such a choice is proposed by Tricomi [tricomi1950sugli,hale2013fast](@cite):

$$x_{i}^{0} = \left[1 - \frac{N-1}{8N^3} - \frac{1}{384N^4}\left(39-\frac{28}{\sin^{2}(\phi_i)}\right)\right]\cos(\phi_i)$$

with

$$\phi_i = \frac{\pi}{2}\left(\frac{4i-1}{2N+1}\right)$$

The derivative of Legendre polynomial, needed for the Newton method, is given by

$$P_{N}'(x) = \frac{N\left(P_{N-1}(x)-xP_{N}(x)\right)}{1-x^2}$$

Finally, the weights are given by

$$\omega_{i} = \frac{2}{(1-x_{i}^2)\left(P_{N}'(x_{i})\right)^2}$$

## 3.2.2 Gauss-Lobatto quadrature

The $N$-order Gauss-Lobatto quadrature is defined over $x \in [-1,1]$. It gives accurate result for polynomial of degree $2N-3$ or less. The $i$-th node $x_{i}$ correspond to the $(i-1)$-th root of the Legendre polynomial derivative $P_{N-1}'(x)$. The first and last node are respectively $x_{1} = -1$ and $x_{N} = 1$. It can be found using Newton iteration method, if a good choice of initial node is done. A good choice of initial value is

$$x_{i} = \cos\left((i-1)\frac{\pi}{N-1}\right)$$

The second derivative of Legendre polynomial, needed for the Newton method, is given by

$$P_{N}''(x) = \frac{2xP_{N}'(x)-N(N+1)P_{N}(x)}{1-x^2}$$

Finally, the weights are given by

$$\omega_{i} = \begin{cases}\dfrac{2}{n(n-1)\left(P_{N-1}(x_{i})\right)^2} & x_{i} \neq \pm 1 \\\dfrac{2}{n(n-1)} & \text{otherwise}\end{cases}$$

## 3.2.3 Product quadrature with Chebychev azimuthal quadrature

The $N$-order product quadrature with Chebychev azimuthal quadrature is defined over the unit sphere. This method is symmetric over 8 octants and has no evaluation point along the reference frame axis. The points distribution over the sphere is uneven. The quadrature has $N$ level containing $2N$ points each. The total number of directions is

$$N_{d} = 2N^{2}$$

For a level $n \in 1,..,N$ and the point $m \in 1,..,2N$ on that level, corresponding to an index $i = 2N(n-1)+m$, with $\mu_{n}$ and $w_{n}$ being the nodes and weights of the $N$-order quadrature over a line segment (e.g. Gauss-Legendre), the direction cosines are given by

$$\begin{split} \mu_{i} &= \mu_{n} \\ \eta_{i} &= \sqrt{1-\mu_{n}^2}\cos(\phi_m) \\ \xi_{i} &= \sqrt{1-\mu_{n}^2}\sin(\phi_m) \end{split}$$

where 

$$\phi_m = \frac{(2m-1)\pi}{2N}$$ 

and the corresponding weight is given by

$$w_{i} = w_{n} \frac{\pi}{N}$$

## 3.2.4 Carlson quadrature

The Carlson quadrature [carlson1976method](@cite) is defined over the unit sphere and require less evaluation point than product quadrature. This method is symmetric over 8 octants and has no evaluation point along the reference frame axis. In the octant with positive cosines, with $N$-order Carlson quadrature, the weight are given by

$$\omega_{m} = \frac{4(N-2n+2)}{N(N+2)}$$

with $m \in \left\{1,N/2\right\}$. The main director cosine is given by

$$\mu_{m} = \bar{\mu}_{m} + f\mu_{m-1/2}$$

with

$$\bar{\mu}_{m} = 1 - \frac{(N-2n+2)^2}{N(N+2)}$$

and

$$\mu_{m-1/2} = 1 - \frac{(N-2n+2)(N-2n+4)}{N(N+2)}$$

The factor $f$ is determined by finding the root of the following equation 

$$\sum_{m=1}^{N/2}\omega_{m}\mu_{m}^2 = 1/3$$

using the Newton-bisection method. The two other director cosines are given

$$\eta_{m,k} = \sqrt{1-\mu_{m}^{2}}\sin\theta_{m,k} \quad \text{and} \quad \xi_{m,k} = \sqrt{1-\mu_{m}^{2}}\cos\theta_{m,k}$$

where

$$\theta_{m,k} = \frac{\pi}{2}\left[\frac{2m-1}{N-2k+2}A_{n}+\frac{1-A_{n}}{2}\right]$$

with $k \in \left\{1,N/2-m+1\right\}$. The $A_{n}$ factors are determined by finding the root of the following equation 

$$\omega_{0}\sum_{m=1}^{N/2}\sum_{m=1}^{N/2-m+1}\eta_{m,k} = \sum_{m=1}^{N/2}\omega_{m}\mu_{m}$$

using the Newton-bisection.

## 3.2.5 Lebedev quadrature

!!! warning 
    Incomplete - Work in progress
