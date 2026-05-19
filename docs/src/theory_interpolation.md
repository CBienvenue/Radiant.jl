# 3.3 Interpolations

Radiant relies on tabulated atomic data (cross-sections, form factors, scattering functions, mean ionization energies, etc.) to build its multigroup libraries. Interpolation between tabulated points is used throughout the cross-section preparation step.

## 3.3.1 Linear and Bilinear Interpolation

Given two abscissae $(x_{i-1},x_{i})$ enclosing a target $x_{i} \le \tilde{x} \le x_{i+1}$ and the associated function values $(f_{i-1},f_{i})$, the linear interpolant is

$$\tilde{f}(\tilde{x}) = \frac{(x_{i}-\tilde{x})\,f_{i-1} + (\tilde{x}-x_{i-1})\,f_{i}}{x_{i}-x_{i-1}} \,.$$

Its bilinear extension on a tensor grid $(x_{i},y_{j})$ is

$$\tilde{f}(\tilde{x},\tilde{y}) = \frac{\sum_{p,q\in\{0,1\}} f_{i-1+p,j-1+q}\,\left(x_{i-p}-\tilde{x}\right)^{1-p}\left(\tilde{x}-x_{i-1+p}\right)^{p}\left(y_{j-q}-\tilde{y}\right)^{1-q}\left(\tilde{y}-y_{j-1+q}\right)^{q}}{(x_{i}-x_{i-1})(y_{j}-y_{j-1})} \,.$$

Linear interpolation is fast and bounded by the surrounding values; it is the default for monotonic and slowly-varying data.

## 3.3.2 Log–Log Interpolation

For cross-section data exhibiting power-law behaviour over several decades — typically photoelectric or pair-production cross-sections — interpolation is performed in $(\log x,\log f)$ space:

$$\log \tilde{f}(\tilde{x}) = \frac{(\log x_{i}-\log\tilde{x})\,\log f_{i-1} + (\log\tilde{x}-\log x_{i-1})\,\log f_{i}}{\log x_{i}-\log x_{i-1}} \,.$$

This is equivalent to linear interpolation of the exponent of a local power law $f(x) = a\,x^{b}$.

## 3.3.3 Natural Cubic Spline

For data where smoothness of the first and second derivative is required (for example mean ionization energies or scattering functions used in differential cross-section integrals), Radiant uses a natural cubic spline. The spline is defined piecewise on each subinterval $[x_{i},x_{i+1}]$ as a cubic polynomial

$$s_{i}(x) = a_{i}\,x^{3} + b_{i}\,x^{2} + c_{i}\,x + d_{i} \,, \quad x \in [x_{i},x_{i+1}] \,,$$

where the $4n$ unknowns (for $n$ intervals) are determined by:

- $2n$ interpolation conditions $s_{i}(x_{i}) = y_{i}$ and $s_{i}(x_{i+1}) = y_{i+1}$ for $i = 1,\ldots,n$,
- $(n-1)$ continuity conditions on the first derivative at the internal nodes $s'_{i-1}(x_{i}) = s'_{i}(x_{i})$,
- $(n-1)$ continuity conditions on the second derivative $s''_{i-1}(x_{i}) = s''_{i}(x_{i})$,
- 2 "natural" end conditions $s''_{1}(x_{1}) = s''_{n}(x_{n+1}) = 0$.

The resulting block-banded system is solved once and reused for every later evaluation of the spline. Splines are used by Radiant when interpolated derivatives are needed (e.g. in tabulated form factors that enter the differential cross-sections).

## 3.3.4 Choice of Interpolant in Radiant

The interpolant is selected at the data-loading stage based on the structure of the table:

| Tabulated quantity                        | Interpolation                |
|-------------------------------------------|------------------------------|
| Subshell binding energies                 | Linear in $E$                |
| Mean excitation energies                  | Linear in $E$                |
| Cross-sections vs. $E$ (over decades)     | Log–log in $(E,\sigma)$      |
| Form factors / incoherent scattering func.| Natural cubic spline         |
| Stopping powers vs. $E$                   | Log–log                      |

These choices reproduce, in deterministic form, the conventions used in the EPDL/EEDL/EADL evaluations and in Monte-Carlo codes such as EGSnrc and Penelope.
