# 2.3 Atomic Relaxation

The ionization of an atomic electron in shell $k$ by an incoming particle (photoelectric absorption, Møller or Bhabha scattering) leaves behind a vacancy in the atom's electronic structure of element $i$. The inner-shell vacancy is filled by an outer-shell electron, leading to the emission of either a fluorescence photon (radiative transition) or an Auger / Coster–Kronig electron (non-radiative transition). The transition leaves additional vacancies in outer shells, which themselves relax, leading to a *cascade* of secondary emissions. A simplified visual diagram of these cascades can be found in Lorence et al. [lorence1989physics](@cite) and Naceur et al. [naceur2024extending](@cite).

## 2.3.1 Cascade Selection

For high-Z atoms a complete relaxation cascade may comprise hundreds of distinct paths. Two filters are applied to retain only the cascades that contribute significantly to the transport calculation:

1. Only cascades with a probability of occurrence following the initial ionization event greater than `set_minimum_probability(ηmin)` (default $\eta_{\min} = 0.1\%$) are included.
2. Among these, only the $N_{t}$ transitions which produce an electron or photon with energy above the cutoff energy $E_{G+1/2}$ are kept. Cascades that deposit all of their energy below the cutoff contribute only to the local energy / charge deposition and do not enter the multigroup scattering matrix.

## 2.3.2 Differential Cross-Section per Cascade

The differential cross-section corresponding to the production of either a fluorescence photon ($p'=\gamma$) or an Auger electron ($p'=\texttt{e-}$) along a specific cascade transition $j$ ($1\le j \le N_{t}$), associated with a produced-particle energy $\Delta E_{i,k,j}^{p'}$ and an occurrence probability $\eta_{i,k,j}^{p'}$, is

$$\sigma_{s}^{i,k,j,p\rightarrow p'}(E \rightarrow E') = \eta_{i,k,j}^{p'}\,\delta(E'-\Delta E_{i,k,j}^{p'})\,\sigma^{i,k,p}_{t}(E) \,,$$

where $\sigma^{i,k,p}_{t}(E)$ is the $k$-shell ionization cross-section associated with one of the following incident channels:

- $p = \texttt{e-}$ : inelastic-collision (Møller) cross-section on subshell $k$ (Section 2.2.1),
- $p = \texttt{e+}$ : inelastic-collision (Bhabha) cross-section on subshell $k$ (Section 2.2.1),
- $p = \gamma$ : photoelectric cross-section on subshell $k$ (Section 2.1.3).

The total cross-section for the production of relaxation radiation is therefore the sum over every cascade $j$ that survives the filtering of Section 2.3.1.

## 2.3.3 Data Source and Element Coverage

The values of $\Delta E_{i,k,j}^{p'}$ and $\eta_{i,k,j}^{p'}$ are precomputed for every element and every initial vacancy using the relaxation data of the JENDL-5 library [iwamoto2023japanese](@cite), itself derived from the EADL library [perkins1991tables](@cite), following the approach of Hébert and Naceur [hebert2023implementation](@cite). The data is available for $Z\in\{6,\ldots,100\}$ and covers the subshells K, L1–L3, M1–M5, N1–N7, O1–O7, P1–P3 and Q1. For $Z\le 5$ the produced relaxation radiation has very low energy and is ignored: the corresponding energy is treated as locally deposited.

## 2.3.4 Angular Distribution

The production of fluorescence photons and Auger electrons is assumed to be isotropic in the laboratory frame. This corresponds to using only the $\ell = 0$ Legendre moment of the differential scattering cross-section: every higher Legendre moment of the relaxation contribution vanishes. The assumption is justified by the random orientation of the residual ion at the moment of relaxation.

## 2.3.5 Coupling with the Triggering Interactions

The relaxation cascade is triggered only when the underlying ionizing interaction is computed in a subshell-dependent way. Concretely:

- `Inelastic_Collision` with `set_is_subshells_dependant(true)` triggers Møller/Bhabha-driven relaxation,
- `Photoelectric` with `model = "epdl97"` triggers photoelectric-driven relaxation.

When `Relaxation()` is added to the interaction list, every consistent combination of triggering interaction and produced particle (`Fluorescence` or `Auger`) is activated. The convenience constructors `Fluorescence()` and `Auger()` restrict the activation to the radiative or non-radiative channels, respectively.
