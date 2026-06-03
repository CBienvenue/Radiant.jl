"""
    DPN()

Build a Double-PN (DPN) solver. The Double-PN method is the generalized-harmonics
(GN) method at the coarsest angular subdivision: one angular patch per base domain
(`Nv = 1`, `L_elem = L`), i.e. two half-spheres in 1D, four quadrants in 2D and eight
octants in 3D. `DPN` is therefore a thin constructor that returns a [`GN`](@ref) solver
pre-configured with these base parameters:

- `subdivision = 1` (one angular patch per base domain),
- `z_fold = true` (the 2D angular domain is folded into four quadrants),
- `tiling = "polar-anchored"` (the DPN-equivalent octant tiling in 3D),
- `legendre_order_local = legendre_order` (the local patch order tracks the global
  Legendre order; the single-argument `set_legendre_order` keeps this coupling).

The returned object is a `GN`; configure it with the usual `GN` setters.

# Examples
```jldoctest
julia> m = DPN()
julia> m.set_solver_type("BTE")
julia> m.set_polynomial_basis("legendre")
julia> m.set_legendre_order(7)
```
"""
function DPN()
    gn = GN()
    gn.legendre_order       = 8
    gn.legendre_order_local = 8                 # L_elem = L (full coupling, one patch)
    gn.subdivision          = 1                 # Nv = 1
    gn.tiling               = "polar-anchored"  # DPN-equivalent octant tiling in 3D
    gn.z_fold               = true              # 2D → four quadrants
    return gn
end
