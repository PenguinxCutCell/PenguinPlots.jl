# Scalar Plots

Use `ScalarSnapshot` for generic scalar fields.

```julia
snap = ScalarSnapshot(x, y, values, mask;
    title="Temperature",
    colorlabel="u",
)
fig = plot_scalar_field(snap; colormap=:viridis)
```

## Interface overlays

- pass `interface_xy=(xγ, yγ)` inside the snapshot
- or pass `phi` to draw the `phi=0` contour
- pass `uγ` to color interface points with interfacial unknowns

```julia
fig = plot_scalar_field(snap; show_interface=true, phi=phi, uγ=uγ)
```

## Diphasic side-by-side panels

```julia
fig = plot_scalar_field([snap_phase1, snap_phase2]; shared_clims=true)
```

## Interface unknown plotting

```julia
isnap = InterfaceSnapshot(xγ, yγ, uγ; label="uγ")
fig = plot_interface_field(isnap; style=:scatter)
```
