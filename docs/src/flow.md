# Flow Plots

Use `VectorSnapshot` to visualize pressure/velocity fields.

```julia
snap = VectorSnapshot(x, y, u, v; p=p, mask=mask)
fig = plot_vector_field(snap; stride=3, velocity_scale=1.0)
```

## Pressure-only

```julia
fig = plot_pressure_field(snap)
```

## Stokes extension entrypoints

When `PenguinStokes.jl` is loaded:

```julia
fig = plot_solution(model, x; field=:pressure)
fig = plot_solution(model, x; field=:velocity)
fig = plot_solution(model, x; field=:pressure_velocity)
```

Two-phase models support `phase=:both` side-by-side plots.
