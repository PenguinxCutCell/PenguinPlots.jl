# PenguinPlots.jl

`PenguinPlots.jl` is a lightweight Makie-based visualization frontend for the PenguinxCutCell ecosystem.

It provides:

- solver-agnostic plotting on normalized snapshot containers
- optional solver extensions via Julia package extensions (`ext/`)
- CairoMakie-first workflows for testing, docs, and batch rendering

## Backend activation

```julia
using PenguinPlots
activate_penguin_backend!(; backend=:cairo)
```

Optional backends:

- `backend=:gl` (requires `GLMakie` in your environment)
- `backend=:wgl` (requires `WGLMakie` in your environment)

## Returning Makie objects

All plotting APIs return Makie `Figure` objects, and can draw into an existing figure/axis:

```julia
fig = Figure(resolution=(900, 500))
ax = Axis(fig[1, 1])
plot_scalar_field(snapshot; fig=fig, ax=ax)
```

## Saving figures and animations

```julia
save_figure("plot.png", fig)
animate_solution(history, "anim.gif"; fps=20)
```

## v0.1 limits

- Cartesian heatmap view of cut-cell solutions (no exact cut-cell polygon rendering yet)
- limited 3D plotting support
- basic FSI trajectory plotting from user-provided time series
