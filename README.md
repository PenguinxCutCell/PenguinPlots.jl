# PenguinPlots.jl

Makie-based plotting frontends for the PenguinxCutCell ecosystem.

`PenguinPlots.jl` keeps a lightweight solver-agnostic core and loads solver-specific methods through Julia package extensions.

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/PenguinxCutCell/PenguinPlots.jl")
```

For local ecosystem development, keep sibling packages in the same workspace and activate this package project.

## Weak deps and extensions

Solver-specific plotting methods are loaded only when the package is present:

- `PenguinDiffusion` -> `PenguinPlotsDiffusionExt`
- `PenguinStokes` -> `PenguinPlotsStokesExt`
- `PenguinStefan` -> `PenguinPlotsStefanExt`
- `PenguinSolverCore` -> `PenguinPlotsSolverCoreExt`

## Quick start

```julia
using PenguinPlots, CairoMakie
activate_penguin_backend!(; backend=:cairo)
```

### Diffusion

```julia
fig = plot_solution(diffusion_model, state; field=:omega)
save_figure("diffusion.png", fig)
```

### Stokes

```julia
fig = plot_solution(stokes_model, state; field=:pressure_velocity, stride=3)
```

### Stefan animation

```julia
result = solve!(solver, (0.0, 0.1); dt=0.01)
animate_solution(result, "stefan.gif"; field=:uω, show_phi=true)
```

### Coupling residuals

```julia
fig = plot_residual_history(history; logy=true, show_blocks=true)
```

## Example outputs

![Diffusion](examples/diffusion_mono_steady.png)
![Stokes](examples/stokes_mono_pressure_velocity.png)
![Coupling](examples/coupling_residuals_twoway.png)

## v0.1 limitations

- Cartesian heatmap rendering (no exact cut-cell polygon meshing yet)
- limited 3D plotting support
- simple FSI history plotting only when time series are available
