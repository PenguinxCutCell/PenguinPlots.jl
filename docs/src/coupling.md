# Coupling Residuals

`PenguinPlots.jl` can visualize coupled Picard histories from `PenguinSolverCore.jl`.

```julia
fig = plot_residual_history(history; logy=true, show_blocks=true)
```

Accepted inputs:

- `ResidualHistorySnapshot`
- `PenguinSolverCore.CouplingHistory` (extension)
- `(problem, history)` tuples returned by `solve_coupled!(...; return_history=true)`

Use logarithmic y-scale by default for outer coupling residuals.
