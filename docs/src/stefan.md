# Stefan Plots and Animation

With `PenguinStefan.jl` loaded, `plot_solution` and `animate_solution` support Stefan solvers/results.

```julia
fig_u = plot_solution(solver; field=:uω)
fig_s = plot_solution(solver; field=:speed)
fig_h = plot_history(result; field=:speedmax)
```

Animation from `solve!` history:

```julia
animate_solution(result, "stefan.gif"; field=:uω, show_phi=true, fps=12)
```

Supported fields include:

- mono: `:uω`, `:uγ`, `:speed`, `:phi`
- diphasic: `:uω1`, `:uω2`, `:uγ1`, `:uγ2`, `:speed`, `:phi`, `:both`
