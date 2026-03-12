# API

## Core functions

- `activate_penguin_backend!(; backend=:cairo, kwargs...)`
- `plot_scalar_field(args...; kwargs...)`
- `plot_interface_field(args...; kwargs...)`
- `plot_vector_field(args...; kwargs...)`
- `plot_pressure_field(args...; kwargs...)`
- `plot_solution(args...; kwargs...)`
- `plot_bulk(args...; kwargs...)`
- `plot_interface(args...; kwargs...)`
- `plot_history(args...; kwargs...)`
- `plot_residual_history(args...; kwargs...)`
- `animate_solution(args...; kwargs...)`
- `save_figure(path, fig; kwargs...)`

## Snapshot containers

- `ScalarSnapshot`
- `InterfaceSnapshot`
- `VectorSnapshot`
- `ResidualHistorySnapshot`
