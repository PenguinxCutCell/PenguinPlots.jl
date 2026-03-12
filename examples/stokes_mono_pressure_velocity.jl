using CairoMakie
using CartesianGrids: CartesianGrid
using PenguinBCs
using PenguinPlots
using PenguinStokes

activate_penguin_backend!(; backend=:cairo)

full_body(args...) = -1.0
grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (25, 21))

bc = BorderConditions(
    ;
    left=Dirichlet(0.0),
    right=Dirichlet(0.0),
    bottom=Dirichlet(0.0),
    top=Dirichlet(0.0),
)

model = StokesModelMono(
    grid,
    full_body,
    1.0,
    1.0;
    bc_u=(bc, bc),
    bc_cut=Dirichlet(0.0),
    force=(0.0, 0.0),
)
sys = solve_steady!(model)

fig = plot_solution(model, sys.x; field=:pressure_velocity, stride=3)
save_figure(joinpath(@__DIR__, "stokes_mono_pressure_velocity.png"), fig)

fig
