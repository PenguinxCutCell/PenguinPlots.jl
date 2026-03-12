using CairoMakie
using CartesianGrids: CartesianGrid
using PenguinBCs
using PenguinPlots
using PenguinStefan

activate_penguin_backend!(; backend=:cairo)

grid = CartesianGrid((0.0,), (1.0,), (81,))
rep = LevelSetRep(grid, x -> x - 0.5)

bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
params = StefanParams(0.5, 2.0; kappa1=1.0, source1=0.0)
opts = StefanOptions(; scheme=:BE, reinit=false, extend_iters=8, interface_band=2.0)
prob = StefanMonoProblem(grid, bc, params, rep, opts)

solver = build_solver(
    prob;
    t0=0.0,
    uω0=(x, t) -> 1.0 - 0.8 * x,
    uγ0=0.5,
    speed0=0.0,
)

result = solve!(solver, (0.0, 0.06); dt=0.01, save_history=true)

fig_u = plot_solution(solver; field=:uω)
save_figure(joinpath(@__DIR__, "stefan_mono_snapshot.png"), fig_u)

fig_h = plot_history(result; field=:speedmax)
save_figure(joinpath(@__DIR__, "stefan_mono_speedmax.png"), fig_h)

animate_solution(result, joinpath(@__DIR__, "stefan_mono_animation.gif"); field=:uω, fps=10)

fig_u
