using CairoMakie
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDiffusion
using PenguinPlots

activate_penguin_backend!(; backend=:cairo)

grid = (range(0.0, 1.0; length=41), range(0.0, 1.0; length=31))
moms = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)

bc = BorderConditions(
    ;
    left=Dirichlet(0.0),
    right=Dirichlet(0.0),
    bottom=Dirichlet(0.0),
    top=Dirichlet(0.0),
)

ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))
source(x, y, t) = 2pi^2 * sin(pi * x) * sin(pi * y)
model = DiffusionModelMono(cap, ops, 1.0; source=source, bc_border=bc)
sys = solve_steady!(model)

fig = plot_solution(model, sys.x; field=:omega, show_interfacial=true)
save_figure(joinpath(@__DIR__, "diffusion_mono_steady.png"), fig)

fig
