using CairoMakie
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDiffusion
using PenguinPlots

activate_penguin_backend!(; backend=:cairo)

grid = (range(0.0, 1.0; length=41), range(0.0, 1.0; length=31))
moms = geometric_moments((x, y) -> x - 0.5, grid, Float64, nan; method=:vofijul)
cap1 = assembled_capacity(moms; bc=0.0)
cap2 = assembled_capacity(geometric_moments((x, y) -> -(x - 0.5), grid, Float64, nan; method=:vofijul); bc=0.0)

bc = BorderConditions(
    ;
    left=Dirichlet(1.0),
    right=Dirichlet(0.0),
    bottom=Dirichlet(0.0),
    top=Dirichlet(0.0),
)

ops1 = DiffusionOps(cap1; periodic=periodic_flags(bc, 2))
ops2 = DiffusionOps(cap2; periodic=periodic_flags(bc, 2))
model = DiffusionModelDiph(cap1, ops1, 1.0, 0.0, cap2, ops2, 0.6, 0.0; bc_border=bc)
sys = solve_steady!(model)

fig = plot_solution(model, sys.x; field=:omega, phase=:both, show_interfacial=true)
save_figure(joinpath(@__DIR__, "diffusion_diph_fixed_interface.png"), fig)

fig
