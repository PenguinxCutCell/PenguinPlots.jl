module PenguinPlots

using CairoMakie
using GeometryBasics
using LinearAlgebra
using Makie
using Statistics

function plot_solution end
function plot_bulk end
function plot_interface end
function plot_history end
function plot_residual_history end
function animate_solution end

include("core/types.jl")
include("core/utils.jl")
include("core/backend.jl")
include("core/grids.jl")
include("core/scalar.jl")
include("core/vector.jl")
include("core/history.jl")
include("core/animation.jl")

plot_solution(args...; kwargs...) = throw(ArgumentError(
    "No `plot_solution` method for arguments $(map(typeof, args)). Load the corresponding solver package so its PenguinPlots extension is activated."
))

plot_bulk(args...; kwargs...) = throw(ArgumentError(
    "No `plot_bulk` method for arguments $(map(typeof, args)). Load PenguinDiffusion to enable diffusion-specific wrappers."
))

plot_interface(args...; kwargs...) = throw(ArgumentError(
    "No `plot_interface` method for arguments $(map(typeof, args)). Load PenguinDiffusion to enable diffusion-specific wrappers."
))

plot_history(args...; kwargs...) = throw(ArgumentError(
    "No `plot_history` method for arguments $(map(typeof, args)). Use `plot_history(times, values; ...)` or load a solver extension."
))

plot_residual_history(args...; kwargs...) = throw(ArgumentError(
    "No `plot_residual_history` method for arguments $(map(typeof, args)). Use `ResidualHistorySnapshot` or load PenguinSolverCore."
))

animate_solution(args...; kwargs...) = throw(ArgumentError(
    "No `animate_solution` method for arguments $(map(typeof, args)). Use snapshot histories or load a solver extension."
))

export ScalarSnapshot
export InterfaceSnapshot
export VectorSnapshot
export ResidualHistorySnapshot

export activate_penguin_backend!
export plot_scalar_field
export plot_interface_field
export plot_vector_field
export plot_pressure_field
export plot_solution
export plot_bulk
export plot_interface
export plot_history
export plot_residual_history
export animate_solution
export save_figure

end # module
