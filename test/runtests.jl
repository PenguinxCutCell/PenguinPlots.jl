using Test
using CairoMakie
using PenguinPlots

PenguinPlots.activate_penguin_backend!(; backend=:cairo)

include("test_scalar_plots.jl")
include("test_vector_plots.jl")
include("test_history_plots.jl")
include("test_animation.jl")
