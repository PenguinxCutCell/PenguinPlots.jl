using CairoMakie
using PenguinPlots
using PenguinSolverCore

activate_penguin_backend!(; backend=:cairo)

struct DriverModel
    value::Float64
end

struct FollowerModel
    α::Float64
end

mutable struct DriverState
    velocity::Float64
end

mutable struct FollowerState
    velocity::Float64
    response::Float64
end

PenguinSolverCore.initialize_state(::DriverModel, init) = DriverState(init === nothing ? 0.0 : Float64(init))
PenguinSolverCore.initialize_state(::FollowerModel, init) = FollowerState(init === nothing ? 0.0 : Float64(init), 0.0)

PenguinSolverCore.copy_state(s::DriverState) = DriverState(s.velocity)
PenguinSolverCore.copy_state(s::FollowerState) = FollowerState(s.velocity, s.response)

function PenguinSolverCore.advance_steady!(block::CoupledBlock{<:DriverModel,<:DriverState}; kwargs...)
    block.state.velocity = block.model.value
    return block
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{<:FollowerModel,<:FollowerState}; kwargs...)
    block.state.response = block.model.α * block.state.velocity
    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{<:DriverModel,<:DriverState}, ::Val{:velocity})
    return block.state.velocity
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{<:FollowerModel,<:FollowerState}, ::Val{:velocity})
    return block.state.velocity
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{<:FollowerModel,<:FollowerState}, ::Val{:response})
    return block.state.response
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{<:DriverModel,<:DriverState}, ::Val{:velocity}, data)
    block.state.velocity = Float64(data)
    return block
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{<:FollowerModel,<:FollowerState}, ::Val{:velocity}, data)
    block.state.velocity = Float64(data)
    return block
end

b1 = CoupledBlock(:driver, DriverModel(1.0), DriverState(0.0))
b2 = CoupledBlock(:follower, FollowerModel(0.8), FollowerState(0.0, 0.0))

maps = [CouplingMap(:driver, :follower, :velocity)]

oneway = CoupledProblem([b1, b2], OneWayCoupling(:driver, :follower); maps=maps)
_, hist_oneway = solve_coupled!(oneway; return_history=true)
fig1 = plot_residual_history(hist_oneway; title="One-way coupling residual")
save_figure(joinpath(@__DIR__, "coupling_residuals_oneway.png"), fig1)

b3 = CoupledBlock(:driver, DriverModel(1.0), DriverState(0.0))
b4 = CoupledBlock(:follower, FollowerModel(0.8), FollowerState(0.0, 0.0))

maps2 = [
    CouplingMap(:driver, :follower, :velocity),
    CouplingMap(:follower, :driver, :velocity; apply! = (data, from, to, problem) -> 0.5 * data),
]

twoway = CoupledProblem([b3, b4], TwoWayCoupling([:driver, :follower]; maxiter=10, atol=1e-10, rtol=1e-8); maps=maps2)
_, hist_twoway = solve_coupled!(twoway; return_history=true, throw_on_nonconvergence=false)
fig2 = plot_residual_history(hist_twoway; title="Two-way coupling residual")
save_figure(joinpath(@__DIR__, "coupling_residuals_twoway.png"), fig2)

fig2
