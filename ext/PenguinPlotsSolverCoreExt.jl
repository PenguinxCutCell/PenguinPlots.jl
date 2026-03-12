module PenguinPlotsSolverCoreExt

using PenguinPlots
using PenguinSolverCore

import PenguinPlots: plot_residual_history

function _residual_snapshot(history::PenguinSolverCore.CouplingHistory; title::AbstractString="Coupling residual history")
    n = length(history.iterations)
    blocks = Dict{Symbol,Vector{Float64}}()

    for rec in history.block_residuals
        for key in keys(rec)
            if !haskey(blocks, key)
                blocks[key] = fill(NaN, n)
            end
        end
    end

    for i in 1:min(n, length(history.block_residuals))
        rec = history.block_residuals[i]
        for (key, val) in rec
            blocks[key][i] = Float64(val)
        end
    end

    block_vectors = Dict{Symbol,AbstractVector}(k => v for (k, v) in blocks)
    return PenguinPlots.ResidualHistorySnapshot(
        history.iterations,
        history.residuals;
        block_residuals=block_vectors,
        title=title,
        metadata=Dict{Symbol,Any}(:source => :coupling_history),
    )
end

function plot_residual_history(
    history::PenguinSolverCore.CouplingHistory;
    logy::Bool=true,
    show_blocks::Bool=true,
    title::AbstractString="Coupling residual history",
    kwargs...,
)
    snap = _residual_snapshot(history; title=title)
    return PenguinPlots.plot_residual_history(snap; logy=logy, show_blocks=show_blocks, kwargs...)
end

plot_residual_history(result::Tuple{Any,PenguinSolverCore.CouplingHistory}; kwargs...) =
    plot_residual_history(result[2]; kwargs...)

plot_residual_history(result::NamedTuple{(:problem, :history),<:Tuple{Any,PenguinSolverCore.CouplingHistory}}; kwargs...) =
    plot_residual_history(result.history; kwargs...)

end # module
