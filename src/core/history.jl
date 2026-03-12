function plot_history(
    times::AbstractVector,
    values::AbstractVector;
    ylabel::AbstractString="value",
    xlabel::AbstractString="t",
    logy::Bool=false,
    title=nothing,
    fig=nothing,
    ax=nothing,
    label=nothing,
    linewidth::Real=2,
    color=:steelblue,
    kwargs...,
)
    length(times) == length(values) || throw(DimensionMismatch("times and values must have same length"))
    fig_local, ax_local, _ = _ensure_fig_ax(; fig=fig, ax=ax)

    lkwargs = (; linewidth=linewidth, color=color, kwargs...)
    if label === nothing
        lines!(ax_local, times, values; lkwargs...)
    else
        lines!(ax_local, times, values; lkwargs..., label=String(label))
        axislegend(ax_local; position=:rt)
    end

    logy && (ax_local.yscale = Makie.log10)
    ax_local.xlabel = xlabel
    ax_local.ylabel = ylabel
    if title !== nothing
        ax_local.title = String(title)
    end
    return fig_local
end

function plot_residual_history(
    snapshot::ResidualHistorySnapshot;
    logy::Bool=true,
    show_blocks::Bool=true,
    fig=nothing,
    ax=nothing,
    title=nothing,
    linewidth::Real=2,
    kwargs...,
)
    fig_local, ax_local, _ = _ensure_fig_ax(; fig=fig, ax=ax)

    lines!(
        ax_local,
        snapshot.iterations,
        snapshot.residuals;
        label="global",
        linewidth=linewidth,
        color=:black,
        kwargs...,
    )

    if show_blocks
        for name in sort(collect(keys(snapshot.block_residuals)); by=String)
            vals = snapshot.block_residuals[name]
            n = min(length(snapshot.iterations), length(vals))
            n == 0 && continue
            lines!(
                ax_local,
                snapshot.iterations[1:n],
                vals[1:n];
                label=String(name),
                linewidth=max(linewidth - 0.3, 1.0),
            )
        end
    end

    logy && (ax_local.yscale = Makie.log10)
    ax_local.xlabel = "iteration"
    ax_local.ylabel = "residual"
    ax_local.title = title === nothing ? snapshot.title : String(title)

    if show_blocks || !isempty(snapshot.residuals)
        axislegend(ax_local; position=:rt)
    end

    return fig_local
end
