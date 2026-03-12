function _scalar_mode_1d(snapshot::ScalarSnapshot)
    snapshot.y === nothing && return true
    return length(snapshot.y) <= 1 || ndims(snapshot.values) == 1
end

function _plot_scalar_1d!(
    snapshot::ScalarSnapshot,
    fig,
    ax;
    show_interface::Bool=true,
    title=nothing,
    phase_label::Bool=true,
    linewidth::Real=2,
    linecolor=:steelblue,
    interface_color=:black,
    uγ=nothing,
    marker_size::Real=8,
    kwargs...,
)
    x = collect(snapshot.x)
    vals = vec(collect(snapshot.values))
    mask = vec(_coerce_mask(snapshot.mask, size(vals)))
    y = similar(vals, Float64)
    @inbounds for i in eachindex(vals)
        y[i] = mask[i] ? Float64(vals[i]) : NaN
    end

    lines!(ax, x, y; linewidth=linewidth, color=linecolor, kwargs...)
    ax.xlabel = "x"
    ax.ylabel = isempty(snapshot.colorlabel) ? "value" : snapshot.colorlabel

    if show_interface && snapshot.interface_xy !== nothing
        xγ = collect(snapshot.interface_xy[1])
        if !isempty(xγ)
            vlines!(ax, xγ; color=interface_color, linewidth=1.5, linestyle=:dash)
        end
        if uγ !== nothing
            uγv = vec(collect(uγ))
            n = min(length(xγ), length(uγv))
            n > 0 && scatter!(ax, xγ[1:n], uγv[1:n]; color=uγv[1:n], colormap=:balance, markersize=marker_size)
        end
    end

    ttl = title === nothing ? _default_title(snapshot.title, snapshot.time, snapshot.phase; phase_label=phase_label) : String(title)
    ax.title = ttl
    return (fig=fig, ax=ax, plot=nothing)
end

function _plot_scalar_2d!(
    snapshot::ScalarSnapshot,
    fig,
    ax;
    colormap=:viridis,
    show_interface::Bool=true,
    show_cutcells::Bool=false,
    show_colorbar::Bool=true,
    clims=nothing,
    phase_label::Bool=true,
    title=nothing,
    phi=nothing,
    uγ=nothing,
    marker_size::Real=8,
    kwargs...,
)
    y = snapshot.y === nothing ? throw(ArgumentError("2D scalar plotting requires snapshot.y")) : collect(snapshot.y)
    x = collect(snapshot.x)

    vals = _coerce_values_to_shape(snapshot.values, (length(x), length(y)))
    mask = _coerce_mask(snapshot.mask, size(vals))
    z = _masked_array(vals, mask)

    hm = if clims === nothing
        heatmap!(ax, x, y, z; colormap=colormap, kwargs...)
    else
        heatmap!(ax, x, y, z; colormap=colormap, colorrange=clims, kwargs...)
    end
    ax.aspect = DataAspect()
    ax.xlabel = "x"
    ax.ylabel = "y"

    if show_cutcells
        xc = Float64[]
        yc = Float64[]
        @inbounds for j in eachindex(y), i in eachindex(x)
            if !mask[i, j]
                push!(xc, Float64(x[i]))
                push!(yc, Float64(y[j]))
            end
        end
        isempty(xc) || scatter!(ax, xc, yc; color=(:black, 0.35), markersize=2)
    end

    if show_interface
        if phi !== nothing
            phi_arr = _coerce_values_to_shape(phi, (length(x), length(y)))
            contour!(ax, x, y, phi_arr; levels=[0.0], color=:black, linewidth=2)
        elseif snapshot.interface_xy !== nothing
            xγ = collect(snapshot.interface_xy[1])
            yγ = snapshot.interface_xy[2]
            if yγ === nothing
                isempty(xγ) || vlines!(ax, xγ; color=:black, linewidth=1.5, linestyle=:dash)
            else
                yγv = collect(yγ)
                n = min(length(xγ), length(yγv))
                n > 0 && scatter!(ax, xγ[1:n], yγv[1:n]; color=:black, markersize=3)
            end
        end
    end

    if uγ !== nothing && snapshot.interface_xy !== nothing && snapshot.interface_xy[2] !== nothing
        xγ = collect(snapshot.interface_xy[1])
        yγ = collect(snapshot.interface_xy[2])
        uγv = vec(collect(uγ))
        n = min(min(length(xγ), length(yγ)), length(uγv))
        n > 0 && scatter!(ax, xγ[1:n], yγ[1:n]; color=uγv[1:n], colormap=:balance, markersize=marker_size)
    end

    if show_colorbar
        _add_colorbar!(fig, hm; label=snapshot.colorlabel)
    end

    ttl = title === nothing ? _default_title(snapshot.title, snapshot.time, snapshot.phase; phase_label=phase_label) : String(title)
    ax.title = ttl
    return (fig=fig, ax=ax, plot=hm)
end

function _plot_scalar_field_impl(snapshot::ScalarSnapshot; fig=nothing, ax=nothing, kwargs...)
    fig_local, ax_local, _ = _ensure_fig_ax(; fig=fig, ax=ax)
    if _scalar_mode_1d(snapshot)
        return _plot_scalar_1d!(snapshot, fig_local, ax_local; kwargs...)
    end
    return _plot_scalar_2d!(snapshot, fig_local, ax_local; kwargs...)
end

function plot_scalar_field(
    snapshot::ScalarSnapshot;
    fig=nothing,
    ax=nothing,
    colormap=:viridis,
    show_interface::Bool=true,
    show_cutcells::Bool=false,
    show_colorbar::Bool=true,
    clims=nothing,
    phase_label::Bool=true,
    title=nothing,
    phi=nothing,
    uγ=nothing,
    marker_size::Real=8,
    kwargs...,
)
    phi_data = phi === nothing ? get(snapshot.metadata, :phi, nothing) : phi
    out = _plot_scalar_field_impl(
        snapshot;
        fig=fig,
        ax=ax,
        colormap=colormap,
        show_interface=show_interface,
        show_cutcells=show_cutcells,
        show_colorbar=show_colorbar,
        clims=clims,
        phase_label=phase_label,
        title=title,
        phi=phi_data,
        uγ=uγ,
        marker_size=marker_size,
        kwargs...,
    )
    return out.fig
end

function plot_scalar_field(
    snapshots::AbstractVector{<:ScalarSnapshot};
    fig=nothing,
    colormap=:viridis,
    show_interface::Bool=true,
    show_cutcells::Bool=false,
    show_colorbar::Bool=true,
    clims=nothing,
    shared_clims::Bool=true,
    phase_label::Bool=true,
    kwargs...,
)
    isempty(snapshots) && throw(ArgumentError("expected at least one snapshot"))

    n = length(snapshots)
    fig_local = fig === nothing ? Figure(size=(500 * n, 500)) : fig

    clims_eff = clims
    if shared_clims && clims === nothing
        arrays = Float64[]
        buffers = Vector{Matrix{Float64}}()
        for snap in snapshots
            if !_scalar_mode_1d(snap)
                x = collect(snap.x)
                y = collect(snap.y)
                vals = _coerce_values_to_shape(snap.values, (length(x), length(y)))
                m = _coerce_mask(snap.mask, size(vals))
                push!(buffers, _masked_array(vals, m))
            end
        end
        !isempty(buffers) && (clims_eff = _auto_clims(buffers))
    end

    first_plot = nothing
    for (i, snap) in enumerate(snapshots)
        ax = Axis(fig_local[1, i])
        ret = _plot_scalar_field_impl(
            snap;
            fig=fig_local,
            ax=ax,
            colormap=colormap,
            show_interface=show_interface,
            show_cutcells=show_cutcells,
            show_colorbar=false,
            clims=clims_eff,
            phase_label=phase_label,
            kwargs...,
        )
        first_plot === nothing && (first_plot = ret.plot)
    end

    if show_colorbar && first_plot !== nothing
        _add_colorbar!(fig_local, first_plot; label=snapshots[1].colorlabel)
    end
    return fig_local
end

plot_scalar_field(snaps::Tuple{Vararg{ScalarSnapshot}}; kwargs...) = plot_scalar_field(collect(snaps); kwargs...)

function plot_interface_field(
    snapshot::InterfaceSnapshot;
    fig=nothing,
    ax=nothing,
    style::Symbol=:scatter,
    signed_coloring::Bool=true,
    markersize::Real=8,
    linewidth::Real=2,
    color=:black,
    title=nothing,
    phase_label::Bool=true,
    kwargs...,
)
    fig_local, ax_local, _ = _ensure_fig_ax(; fig=fig, ax=ax)
    x = collect(snapshot.xγ)
    vals = collect(snapshot.valuesγ)

    if snapshot.yγ === nothing
        n = min(length(x), length(vals))
        if style === :line
            lines!(ax_local, x[1:n], vals[1:n]; color=color, linewidth=linewidth, kwargs...)
        else
            c = signed_coloring ? vals[1:n] : color
            scatter!(ax_local, x[1:n], vals[1:n]; color=c, colormap=:balance, markersize=markersize, kwargs...)
        end
        ax_local.xlabel = "x"
        ax_local.ylabel = snapshot.label
    else
        y = collect(snapshot.yγ)
        n = min(min(length(x), length(y)), length(vals))
        if style === :line
            lines!(ax_local, x[1:n], y[1:n]; color=color, linewidth=linewidth, kwargs...)
        else
            c = signed_coloring ? vals[1:n] : color
            scatter!(ax_local, x[1:n], y[1:n]; color=c, colormap=:balance, markersize=markersize, kwargs...)
        end
        ax_local.xlabel = "x"
        ax_local.ylabel = "y"
        ax_local.aspect = DataAspect()
    end

    ttl = title === nothing ? _default_title(snapshot.label, snapshot.time, snapshot.phase; phase_label=phase_label) : String(title)
    ax_local.title = ttl
    return fig_local
end
