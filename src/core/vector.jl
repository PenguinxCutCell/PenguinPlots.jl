function _vector_mode_1d(snapshot::VectorSnapshot)
    return length(snapshot.y) <= 1 || ndims(snapshot.u) == 1
end

function _vector_title(snapshot::VectorSnapshot, title, phase_label::Bool)
    if title !== nothing
        return String(title)
    end
    return _default_title(snapshot.title, snapshot.time, snapshot.phase; phase_label=phase_label)
end

function _coerce_vector_data(snapshot::VectorSnapshot)
    nx = length(snapshot.x)
    ny = max(length(snapshot.y), 1)

    if _vector_mode_1d(snapshot)
        u = vec(collect(snapshot.u))
        v = vec(collect(snapshot.v))
        mask = vec(_coerce_mask(snapshot.mask, size(u)))
        p = snapshot.p === nothing ? nothing : vec(collect(snapshot.p))
        return u, v, p, mask
    end

    dims = (nx, ny)
    u = _coerce_values_to_shape(snapshot.u, dims)
    v = _coerce_values_to_shape(snapshot.v, dims)
    p = snapshot.p === nothing ? nothing : _coerce_values_to_shape(snapshot.p, dims)
    mask = _coerce_mask(snapshot.mask, dims)
    return u, v, p, mask
end

function _plot_vector_1d!(
    snapshot::VectorSnapshot,
    fig,
    ax;
    title=nothing,
    phase_label::Bool=true,
    linewidth::Real=2,
    kwargs...,
)
    x = collect(snapshot.x)
    u, _, p, mask = _coerce_vector_data(snapshot)

    y = similar(u, Float64)
    @inbounds for i in eachindex(u)
        y[i] = mask[i] ? Float64(u[i]) : NaN
    end
    lines!(ax, x, y; linewidth=linewidth, color=:steelblue, label="u", kwargs...)

    if p !== nothing
        yp = similar(p, Float64)
        @inbounds for i in eachindex(p)
            yp[i] = mask[i] ? Float64(p[i]) : NaN
        end
        lines!(ax, x, yp; linewidth=linewidth, color=:firebrick, linestyle=:dash, label="p")
        axislegend(ax; position=:rt)
    end

    ax.xlabel = "x"
    ax.ylabel = "value"
    ax.title = _vector_title(snapshot, title, phase_label)
    return (fig=fig, ax=ax, plot=nothing)
end

function _plot_vector_2d!(
    snapshot::VectorSnapshot,
    fig,
    ax;
    stride::Int=2,
    velocity_scale::Real=1.0,
    show_pressure::Bool=true,
    show_colorbar::Bool=true,
    colormap=:viridis,
    clims=nothing,
    show_interface::Bool=true,
    title=nothing,
    phase_label::Bool=true,
    arrow_color=:white,
    marker_size::Real=4,
    kwargs...,
)
    stride >= 1 || throw(ArgumentError("stride must be >= 1"))

    x = collect(snapshot.x)
    y = collect(snapshot.y)
    u, v, p, mask = _coerce_vector_data(snapshot)

    pplot = nothing
    if show_pressure && p !== nothing
        pm = _masked_array(p, mask)
        pplot = if clims === nothing
            heatmap!(ax, x, y, pm; colormap=colormap)
        else
            heatmap!(ax, x, y, pm; colormap=colormap, colorrange=clims)
        end
        if show_colorbar
            _add_colorbar!(fig, pplot; label="pressure")
        end
    end

    pts = Point2f[]
    vecs = Vec2f[]
    @inbounds for j in 1:stride:length(y), i in 1:stride:length(x)
        if mask[i, j]
            ui = Float32(velocity_scale * u[i, j])
            vi = Float32(velocity_scale * v[i, j])
            if isfinite(ui) && isfinite(vi)
                push!(pts, Point2f(Float32(x[i]), Float32(y[j])))
                push!(vecs, Vec2f(ui, vi))
            end
        end
    end
    isempty(pts) || arrows!(ax, pts, vecs; color=arrow_color)

    if show_interface && snapshot.interface_xy !== nothing
        xγ = collect(snapshot.interface_xy[1])
        yγ = snapshot.interface_xy[2]
        if yγ === nothing
            isempty(xγ) || vlines!(ax, xγ; color=:black, linewidth=1.5, linestyle=:dash)
        else
            yγv = collect(yγ)
            n = min(length(xγ), length(yγv))
            n > 0 && scatter!(ax, xγ[1:n], yγv[1:n]; color=:black, markersize=marker_size)
        end
    end

    ax.aspect = DataAspect()
    ax.xlabel = "x"
    ax.ylabel = "y"
    ax.title = _vector_title(snapshot, title, phase_label)

    return (fig=fig, ax=ax, plot=pplot)
end

function _plot_vector_field_impl(snapshot::VectorSnapshot; fig=nothing, ax=nothing, kwargs...)
    fig_local, ax_local, _ = _ensure_fig_ax(; fig=fig, ax=ax)
    if _vector_mode_1d(snapshot)
        return _plot_vector_1d!(snapshot, fig_local, ax_local; kwargs...)
    end
    return _plot_vector_2d!(snapshot, fig_local, ax_local; kwargs...)
end

function plot_vector_field(
    snapshot::VectorSnapshot;
    fig=nothing,
    ax=nothing,
    stride::Int=2,
    velocity_scale::Real=1.0,
    show_pressure::Bool=true,
    show_colorbar::Bool=true,
    colormap=:viridis,
    clims=nothing,
    show_interface::Bool=true,
    title=nothing,
    phase_label::Bool=true,
    kwargs...,
)
    out = _plot_vector_field_impl(
        snapshot;
        fig=fig,
        ax=ax,
        stride=stride,
        velocity_scale=velocity_scale,
        show_pressure=show_pressure,
        show_colorbar=show_colorbar,
        colormap=colormap,
        clims=clims,
        show_interface=show_interface,
        title=title,
        phase_label=phase_label,
        kwargs...,
    )
    return out.fig
end

function plot_pressure_field(
    snapshot::VectorSnapshot;
    fig=nothing,
    ax=nothing,
    colormap=:viridis,
    show_interface::Bool=true,
    show_colorbar::Bool=true,
    clims=nothing,
    title=nothing,
    phase_label::Bool=true,
    kwargs...,
)
    snapshot.p === nothing && throw(ArgumentError("pressure data is not available in this vector snapshot"))
    snap = ScalarSnapshot(
        snapshot.x,
        _vector_mode_1d(snapshot) ? nothing : snapshot.y,
        snapshot.p,
        snapshot.mask;
        title=isempty(snapshot.title) ? "pressure" : snapshot.title,
        colorlabel="pressure",
        time=snapshot.time,
        phase=snapshot.phase,
        interface_xy=snapshot.interface_xy,
        metadata=copy(snapshot.metadata),
    )
    return plot_scalar_field(
        snap;
        fig=fig,
        ax=ax,
        colormap=colormap,
        show_interface=show_interface,
        show_colorbar=show_colorbar,
        clims=clims,
        title=title,
        phase_label=phase_label,
        kwargs...,
    )
end

function plot_pressure_field(
    x::AbstractVector,
    y::Union{Nothing,AbstractVector},
    p::AbstractArray;
    mask::AbstractArray{Bool}=trues(size(p)),
    kwargs...,
)
    snap = ScalarSnapshot(x, y, p, mask; title="pressure", colorlabel="pressure")
    return plot_scalar_field(snap; kwargs...)
end
