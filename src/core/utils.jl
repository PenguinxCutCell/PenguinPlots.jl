const _DEFAULT_FIGSIZE = (900, 500)

function _figure_from_axis(ax)
    if hasproperty(ax, :figure)
        return getproperty(ax, :figure)
    end
    try
        return Makie.parent_figure(ax)
    catch
        return nothing
    end
end

function _ensure_fig_ax(
    ;
    fig=nothing,
    ax=nothing,
    figsize=_DEFAULT_FIGSIZE,
    axis_kwargs...,
)
    if ax !== nothing
        fig_local = fig === nothing ? _figure_from_axis(ax) : fig
        return fig_local, ax, false
    end
    fig_local = fig === nothing ? Figure(size=figsize) : fig
    ax_local = Axis(fig_local[1, 1]; axis_kwargs...)
    return fig_local, ax_local, true
end

function _coerce_values_to_shape(values::AbstractArray, dims::NTuple)
    size(values) == dims && return Array(values)
    length(values) == prod(dims) || throw(DimensionMismatch("cannot reshape values of size $(size(values)) into $dims"))
    return reshape(collect(values), dims)
end

function _coerce_mask(mask::AbstractArray{Bool}, dims::NTuple)
    size(mask) == dims && return Array(mask)
    length(mask) == prod(dims) || throw(DimensionMismatch("cannot reshape mask of size $(size(mask)) into $dims"))
    return reshape(collect(mask), dims)
end

function _masked_array(values::AbstractArray, mask::AbstractArray{Bool})
    dims = size(values)
    m = _coerce_mask(mask, dims)
    out = Array{Float64}(undef, dims)
    @inbounds for i in eachindex(values)
        out[i] = m[i] ? Float64(values[i]) : NaN
    end
    return out
end

function _auto_clims(arrays::AbstractVector{<:AbstractArray})
    minv = Inf
    maxv = -Inf
    for arr in arrays
        @inbounds for v in arr
            isfinite(v) || continue
            minv = min(minv, v)
            maxv = max(maxv, v)
        end
    end
    if !isfinite(minv) || !isfinite(maxv)
        return nothing
    end
    if minv == maxv
        δ = max(abs(minv), 1.0) * 1e-6
        return (minv - δ, maxv + δ)
    end
    return (minv, maxv)
end

function _default_title(base::String, time, phase::Union{Nothing,Int}; phase_label::Bool=true)
    parts = String[]
    !isempty(base) && push!(parts, base)
    if phase_label && phase !== nothing
        push!(parts, "phase $phase")
    end
    if time !== nothing
        push!(parts, "t=$(time)")
    end
    return join(parts, " | ")
end

function _add_colorbar!(fig, hm; label::AbstractString="")
    fig === nothing && return nothing
    if hasproperty(fig, :layout)
        return Colorbar(fig[1, end + 1], hm; label=label)
    end
    return Colorbar(fig, hm; label=label)
end

function _clear_axis!(ax)
    try
        Makie.empty!(ax)
    catch
        if hasproperty(ax, :scene) && hasproperty(ax.scene, :plots)
            for p in reverse(copy(ax.scene.plots))
                try
                    delete!(ax.scene, p)
                catch
                end
            end
        end
    end
    return ax
end
