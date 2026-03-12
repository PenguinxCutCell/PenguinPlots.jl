module PenguinPlotsStokesExt

using Makie
using PenguinPlots
using PenguinStokes

import PenguinPlots: plot_history, plot_solution

function _statevec(x)
    if x isa AbstractVector
        return x
    elseif x isa NamedTuple
        if haskey(x, :states)
            states = getproperty(x, :states)
            isempty(states) && throw(ArgumentError("states history is empty"))
            return states[end]
        elseif haskey(x, :system)
            return getproperty(getproperty(x, :system), :x)
        elseif haskey(x, :x)
            return getproperty(x, :x)
        end
    elseif hasproperty(x, :x)
        return getproperty(x, :x)
    end
    throw(ArgumentError("cannot extract state vector from $(typeof(x))"))
end

function _expand_mono_state(layout, xvec::AbstractVector)
    nsys = PenguinStokes.nunknowns(layout)
    nt = layout.nt
    N = length(layout.uomega)

    if length(xvec) == nsys
        return collect(xvec)
    elseif length(xvec) == N * nt
        out = zeros(eltype(xvec), nsys)
        for d in 1:N
            src = ((d - 1) * nt + 1):(d * nt)
            out[layout.uomega[d]] .= xvec[src]
        end
        return out
    end
    throw(DimensionMismatch("state length must be $nsys (full) or $(N * nt) (uomega blocks)"))
end

function _expand_twophase_state(layout, xvec::AbstractVector)
    nsys = PenguinStokes.nunknowns(layout)
    length(xvec) == nsys || throw(DimensionMismatch("two-phase state length must be $nsys"))
    return collect(xvec)
end

function _average_staggered_to_center(uω::NTuple{N,AbstractVector}, dims::NTuple{N,Int}) where {N}
    arrays = ntuple(d -> reshape(collect(uω[d]), dims), N)

    if N == 1
        u = arrays[1]
        v = zeros(eltype(u), size(u))
        return u, v
    elseif N >= 2
        nx, ny = dims[1], dims[2]
        u = Array{Float64}(undef, nx, ny)
        v = Array{Float64}(undef, nx, ny)

        ux = arrays[1]
        vy = arrays[2]
        @inbounds for j in 1:ny, i in 1:nx
            u[i, j] = i < nx ? 0.5 * (ux[i, j] + ux[i + 1, j]) : Float64(ux[i, j])
            v[i, j] = j < ny ? 0.5 * (vy[i, j] + vy[i, j + 1]) : Float64(vy[i, j])
        end
        return u, v
    end

    throw(ArgumentError("unsupported dimension N=$N"))
end

function _mono_snapshot(model::PenguinStokes.StokesModelMono, xvec::AbstractVector; title::AbstractString="Stokes mono")
    full = _expand_mono_state(model.layout, xvec)
    N = length(model.layout.uomega)

    pω = full[model.layout.pomega]
    uω = ntuple(d -> full[model.layout.uomega[d]], N)

    dims = PenguinPlots.capacity_dims(model.cap_p)
    p = PenguinPlots.reshape_on_capacity(pω, model.cap_p)
    u, v = _average_staggered_to_center(uω, dims)

    axes = PenguinPlots.capacity_axes(model.cap_p)
    x = axes[1]
    y = length(axes) >= 2 ? axes[2] : [0.0]
    mask = PenguinPlots.active_mask_from_capacity(model.cap_p)
    iface = PenguinPlots.interface_coordinates_from_capacity(model.cap_p)

    return PenguinPlots.VectorSnapshot(
        x,
        y,
        u,
        v;
        p=p,
        mask=mask,
        interface_xy=iface,
        title=title,
        metadata=Dict{Symbol,Any}(:dimension => N),
    )
end

function _moving_mono_snapshot(model::PenguinStokes.MovingStokesModelMono, xvec::AbstractVector; title::AbstractString="Moving Stokes mono")
    full = _expand_mono_state(model.layout, xvec)
    N = length(model.layout.uomega)

    cap_p = something(model.cap_p_end, model.cap_p_slab, nothing)
    cap_p === nothing && throw(ArgumentError("moving Stokes plotting requires available slab/end capacities; run an assembly/solve first"))

    pω = full[model.layout.pomega]
    uω = ntuple(d -> full[model.layout.uomega[d]], N)

    dims = PenguinPlots.capacity_dims(cap_p)
    p = PenguinPlots.reshape_on_capacity(pω, cap_p)
    u, v = _average_staggered_to_center(uω, dims)

    axes = PenguinPlots.capacity_axes(cap_p)
    x = axes[1]
    y = length(axes) >= 2 ? axes[2] : [0.0]
    mask = PenguinPlots.active_mask_from_capacity(cap_p)
    iface = PenguinPlots.interface_coordinates_from_capacity(cap_p)

    return PenguinPlots.VectorSnapshot(
        x,
        y,
        u,
        v;
        p=p,
        mask=mask,
        interface_xy=iface,
        title=title,
        metadata=Dict{Symbol,Any}(:dimension => N),
    )
end

function _twophase_snapshot(
    model::PenguinStokes.StokesModelTwoPhase,
    xvec::AbstractVector;
    phase::Int=1,
    title::AbstractString="Stokes two-phase",
)
    phase in (1, 2) || throw(ArgumentError("phase must be 1 or 2"))
    full = _expand_twophase_state(model.layout, xvec)
    N = length(model.layout.uomega1)

    if phase == 1
        pω = full[model.layout.pomega1]
        uω = ntuple(d -> full[model.layout.uomega1[d]], N)
        cap = model.cap_p1
    else
        pω = full[model.layout.pomega2]
        uω = ntuple(d -> full[model.layout.uomega2[d]], N)
        cap = model.cap_p2
    end

    dims = PenguinPlots.capacity_dims(cap)
    p = PenguinPlots.reshape_on_capacity(pω, cap)
    u, v = _average_staggered_to_center(uω, dims)

    axes = PenguinPlots.capacity_axes(cap)
    x = axes[1]
    y = length(axes) >= 2 ? axes[2] : [0.0]
    mask = PenguinPlots.active_mask_from_capacity(cap)
    iface = PenguinPlots.interface_coordinates_from_capacity(cap)

    return PenguinPlots.VectorSnapshot(
        x,
        y,
        u,
        v;
        p=p,
        mask=mask,
        interface_xy=iface,
        title=title,
        phase=phase,
        metadata=Dict{Symbol,Any}(:dimension => N),
    )
end

function _plot_stokes_snapshot(snapshot::PenguinPlots.VectorSnapshot; field::Symbol=:pressure_velocity, kwargs...)
    if field === :pressure || field === :p
        return PenguinPlots.plot_pressure_field(snapshot; kwargs...)
    elseif field === :velocity || field === :u
        return PenguinPlots.plot_vector_field(snapshot; show_pressure=false, kwargs...)
    elseif field === :pressure_velocity || field === :combined
        return PenguinPlots.plot_vector_field(snapshot; show_pressure=true, kwargs...)
    end
    throw(ArgumentError("unsupported Stokes field `$field`"))
end

function plot_solution(
    model::PenguinStokes.StokesModelMono,
    x;
    field::Symbol=:pressure,
    kwargs...,
)
    snap = _mono_snapshot(model, _statevec(x))
    return _plot_stokes_snapshot(snap; field=field, kwargs...)
end

function plot_solution(
    model::PenguinStokes.MovingStokesModelMono,
    x;
    field::Symbol=:pressure_velocity,
    kwargs...,
)
    snap = _moving_mono_snapshot(model, _statevec(x))
    return _plot_stokes_snapshot(snap; field=field, kwargs...)
end

function plot_solution(
    model::PenguinStokes.StokesModelTwoPhase,
    x;
    field::Symbol=:pressure_velocity,
    phase::Union{Int,Symbol}=:both,
    shared_clims::Bool=true,
    kwargs...,
)
    xvec = _statevec(x)

    if phase === :both
        s1 = _twophase_snapshot(model, xvec; phase=1, title="Stokes phase 1")
        s2 = _twophase_snapshot(model, xvec; phase=2, title="Stokes phase 2")

        if field === :pressure || field === :p
            fig = Figure(size=(1200, 550))
            clims = nothing
            if shared_clims
                m1 = PenguinPlots.nanfilled(s1.p, s1.mask)
                m2 = PenguinPlots.nanfilled(s2.p, s2.mask)
                clims = PenguinPlots._auto_clims([m1, m2])
            end
            PenguinPlots.plot_pressure_field(s1; fig=fig, ax=Axis(fig[1, 1]), show_colorbar=false, clims=clims, title="phase 1", kwargs...)
            PenguinPlots.plot_pressure_field(s2; fig=fig, ax=Axis(fig[1, 2]), show_colorbar=false, clims=clims, title="phase 2", kwargs...)
            return fig
        elseif field === :velocity || field === :u
            fig = Figure(size=(1200, 550))
            PenguinPlots.plot_vector_field(s1; fig=fig, ax=Axis(fig[1, 1]), show_pressure=false, title="phase 1", kwargs...)
            PenguinPlots.plot_vector_field(s2; fig=fig, ax=Axis(fig[1, 2]), show_pressure=false, title="phase 2", kwargs...)
            return fig
        elseif field === :pressure_velocity || field === :combined
            fig = Figure(size=(1300, 560))
            PenguinPlots.plot_vector_field(s1; fig=fig, ax=Axis(fig[1, 1]), show_pressure=true, title="phase 1", kwargs...)
            PenguinPlots.plot_vector_field(s2; fig=fig, ax=Axis(fig[1, 2]), show_pressure=true, title="phase 2", kwargs...)
            return fig
        end
        throw(ArgumentError("unsupported Stokes two-phase field `$field`"))
    end

    snap = _twophase_snapshot(model, xvec; phase=Int(phase))
    return _plot_stokes_snapshot(snap; field=field, kwargs...)
end

function plot_history(
    history::AbstractVector{<:NamedTuple};
    field::Symbol=:residual,
    component::Int=1,
    kwargs...,
)
    isempty(history) && throw(ArgumentError("history is empty"))
    haskey(history[1], :t) || throw(ArgumentError("history entries must contain `t`"))
    haskey(history[1], field) || throw(ArgumentError("history entries do not contain `$field`"))

    t = [getproperty(h, :t) for h in history]
    raw = [getproperty(h, field) for h in history]

    if raw[1] isa Number
        values = raw
    elseif raw[1] isa AbstractVector || raw[1] isa Tuple
        values = [ri[component] for ri in raw]
    else
        throw(ArgumentError("unsupported history field value type $(typeof(raw[1]))"))
    end

    ylabel = component == 1 ? String(field) : "$(field)[$component]"
    return PenguinPlots.plot_history(t, values; ylabel=ylabel, title="Stokes history: $(field)", kwargs...)
end

end # module
