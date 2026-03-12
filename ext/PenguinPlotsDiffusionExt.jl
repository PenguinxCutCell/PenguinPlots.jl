module PenguinPlotsDiffusionExt

using PenguinDiffusion
using PenguinPlots

import PenguinPlots: plot_solution

function _statevec(x)
    if x isa AbstractVector
        return x
    elseif x isa NamedTuple
        if haskey(x, :states)
            states = getproperty(x, :states)
            isempty(states) && throw(ArgumentError("states history is empty"))
            return states[end]
        elseif haskey(x, :system)
            sys = getproperty(x, :system)
            hasproperty(sys, :x) || throw(ArgumentError("`system` does not expose field `x`"))
            return getproperty(sys, :x)
        elseif haskey(x, :x)
            return getproperty(x, :x)
        end
    elseif hasproperty(x, :x)
        return getproperty(x, :x)
    end
    throw(ArgumentError("cannot extract state vector from object of type $(typeof(x))"))
end

function _expand_mono_state(model, xvec::AbstractVector)
    lay = model.layout.offsets
    nt = length(lay.ω)
    nsys = max(last(lay.ω), last(lay.γ))
    if length(xvec) == nsys
        return collect(xvec)
    elseif length(xvec) == nt
        out = zeros(eltype(xvec), nsys)
        out[lay.ω] .= xvec
        out[lay.γ] .= xvec
        return out
    end
    throw(DimensionMismatch("state length must be $nt (ω) or $nsys (ω+γ)"))
end

function _expand_diph_state(model, xvec::AbstractVector)
    lay = model.layout.offsets
    nt = length(lay.ω1)
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    if length(xvec) == nsys
        return collect(xvec)
    elseif length(xvec) == 2 * nt
        out = zeros(eltype(xvec), nsys)
        out[lay.ω1] .= xvec[1:nt]
        out[lay.ω2] .= xvec[(nt + 1):(2 * nt)]
        out[lay.γ1] .= out[lay.ω1]
        out[lay.γ2] .= out[lay.ω2]
        return out
    end
    throw(DimensionMismatch("state length must be $(2 * nt) (ω1+ω2) or $nsys (full system)"))
end

function _mask_and_coords(cap)
    axes = PenguinPlots.capacity_axes(cap)
    x = axes[1]
    y = length(axes) >= 2 ? axes[2] : nothing
    mask = PenguinPlots.active_mask_from_capacity(cap)
    iface = PenguinPlots.interface_coordinates_from_capacity(cap)
    return x, y, mask, iface
end

function _interface_values(cap, γ::AbstractVector)
    im = vec(PenguinPlots.interface_mask_from_capacity(cap))
    γv = vec(collect(γ))
    vals = Float64[]
    n = min(length(im), length(γv))
    @inbounds for i in 1:n
        im[i] || continue
        push!(vals, Float64(γv[i]))
    end
    return vals
end

function _mono_snapshot(model::PenguinDiffusion.DiffusionModelMono, xvec::AbstractVector; field::Symbol=:omega, title::AbstractString="Diffusion mono")
    full = _expand_mono_state(model, xvec)
    lay = model.layout.offsets
    cap = model.cap

    ω = full[lay.ω]
    x, y, mask, iface = _mask_and_coords(cap)
    vals = PenguinPlots.reshape_on_capacity(ω, cap)

    return PenguinPlots.ScalarSnapshot(
        x,
        y,
        vals,
        mask;
        title=title,
        colorlabel=field === :omega ? "ω" : String(field),
        interface_xy=iface,
        metadata=Dict{Symbol,Any}(),
    ), full[lay.γ], cap
end

function _diph_snapshot(
    model::PenguinDiffusion.DiffusionModelDiph,
    xvec::AbstractVector;
    phase::Int=1,
    title::AbstractString="Diffusion diph",
)
    phase in (1, 2) || throw(ArgumentError("phase must be 1 or 2"))
    full = _expand_diph_state(model, xvec)
    lay = model.layout.offsets

    if phase == 1
        cap = model.cap1
        ω = full[lay.ω1]
        γ = full[lay.γ1]
    else
        cap = model.cap2
        ω = full[lay.ω2]
        γ = full[lay.γ2]
    end

    x, y, mask, iface = _mask_and_coords(cap)
    vals = PenguinPlots.reshape_on_capacity(ω, cap)

    snap = PenguinPlots.ScalarSnapshot(
        x,
        y,
        vals,
        mask;
        title=title,
        colorlabel="ω$phase",
        phase=phase,
        interface_xy=iface,
        metadata=Dict{Symbol,Any}(),
    )
    return snap, γ, cap
end

function _moving_axes_mask(model)
    axes = PenguinPlots.cell_center_axes_from_grid(model.grid)
    dims = Tuple(getproperty(model.grid, :n))
    x = axes[1]
    y = length(axes) >= 2 ? axes[2] : nothing
    return x, y, trues(dims), nothing
end

function _moving_mono_snapshot(
    model::PenguinDiffusion.MovingDiffusionModelMono,
    xvec::AbstractVector;
    title::AbstractString="Moving diffusion mono",
)
    full = _expand_mono_state(model, xvec)
    lay = model.layout.offsets
    ω = full[lay.ω]
    γ = full[lay.γ]

    cap = something(model.cap_slab, nothing)
    if cap === nothing
        x, y, mask, iface = _moving_axes_mask(model)
        vals = PenguinPlots._coerce_values_to_shape(ω, Tuple(getproperty(model.grid, :n)))
    else
        x, y, mask, iface = _mask_and_coords(cap)
        vals = PenguinPlots.reshape_on_capacity(ω, cap)
    end

    snap = PenguinPlots.ScalarSnapshot(
        x,
        y,
        vals,
        mask;
        title=title,
        colorlabel="ω",
        interface_xy=iface,
        metadata=Dict{Symbol,Any}(),
    )
    return snap, γ, cap
end

function _moving_diph_snapshot(
    model::PenguinDiffusion.MovingDiffusionModelDiph,
    xvec::AbstractVector;
    phase::Int=1,
    title::AbstractString="Moving diffusion diph",
)
    phase in (1, 2) || throw(ArgumentError("phase must be 1 or 2"))
    full = _expand_diph_state(model, xvec)
    lay = model.layout.offsets

    if phase == 1
        ω = full[lay.ω1]
        γ = full[lay.γ1]
        cap = something(model.cap1_slab, nothing)
    else
        ω = full[lay.ω2]
        γ = full[lay.γ2]
        cap = something(model.cap2_slab, nothing)
    end

    if cap === nothing
        x, y, mask, iface = _moving_axes_mask(model)
        vals = PenguinPlots._coerce_values_to_shape(ω, Tuple(getproperty(model.grid, :n)))
    else
        x, y, mask, iface = _mask_and_coords(cap)
        vals = PenguinPlots.reshape_on_capacity(ω, cap)
    end

    snap = PenguinPlots.ScalarSnapshot(
        x,
        y,
        vals,
        mask;
        title=title,
        colorlabel="ω$phase",
        phase=phase,
        interface_xy=iface,
        metadata=Dict{Symbol,Any}(),
    )
    return snap, γ, cap
end

function _interface_snapshot(cap, γ; label::AbstractString="uγ", phase=nothing)
    cap === nothing && return PenguinPlots.InterfaceSnapshot(Float64[], nothing, Float64[]; label=label, phase=phase)

    im = vec(PenguinPlots.interface_mask_from_capacity(cap))
    Cγ = getproperty(cap, :C_γ)
    γv = vec(collect(γ))
    n = min(min(length(im), length(Cγ)), length(γv))

    xs = Float64[]
    ys = Float64[]
    vals = Float64[]
    for i in 1:n
        im[i] || continue
        push!(xs, Float64(Cγ[i][1]))
        if length(Cγ[i]) >= 2
            push!(ys, Float64(Cγ[i][2]))
        end
        push!(vals, Float64(γv[i]))
    end

    if isempty(ys)
        return PenguinPlots.InterfaceSnapshot(xs, nothing, vals; label=label, phase=phase)
    end
    return PenguinPlots.InterfaceSnapshot(xs, ys, vals; label=label, phase=phase)
end

function plot_solution(
    model::PenguinDiffusion.DiffusionModelMono,
    x;
    field::Symbol=:omega,
    show_interfacial::Bool=false,
    kwargs...,
)
    xvec = _statevec(x)
    snap, γ, cap = _mono_snapshot(model, xvec; field=field)

    if field === :omega || field === :bulk
        uγ = show_interfacial ? _interface_values(cap, γ) : nothing
        return PenguinPlots.plot_scalar_field(snap; uγ=uγ, kwargs...)
    elseif field === :gamma || field === :interface || field === :uγ
        isnap = _interface_snapshot(cap, γ; label="uγ")
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    end

    throw(ArgumentError("unsupported diffusion field `$field`"))
end

function plot_solution(
    model::PenguinDiffusion.DiffusionModelDiph,
    x;
    phase::Union{Int,Symbol}=1,
    field::Symbol=:omega,
    show_interfacial::Bool=false,
    kwargs...,
)
    xvec = _statevec(x)

    if phase === :both
        field === :omega || field === :bulk || throw(ArgumentError("phase=:both currently supports field=:omega/:bulk"))
        s1, _, _ = _diph_snapshot(model, xvec; phase=1)
        s2, _, _ = _diph_snapshot(model, xvec; phase=2)
        return PenguinPlots.plot_scalar_field([s1, s2]; kwargs...)
    end

    phase_i = Int(phase)
    snap, γ, cap = _diph_snapshot(model, xvec; phase=phase_i)

    if field === :omega || field === :bulk
        uγ = show_interfacial ? _interface_values(cap, γ) : nothing
        return PenguinPlots.plot_scalar_field(snap; uγ=uγ, kwargs...)
    elseif field === :gamma || field === :interface || field === :uγ
        isnap = _interface_snapshot(cap, γ; label="uγ$phase_i", phase=phase_i)
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    end

    throw(ArgumentError("unsupported diffusion field `$field`"))
end

function plot_solution(
    model::PenguinDiffusion.MovingDiffusionModelMono,
    x;
    field::Symbol=:omega,
    show_interfacial::Bool=false,
    kwargs...,
)
    xvec = _statevec(x)
    snap, γ, cap = _moving_mono_snapshot(model, xvec)

    if field === :omega || field === :bulk
        uγ = (show_interfacial && cap !== nothing) ? _interface_values(cap, γ) : nothing
        return PenguinPlots.plot_scalar_field(snap; uγ=uγ, kwargs...)
    elseif field === :gamma || field === :interface || field === :uγ
        isnap = _interface_snapshot(cap, γ; label="uγ")
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    end

    throw(ArgumentError("unsupported diffusion field `$field`"))
end

function plot_solution(
    model::PenguinDiffusion.MovingDiffusionModelDiph,
    x;
    phase::Union{Int,Symbol}=1,
    field::Symbol=:omega,
    show_interfacial::Bool=false,
    kwargs...,
)
    xvec = _statevec(x)

    if phase === :both
        field === :omega || field === :bulk || throw(ArgumentError("phase=:both currently supports field=:omega/:bulk"))
        s1, _, _ = _moving_diph_snapshot(model, xvec; phase=1)
        s2, _, _ = _moving_diph_snapshot(model, xvec; phase=2)
        return PenguinPlots.plot_scalar_field([s1, s2]; kwargs...)
    end

    phase_i = Int(phase)
    snap, γ, cap = _moving_diph_snapshot(model, xvec; phase=phase_i)

    if field === :omega || field === :bulk
        uγ = (show_interfacial && cap !== nothing) ? _interface_values(cap, γ) : nothing
        return PenguinPlots.plot_scalar_field(snap; uγ=uγ, kwargs...)
    elseif field === :gamma || field === :interface || field === :uγ
        isnap = _interface_snapshot(cap, γ; label="uγ$phase_i", phase=phase_i)
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    end

    throw(ArgumentError("unsupported diffusion field `$field`"))
end

PenguinPlots.plot_bulk(model::PenguinDiffusion.DiffusionModelMono, x; kwargs...) =
    plot_solution(model, x; field=:bulk, kwargs...)

PenguinPlots.plot_bulk(model::PenguinDiffusion.MovingDiffusionModelMono, x; kwargs...) =
    plot_solution(model, x; field=:bulk, kwargs...)

PenguinPlots.plot_bulk(model::PenguinDiffusion.DiffusionModelDiph, x; phase=1, kwargs...) =
    plot_solution(model, x; phase=phase, field=:bulk, kwargs...)

PenguinPlots.plot_bulk(model::PenguinDiffusion.MovingDiffusionModelDiph, x; phase=1, kwargs...) =
    plot_solution(model, x; phase=phase, field=:bulk, kwargs...)

PenguinPlots.plot_interface(model::PenguinDiffusion.DiffusionModelMono, x; kwargs...) =
    plot_solution(model, x; field=:interface, kwargs...)

PenguinPlots.plot_interface(model::PenguinDiffusion.MovingDiffusionModelMono, x; kwargs...) =
    plot_solution(model, x; field=:interface, kwargs...)

PenguinPlots.plot_interface(model::PenguinDiffusion.DiffusionModelDiph, x; phase=1, kwargs...) =
    plot_solution(model, x; phase=phase, field=:interface, kwargs...)

PenguinPlots.plot_interface(model::PenguinDiffusion.MovingDiffusionModelDiph, x; phase=1, kwargs...) =
    plot_solution(model, x; phase=phase, field=:interface, kwargs...)

end # module
