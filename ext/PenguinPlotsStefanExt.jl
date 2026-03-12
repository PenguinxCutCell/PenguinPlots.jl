module PenguinPlotsStefanExt

using PenguinPlots
using PenguinStefan

import PenguinPlots: animate_solution, plot_history, plot_solution

function _rep_grid(rep)
    if hasproperty(rep, :grid)
        return getproperty(rep, :grid)
    elseif hasproperty(rep, :space_grid)
        return getproperty(rep, :space_grid)
    end
    throw(ArgumentError("cannot infer grid from interface representation $(typeof(rep))"))
end

function _grid_axes_for_scalar(grid)
    axes = PenguinPlots.cell_center_axes_from_grid(grid)
    if length(axes) == 1
        return axes[1], nothing
    end
    return axes[1], axes[2]
end

function _reshape_state_field(field_values, dims::NTuple)
    if field_values isa AbstractArray && size(field_values) == dims
        return field_values
    end
    length(field_values) == prod(dims) || throw(DimensionMismatch("field size does not match grid dimensions $dims"))
    return reshape(collect(field_values), dims)
end

function _mono_scalar_snapshot(
    state::PenguinStefan.StefanMonoState,
    rep;
    field::Symbol=:uω,
    title::AbstractString="Stefan mono",
)
    grid = _rep_grid(rep)
    dims = Tuple(getproperty(grid, :n))
    x, y = _grid_axes_for_scalar(grid)
    phi = PenguinStefan.phi_values(rep)

    if field === :uω || field === :temperature
        vals = _reshape_state_field(state.uω, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="uω",
            time=state.t,
            metadata=Dict{Symbol,Any}(:phi => phi),
        )
    elseif field === :speed
        vals = _reshape_state_field(state.speed_full, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="speed",
            time=state.t,
            metadata=Dict{Symbol,Any}(:phi => phi),
        )
    elseif field === :phi
        vals = _reshape_state_field(phi, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="phi",
            time=state.t,
            metadata=Dict{Symbol,Any}(),
        )
    end

    throw(ArgumentError("unsupported Stefan mono field `$field`"))
end

function _diph_scalar_snapshot(
    state::PenguinStefan.StefanDiphState,
    rep;
    field::Symbol=:uω1,
    title::AbstractString="Stefan diph",
)
    grid = _rep_grid(rep)
    dims = Tuple(getproperty(grid, :n))
    x, y = _grid_axes_for_scalar(grid)
    phi = PenguinStefan.phi_values(rep)

    if field === :uω1 || field === :temperature1
        vals = _reshape_state_field(state.uω1, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="uω1",
            time=state.t,
            phase=1,
            metadata=Dict{Symbol,Any}(:phi => phi),
        )
    elseif field === :uω2 || field === :temperature2
        vals = _reshape_state_field(state.uω2, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="uω2",
            time=state.t,
            phase=2,
            metadata=Dict{Symbol,Any}(:phi => phi),
        )
    elseif field === :speed
        vals = _reshape_state_field(state.speed_full, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="speed",
            time=state.t,
            metadata=Dict{Symbol,Any}(:phi => phi),
        )
    elseif field === :phi
        vals = _reshape_state_field(phi, dims)
        return PenguinPlots.ScalarSnapshot(
            x,
            y,
            vals,
            trues(size(vals));
            title=title,
            colorlabel="phi",
            time=state.t,
        )
    end

    throw(ArgumentError("unsupported Stefan diph field `$field`"))
end

function _interface_snapshot_from_cap(cap, uγ; label::AbstractString="uγ", time=nothing, phase=nothing)
    mask = vec(PenguinPlots.interface_mask_from_capacity(cap))
    xγ = Float64[]
    yγ = Float64[]
    vals = Float64[]

    Cγ = getproperty(cap, :C_γ)
    uγv = vec(collect(uγ))
    n = min(length(mask), min(length(Cγ), length(uγv)))
    for i in 1:n
        mask[i] || continue
        push!(xγ, Float64(Cγ[i][1]))
        if length(Cγ[i]) >= 2
            push!(yγ, Float64(Cγ[i][2]))
        end
        push!(vals, Float64(uγv[i]))
    end

    if isempty(xγ)
        return PenguinPlots.InterfaceSnapshot(Float64[], nothing, Float64[]; label=label, time=time, phase=phase)
    end

    if isempty(yγ)
        return PenguinPlots.InterfaceSnapshot(xγ, nothing, vals; label=label, time=time, phase=phase)
    end
    return PenguinPlots.InterfaceSnapshot(xγ, yγ, vals; label=label, time=time, phase=phase)
end

function plot_solution(
    solver::PenguinStefan.StefanMonoSolver;
    field::Symbol=:uω,
    show_phi::Bool=true,
    kwargs...,
)
    if field === :uγ
        cap = something(solver.cache.model.cap_slab, nothing)
        cap === nothing && throw(ArgumentError("interface plotting requires slab capacity; step or solve the solver first"))
        isnap = _interface_snapshot_from_cap(cap, solver.state.uγ; label="uγ", time=solver.state.t)
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    end

    snap = _mono_scalar_snapshot(solver.state, solver.problem.interface_rep; field=field)
    phi_kw = show_phi ? get(snap.metadata, :phi, nothing) : nothing
    return PenguinPlots.plot_scalar_field(snap; phi=phi_kw, kwargs...)
end

function plot_solution(
    solver::PenguinStefan.StefanDiphSolver;
    field::Symbol=:uω1,
    show_phi::Bool=true,
    kwargs...,
)
    if field === :uγ1 || field === :uγ2
        cap = field === :uγ1 ? something(solver.cache.model.cap1_slab, nothing) : something(solver.cache.model.cap2_slab, nothing)
        cap === nothing && throw(ArgumentError("interface plotting requires slab capacity; step or solve the solver first"))
        uγ = field === :uγ1 ? solver.state.uγ1 : solver.state.uγ2
        phase = field === :uγ1 ? 1 : 2
        isnap = _interface_snapshot_from_cap(cap, uγ; label=String(field), time=solver.state.t, phase=phase)
        return PenguinPlots.plot_interface_field(isnap; kwargs...)
    elseif field === :both
        s1 = _diph_scalar_snapshot(solver.state, solver.problem.interface_rep; field=:uω1)
        s2 = _diph_scalar_snapshot(solver.state, solver.problem.interface_rep; field=:uω2)
        if show_phi
            s1.metadata[:phi] = PenguinStefan.phi_values(solver.problem.interface_rep)
            s2.metadata[:phi] = PenguinStefan.phi_values(solver.problem.interface_rep)
        end
        return PenguinPlots.plot_scalar_field([s1, s2]; kwargs...)
    end

    snap = _diph_scalar_snapshot(solver.state, solver.problem.interface_rep; field=field)
    phi_kw = show_phi ? get(snap.metadata, :phi, nothing) : nothing
    return PenguinPlots.plot_scalar_field(snap; phi=phi_kw, kwargs...)
end

function plot_solution(
    state::PenguinStefan.StefanMonoState,
    rep;
    field::Symbol=:uω,
    show_phi::Bool=true,
    kwargs...,
)
    field === :uγ && throw(ArgumentError("`plot_solution(state, rep; field=:uγ)` needs solver cache; pass solver instead"))
    snap = _mono_scalar_snapshot(state, rep; field=field)
    phi_kw = show_phi ? get(snap.metadata, :phi, nothing) : nothing
    return PenguinPlots.plot_scalar_field(snap; phi=phi_kw, kwargs...)
end

function plot_solution(
    state::PenguinStefan.StefanDiphState,
    rep;
    field::Symbol=:uω1,
    show_phi::Bool=true,
    kwargs...,
)
    if field === :uγ1 || field === :uγ2
        throw(ArgumentError("`plot_solution(state, rep; field=$(field))` needs solver cache; pass solver instead"))
    elseif field === :both
        s1 = _diph_scalar_snapshot(state, rep; field=:uω1)
        s2 = _diph_scalar_snapshot(state, rep; field=:uω2)
        return PenguinPlots.plot_scalar_field([s1, s2]; kwargs...)
    end

    snap = _diph_scalar_snapshot(state, rep; field=field)
    phi_kw = show_phi ? get(snap.metadata, :phi, nothing) : nothing
    return PenguinPlots.plot_scalar_field(snap; phi=phi_kw, kwargs...)
end

function _history_times_and_values(result::NamedTuple, field::Symbol)
    if field === :speedmax
        if haskey(result, :history) && !isempty(result.history) && haskey(result.history[1], :speed)
            t = [fr.t for fr in result.history]
            v = [maximum(abs.(vec(fr.speed))) for fr in result.history]
            return t, v, "speed_max"
        elseif haskey(result, :solver)
            logs = result.solver.state.logs
            if haskey(logs, :times) && haskey(logs, :speed_max)
                t = collect(logs[:times])
                v = collect(logs[:speed_max])
                if length(v) == length(t) - 1
                    return t[2:end], v, "speed_max"
                end
                return t, v, "speed_max"
            end
        end
        throw(ArgumentError("result does not expose speed_max history"))
    end

    if haskey(result, :solver)
        logs = result.solver.state.logs
        if haskey(logs, field) && haskey(logs, :times)
            t = collect(logs[:times])
            v = collect(logs[field])
            if length(v) == length(t) - 1
                return t[2:end], v, String(field)
            end
            return t, v, String(field)
        end
    end

    throw(ArgumentError("history field `$field` is not available"))
end

function plot_history(
    result::NamedTuple;
    field::Symbol=:speedmax,
    kwargs...,
)
    t, v, ylabel = _history_times_and_values(result, field)
    ttl = "Stefan history: $(field)"
    return PenguinPlots.plot_history(t, v; ylabel=ylabel, title=ttl, kwargs...)
end

function _frame_to_scalar_snapshot(frame::NamedTuple, grid; field::Symbol=:uω, show_phi::Bool=true)
    dims = Tuple(getproperty(grid, :n))
    xaxes = PenguinPlots.cell_center_axes_from_grid(grid)
    x = xaxes[1]
    y = length(xaxes) >= 2 ? xaxes[2] : nothing

    if !haskey(frame, field)
        throw(ArgumentError("history frame does not contain `$field`"))
    end

    vals_raw = getproperty(frame, field)
    vals = _reshape_state_field(vals_raw, dims)
    md = Dict{Symbol,Any}()
    if show_phi && haskey(frame, :phi)
        md[:phi] = getproperty(frame, :phi)
    end

    return PenguinPlots.ScalarSnapshot(
        x,
        y,
        vals,
        trues(size(vals));
        title="Stefan $(field)",
        colorlabel=String(field),
        time=getproperty(frame, :t),
        metadata=md,
    )
end

function animate_solution(
    result::NamedTuple,
    path;
    field::Symbol=:uω,
    show_phi::Bool=true,
    fps::Integer=20,
    backend::Symbol=:cairo,
    kwargs...,
)
    haskey(result, :history) || throw(ArgumentError("result must contain `history`"))
    haskey(result, :solver) || throw(ArgumentError("result must contain `solver`"))

    grid = result.solver.problem.grid
    snaps = PenguinPlots.ScalarSnapshot[]
    for frame in result.history
        push!(snaps, _frame_to_scalar_snapshot(frame, grid; field=field, show_phi=show_phi))
    end

    return PenguinPlots.animate_solution(snaps, path; fps=fps, field=field, backend=backend, kwargs...)
end

end # module
