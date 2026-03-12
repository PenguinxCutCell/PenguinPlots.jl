function animate_solution(
    history::AbstractVector{<:ScalarSnapshot},
    path;
    fps::Integer=20,
    field::Symbol=:u,
    backend::Symbol=:cairo,
    kwargs...,
)
    isempty(history) && throw(ArgumentError("history cannot be empty"))
    activate_penguin_backend!(; backend=backend)

    fig = Figure(size=(900, 500))
    ax = Axis(fig[1, 1])

    record(fig, path, eachindex(history); framerate=fps) do i
        _clear_axis!(ax)
        snap = history[i]
        phi = get(snap.metadata, :phi, nothing)
        plot_scalar_field(snap; fig=fig, ax=ax, show_colorbar=false, phi=phi, kwargs...)
    end

    return (path=path, fig=fig, field=field)
end

function animate_solution(
    history::AbstractVector{<:VectorSnapshot},
    path;
    fps::Integer=20,
    field::Symbol=:velocity,
    backend::Symbol=:cairo,
    kwargs...,
)
    isempty(history) && throw(ArgumentError("history cannot be empty"))
    activate_penguin_backend!(; backend=backend)

    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])

    record(fig, path, eachindex(history); framerate=fps) do i
        _clear_axis!(ax)
        plot_vector_field(history[i]; fig=fig, ax=ax, show_colorbar=false, kwargs...)
    end

    return (path=path, fig=fig, field=field)
end
