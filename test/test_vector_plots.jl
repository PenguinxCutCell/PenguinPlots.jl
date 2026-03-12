@testset "Vector plots" begin
    x = collect(range(0.0, 1.0; length=30))
    y = collect(range(0.0, 1.0; length=20))

    u = [sin(pi * xi) for xi in x, _ in y]
    v = [cos(pi * yj) for _ in x, yj in y]
    p = [xi - yj for xi in x, yj in y]
    mask = trues(size(u))
    mask[2:4, 2:4] .= false

    snap = VectorSnapshot(
        x,
        y,
        u,
        v;
        p=p,
        mask=mask,
        interface_xy=([0.2, 0.5, 0.8], [0.2, 0.5, 0.8]),
        title="Synthetic flow",
    )

    fig = plot_vector_field(snap; stride=3)
    @test fig isa Figure

    figp = plot_pressure_field(snap)
    @test figp isa Figure

    mktempdir() do tmp
        path = joinpath(tmp, "vector.png")
        save_figure(path, fig)
        @test isfile(path)
        @test filesize(path) > 0
    end
end

@testset "Stokes extension smoke" begin
    required = ("PenguinStokes", "CartesianGrids", "PenguinBCs")
    if any(Base.find_package(pkg) === nothing for pkg in required)
        @info "Stokes plotting smoke skipped: missing one of $required"
        @test true
    else
        @eval begin
            using CartesianGrids: CartesianGrid
            using PenguinBCs
            using PenguinStokes
        end

        full_body(args...) = -1.0
        bc = BorderConditions(
            ;
            left=Dirichlet(0.0),
            right=Dirichlet(0.0),
            bottom=Dirichlet(0.0),
            top=Dirichlet(0.0),
        )

        grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (9, 9))
        model = StokesModelMono(
            grid,
            full_body,
            1.0,
            1.0;
            bc_u=(bc, bc),
            bc_cut=Dirichlet(0.0),
            force=(0.0, 0.0),
        )

        sys = solve_steady!(model)

        fig_p = plot_solution(model, sys.x; field=:pressure)
        @test fig_p isa Figure

        fig_u = plot_solution(model, sys.x; field=:velocity, stride=2)
        @test fig_u isa Figure

        fig_all = plot_solution(model, sys.x; field=:pressure_velocity, stride=2)
        @test fig_all isa Figure
    end
end
