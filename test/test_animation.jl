@testset "Animation" begin
    x = collect(range(0.0, 1.0; length=24))
    y = collect(range(0.0, 1.0; length=18))
    mask = trues(length(x), length(y))

    snaps = ScalarSnapshot[]
    for k in 1:3
        vals = [sin(2pi * xi + 0.3 * k) * cos(2pi * yj) for xi in x, yj in y]
        phi = [yi - 0.5 for xi in x, yi in y]
        push!(snaps, ScalarSnapshot(x, y, vals, mask; title="frame $k", time=0.1 * (k - 1), metadata=Dict(:phi => phi)))
    end

    mktempdir() do tmp
        path = joinpath(tmp, "scalar_anim.gif")
        out = animate_solution(snaps, path; fps=4)
        @test out.path == path
        @test isfile(path)
        @test filesize(path) > 0
    end
end

@testset "Stefan extension smoke" begin
    required = ("PenguinStefan", "CartesianGrids", "PenguinBCs")
    if any(Base.find_package(pkg) === nothing for pkg in required)
        @info "Stefan plotting smoke skipped: missing one of $required"
        @test true
    else
        @eval begin
            using CartesianGrids: CartesianGrid
            using PenguinBCs
            using PenguinStefan
        end

        grid = CartesianGrid((0.0,), (1.0,), (41,))
        rep = LevelSetRep(grid, x -> x - 0.5)
        bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
        params = StefanParams(0.5, 2.0; kappa1=1.0, source1=0.0)
        opts = StefanOptions(; scheme=:BE, reinit=false, extend_iters=6, interface_band=2.0)
        prob = StefanMonoProblem(grid, bc, params, rep, opts)

        solver = build_solver(
            prob;
            t0=0.0,
            uω0=(x, t) -> 1.0 - 0.8 * x,
            uγ0=0.5,
            speed0=0.0,
        )

        result = solve!(solver, (0.0, 0.03); dt=0.01, save_history=true)

        fig_sol = plot_solution(solver; field=:uω)
        @test fig_sol isa Figure

        fig_hist = plot_history(result; field=:speedmax)
        @test fig_hist isa Figure

        mktempdir() do tmp
            path = joinpath(tmp, "stefan_anim.gif")
            out = animate_solution(result, path; field=:uω, fps=4)
            @test out.path == path
            @test isfile(path)
            @test filesize(path) > 0
        end
    end
end
