@testset "History plots" begin
    t = collect(range(0.0, 1.0; length=20))
    v = exp.(-3 .* t)

    fig = plot_history(t, v; ylabel="energy", title="Decay")
    @test fig isa Figure

    res = ResidualHistorySnapshot(
        1:6,
        [1.0, 0.4, 0.1, 0.03, 0.01, 0.003];
        block_residuals=Dict{Symbol,AbstractVector}(
            :diffusion => [0.8, 0.35, 0.09, 0.02, 0.009, 0.002],
            :stokes => [0.9, 0.42, 0.12, 0.04, 0.014, 0.004],
        ),
        title="Coupling residuals",
    )

    figr = plot_residual_history(res)
    @test figr isa Figure
end

@testset "SolverCore extension smoke" begin
    if Base.find_package("PenguinSolverCore") === nothing
        @test_skip "PenguinSolverCore not available"
    else
        @eval using PenguinSolverCore

        h = PenguinSolverCore.CouplingHistory()
        push!(h.iterations, 1)
        push!(h.residuals, 1.0)
        push!(h.block_residuals, Dict(:a => 1.0, :b => 0.5))
        push!(h.iterations, 2)
        push!(h.residuals, 0.2)
        push!(h.block_residuals, Dict(:a => 0.2, :b => 0.1))
        push!(h.iterations, 3)
        push!(h.residuals, 0.04)
        push!(h.block_residuals, Dict(:a => 0.04, :b => 0.02))

        fig = plot_residual_history(h)
        @test fig isa Figure

        fig_tuple = plot_residual_history((nothing, h))
        @test fig_tuple isa Figure
    end
end
