@testset "Scalar plots" begin
    x = collect(range(0.0, 1.0; length=32))
    y = collect(range(0.0, 1.0; length=24))
    vals = [sin(2pi * xi) * cos(2pi * yj) for xi in x, yj in y]
    mask = trues(size(vals))
    mask[1:3, :] .= false

    snap = ScalarSnapshot(
        x,
        y,
        vals,
        mask;
        title="Synthetic scalar",
        colorlabel="u",
        interface_xy=([0.25, 0.5, 0.75], [0.2, 0.5, 0.8]),
    )

    fig = plot_scalar_field(snap)
    @test fig isa Figure

    fig_phase = plot_scalar_field([
        ScalarSnapshot(x, y, vals, mask; title="phase 1", colorlabel="u", phase=1),
        ScalarSnapshot(x, y, 0.5 .* vals, mask; title="phase 2", colorlabel="u", phase=2),
    ]; shared_clims=true)
    @test fig_phase isa Figure

    itf = InterfaceSnapshot([0.2, 0.4, 0.6], [0.3, 0.5, 0.7], [1.0, -0.5, 0.2]; label="uγ")
    fig_i = plot_interface_field(itf; style=:scatter)
    @test fig_i isa Figure

    itf_empty = InterfaceSnapshot(Float64[], nothing, Float64[]; label="uγ")
    fig_ie = plot_interface_field(itf_empty)
    @test fig_ie isa Figure

    mktempdir() do tmp
        path = joinpath(tmp, "scalar.png")
        save_figure(path, fig)
        @test isfile(path)
        @test filesize(path) > 0
    end
end

@testset "Diffusion extension smoke" begin
    required = ("PenguinDiffusion", "CartesianGeometry", "CartesianOperators", "PenguinBCs")
    if any(Base.find_package(pkg) === nothing for pkg in required)
        @info "Diffusion plotting smoke skipped: missing one of $required"
        @test true
    else
        @eval begin
            using CartesianGeometry: geometric_moments, nan
            using CartesianOperators
            using PenguinBCs
            using PenguinDiffusion
        end

        grid = (range(0.0, 1.0; length=21), range(0.0, 1.0; length=17))
        moms = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)
        cap = assembled_capacity(moms; bc=0.0)
        bc = BorderConditions(
            ;
            left=Dirichlet(0.0),
            right=Dirichlet(0.0),
            bottom=Dirichlet(0.0),
            top=Dirichlet(0.0),
        )
        ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))

        mono = DiffusionModelMono(cap, ops, 1.0; source=0.0, bc_border=bc)
        mono_sys = solve_steady!(mono)
        fig_mono = plot_solution(mono, mono_sys.x; field=:omega)
        @test fig_mono isa Figure

        diph = DiffusionModelDiph(cap, ops, 1.0, 0.6; source=(0.0, 0.0), bc_border=bc)
        diph_sys = solve_steady!(diph)
        fig_diph = plot_solution(diph, diph_sys.x; field=:omega, phase=:both)
        @test fig_diph isa Figure

        fig_bulk = plot_bulk(mono, mono_sys.x)
        @test fig_bulk isa Figure
    end
end
