function _activate_optional_backend!(pkg::Symbol, activate_name::Symbol; kwargs...)
    if Base.find_package(String(pkg)) === nothing
        throw(ArgumentError("backend `$(lowercase(String(pkg))[1:(end - 5)])` requires package $(pkg). Add it to your environment and retry."))
    end
    @eval import $pkg
    mod = getproperty(@__MODULE__, pkg)
    getproperty(mod, activate_name)(; kwargs...)
    return pkg
end

function activate_penguin_backend!(; backend::Symbol=:cairo, kwargs...)
    if backend === :cairo
        CairoMakie.activate!(; kwargs...)
        return :cairo
    elseif backend === :gl
        _activate_optional_backend!(:GLMakie, :activate!; kwargs...)
        return :gl
    elseif backend === :wgl
        _activate_optional_backend!(:WGLMakie, :activate!; kwargs...)
        return :wgl
    end
    throw(ArgumentError("unsupported backend `$backend`; expected :cairo, :gl, or :wgl"))
end

function save_figure(path, fig; kwargs...)
    Makie.save(path, fig; kwargs...)
    return path
end
