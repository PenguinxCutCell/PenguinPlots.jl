struct ScalarSnapshot
    x::AbstractVector
    y::Union{Nothing,AbstractVector}
    values::AbstractArray
    mask::AbstractArray{Bool}
    title::String
    colorlabel::String
    time
    phase::Union{Nothing,Int}
    interface_xy::Union{Nothing,Tuple}
    metadata::Dict{Symbol,Any}
end

function ScalarSnapshot(
    x::AbstractVector,
    y::Union{Nothing,AbstractVector},
    values::AbstractArray,
    mask::AbstractArray{Bool};
    title::AbstractString="",
    colorlabel::AbstractString="",
    time=nothing,
    phase::Union{Nothing,Int}=nothing,
    interface_xy::Union{Nothing,Tuple}=nothing,
    metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}(),
)
    return ScalarSnapshot(
        x,
        y,
        values,
        mask,
        String(title),
        String(colorlabel),
        time,
        phase,
        interface_xy,
        Dict{Symbol,Any}(metadata),
    )
end

struct InterfaceSnapshot
    xγ::AbstractVector
    yγ::Union{Nothing,AbstractVector}
    valuesγ::AbstractVector
    label::String
    time
    phase::Union{Nothing,Int}
    metadata::Dict{Symbol,Any}
end

function InterfaceSnapshot(
    xγ::AbstractVector,
    yγ::Union{Nothing,AbstractVector},
    valuesγ::AbstractVector;
    label::AbstractString="uγ",
    time=nothing,
    phase::Union{Nothing,Int}=nothing,
    metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}(),
)
    return InterfaceSnapshot(
        xγ,
        yγ,
        valuesγ,
        String(label),
        time,
        phase,
        Dict{Symbol,Any}(metadata),
    )
end

struct VectorSnapshot
    x::AbstractVector
    y::AbstractVector
    u::AbstractArray
    v::AbstractArray
    p::Union{Nothing,AbstractArray}
    mask::AbstractArray{Bool}
    interface_xy::Union{Nothing,Tuple}
    title::String
    time
    phase::Union{Nothing,Int}
    metadata::Dict{Symbol,Any}
end

function VectorSnapshot(
    x::AbstractVector,
    y::AbstractVector,
    u::AbstractArray,
    v::AbstractArray;
    p::Union{Nothing,AbstractArray}=nothing,
    mask::AbstractArray{Bool}=trues(size(u)),
    interface_xy::Union{Nothing,Tuple}=nothing,
    title::AbstractString="",
    time=nothing,
    phase::Union{Nothing,Int}=nothing,
    metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}(),
)
    return VectorSnapshot(
        x,
        y,
        u,
        v,
        p,
        mask,
        interface_xy,
        String(title),
        time,
        phase,
        Dict{Symbol,Any}(metadata),
    )
end

struct ResidualHistorySnapshot
    iterations::AbstractVector
    residuals::AbstractVector
    block_residuals::Dict{Symbol,AbstractVector}
    title::String
    metadata::Dict{Symbol,Any}
end

function ResidualHistorySnapshot(
    iterations::AbstractVector,
    residuals::AbstractVector;
    block_residuals::AbstractDict{Symbol,<:AbstractVector}=Dict{Symbol,Vector{Float64}}(),
    title::AbstractString="Residual history",
    metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}(),
)
    br = Dict{Symbol,AbstractVector}()
    for (k, v) in pairs(block_residuals)
        br[k] = v
    end
    return ResidualHistorySnapshot(
        iterations,
        residuals,
        br,
        String(title),
        Dict{Symbol,Any}(metadata),
    )
end
