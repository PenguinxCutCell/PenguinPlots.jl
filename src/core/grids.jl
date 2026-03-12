function capacity_dims(cap)
    return Tuple(getproperty(cap, :nnodes))
end

function capacity_axes(cap)
    xyz = getproperty(cap, :xyz)
    return ntuple(d -> collect(xyz[d]), length(xyz))
end

function active_mask_from_capacity(cap)
    V = getproperty(getproperty(cap, :buf), :V)
    out = Vector{Bool}(undef, length(V))
    @inbounds for i in eachindex(V)
        vi = V[i]
        out[i] = isfinite(vi) && vi > 0
    end
    return reshape(out, capacity_dims(cap))
end

function interface_mask_from_capacity(cap)
    Γ = getproperty(getproperty(cap, :buf), :Γ)
    out = Vector{Bool}(undef, length(Γ))
    @inbounds for i in eachindex(Γ)
        gi = Γ[i]
        out[i] = isfinite(gi) && gi > 0
    end
    return reshape(out, capacity_dims(cap))
end

function nanfilled(values::AbstractArray, mask::AbstractArray{Bool})
    return _masked_array(values, mask)
end

function reshape_on_capacity(values::AbstractArray, cap)
    return _coerce_values_to_shape(values, capacity_dims(cap))
end

function interface_coordinates_from_capacity(cap)
    Γ = getproperty(getproperty(cap, :buf), :Γ)
    Cγ = getproperty(cap, :C_γ)
    N = length(first(Cγ))
    xγ = Float64[]
    yγ = Float64[]

    @inbounds for i in eachindex(Γ)
        gi = Γ[i]
        if !(isfinite(gi) && gi > 0)
            continue
        end
        xi = Cγ[i]
        if N >= 1 && isfinite(xi[1])
            push!(xγ, Float64(xi[1]))
        end
        if N >= 2 && isfinite(xi[2])
            push!(yγ, Float64(xi[2]))
        end
    end

    isempty(xγ) && return nothing
    if N == 1
        return (xγ, nothing)
    end
    return (xγ, yγ)
end

function cell_center_axes_from_grid(grid)
    n = Tuple(getproperty(grid, :n))
    lc = Tuple(getproperty(grid, :lc))
    hc = Tuple(getproperty(grid, :hc))
    N = length(n)
    Δ = ntuple(d -> (hc[d] - lc[d]) / n[d], N)
    return ntuple(d -> [lc[d] + (i - 0.5) * Δ[d] for i in 1:n[d]], N)
end

function full_mask_from_dims(dims::NTuple)
    return trues(dims)
end
