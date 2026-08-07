residual(C::BilinearLevelSet, point::Tuple{<:AbstractArray,<:AbstractArray}) = begin
    x, y = point
    _check_pair(x, y)
    real(dot(x, y)) - C.level
end

_residual_scale(C::BilinearLevelSet, point) = begin
    x, y = point
    norm(x) * norm(y) + abs(C.level)
end

function projection_set(C::BilinearLevelSet,
                        point::Tuple{<:AbstractArray,<:AbstractArray};
                        atol=nothing, rtol=nothing, maxiters::Int=200)
    x0, y0 = point
    _check_pair(x0, y0)
    T = _float_type(x0, y0, C.level)
    gamma = T(C.level)
    xT = _asarray(x0, T)
    yT = _asarray(y0, T)

    if iszero(gamma)
        P = projection_set(CrossConstraint(), (xT, yT))
        return P
    end

    u0, v0 = _canonical_coordinates(xT, yT, T)
    raw = _canonical_level_raw(u0, v0, gamma, T;
        atol, rtol, maxiters)
    mapped = _map_raw(raw, canonical_point ->
        _bilinear_coordinates(canonical_point[1], canonical_point[2]))
    input = (xT, yT)
    CT = BilinearLevelSet(gamma)
    return _public_projection(mapped) do projected, multiplicity, branch, iterations, converged
        _pair_info(CT, input, projected, multiplicity, branch, iterations, converged)
    end
end
