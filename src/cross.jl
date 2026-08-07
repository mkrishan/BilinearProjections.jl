residual(::CrossConstraint, point::Tuple{<:AbstractArray,<:AbstractArray}) = begin
    x, y = point
    _check_pair(x, y)
    real(dot(x, y))
end

_residual_scale(::CrossConstraint, point) = begin
    x, y = point
    norm(x) * norm(y)
end

function _pair_info(C::Union{CrossConstraint,BilinearLevelSet}, input, point,
                    multiplicity::Symbol, branch::Symbol,
                    iterations::Int, converged::Bool)
    x0, y0 = input
    x, y = point
    T = promote_type(eltype(x), eltype(y))
    dx = x - x0
    dy = y - y0
    denominator = T(_sumsq(x) + _sumsq(y))
    lambda = if iszero(denominator)
        zero(T)
    else
        -T(real(dot(dx, y) + dot(dy, x))) / denominator
    end
    stationarity = sqrt(T(_sumsq(dx + lambda * y) +
                              _sumsq(dy + lambda * x)))
    return ProjectionInfo(lambda,
                          T(_distance2_pair(input, point)),
                          T(abs(residual(C, point))),
                          stationarity,
                          branch,
                          multiplicity,
                          iterations,
                          converged)
end

function projection_set(C::CrossConstraint,
                        point::Tuple{<:AbstractArray,<:AbstractArray}; kwargs...)
    x0, y0 = point
    _check_pair(x0, y0)
    T = _float_type(x0, y0)
    xT = _asarray(x0, T)
    yT = _asarray(y0, T)
    u0, v0 = _canonical_coordinates(xT, yT, T)
    raw = _canonical_cross_raw(u0, v0, T)
    mapped = _map_raw(raw, canonical_point ->
        _bilinear_coordinates(canonical_point[1], canonical_point[2]))
    input = (xT, yT)
    return _public_projection(mapped) do projected, multiplicity, branch, iterations, converged
        _pair_info(C, input, projected, multiplicity, branch, iterations, converged)
    end
end
