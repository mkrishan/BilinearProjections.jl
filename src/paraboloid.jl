residual(C::CanonicalHyperbolicParaboloid,
         point::Tuple{<:AbstractArray,<:AbstractArray,<:Real}) = begin
    u, v, z = point
    _check_pair(u, v)
    _sumsq(u) - _sumsq(v) - 2 * C.alpha * z
end

residual(C::HyperbolicParaboloid,
         point::Tuple{<:AbstractArray,<:AbstractArray,<:Real}) = begin
    x, y, z = point
    _check_pair(x, y)
    real(dot(x, y)) - C.alpha * z
end

_residual_scale(C::CanonicalHyperbolicParaboloid, point) = begin
    u, v, z = point
    _sumsq(u) + _sumsq(v) + 2 * abs(C.alpha * z)
end

_residual_scale(C::HyperbolicParaboloid, point) = begin
    x, y, z = point
    norm(x) * norm(y) + abs(C.alpha * z)
end

function _canonical_paraboloid_info(C::CanonicalHyperbolicParaboloid,
                                    input, point, multiplicity::Symbol,
                                    branch::Symbol, iterations::Int,
                                    converged::Bool)
    u0, v0, z0 = input
    u, v, z = point
    T = promote_type(eltype(u), eltype(v), typeof(z))
    du = u - u0
    dv = v - v0
    dz = z - z0
    denominator = T(_sumsq(u) + _sumsq(v) + C.alpha^2)
    numerator = T(real(dot(du, u) - dot(dv, v)) - C.alpha * C.beta^2 * dz)
    lambda = iszero(denominator) ? zero(T) : -numerator / denominator
    stationarity = sqrt(T(_sumsq(du + lambda * u) +
                              _sumsq(dv - lambda * v) +
                              abs2(C.beta^2 * dz - lambda * C.alpha)))
    return ProjectionInfo(lambda,
                          T(_distance2_triple(input, point, C.beta)),
                          T(abs(residual(C, point))),
                          stationarity,
                          branch,
                          multiplicity,
                          iterations,
                          converged)
end

function _bilinear_paraboloid_info(C::HyperbolicParaboloid,
                                   input, point, multiplicity::Symbol,
                                   branch::Symbol, iterations::Int,
                                   converged::Bool)
    x0, y0, z0 = input
    x, y, z = point
    T = promote_type(eltype(x), eltype(y), typeof(z))
    dx = x - x0
    dy = y - y0
    dz = z - z0
    denominator = T(_sumsq(x) + _sumsq(y) + C.alpha^2)
    numerator = T(real(dot(dx, y) + dot(dy, x)) - C.alpha * C.beta^2 * dz)
    lambda = iszero(denominator) ? zero(T) : -numerator / denominator
    stationarity = sqrt(T(_sumsq(dx + lambda * y) +
                              _sumsq(dy + lambda * x) +
                              abs2(C.beta^2 * dz - lambda * C.alpha)))
    return ProjectionInfo(lambda,
                          T(_distance2_triple(input, point, C.beta)),
                          T(abs(residual(C, point))),
                          stationarity,
                          branch,
                          multiplicity,
                          iterations,
                          converged)
end

function projection_set(C::CanonicalHyperbolicParaboloid,
                        point::Tuple{<:AbstractArray,<:AbstractArray,<:Real};
                        atol=nothing, rtol=nothing, maxiters::Int=200)
    u0, v0, z0 = point
    _check_pair(u0, v0)
    T = _float_type(u0, v0, z0, C.alpha, C.beta)
    uT = _asarray(u0, T)
    vT = _asarray(v0, T)
    zT = T(z0)
    CT = CanonicalHyperbolicParaboloid(T(C.alpha), T(C.beta))
    raw = _canonical_paraboloid_raw(uT, vT, zT, CT.alpha, CT.beta, T;
        atol, rtol, maxiters)
    input = (uT, vT, zT)
    return _public_projection(raw) do projected, multiplicity, branch, iterations, converged
        _canonical_paraboloid_info(CT, input, projected, multiplicity,
                                   branch, iterations, converged)
    end
end

function projection_set(C::HyperbolicParaboloid,
                        point::Tuple{<:AbstractArray,<:AbstractArray,<:Real};
                        atol=nothing, rtol=nothing, maxiters::Int=200)
    x0, y0, z0 = point
    _check_pair(x0, y0)
    T = _float_type(x0, y0, z0, C.alpha, C.beta)
    xT = _asarray(x0, T)
    yT = _asarray(y0, T)
    zT = T(z0)
    u0, v0 = _canonical_coordinates(xT, yT, T)
    CTcanonical = CanonicalHyperbolicParaboloid(T(C.alpha), T(C.beta))
    raw = _canonical_paraboloid_raw(u0, v0, zT,
                                    CTcanonical.alpha, CTcanonical.beta, T;
                                    atol, rtol, maxiters)
    mapped = _map_raw(raw, canonical_point -> begin
        x, y = _bilinear_coordinates(canonical_point[1], canonical_point[2])
        (x, y, canonical_point[3])
    end)
    input = (xT, yT, zT)
    CT = HyperbolicParaboloid(T(C.alpha), T(C.beta))
    return _public_projection(mapped) do projected, multiplicity, branch, iterations, converged
        _bilinear_paraboloid_info(CT, input, projected, multiplicity,
                                  branch, iterations, converged)
    end
end
