function _canonical_cross_raw(u0::AbstractArray, v0::AbstractArray, ::Type{T}) where {T<:AbstractFloat}
    a = T(norm(u0))
    b = T(norm(v0))
    r = (a + b) / T(2)

    if iszero(a) && iszero(b)
        return _raw_singleton((copy(u0), copy(v0)), :already_feasible)
    elseif a > zero(T) && b > zero(T)
        u = (r / a) * u0
        v = (r / b) * v0
        return _raw_singleton((u, v), :generic_closed_form)
    elseif iszero(a)
        v = v0 / T(2)
        generator = direction -> begin
            d = _normalize_direction(direction, u0, T)
            (r * d, copy(v))
        end
        default = generator(_default_direction(v0, T))
        return _raw_sphere(generator, default, :zero_first_canonical_component;
            description="all unit directions for the first canonical component")
    else
        u = u0 / T(2)
        generator = direction -> begin
            d = _normalize_direction(direction, v0, T)
            (copy(u), r * d)
        end
        default = generator(_default_direction(u0, T))
        return _raw_sphere(generator, default, :zero_second_canonical_component;
            description="all unit directions for the second canonical component")
    end
end

function _canonical_positive_level_raw(u0::AbstractArray, v0::AbstractArray,
                                       gamma::T, ::Type{T};
                                       atol=nothing, rtol=nothing,
                                       maxiters::Int=200) where {T<:AbstractFloat}
    gamma > zero(T) || throw(ArgumentError("gamma must be positive"))
    a2 = T(_sumsq(u0))
    b2 = T(_sumsq(v0))
    a = sqrt(a2)
    b = sqrt(b2)

    if a > zero(T) && b > zero(T)
        f = lambda -> a2 / (one(T) + lambda)^2 -
                      b2 / (one(T) - lambda)^2 - T(2) * gamma
        lambda, iters, _, converged = _bisect_decreasing(f, T;
            atol, rtol, maxiters)
        u = u0 / (one(T) + lambda)
        v = v0 / (one(T) - lambda)
        return _raw_singleton((u, v), :generic_root;
            iterations=iters, converged=converged)
    elseif iszero(a)
        v = v0 / T(2)
        r2 = T(2) * gamma + b2 / T(4)
        r = _radius(r2, T(2) * gamma + b2)
        generator = direction -> begin
            d = _normalize_direction(direction, u0, T)
            (r * d, copy(v))
        end
        default = generator(_default_direction(v0, T))
        return _raw_sphere(generator, default, :zero_first_canonical_component;
            description="all unit directions for the first canonical component")
    else
        threshold = T(2) * sqrt(T(2) * gamma)
        if a < threshold
            u = sqrt(T(2) * gamma) * u0 / a
            v = _zero_like(v0, T)
            return _raw_singleton((u, v), :radial_singleton)
        end
        u = u0 / T(2)
        r2 = a2 / T(4) - T(2) * gamma
        r = _radius(r2, a2 + T(2) * gamma)
        if iszero(r)
            return _raw_singleton((u, _zero_like(v0, T)), :threshold_singleton)
        end
        generator = direction -> begin
            d = _normalize_direction(direction, v0, T)
            (copy(u), r * d)
        end
        default = generator(_default_direction(u0, T))
        return _raw_sphere(generator, default, :zero_second_canonical_component;
            description="all unit directions for the second canonical component")
    end
end

function _canonical_level_raw(u0::AbstractArray, v0::AbstractArray,
                              gamma::T, ::Type{T}; kwargs...) where {T<:AbstractFloat}
    if iszero(gamma)
        return _canonical_cross_raw(u0, v0, T)
    elseif gamma > zero(T)
        return _canonical_positive_level_raw(u0, v0, gamma, T; kwargs...)
    end

    raw = _canonical_positive_level_raw(v0, u0, -gamma, T; kwargs...)
    mapped = _map_raw(raw, point -> (point[2], point[1]))
    branch = Symbol("negative_level_", String(mapped.branch))
    return merge(mapped, (branch=branch,))
end

function _canonical_paraboloid_raw(u0::AbstractArray, v0::AbstractArray, z0::T,
                                   alpha::T, beta::T, ::Type{T};
                                   atol=nothing, rtol=nothing,
                                   maxiters::Int=200) where {T<:AbstractFloat}
    a2 = T(_sumsq(u0))
    b2 = T(_sumsq(v0))
    a = sqrt(a2)
    b = sqrt(b2)
    beta2 = beta^2
    alpha2_over_beta2 = alpha^2 / beta2

    if iszero(a) && iszero(b)
        az = alpha * z0
        if az > alpha2_over_beta2
            z = z0 - alpha / beta2
            r = _radius(T(2) * alpha * z, abs(T(2) * alpha * z))
            generator = direction -> begin
                d = _normalize_direction(direction, u0, T)
                (r * d, _zero_like(v0, T), z)
            end
            default = generator(_first_unit(u0, T))
            return _raw_sphere(generator, default, :origin_positive_branch;
                description="all unit directions for the first canonical component")
        elseif abs(az) <= alpha2_over_beta2
            return _raw_singleton((_zero_like(u0, T), _zero_like(v0, T), zero(T)),
                                  :origin_zero_branch)
        else
            z = z0 + alpha / beta2
            r = _radius(-T(2) * alpha * z, abs(T(2) * alpha * z))
            generator = direction -> begin
                d = _normalize_direction(direction, v0, T)
                (_zero_like(u0, T), r * d, z)
            end
            default = generator(_first_unit(v0, T))
            return _raw_sphere(generator, default, :origin_negative_branch;
                description="all unit directions for the second canonical component")
        end
    elseif a > zero(T) && b > zero(T)
        f = lambda -> a2 / (one(T) + lambda)^2 -
                      b2 / (one(T) - lambda)^2 -
                      T(2) * alpha * (z0 + lambda * alpha / beta2)
        lambda, iters, _, converged = _bisect_decreasing(f, T;
            atol, rtol, maxiters)
        u = u0 / (one(T) + lambda)
        v = v0 / (one(T) - lambda)
        z = z0 + lambda * alpha / beta2
        return _raw_singleton((u, v, z), :generic_root;
            iterations=iters, converged=converged)
    elseif iszero(a)
        threshold_value = alpha * (z0 - alpha / beta2)
        if threshold_value < -b2 / T(8)
            f = lambda -> b2 / (one(T) - lambda)^2 +
                          T(2) * lambda * alpha^2 / beta2 +
                          T(2) * alpha * z0
            lambda, iters, _, converged = _bisect_decreasing(lambda -> -f(lambda), T;
                atol, rtol, maxiters)
            # `f` is strictly increasing, hence `-f` is strictly decreasing.
            v = v0 / (one(T) - lambda)
            z = z0 + lambda * alpha / beta2
            return _raw_singleton((_zero_like(u0, T), v, z), :zero_first_root;
                iterations=iters, converged=converged)
        end
        v = v0 / T(2)
        z = z0 - alpha / beta2
        r2 = T(2) * alpha * z + b2 / T(4)
        r = _radius(r2, abs(T(2) * alpha * z) + b2)
        if iszero(r)
            return _raw_singleton((_zero_like(u0, T), v, z), :zero_first_threshold)
        end
        generator = direction -> begin
            d = _normalize_direction(direction, u0, T)
            (r * d, copy(v), z)
        end
        default = generator(_default_direction(v0, T))
        return _raw_sphere(generator, default, :zero_first_family;
            description="all unit directions for the first canonical component")
    else
        threshold_value = alpha * (z0 + alpha / beta2)
        if threshold_value > a2 / T(8)
            f = lambda -> a2 / (one(T) + lambda)^2 -
                          T(2) * lambda * alpha^2 / beta2 -
                          T(2) * alpha * z0
            lambda, iters, _, converged = _bisect_decreasing(f, T;
                atol, rtol, maxiters)
            u = u0 / (one(T) + lambda)
            z = z0 + lambda * alpha / beta2
            return _raw_singleton((u, _zero_like(v0, T), z), :zero_second_root;
                iterations=iters, converged=converged)
        end
        u = u0 / T(2)
        z = z0 + alpha / beta2
        r2 = -T(2) * alpha * z + a2 / T(4)
        r = _radius(r2, abs(T(2) * alpha * z) + a2)
        if iszero(r)
            return _raw_singleton((u, _zero_like(v0, T), z), :zero_second_threshold)
        end
        generator = direction -> begin
            d = _normalize_direction(direction, v0, T)
            (copy(u), r * d, z)
        end
        default = generator(_default_direction(u0, T))
        return _raw_sphere(generator, default, :zero_second_family;
            description="all unit directions for the second canonical component")
    end
end
