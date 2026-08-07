function _check_array(x::AbstractArray)
    eltype(x) <: Real || throw(ArgumentError("only real arrays are supported"))
    isempty(x) && throw(ArgumentError("Hilbert-space variables must be nonempty"))
    return nothing
end

function _check_pair(x::AbstractArray, y::AbstractArray)
    _check_array(x)
    _check_array(y)
    axes(x) == axes(y) || throw(DimensionMismatch("x and y must have the same axes"))
    return nothing
end

function _float_type(x::AbstractArray, y::AbstractArray, scalars::Real...)
    T = promote_type(typeof(float(zero(eltype(x)))),
                     typeof(float(zero(eltype(y)))))
    for s in scalars
        T = promote_type(T, typeof(float(s)))
    end
    return T
end

_asarray(x::AbstractArray, ::Type{T}) where {T<:Real} = T.(x)
_sumsq(x::AbstractArray) = real(dot(x, x))

function _zero_like(x::AbstractArray, ::Type{T}) where {T<:Real}
    z = similar(x, T)
    fill!(z, zero(T))
    return z
end

function _first_unit(x::AbstractArray, ::Type{T}) where {T<:Real}
    e = _zero_like(x, T)
    e[firstindex(e)] = one(T)
    return e
end

function _normalize_direction(direction::AbstractArray, template::AbstractArray, ::Type{T}) where {T<:Real}
    axes(direction) == axes(template) ||
        throw(DimensionMismatch("direction must have the same axes as the projected variable"))
    eltype(direction) <: Real || throw(ArgumentError("direction must be real"))
    d = T.(direction)
    nd = norm(d)
    nd > zero(T) || throw(ArgumentError("direction must be nonzero"))
    return d / nd
end

function _default_direction(reference::AbstractArray, ::Type{T}) where {T<:Real}
    r = T.(reference)
    nr = norm(r)
    return nr > zero(T) ? r / nr : _first_unit(reference, T)
end

function _canonical_coordinates(x::AbstractArray, y::AbstractArray, ::Type{T}) where {T<:Real}
    s = sqrt(T(2))
    return ((T.(x) + T.(y)) / s, (-T.(x) + T.(y)) / s)
end

function _bilinear_coordinates(u::AbstractArray, v::AbstractArray)
    T = promote_type(eltype(u), eltype(v))
    s = sqrt(T(2))
    return ((u - v) / s, (u + v) / s)
end

function _radius(r2::T, scale::T) where {T<:Real}
    tol = T(64) * eps(T) * max(one(T), abs(scale))
    r2 < -tol && throw(ErrorException("internal branch produced a negative squared radius: $r2"))
    abs(r2) <= tol && return zero(T)
    return sqrt(max(zero(T), r2))
end

function _distance2_pair(input, point)
    x0, y0 = input
    x, y = point
    return _sumsq(x - x0) + _sumsq(y - y0)
end

function _distance2_triple(input, point, beta)
    x0, y0, z0 = input
    x, y, z = point
    return _sumsq(x - x0) + _sumsq(y - y0) + beta^2 * abs2(z - z0)
end

function _map_raw(raw, mapper)
    if raw.kind === :singleton
        return (; kind=:singleton,
                point=mapper(raw.point),
                branch=raw.branch,
                iterations=raw.iterations,
                converged=raw.converged,
                description=raw.description)
    end
    generator = direction -> mapper(raw.generator(direction))
    extras = map(mapper, raw.extras)
    return (; kind=:sphere,
            generator,
            default=mapper(raw.default),
            extras,
            branch=raw.branch,
            iterations=raw.iterations,
            converged=raw.converged,
            description=raw.description)
end

_raw_singleton(point, branch; iterations=0, converged=true, description="unique projection") =
    (; kind=:singleton, point, branch, iterations, converged, description)

_raw_sphere(generator, default, branch;
            extras=(), iterations=0, converged=true,
            description="unit-sphere parameterized projection family") =
    (; kind=:sphere, generator, default, extras, branch, iterations, converged, description)
