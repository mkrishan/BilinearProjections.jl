"""
    projection_set(C, point; kwargs...)

Return the complete metric projection. The result is a `SingletonProjection`
or a `SphereProjection`.
"""
function projection_set end

"""
    residual(C, point)

Signed constraint residual at `point`.
"""
function residual end

"""
    select(P, selector=DefaultSelector())

Select one point from a projection-set representation. `UnitDirection(d)` is
available for `SphereProjection` objects.
"""
select(P::SingletonProjection, ::DefaultSelector=DefaultSelector()) = P.point
select(P::SphereProjection, ::DefaultSelector=DefaultSelector()) = P.default
select(P::SphereProjection, selector::UnitDirection) = P.generator(selector.direction)
select(::SingletonProjection, ::UnitDirection) =
    throw(ArgumentError("UnitDirection is only meaningful for a non-singleton projection set"))

"""
    project(C, point; selector=DefaultSelector(), return_info=false, kwargs...)

Return a deterministic nearest point in the constraint set. Set
`return_info=true` to also receive `ProjectionInfo`. For a complete nonunique
answer, use `projection_set`.
"""
function project(C::AbstractBilinearConstraint, point;
                 selector::AbstractSelector=DefaultSelector(),
                 return_info::Bool=false, kwargs...)
    P = projection_set(C, point; kwargs...)
    projected = select(P, selector)
    return return_info ? (projected, P.info) : projected
end

"""
    project!(output, C, point; kwargs...)

Write a selected projection into preallocated arrays. Pair constraints accept
`(xout,yout)`. Paraboloids accept `(xout,yout,Ref(zout))`. The method returns
`ProjectionInfo`.
"""
function project!(output::Tuple{<:AbstractArray,<:AbstractArray},
                  C::Union{CrossConstraint,BilinearLevelSet}, point; kwargs...)
    projected, info = project(C, point; return_info=true, kwargs...)
    copyto!(output[1], projected[1])
    copyto!(output[2], projected[2])
    return info
end

function project!(output::Tuple{<:AbstractArray,<:AbstractArray,<:Base.RefValue},
                  C::Union{CanonicalHyperbolicParaboloid,HyperbolicParaboloid}, point; kwargs...)
    projected, info = project(C, point; return_info=true, kwargs...)
    copyto!(output[1], projected[1])
    copyto!(output[2], projected[2])
    output[3][] = projected[3]
    return info
end

"""
    distance2(C, point; kwargs...)

Squared distance from `point` to the constraint set.
"""
distance2(C::AbstractBilinearConstraint, point; kwargs...) =
    projection_set(C, point; kwargs...).info.distance2

"""
    isfeasible(C, point; atol=0, rtol=sqrt(eps(Float64)))

Test feasibility using a scale-aware residual tolerance.
"""
function isfeasible(C::AbstractBilinearConstraint, point;
                    atol::Real=0, rtol::Real=sqrt(eps(Float64)))
    r = abs(residual(C, point))
    scale = one(float(r)) + _residual_scale(C, point)
    return r <= atol + rtol * scale
end

function _public_projection(info_builder, raw)
    multiplicity = raw.kind === :singleton ? :singleton : :continuum
    default_point = raw.kind === :singleton ? raw.point : raw.default
    info = info_builder(default_point, multiplicity, raw.branch,
                        raw.iterations, raw.converged)
    if raw.kind === :singleton
        return SingletonProjection(raw.point, info)
    end
    return SphereProjection(raw.generator, raw.default, raw.extras,
                            raw.description, info)
end
