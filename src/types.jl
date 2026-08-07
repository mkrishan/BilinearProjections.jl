abstract type AbstractBilinearConstraint end

"""
    CrossConstraint()

The nonconvex cross

```math
C_0 = \{(x,y): \langle x,y\rangle = 0\}.
```
"""
struct CrossConstraint <: AbstractBilinearConstraint end

"""
    BilinearLevelSet(level)

The bilinear level set

```math
C_\gamma = \{(x,y): \langle x,y\rangle = \gamma\}.
```
"""
struct BilinearLevelSet{T<:Real} <: AbstractBilinearConstraint
    level::T
end

"""
    CanonicalHyperbolicParaboloid(alpha; beta=1)

The canonical rectangular hyperbolic paraboloid

```math
\widetilde C_\alpha
= \{(u,v,z): \|u\|^2-\|v\|^2 = 2\alpha z\},
```
with weighted product norm

```math
\|(u,v,z)\|_\beta^2 = \|u\|^2+\|v\|^2+\beta^2|z|^2.
```
"""
struct CanonicalHyperbolicParaboloid{T<:Real} <: AbstractBilinearConstraint
    alpha::T
    beta::T
    function CanonicalHyperbolicParaboloid(alpha::T, beta::T) where {T<:Real}
        iszero(alpha) && throw(ArgumentError("alpha must be nonzero"))
        beta > zero(T) || throw(ArgumentError("beta must be positive"))
        new{T}(alpha, beta)
    end
end

function CanonicalHyperbolicParaboloid(alpha::Real; beta::Real=1)
    T = promote_type(typeof(float(alpha)), typeof(float(beta)))
    return CanonicalHyperbolicParaboloid(T(alpha), T(beta))
end

"""
    HyperbolicParaboloid(alpha; beta=1)

The bilinear-coordinate rectangular hyperbolic paraboloid

```math
C_\alpha = \{(x,y,z): \langle x,y\rangle = \alpha z\},
```
with weighted product norm

```math
\|(x,y,z)\|_\beta^2 = \|x\|^2+\|y\|^2+\beta^2|z|^2.
```
"""
struct HyperbolicParaboloid{T<:Real} <: AbstractBilinearConstraint
    alpha::T
    beta::T
    function HyperbolicParaboloid(alpha::T, beta::T) where {T<:Real}
        iszero(alpha) && throw(ArgumentError("alpha must be nonzero"))
        beta > zero(T) || throw(ArgumentError("beta must be positive"))
        new{T}(alpha, beta)
    end
end

function HyperbolicParaboloid(alpha::Real; beta::Real=1)
    T = promote_type(typeof(float(alpha)), typeof(float(beta)))
    return HyperbolicParaboloid(T(alpha), T(beta))
end

"""
Metadata associated with a projection computation.
"""
struct ProjectionInfo{T<:Real}
    multiplier::Union{Nothing,T}
    distance2::T
    feasibility_residual::T
    stationarity_residual::T
    branch::Symbol
    multiplicity::Symbol
    iterations::Int
    converged::Bool
end

abstract type AbstractProjectionSet end

struct SingletonProjection{P,I<:ProjectionInfo} <: AbstractProjectionSet
    point::P
    info::I
end

"""
A complete parameterization of a non-singleton projection set.

`generator(direction)` returns a projected point. The direction may be any nonzero
array with the same axes as the relevant Hilbert-space variable; it is normalized
internally. `default` is the deterministic point returned by `project` unless another
selector is supplied.
"""
struct SphereProjection{F,P,E,I<:ProjectionInfo} <: AbstractProjectionSet
    generator::F
    default::P
    extras::E
    description::String
    info::I
end

abstract type AbstractSelector end
struct DefaultSelector <: AbstractSelector end
struct UnitDirection{V} <: AbstractSelector
    direction::V
end

Base.length(::SingletonProjection) = 1
Base.isempty(::AbstractProjectionSet) = false
(P::SphereProjection)(direction) = P.generator(direction)

function Base.show(io::IO, P::SingletonProjection)
    print(io, "SingletonProjection(branch=", P.info.branch,
          ", residual=", P.info.feasibility_residual, ")")
end

function Base.show(io::IO, P::SphereProjection)
    print(io, "SphereProjection(branch=", P.info.branch,
          ", description=\"", P.description, "\")")
end
