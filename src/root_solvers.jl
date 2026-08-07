"""
Bisection for a continuous strictly decreasing function on `(-1,1)`.
The caller must provide a problem for which the left limit is positive and the
right limit is negative.
"""
function _bisect_decreasing(f, ::Type{T}; atol=nothing, rtol=nothing, maxiters::Int=200) where {T<:AbstractFloat}
    maxiters > 0 || throw(ArgumentError("maxiters must be positive"))
    atolT = isnothing(atol) ? T(16) * eps(T) : T(atol)
    rtolT = isnothing(rtol) ? T(32) * eps(T) : T(rtol)
    atolT >= zero(T) || throw(ArgumentError("atol must be nonnegative"))
    rtolT >= zero(T) || throw(ArgumentError("rtol must be nonnegative"))

    lo = nextfloat(-one(T))
    hi = prevfloat(one(T))
    flo = f(lo)
    fhi = f(hi)

    (flo > zero(T) || isinf(flo) && flo > 0) ||
        throw(ErrorException("root was not bracketed at the left endpoint; f(lo)=$flo"))
    (fhi < zero(T) || isinf(fhi) && fhi < 0) ||
        throw(ErrorException("root was not bracketed at the right endpoint; f(hi)=$fhi"))

    mid = (lo + hi) / T(2)
    fmid = f(mid)
    for k in 1:maxiters
        mid = (lo + hi) / T(2)
        fmid = f(mid)
        if iszero(fmid) || (hi - lo) <= atolT + rtolT * max(one(T), abs(mid))
            return mid, k, abs(fmid), true
        end
        if fmid > zero(T)
            lo = mid
        else
            hi = mid
        end
    end
    return mid, maxiters, abs(fmid), false
end
