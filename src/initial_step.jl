"""
    InitialStatic(; alpha=1.0)

Propose the same positive, finite step `alpha` on every search. A per-call
consumer bound may restrict that initial trial; no previous step is remembered.
"""
struct InitialStatic{T<:Real}
    alpha::T
end

InitialStatic(; alpha::Real=1.0) = InitialStatic(alpha)
