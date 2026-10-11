"""
    InitialStatic(; alpha=1.0)

Start every search at the same positive, finite step `alpha`.
"""
struct InitialStatic{T<:Real}
    alpha::T
end

InitialStatic(; alpha::Real=1.0) = InitialStatic(alpha)
