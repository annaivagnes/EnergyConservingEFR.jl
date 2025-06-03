export AbstractFilter, apply_filter

abstract type AbstractFilter end

# Default call overload for filters
function (F::AbstractFilter)(u, setup; kwargs...)
    apply_filter(F, u, setup; kwargs...)
end