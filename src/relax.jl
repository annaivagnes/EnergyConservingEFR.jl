export AbstractRelax, apply_relax, ConstantRelax, AdaptiveRelax

abstract type AbstractRelax end

# Default call overload for relax step
function (r::AbstractRelax)(u, ufilt; kwargs...)
    apply_relax(r, u, ufilt; kwargs...)
end

struct ConstantRelax <: AbstractRelax
    χ::Float64
end

function apply_relax(r::ConstantRelax, u, ufilt)
    return (1 - r.χ) .* u .+ r.χ .* ufilt
end

struct AdaptiveRelax <: AbstractRelax
    compute_chi::Function
end

function apply_relax(r::AdaptiveRelax, u, ufilt)
    χ = r.compute_chi(u, ufilt, setup)
    return (1 - χ) .* u .+ χ .* ufilt
end