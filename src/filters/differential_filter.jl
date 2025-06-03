export DifferentialFilter

struct DifferentialFilter <: AbstractFilter
    filter_radius::Float64
end

function decompose_filter_mat(setup, filter_radius)
	D = diffusion_mat(setup)
	Id = sparse(I, size(D))
	filter_mat = Id-2*(filter_radius^2)*D
	lu(filter_mat)
end

function apply_filter(f::DifferentialFilter, u, setup, lu_filter_mats)
    ufilt = copy(u)

    D = diffusion_mat(setup)
    yu = apply_bc_u(zero(ufilt), 0.0, setup)
    B = bc_u_mat(setup)

    rhs = B * u[:] + 2 * (f.filter_radius^2) * D * yu[:]
    ufilt .= reshape(lu_filter_mats \ rhs, size(u))

    return ufilt
end