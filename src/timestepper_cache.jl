
function efr_ode_method_cache(method::ExplicitRungeKuttaMethod, setup, filter = nothing)
    ustart = vectorfield(setup)
    ns = length(method.b)
    ku = map(i -> vectorfield(setup), 1:ns)
    p = scalarfield(setup)
    if isnothing(setup.temperature)
        tempstart = nothing
        ktemp = nothing
        diff = nothing
    else
        tempstart = scalarfield(setup)
        ktemp = map(i -> scalarfield(setup), 1:ns)
        diff = vectorfield(setup)
    end
    (; ustart, ku, p, tempstart, ktemp, diff)

    # Pre-compute LU decomposition before time-stepping
    println(filter)
    lu_filter_mats = if (isnothing(filter))
        nothing
    elseif (!isnothing(filter) && !isnothing(filter.filter_radius))
        decompose_filter_mat(setup, filter.filter_radius)
    end
    (; ustart, ku, p, tempstart, ktemp, diff, lu_filter_mats)
end