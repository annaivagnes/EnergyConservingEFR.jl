
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
    lu_filter_mats = nothing 
    if (!isnothing(filter) && hasproperty(filter, :filter_radius))
        lu_filter_mats = decompose_filter_mat(setup, filter.filter_radius)
    end
    # if the filter is a FrequenciesFilter, we need to compute the filter matrix
    if (!isnothing(filter) && isa(filter, FrequenciesFilter))
        if !isfile(filter.filename)
            f_star = find_f_star(setup, filter.time_train, filter.Δt_train, filter.every_train, filter.nseeds)
            f_star = reshape(f_star, setup.grid.Nu[1][1]+2, setup.grid.Nu[1][2]+2, 2)
            CSV.write(filter.filename, array_to_dataframe(f_star))
        else
            f_star = CSV.read(filter.filename, DataFrame)
            f_star = complex_dataframe_to_array(f_star, setup.grid.Nu[1][1]+2, setup.grid.Nu[1][1]+2, 2)
        end
    end
    (; ustart, ku, p, tempstart, ktemp, diff, lu_filter_mats, f_star)
end