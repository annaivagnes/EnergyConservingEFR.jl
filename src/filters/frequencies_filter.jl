export FrequenciesFilter
export find_f_star

struct FrequenciesFilter <: AbstractFilter
    time_train::Float64
    Δt_train::Float64
    every_train::Int
    nseeds::Int
    filename::String
end

function find_f_star(setup, time_train, Δt_train, every_train, nseeds)
    nsnaps_seed = Int(time_train / (Δt_train * every_train))
    nsnaps = nsnaps_seed*nseeds
    (; grid, boundary_conditions) = setup
    (; dimension, x, N, Np, Nu, Ip, Iu, Δ, Δu) = grid
    n = Nu[1][1]

    # if files are not found: run DNS, filter it (indipendently at each time step) and save
    if !isfile("csv_files/csv_ref_filtered_long/seed_1_uref_0.0000.csv")
        println("Files of filtered DNS not found, run DNS and filter it through FaceAverage at first!")
        # TODO: run DNS and filter it through FaceAverage at first
    else
        println("Files of filtered DNS found!")
    end

    uref_end = zeros(n+2, n+2, 2, Int(nsnaps))
    wevolved = zeros(n+2, n+2, 2, Int(nsnaps))
    global i_index = 1
    for seed_ in 1:nseeds
        println("Processing seed: $seed_")
        for i in 1:every_train:nsnaps-1

            # read the data from the csv files at time step i and i+1
            local formatted_t_start = @sprintf("%.4f", i * Δt_train)
            local formatted_t_end = @sprintf("%.4f", (i + 1) * Δt_train)
            local uread_start = CSV.read("csv_files/csv_ref_filtered_long/seed_$(seed_)_uref_$(formatted_t_start).csv", DataFrame)
            local uread_end = CSV.read("csv_files/csv_ref_filtered_long/seed_$(seed_)_uref_$(formatted_t_end).csv", DataFrame)
            uref_end[:, :, :, i_index] .= dataframe_to_array(uread_end, n+2, n+2, 2)
            
            # evolve the filtered DNS data from time step i to i+1
            local ustart_i = reshape(dataframe_to_array(uread_start, n+2, n+2, 2), (n+2, n+2, 2))
            local state_, outputs_ = solve_unsteady(;
                                    setup=setup,
                                    tlims=(0.0, Δt_train),
                                    ustart=ustart_i,
                                    Δt=Δt_train,
                                    )
            wevolved[:, :, :, i_index] .= state_.u
            global i_index += 1
        end
    end

    # Compute the Fourier coefficients for the difference between the evolved and reference data
    dims = length(size(wevolved)) - 2 #should be 2 for 2D
    W_hat = fft(wevolved, collect(1:dims)) 
    res_hat = fft(uref_end, collect(1:dims))

    # Least-squares solution for f_star
    f_star_num = sum((conj.(W_hat) .* res_hat), dims = [dims+2])
    f_star_den = sum((W_hat .* conj.(W_hat)), dims = [dims+2])
    f_star = f_star_num ./ ( f_star_den)[[(:) for i in 1:dims+1]...,1]
    f_star
end

function apply_filter(f::FrequenciesFilter, u, setup, f_star)
    # Extract grid and boundary information from the setup
    (; grid, boundary_conditions) = setup
    (; dimension, x, N, Np, Nu, Ip, Iu, Δ, Δu) = grid
	dims = dimension()
	local u_reshaped = reshape(u, (size(u)[1], size(u)[2], dims, 1))
	u_new = real.(ifft(f_star .* fft(u_reshaped, collect(1:dims)), collect(1:dims)))
	ufilt = reshape(u_new, size(u))
	return ufilt
end