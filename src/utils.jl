function dataframe_to_array(df, Nx, Ny, Nz)
    A = zeros(Nx, Ny, Nz)  # Create an empty array
    
    for row in eachrow(df)
        i, j = row.i, row.j  # Extract indices
        A[i, j, 1] = row.value1
        A[i, j, 2] = row.value2
    end
    return A
end

function complex_dataframe_to_array(df, Nx, Ny, Nz)
    A = zeros(ComplexF64, Nx, Ny, Nz)  # Create an empty array
    
    for row in eachrow(df)
        i, j = row.i, row.j  # Extract indices
        A[i, j, 1] = parse(ComplexF64, row.value1) 
        A[i, j, 2] = parse(ComplexF64, row.value2)
    end
    return A
end

function array_to_dataframe(A)
    Nx, Ny, Nz = size(A)  # Get the size of the array
    df = DataFrame(i=Int[], j=Int[], value1=Float64[], value2=Float64[])  # Empty DataFrame

    for i in 1:Nx
        for j in 1:Ny
            push!(df, (i, j, A[i, j, 1], A[i, j, 2]), promote=true)  # Flatten into rows
        end
    end

    return df
end