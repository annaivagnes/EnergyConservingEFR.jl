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

# computation of enstrophy
enstrophy(u, setup; kwargs...) =
    enstrophy!(scalarfield(setup), u, setup; kwargs...)

ChainRulesCore.rrule(::typeof(enstrophy), u, setup; kwargs...) =
    (enstrophy(u, setup; kwargs...), φ -> error("Not yet implemented"))


function enstrophy!(ens_e, u, setup)
    (; grid, backend, workgroupsize) = setup
    (; dimension, Np, Ip) = grid
    D = dimension()
    e = Offset(D)
    ω = vorticity(u, setup) 
    @kernel function efirst!(ens_e, u, I0)
        I = @index(Global, Cartesian)
        I = I + I0
        ens_e[I] = ω[I]^2
    end
    ens_e! = efirst!
    I0 = getoffset(Ip)
    ens_e!(backend, workgroupsize)(ens_e, u, I0; ndrange = Np)
    ens_e
end

function total_enstrophy(u, setup; kwargs...)
    (; Ip) = setup.grid
    ens = enstrophy(u, setup; kwargs...)
    ens = scalewithvolume(ens, setup)
    sum(view(ens, Ip))
end