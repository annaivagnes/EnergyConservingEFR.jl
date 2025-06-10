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
    df = DataFrame(i=Int[], j=Int[], value1=ComplexF64[], value2=ComplexF64[])  # Empty DataFrame

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


#function enstrophy!(ens_e, u, setup)
#    (; grid, backend, workgroupsize) = setup
#    (; dimension, Np, Ip) = grid
#    D = dimension()
#    e = IncompressibleNavierStokes.Offset(D)
#    ω = vorticity(u, setup) 
#    @kernel function efirst!(ens_e, ω, I0)
#        I = @index(Global, Cartesian)
#        I = I + I0
#        ens_e[I] = ω[I]^2
#    end
#    ens_e! = efirst!
#    I0 = getoffset(Ip)
#    ens_e!(backend, workgroupsize)(ens_e, ω, I0; ndrange = Np)
#    ens_e
#end

function enstrophy!(ens, u, setup; interpolate_first = false)
    (; grid, backend, workgroupsize) = setup
    (; dimension, Np, Ip) = grid
    D = dimension()
    e = IncompressibleNavierStokes.Offset(D)
    ω = vorticity(u, setup)

    @kernel function efirst!(ens, ω, I0)
        I = @index(Global, Cartesian)
        I = I + I0
        v = zero(eltype(ens))
        v += (ω[I] + ω[I-e(1)])^2
        ens[I] = v / 8
    end

    @kernel function elast!(ens, ω, I0)
        I = @index(Global, Cartesian)
        I = I + I0
        v = zero(eltype(ens))
        v += ω[I]^2 + ω[I-e(1)]^2
        ens[I] = v / 4
    end

    ens! = interpolate_first ? efirst! : elast!
    I0 = getoffset(Ip)
    ens!(backend, workgroupsize)(ens, ω, I0; ndrange = Np)
    ens
end


function total_enstrophy(u, setup; kwargs...)
    (; Ip) = setup.grid
    ens = enstrophy(u, setup; kwargs...)
    ens = scalewithvolume(ens, setup)
    0.5 * sum(view(ens, Ip))
end