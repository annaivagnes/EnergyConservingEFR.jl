
export FaceAverage, VolumeAverage

struct FaceAverage <: AbstractFilter
    compression::Int
end

struct VolumeAverage <: AbstractFilter
    compression::Int
end

function apply_filter(f::FaceAverage, u, setup_les)
    compression = f.compression
    (; grid, backend, workgroupsize) = setup_les
    (; dimension, Nu, Iu) = grid
    D = dimension()

    v = vectorfield(setup_les)

    @kernel function Φ!(v, u, ::Val{α}, face, I0) where {α}
        I = @index(Global, Cartesian)
        J = I0 + compression * (I - oneunit(I))
        s = zero(eltype(v))
        for i in face
            s += u[J + i, α]
        end
        v[I0 + I, α] = s / compression^(D - 1)
    end

    for α = 1:D
        ndrange = Nu[α]
        I0 = getoffset(Iu[α])
        face = CartesianIndices(ntuple(β -> β == α ? (compression:compression) : (1:compression), D))
        Φ!(backend, workgroupsize)(v, u, Val(α), face, I0; ndrange)
    end

    return v
end

function apply_filter(f::VolumeAverage, u, setup_les)
    compression = f.compression
    (; grid, boundary_conditions, backend, workgroupsize) = setup_les
    (; dimension, N, Nu, Iu) = grid
    D = dimension()

    v = vectorfield(setup_les)

    @assert all(bc -> bc[1] isa PeriodicBC && bc[2] isa PeriodicBC, boundary_conditions)
    @kernel function Φ!(v, u, ::Val{α}, volume, I0) where {α}
        I = @index(Global, Cartesian)
        J = I0 + compression * (I - oneunit(I))
        s = zero(eltype(v))
        # n = 0
        for i in volume
            # Periodic extension
            K = J + i
            K = mod1.(K.I, compression .* (N .- 2))
            K = CartesianIndex(K)
            s += u[K, α]
            # n += 1
        end
        n = (iseven(comp) ? compression + 1 : comp) * comp^(D - 1)
        v[I0+I, α] = s / n
    end
    for α = 1:D
        ndrange = Nu[α]
        I0 = getoffset(Iu[α])
        volume = CartesianIndices(
            ntuple(
                β ->
                    α == β ?
                    iseven(compression) ? (div(compression, 2):div(compression, 2)+compression) :
                    (div(compression, 2)+1:div(compression, 2)+compression) : (1:compression),
                D,
            ),
        )
        Φ!(backend, workgroupsize)(v, u, Val(α), volume, I0; ndrange)
    end
    v
end