module EnergyConservingEFRMakieExt
using IncompressibleNavierStokes: Dimension, scalewithvolume!, scalarfield, vorticity
using EnergyConservingEFR: enstrophy!, total_enstrophy
using Makie
using Observables
import EnergyConservingEFR: enstrophy_history_plot

function enstrophy_history_plot(state; setup)
    @assert state isa Observable "Enstrophy history requires observable state."
    (; Ip) = setup.grid
    e = scalarfield(setup)
    _points = Point2f[]
    points = lift(state) do (; u, t)
        E = total_enstrophy(u, setup)
        push!(_points, Point2f(t, E))
    end
    fig = lines(points; axis = (; xlabel = "t", ylabel = "Total enstrophy"))
    on(_ -> autolimits!(fig.axis), points)
    fig
end
export enstrophy_history_plot

end