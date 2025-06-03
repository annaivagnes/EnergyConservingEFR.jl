# EnergyConservingEFR

This repository implements the Evolve-Filter-Relax (**EFR**) formulation in the context of
incompressible Navier-Stokes simulations on a staggered Cartesian grid.

The main dependency is the package [IncompressibleNavierStokes.jl](https://github.com/agdestein/IncompressibleNavierStokes.jl), where the Navier-Stokes framework on staggered
grids is implemented.

## The Evolve-Filter-Relax formulation

The EFR formulation consists of three different steps, that can be implemented
in a succession, one after the other. At the time step $t_n$, the three steps
can be written as:

1. **Evolve**: the solution is evolved in time using a time-stepping scheme. 

2. **Filter**: filter the intermediate velocity to
obtain the filtered one. Different filters
may be defined, like for instance:

    - a standard differential elliptic filter first presented in
[Differential filters of elliptic type, by Germano](https://pubs.aip.org/aip/pfl/article/29/6/1757/943987).
    - a simple face-average;
    - a Smagorinsky closure model. Indeed, a filter can be re-written as a closure additional term into the *evolve* equations, as showed in [Bridging Large Eddy Simulation and
Reduced Order Modeling of Convection-Dominated
Flows through Spatial Filtering: Review and Perspectives by Quaini et al.](https://arxiv.org/pdf/2407.00231),
where also many other filters can be found.

3. **Relax**: combine filtered and unfiltered fields in a convex combination, obtaining the final velocity.

An example of EFR formulation with the elliptic filter follows here:
![efr formulation](images/efr_readme.svg)