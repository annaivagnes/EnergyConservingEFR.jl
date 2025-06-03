# EnergyConservingEFR

This repository implements the Evolve-Filter-Relax (**EFR**) formulation in the context of
incompressible Navier-Stokes simulations on a staggered Cartesian grid.

The main dependency is the package [IncompressibleNavierStokes.jl](https://github.com/agdestein/IncompressibleNavierStokes.jl), where the Navier-Stokes framework on staggered
grids is implemented.

## The Evolve-Filter-Relax formulation

The EFR formulation consists of three different steps, that can be implemented
in a succession, one after the other. At the time step $t_n$, the three steps
can be written as:

- **Evolve**: the solution is evolved in time using a time-stepping scheme. 

**Evolve equations**:
Let  $\mathbf{w}_{n+1}$ be the intermediate velocity. Then the scheme is: $$\frac{\mathbf{w}_{n+1} - \mathbf{u}_n}{\Delta t}+ (\mathbf{w}_{n+1} \cdot \nabla) \mathbf{w}_{n+1}- \nu \Delta \mathbf{w}_{n+1}+ \nabla p_{n+1} = 0,$$ with $\nabla \cdot \mathbf{w}_{n+1} = 0$, boundary conditions applied as needed.

- **Filter**: filter the intermediate velocity $\mathbf{w}_{n+1}$, to
obtain the filtered one, $\overline{\mathbf{w}}_{n+1}$. Different filters
may be defined, like for instance:

    - a standard differential elliptic filter first presented in
[Differential filters of elliptic type, by Germano](https://pubs.aip.org/aip/pfl/article/29/6/1757/943987).
For example, it can be written as: $$-2\delta^2 \Delta \overline{\mathbf{w}}_{n+1} + \overline{\mathbf{w}}_{n+1} = \mathbf{w}_{n+1},$$

with the corresponding velocity boundary conditions.
    - a simple face-average;
    - a Smagorinsky closure model. Indeed, a filter can be re-written as a closure additional term into the
    *evolve* equations, as showed in [Bridging Large Eddy Simulation and
Reduced Order Modeling of Convection-Dominated
Flows through Spatial Filtering: Review and Perspectives by Quaini et al.](https://arxiv.org/pdf/2407.00231),
where also many other filters can be found.
    - a filter defined in the frequencies domain, the novelty presented in this repository.

- **Relax**: combine filtered and unfiltered fields in a convex combination, obtaining the final velocity $\mathbf{u}_{n+1}$: $$\mathbf{u}_{n+1} = (1 - \chi) \mathbf{w}_{n+1} + \chi \overline{\mathbf{w}}_{n+1}.$$.