# # Decaying Homogeneous Isotropic Turbulence - 2D
#
# In this example we consider decaying homogeneous isotropic turbulence,
# similar to the cases considered in [Kochkov2021](@cite) and
# [Kurz2022](@cite). The initial velocity field is created randomly, but with a
# specific energy spectrum. Due to viscous dissipation, the turbulent features
# eventually group to form larger visible eddies.

# ## Packages
#
# We just need IncompressibleNavierStokes and a Makie plotting backend.
                         #src

#md using CairoMakie
using GLMakie #!md
using IncompressibleNavierStokes
using EnergyConservingEFR
using Random
using CSV, DataFrames
using FFTW
Random.seed!(1234)

# Setup
n_dns = 512
n = 128
compression = 4
ax_dns = LinRange(0.0, 1.0, n_dns + 1)
ax = LinRange(0.0, 1.0, n + 1)
Re = 4e4
kolm = Re^(-3/4)
setup_dns = Setup(; x = (ax_dns, ax_dns), Re = Re,)
setup = Setup(; x = (ax, ax), Re = Re,)

# Filter (via face average) the initial velocity of the DNS simulation
filter_dns = FaceAverage(compression)
ustart_dns = random_field(setup_dns, 0.0);
ustart = apply_filter(filter_dns, ustart_dns, setup);

# Define Filter and relax
filter = DifferentialFilter(kolm)
relax = ConstantRelax(1.)

# Solve unsteady problem
state, outputs = efr_solve_unsteady(;
    setup,
    ustart,
	tlims = (0.0, 0.35),
	Δt = 0.0015,
    filter = filter,
    relax = relax,
    project_after_filter = true,
    processors = (
        rtp = realtimeplotter(; setup, nupdate = 1),
        ehist = realtimeplotter(;
            setup,
            plot = energy_history_plot,
            nupdate = 2,
            displayfig = false,
        ),
        espec = realtimeplotter(;
            setup,
            plot = energy_spectrum_plot,
            nupdate = 2,
            displayfig = false,
        ),
		anim = animator(;
						setup,
						path="examples/animations/standard_efr.mp4",
						fieldname = :vorticity,
						size=(500, 500),
						nupdate=2,
						),
        log = timelogger(; nupdate = 2),
    ),
);

#md # ```@raw html
#md # <video src="/DecayingTurbulence2D.mp4" controls="controls" autoplay="autoplay" loop="loop"></video>
#md # ```

# ## Post-process
#
# We may visualize or export the computed fields
u = state.u
