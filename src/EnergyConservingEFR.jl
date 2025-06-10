
module EnergyConservingEFR
using IncompressibleNavierStokes
using EnergyConservingEFR
using GLMakie
using Adapt
using Atomix: @atomic
using ChainRulesCore
using DocStringExtensions
using FFTW
using CSV
using DataFrames
using IterativeSolvers
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using LinearAlgebra
using Observables
using PrecompileTools
using Revise
using Printf
using Random
using SparseArrays
using StaticArrays
using Statistics
using WriteVTK: CollectionFile, paraview_collection, vtk_grid, vtk_save

# Docstring templates
@template MODULES = """
                   $(DOCSTRING)

                   ## Exports

                   $(EXPORTS)
                   """
@template (FUNCTIONS, METHODS) = """
                                 $TYPEDSIGNATURES

                                 $DOCSTRING
                                 """
@template TYPES = """
                  $TYPEDEF

                  $DOCSTRING

                  ## Fields

                  $FIELDS
                  """

"$LICENSE"
license = "MIT"

# # Easily retrieve value from Val
# (::Val{x})() where {x} = x

include("filters/abstract_filter.jl")
include("filters/average_filter.jl")
include("filters/differential_filter.jl")
include("filters/frequencies_filter.jl")
include("relax.jl")
include("utils.jl")
include("timestepper.jl")
include("solver.jl")
include("timestepper_cache.jl")
include("processors.jl")

export AbstractFilter
export FaceAverage, VolumeAverage
export DifferentialFilter
export FrequenciesFilter, find_f_star
export efr_timestep!
export efr_solve_unsteady
export efr_ode_method_cache
export AbstractRelax
export ConstantRelax, AdaptiveRelax
export dataframe_to_array, complex_dataframe_to_array, array_to_dataframe, total_enstrophy
export enstrophy!, total_enstrophy
export enstrophy_history_plot
end