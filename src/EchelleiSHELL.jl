module EchelleiSHELL

using FITSIO, JLD2
using Infiltrator
using Polynomials
using OrderedCollections
using AstroAngles, SkyCoords
using NaNStatistics
using Glob
using Distributed
using Reexport

@reexport using Echelle
using EchelleReduce


# Exports
export iSHELLL0, iSHELLL1


# Basic info
const NAME = "iSHELL"
const OBSERVATORY = "IRTF"
const TIMEZONE = -10


# Gas Cell info
const GASCELL_FILE = "methane_gas_cell_ishell_kgas.jld"
const GASCELL_τ = 0.97


# LSF Info
const LSFσ_BOUNDS_KGAS_0375 = [0.011, 0.014] # units of nm


# Data types for ishell -> L0 and L1 for all images and extract spectra, respectively.
# No formal L2 support.
abstract type AnyiSHELL{L} <: SpecData{Symbol(lowercase(NAME)), L} end
struct iSHELLL0 <: AnyiSHELL{0}
    filename::String
end
struct iSHELLL1 <: AnyiSHELL{1}
    filename::String
end


# Read spectrum and header
include("parsing.jl")

# Starting wavelengths for iSHELL
include("wavelength.jl")

# Reduction
include("reduction.jl")

end
