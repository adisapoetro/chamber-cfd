module ChamberWaterLily

using WaterLily, StaticArrays, LinearAlgebra, Statistics, SparseArrays
using TOML, JSON3, SHA, Dates, Printf

# Grid-sized vector reductions do not benefit from a second thread pool.
function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
end

include("Config.jl")
include("Geometry.jl")
include("Forcing.jl")
include("Flow.jl")
include("CutCells.jl")
include("Transport.jl")
include("LeafSource.jl")
include("Coupled.jl")

export load_config, validate_config, cycle_state, fan_configuration, nominal_volume
export chamber_body, build_flow, advance_flow!, run_airflow
export run_coupled

end
