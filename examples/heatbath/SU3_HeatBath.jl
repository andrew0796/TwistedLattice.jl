using TwistedLattice
using Random

N = 3 # number of colours for SU(N) Yang-Mills
twists_14_23 = [1*(i==1 && j==4) + 1*(i==2 && j==3) for i=1:4, j=1:4]
twists_14_23 .= twists_14_23 - transpose(twists_14_23)

dims = (10,6,6,6)

nSweeps = 1
heatBathRate = 1
overrelaxationRate = 0
coolingRate = 1
restartRate = 0

algMinTemp = 0.0
algMaxTemp = Inf

initialTemp = 2.0
finalTemp = 0.0 # will never reach this temperature, as long as the minimum convergence is reached

temperatureUpdateCoef = 0.3

nIterations_thermalization = 100

defaultSchedule(t::Int, L::Lattice, params::MCParameters) = params.initialTemperature*(1-params.temperatureUpdateCoef)^t
thermalizeSchedule(t::Int, L::Lattice, params::MCParameters) = t < nIterations_thermalization ? params.initialTemperature : params.finalTemperature

params = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef))
params_thermalization = readMCparameters(Dict{String,Any}("nSweeps"=>1, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef))

stoppingConditionInitial = defaultstoppingcondition(0.0, 0.0)
stoppingConditionFinal = defaultstoppingcondition(N-0.5, 2e-8)

directory = "./"
if !isdir(directory)
    println("Directory doesn't exist, creating directory $directory")
    mkpath(directory)
end
for seed = 1:1
    # set the random seed (for reproducability)
    println("\nAt seed $seed")
    Random.seed!(seed)
    
    L = Lattice(N, dims, twists_14_23)

    output = "SU$(N)_$(dims[1])_$(dims[2])_$(dims[3])_$(dims[4])_twists_14_23_seed_$seed.h5"
    f = createdatafile(output)

    # thermalize at the initial temperature
    @time results = minimizeyangmills!(L, params_thermalization, thermalizeSchedule, stoppingConditionInitial)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)

    # minimize until convergence is sufficiently small    
    @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionFinal)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)
    close(f)
end