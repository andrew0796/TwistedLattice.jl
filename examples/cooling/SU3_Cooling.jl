using TwistedLattice
using Random

N = 3 # number of colours for SU(N) Yang-Mills
twists_14_23 = [1*(i==1 && j==4) + 1*(i==2 && j==3) for i=1:4, j=1:4]
twists_14_23 .= twists_14_23 - transpose(twists_14_23)

dims = (24,6,6,6)

nSweeps = 1
heatBathRate = 0
overrelaxationRate = 0
coolingRate = 1
restartRate = 0

algMinTemp = 0.0
algMaxTemp = Inf

initialTemp = 1.0
finalTemp = 0.0

temperatureUpdateCoef = 0.01

mod_twists = zeros(Int64, 4,4)
nIterations_modTwists = 70

defaultSchedule(t::Int, L::Lattice, params::MCParameters) = params.initialTemperature*(1-params.temperatureUpdateCoef)^t

params = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef))
initialParams = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>:polar))

stoppingConditionInitial = defaultstoppingcondition(Inf, 1e-4)
stoppingConditionMod(L::Lattice, actions::Array{Float64}, T::Float64, params::MCParameters) = length(actions)-1 >= nIterations_modTwists
stoppingConditionFinal = defaultstoppingcondition(3.0, 2e-8)

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

    # minimize using polar cooling method until convergence is < 10⁻⁴
    @time results = minimizeyangmills!(L, initialParams, defaultSchedule, stoppingConditionInitial; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)

    # parameters for doing sweeps with modified twists
    mod_params = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef))

    while L.action/(8*pi^2/N) > N-0.5
        # while it's stuck, reset twists to zero, try doing some iterations, and then try minimizing again with original twists
        L.twists = mod_twists
        
        @time results = minimizeyangmills!(L, mod_params, defaultSchedule, stoppingConditionMod; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

        savelattice!(L, f)
        dumpmetadata!(f, results)
        dumpMCparams!(f, mod_params)

        # switch twists back to the original ones
        L.twists = twists_14_23
        @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionInitial; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

        savelattice!(L, f)
        dumpmetadata!(f, results)
        dumpMCparams!(f, params)
    end

    # finish off the minimization
    @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionFinal; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)
    close(f)
end