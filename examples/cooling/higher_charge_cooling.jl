using TwistedLattice
using Random

N = 3 # number of colours for SU(N) Yang-Mills
twists_14_23 = [1*(i==1 && j==4) + 1*(i==2 && j==3) for i=1:4, j=1:4]
twists_14_23 .= twists_14_23 - transpose(twists_14_23)

dims = (20,6,6,6)

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

defaultSchedule(t::Int, L::Lattice, params::MCParameters) = params.initialTemperature*(1-params.temperatureUpdateCoef)^t

# parameter for improved action, less than 0 to stabilize higher charge instantons
ϵ_improved = -1.0

params = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>:polar))
finalParams = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>:polar, "improvedAction"=>ϵ_improved))

stoppingConditionInitial = defaultstoppingcondition(2.0*N, 1e-3) # stop when the action converges and is less than the two-instanton action
stoppingConditionFinal = defaultstoppingcondition(2.0*N, 3e-6)

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

    output = "SU$(N)_$(dims[1])_$(dims[2])_$(dims[3])_$(dims[4])_twists_14_23_higher_charge_seed_$seed.h5"
    f = createdatafile(output)

    # minimize using polar cooling method until convergence is < 10⁻³
    @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionInitial; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)
    

    # finish off the minimization with the improved action to stabilize higher charges
    improvedaction!(L, ϵ_improved)
    setimprovedactiondensity!(L, ϵ_improved)
    @time results = minimizeyangmills!(L, finalParams, defaultSchedule, stoppingConditionFinal; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)
    close(f)
end