using TwistedLattice
using Random
using HDF5

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

# schedule for simulated annealing, not really needed here
defaultSchedule(t::Int, L::Lattice, params::MCParameters) = params.initialTemperature*(1-params.temperatureUpdateCoef)^t

# parameter for improved action
ϵ_improved = 0.0

# parameters for different cooling methods and actions (Wilson vs improved)
params_SU2 = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef))
params_polar = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>:polar))
params_SU2_improved = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "improvedAction"=>0.0))
params_polar_improved = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>:polar, "improvedAction"=>0.0))

# stopping conditions
stoppingConditionInitial = defaultstoppingcondition(Inf, 1e-4) # stop when convergence falls below 1e-4, at any action
stoppingConditionMod(L::Lattice, actions::Array{Float64}, T::Float64, params::MCParameters) = length(actions)-1 >= nIterations_modTwists # stop when nIterations_modTwists have been reached
stoppingConditionFinal = defaultstoppingcondition(3.0, 2e-8) # stop when the action (in units of 8π²/N) falls below 3.0 and the convergence is less than 2e-8

# for improved action, stop when the action converges OR the first derivative converges
minimumAction = 3.0
minimumConvergence = 1e-10
minimumSecondConvergence = 7e-12
stoppingCondition_improved(L::Lattice, action::Array{Float64}, T::Float64, params::MCParameters) = length(action) > 3 && action[end] < minimumAction && (action[end-1]-action[end] < minimumConvergence || action[end]-2*action[end-1]+action[end-2] < minimumSecondConvergence)


# clean up, this just makes sure that any open files are closed nicely
open_files = Union{HDF5.File, IOStream}[]
function cleanup()
    println("running cleanup!")
    for file in open_files
        close(file)
    end
end
atexit(cleanup)

# for randomizing twists
function newtwists(oldtwists::Matrix{Int}, N::Int)::Matrix{Int}
    if size(oldtwists) != (4,4)
        error("old twists must be 4x4")
    end
    T = rand(Int, size(oldtwists))
    T = mod.(T - transpose(T), N)
    while all(mod.(T - oldtwists, N) .== 0)
        T = rand(Int, size(oldtwists))
        T = mod.(T - transpose(T), N)
    end
    return T
end

function cool_lattice(params::MCParameters, params_improved::MCParameters, label::String)
    output = "SU$(N)_$(dims[1])_$(dims[2])_$(dims[3])_$(dims[4])_twists_14_23_$label.h5"
    progressfile = replace(output, ".h5"=>"_progress.txt")

    if isfile(output)
        println("Deleting file $output and recreating")
        rm(output)
    end

    # create datafile
    close(createdatafile(output))

    # initialize file with SWMR read/write
    f = h5open(output, "r+")
    push!(open_files, f)
    try
        HDF5.start_swmr_write(f) # use swmr to protect from corruption
    catch e
        if !isa(e, HDF5.API.H5Error)
            rethrow(e)
        end
    end

    progressoutput = open(progressfile, "w")
    push!(open_files, progressoutput)

    # set random seed
    Random.seed!(1)
    L = Lattice(N, dims, twists_14_23)

    # minimize using cooling method until convergence is < 10⁻⁴
    @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionInitial; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true, progressBarFrequency=20, printOutput=progressoutput)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)

    while L.action/(8*pi^2/N) > N-0.5
        # while it's stuck, randomize twists, try doing some iterations, and then try minimizing again with original twists
        L.twists = newtwists(L.twists, L.N)
        @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionMod; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true, progressBarFrequency=20, printOutput=progressoutput)

        savelattice!(L, f)
        dumpmetadata!(f, results)
        dumpMCparams!(f, params)

        # switch twists back to the original twists
        L.twists = twists_14_23
        @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionInitial; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true, progressBarFrequency=20, printOutput=progressoutput)

        savelattice!(L, f)
        dumpmetadata!(f, results)
        dumpMCparams!(f, params)
    end

    # finish off the minimization
    @time results = minimizeyangmills!(L, params, defaultSchedule, stoppingConditionFinal; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true, progressBarFrequency=20, printOutput=progressoutput)

    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)

    # improve action
    improvedaction!(L, ϵ_improved)
    setimprovedactiondensity!(L, ϵ_improved)
    @time results = minimizeyangmills!(L, params_improved, defaultSchedule, stoppingCondition_improved; snapshotFrequency=50, snapshotFile=f, overwriteSnapshots=true, progressBarFrequency=20, printOutput=progressoutput)
    
    savelattice!(L, f)
    dumpmetadata!(f, results)
    dumpMCparams!(f, params)

    close(progressoutput)
    pop!(open_files)

    close(f)
    pop!(open_files)
end



################################################################################################
# Test different cooling methods
################################################################################################

# SU(2) subgroups
cool_lattice(params_SU2, params_SU2_improved, "SU2_subgroups")

# polar decomposition
cool_lattice(params_polar, params_polar_improved, "polar_decomposition")