using Printf
using HDF5

function acceptanceprobability(oldAction::Float64, newAction::Float64, temperature::Float64)::Float64
    if newAction < oldAction
        return 1.0
    else   
        return exp(-(newAction-oldAction)/temperature)
    end
end

struct MCParameters
    nSweeps::Int64            # number of sweeps per minimization
    heatBathRate::Int64       # number of heat bath sweeps to do per algorithm sweep
    overrelaxationRate::Int64 # number of overrelaxation sweeps to do per algorithm sweep
    coolingRate::Int64        # number of cooling sweeps to do per algorithm sweep
    restartRate::Int64        # restart rate for simulated annealing

    parallelCooling::Bool     # whether or not to use parallelization for cooling
    coolingMethod::Symbol     # cooling method to use, either :default (which is SU2), :SU2, or :polar

    improvedAction::Float64   # ϵ for improved action, 0=improved, 1=Wilson action

    heatBathMinTemp::Float64       # minimum temperature to use heat bath
    heatBathMaxTemp::Float64       # maximum temperature to use heat bath
    overrelaxationMinTemp::Float64 # minimum temperature to use overrelaxation
    overrelaxationMaxTemp::Float64 # maximum temperature to use overrelaxation
    coolingMinTemp::Float64        # minimum temperature to use cooling
    coolingMaxTemp::Float64        # maximum temperature to use cooling

    initialTemperature::Float64    # initial simulated annealing temperature
    finalTemperature::Float64      # final simulated annealing temperature
    temperatureUpdateCoef::Float64 # coefficient for updating T -> T*(1-temperatureUpdateCoef)
end

function readMCparameters(params::Dict{String,Any})::MCParameters
    allowedKeys = ["nSweeps", "heatBathRate", "overrelaxationRate", "coolingRate", "restartRate", "parallelCooling", "coolingMethod", "improvedAction",
                   "heatBathMinTemp", "heatBathMaxTemp", "overrelaxationMinTemp", "overrelaxationMaxTemp", "coolingMinTemp", "coolingMaxTemp", 
                   "initialTemperature", "finalTemperature", "temperatureUpdateCoef"]
    defaultValues = [2, 1,0,3, 0, true, :default, 1.0, 0.0,Inf, 0.0,Inf, 0.0,Inf, 1.0,1e-7,0.01]
    orderedValues = []
    allowedValues = Dict{String,Any}("coolingMethod"=>[:default, :SU2, :polar])
    for (i,key) in enumerate(allowedKeys)
        if key ∉ keys(params)
            params[key] = defaultValues[i]
        end
        if key in keys(allowedValues) && params[key] ∉ allowedValues[key]
            @warn "value given ($(params[key])) for $key is not in the allowed values: $(allowedValues[key]), defaulting to $(defaultValues[i])"
            params[key] = defaultValues[i]
        end
        push!(orderedValues, params[key])
    end

    for key in keys(params)
        if key ∉ allowedKeys
            @warn "key \"$key\" not an accepted MC parameter"
        end
    end
    return MCParameters(orderedValues...)
end
    


function generatenewsite_heatbath(N::Int64, site::Matrix{ComplexF64}, staple::Matrix{ComplexF64}, temperature::Float64)::Matrix{ComplexF64}
    newSite = copy(site)
    U = staple*site
    temp = zeros(ComplexF64, 2,2)
    w = zeros(ComplexF64, 2,2)

    beta = 1/temperature

    for i=1:N-1, j=i+1:N
        temp .= [U[i,i] U[i,j] ; U[j,i] U[j,j]]
        temp .= (adjoint(temp).-temp.+I(2)*tr(temp))
        xi = sqrt(real(det(temp)))/2
        temp .= temp/(2*xi) # u^dagger

        xi *= 2/N

        # use rejection sampling
        a0 = 1 + 1/(beta*xi)*log(rand(exp(-2*beta*xi):eps():1))
        while 1 - sqrt(1-a0^2) < rand()
            a0 = 1 + 1/(beta*xi)*log(rand(exp(-2*beta*xi):eps():1))
        end
        theta = acos(2*rand()-1)
        phi = 2*pi*rand()
        sina = (-1)^(rand(0:1))*sqrt(1-a0^2)
        w = adjoint(temp)*(a0*I(2) + 1im*sina*[cos(theta) exp(-1im*phi)*sin(theta) ; exp(1im*phi)*sin(theta) -cos(theta)])

        for k=1:N
            temp0 = newSite[k,i]*w[1,1] + newSite[k,j]*w[2,1]
            temp1 = newSite[k,i]*w[1,2] + newSite[k,j]*w[2,2]
            newSite[k,i] = temp0
            newSite[k,j] = temp1

            temp0 = U[k,i]*w[1,1] + U[k,j]*w[2,1]
            temp1 = U[k,i]*w[1,2] + U[k,j]*w[2,2]
            U[k,i] = temp0
            U[k,j] = temp1
        end
    end

    enforceunitarity!(newSite)

    return newSite
end

function generatenewsite_cooling(N::Int64, site::Matrix{ComplexF64}, staple::Matrix{ComplexF64})::Matrix{ComplexF64}
    if N == 2
        newSite = staple/sqrt(abs(det(staple)))
    else
        U = staple*site
        newSite = copy(site)
        for i=1:N-1, j=i+1:N
            temp = [U[i,i] U[i,j] ; U[j,i] U[j,j]]
            xi = sqrt(real(det(temp - adjoint(temp) + I(2)*tr(adjoint(temp)))))/2
            temp = (adjoint(temp)-temp+I(2)*tr(temp))/(2*xi) # u^dagger

            for k=1:N
                temp0 = newSite[k,i]*temp[1,1] + newSite[k,j]*temp[2,1]
                temp1 = newSite[k,i]*temp[1,2] + newSite[k,j]*temp[2,2]
                newSite[k,i] = temp0
                newSite[k,j] = temp1

                temp0 = U[k,i]*temp[1,1] + U[k,j]*temp[2,1]
                temp1 = U[k,i]*temp[1,2] + U[k,j]*temp[2,2]
                U[k,i] = temp0
                U[k,j] = temp1
            end
        end
    end

    enforceunitarity!(newSite)

    return newSite
end

function generatenewsite_overrelaxation(N::Int64, site::Matrix{ComplexF64}, staple::Matrix{ComplexF64})::Matrix{ComplexF64}
    X = -adjoint(staple)
    SVD = svd(X) # X = SVD.U * Diagonal(SVD.S) * SVD.Vt
    phi = angle(det(X))
    Xhat = exp(-im*phi/N)*(SVD.U)*(SVD.Vt)
    newSite = Xhat * adjoint(site) * Xhat
    
    enforceunitarity!(newSite)

    return newSite
end


function heatbathupdatesite!(ndx::Int, mu::Int, L::Lattice, temperature::Float64)::Int64
    M = calculatestaple(ndx, mu, L)

    currentAction = 2*real(3*L.N-tr(L[ndx, mu]*M))
    newSite = generatenewsite_heatbath(L.N, L[ndx, mu], M, temperature)
    newAction = 2*real(3*L.N-tr(newSite*M))

    if acceptanceprobability(currentAction, newAction, temperature) >= rand()
        L[ndx, mu] = newSite
        return 1
    else
        return 0
    end
end


function heatbathsweep!(L::Lattice, temperature::Float64)::Int64
    nAccepted = 0
    for i=1:prod(L.size), mu=1:4
        ndx = rand(1:prod(L.size))
        nu = rand(1:4)
        nAccepted += heatbathupdatesite!(ndx, nu, L, temperature)
    end
    return nAccepted
end

function overrelaxationupdatesite!(ndx::Int, mu::Int, L::Lattice)
    M = calculatestaple(ndx, mu, L)
    newSite = generatenewsite_overrelaxation(L.N, L[ndx, mu], M)
    L[ndx, mu] = newSite
end

function overrelaxationsweep!(L::Lattice)
    for i=1:prod(L.size), mu=1:4
        overrelaxationupdatesite!(i, mu, L)
    end
end

function coolingupdatesite_SU2subgroups!(ndx::Int, mu::Int, L::Lattice, ϵ::Float64)
    if L.N == 2
        newSite = zeros(ComplexF64, L.N,L.N)
        @simd for nu in 1:4
            if nu != mu
                mul!(newSite, L[ndx,nu] * L[L._neighbours_pos[ndx,nu],mu], adjoint(L[L._neighbours_pos[ndx,mu],nu]), conj(centerbackground(L, ndx, mu, nu)), 1.0)
                mul!(newSite, adjoint(L[L._neighbours_neg[ndx,nu],nu]) * L[L._neighbours_neg[ndx,nu],mu], L[L._neighbours_neg[L._neighbours_pos[ndx,mu],nu],nu], centerbackground(L, L._neighbours_neg[ndx,nu], mu, nu), 1.0)
            end
        end
        newSite = newSite/sqrt(abs(det(newSite)))
        L[ndx, mu] = newSite
    else
        staple = calculateimprovedstaple(ndx, mu, L, ϵ)
        U = staple*L[ndx, mu]
        for i=1:L.N-1, j=i+1:L.N
            temp = [U[i,i] U[i,j] ; U[j,i] U[j,j]]
            xi = sqrt(real(det(temp - adjoint(temp) + I(2)*tr(adjoint(temp)))))/2
            temp = (adjoint(temp)-temp+I(2)*tr(temp))/(2*xi) # u^dagger

            for k=1:L.N
                temp0 = L.sites[ndx,mu,k,i]*temp[1,1] + L.sites[ndx,mu,k,j]*temp[2,1]
                temp1 = L.sites[ndx,mu,k,i]*temp[1,2] + L.sites[ndx,mu,k,j]*temp[2,2]
                L.sites[ndx,mu,k,i] = temp0
                L.sites[ndx,mu,k,j] = temp1

                temp0 = U[k,i]*temp[1,1] + U[k,j]*temp[2,1]
                temp1 = U[k,i]*temp[1,2] + U[k,j]*temp[2,2]
                U[k,i] = temp0
                U[k,j] = temp1
            end
        end
    end

    enforceunitarity!(ndx,mu, L)
end

#= function coolingupdatesite_polar!(x::Vector{Int64}, L::Lattice, ϵ::Float64)
    staple = calculateimprovedstaple(x, L, ϵ)
    F = svd(staple)
    L[x...] = F.V*adjoint(F.U)

    enforceunitarity!(x, L)
end
 =#
function coolingupdatesite_polar!(ndx::Int, mu::Int, L::Lattice, ϵ::Float64)
    staple = calculateimprovedstaple(ndx, mu, L, ϵ)
    F = svd(staple)
    L[ndx, mu] = F.V*adjoint(F.U)

    enforceunitarity!(ndx, mu, L)
end

function coolingsweep!(L::Lattice, method::Symbol, ϵ::Float64)
    for i=1:prod(L.size), mu=1:4
        if method == :polar
            coolingupdatesite_polar!(i, mu, L, ϵ)
        else
            coolingupdatesite_SU2subgroups!(i, mu, L, ϵ)
        end
    end
end

function parallelcoolingsweep!(L::Lattice, method::Symbol, ϵ::Float64)
    for mu=1:4
        @sync begin
            for indices in Iterators.partition(L._checkerboard_red, div(length(L._checkerboard_red), Threads.nthreads()))
                if method == :polar
                    Threads.@spawn for ndx in indices
                        coolingupdatesite_polar!(ndx, mu, L, ϵ)
                    end
                else
                    Threads.@spawn for ndx in indices
                        coolingupdatesite_SU2subgroups!(ndx, mu, L, ϵ)
                    end
                end
            end
        end

        @sync begin
            for indices in Iterators.partition(L._checkerboard_blk, div(length(L._checkerboard_blk), Threads.nthreads()))
                if method == :polar
                    Threads.@spawn for ndx in indices
                        coolingupdatesite_polar!(ndx, mu, L, ϵ)
                    end
                else
                    Threads.@spawn for ndx in indices
                        coolingupdatesite_SU2subgroups!(ndx, mu, L, ϵ)
                    end
                end
            end
        end
    end
end

function improvedparallelcoolingsweep!(L::Lattice, method::Symbol, ϵ::Float64)
    for mu=1:4
        for i=1:length(L._ranges_improved[mu])
            Threads.@threads for ndx in L._ranges_improved[mu][i]
                if method == :polar
                    coolingupdatesite_polar!(ndx, mu, L, ϵ)
                else
                    coolingupdatesite_SU2subgroups!(ndx, mu, L, ϵ)
                end
            end
        end
    end
end

function defaultstoppingcondition(stoppingAction::Float64, minimumConvergence::Float64)::Function
    f(L::Lattice, actions::Array{Float64}, T::Float64, params::MCParameters) = !(T > params.finalTemperature && (actions[end] > stoppingAction  || length(actions) == 1 || (length(actions) > 1 && (actions[end-1] - actions[end]) > minimumConvergence)))
    return f
end

function minimizeyangmills!(L::Lattice, params::MCParameters, schedule::Function, stoppingCondition::Function; quiet::Bool=false, progressBar::Bool=true, progressBarFrequency::Int64=1, snapshotFrequency::Int64=0, snapshotFile::Union{HDF5.File,Nothing}=nothing, overwriteSnapshots::Bool=false, printOutput::IO=stdout)::Dict{String,<:Any}
    """ 
    Minimize the Yang-Mills action for lattice L, returning a dictionary of MC metadata
    
    L: Lattice object to be minimized
    params: MCParameters object, likely created by calling ReadMCParameters
    schedule: simulated annealing schedule, which is a function of time (number of iterations completed), the lattice itself (as a Lattice), and params. Typically will be T0*(1-ε)^time where ε=params.temperatureUpdateCoef and T0=params.initialTemperature
    stoppingCondition::(Lattice, actions, temperature, params) -> Bool: when to stop the algorithm. An example would be f(L, actions, T, params) = !(T > params.finalTemperature && (actions[end] > stoppingAction  || length(actions) == 1 || (length(actions) > 1 && (actions[end-1] - actions[end]) > minimumConvergence)))

    quiet: whether or not to print information
    progressBar: whether or not to display a progress bar
    progressBarFrequency: how often to update the progressBar
    snapshotFrequency: how often to save snapshots
    snapshotFile: file to save snapshots, or nothing
    overwriteSnapshots: whether to overwrite the previous snapshot or not
    printOutput: IO to print output to
    """
    if !hasmethod(schedule, (Int, Lattice, MCParameters))
        error("schedule given does not have the correct method, should be a function taking an integer (time), the lattice, and the MCParameters")
    end
    if !hasmethod(stoppingCondition, (Lattice, Array{Float64}, Float64, MCParameters))
        error("stoppingFunction given does not have the correct method, should be a function taking a Lattice, Array{Float64}, Float64, MCParameters")
    end
    if snapshotFrequency > 0 && isnothing(snapshotFile)
        error("snapshotFrequency is positive, but there is no file to save snapshots to")
    end
    if progressBarFrequency <= 0 && progressBar
        @warn "progressBarFrequency non-positive, defaulting to 1"
        progressBarFrequency = 1
    end

    action = zeros(Float64, 0)
    push!(action, improvedaction!(L, params.improvedAction)/(8*pi^2/L.N))

    if !quiet
        println(printOutput, "Initial action/(8π²/$(L.N)): $(L.electricAction/(8*pi^2/L.N)) + $(L.magneticAction/(8*pi^2/L.N)) = $(action[1])")
        flush(printOutput)
    end

    acceptRate_heatBath = zeros(Float64, 0)

    T = params.initialTemperature
    temps = [T]

    nSites = prod(L.size)*4

    bestAction = Inf
    bestLattice = similar(L.sites)

    while !stoppingCondition(L, action, T, params)
        nAccepted_heatBath = 0

        for n = 1:params.nSweeps
            if params.overrelaxationRate > 0 && (params.overrelaxationMinTemp <= T <= params.overrelaxationMaxTemp)
                for m = 1:params.overrelaxationRate
                    overrelaxationsweep!(L)
                end
            end

            if params.coolingRate > 0 && (params.coolingMinTemp <= T <= params.coolingMaxTemp)
                for m = 1:params.coolingRate
                    if params.parallelCooling
                        if params.improvedAction == 1.0
                            parallelcoolingsweep!(L, params.coolingMethod, params.improvedAction)
                        else
                            improvedparallelcoolingsweep!(L, params.coolingMethod, params.improvedAction)
                        end
                    else
                        coolingsweep!(L, params.coolingMethod, params.improvedAction)
                    end
                end
            end

            if params.heatBathRate > 0 && (params.heatBathMinTemp <= T <= params.heatBathMaxTemp)
                for m = 1:params.heatBathRate
                    nAccepted_heatBath += heatbathsweep!(L, T)
                end
            end
        end

        params.heatBathRate > 0 && push!(acceptRate_heatBath, nAccepted_heatBath/(nSites*params.nSweeps*params.heatBathRate))

        push!(action, improvedaction!(L, params.improvedAction)/(8*pi^2/L.N))
        push!(temps, T)

        T = schedule(length(action)-1, L, params)

        if progressBar && (length(action)-1)%progressBarFrequency == 0
            if params.heatBathRate > 0 
                @printf "\r\033[K\r%i iterations, action/(8π^2/%i) = %.14f, temperature = %.3e, metropolis acceptance rate = %.4f%%, log(convergence) = (%c)%.4e" (length(action)-1) L.N action[end] T (acceptRate_heatBath[end]*100) charactersign(action[end]-action[end-1]) log10(abs(action[end]-action[end-1]))
            else
                printprogressbar(printOutput, L, action)
            end
        end

        if params.restartRate > 0 && action[end] < bestAction
            bestAction = action[end]
            copy!(bestLattice, L.sites)
        end
        if params.restartRate > 0 && length(action) % params.restartRate == 0
            copy!(L.sites, bestLattice)
        end

        if snapshotFrequency > 0 && mod(length(action)-1, snapshotFrequency) == 0
            savesnapshot!(L, snapshotFile; overwriteprevious=overwriteSnapshots)
        end
    end
    if progressBar
        print(printOutput, "\n")
        flush(printOutput)
    end

    setimprovedactiondensity!(L, params.improvedAction)

    if !quiet
        println(printOutput, "Final action/(8π²/$(L.N)): $(L.electricAction/(8*pi^2/L.N)) + $(L.magneticAction/(8*pi^2/L.N)) = $(action[end])")
        flush(printOutput)
    end

    return Dict("action"=>action, "temperature"=>temps, "acceptRate_heatBath"=>acceptRate_heatBath)
end