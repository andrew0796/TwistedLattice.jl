using TwistedLattice
using Test

tolerance = 1e-1

function test_cooling(dims::NTuple{4, Int}, method::Symbol, N::Int, twists::Matrix{Int}, improvement_parameter::Float64)
    nSweeps = 1
    heatBathRate = 0
    coolingRate = 1

    initialTemp = 1.0
    finalTemp = 0.0

    temperatureUpdateCoef = 0.01

    defaultSchedule(t::Int, L::Lattice, params::MCParameters) = t*0.1

    params = readMCparameters(Dict{String,Any}("nSweeps"=>nSweeps, "coolingRate"=>coolingRate, "heatBathRate"=>heatBathRate, "initialTemperature"=>initialTemp, "finalTemperature"=>finalTemp, "temperatureUpdateCoef"=>temperatureUpdateCoef, "coolingMethod"=>method, "improvedAction"=>improvement_parameter))

    stoppingCondition = defaultstoppingcondition(1.0*N, 1e-5)

    L = Lattice(N, dims, twists)

    results = minimizeyangmills!(L, params, defaultSchedule, stoppingCondition)

    improvedaction!(L, improvement_parameter) 
    S = L.action/(8*pi^2/N)

    return S < N && ≈(S-round(S), 0, atol = tolerance)
end

@testset "TwistedLattice.jl" begin
    twists_14_23 = [1*(i==1 && j==4) + 1*(i==2 && j==3) for i=1:4, j=1:4]
    twists_14_23 .= twists_14_23 - transpose(twists_14_23)

    twists_zero = zeros(Int64, 4,4)

    dims = (4,4,4,4)

    # test Wilson action cooling
    @test test_cooling(dims, :SU2, 2, twists_zero, 1.0)
    @test test_cooling(dims, :SU2, 3, twists_zero, 1.0)
    @test test_cooling(dims, :polar, 2, twists_zero, 1.0)
    @test test_cooling(dims, :polar, 3, twists_zero, 1.0)
    @test test_cooling(dims, :SU2, 2, twists_14_23, 1.0)
    @test test_cooling(dims, :SU2, 3, twists_14_23, 1.0)
    @test test_cooling(dims, :polar, 2, twists_14_23, 1.0)
    @test test_cooling(dims, :polar, 3, twists_14_23, 1.0)

    # test improved action cooling
    @test test_cooling(dims, :SU2, 2, twists_zero, 0.0)
    @test test_cooling(dims, :SU2, 3, twists_zero, 0.0)
    @test test_cooling(dims, :polar, 2, twists_zero, 0.0)
    @test test_cooling(dims, :polar, 3, twists_zero, 0.0)
    @test test_cooling(dims, :SU2, 2, twists_14_23, 0.0)
    @test test_cooling(dims, :SU2, 3, twists_14_23, 0.0)
    @test test_cooling(dims, :polar, 2, twists_14_23, 0.0)
    @test test_cooling(dims, :polar, 3, twists_14_23, 0.0)
end
