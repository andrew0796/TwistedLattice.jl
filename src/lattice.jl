using LinearAlgebra

mutable struct Lattice
    """ Four dimensional Lattice, using lexicographic indices. This is quite a bit faster than the Cartesian indexing, though less intuitive """
    N::Int64
    sites::Array{ComplexF64, 4} # (N1 x N2 x N3 x N4) x 4 x NxN
    size::NTuple{4, Int64}
    twists::Matrix{Int64} # 4x4 antisymmetric matrix of twists
    
    action::Float64
    electricAction::Float64
    magneticAction::Float64
    actionDensity::Array{Float64}

    _neighbours_pos::Array{Int, 2} # V x 4 array of neighbouring indices in position directions
    _neighbours_neg::Array{Int, 2} # V x 4 array of neighbouring indices in negative directions

    _checkerboard_red::Array{Int} # V/2 array of "red" sites (where the sum of all coordinates is even)
    _checkerboard_blk::Array{Int} # V/2 array of "black" sites (where the sum of all coordinates is odd)

    _ranges_improved::Array{Vector{Vector{Int}}} # array of sites for iterating over safely in parallel with improved cooling

    Lattice() = new()
end


function singleindextocartesian(i::Int, L::Lattice)::Vector{Int}
    ndx = zeros(Int, 4)
    i_temp = i
    for mu = 4:-1:1
        ndx[mu] = div(i_temp-1, prod(L.size[1:mu-1])) + 1
        i_temp = i_temp - (ndx[mu]-1)*prod(L.size[1:mu-1])
    end
    return ndx
end

function cartesiantosingleindex(x::Vector{Int}, L::Lattice)::Int
    ndx = 1
    for mu=1:4
        ndx += (x[mu]-1)*prod(L.size[1:mu-1])
    end
    return ndx
end

function cartesiantosingleindex(x::Vector{Int}, dimensions::NTuple{4, Int})::Int
    ndx = 1
    for mu=1:4
        ndx += (x[mu]-1)*prod(dimensions[1:mu-1])
    end
    return ndx
end

function cartesiantosingleindex(L::Lattice, x::Vararg{Int, 5})::Int
    ndx = 1
    for mu=1:4
        ndx += (x[mu]-1)*prod(L.size[1:mu-1])
    end
    return ndx
end

@inline Base.getindex(L::Lattice, ndx::Vararg{Int,2}) = getindex(L.sites, ndx..., :,:)
@inline Base.setindex!(L::Lattice, M::Matrix{ComplexF64}, ndx::Vararg{Int,2}) = setindex!(L.sites, M, ndx..., :,:)

@inline Base.getindex(L::Lattice, ndx::Vararg{Int, 5}) = getindex(L.sites, cartesiantosingleindex(L, ndx...), ndx[5], :,:)
@inline Base.setindex!(L::Lattice, M::Matrix{ComplexF64}, ndx::Vararg{Int,5}) = setindex!(L.sites, M, cartesiantosingleindex(L, ndx...), ndx[5], :,:)


function setneighbours!(SIL::Lattice)
    SIL._neighbours_neg = Array{Int, 2}(undef, prod(SIL.size), 4)
    SIL._neighbours_pos = Array{Int, 2}(undef, prod(SIL.size), 4)
    for i1=1:SIL.size[1], i2=1:SIL.size[2], i3=1:SIL.size[3], i4=1:SIL.size[4], mu=1:4
        site = [i1,i2,i3,i4]
        ndx = cartesiantosingleindex(site, SIL.size)
        if site[mu] == SIL.size[mu]
            SIL._neighbours_pos[ndx, mu] = ndx - (SIL.size[mu]-1)*prod(SIL.size[1:mu-1])
        else
            SIL._neighbours_pos[ndx, mu] = ndx + prod(SIL.size[1:mu-1])
        end
        if site[mu] == 1
            SIL._neighbours_neg[ndx, mu] = ndx + (SIL.size[mu]-1)*prod(SIL.size[1:mu-1])
        else
            SIL._neighbours_neg[ndx, mu] = ndx - prod(SIL.size[1:mu-1])
        end
    end
end

function setcheckerboards!(SIL::Lattice)
    SIL._checkerboard_red = Array{Int}(undef, div(prod(SIL.size), 2))
    SIL._checkerboard_blk = Array{Int}(undef, div(prod(SIL.size), 2))
    ndx_red = 1
    ndx_blk = 1
    for i1=1:SIL.size[1], i2=1:SIL.size[2], i3=1:SIL.size[3], i4=1:SIL.size[4]
        site = [i1,i2,i3,i4]
        ndx = cartesiantosingleindex(site, SIL.size)
        if iseven(sum(site))
            SIL._checkerboard_red[ndx_red] = ndx
            ndx_red += 1
        else
            SIL._checkerboard_blk[ndx_blk] = ndx
            ndx_blk += 1
        end
    end
    if ndx_red != div(prod(SIL.size), 2) + 1 || ndx_blk != div(prod(SIL.size), 2) + 1
        error("something went wrong with constructing checkerboards! got to $ndx_red and $ndx_blk for red and black")
    end
    sort!(SIL._checkerboard_red)
    sort!(SIL._checkerboard_blk)
end

function createranges_improvedparallel(dims::NTuple{4,Int64}, indices::Vector{Int64})::Vector{Vector{UnitRange{Int64}}}
    """ create the ranges needed for parallelizing with the improved action """
    ranges = Vector{UnitRange{Int64}}[]

    for i=1:3
        push!(ranges, [0:2*fld(dims[indices[i]],4)-1])
        if iseven(div(dims[indices[i]],2)) # don't need the extra range at the end, just add an empty range
            push!(ranges[end], 1:0)
        else
            push!(ranges[end], div(dims[indices[i]]-2,2):div(dims[indices[i]]-2,2))
        end
    end
    return ranges
end

function setimprovedranges!(L::Lattice)
    L._ranges_improved = Array{Vector{Vector{Int}}}(undef, 4)

    for mu=1:4
        indices = [1,2,3,4]
        filter!(e->e!=mu, indices)

        ranges = createranges_improvedparallel(L.size, indices)
        L._ranges_improved[mu] = Vector{Int}[]
        
        for a=0:3, b1=0:1, b2=0:1, b3=0:1
            for i=1:2, j=1:2, k=1:2
                range = []
                for i1=ranges[1][i], i2=ranges[2][j], i3=ranges[3][k]
                    for i4=cld(-2*(i1+i2+i3)-a,4):fld(L.size[mu]-2*(i1+i2+i3)-a-1,4)
                        x = (2*i1+b1+1)*I[1:4,indices[1]] + (2*i2+b2+1)*I[1:4,indices[2]] + (2*i3+b3+1)*I[1:4,indices[3]] + (4*i4 + 2*(i1+i2+i3) + a + 1)*I[1:4,mu]
                        ndx = cartesiantosingleindex(x, L)
                        push!(range, ndx)
                    end
                end
                if length(range) > 0
                    sort!(range)
                    push!(L._ranges_improved[mu], range)
                end
            end
        end
    end
end


function randominitializelattice!(L::Lattice)
    for i in eachindex(IndexCartesian(), L.sites[:,:,1,1])
        L.sites[i,:,:] = randomSU(L.N)
    end
end

function Lattice(N::Int64, dimensions::NTuple{4, Int64}, twists::Matrix{Int64}; initialization::Symbol=:random)
    L = Lattice()
    N < 2 && error("N must be at least 2")
    L.N = N
    L.size = dimensions
    L.sites = Array{ComplexF64, 4}(undef, prod(dimensions), 4, N,N)

    L.action = 0.0
    L.electricAction = 0.0
    L.magneticAction = 0.0
    L.actionDensity = zeros(Float64, prod(dimensions))

    if initialization == :random
        randominitializelattice!(L)
    end

    if size(twists) != (4,4)
        error("twists must be a 4x4 (antisymmetric) matrix, given a matrix with size $(size(twists))")
    end
    L.twists = twists
    for i = 1:3, j= i+1:4
        L.twists[j,i] = -L.twists[i,j]
    end

    setneighbours!(L)
    setcheckerboards!(L)
    setimprovedranges!(L)

    return L
end

function neighbour_pos_index(ndx::Int, mu::Int, L::Lattice)::Int
    return L._neighbours_pos[ndx, mu]
end

function neighbour_neg_index(ndx::Int, mu::Int, L::Lattice)::Int
    return L._neighbours_neg[ndx, mu]
end

function centerbackground(L::Lattice, ndx::Int, mu::Int, nu::Int)::ComplexF64
    if mu == nu
        return 0.0+0.0im
    end
    x = singleindextocartesian(ndx, L) 
    if x[mu] == x[nu] == 1
        return exp(-2im*pi*L.twists[mu,nu]/L.N)
    else
        return 1.0+0.0im
    end
end

function doubleplaquettecenterbackground(L::Lattice, ndx::Int, mu::Int, nu::Int)::ComplexF64
    if mu == nu
        return 0.0+0.0im
    end
    x = singleindextocartesian(ndx, L) 
    if (x[mu] == 1 || x[mu] == L.size[mu]) && (x[nu] == 1 || x[nu] == L.size[nu])
        return exp(-2im*pi*L.twists[mu,nu]/L.N)
    else
        return 1.0+0.0im
    end
end

function enforceunitarity!(ndx::Int, mu::Int, L::Lattice)
    # use Gram-Schmidt on the rows to guarantee unitarity, then add an appropriate phase to make the determinant unity
    L.sites[ndx,mu,1,:] = L.sites[ndx,mu,1,:]/norm(L.sites[ndx,mu,1,:])
    for i = 2:L.N
        L.sites[ndx,mu,i,:] = L.sites[ndx,mu,i,:] - sum(dot(L.sites[ndx,mu,j,:],L.sites[ndx,mu,i,:])*L.sites[ndx,mu,j,:] for j=1:i-1)
        L.sites[ndx,mu,i,:] = L.sites[ndx,mu,i,:]/norm(L.sites[ndx,mu,i,:])
    end
    phase = exp(-im*angle(det(L[ndx,mu]))/L.N)
    L[ndx,mu] = phase*L[ndx,mu]
end

function shiftlattice!(L::Lattice, distance::Int64, direction::Int64)
    """ Shift lattice a certain distance along a particular direction by incrementally moving the twists and shifting the lattice sites """
    if direction < 1 || direction > 4
        error("direction must be between 1 and 4, given $direction")
    end

    if distance == 0
        return
    end

    indices = [1,2,3,4]
    filter!(e->e!=direction, indices)

    new_sites = copy(L.sites)
    new_action_density = copy(L.actionDensity)

    for n=1:abs(distance)
        for j in indices
            if mod(L.twists[direction, j], L.N) != 0
                otherIndices = filter(e->e!=j, indices)
                for i1=1:L.size[otherIndices[1]], i2=1:L.size[otherIndices[2]]
                    if distance > 0
                        L[(I[1:5,direction] + I[1:5,j] + i1*I[1:5,otherIndices[1]]+ i2*I[1:5,otherIndices[2]] + j*I[1:5,5])...] *= exp(2im*pi*L.twists[direction,j]/L.N)
                    else
                        L[(2*I[1:5,direction] + I[1:5,j] + i1*I[1:5,otherIndices[1]]+ i2*I[1:5,otherIndices[2]] + j*I[1:5,5])...] *= exp(-2im*pi*L.twists[direction,j]/L.N)
                    end
                end
            end
        end

        for ndx = 1:prod(L.size), mu=1:4
            if distance > 0
                new_sites[ndx, mu, :,:] = L[L._neighbours_neg[ndx, direction], mu]
                new_action_density[ndx] = L.actionDensity[L._neighbours_neg[ndx, direction]]
            else
                new_sites[ndx, mu, :,:] = L[L._neighbours_pos[ndx, direction], mu]
                new_action_density[ndx] = L.actionDensity[L._neighbours_pos[ndx, direction]]
            end
        end
        L.sites = copy(new_sites)
        L.actionDensity = copy(new_action_density)
    end
end

function shiftlattice!(L::Lattice, distance::Vector{Int64})
    if length(distance) != 4
        error("distance should be a 4-dimensional vector, given $distance")
    end
    for i=1:4
        if distance[i] != 0
            shiftlattice!(L, distance[i], i)
        end
    end
end

function centerlattice!(L::Lattice)
    """ Center lattice around the maxima of the action density in each direction, calculate action density if necessary """
    maxima = argmax(L.actionDensity)

    shiftlattice!(L, [fld(L.size[i],2)-maxima[i] for i=1:4])
end

function settwist!(L::Lattice, mu::Int, nu::Int, twist::Int)
    """ Set the mu-nu component of the twists of L to be twist (mod N), automatically setting the nu-mu component so L.twists is antisymmetric """
    if mu == nu
        @warn "mu and nu are the same component, not setting twist"
        return
    end
    if mu < 1 || mu > 4 || nu < 1 || nu > 4
        error("mu and nu must be between 1 and 4, given $mu and $nu")
    end
    L.twists[mu,nu] = mod(twist, L.N)
    L.twists[nu,mu] = -L.twists[mu,nu]
    nothing
end

function settwists!(L::Lattice, twists::Matrix{Int})
    """ Set the twists of the lattice L. twists must be a 4x4 antisymmetric (up to mod L.N) matrix """
    if size(twists) != (4,4)
        error("twists must be a 4x4 matrix, given $(size(twists))")
    end
    if all(mod.(twists+transpose(twists), L.N) .!= 0)
        error("twists must be an antisymmetric matrix, up to mod(L.N)")
    end
    L.twists = twists
    nothing
end