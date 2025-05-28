using LinearAlgebra
using Folds

function wilsonactiondensity(SL::Lattice, ndx::Int, mu::Int, nu::Int)::Float64
    return real(SL.N - centerbackground(SL, ndx, mu,nu) * tr(SL[ndx,mu] * SL[SL._neighbours_pos[ndx,mu],nu] * adjoint(SL[SL._neighbours_pos[ndx,nu],mu]) * adjoint(SL[ndx,nu])))
end

function wilsonaction(SL::Lattice)::Float64
    return 2*Folds.sum(wilsonactiondensity(SL, ndx,mu,nu) for ndx in 1:prod(SL.size), mu in 1:4 for nu in mu+1:4)
end

function electricmagneticwilsonaction(L::Lattice)::NTuple{2, Float64}
    SE = 2*Folds.sum(wilsonactiondensity(L, ndx,4,i) for ndx in 1:prod(L.size), i in 1:3)
    SB = 2*Folds.sum(wilsonactiondensity(L, ndx,i,j) for ndx in 1:prod(L.size), i in 1:2 for j in i+1:3)
    return (SE, SB)
end


function wilsonaction!(L::Lattice)::Float64
    SE, SB = electricmagneticwilsonaction(L)
    L.action = SE + SB
    L.electricAction = SE
    L.magneticAction = SB
    return L.action
end


function setactiondensity!(L::Lattice)
    for ndx = 1:prod(L.size)
        L.actionDensity[ndx] = 0.0
        for mu=1:3, nu=mu+1:4
            L.actionDensity[ndx] += 2*wilsonactiondensity(L, ndx, mu,nu)
        end
    end
end

function improvedactiondensity(L::Lattice, ϵ::Float64, ndx::Int, mu::Int, nu::Int)::Float64
    """ Improved action with minimal subleading terms at ϵ=0, and Wilson action at ϵ=1 """
    if ϵ == 1.0
        return wilsonactiondensity(L, ndx,mu,nu)
    else
        return (4-ϵ)/3*wilsonactiondensity(L, ndx,mu,nu) + (ϵ-1)/48*real(L.N - doubleplaquettecenterbackground(L,ndx,mu,nu)*tr(L[ndx,mu] * L[L._neighbours_pos[ndx,mu],mu] * L[L._neighbours_pos[L._neighbours_pos[ndx,mu],mu],nu] * L[L._neighbours_pos[L._neighbours_pos[L._neighbours_pos[ndx,mu],mu],nu],nu] * adjoint(L[L._neighbours_pos[L._neighbours_pos[L._neighbours_pos[ndx,mu],nu],nu],mu]) * adjoint(L[L._neighbours_pos[L._neighbours_pos[ndx,nu],nu],mu]) * adjoint(L[L._neighbours_pos[ndx,nu],nu]) * adjoint(L[ndx,nu])))
    end
end

function electricmagneticimprovedaction(L::Lattice, ϵ::Float64)::NTuple{2, Float64}
    SE = 2*Folds.sum(improvedactiondensity(L,ϵ, ndx,4,i) for ndx in 1:prod(L.size), i in 1:3)
    SB = 2*Folds.sum(improvedactiondensity(L,ϵ, ndx,i,j) for ndx in 1:prod(L.size), i in 1:2 for j in i+1:3)
    return (SE, SB)
end

function improvedaction!(L::Lattice, ϵ::Float64)::Float64
    SE, SB = electricmagneticimprovedaction(L, ϵ)
    L.action = SE + SB
    L.electricAction = SE
    L.magneticAction = SB
    return L.action
end

function setimprovedactiondensity!(L::Lattice, ϵ::Float64)
    for ndx = 1:prod(L.size)
        L.actionDensity[ndx] = 0.0
        for mu=1:3, nu=mu+1:4
            L.actionDensity[ndx] += 2*improvedactiondensity(L,ϵ, ndx, mu,nu)
        end
    end
end

function calculatestaple(ndx::Int, mu::Int, L::Lattice)::Matrix{ComplexF64}
    M = zeros(ComplexF64, L.N,L.N)

    @simd for nu=1:4
        if nu != mu
            mul!(M, L[L._neighbours_pos[ndx,mu],nu], adjoint(L[ndx,nu] * L[L._neighbours_pos[ndx,nu],mu]), centerbackground(L, ndx, mu,nu), 1.0)
            mul!(M, adjoint(L[L._neighbours_neg[ndx,nu],mu] * L[L._neighbours_pos[L._neighbours_neg[ndx, nu],mu], nu]), L[L._neighbours_neg[ndx,nu],nu], conj(centerbackground(L, L._neighbours_neg[ndx,nu], mu,nu)), 1.0)
        end
    end
    return M
end

function calculatedoublestaple(ndx::Int, mu::Int, L::Lattice)::Matrix{ComplexF64}
    M = zeros(ComplexF64, L.N,L.N)

    @simd for nu=1:4
        if nu != mu
            mul!(M, L[L._neighbours_pos[ndx,mu],mu] * L[L._neighbours_pos[L._neighbours_pos[ndx,mu],mu],nu] * L[L._neighbours_pos[L._neighbours_pos[L._neighbours_pos[ndx,nu],mu],mu],nu], adjoint(L[ndx,nu] * L[L._neighbours_pos[ndx,nu],nu] * L[L._neighbours_pos[L._neighbours_pos[ndx,nu],nu],mu] * L[L._neighbours_pos[L._neighbours_pos[L._neighbours_pos[ndx,mu],nu],nu],mu]), doubleplaquettecenterbackground(L, ndx, mu, nu), 1.0)
            mul!(M, L[L._neighbours_pos[ndx,mu],nu] * L[L._neighbours_pos[L._neighbours_pos[ndx,mu],nu],nu] * adjoint(L[L._neighbours_neg[ndx,mu],nu] * L[L._neighbours_neg[L._neighbours_pos[ndx,nu],mu],nu] * L[L._neighbours_neg[L._neighbours_pos[L._neighbours_pos[ndx,nu],nu],mu],mu] * L[L._neighbours_pos[L._neighbours_pos[ndx,nu],nu],mu]), L[L._neighbours_neg[ndx,mu],mu], doubleplaquettecenterbackground(L, L._neighbours_neg[ndx,mu], mu, nu), 1.0)
            mul!(M, L[L._neighbours_pos[ndx,mu],mu] * adjoint(L[L._neighbours_neg[L._neighbours_neg[ndx,nu],nu],mu] * L[L._neighbours_neg[L._neighbours_neg[L._neighbours_pos[ndx,mu],nu],nu],mu] * L[L._neighbours_neg[L._neighbours_neg[L._neighbours_pos[L._neighbours_pos[ndx,mu],mu],nu],nu],nu] * L[L._neighbours_neg[L._neighbours_pos[L._neighbours_pos[ndx,mu],mu],nu],nu]) * L[L._neighbours_neg[L._neighbours_neg[ndx,nu],nu],nu], L[L._neighbours_neg[ndx,nu],nu], conj(doubleplaquettecenterbackground(L, L._neighbours_neg[L._neighbours_neg[ndx,nu],nu], mu,nu)), 1.0)
            mul!(M, adjoint(L[L._neighbours_neg[L._neighbours_neg[L._neighbours_neg[ndx,mu],nu],nu],mu] * L[L._neighbours_neg[L._neighbours_neg[ndx,nu],nu],mu] * L[L._neighbours_neg[L._neighbours_neg[L._neighbours_pos[ndx,mu],nu],nu],nu] * L[L._neighbours_pos[L._neighbours_neg[ndx,nu],mu],nu]) * L[L._neighbours_neg[L._neighbours_neg[L._neighbours_neg[ndx,mu],nu],nu],nu] * L[L._neighbours_neg[L._neighbours_neg[ndx,nu],mu],nu], L[L._neighbours_neg[ndx,mu],mu], conj(doubleplaquettecenterbackground(L, L._neighbours_neg[L._neighbours_neg[L._neighbours_neg[ndx,mu],nu],nu], mu, nu)), 1.0)
        end
    end

    return M
end

function calculateimprovedstaple(ndx::Int, mu::Int, L::Lattice, ϵ::Float64)::Matrix{ComplexF64}
    if ϵ == 1.0
        return calculatestaple(ndx, mu, L)
    else
        return lmul!((4-ϵ)/3, calculatestaple(ndx, mu, L)) .+ lmul!((ϵ-1)/48,calculatedoublestaple(ndx, mu, L))
    end
end