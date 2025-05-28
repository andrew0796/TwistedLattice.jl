using LinearAlgebra

function randomSU2()::Matrix{ComplexF64}
    a = 2*pi*rand()
    theta = acos(2*rand()-1)
    phi = 2*pi*rand()

    return cos(a)*I(2) + 1im*sin(a)*[cos(theta) exp(-1im*phi)*sin(theta) ; exp(1im*phi)*sin(theta) -cos(theta)]
end

function randomembeddedSU2(i::Int64, j::Int64, N::Int64)::Matrix{ComplexF64}
    N < 2 && error("N must be at least 2")
    (i <= 0 || i >= N) && error("i must be between 1 and N-1")
    (j <= 1 || j > N)  && error("j must be between 2 and N")
    i >= j && error("i must be strictly less than j")

    a = 2*pi*rand()
    theta = acos(2*rand()-1)
    phi = 2*pi*rand()

    R::Matrix{ComplexF64} = I(N)

    R[i,i] = cos(a) + im*sin(a)*cos(theta)
    R[i,j] = im*sin(a)*sin(theta)*exp(-im*phi)
    R[j,i] = im*sin(a)*sin(theta)*exp(im*phi)
    R[j,j] = cos(a) - im*sin(a)*cos(theta)

    return R
end

function randomSU(N::Int64)::Matrix{ComplexF64}
    N < 2 && error("N must be at least 2")

    M::Matrix{ComplexF64} = I(N)
    for i = 1:N-1, j=i+1:N
        M = M*randomembeddedSU2(i,j,N)
    end
    return M
end

function enforceunitarity!(M::Matrix{ComplexF64})
    # use Gram-Schmidt on the columns to guarantee unitarity, then add an appropriate phase to make the determinant unity
    M[1,:] = M[1,:]/norm(M[1,:])
    for i = 2:size(M)[1]
        for j = 1:i-1
            M[i,:] = M[i,:] - dot(M[j,:],M[i,:])*M[j,:]
        end
        M[i,:] = M[i,:]/norm(M[i,:])
    end
    phase = exp(-im*angle(det(M))/size(M)[2])
    lmul!(phase, M)
end