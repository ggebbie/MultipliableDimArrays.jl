module MultipliableDimArrays

using DimensionalData
using DimensionalData: @dim
using LinearAlgebra
using Unitful
using UnitfulLinearAlgebra

export MultipliableDimArray
export DiagonalDimArray
export Matrix
export uniform
export \, /
export eigen
export exp

import Base: Matrix
import Base: \, /
import Base: exp 
import LinearAlgebra: eigen 

# an alias
#MultiDimArray{T} = DimArray{T} where T <: AbstractDimArray
@dim Eigenmode "eigenmode"

"""
function Matrix(P::DimArray{T}) where T <: AbstractDimArray
"""
function Matrix(P::DimArray{T}) where T <: DimArray
    # number of columns/ outer dims
    N = length(P)
    # number of rows, take first inner element as example
    M = length(first(P))
    A = Array{eltype(first(P))}(undef,M,N)
    if N > 1  
        for j in eachindex(P)
            A[:,j] = P[j][:]
        end
    elseif N == 1
        for i in eachindex(first(P))
            A[i,1] = first(P)[i]
        end
    end
    return A 
end

"""
function MultipliableDimArray(A::AbstractMatrix,rangedims,domaindims)
"""
function MultipliableDimArray(A::AbstractMatrix,rangedims,domaindims)

    # extra step for type stability
    Q1 = reshape(A[:,1],size(rangedims))
    P1 = DimArray(Q1, rangedims)

    P = Array{typeof(P1)}(undef,size(domaindims))
    for j in eachindex(P)
        Q = reshape(A[:,j],size(rangedims))
        P[j] = DimArray(Q, rangedims)
    end
    return DimArray(P, domaindims)
end

function MultipliableDimArray(A::AbstractVector, rdims)
    Q1 = reshape(A,size(rdims))
    return DimArray(Q1, rdims)
end

# Note: there is no equivalent function for Multipliable Vectors,
# because it would conflict with default functions in DimensionalData.jl.
# Use a Matrix to represent a Vector if this functionality is desired.
function Base.transpose(P::DimArray{T}) where T <: AbstractDimArray
    ddims = dims(P)
    rdims = dims(first(P))
    A = Matrix(P)
    return MultipliableDimArray( transpose(A), ddims, rdims)
end

function Base.adjoint(P::DimArray{T}) where T <: AbstractDimArray
    ddims = dims(P)
    rdims = dims(first(P))
    A = Matrix(P)
    return MultipliableDimArray( adjoint(A), ddims, rdims)
end

function Base.:*(A::DimArray{T1}, b::DimArray{T2}) where T1 <: AbstractDimArray where T2 <: Number
    Amat = Matrix(A) * vec(b)
    (Amat isa Number) && (Amat = [Amat])
    rdims = dims(first(A))
    return MultipliableDimArray(Amat, rdims)
end
function Base.:*(A::DimArray{T1}, B::DimArray{T2}) where T1 <: AbstractDimArray where T2 <: AbstractDimArray
    Amat = Matrix(A) * Matrix(B)
    (Amat isa Number) && (Amat = [Amat])
    ddims = dims(B)
    rdims = dims(first(A))
    return MultipliableDimArray(Amat, rdims, ddims)
end

function Base.:(\)(A::DimArray{T1}, b::DimArray{T2}) where T1<: AbstractDimArray where T2 <: Number
    Amat = Matrix(A) \ vec(b)
    (Amat isa Number) && (Amat = [Amat])
    return DimArray(reshape(Amat, size(dims(A))), dims(A))
end

function Base.:(\)(A::DimArray{T1}, B::DimArray{T2})  where T1 <: AbstractDimArray where T2 <: AbstractDimArray 
    Amat = Matrix(A) \ Matrix(B)
    (Amat isa Number) && (Amat = [Amat])
    ddims = dims(B)
    rdims = dims(A)
    return MultipliableDimArray(Amat, rdims, ddims)
end
# Unitful doesn't handle complex quantities
function Base.:(\)(A::Matrix{T}, B::Matrix{T2}) where T <: Complex where T2 <: Quantity
    if uniform(B)
        Bunit = unit(first(B))
        return (1/Bunit) .* (A \ ustrip.(B))
    else
        error("matrix left divide not handled for non-uniform matrices")
    end
end

# Unitful doesn't handle matrix left divide between Quantity and non-Quantity
function Base.:(\)(A::Matrix{T}, b::Vector{T2}) where T <: Quantity where T2 <: Number
    if uniform(A)
        Aunit = unit(first(A))
        return (1/Aunit) .* (ustrip.(A) \ b)
    else
        error("matrix left divide not handled for non-uniform matrices")
    end
end

"""
     function left divide

     Left divide of Multipliable Matrix.
     Reverse mapping from unitdomain to range.
     Is `exact` if input is exact.
"""
function Base.:(\)(A::AbstractDimMatrix,b::AbstractDimVector)
    DimensionalData.comparedims(first(dims(A)), first(dims(b)); val=true)
    return rebuild(A,parent(A)\parent(b),(last(dims(A)),)) 
end

"""
function matrix right divide

`A/B = ( B'\\A')'
"""
function Base.:(/)(A::DimArray{T1}, B::DimArray{T2})  where T1 <: AbstractDimArray where T2 <: AbstractDimArray 
    Amat = Matrix(A) / Matrix(B)
    (Amat isa Number) && (Amat = [Amat])
    ddims = dims(first(B))
    rdims = dims(first(A))
    return MultipliableDimArray(Amat, rdims, ddims)
end

# Unitful doesn't handle complex quantities
function Base.:(/)(A::Matrix{T}, B::Matrix{T2}) where T <: Quantity where T2 <: Complex
    if uniform(A)
        Aunit = unit(first(A))
        return Aunit.* (ustrip.(A) / B)
    else
        error("matrix right divide not handled for non-uniform matrices")
    end
end
function Base.:(/)(A::Matrix{T2}, B::Matrix{T}) where T <: Quantity where T2 <: Complex
    if uniform(B)
        Bunit = unit(first(B))
        return (1/Bunit) .* (A / ustrip.(B)) # needs to be tested for proper units
    else
        error("matrix right divide not handled for non-uniform matrices")
    end
end

# eigenstructure only exists if A is uniform
# should be a better way by reading type
uniform(A::DimArray{<:DimArray}) = uniform(Matrix(A))
function uniform(A::Matrix)
    ulist = unit.(A)
    return allequal(ulist)
end

allequal(x) = all(y -> y == first(x), x)

function DiagonalDimArray(v::AbstractVector{T},Pdims::Tuple) where T

    tmp = zeros(T,Pdims)
    typetmp = typeof(tmp)

    P = Array{typetmp}(undef,size(tmp))

    for i in eachindex(P)
        P[i] = zeros(T, Pdims)
        P[i][i] += v[i]
    end
    return DimArray(P,Pdims)
end

# reverse the order of DiagonalDimArray
function LinearAlgebra.diag(P::DimArray{T}) where T <: DimArray{T2} where T2 <: Number

    d = zeros(T2,length(P))
    for i in eachindex(P)
        d[i] = P[i][i]
    end
    return d
end

function LinearAlgebra.eigen(A::DimArray{<:DimArray})

    uniform(A) ? uA = unit(first(first(A))) : error("No eigendecomposition for a non-uniform matrix")
    A_matrix = MultipliableDimArrays.Matrix(A)
    F = eigen(ustrip.(A_matrix))

    eigen_dims = Eigenmode(1:size(A_matrix,2))
    model_dims = dims(A)

    values = MultipliableDimArray(uA * F.values, eigen_dims)
    μ = DiagonalDimArray(values, dims(values))

    vectors = MultipliableDimArray(F.vectors,
            model_dims, eigen_dims)    

    return μ, vectors
    # ideally, would return an Eigen factorization, in spirit like:
    #    return Eigen(QuantityArray(F.values, dimension(A)), F.vectors)
end

endomorphic(A::DimArray{T}) where T<: AbstractDimArray = isequal(dims(A), dims(first(A)))

function exp(A::DimArray{T}) where T<: AbstractDimArray
    # A must be endomorphic (check type signature someday)
    !endomorphic(A) && error("A must be endomorphic to be consistent with matrix exponential")
    eA = exp(Matrix(A)) # move upstream to MultipliableDimArrays eventually
    return MultipliableDimArray(eA,dims(A),dims(A)) # wrap with same labels and format as A
end

include("unitful_linear_algebra.jl")

end
