module MultipliableDimArrays

using DimensionalData
using LinearAlgebra
using Unitful
using UnitfulLinearAlgebra
export MultipliableDimArray
export DiagonalDimArray
export Matrix

import Base: Matrix

# an alias
#MultiDimArray{T} = DimArray{T} where T <: AbstractDimArray

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

function Base.transpose(P::DimArray{T}) where T <: AbstractDimArray
    ddims = dims(P)
    rdims = dims(first(P))
    # if T isa Number
    #     A = vec(P)
    # elseif T isa AbstractDimArray

    A = Matrix(P)
    #end
    return MultipliableDimArray( transpose(A), ddims, rdims)
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

function Base.:\(A::DimArray{T1}, b::DimArray{T2}) where T1<: AbstractDimArray where T2 <: Number
    Amat = Matrix(A) \ vec(b)
    (Amat isa Number) && (Amat = [Amat])
    return DimArray(reshape(Amat, size(dims(A))), dims(A))
end
function Base.:\(A::DimArray{T1}, B::DimArray{T2})  where T1 <: AbstractDimArray where T2 <: AbstractDimArray 
    Amat = Matrix(A) \ Matrix(B)
    (Amat isa Number) && (Amat = [Amat])
    ddims = dims(B)
    rdims = dims(A)
    return MultipliableDimArray(Amat, rdims, ddims)
end

"""
     function left divide

     Left divide of Multipliable Matrix.
     Reverse mapping from unitdomain to range.
     Is `exact` if input is exact.
"""
function Base.:\(A::AbstractDimMatrix,b::AbstractDimVector)
    DimensionalData.comparedims(first(dims(A)), first(dims(b)); val=true)
    return rebuild(A,parent(A)\parent(b),(last(dims(A)),)) 
end

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

# vec works just as well (maybe an issue when units appear)
# """
# function algebraic_object(P::DimArray{Number})
# """
# function algebraic_object(P::DimArray{T}) where T <: Number
#     M = length(P)
#     A = Vector{T}(undef,M)
#     for i in eachindex(P)
#         A[i] = P[i]
#     end
#     return A 
# end

include("unitful_linear_algebra.jl")

end
