module MultipliableDimArrays

using DimensionalData
using LinearAlgebra

export algebraic_object
export MultipliableDimArray
export matrix_to_dimarray
    
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

"""
function algebraic_object(P::DimArray{T}) where T <: AbstractDimArray
"""
function algebraic_object(P::DimArray{T}) where T <: AbstractDimArray
    # number of columns/ outer dims
    N = length(P)
    # number of rows, take first inner element as example
    M = length(first(P))
    A = Matrix{eltype(first(P))}(undef,M,N)
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
function matrix_to_dimarray(A,rangedims,domaindims)
"""
function matrix_to_dimarray(A,rangedims,domaindims)

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

function MultipliableDimArray(A::AbstractVector,
                              rdims)
    Q1 = reshape(A,size(rdims))
    return DimArray(Q1, rdims)
end

function algebraic_transpose(P::DimArray)
    ddims = dims(P)
    rdims = dims(first(P))
    A = algebraic_object(P)
    return matrix_to_dimarray( transpose(A), ddims, rdims)
end

function ldiv(A::DimArray, b::DimArray) 
    Amat = algebraic_object(A) \ algebraic_object(b)
    (Amat isa Number) && (Amat = [Amat])
    ddims = dims(b)
    rdims = dims(A)
    return matrix_to_dimarray(Amat, rdims, ddims)
#    return DimArray(reshape(Amat, size(dims(A))), dims(A))
end
function ldiv(A::DimArray{T1}, b::DimArray{T2}) where T1<: AbstractDimArray where T2 <: Number
    Amat = algebraic_object(A) \ algebraic_object(b)
    (Amat isa Number) && (Amat = [Amat])
    return DimArray(reshape(Amat, size(dims(A))), dims(A))
end

function matmul(A::DimArray, b::DimArray{T}) where T <: Number
    Amat = algebraic_object(A) * algebraic_object(b)
    (Amat isa Number) && (Amat = [Amat])
    rdims = dims(first(A))
    return vector_to_dimarray(Amat, rdims)
  # return DimArray( reshape(Amat, size(rdims)), rdims)
end
function matmul(A::DimArray, b::DimArray)
    Amat = algebraic_object(A) * algebraic_object(b)
    (Amat isa Number) && (Amat = [Amat])
#    rdims = dims(first(A))
 #   return DimArray( reshape(Amat, size(rdims)), rdims)
    ddims = dims(b)
    rdims = dims(first(A))
    return matrix_to_dimarray(Amat, rdims, ddims)

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

end
