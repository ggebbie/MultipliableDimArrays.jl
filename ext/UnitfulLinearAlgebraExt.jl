module UnitfulLinearAlgebraExt

using MultipliableDimArrays, UnitfulLinearAlgebra, LinearAlgebra, DimensionalData, Unitful

println("loading extension")

# extended constructor, just save it for now, not much type checking
function UnitfulLinearAlgebra.UnitfulMatrix(data::DimArray{T}, dims::Tuple; exact = false) where T <: DimArray
# should also check that dims(dims) == dims(data)
    return UnitfulLinearAlgebra.UnitfulMatrix(data,
        dims, exact)
#                DimensionalData.format(dims, data), exact) # formatting not working
end

"""
function Matrix(P::UnitfulMatrix{T}) where T <: AbstractDimArray
"""
function Matrix(P::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: AbstractDimArray

    # turn the data into a Matrix
    D = MultipliableDimArrays.Matrix(parent(P))

    # vectorize the unit dimensions
    # make a UnitfulMatrix without axis label Dimensions
    P2 = UnitfulLinearAlgebra.UnitfulMatrix(D,vec(unitrange(P)), vec(unitdomain(P)))

    # Turn into a Matrix that can be used for algebraic operations
    return Matrix(P2)
end

UnitfulLinearAlgebra.unitdomain(A::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: AbstractDimArray = last(unitdims(A))
#UnitfulLinearAlgebra.unitdomain(A::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: Number =  Unitful.Units([Unitful.NoUnits]) # kludge for a nondimensional scalar 

# Write your package code here.
# function LinearAlgebra.diag(A::AbstractUnitfulMatrix{T}) where T <: AbstractDimArray

#     if uniform(A)
#         inside_type = typeof(getindexqty(A,1,1))
#         sigma = Array{inside_type}(undef,size(A))
#     else
#         sigma = Array{Quantity}(undef,size(A))
#     end
    
#     for i in eachindex(A)
#         sigma[i] = getindexqty(A,i,i)
#     end
#     return DimArray(sigma,dims(A))
# end

end
