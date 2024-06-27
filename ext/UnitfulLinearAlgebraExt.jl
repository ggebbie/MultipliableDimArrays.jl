module UnitfulLinearAlgebraExt

using MultipliableDimArrays, UnitfulLinearAlgebra, LinearAlgebra

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
