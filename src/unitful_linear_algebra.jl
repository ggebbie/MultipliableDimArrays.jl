
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
function Matrix(P::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: DimArray

    # turn the data into a Matrix
    D = MultipliableDimArrays.Matrix(parent(P))

    # vectorize the unit dimensions
    # make a UnitfulMatrix without axis label Dimensions
    P2 = UnitfulLinearAlgebra.UnitfulMatrix(D,vec(unitrange(P)), vec(unitdomain(P)))

    # Turn into a Matrix that can be used for algebraic operations
    return UnitfulLinearAlgebra.Matrix(P2)
end

UnitfulLinearAlgebra.unitdomain(A::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: DimArray = last(unitdims(A))
#UnitfulLinearAlgebra.unitdomain(A::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: Number =  Unitful.Units([Unitful.NoUnits]) # kludge for a nondimensional scalar 

function matrix_fancy_print(P::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: DimArray
    # unwrap into algebraic form
    P2 = MultipliableDimArrays.Matrix(P)

    # re-label with axis Dimensions
    ddims = dims(parent(P))
    rdims = dims(first(parent(P)))

    return MultipliableDimArray(P2, rdims, ddims)

end

function Base.transpose(P::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: DimArray
    # unwrap into algebraic form
    P2 = MultipliableDimArrays.Matrix(P)

    # re-label with axis Dimensions
    ddims = dims(parent(P))
    rdims = dims(first(parent(P)))

    return MultipliableDimArray(transpose(P2), ddims, rdims)
end

function Base.:*(A::UnitfulMatrix{T1}, b::DimArray{T2}) where T1 <: DimArray where T2 <: Number

    # make a lowered UnitfulMatrix version
    # lowered = drops label Dimensions
    parent_lower = MultipliableDimArrays.Matrix(parent(A))
    range_lower = vec(unitrange(A))
    domain_lower = vec(unitdomain(A))
    exact_lower = exact(A)
    Alower = UnitfulMatrix(parent_lower,
        range_lower,domain_lower,exact=exact_lower)
    
    Amat = Alower * vec(b)
    (Amat isa Number) && (Amat = [Amat])
    rdims = dims(first(parent(A)))
    return MultipliableDimArray(Amat, rdims)
end
# function Base.:*(A::DimArray{T1}, B::DimArray{T2}) where T1 <: AbstractDimArray where T2 <: AbstractDimArray
#     Amat = Matrix(A) * Matrix(B)
#     (Amat isa Number) && (Amat = [Amat])
#     ddims = dims(B)
#     rdims = dims(first(A))
#     return MultipliableDimArray(Amat, rdims, ddims)
# end

Base.show(io::IO, mime::MIME"text/plain", A::UnitfulLinearAlgebra.UnitfulMatrix{T}) where T <: DimArray = Base.show(io::IO, mime::MIME"text/plain", matrix_fancy_print(A)) 
