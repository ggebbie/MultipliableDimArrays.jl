using Revise
using MultipliableDimArrays
using DimensionalData
using DimensionalData:@dim
using Test

# fixed paramete
@dim YearCE "years Common Era"
@dim SurfaceRegion "surface location"
@dim InteriorLocation "interior location"
@dim StateVariable "state variable"
surfaceregions = [:NATL,:ANT,:SUBANT]
N = length(surfaceregions)
years = (1990:2000)
statevariables = [:θ, :δ¹⁸O] 
M = 5 # Interior Locations with obs

function source_water_solution(surfaceregions,years)
    m = length(years)
    n = length(surfaceregions)
    x = DimArray(randn(m,n),(Ti(years),SurfaceRegion(surfaceregions)))
    return x
end

function source_water_solution(surfaceregions, years, statevar)
    m = length(years)
    n = length(surfaceregions)
    mat = cat(randn(m, n, 1), randn(m, n, 1); dims = 3)
    x = DimArray(mat, (Ti(years), SurfaceRegion(surfaceregions), StateVariable(statevar)))
    return x
end

@testset "MultipliableDimArrays.jl" begin

    x = source_water_solution(surfaceregions,
        years,
        statevariables)

    # test that these vectors;matrices can be used in algebraic expressions
    xvec = vec(x)
    x2 = vector_to_dimarray(xvec, dims(x))
    x3 = MultipliableDimArray(xvec, dims(x))
    @test x == x2
    @test x == x3
    
    yPmat = BLUEs.algebraic_object(y.P)
    yPda = BLUEs.matrix_to_dimarray(yPmat, dims(y.P), dims(y.P))
    @test y.P == yPda



end
