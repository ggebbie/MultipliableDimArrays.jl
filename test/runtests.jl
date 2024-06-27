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
    x3 = MultipliableDimArray(xvec, dims(x))
    @test x == x3

    v = ones(length(x))
    Px = DiagonalDimArray(v,dims(x))
    Pxmat = Matrix(Px)
    Px2 = MultipliableDimArray(Pxmat, dims(x), dims(x))
    @test Px == Px2

    PxT = transpose(Px)
    PxTT = transpose(PxT)
    @test Px == PxTT

    Rx = MultipliableDimArray(rand(length(x),length(x)),
        dims(x), dims(x))    
    RxT = transpose(Rx)
    RxTT = transpose(RxT)
    @test Rx == RxTT

    q = Rx * x
    x2 = Rx \ q
    @test isapprox(x, x2, atol = 1e-8)

    Sx = MultipliableDimArray(rand(length(x),length(x)),
        dims(x), dims(x))    
    Q = Rx * Sx
    Sx2 = Rx \ Q 
    @test isapprox(Sx2, Sx, atol = 1e-8)
    
end

