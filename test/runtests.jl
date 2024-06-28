using Revise
using MultipliableDimArrays
using Unitful
using UnitfulLinearAlgebra
using DimensionalData
using DimensionalData:@dim
using Test

# fixed parameters
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

    @testset "no units" begin
        x = source_water_solution(surfaceregions,
            years,
            statevariables)
            
        # test that these vectors;matrices can be used in algebraic expressions
        xvec = vec(x)
        x3 = MultipliableDimArray(xvec, dims(x))
        @test x == x3

        v = ones(length(x))
        Px = DiagonalDimArray(v,dims(x))
        # try to remove reference 
        Pxmat = MultipliableDimArrays.Matrix(Px)
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

    @testset "UnitfulLinearAlgebra extension" begin

        using Unitful
        
        function source_water_solution_with_units(surfaceregions, years, statevar)
            yr = u"yr"
            K = u"K"
            permil = u"permille"
            m = length(years)
            n = length(surfaceregions)
            mat = cat(randn(m, n, 1)K, randn(m, n, 1)permil; dims = 3)
            x = DimArray(mat, (Ti(years), SurfaceRegion(surfaceregions), StateVariable(statevar)))
            return x
        end

        x = source_water_solution_with_units(surfaceregions,
            years,
            statevariables)

        # test that these vectors;matrices can be used in algebraic expressions
        xvec = vec(x)
        x3 = MultipliableDimArray(xvec, dims(x))
        @test x == x3

#        using UnitfulLinearAlgebra

        D = MultipliableDimArray(randn(length(x),length(x)),dims(x),dims(x))
        urange = unit.(x)
        udomain = unit.(inv.(x))
        Px = UnitfulMatrix(D,(urange,udomain))

        # try to remove reference 
        Pxmat = MultipliableDimArrays.Matrix(Px)
        Px2 = MultipliableDimArray(Pxmat, dims(x), dims(x))
        @test Px[2][3] == ustrip(Px2[2][3]) # weak test, but function to revert operation is not completed, one is DimArray, one is UnitfulMatrix

        PxT = transpose(Px)
        PxTT = transpose(PxT)
#        @test Px == PxTT
        @test Px[2][3] == ustrip(PxTT[2][3]) # weak test, but function to revert operation is not completed, one is DimArr    
        # make quick inverse of Px for unit compatibility
        iPx = UnitfulMatrix(D,(udomain,urange))

        # matrix-vector multiplication
        q = iPx * x

        # x2 = iPx \ q
        # @test isapprox(x, x2, atol = 1e-8)

        Ialmost = iPx * Px # matrix-matrix multiplication, except iPx is not actually the inverse of Px
        MultipliableDimArrays.Matrix(Ialmost)  # visually reasonable
end
