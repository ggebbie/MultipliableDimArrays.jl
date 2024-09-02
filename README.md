# MultipliableDimArrays

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/MultipliableDimArrays.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/MultipliableDimArrays.jl/dev/)
[![Build Status](https://github.com/ggebbie/MultipliableDimArrays.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/MultipliableDimArrays.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/MultipliableDimArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/MultipliableDimArrays.jl)

# The issue

Box models and other gridded datasets can be stored in concise arrays that are organized according to physical space or some other organizational system. Geo-scientists are familiar with storing output data to file with meta-data in a self-describing format. For example, N-dimensional gridded data is naturally stored in N-dimensional arrays, but N-dimensional arrays are not in the right format to perform linear algebra operations with matrices and vectors. 

# A proposed solution

The Julia package `MultipliableDimArrays.jl` aims to do the right thing to permit linear algebra operations and then returns the output in the same human-readable output that the investigator originally provided. Cryptic references to boxes by a given sequential number are minimized in the code, and instead values can be looked up from more easily-interpretable names. Here we build upon `DimensionalData.jl` to do the heavy lifting of keeping data organized during the computational phase of the work.

# Implementation

An N-dimensional state vector, `x`, is stored as a `DimArray` where the dimensions of the array, i.e., `dims(x)` give the spatial or temporal locations.  These dimensions could also be used to denote which variable type is being referred to. For an array `x` that is meant to be handled as a vector during linear algebra operation, the Julia statement `A isa DimArray{T} where T <: Number` is `true`. The word, dimensions, unfortunately, has many different meanings. Do not confuse these dimensions with the units or physical quantities of the numerical entries in the array. The book, "Multidimensional Analysis" by George Hart deals with these physical units, instead. One may still input numerical values with units in the parent array of the `DimArray` by using the `Unitful` or `DynamicQuantities` Julia packages.a

A matrix is stored as a nested `DimArray`, that is a `DimArray` with element type `DimArray`. For a matrix `A`, the Julia statement `A isa DimArray{T} where T <: AbstractDimArray` is `true`. In accordance with the CR decomposition of a matrix, here we choose to store the outer dimension as the columns of the matrix (i.e., the range space), and the inner dimensions correspond to the rows of the matrix (i.e., the domain space). With this choice, a word of caution is necessary. Accessing the ith row and jth column of `A` cannot be done with `A[i,j]`. Instead one must use nested indices that are in the reverse order of typical mathematical notation `A[j][i]`.

# Type piracy

This package extends linear algebra operations for `DimArray`s. This package owns neither the linear algebra methods nor the `DimArray` struct. Thus the situation is one of type piracy, and should be avoided. Suggestions or discussions to handle this situation are enthusiastically sought.
