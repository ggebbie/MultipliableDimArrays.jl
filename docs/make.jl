using MultipliableDimArrays
using Documenter

DocMeta.setdocmeta!(MultipliableDimArrays, :DocTestSetup, :(using MultipliableDimArrays); recursive=true)

makedocs(;
    modules=[MultipliableDimArrays],
    authors="G Jake Gebbie <ggebbie@whoi.edu>",
    sitename="MultipliableDimArrays.jl",
    format=Documenter.HTML(;
        canonical="https://ggebbie.github.io/MultipliableDimArrays.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/MultipliableDimArrays.jl",
    devbranch="main",
)
