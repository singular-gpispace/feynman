#=
To build the documentation, do this from within the GromovWitten.jl directory:

julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl

To run the tests again, it suffices to repeat the second command.
=#

# Import required packages
using feynman
using Documenter

# Set up documentation test setup
DocMeta.setdocmeta!(feynman, :DocTestSetup, :(using feynman; recursive=true))

makedocs(;
    modules=[feynman],
    clean = true,
    checkdocs = :none,
    doctest = true,
    #strict = true,
    sitename="feynman.jl",
    expandfirst = ["Overview.md"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/feynman.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" =>[
            "index.md",
            "Installation.md",
        ],
        
        "Examples" =>[
            
        ],
        "Functions" => "Overview.md",
        
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/feynman.git",
    devbranch="main",
)
