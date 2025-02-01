##-----direct to the julia working enviorenment
#cd /hoome/dushan/julia_test

#JULIA_PROJECT=docs/ julia -e 'using Pkg; Pkg.develop(path="/home/dushan/FM/feynman"); Pkg.instantiate()'
#JULIA_PROJECT=docs/ julia /home/dushan/FM/feynman/docs/make.jl
using Feynman
using Documenter

DocMeta.setdocmeta!(Feynman, :DocTestSetup, :(using Feynman); recursive=true)

makedocs(;
    modules=[Feynman],
    clean = true,
    checkdocs = :none,
    doctest = true,
    authors = "Dushan Priyasad, Janko Boehm",
    sitename = "Feynman.jl",
    expandfirst = ["Overview.md"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/Feynman.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => [
            "src/index.md",
            "src/Installation.md",
            "src/Example.md",
         #   "SmallExample.md",
        ],

        "Examples" =>[
            "src/Example.md",

        ],

        "Functions" =>[ 
            "src/Overview.md",
        ],

    ],


)

deploydocs(;
    repo="github.com/singular-gpispace/Feynman.jl",
    devbranch="main",
    branch="gh-pages",

)
