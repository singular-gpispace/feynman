using Feynman
using Documenter

DocMeta.setdocmeta!(Feynman, :DocTestSetup, :(using Feynman); recursive=true)
# Run Singular and capture the output
example_output = readchomp(`Singular run_examples.sing`)

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
            "index.md",
            "Installation.md",
        ],
        "Examples" => [
            "Example.md",
        ],
        "Functions" => [
            "Overview.md",
        ],
        "Singular version of Feynman" => [
            "feynman_lib.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/singular-gpispace/Feynman.jl",
    devbranch="main",
    branch="gh-pages",

)

