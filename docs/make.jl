using Feynman
using Documenter

DocMeta.setdocmeta!(Feynman, :DocTestSetup, :(using Feynman); recursive=true)

makedocs(;
    modules=[Feynman],
    authors="Dushan Priyasad, Janko Boehm",
    sitename="Feynman.jl",
    format=Documenter.HTML(;
        canonical="https://honnatht.github.io/Feynman.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/honnatht/Feynman.jl",
    devbranch="main",
)
