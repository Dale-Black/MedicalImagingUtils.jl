using MedicalImagingUtils
using Documenter

DocMeta.setdocmeta!(MedicalImagingUtils, :DocTestSetup, :(using MedicalImagingUtils); recursive=true)

makedocs(;
    modules=[MedicalImagingUtils],
    authors="Dale <djblack@uci.edu> and contributors",
    repo="https://github.com/Dale-Black/MedicalImagingUtils.jl/blob/{commit}{path}#{line}",
    sitename="MedicalImagingUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dale-Black.github.io/MedicalImagingUtils.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Dale-Black/MedicalImagingUtils.jl",
)
