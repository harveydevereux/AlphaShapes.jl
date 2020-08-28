using Documenter, AlphaShapes

makedocs(sitename="AlphaShapes.jl",format = Documenter.HTML(prettyurls = true))

deploydocs(
    repo = "github.com/harveydevereux/AlphaShapes.jl.git"
)
