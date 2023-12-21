using Documenter
using DocumenterTools
using DocStringExtensions
push!(LOAD_PATH, "../src/")
using Gears

@show pkgPath = joinpath(dirname(pathof(Gears)), "..") # this includes Gears/src, so ..
makedocs(
  sitename="Gears.jl",
  modules=[Gears],
  root = joinpath(pkgPath, "docs"),
  source = "src",
  build = "build",
  clean=true,
  doctest=true,
  draft=false,
  checkdocs=:all,
  # linkcheck=true, fails to find internal links to bookmarks..
  )

# compile custom theme scss in to css, copying over the default themes
DocumenterTools.Themes.compile(joinpath(pkgPath,"docs","src","assets","themes","documenter-mechanomy.scss"), joinpath(pkgPath,"docs","build","assets","themes","documenter-dark.css") )
DocumenterTools.Themes.compile(joinpath(pkgPath,"docs","src","assets","themes","documenter-mechanomy.scss"), joinpath(pkgPath,"docs","build","assets","themes","documenter-light.css") )
deploydocs(
  root = joinpath(pkgPath, "docs"),
  target = "build",
  dirname = "",
  repo = "github.com/mechanomy/Gears.jl.git",
  branch = "gh-pages",
  deps = nothing, 
  make = nothing,
  devbranch = "main",
  devurl = "dev",
  versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
  forcepush = false,
  deploy_config = Documenter.auto_detect_deploy_system(),
  push_preview = false,
)
