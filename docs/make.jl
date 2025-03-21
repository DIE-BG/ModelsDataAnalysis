CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using DrWatson
@quickactivate "ModelsDataAnalysis"
using Documenter

# Here you may include files from the source directory
# include(srcdir("dummy_src_file.jl"))

# Copy plots folder 
# SRC_PLOTS_PATH = mkpath(projectdir("docs", "src", "images", "periodograms"))
# cp(plotsdir("periodograms"), SRC_PLOTS_PATH; force=true)

@info "Building Documentation"
makedocs(;
    sitename = "ModelsDataAnalysis",
    # This argument is only so that the sequence of pages in the sidebar is configured
    # By default all markdown files in `docs/src` are expanded and included.
    pages = [
        "index.md",
        "Frequency Analysis" => "periodogram.md",
        "Block Bootstrap" => [
            "Simulation Study" => "block_bootstrap_exercise.md"
            "Detailed Results" => [
                "mean.md",
                "std.md",
                "acf.md",
                "covariance.md"
            ]
        ]
    ],
    # Don't worry about what `CI` does in this line.
    format = Documenter.HTML(prettyurls = CI),
)

@info "Deploying Documentation"
if CI
    deploydocs(
        # `repo` MUST be set correctly. Once your GitHub name is set
        # the auto-generated documentation will be hosted at:
        # https://PutYourGitHubNameHere.github.io/ModelsDataAnalysis/dev/
        # (assuming you have enabled `gh-pages` deployment)
        repo = "github.com/DIE-BG/ModelsDataAnalysis.git",
        target = "build",
        push_preview = true,
        devbranch = "main",
    )
end

@info "Finished with Documentation"
