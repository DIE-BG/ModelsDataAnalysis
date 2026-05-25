using DrWatson
@quickactivate "ModelsDataAnalysis"

using CairoMakie
import Printf
using DataFramesMeta
using DependentBootstrap: optblocklength
using Statistics: mean, median, std

# Helper functions
ft1(x) = Printf.@sprintf("%.1f", x)
ft2(x) = Printf.@sprintf("%.2f", x)


## Load data
WITH_REMITTANCES = true  # Configuration: choose whether to include remittances in the analysis
DLA_DATA = load_dla_data(WITH_REMITTANCES)
# Set output directory based on dataset choice
PLOTSDIR = mkpath(plotsdir("dla_analysis", "optblocklength", WITH_REMITTANCES ? "with_rem_gdp" : "without_rem_gdp"))


# Matrix data
mat_dla_data = Matrix(DLA_DATA[:, 2:end]) |> disallowmissing

K = size(mat_dla_data, 2)
varnames = string.(propertynames(DLA_DATA))[(begin + 1):end]

## Analysis of DLA data

# Median optimal block length
# Default is bootmethod=:stationary
optblocklength(mat_dla_data, blocklengthmethod = :ppw2009)
# optblocklength(mat_dla_data, blocklengthmethod=:ppw2009, bootmethod=:stationary)
# optblocklength(mat_dla_data, blocklengthmethod=:ppw2009, bootmethod=:moving)

# Optimal block length for each variable
opt_block_sizes = map(col -> optblocklength(col, bootmethod = :stationary), eachcol(mat_dla_data))
median(opt_block_sizes)

## Makie figure
fig = Figure(size = (950, 550))
lt = Label(fig[1, 1], "Patton, Politis, and White's (2009) Optimal Block Length", fontsize = 18, font = :bold, tellwidth = false)
ax = Axis(fig[2, 1])

barplot!(
    ax, 1:K, opt_block_sizes,
    bar_labels = ft1.(opt_block_sizes),
)
ax.xticks = (1:K, varnames)
ax.xticklabelrotation = pi / 4
ylims!(ax, 0, 30)

mean_block_size = mean(opt_block_sizes)
median_block_size = median(opt_block_sizes)

colors = Makie.wong_colors(0.8)
hlines!(mean_block_size, linewidth = 2, linestyle = :dot, color = colors[3])
hlines!(median(opt_block_sizes), linewidth = 2, linestyle = :dash, color = colors[2])
text!(ax, 9.5, mean_block_size, text = "Average = $(ft1(mean_block_size))", fontsize = 13)
text!(ax, 9.5, median_block_size + 1, text = "Median = $(ft1(median_block_size))", fontsize = 13)

filename = savename("ppw2009_optimal_block_length", (; data = "dla"), "png")
# save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig

## Analysis of DLA data

mat_dla_data = DLA_DATA[(begin + 1):(end - 1), (begin + 1):end] |> disallowmissing! |> Matrix
dla_varnames = string.(propertynames(DLA_DATA))[(begin + 1):end]

# Optimal block length for each variable
opt_block_sizes = map(col -> optblocklength(col), eachcol(mat_dla_data))

## Makie figure
fig = Figure(size = (950, 550))
lt = Label(fig[1, 1], "Patton, Politis, and White's (2009) Optimal Block Length", fontsize = 18, font = :bold, tellwidth = false)
ax = Axis(fig[2, 1])

barplot!(
    ax, 1:K, opt_block_sizes,
    bar_labels = ft1.(opt_block_sizes),
)
ax.xticks = (1:K, dla_varnames)
ax.xticklabelrotation = pi / 4
ylims!(ax, 0, 32)

mean_block_size = mean(opt_block_sizes)
median_block_size = median(opt_block_sizes)

colors = Makie.wong_colors(0.8)
hlines!(mean_block_size, linewidth = 2, linestyle = :dot, color = colors[3])
hlines!(median(opt_block_sizes), linewidth = 2, linestyle = :dash, color = colors[2])
text!(ax, 0.5, mean_block_size, text = "Average = $(ft1(mean_block_size))", fontsize = 13)
text!(ax, 0.5, median_block_size, text = "Median = $(ft1(median_block_size))", fontsize = 13)

filename = savename("ppw2009_optimal_block_length", (; data = "dla"), "png")
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit = 2.0)
fig
