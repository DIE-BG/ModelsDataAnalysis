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

## Setup directories
PLOTSDIR = mkpath(plotsdir("optblocklength"))

## Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")
dla_data = load(datadir("data_QPM.jld2"), "dla_data")
disallowmissing!(d4l_data)

# Matrix data
const mat_d4l_data = Matrix(d4l_data[:, 2:end]) |> disallowmissing 
const K = size(mat_d4l_data, 2)
const varnames = string.(propertynames(d4l_data))[begin+1:end]

## Analysis of D4L data

# Median optimal block length
optblocklength(mat_d4l_data, blocklengthmethod=:ppw2009)

# Optimal block length for each variable
opt_block_sizes = map(col -> optblocklength(col), eachcol(mat_d4l_data))
median(opt_block_sizes)

# Makie figure 

fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Patton, Politis, and White's (2009) Optimal Block Length", fontsize=18, font=:bold, tellwidth=false)
ax  = Axis(fig[2,1])

barplot!(ax, 1:K, opt_block_sizes, 
    bar_labels = ft1.(opt_block_sizes), 
)
ax.xticks = (1:K, varnames)
ax.xticklabelrotation = pi/4
ylims!(ax, 0, 25)

mean_block_size = mean(opt_block_sizes)
median_block_size = median(opt_block_sizes)

colors = Makie.wong_colors(0.8)
hlines!(mean_block_size, linewidth=2, linestyle=:dot, color=colors[3])
hlines!(median(opt_block_sizes), linewidth=2, linestyle=:dash, color=colors[2])
text!(ax, 9.5, mean_block_size, text="Average = $(ft1(mean_block_size))", fontsize=13)
text!(ax, 9.5, median_block_size, text="Median = $(ft1(median_block_size))", fontsize=13)

filename = savename("ppw2009_optimal_block_length", (;data="d4l"), "png")
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig

## Analysis of DLA data

mat_dla_data = dla_data[begin+1:end-1, begin+1:end] |> disallowmissing! |> Matrix
const dla_varnames = string.(propertynames(dla_data))[begin+1:end]

# Optimal block length for each variable
opt_block_sizes = map(col -> optblocklength(col), eachcol(mat_dla_data))

# Makie figure 
fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Patton, Politis, and White's (2009) Optimal Block Length", fontsize=18, font=:bold, tellwidth=false)
ax  = Axis(fig[2,1])

barplot!(ax, 1:K, opt_block_sizes, 
    bar_labels = ft1.(opt_block_sizes), 
)
ax.xticks = (1:K, dla_varnames)
ax.xticklabelrotation = pi/4
ylims!(ax, 0, 32)

mean_block_size = mean(opt_block_sizes)
median_block_size = median(opt_block_sizes)

colors = Makie.wong_colors(0.8)
hlines!(mean_block_size, linewidth=2, linestyle=:dot, color=colors[3])
hlines!(median(opt_block_sizes), linewidth=2, linestyle=:dash, color=colors[2])
text!(ax, 0.5, mean_block_size, text="Average = $(ft1(mean_block_size))", fontsize=13)
text!(ax, 0.5, median_block_size, text="Median = $(ft1(median_block_size))", fontsize=13)

filename = savename("ppw2009_optimal_block_length", (; data="dla"), "png")
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig
