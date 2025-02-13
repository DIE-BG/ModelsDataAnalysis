using DrWatson
@quickactivate "ModelsDataAnalysis"

using Dates
using CairoMakie
using DataFrames
using DataFramesMeta
import Printf
using DependentBootstrap: dbootinds
using Statistics: mean, cov, cor
using LinearAlgebra: tril

# Helper functions
ft1(x) = Printf.@sprintf("%.1f", x)
ft2(x) = Printf.@sprintf("%.2f", x)

## Setup directories
PLOTSDIR = mkpath(plotsdir("covariance"))

## Load data 
const d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")
const dla_data = load(datadir("data_QPM.jld2"), "dla_data")
disallowmissing!(d4l_data)

# Matrix data
const mat_d4l_data = Matrix(d4l_data[:, 2:end]) |> disallowmissing 

const K = size(mat_d4l_data, 2)
const B = 10_000
const L_block = 40

# Observed correlation between covariates
const obs_cor = cor(mat_d4l_data)

# Sample using block bootstrap method

# Function to fill array of resamples
function resample_stats!(bootstrap_samples_cor, method=:stationary, l=10)
    # Compute indices for block bootstrap
    inds = dbootinds(mat_d4l_data, bootmethod=method, blocklength=l, numresample=B)
    # Perform the block bootstrap with length l and compute the correlation matrix
    for b in 1:B 
        bootstrap_sample = mat_d4l_data[inds[b], :]
        bootstrap_samples_cor[:, :, l, b] = cor(bootstrap_sample)
    end
end

# Array of (variable, variable, block length, bootstrap realization)
sbb_bootstrap_samples_cor = Array{eltype(mat_d4l_data)}(undef, K, K, L_block, B)
mbb_bootstrap_samples_cor = Array{eltype(mat_d4l_data)}(undef, K, K, L_block, B)

# Fill array of statistics with resamples for each method
map(l -> resample_stats!(sbb_bootstrap_samples_cor, :stationary, l), 1:L_block) 
map(l -> resample_stats!(mbb_bootstrap_samples_cor, :moving, l), 1:L_block) 

## Prototype a distribution figure 

l = 5

# Generate a figure for each block length
map(1:L_block) do l 

    fig = Figure(size=(1050, 1050))
    lt  = Label(fig[1,1:K+1], "Stationary Block Bootstrap correlation distribution between covariates", fontsize=18, font=:bold, tellwidth=false)
    st  = Label(fig[2,1:K+1], L"\ell=%$(l)", fontsize=18, tellwidth=false)

    for j in 1:K 
        Label(fig[3, j+1], varnames[j], tellwidth=false, fontsize=10)
    end
    for i in 1:K 
        Label(fig[i+3, 1], varnames[i], tellheight=false, rotation=pi/2, fontsize=10)
    end
    for i in 1:K, j in 1:K
        if i >= j 
            ax = Axis(fig[i+3, j+1])
            cor_dist = @view sbb_bootstrap_samples_cor[i, j, l, :] 
            if i == j 
                vlines!(ax, cor_dist[1], linewidth=4)
                ax.xticks = 0.:0.2:1.0
            else
                hist!(ax, cor_dist, normalization=:probability)
                # density!(ax, cor_dist, bandwidth=0.05)
                ax.xticks = -1.0:0.2:1.0
            end
            vlines!(obs_cor[i, j], linewidth=2, color=:red)
            ylims!(ax, 0, 0.3)
            ax.yticks = 0.0:0.1:0.3
            ax.xticklabelrotation = pi/4
            ax.xticklabelsize=10
        end

    end

    fig
    filename = savename("cor", (l=l,), "png", sort=false)
    save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)

end;


## MSE analysis

sqr(x) = x^2

# Eval mask to obtain lower triangular elements of the correlation matrices
eval_mask = tril(repeat([true], 10, 10), -1)

# Plot the overall MSE of the correlation matrix vs block length

function mse_cor(bootstrap_samples_cor, l)
    # Get bootstrap samples for each block length l
    bootstrap_samples = @view bootstrap_samples_cor[:, :, l, :]
    # Compute MSE
    mse = sum(mean(sqr, eval_mask .* (bootstrap_samples .- obs_cor), dims=3))::Float64
    mse
end

sbb_mse_cor = map(l -> mse_cor(sbb_bootstrap_samples_cor, l), 1:L_block)
mbb_mse_cor = map(l -> mse_cor(mbb_bootstrap_samples_cor, l), 1:L_block)

fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Block Bootstrap MSE of lower triangular covariance matrix", fontsize=18, font=:bold, tellwidth=false)
st  = Label(fig[2,1], "", fontsize=18, tellwidth=false)
ax  = Axis(fig[3,1], xlabel="Block length", ylabel="Mean Squared Error")

lines!(ax, 1:L_block, sbb_mse_cor, linewidth=2, label="Stationary")
lines!(ax, 1:L_block, mbb_mse_cor, linewidth=2, label="Moving")
axislegend(framevisible=false)
filename = savename("cor_mse", (method="stationary_moving",), "png", sort=false)
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig
