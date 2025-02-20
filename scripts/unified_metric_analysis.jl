using DrWatson
@quickactivate "ModelsDataAnalysis"

using Dates
using CairoMakie
using DataFrames
using DataFramesMeta
import Printf
using DependentBootstrap: dbootinds
using Statistics: mean, std, var, cor
using StatsBase: autocor, autocor!
using LinearAlgebra: tril

# Helper functions
ft1(x) = Printf.@sprintf("%.1f", x)
ft2(x) = Printf.@sprintf("%.2f", x)
sqr(x) = x^2

## Setup directories
PLOTSDIR = mkpath(plotsdir("unified"))

## Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")
dla_data = load(datadir("data_QPM.jld2"), "dla_data")
disallowmissing!(d4l_data)

# Matrix data
const mat_d4l_data = Matrix(d4l_data[:, 2:end]) |> disallowmissing 

const K = size(mat_d4l_data, 2)
const B = 10_000
const L_autocor = 12
const L_block   = 40

# Historically observed statistics
const obs_mean      = mean(mat_d4l_data, dims=1)
const obs_var       = var(mat_d4l_data, dims=1)
const obs_cor       = cor(mat_d4l_data)
const obs_autocor   = autocor(mat_d4l_data, 1:L_autocor)

# Sample using block bootstrap method

TT = eltype(mat_d4l_data)
# ACF array (lags of ACF, variable, block length, bootstrap realization)
sbb_samples_autocor = Array{TT}(undef, L_autocor, K, L_block, B)
mbb_samples_autocor = Array{TT}(undef, L_autocor, K, L_block, B)
# Correlation array (variable, variable, block length, bootstrap realization)
sbb_samples_cor = Array{TT}(undef, K, K, L_block, B)
mbb_samples_cor = Array{TT}(undef, K, K, L_block, B)
# Mean array (1, variable, block length, bootstrap realization)
sbb_samples_mean = Array{TT}(undef, 1, K, L_block, B)
mbb_samples_mean = Array{TT}(undef, 1, K, L_block, B)
# Variance array (1, variable, block length, bootstrap realization)
sbb_samples_var = Array{TT}(undef, 1, K, L_block, B)
mbb_samples_var = Array{TT}(undef, 1, K, L_block, B)

sbb_dict = @dict mean=sbb_samples_mean var=sbb_samples_var autocor=sbb_samples_autocor cor=sbb_samples_cor
mbb_dict = @dict mean=mbb_samples_mean var=mbb_samples_var autocor=mbb_samples_autocor cor=mbb_samples_cor

# Function to fill array of resamples
function resample_stats!(dict_bb_samples, method=:stationary, l=10)

    mean_array = dict_bb_samples[:mean]
    var_array = dict_bb_samples[:var]
    autocor_array = dict_bb_samples[:autocor]
    cor_array = dict_bb_samples[:cor]

    # Compute indices for block bootstrap
    inds = dbootinds(mat_d4l_data, bootmethod=method, blocklength=l, numresample=B)
    # Perform the block bootstrap with length l
    @info "Performing bootstrap resampling for method $method and block length $l"
    for b in 1:B 
        # Get bootstrap sample
        bootstrap_sample = mat_d4l_data[inds[b], :]

        # Mean 
        mean_array[:, :, l, b] = mean(bootstrap_sample, dims=1)
        # Variance
        var_array[:, :, l, b]  = var(bootstrap_sample, dims=1)
        # Correlation
        cor_array[:, :, l, b] = cor(bootstrap_sample)
        # Autocorrelation 
        autocor_mat = @view autocor_array[:, :, l, b]
        autocor!(autocor_mat, bootstrap_sample, 1:L_autocor)
    end
end

# Stationary block
map(l -> resample_stats!(sbb_dict, :stationary, l), 1:L_block);
# Moving block
map(l -> resample_stats!(mbb_dict, :moving, l), 1:L_block);


## Compute unified metric per block length 

# Aggregates the ACF lags of each variable
function autocor_agg_exp_fn(mse_autocor_lags_variable, alpha=0.9)
    # Weighted average of the MSE of each of the ACF lags
    w = [alpha^i for i in 0:L_autocor-1]
    # w = ones(L_autocor)
    w = w / sum(w)
    mse_autocor_variable = sum(mse_autocor_lags_variable .* w, dims=1) 
    mse_autocor_variable
end

function unified_metric(
    dict_bb_samples::Dict, l::Int; 
    variable_agg_fn=mean, 
    autocor_agg_fn=autocor_agg_exp_fn, 
    cor_agg_fn=(m -> 2 * sum(m) / (K * (K-1)))
    ) 

    # Get the arrays for block length l
    mean_array      = @view dict_bb_samples[:mean][:, :, l, :]
    var_array       = @view dict_bb_samples[:var][:, :, l, :]
    autocor_array  = @view dict_bb_samples[:autocor][:, :, l, :]
    cor_array       = @view dict_bb_samples[:cor][:, :, l, :]

    # Metric for the mean 
    std_mean = std(mean_array, dims=3)
    # Normalized MSE for the mean estimator each variable
    mse_mean_variable = mean(sqr, ((mean_array .- obs_mean ) ./ std_mean), dims=3)
    # Normalized MSE of method 
    mse_mean = variable_agg_fn(mse_mean_variable)

    # Metric for the variance 
    std_var = std(var_array, dims=3)
    # Normalized MSE for the variance estimator each variable
    mse_var_variable = mean(sqr, ((var_array .- obs_var) ./ std_var), dims=3)
    # Normalized MSE of method 
    mse_var = variable_agg_fn(mse_var_variable)

    # Metric for autocorrelation 
    std_autocor = std(autocor_array, dims=3)
    # Normalized MSE for the autocorrelation estimator of each variable, from 1 to L_autocor lags
    mse_autocor_lags_variable = mean(sqr, ((autocor_array .- obs_autocor) ./ std_autocor), dims=3)
    # Normalized MSE for the autocorrelation estimator of each variable
    mse_autocor_variable = autocor_agg_fn(mse_autocor_lags_variable)
    # Normalized MSE of method 
    mse_autocor = variable_agg_fn(mse_autocor_variable)

    # Metric for correlation
    std_cor = std(cor_array, dims=3)
    # Evaluation mask to obtain lower triangular elements of the correlation matrices
    eval_mask = tril(repeat([true], K, K), -1)
    # Normalized MSE for the correlation estimator of each variable, from 1 to L_autocor lags
    mse_cor_matrix = mean(sqr, eval_mask .* ((cor_array .- obs_cor) ./ std_cor), dims=3)
    # Normalized MSE of method 
    mse_cor = cor_agg_fn(mse_cor_matrix)

    mse_mean, mse_var, mse_autocor, mse_cor
end

# tt = unified_metric(mbb_dict, 10)

# Compute statistics for each block length
sbb_stats_metrics = zeros(L_block, 4)
for l in 1:L_block 
    sbb_stats_metrics[l, :] .= unified_metric(sbb_dict, l)
end

mbb_stats_metrics = zeros(L_block, 4)
for l in 1:L_block 
    mbb_stats_metrics[l, :] .= unified_metric(mbb_dict, l)
end


## Plot unified metrics 

fig = Figure(size=(950,550))
lt  = Label(fig[1,1:2], "Components of the Normalized MSE for Block Bootstrap Methods", fontsize=18, font=:bold, tellwidth=false)
ax1  = Axis(fig[2,1], title="Stationary", xlabel="Block length")
ax2  = Axis(fig[2,2], title="Moving", xlabel="Block length")

lines!(ax1, 1:L_block, sbb_stats_metrics[:, 1], label="Mean")
lines!(ax1, 1:L_block, sbb_stats_metrics[:, 2], label="Variance")
lines!(ax1, 1:L_block, sbb_stats_metrics[:, 3], linewidth=2, label="Autocorrelation")
lines!(ax1, 1:L_block, sbb_stats_metrics[:, 4], linewidth=3, linestyle=:dash, label="Correlation")
axislegend(ax1)

lines!(ax2, 1:L_block, mbb_stats_metrics[:, 1], label="Mean")
lines!(ax2, 1:L_block, mbb_stats_metrics[:, 2], label="Variance")
lines!(ax2, 1:L_block, mbb_stats_metrics[:, 3], linewidth=2, label="Autocorrelation")
lines!(ax2, 1:L_block, mbb_stats_metrics[:, 4], linewidth=3, linestyle=:dash, label="Correlation")
axislegend(ax2)

# Synchronize the axes
linkyaxes!(ax1, ax2)
# hideydecorations!(ax2)

filename = savename("unified_metrics_components", (method="stationary_moving",), "png")
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig


## Energy analysis, block size that corresponds to beta % of the decrease in the MSE

beta = 0.95

sbb_mse = sum(sbb_stats_metrics, dims=2) |> vec
sbb_min_mse, sbb_max_mse = extrema(sbb_mse)
threshold_sbb_mse = sbb_min_mse + (1-beta) * (sbb_max_mse - sbb_min_mse)

sbb_l = findfirst(<=(threshold_sbb_mse), sbb_mse)
@info "SBB block size for β=$beta is L=$sbb_l"

mbb_mse = sum(mbb_stats_metrics, dims=2) |> vec
mbb_min_mse, mbb_max_mse = extrema(mbb_mse)
threshold_mbb_mse = mbb_min_mse + (1-beta) * (mbb_max_mse - mbb_min_mse)

mbb_l = findfirst(<=(threshold_mbb_mse), mbb_mse)
@info "MBB block size for β=$beta is L=$mbb_l"

# [ Info: SBB block size for β=0.9 is L=6
# [ Info: SBB block size for β=0.95 is L=10
# [ Info: SBB block size for β=0.99 is L=24
# [ Info: MBB block size for β=0.9 is L=6
# [ Info: MBB block size for β=0.95 is L=10
# [ Info: MBB block size for β=0.99 is L=24

mbb_l_star = findfirst(<=(mbb_min_mse), mbb_mse)

## Compare the sum of the metrics

fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Unified Normalized MSE for Block Bootstrap Methods", fontsize=18, font=:bold, tellwidth=false)
ax  = Axis(fig[2,1], xlabel="Block length")

lines!(ax, 1:L_block, sbb_mse, label="Stationary", linewidth=2)
lines!(ax, 1:L_block, mbb_mse, label="Moving", linewidth=2)
hlines!(ax, sbb_mse[sbb_l], color=(:gray, 0.6), linestyle=:dash)
# text!(ax, 0, sbb_mse[sbb_l], text="SBB\n$(ft1(100*beta))% of total decrease", fontsize=12)
scatter!(ax, [sbb_l], [sbb_mse[sbb_l]], 
    label="SBB $(ft1(100*beta))% of total decrease", 
    markersize=14, 
    color=:blue,
)
scatter!(ax, [mbb_l_star], [mbb_min_mse], 
    label="MBB minimum",
    markersize=14, 
    color=:orange
)
axislegend(ax)

ylims!(ax, 0, 27)
ax.yticks = 0:2:28
ax.xticks = 0:2:L_block

filename = savename("unified_metrics", (method="stationary_moving",), "png")
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig

