using DrWatson
@quickactivate "ModelsDataAnalysis"

using Dates
using CairoMakie
using DataFrames
using DataFramesMeta
import Printf
using DependentBootstrap: dbootinds
using Statistics: mean
using StatsBase: autocor, autocor!

# Helper functions
ft1(x) = Printf.@sprintf("%.1f", x)
ft2(x) = Printf.@sprintf("%.2f", x)


## Setup directories
PLOTSDIR = mkpath(plotsdir("autocorrelation"))

## Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")
dla_data = load(datadir("data_QPM.jld2"), "dla_data")
disallowmissing!(d4l_data)

# Matrix data
mat_d4l_data = Matrix(d4l_data[:, 2:end]) |> disallowmissing 

K = size(mat_d4l_data, 2)
B = 10_000
L_autocor = 12
L_block   = 40

# Observed autocorrelation 
obs_autocor = autocor(mat_d4l_data, 0:L_autocor)

# Sample using block bootstrap method

# Array of (lags of ACF, variable, block length, bootstrap realization)
sbb_bootstrap_samples_autocor = Array{eltype(mat_d4l_data)}(undef, L_autocor+1, K, L_block, B)
mbb_bootstrap_samples_autocor = Array{eltype(mat_d4l_data)}(undef, L_autocor+1, K, L_block, B)

# Function to fill array of resamples
function resample_stats!(bootstrap_samples_autocor, method=:stationary, l=10)
    # Compute indices for block bootstrap
    inds = dbootinds(mat_d4l_data, bootmethod=method, blocklength=l, numresample=B)
    # Perform the block bootstrap with length l
    for b in 1:B 
        bootstrap_sample = mat_d4l_data[inds[b], :]
        output_mat = @view bootstrap_samples_autocor[:, :, l, b]
        autocor!(output_mat, bootstrap_sample, 0:L_autocor)
    end
end

# Stationary block
map(1:L_block) do l
    resample_stats!(sbb_bootstrap_samples_autocor, :stationary, l)
end; 

# Moving block
map(1:L_block) do l
    resample_stats!(mbb_bootstrap_samples_autocor, :moving, l)
end; 


## Prototype a distribution figure 

varnames = string.(propertynames(d4l_data))[begin+1:end]
k = 1
l = 10


function plot_acf_distribution(k, l, bootstrap_samples_autocor=bootstrap_samples_autocor)

    fig = Figure(size=(950,550))
    title = varnames[k] * ", blocklength=$l"
    lt  = Label(fig[1,1], varnames[k], fontsize=18, font=:bold, tellwidth=false)
    st  = Label(fig[2,1], L"\ell=%$(l)", fontsize=18, tellwidth=false)
    ax  = Axis(fig[3,1]) #, xlabel="Lag", ylabel="Autocorrelation")

    # for j in 1:L_autocor+1
    #     density!(
    #         ax,
    #         bootstrap_samples_autocor[j, k, l, :],
    #         direction   = :y, 
    #         offset  = j - 1, 
    #         color   = (:slategray, 0.4), 
    #         bandwidth = 0.05, 
    #     )
    # end

    # violin!(
    #     ax, 
    #     1:L_autocor, 
    #     bootstrap_samples_autocor[begin+1:end, k, l, :],
    #     color   = (:slategray, 0.4), 
    #     bandwidth = 0.05, 
    #     datalimits = (-1.0, 1.0),
    #     show_median = true,
    #     mediancolor = (:slategray, 0.8),
    # )

    for j in 1:L_autocor
        violin!(
            ax, 
            Iterators.repeated(j, B),
            bootstrap_samples_autocor[begin+j, k, l, :],
            color   = (:slategray, 0.4), 
            bandwidth = 0.0125, 
            datalimits = (-1.0, 1.0),
            show_median = true,
            mediancolor = (:slategray, 0.8),
        )
    end

    stl = stem!(0:L_autocor, obs_autocor[:, k], trunkwidth=2, label="ACF")
    ax.xticks = 0:L_autocor
    ax.yticks = -1.0:0.2:1.

    dist_elem = PolyElement(color=(:slategray, 0.4))
    Legend(fig[4,1], [stl, dist_elem], ["Autocorrelation function", "Bootstrap distribution"], orientation=:horizontal, framevisible=false, tellwidth=false)
    fig
end

plot_acf_distribution(1, 10, sbb_bootstrap_samples_autocor)
plot_acf_distribution(1, 10, mbb_bootstrap_samples_autocor)

# # Save plots for all variables
# for k in 1:10, l in 1:40
#     fig = plot_acf_distribution(1, 2, sbb_bootstrap_samples_autocor)
#     filename = savename("acf", (variable=varnames[k], l=l, method="stationary"), "png", sort=false)
#     save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
# end


## MSE analysis per variable

k = 1
l = 1

sqr(x) = x^2

# Output folder
mkpath(plotsdir(PLOTSDIR, "variable"))

# Plot the overall MSE of the ACF per variable
fig_acf_variable = map(1:K) do k 

    # Compute the MSE for each block length l
    mse_acf = map(1:L_block) do l 
        bootstrap_samples = @view sbb_bootstrap_samples_autocor[begin+1:end, k, l, :]
        obs_autocor_var   = obs_autocor[begin+1:end, k]
        mse = mean(sqr, bootstrap_samples .- obs_autocor_var, dims=2) |> vec
        # mse = mean(sqr, bootstrap_samples .- obs_autocor_var)

        # weights = collect(L_autocor:-1:1) / sum(1:12)
        weights = ones(12) / 12
        sum(weights .* mse)
    end

    # Plot MSE vs block length 
    fig = Figure(size=(950,550))
    title = varnames[k] * ", blocklength=$l"
    lt  = Label(fig[1,1], varnames[k], fontsize=18, font=:bold, tellwidth=false)
    st  = Label(fig[2,1], "", fontsize=18, tellwidth=false)
    ax  = Axis(fig[3,1], xlabel="Block length", ylabel="Mean Squared Error")
    lines!(ax, 1:L_block, mse_acf, linewidth=2)
    filename = savename("acf_mse", (variable=varnames[k], method="stationary"), "png", sort=false)
    save(plotsdir(PLOTSDIR, "variable", filename), fig, px_per_unit=2.0)
    fig
end

fig_acf_variable[1]


## MSE aggregated analysis 

# Computes the MSE of the overall ACF plot (for all lags) for variable k of bootstrap_samples_autocor
function mse_acf_variable(bootstrap_samples_autocor, k) 
    # Compute the MSE for each block length l
    mse_acf = map(1:L_block) do l 
        # Get bootstrap samples for variable k, blocklength l, only lags 1-L_autocor
        bootstrap_samples = @view bootstrap_samples_autocor[begin+1:end, k, l, :]
        obs_autocor_var   = obs_autocor[begin+1:end, k]

        # MSE is the mean over the realizations' dimension (dims=2)
        mse = mean(sqr, bootstrap_samples .- obs_autocor_var, dims=2) |> vec
        # Average the MSE values over the ACF lags
        # weights = collect(L_autocor:-1:1) / sum(1:12)     # More weight to the first lags
        weights = ones(L_autocor) / L_autocor        # Equal weights
        sum(weights .* mse)
    end
end

# Normalized MSE aggregated analysis
function norm_mse_acf_variable(bootstrap_samples_autocor, k) 
    # Compute the MSE for each block length l
    mse_acf = map(1:L_block) do l 
        # Get bootstrap samples for variable k, blocklength l, only lags 1-L_autocor
        bootstrap_samples = @view bootstrap_samples_autocor[begin+1:end, k, l, :]
        obs_autocor_var   = obs_autocor[begin+1:end, k]

        # MSE is the mean over the realizations' dimension (dims=2)
        std_autocor = std(bootstrap_samples, dims=2)
        mse = mean(sqr, ((bootstrap_samples .- obs_autocor_var)./std_autocor), dims=2) |> vec
        # Average the MSE values over the ACF lags
        # weights = collect(L_autocor:-1:1) / sum(1:12)     # More weight to the first lags
        weights = ones(L_autocor) / L_autocor        # Equal weights
        sum(weights .* mse)
    end
end

mse_acf_variable(sbb_bootstrap_samples_autocor, 1)
norm_mse_acf_variable(sbb_bootstrap_samples_autocor, 1)


# Plot the overall MSE of the ACF per variable
sbb_mse_per_variable = mapreduce(hcat, 1:K) do k 
    mse_acf_variable(sbb_bootstrap_samples_autocor, k)
end
mbb_mse_per_variable = mapreduce(hcat, 1:K) do k 
    mse_acf_variable(mbb_bootstrap_samples_autocor, k)
end

# Plot the overall normalized MSE of the ACF per variable
sbb_norm_mse_per_variable = mapreduce(hcat, 1:K) do k 
    norm_mse_acf_variable(sbb_bootstrap_samples_autocor, k)
end
mbb_norm_mse_per_variable = mapreduce(hcat, 1:K) do k 
    norm_mse_acf_variable(mbb_bootstrap_samples_autocor, k)
end

# Overall MSE is the sum of the per_variable MSEs 
sbb_agg_mse_acf = mean(sbb_mse_per_variable, dims=2) |> vec
mbb_agg_mse_acf = mean(mbb_mse_per_variable, dims=2) |> vec

# Overall MSE is the sum of the per_variable MSEs 
norm_sbb_agg_mse_acf = mean(sbb_norm_mse_per_variable, dims=2) |> vec
norm_mbb_agg_mse_acf = mean(mbb_norm_mse_per_variable, dims=2) |> vec



# Plot aggregated MSE vs block length 
fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Stationary & Moving Block Bootstrap aggregated MSE of autocorrelation functions", fontsize=18, font=:bold, tellwidth=false)
st  = Label(fig[2,1], "(ACF lags 1 - $L_autocor)", fontsize=18, tellwidth=false)
ax  = Axis(fig[3,1], xlabel="Block length", ylabel="Mean Squared Error")
lines!(ax, 1:L_block, sbb_agg_mse_acf, linewidth=2, label="Stationary")
lines!(ax, 1:L_block, mbb_agg_mse_acf, linewidth=2, label="Moving")
filename = savename("agg_acf_mse", (method="stationary_moving",), "png", sort=false)
axislegend(framevisible=false)
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig

# Plot aggregated normalized MSE vs block length 
fig = Figure(size=(950,550))
lt  = Label(fig[1,1], "Stationary & Moving Block Bootstrap aggregated normalized MSE of autocorrelation functions", fontsize=18, font=:bold, tellwidth=false)
st  = Label(fig[2,1], "(ACF lags 1 - $L_autocor)", fontsize=18, tellwidth=false)
ax  = Axis(fig[3,1], xlabel="Block length", ylabel="Mean Squared Error")
lines!(ax, 1:L_block, norm_sbb_agg_mse_acf, linewidth=2, label="Stationary")
lines!(ax, 1:L_block, norm_mbb_agg_mse_acf, linewidth=2, label="Moving")
filename = savename("agg_acf_norm_mse", (method="stationary_moving",), "png", sort=false)
axislegend(framevisible=false)
save(plotsdir(PLOTSDIR, filename), fig, px_per_unit=2.0)
fig


norm_sbb_agg_mse_acf

1 - norm_sbb_agg_mse_acf[6]/maximum(norm_sbb_agg_mse_acf)



