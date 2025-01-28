using DrWatson
@quickactivate "ModelsDataAnalysis"

using Statistics
using DependentBootstrap
using CairoMakie
using DataFramesMeta

include(srcdir("mean_helpers.jl"))

# Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")

mat_d4l_data = Matrix(d4l_data[:, 2:end] |> disallowmissing)
# dboot(mat_d4l_data, blocklength=6, bootmethod=:moving, flevel1=mean)
# optblocklength(mat_d4l_data, blmethod=:ppw2009, fblocklengthcombine=mean)
# optblocklength(mat_d4l_data, BootInput(bootmethod=:moving, blmethod=:ppw2009))

K = size(mat_d4l_data, 2)
B = 1_000

# Get the distribution of the sample mean for each variable
inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=5, numresample=B)

sample_mean = Matrix{Float64}(undef, B, K)

for b in 1:B 
    sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
end

d4l_data.D4L_CPI |> mean

fig, ax = CairoMakie.hist(sample_mean[:, 1], normalization=:probability)
vlines!(ax, d4l_data.D4L_CPI |> mean)
fig

## Plot the mean estimate for different block sizes

L = 91

l2means = mapreduce(hcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean(sample_mean, dims=1) |> vec
end 

## Plot means vs block size
nvar = 1
actual_means = mean(mat_d4l_data, dims=1)
varnames = propertynames(d4l_data)[2:end]

fig = Figure()
ax  = Axis(fig[1,1], title=string(varnames[nvar]))


lines!(ax, 1:L, l2means[nvar, :])
hlines!(ax, actual_means[nvar], color=:red, linewidth=2, linestyle=:dash)
fig


## Plot the variance of the sample mean for different block size

l2std_moving = mapreduce(hcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    std(sample_mean, dims=1) |> vec
end 

## Plot means vs block size
nvar = 1
actual_std = std(mat_d4l_data, dims=1)
varnames = propertynames(d4l_data)[2:end]

fig = Figure()
ax  = Axis(fig[1,1], title=string(varnames[nvar]))

lines!(ax, 1:L, l2std_stationary[nvar, :], label="Stationary")
lines!(ax, 1:L, l2std_moving[nvar, :], label="Moving")
axislegend()
# hlines!(ax, actual_std[nvar], color=:red, linewidth=2, linestyle=:dash)
fig



## Plot the MSE of the sample mean for different block size

actual_means = mean(mat_d4l_data, dims=1)

l2mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean((sample_mean .- actual_means).^2, dims=1)
end 

## Plot means vs block size
nvar = 3
actual_std = std(mat_d4l_data, dims=1)
varnames = propertynames(d4l_data)[2:end]

fig = Figure()
ax  = Axis(fig[1,1], title=string(varnames[nvar]))

lines!(ax, 1:L, l2mse_moving[:, nvar], label="Moving")
lines!(ax, 1:L, l2mse_stationary[:, nvar], label="Stationary")

axislegend()
fig

## 

savepath = mkpath(plotsdir("mean"))
save(plotsdir(savepath, "a.png"), fig, px_per_unit=2.0)


## To-do 
# For each variable compare between :moving, :stationary, :nooverlap
# - Plots of mean, standard deviation and MSE of the vs. block size 