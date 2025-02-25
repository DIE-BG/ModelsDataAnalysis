using DrWatson
@quickactivate "ModelsDataAnalysis"

using Statistics
using DependentBootstrap
using CairoMakie
using DataFramesMeta

#include(srcdir("mean_helpers.jl"))

# Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")

mat_d4l_data = Matrix(d4l_data[:, 2:end] |> disallowmissing)

K = size(mat_d4l_data, 2)
B = 10_000

sample_mean = Matrix{Float64}(undef, B, K)

## Plot the mean estimate for different block sizes

L = 91

# Mean of moving block bootstrap
l2means_moving = mapreduce(hcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean(sample_mean, dims=1) |> vec
end 

# Mean of stationary block bootstrap
l2means_stationary = mapreduce(hcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean(sample_mean, dims=1) |> vec
end 

# standar deviation for MBB
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

# standar deviation for SBB
l2std_stationary = mapreduce(hcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    std(sample_mean, dims=1) |> vec
end 

## MSE
actual_means = mean(mat_d4l_data, dims=1)
# Normalized MSE for MBB
l2mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean((sample_mean .- actual_means).^2, dims=1)
end

# MSE for SBB
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


# Normalized MSE for MBB
l2norm_mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end
    
    sample_std = std(sample_mean, dims=1)

    # Compute summary statistic 
    mean(((sample_mean .- actual_means)./sample_std).^2, dims=1)
end

# Normalized MSE for SBB
l2norm_mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    sample_std = std(sample_mean, dims=1)

    # Compute summary statistic 
    mean(((sample_mean .- actual_means)./sample_std).^2, dims=1)
end

# Average MSE of al macroeconomic series
l2mse_moving_all = mean(l2mse_moving, dims = 2)
l2mse_stationary_all = mean(l2mse_stationary, dims = 2)

# Average normalized MSE of al macroeconomic series
l2norm_mse_moving_all = mean(l2norm_mse_moving, dims = 2)
l2norm_mse_stationary_all = mean(l2norm_mse_stationary, dims = 2)

## Plots
varnames = ["total inflation", "core inflation", "import prices", "exchange rate",
            "Monetary base", "external inflation", "policy rate", "external policy rate", "Domestic product",
            "external product"]

var_cod = propertynames(d4l_data)[2:end]

nvar = length(varnames)

mkdir(plotsdir())
mkdir(plotsdir()*"\\mean")

map(1:nvar) do nvar

            fig = Figure(size = (2750,1500), fontsize = 25)

Label(fig[1,1:2],  string(varnames[nvar]), fontsize = 60, tellwidth = false, halign = :center)

# Mean
ax = Axis(fig[2,1], title = "Average (L2) of the historical mean estimator")

lines!(ax, 1:L, l2means_moving[nvar, :], linewidth=2, label = "Moving")
lines!(ax, 1:L, l2means_stationary[nvar, :], linewidth=2, label = "Stationary")
hlines!(ax, actual_means[nvar], color=:red, linewidth=2, linestyle=:dash, label = "Media historica")

axislegend(position = :rb, framevisible = false)

# Standard deviation
ax = Axis(fig[2,2], title = "Average (L2) \nof the historical mean estimator")

lines!(ax, 1:L, l2std_moving[nvar, :], label="Moving")
lines!(ax, 1:L, l2std_stationary[nvar, :], label="Stationary")

axislegend(position = :rt, framevisible = false)

# MSE moving
ax = Axis(fig[3,1], title = "Mean square error \nof the historial mean estimator",
            subtitle = "Moving")

lines!(ax, 1:L, l2mse_moving[:, nvar])

# MSE stationary
ax = Axis(fig[3,2], title = "Mean square error \nof the historical mean estimator",
            subtitle = "Stationary")

lines!(ax, 1:L, l2mse_stationary[:,nvar])

save(plotsdir()*"\\mean\\"*string(nvar)*"_"*string(var_cod[nvar])*".png", fig, px_per_unit=2.0)

end

# Average MSE of all variables
fig = Figure(size = (1200, 600))
ax = Axis(fig[1,1], title = "Average MSE of the historial mean estimator", 
         subtitle = "Moving",
         xlabel = L"\text{Block length } l")

lines!(ax, 1:L, l2mse_moving_all[:,1])

ax = Axis(fig[1,2], title = "Average MSE of the historial mean estimator",
          subtitle = "Stationary",
          xlabel = L"\text{Block length } l")
lines!(ax, 1:L, l2mse_stationary_all[:,1])
fig

# Average Normalized MSE of all variables
fig = Figure(size = (1200, 600))
ax = Axis(fig[1,1], title = "Average MSE of the historial mean estimator", 
         subtitle = "Moving",
         xlabel = L"\text{Block length } l")

lines!(ax, 1:L, l2mse_moving_all[:,1])

ax = Axis(fig[1,2], title = "Average MSE of the historial mean estimator",
          subtitle = "Stationary",
          xlabel = L"\text{Block length } l")
lines!(ax, 1:L, l2mse_stationary_all[:,1])
fig

save(plotsdir()*"\\mean\\"*"MSE_stationary.png", fig, px_per_unit=2.0)

