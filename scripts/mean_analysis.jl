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

## Plot the mean estimate for different block sizes

L = 91

# Mean of moving blocks bootstrap
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

# Mean of stationary blocks bootstrap
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

## Plot the MSE of the sample mean for different block size

actual_means = mean(mat_d4l_data, dims=1)
# MSE for MBB
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

## Plots
#varnames = propertynames(d4l_data)[2:end]
varnames = ["Inflación total", "Inflación subyacente", "Precios de Importados", "Tipo de Cambio",
            "Base Monetaria", "Inflación externa", "Tasa líder", "Tasa externa", "Producto doméstico",
            "Producto Externo"]

var = propertynames(d4l_data)[2:end]

nvar = length(varnames)

map(1:nvar) do nvar

            fig = Figure(size = (2500,1000), fontsize = 25)

Label(fig[1, 2],  string(varnames[nvar]), fontsize = 60, tellwidth = false, halign = :center)

# Mean
ax = Axis(fig[2,1], title = "Media de cada realización")

lines!(ax, 1:L, l2means_moving[nvar, :], linewidth=2, label = "Moving")
lines!(ax, 1:L, l2means_stationary[nvar, :], linewidth=2, label = "Stationary")
hlines!(ax, actual_means[nvar], color=:red, linewidth=2, linestyle=:dash, label = "Media historica")

axislegend(position = :rb, framevisible = false)

# Standar deviation
ax = Axis(fig[2,2], title = "Desviación estandar de las medias \nen cada realización")

lines!(ax, 1:L, l2std_moving[nvar, :], label="Moving")
lines!(ax, 1:L, l2std_stationary[nvar, :], label="Stationary")

axislegend(position = :rt, framevisible = false)

# MSE
ax = Axis(fig[2,3], title = "Error Cuadrático Medio \nen cada realización")

lines!(ax, 1:L, l2mse_moving[:, nvar], label="Moving")
lines!(ax, 1:L, l2mse_stationary[:, nvar], label="Stationary")

axislegend(position = :rt, framevisible = false)


fig

save(plotsdir(string(var[nvar])*".png"), fig, px_per_unit=2.0)

end
