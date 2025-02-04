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

# MSE for MBB
actual_means = mean(mat_d4l_data, dims=1)
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

## Bias and variance for MBB and SBB
# Bias MBB
bias_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    (mean(sample_mean, dims = 1) .- actual_means).^2
end 

variance_moving = l2mse_moving .- bias_moving

# Bias SBB
bias_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_mean[b, :] = mean(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    (mean(sample_mean, dims = 1) .- actual_means).^2
end 

# Variance MBB
variance_stationary = l2mse_stationary .- bias_stationary

## Plots
varnames = ["Inflación total", "Inflación subyacente", "Precios de Importados", "Tipo de Cambio",
            "Base Monetaria", "Inflación externa", "Tasa líder", "Tasa externa", "Producto doméstico",
            "Producto Externo"]

Bvar_cod = propertynames(d4l_data)[2:end]

nvar = length(varnames)

map(1:nvar) do nvar

            fig = Figure(size = (2750,1000), fontsize = 25)

Label(fig[1, 2],  string(varnames[nvar]), fontsize = 60, tellwidth = false, halign = :center)

# Mean
ax = Axis(fig[2,1], title = "Promedio (L2) del estimador de la media histórica")

lines!(ax, 1:L, l2means_moving[nvar, :], linewidth=2, label = "Moving")
lines!(ax, 1:L, l2means_stationary[nvar, :], linewidth=2, label = "Stationary")
hlines!(ax, actual_means[nvar], color=:red, linewidth=2, linestyle=:dash, label = "Media historica")

axislegend(position = :rb, framevisible = false)

# Standar deviation
ax = Axis(fig[2,2], title = "Desviación estándar (L2) \ndel estimador de la media histórica")

lines!(ax, 1:L, l2std_moving[nvar, :], label="Moving")
lines!(ax, 1:L, l2std_stationary[nvar, :], label="Stationary")

axislegend(position = :rt, framevisible = false)

# MSE
ax = Axis(fig[2,3], title = "Error Cuadrático Medio \ndel estimador de la media histórica")

lines!(ax, 1:L, l2mse_moving[:, nvar], label="Moving")
lines!(ax, 1:L, l2mse_stationary[:, nvar], label="Stationary")
band!(1:L, repeat([0], L), bias_stationary[:,nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35, label = "sesgo")
band!(1:L, bias_stationary[:,nvar], variance_stationary[:,nvar].+ bias_stationary[:,nvar], color = :red, alpha = 0.35, label = "varianza")

axislegend(position = :rt, framevisible = false)

save(plotsdir("mean//"*string(nvar)*"_"*string(var_cod[nvar])*".png"), fig, px_per_unit=2.0)

end

## Plots bias and variance
map(1:nvar) do nvar

   fig = Figure(size = (950, 600))
   
   ax = Axis(fig[1,1], title = "Sesgo y varianza \nStationary Block Bootstrap (SBB)",
             subtitle = string(var_cod[nvar]))
   lines!(1:L, l2mse_stationary[:,nvar], linewidth = 2, label = "MSE SBB")
   band!(1:L, repeat([0], L), bias_stationary[:,nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35, 
        label = "Sesgo")
   band!(1:L, bias_stationary[:,nvar], variance_stationary[:,nvar].+ bias_stationary[:,nvar], color = :red, alpha = 0.35,
        label = "Varianza")

    axislegend(position = :rt, framevisible = false)
    
    fig

    save(plotsdir("mean\\bias and variance"*"\\"*string(nvar)*"_"*string(var_cod[nvar])*"_stationary.png"), fig, px_per_unit = 2.0)
end
