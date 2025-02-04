using DrWatson
@quickactivate "ModelsDataAnalysis"

using Statistics
using DependentBootstrap
using CairoMakie
using DataFramesMeta

# Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")

mat_d4l_data = Matrix(d4l_data[:, 2:end] |> disallowmissing)

K = size(mat_d4l_data, 2)
B = 10_000

sample_variance = Matrix{Float64}(undef, B, K)

## Plot the variance estimate for different block sizes
L = 91

## Mean
# mean of variance for moving block bootstrap
actual_variance = var(mat_d4l_data, dims = 1)
l2means_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean(sample_variance, dims=1)

end 

# mean of variance for stationary block bootstrap
l2means_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean(sample_variance, dims=1) 

end 

## Standar Deviation
# standar deviation of variance for moving block bootstrap
l2std_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    std(sample_variance, dims=1)

end 

# standar deviation of variance for stationary block bootstrap
l2std_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    std(sample_variance, dims=1) 

end 

## MSE
# MSE for MBB
l2mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean((sample_variance .- actual_variance).^2, dims=1)
end 

# MSE for SBB
l2mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    mean((sample_variance .- actual_variance).^2, dims=1)
end 

## Bias and variance for MBB and SBB
# Bias MBB
bias_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    (mean(sample_variance, dims = 1) .- actual_variance).^2
end 

variance_moving = l2mse_moving .- bias_moving

# Bias SBB
bias_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    # Compute summary statistic 
    (mean(sample_variance, dims = 1) .- actual_variance).^2
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
    ax = Axis(fig[2,1], title = "Promedio (L2) del estimador de la varianza histórica")

    lines!(ax, 1:L, l2means_moving[:,nvar], linewidth=2, label = "Moving")
    lines!(ax, 1:L, l2means_stationary[:,nvar], linewidth=2, label = "Stationary")
    hlines!(ax, actual_variance[nvar], color=:red, linewidth=2, linestyle=:dash, label = "varianza historica")

    axislegend(position = :rb, framevisible = false)

    # Standar deviation
    ax = Axis(fig[2,2], title = "Desviación estándar (L2) \ndel estimador de la varianza histórica")

    lines!(ax, 1:L, l2std_moving[:,nvar], label="Moving")
    lines!(ax, 1:L, l2std_stationary[:,nvar], label="Stationary")

    axislegend(position = :rt, framevisible = false)

    # MSE
    ax = Axis(fig[2,3], title = "Error Cuadrático Medio \ndel estimador de la varianza histórica")

    lines!(ax, 1:L, l2mse_moving[:,nvar], label="Moving")
    lines!(ax, 1:L, l2mse_stationary[:,nvar], label="Stationary")
    band!(1:L, repeat([0], L), bias_stationary[:,nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35, label = "sesgo")
    band!(1:L, bias_stationary[:,nvar], variance_stationary[:,nvar].+ bias_stationary[:,nvar], color = :red, alpha = 0.35, label = "varianza")

    axislegend(position = :rt, framevisible = false)

    save(plotsdir("variance\\"*string(nvar)*"_"*string(var_cod[nvar])*".png"), fig, px_per_unit=2.0)

end

## Plots bias and variance
map(1:nvar) do nvar

    fig = Figure(size = (1400, 600))
    
    ax = Axis(fig[1,1], title = "Sesgo y varianza \nStationary Block Bootstrap (SBB)",
              subtitle = string(var_cod[nvar]))
    lines!(1:L, l2mse_stationary[:,nvar], linewidth = 2, label = "MSE SBB")
    band!(1:L, repeat([0], L), bias_stationary[:,nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35, 
         label = "Sesgo")
    band!(1:L, bias_stationary[:,nvar], variance_stationary[:,nvar].+ bias_stationary[:,nvar], color = :red, alpha = 0.35,
         label = "Varianza")
 
     axislegend(position = :rt, framevisible = false)
     
    ax = Axis(fig[1,2], title = "Sesgo y varianza \nMoving Block Bootstrap (MBB)",
              subtitle = string(var_cod[nvar]))
    lines!(1:L, l2mse_moving[:,nvar], linewidth = 2, label = "MSE MBB")
    band!(1:L, repeat([0], L), bias_moving[:,nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35, 
         label = "Sesgo")
    band!(1:L, bias_moving[:,nvar], variance_moving[:,nvar].+ bias_moving[:,nvar], color = :red, alpha = 0.35,
         label = "Varianza")
 
     axislegend(position = :rt, framevisible = false)


     fig

     save(plotsdir("variance\\bias and variance"*"\\"*string(nvar)*"_"*string(var_cod[nvar])*"_bias_variance.png"), fig, px_per_unit = 2.0)
 end