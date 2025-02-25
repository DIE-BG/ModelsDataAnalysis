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
sample_std = Matrix{Float64}(undef, B, K)

## Plot the variance estimate for different block sizes
L = 91

## Mean
# mean of variance for moving block bootstrap
actual_variance = var(mat_d4l_data, dims = 1)
l2vars_moving = mapreduce(vcat, 1:L) do l

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
l2vars_stationary = mapreduce(vcat, 1:L) do l

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

# Normalized MSE for MBB
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

# Normalized MSE for SBB
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

# Normalized MSE for MBB
l2norm_mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:moving, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    sample_std = std(sample_variance, dims=1)

    # Compute summary statistic 
    mean(((sample_variance .- actual_variance)./sample_std).^2, dims=1)
end

# Normalized MSE for SBB
l2norm_mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_d4l_data, bootmethod=:stationary, blocklength=l, numresample=B)
    
    # Perform the block bootstrap with length l
    for b in 1:B 
        sample_variance[b, :] = var(mat_d4l_data[inds[b], :], dims=1)
    end

    sample_std = std(sample_variance, dims=1)

    # Compute summary statistic 
    mean(((sample_variance .- actual_variance)./sample_std).^2, dims=1)
end

# Average MSE of al macroeconomic series
l2mse_moving_all = mean(l2mse_moving, dims = 2)
l2mse_stationary_all = mean(l2mse_stationary, dims = 2)

# Average MSE of al macroeconomic series
l2norm_mse_moving_all = mean(l2norm_mse_moving, dims = 2)
l2norm_mse_stationary_all = mean(l2norm_mse_stationary, dims = 2)

## Plots
varnames = ["total inflation", "core inflation", "import prices", "exchange rate",
            "Monetary base", "external inflation", "policy rate", "external policy rate", "Domestic product",
            "external product"]

Bvar_cod = propertynames(d4l_data)[2:end]

nvar = length(varnames)

mkdir(plotsdir())
mkdir(plotsdir()*"\\variance")

map(1:nvar) do nvar

    fig = Figure(size = (2750,1000), fontsize = 25)

    Label(fig[1,1:2],  string(varnames[nvar]), fontsize = 60, tellwidth = false, halign = :center)

    # variance
    ax = Axis(fig[2,1], title = "Average (L2) of the historical variance estimator")

    lines!(ax, 1:L, l2vars_moving[:,nvar], linewidth=2, label = "Moving")
    lines!(ax, 1:L, l2vars_stationary[:,nvar], linewidth=2, label = "Stationary")
    hlines!(ax, actual_variance[nvar], color=:red, linewidth=2, linestyle=:dash, label = "historical variance")

    axislegend(position = :rb, framevisible = false)

    # Standar deviation
    ax = Axis(fig[2,2], title = "Standard deviation (L2) \nof the historical variance estimator")

    lines!(ax, 1:L, l2std_moving[:,nvar], label="Moving")
    lines!(ax, 1:L, l2std_stationary[:,nvar], label="Stationary")

    axislegend(position = :rt, framevisible = false)

    # MSE
    ax = Axis(fig[2,3], title = "square mean error \nof the historical variance estimator")

    lines!(ax, 1:L, l2mse_moving[:,nvar], label="Moving")
    lines!(ax, 1:L, l2mse_stationary[:,nvar], label="Stationary")

    axislegend(position = :rt, framevisible = false)

    save(plotsdir("variance\\"*string(nvar)*"_"*string(var_cod[nvar])*".png"), fig, px_per_unit=2.0)

end

# Average MSE of all variables

fig = Figure(size = (950, 600))
ax = Axis(fig[1,1], title = "Average MSE of the historial variance estimator", subtitle = "Moving")
lines!(ax, 1:L, l2mse_moving_all[:,1])
save(plotsdir()*"\\variance\\"*"MSE_moving.png", fig, px_per_unit=2.0)

fig = Figure(size = (950, 600))
ax = Axis(fig[1,1], title = "Average MSE of the historial variance estimator", subtitle = "Moving")
lines!(ax, 1:L, l2mse_stationary_all[:,1])
save(plotsdir()*"\\variance\\"*"MSE_stationary.png", fig, px_per_unit=2.0)



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

 # Average MSE of all variables
fig = Figure(size = (1200, 600))
ax = Axis(fig[1,1], title = "Average MSE of the historial variance estimator", 
         subtitle = "Moving",
         xlabel = L"\text{Block length } l")

lines!(ax, 1:L, l2mse_moving_all[:,1])
ylims!(ax, 0, 16)

ax = Axis(fig[1,2], title = "Average MSE of the historial variance estimator",
          subtitle = "Stationary",
          xlabel = L"\text{Block length } l")
lines!(ax, 1:L, l2mse_stationary_all[:,1])
ylims!(ax, 0, 16)
save(plotsdir()*"\\variance\\"*"all_MSE.png", fig, px_per_unit=2.0)

# Average Normalized MSE of all variables
fig = Figure(size = (1200, 600))
ax = Axis(fig[1,1], title = "Average normalized MSE of the historial mean estimator", 
         subtitle = "Moving",
         xlabel = L"\text{Block length } l")

lines!(ax, 1:L, l2norm_mse_moving_all[:,1])

ax = Axis(fig[1,2], title = "Average normalized MSE of the historial mean estimator",
          subtitle = "Stationary",
          xlabel = L"\text{Block length } l")
lines!(ax, 1:L, l2norm_mse_stationary_all[:,1])

save(plotsdir()*"\\variance\\"*"all_normalized_MSE.png", fig, px_per_unit=2.0)