using DrWatson
@quickactivate "ModelsDataAnalysis"

using Statistics
using DependentBootstrap
using CairoMakie
using DataFramesMeta

## Load data
WITH_REMITTANCES = false  # Configuration: choose whether to include remittances in the analysis
DLA_DATA = load_dla_data(WITH_REMITTANCES)
# Set output directory based on dataset choice
PLOTSDIR = mkpath(plotsdir("dla_analysis", "variance", WITH_REMITTANCES ? "with_rem_gdp" : "without_rem_gdp"))


# Matrix data
mat_dla_data = Matrix(DLA_DATA[:, 2:end]) |> disallowmissing


K = size(mat_dla_data, 2)
B = 10_000

sample_variance = Matrix{Float64}(undef, B, K)
sample_std = Matrix{Float64}(undef, B, K)

## Plot the variance estimate for different block sizes
L = 91

## Mean
# mean of variance for moving block bootstrap
actual_variance = var(mat_dla_data, dims = 1)
l2vars_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :moving, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    mean(sample_variance, dims = 1)

end

# mean of variance for stationary block bootstrap
l2vars_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :stationary, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    mean(sample_variance, dims = 1)

end

## Standar Deviation
# standar deviation of variance for moving block bootstrap
l2std_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :moving, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    std(sample_variance, dims = 1)

end

# standar deviation of variance for stationary block bootstrap
l2std_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :stationary, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    std(sample_variance, dims = 1)

end

# Normalized MSE for MBB
l2mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :moving, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    mean((sample_variance .- actual_variance) .^ 2, dims = 1)
end

# Normalized MSE for SBB
l2mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :stationary, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    # Compute summary statistic
    mean((sample_variance .- actual_variance) .^ 2, dims = 1)
end

# Normalized MSE for MBB
l2norm_mse_moving = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :moving, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    sample_std = std(sample_variance, dims = 1)

    # Compute summary statistic
    mean(((sample_variance .- actual_variance) ./ sample_std) .^ 2, dims = 1)
end

# Normalized MSE for SBB
l2norm_mse_stationary = mapreduce(vcat, 1:L) do l

    # Compute indices
    inds = dbootinds(mat_dla_data, bootmethod = :stationary, blocklength = l, numresample = B)

    # Perform the block bootstrap with length l
    for b in 1:B
        sample_variance[b, :] = var(mat_dla_data[inds[b], :], dims = 1)
    end

    sample_std = std(sample_variance, dims = 1)

    # Compute summary statistic
    mean(((sample_variance .- actual_variance) ./ sample_std) .^ 2, dims = 1)
end

# Average MSE of al macroeconomic series
l2mse_moving_all = mean(l2mse_moving, dims = 2)
l2mse_stationary_all = mean(l2mse_stationary, dims = 2)

# Average MSE of al macroeconomic series
l2norm_mse_moving_all = mean(l2norm_mse_moving, dims = 2)
l2norm_mse_stationary_all = mean(l2norm_mse_stationary, dims = 2)

## Plots
varnames = [
    "total inflation", "core inflation", "import prices", "exchange rate",
    "Monetary base", "external inflation", "Domestic product", "external product",
    "policy rate", "external policy rate",
]

Bvar_cod = propertynames(DLA_DATA)[(begin + 1):end]

nvar = length(varnames)


map(1:nvar) do nvar

    fig = Figure(size = (2750, 1000), fontsize = 25)

    Label(fig[1, 1:2], string(varnames[nvar]), fontsize = 60, tellwidth = false, halign = :center)

    # variance
    ax = Axis(fig[2, 1], title = "Average (L2) of the historical variance estimator")

    lines!(ax, 1:L, l2vars_moving[:, nvar], linewidth = 2, label = "Moving")
    lines!(ax, 1:L, l2vars_stationary[:, nvar], linewidth = 2, label = "Stationary")
    hlines!(ax, actual_variance[nvar], color = :red, linewidth = 2, linestyle = :dash, label = "historical variance")

    axislegend(position = :rb, framevisible = false)

    # Standar deviation
    ax = Axis(fig[2, 2], title = "Standard deviation (L2) \nof the historical variance estimator")

    lines!(ax, 1:L, l2std_moving[:, nvar], label = "Moving")
    lines!(ax, 1:L, l2std_stationary[:, nvar], label = "Stationary")

    axislegend(position = :rt, framevisible = false)

    # MSE
    ax = Axis(fig[2, 3], title = "square mean error \nof the historical variance estimator")

    lines!(ax, 1:L, l2mse_moving[:, nvar], label = "Moving")
    lines!(ax, 1:L, l2mse_stationary[:, nvar], label = "Stationary")

    axislegend(position = :rt, framevisible = false)

    save(plotsdir(PLOTSDIR, "$(nvar)_$(Bvar_cod[nvar]).png"), fig, px_per_unit = 2.0)

end

## Average MSE of all variables

fig = Figure(size = (950, 600))
ax = Axis(fig[1, 1], title = "Average MSE of the historial variance estimator", subtitle = "Moving")
lines!(ax, 1:L, l2mse_moving_all[:, 1])
save(plotsdir(PLOTSDIR, "MSE_moving.png"), fig, px_per_unit = 2.0)
fig

##
fig = Figure(size = (950, 600))
ax = Axis(fig[1, 1], title = "Average MSE of the historial variance estimator", subtitle = "Moving")
lines!(ax, 1:L, l2mse_stationary_all[:, 1])
save(plotsdir(PLOTSDIR, "MSE_stationary.png"), fig, px_per_unit = 2.0)
fig

## Plots bias and variance
#=
map(1:nvar) do nvar

    fig = Figure(size = (1400, 600))

    ax = Axis(
        fig[1, 1], title = "Sesgo y varianza \nStationary Block Bootstrap (SBB)",
        subtitle = string(Bvar_cod[nvar])
    )
    lines!(1:L, l2mse_stationary[:, nvar], linewidth = 2, label = "MSE SBB")
    band!(
        1:L, repeat([0], L), bias_stationary[:, nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35,
        label = "Sesgo"
    )
    band!(
        1:L, bias_stationary[:, nvar], variance_stationary[:, nvar] .+ bias_stationary[:, nvar], color = :red, alpha = 0.35,
        label = "Varianza"
    )

    axislegend(position = :rt, framevisible = false)

    ax = Axis(
        fig[1, 2], title = "Sesgo y varianza \nMoving Block Bootstrap (MBB)",
        subtitle = string(Bvar_cod[nvar])
    )
    lines!(1:L, l2mse_moving[:, nvar], linewidth = 2, label = "MSE MBB")
    band!(
        1:L, repeat([0], L), bias_moving[:, nvar], color = RGBf(0.008, 0.467, 0.878), alpha = 0.35,
        label = "Sesgo"
    )
    band!(
        1:L, bias_moving[:, nvar], variance_moving[:, nvar] .+ bias_moving[:, nvar], color = :red, alpha = 0.35,
        label = "Varianza"
    )

    axislegend(position = :rt, framevisible = false)
    save(plotsdir(PLOTSDIR, "bias_variance_$(nvar)_$(Bvar_cod[nvar]).png"), fig, px_per_unit = 2.0)
    fig
end
=#

## Average MSE of all variables
fig = Figure(size = (1500, 600), fontsize = 20)
ax = Axis(
    fig[1, 1], title = "Average MSE of the historial variance estimator",
    subtitle = "Moving",
    xlabel = L"\text{Block length } l"
)

lines!(ax, 1:L, l2mse_moving_all[:, 1])
#ylims!(ax, 0, 16.5)

ax = Axis(
    fig[1, 2], title = "Average MSE of the historial variance estimator",
    subtitle = "Stationary",
    xlabel = L"\text{Block length } l"
)
lines!(ax, 1:L, l2mse_stationary_all[:, 1])
#ylims!(ax, 0, 16.5)
save(plotsdir(PLOTSDIR, "all_MSE.png"), fig, px_per_unit = 2.0)
fig


## Average Normalized MSE of all variables
fig = Figure(size = (1500, 600), fontsize = 20)
ax = Axis(
    fig[1, 1], title = "Average normalized MSE of the historial mean estimator",
    subtitle = "Moving",
    xlabel = L"\text{Block length } l"
)

lines!(ax, 1:L, l2norm_mse_moving_all[:, 1])

ax = Axis(
    fig[1, 2], title = "Average normalized MSE of the historial mean estimator",
    subtitle = "Stationary",
    xlabel = L"\text{Block length } l"
)
lines!(ax, 1:L, l2norm_mse_stationary_all[:, 1])
save(plotsdir(PLOTSDIR, "all_normalized_MSE.png"), fig, px_per_unit = 2.0)
fig
