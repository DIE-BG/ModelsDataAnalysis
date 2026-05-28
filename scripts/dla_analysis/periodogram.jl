using DrWatson
@quickactivate "ModelsDataAnalysis"

using Dates
using CairoMakie
using DSP
using DataFrames
import Printf
using Statistics: mean
using DataFramesMeta: disallowmissing!

# Helper functions
ft1(x) = Printf.@sprintf("%.1f", x)
ft2(x) = Printf.@sprintf("%.2f", x)


## Load data
WITH_REMITTANCES = false  # Configuration: choose whether to include remittances in the analysis
DLA_DATA = load_dla_data(WITH_REMITTANCES)
# Set output directory based on dataset choice
PLOTSDIR = mkpath(plotsdir("dla_analysis", "periodograms", WITH_REMITTANCES ? "with_rem_gdp" : "without_rem_gdp"))


## Periodogram function

function plot_periodogram(y, varname, dates; centered = true)
    ## Centered variable
    centered_y = y .- mean(y)
    if centered
        pg = periodogram(centered_y)
    else
        pg = periodogram(y)
    end
    norm_freq = freq(pg)
    periodicity_str = @. ft2(norm_freq) * " (" * ft1(1 / norm_freq) * ")"
    periodicity_str[1] = "0 (∞)"


    # Periodogram Makie figure
    fig = Figure(size = (950, 650))
    lt = Label(fig[1, 1:2], varname, fontsize = 18, font = :bold)

    # Time series plot
    ratadates = datetime2rata.(dates)
    ax1 = Axis(fig[2, 1:2], title = "Time series")
    lines!(ax1, ratadates, y)
    dates_str = Dates.format.(dates, dateformat"yyyy:mm")
    ax1.xticks = (ratadates[begin:6:end], dates_str[begin:6:end])
    ax1.xticklabelrotation = pi / 4

    # Periodogram
    ax2_title = centered ? "Periodogram (centered variable)" : "Periodogram"
    ax2 = Axis(fig[3, 1:2], title = ax2_title, xlabel = "Normalized frequency (Period in quarters)")
    stem!(ax2, norm_freq, power(pg))
    elems = [1:5..., 6:2:length(norm_freq)...]
    ax2.xticks = (norm_freq[elems], periodicity_str[elems])
    ax2.xticklabelrotation = pi / 4
    ax2.xticklabelsize = 12

    return fig
end

## Example periodogram

t = 1:size(DLA_DATA, 1)
y = cos.(2 * pi * t / 8)

# Compute the periodogram
pg = periodogram(y)
var(y)  # = 0.5

# This function approximates the area under the periodogram pg = periodogram(y)
# The result should be close to the variance of y (var(y))
function periodogram_area(pg)
    pow = power(pg)
    f = freq(pg)
    area = 0.0
    L = length(pow)
    for l in 2:L
        area += pow[l] * (f[l] - f[l - 1])
    end
    return area
end

periodogram_area(pg)    # 0.4945

# Plot the example periodogram
example_fig = plot_periodogram(y, L"Cosine function $y = \cos\left(\frac{2\pi t}{8} \right)$", DLA_DATA.Date)
save(plotsdir(PLOTSDIR, "example_periodogram.png"), example_fig, px_per_unit = 2.0)


## Compute the periodogram for all variables in DLA_DATA

varnames = string.(propertynames(DLA_DATA))[(begin + 1):end]
dates = DLA_DATA.Date
centered = true

map(varnames) do varname
    # Get the periodogram figure
    y = DLA_DATA[!, varname]
    mask = @. !ismissing(y)
    fig = plot_periodogram(disallowmissing(y[mask]), varname, dates[mask]; centered)

    # Save plot
    plot_name = savename("periodogram", (; variable = varname, centered), "png")
    save(plotsdir(PLOTSDIR, plot_name), fig)
end


## Compute the periodogram for all variables in DLA_DATA

varnames = string.(propertynames(DLA_DATA))[(begin + 1):end]
dates = DLA_DATA.Date
centered = false

map(varnames) do varname
    # Get the periodogram figure
    y = DLA_DATA[!, varname]
    mask = @. !ismissing(y)
    fig = plot_periodogram(disallowmissing(y[mask]), varname, dates[mask]; centered)

    # Save plot
    plot_name = savename("periodogram", (; variable = varname, centered), "png")
    save(plotsdir(PLOTSDIR, plot_name), fig)
end
