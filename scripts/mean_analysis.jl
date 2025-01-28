using DrWatson
@quickactivate "ModelsDataAnalysis"

using Statistics
using DependentBootstrap
using CairoMakie

# Load data 
d4l_data = load(datadir("data_QPM.jld2"), "d4l_data")

mat_d4l_data = Matrix(d4l_data[:, 2:end] |> disallowmissing)
dboot(mat_d4l_data, blocklength=6, bootmethod=:moving, flevel1=mean)


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




CairoMakie.hist(sample_mean[:, 2], normalization=:probability)
