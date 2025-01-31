using DrWatson
@quickactivate "ModelsDataAnalysis"

using Dates
using Statistics
using CSV
using DataFrames
using DataFramesMeta
using TimeSeriesEcon

d4lfn(x) = vcat(missings(4), x[5:end] .- x[1:end-4])
dlafn(x) = vcat(missings(1), 4 * (x[begin+1:end] .- x[begin:end-1]))

monthly_data = CSV.read(datadir("monthly.csv"), DataFrame, ignoreemptyrows=true, comment=".")

# Get quarterly data from monthly data
quarterly_from_monthly_data = @chain monthly_data begin 
    @rtransform(
        :Date = Date(parse.(Int, split(:DateString, "M"))...), 
        :IPEI = :ind_prec_expus * :A_prom + :ind_prec_impus * (1 - :A_prom),
    )
    # Combine monthly data to form quarterly data
    @rtransform(
        :Year    = Dates.year(:Date), 
        :Quarter = quarterofyear(:Date), 
    )
    groupby([:Year, :Quarter])
    @combine(
        :CPI    = mean(:CPI), 
        :CPIXFE = mean(:CPIXFE), 
        :IPEI   = mean(:IPEI), 
        :S      = mean(:S), 
        :MB     = mean(:MB),
        :CPI_RW = mean(:CPI_RW),  
        :RS     = mean(:RS),
        :RS_RW  = mean(:RS_RW),
        :REM    = sum(:REM)
    )
    @transform(:Date = Date.(:Year, :Quarter*3))
    @select(:Date, Not([:Year, :Quarter]))
end


## Read quarterly data

raw_quarterly_data = CSV.read(datadir("quarterly.csv"), DataFrame, ignoreemptyrows=true, comment=".")

# Pre-process quarterly data
quarterly_data = @chain raw_quarterly_data begin
    @rtransform _ @astable begin 
        vals = parse.(Int, split(:DateString, "Q"))
        :Year = first(vals)
        :Quarter = last(vals)
    end
    @transform(:Date   = Date.(:Year, :Quarter*3))
    @select(:Date, Not(:DateString, :Year, :Quarter))
end

## REM/GDP data
rem_gdp_data = @chain quarterly_data begin 
    leftjoin(_, quarterly_from_monthly_data, on=:Date)
    @rtransform(:REM_GDP = (:REM * :S) / :NGDP)
    @select(:Date, :REM_GDP)
end

## D4L data (approximated year-on-year changes by difflog(., 4))

d4l_quarterly_from_monthly_data = @chain quarterly_from_monthly_data begin 
    @transform(
        :L_CPI = 100log.(:CPI), 
        :L_CPIXFE = 100log.(:CPIXFE), 
        :L_IPEI = 100log.(:IPEI), 
        :L_S = 100log.(:S), 
        :L_MB = 100log.(:MB),
        :L_CPI_RW = 100log.(:CPI_RW),  
    )
    @transform(
        :D4L_CPI = d4lfn(:L_CPI), 
        :D4L_CPIXFE = d4lfn(:L_CPIXFE), 
        :D4L_IPEI = d4lfn(:L_IPEI), 
        :D4L_S = d4lfn(:L_S), 
        :D4L_MB = d4lfn(:L_MB),
        :D4L_CPI_RW = d4lfn(:L_CPI_RW),  
    )
    @select(:Date, Cols(r"D4L_"), :RS, :RS_RW)
    @rsubset(:Date >= Date(2002,1) && :Date <= Date(2024,11))
    disallowmissing
end

d4l_quarterly_data = @chain quarterly_data begin 
    # Transform the data
    @transform(
        :L_GDP = 100log.(:GDP), 
        :L_NGDP = 100log.(:NGDP), 
        :L_GDP_RW = 100log.(:GDP_RW), 
    )
    @transform(
        :D4L_GDP = d4lfn(:L_GDP), 
        :D4L_NGDP = d4lfn(:L_NGDP), 
        :D4L_GDP_RW = d4lfn(:L_GDP_RW), 
    )
    @select(:Date, Cols(r"D4L_"))
    @rsubset(:Date >= Date(2002,1) && :Date <= Date(2024,11))
    @select(Not(:D4L_NGDP))
    disallowmissing
end

d4l_data = leftjoin(d4l_quarterly_from_monthly_data, d4l_quarterly_data, on=:Date)


## DLA data (approximated quarter-on-quarter changes by difflog, annualized)

dla_data = @chain quarterly_from_monthly_data begin 
    leftjoin(_, quarterly_data, on=:Date)
    @transform(
        :L_CPI = 100log.(:CPI), 
        :L_CPIXFE = 100log.(:CPIXFE), 
        :L_IPEI = 100log.(:IPEI), 
        :L_S = 100log.(:S), 
        :L_MB = 100log.(:MB),
        :L_CPI_RW = 100log.(:CPI_RW),  
        :L_GDP = 100log.(:GDP),  
        :L_GDP_RW = 100log.(:GDP_RW),  
    )
    @transform(
        :DLA_CPI    = dlafn(:L_CPI), 
        :DLA_CPIXFE = dlafn(:L_CPIXFE), 
        :DLA_IPEI   = dlafn(:L_IPEI), 
        :DLA_S      = dlafn(:L_S), 
        :DLA_MB     = dlafn(:L_MB),
        :DLA_CPI_RW = dlafn(:L_CPI_RW), 
        :DLA_GDP    = dlafn(:L_GDP),
        :DLA_GDP_RW = dlafn(:L_GDP_RW),
    )
    @select(:Date, Cols(r"DLA_"), Cols(r"RS"))
end


## Save preprocessed data

save(datadir("data_QPM.jld2"), "d4l_data", d4l_data, "dla_data", dla_data)

