using CSV
using Distributions
using DataFrames
using LinearAlgebra
using Measures
using Plots
using Random
using SpecialFunctions
using StatsPlots
using StatsBase
using ProgressBars
include("heavytail.jl")
include("outbreaklikelihood.jl")
gr(fontfamily="Helvetica",foreground_color_legend = nothing,background_color_legend = nothing, titlefontsize=11, tickfontsize=10, legendfontsize=10,labelfontsize=10,grid=true, tick_direction=:out,size=(400,300),fmt = :svg)
    
# Util function
Base.:-(x::Tuple{Float64,Float64,Float64},y::Int)=x.-y
TruncatedNormal(mu, sigma, l, u)=Truncated(Normal(mu, sigma), l, u)

# Parameters
## Weibull parameters
if !isdefined(@__MODULE__,:κ_msm) κ_msm=0.6 end # Pareto exponent
if !isdefined(@__MODULE__,:α_msm) α_msm=0.15 end
#if !isdefined(@__MODULE__,:κ_other) κ_other=1.66 end # Pareto exponent
#if !isdefined(@__MODULE__,:κ_other) α_other = 0.04 end
## Epi parameters
if !isdefined(@__MODULE__,:inf_period) inf_period=21 end
## Demography parameters
if !isdefined(@__MODULE__,:p_msm) p_msm = 0.02 end# proportion of sexually active MSM
#if !isdefined(@__MODULE__,:p_other) p_other = 0.86 end # proportion of sexually active non-MSM
if !isdefined(@__MODULE__,:p_hem) p_hem = 0.86 end # proportion of sexually active heterosexual men
if !isdefined(@__MODULE__,:p_hew) p_hew = 0.89 end # proportion of sexually active heterosexual women

## outbreak size specification
if !isdefined(@__MODULE__,:currentcases) currentcases = 1000 end

## crontrol variables
if !isdefined(@__MODULE__,:printtable) printtable = true end
if !isdefined(@__MODULE__,:outbreak_iters) outbreak_iters = 100000 end
if !isdefined(@__MODULE__,:hist_upp) hist_upp = :auto end

# Degree distribution over 1 year
deg1y_msm = truncated(Weibull2(α_msm, κ_msm), lower=1)
#deg1y_other = truncated(Weibull2(α_other, κ_other), lower=1)
deg1y_hem = truncated(Weibull2(α_hem, κ_hem), lower=1)
deg1y_hew = truncated(Weibull2(α_hew, κ_hew), lower=1)

# Preparation for draws from ~ x*pdf(x)
ED_3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull 
ED_3w_hem = truncated(Weibull(params(deg1y_hem)[1], params(deg1y_hem)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull
ED_3w_hew = truncated(Weibull(params(deg1y_hew)[1], params(deg1y_hew)[2]*inf_period/365), lower=inf_period/365) |> EDTWeibull

# Pre-draws from excess degree +1 distribution for efficiency
Random.seed!(2022)
drawsize = 1000000
edeg3w_msm = rand(ED_3w_msm, drawsize)
#edeg3w_other = rand(ED_3w_other, drawsize)
edeg3w_hem = rand(ED_3w_hem, drawsize)
edeg3w_hew = rand(ED_3w_hew, drawsize)

# Degree distributions for different types
deg3w_msm = truncated(Weibull(params(deg1y_msm)[1], params(deg1y_msm)[2]*inf_period/365), lower=inf_period/365)
#deg3w_other = truncated(Weibull(params(deg1y_other)[1], params(deg1y_other)[2]*inf_period/365), lower=inf_period/365)
deg3w_hem = truncated(Weibull(params(deg1y_hem)[1], params(deg1y_hem)[2]*inf_period/365), lower=inf_period/365)
deg3w_hew = truncated(Weibull(params(deg1y_hew)[1], params(deg1y_hew)[2]*inf_period/365), lower=inf_period/365)

contDirac = TruncatedNormal(0,1e-10,0,0.1)
deg3w_msmavg = MixtureModel([contDirac,deg3w_msm],[1-p_msm ,p_msm])
deg3w_hemavg = MixtureModel([contDirac,deg3w_hem],[1-p_hem ,p_hem])
deg3w_hewavg = MixtureModel([contDirac,deg3w_hew],[1-p_hew ,p_hew])


### Outputs
if !isdefined(@__MODULE__,:R0only) R0only = false end
if !isdefined(@__MODULE__,:R0finite) R0finite = false end
if !isdefined(@__MODULE__,:R0deplete) R0deplete = false end
if !isdefined(@__MODULE__,:R0assortative) R0assortative = false end


# Likelihood of an outbreak of current size or greater
if printtable
for SAR in [0.05, 0.1, 0.2, 0.4]
    outbreaklikelihood1(edeg3w_msm, edeg3w_hem, edeg3w_hew, currentcases, SAR, iter=outbreak_iters)
    print("\n\n")
    outbreaklikelihood2(edeg3w_msm, edeg3w_hem, edeg3w_hew, deg3w_msm, deg3w_msmavg, deg3w_hemavg,deg3w_hewavg, currentcases, SAR, iter=outbreak_iters)
    print("\n\n\n")
end
end

