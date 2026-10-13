# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Julia 1.12
#     language: julia
#     name: julia-1.12
# ---

# %%
using Pkg
Pkg.activate(".")
Pkg.instantiate() 

# %%
params_main = (
    κ_msm=0.77, α_msm=0.10,  # 2000–2020
    κ_hem=1.68, α_hem=0.01,  # 2000–2010
    κ_hew=2.14, α_hew=0.01   # 2010
)

# %%
module Baseline
κ_msm, α_msm, κ_hem, α_hem, κ_hew, α_hew = 
    getindex.(Ref(Main.params_main), [:κ_msm, :α_msm, :κ_hem, :α_hem, :κ_hew, :α_hew])

printtable=true
outbreak_iters=100000
include("src/outputs.jl")
end

# %%
module InfPeriodSAR_genpop
using Printf

κ_msm, α_msm, κ_hem, α_hem, κ_hew, α_hew = 
    getindex.(Ref(Main.params_main), [:κ_msm, :α_msm, :κ_hem, :α_hem, :κ_hew, :α_hew])

currentcases = 1000
printtable = false
outbreak_iters = 100000
inf_period = 0

sarlist = [0.05, 0.1, 0.2, 0.4]
results = Dict{Int, Vector{Float64}}()

for inf_p in [3, 7, 14, 21]
    global inf_period = inf_p
    include("src/outputs.jl") 
    
    probs = Float64[]
    for SAR in sarlist
        prob = simulatetillsize(edeg3w_msm, currentcases, SAR;
                                init_deg=deg3w_msmavg, iter=outbreak_iters)
        push!(probs, prob)
    end
    results[inf_p] = probs
end

@printf("%-12s", "inf_period")
for SAR in sarlist; @printf(" | SAR=%-4.0f%%", SAR*100) end
println()
println("-"^55)
for inf_p in [3, 7, 14, 21]
    @printf("%-12s", "$(inf_p) days")
    for p in results[inf_p]; @printf(" | %7.4f%%", p*100) end
    println()
end
end

# %%
module InfPeriodSAR_genpop_n
using Printf

sarlist = [0.05, 0.1, 0.2, 0.4]
results = Main.InfPeriodSAR_genpop.results

for inf_p in [3, 7, 14, 21]
    println("\n=== inf_period = $inf_p days ===")
    @printf("%-10s", "n cases")
    for SAR in sarlist; @printf(" | SAR=%-4.0f%%", SAR*100) end
    println()
    println("-"^55)
    for n in [100, 150, 200]
        @printf("%-10s", "n=$n")
        for p in results[inf_p]
            q = 1 - (1 - p)^n
            @printf(" | %7.4f%%", q*100)
        end
        println()
    end
end
end
