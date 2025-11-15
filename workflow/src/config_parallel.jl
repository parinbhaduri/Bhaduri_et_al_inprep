using Distributed
addprocs(12, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
end

@everywhere begin
    include("data_include.jl")
    include("functions.jl")
    include("data_collect.jl")
    
    phil_cbsa_base_pop = load_pop(0)
end 
