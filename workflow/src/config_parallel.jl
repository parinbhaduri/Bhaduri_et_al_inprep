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

@everywhere include("data_collect.jl")
@everywhere include("config.jl")