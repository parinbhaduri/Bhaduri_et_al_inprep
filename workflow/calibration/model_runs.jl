#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


### PARALLEL ENSEMBLE RUN ###

calib_params = Dict(
    :risk_averse=>collect(range(0.1,0.9,step=0.2)),
    :flood_mem=>10,#[5,10,15],
    :build_perc=>[0.05,0.1,0.25, 0.4],
    :price_perc=>[0.1,0.15,0.2],
    :perc_growth=>collect(range(0.0, 0.05, step = 0.01)),
    :base_move=>collect(range(0.0, 0.05, step = 0.005)),
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>collect(range(50000,500000,step=5000)), 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :penalty=>push!(collect(range(0.0,100000.0,step=10000.0)), 10000000.0),
    :flood_coefficient=>collect(range(100000.0,1000000.0,step=50000.0)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500 
)

#Calculate number of parameter combinations
length(Iterators.product(values(calib_params)...))

count = 0
for iter in Iterators.product(values(calib_params)...)
    if count == 5
        break
    else
        println(iter)
        count += 1
    end
end


#Run Models. Collect Data
adf_calib,mdf_calib = paramscan(calib_params, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

#Save Data as df
CSV.write(joinpath(@__DIR__,"data/model_run_adf.csv"), adf_calib)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf.csv"), mdf_calib)