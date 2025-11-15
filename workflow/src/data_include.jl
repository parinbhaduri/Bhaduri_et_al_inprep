###import input data 
data_location = "philadelphia-data/model_inputs"
bg_file = "phil_flood_bg_2019_nomiss_v1.csv"
base_pop_file = "philly_cbsa_pop_0.csv"
flood_file = "phil_flood_hist_year.csv"

pop_dir = joinpath(dirname(pwd()), data_location, "pop_files")

##Read in Demographic Data
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, bg_file)))

##Read in flood data
phil_flood = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, flood_file)))
#transform df to correct format
phil_flood_record = unstack(phil_flood, :GEOID, :year, :perc_flood_extent)
#Extra edits
phil_flood_record[!,"1982"] = zeros(size(phil_flood_record)[1])
select!(phil_flood_record, "GEOID", "1981", "1982", Not(["1982", "2019"]), "2019")

#Create synthetic flood record of no floods 
synth_mat = zeros(size(phil_flood_record[:,2:end]))
synth_flood_record = DataFrame(hcat(phil_flood_record.GEOID,synth_mat), names(phil_flood_record))
synth_flood_record[!,:GEOID] = convert.(Int64, synth_flood_record[!,:GEOID])


