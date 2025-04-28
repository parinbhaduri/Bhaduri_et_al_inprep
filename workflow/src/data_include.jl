###import input data 
data_location = "philadelphia-data/model_inputs"
bg_file = "phil_flood_bg_2019_v2.csv"
pop_file = "philly_cbsa_pop_0.csv"
flood_file = "phil_flood_hist_year.csv"

##Read in Demographic Data
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, bg_file)))
phil_cbsa_base_pop = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, "pop_files", pop_file)))
#drop missing values
dropmissing!(phil_cbsa_base_pop, :NP)
#For rows with people and negative income, set income to bottom 10%
inc_bot_10 = quantile(subset(phil_cbsa_base_pop, [:NP .=> ByRow(>(0)), :adj_income_2019 .=> ByRow(>(0))]).adj_income_2019, [0.10])[1]
@. phil_cbsa_base_pop.adj_income_2019 = ifelse.(phil_cbsa_base_pop.NP > 0 && phil_cbsa_base_pop.adj_income_2019 <= 0, inc_bot_10, phil_cbsa_base_pop.adj_income_2019)

##Read in flood data
phil_flood = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, flood_file)))
#transform df to correct format
phil_flood_record = unstack(phil_flood, :GEOID, :year, :perc_flood_extent)
#Extra edits
phil_flood_record[!,"1982"] = zeros(size(phil_flood_record)[1])
select!(phil_flood_record, "GEOID", "1981", "1982", Not(["1982", "2019"]), "2019")

synth_flood_record = DataFrame(CSV.File(joinpath(dirname(pwd()), "CHANCE_C.jl","data","synth_flood_phil.csv")))