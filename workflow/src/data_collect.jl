
##Calculating Population Characteristics
BG(agent) = agent isa BlockGroup
HH(agent) = agent isa HHAgent && agent.bg_id >= 1
#Pop. Category Pop
hh_low(agent) = agent.group == 1 ? 1 : 0
hh_med(agent) = agent.group == 2 ? 1 : 0
hh_high(agent) = agent.group == 3 ? 1 : 0
#Occupancy Pop
occ_low(agent) = agent.occupied_units[1]
occ_med(agent) = agent.occupied_units[2]
occ_high(agent) = agent.occupied_units[3]

##Calculating Transaction Characteristics
function moved(model, cat)
    count = 0
    cat_ids = [id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[id].group == cat]
    for id in cat_ids
        count += model[id].year_of_residence == model.tick ? 1.0 : 0.0
    end
    return count  
end

moved_low(model) = moved(model, 1)
moved_med(model) = moved(model, 2)
moved_high(model) = moved(model, 3)

price_low(agent) = agent.new_price[1]
price_med(agent) = agent.new_price[2]
price_high(agent) = agent.new_price[3]

simul_adata = [(hh_low, sum, HH), (hh_med, sum, HH), (hh_high, sum, HH) , (occ_low, sum, BG), (occ_med, sum, BG), (occ_high, sum, BG), 
    (price_low, mean, BG), (price_med, mean, BG), (price_high, mean, BG)
]

simul_mdata = [moved_low, moved_med, moved_high]
