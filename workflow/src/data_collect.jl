
##Calculating Population Characteristics
BG(agent) = agent isa BlockGroup
FL_BG(agent) = agent isa BlockGroup && agent.perc_fld_area > 0
HH(agent) = agent isa HHAgent && agent.bg_id >= 1
#Pop. Category Pop
hh_low(agent) = agent.group == 1 ? 1 : 0
hh_med(agent) = agent.group == 2 ? 1 : 0
hh_high(agent) = agent.group == 3 ? 1 : 0

#Pop. Category Pop
hh_diff_low(agent) = agent.group == 1 && agent.occ_cat != 1 ? 1 : 0 #Low Income Agents not living in Low Income Housing
hh_diff_med(agent) = agent.group == 2 && agent.occ_cat != 2 ? 1 : 0
hh_diff_high(agent) = agent.group == 3 && agent.occ_cat != 3 ? 1 : 0
#Occupancy Pop
occ_low(agent) = agent.occupied_units[1]
occ_med(agent) = agent.occupied_units[2]
occ_high(agent) = agent.occupied_units[3]
#Housing Prices
price_low(agent) = agent.new_price[1]
price_med(agent) = agent.new_price[2]
price_high(agent) = agent.new_price[3]
#Housing Capacities
cap_low(agent) = agent.occupied_units[1] + agent.available_units[1]
cap_med(agent) = agent.occupied_units[2] + agent.available_units[2]
cap_high(agent) = agent.occupied_units[3] + agent.available_units[3]


##Calculating Transaction Characteristics
function moved(model, cat)
    count = 0
    cat_ids = [id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[id].occ_cat == cat]
    for id in cat_ids
        count += model[id].year_of_residence == model.tick ? 1.0 : 0.0
    end
    return count  
end
moved_low(model) = moved(model, 1)
moved_med(model) = moved(model, 2)
moved_high(model) = moved(model, 3)

#Calculate agents moving into floodplain at each timestep
function flpn_moved(model, cat)
    count = 0
    cat_ids = [id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[id].occ_cat == cat && model[model[id].bg_id].perc_fld_area > 0]
    for id in cat_ids
        count += model[id].year_of_residence == model.tick ? 1.0 : 0.0
    end
    return count  
end
flpn_moved_low(model) = flpn_moved(model, 1)
flpn_moved_med(model) = flpn_moved(model, 2)
flpn_moved_high(model) = flpn_moved(model, 3)

#Calculate agents moving to different categories at each time step 
function moved_diff_cat(model, cat)
    count = 0
    cat_ids = [id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[id].occ_cat == cat && model[id].group != cat]
    for id in cat_ids
        count += model[id].year_of_residence == model.tick ? 1.0 : 0.0
    end
    return count  
end
moved_diff_low(model) = moved_diff_cat(model, 1)
moved_diff_med(model) = moved_diff_cat(model, 2)
moved_diff_high(model) = moved_diff_cat(model, 3)

##Calculating Population in floodplain
function flpn_pop(model, cat)
    return length([id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[model[id].bg_id].perc_fld_area > 0 && model[id].group == cat])
end

flpn_pop_low(model) = flpn_pop(model, 1)
flpn_pop_med(model) = flpn_pop(model, 2)
flpn_pop_high(model) = flpn_pop(model, 3)

##Calculating Population exposed to flooding at each time step
function pop_exposed(model, cat)
    return length([id for id in allids(model) if model[id] isa HHAgent && model[id].bg_id >= 1 && model[id].group == cat && model[id].flood_hazard[model.tick] > 0])
end

pop_exposed_low(model) = pop_exposed(model, 1)
pop_exposed_med(model) = pop_exposed(model, 2)
pop_exposed_high(model) = pop_exposed(model, 3)


#Simulation Studies data collection
simul_adata = [(hh_low, sum, HH), (hh_med, sum, HH), (hh_high, sum, HH), (hh_diff_low, sum, HH), (hh_diff_med, sum, HH), (hh_diff_high, sum, HH), 
    (occ_low, sum, BG), (occ_med, sum, BG), (occ_high, sum, BG), (price_low, mean, BG), (price_med, mean, BG), (price_high, mean, BG), (:flood_hazard,mean,FL_BG)
]
simul_mdata = [moved_low, moved_med, moved_high, moved_diff_low, moved_diff_med, moved_diff_high]

#Simulation Studies data collection for Floodplain
flpn_adata = [(occ_low, sum, FL_BG), (occ_med, sum, FL_BG), (occ_high, sum, FL_BG), 
    (price_low, mean, FL_BG), (price_med, mean, FL_BG), (price_high, mean, FL_BG)
]
flpn_mdata = [flpn_pop_low, flpn_pop_med, flpn_pop_high, flpn_moved_low, flpn_moved_med, flpn_moved_high]

#Calibration data collection
calib_adata = [(hh_low, sum, HH), (hh_med, sum, HH), (hh_high, sum, HH), (price_low, mean, BG), (price_med, mean, BG), (price_high, mean, BG)]
calib_mdata = [moved_low, moved_med, moved_high, flpn_pop_low, flpn_pop_med, flpn_pop_high]

#Shapley data collection
shap_adata = [:pos, :GEOID, hh_low, hh_med, hh_high, price_low, price_med, price_high, cap_low, cap_med, cap_high]
