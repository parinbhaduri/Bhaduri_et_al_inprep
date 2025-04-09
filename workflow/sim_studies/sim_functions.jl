


function simul_plot(adf, group_col; leg = :outertopright, color = palette(:BrBG_11))
    plots = []
    push!(plots, plot(adf.time, adf.sum_hh_low_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Low Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_low_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Low Income Residents"))

    push!(plots, plot(adf.time, adf.sum_hh_med_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Medium Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_med_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Medium Income Residents"))

    push!(plots, plot(adf.time, adf.sum_hh_high_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="High Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_high_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="High Income Residents"))
    
    return plots
end

function simul_market(adf, mdf, group_col; leg = :outertopright, color = palette(:BrBG_11))
    plots = []
    #Ignore time=0 for moving plots
    move_ind = mdf.time .!= 0 
    push!(plots, plot(mdf.time[move_ind], mdf.moved_low[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, title="Low Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_low_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Low Income House Price"))
    
    push!(plots, plot(mdf.time[move_ind], mdf.moved_med[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, title="Medium Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_med_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="Medium Income House Price"))
    
    push!(plots, plot(mdf.time[move_ind], mdf.moved_high[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, title="High Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_high_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, title="High Income House Price"))
    
    
    return plots
end