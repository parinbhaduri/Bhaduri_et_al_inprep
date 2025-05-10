


function simul_plot(adf, group_col; leg = :outertopright, color = palette(:BrBg), lim = (7000,15000))
    plots = []
    push!(plots, plot(adf.time, adf.sum_hh_low_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Low Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_low_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Low Income Residents"))

    push!(plots, plot(adf.time, adf.sum_hh_med_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Medium Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_med_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Medium Income Residents"))

    push!(plots, plot(adf.time, adf.sum_hh_high_HH, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="High Income Pop."))
    push!(plots, plot(adf.time, adf.sum_occ_high_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=lim, title="High Income Residents"))
    
    return plots
end

function simul_market(adf, mdf, group_col; leg = :outertopright, color = palette(:BrBg), lim =(0,600), price_lim =(1e5,1.5e6))
    plots = []
    #Ignore time=0 for moving plots
    move_ind = mdf.time .!= 0 
    push!(plots, plot(mdf.time[move_ind], mdf.moved_low[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Low Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_low_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=price_lim, title="Low Income House Price"))
    
    push!(plots, plot(mdf.time[move_ind], mdf.moved_med[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Medium Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_med_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=price_lim, title="Medium Income House Price"))
    
    push!(plots, plot(mdf.time[move_ind], mdf.moved_high[move_ind], group = mdf[:, group_col][move_ind], legend=leg, palette = color, linewidth=2, ylimits=lim, title="High Income Movers"))
    push!(plots, plot(adf.time, adf.mean_price_high_BG, group = adf[:, group_col], legend=leg, palette = color, linewidth=2, ylimits=price_lim, title="High Income House Price"))
    
    
    return plots
end

##Repeat plot functions for utility coefficients
function util_plot(adf_low,adf_med,adf_high, group_col; leg = :outertopright, color = palette(:BrBg), lims = (8000,15000))
    plots = []
    #Plot Area results
    push!(plots, plot(adf_low.time, adf_low.sum_hh_low_HH, group = select(adf_low,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="Low Income Pop."))
    push!(plots, plot(adf_low.time, adf_low.sum_occ_low_BG, group = select(adf_low,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="Low Income Residents"))

    push!(plots, plot(adf_med.time, adf_med.sum_hh_med_HH, group = select(adf_med,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="Medium Income Pop."))
    push!(plots, plot(adf_med.time, adf_med.sum_occ_med_BG, group = select(adf_med,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="Medium Income Residents"))

    push!(plots, plot(adf_high.time, adf_high.sum_hh_high_HH, group = select(adf_high,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="High Income Pop."))
    push!(plots, plot(adf_high.time, adf_high.sum_occ_high_BG, group = select(adf_high,Regex("$(group_col)"))[:,1], legend=:outertopright, palette = color, linewidth=2, ylimits=lims, title="High Income Residents"))
    
    return plots
end


function util_market(adf_low,adf_med,adf_high, mdf_low,mdf_med,mdf_high, group_col; leg = :outertopright, color = palette(:BrBg), lims = (0,800))
    plots = []
    #Ignore time=0 for moving plots
    move_ind = mdf_low.time .!= 0 
    push!(plots, plot(mdf_low.time[move_ind], mdf_low.moved_low[move_ind], group = select(mdf_low,Regex("$(group_col)"))[move_ind,1], legend=leg, palette = color, linewidth=2, ylimits=lims, title="Low Income Movers"))
    push!(plots, plot(adf_low.time, adf_low.mean_price_low_BG, group = select(adf_low,Regex("$(group_col)"))[:,1], legend=leg, palette = color, linewidth=2, ylimits=(1e5,1.5e6), title="Low Income House Price"))
    
    push!(plots, plot(mdf_med.time[move_ind], mdf_med.moved_med[move_ind], group = select(mdf_med,Regex("$(group_col)"))[move_ind,1], legend=leg, palette = color, linewidth=2, ylimits=lims, title="Medium Income Movers"))
    push!(plots, plot(adf_med.time, adf_med.mean_price_med_BG, group = select(adf_med,Regex("$(group_col)"))[:,1], legend=leg, palette = color, linewidth=2, ylimits=(1e5,1.5e6), title="Medium Income House Price"))
    
    push!(plots, plot(mdf_high.time[move_ind], mdf_high.moved_high[move_ind], group = select(mdf_high,Regex("$(group_col)"))[move_ind,1], legend=leg, palette = color, linewidth=2, ylimits=lims, title="High Income Movers"))
    push!(plots, plot(adf_high.time, adf_high.mean_price_high_BG, group = select(adf_high,Regex("$(group_col)"))[:,1], legend=leg, palette = color, linewidth=2, ylimits=(1e5,1.5e6), title="High Income House Price"))
    
    
    return plots
end