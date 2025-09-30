### File holds functions needed for workflow analysis

function obj_err(f_x, z)
    return (f_x - z)^2    
end

function model_error(outputs::Vector, obs::Vector)
    return sum(obj_err.(outputs, obs))
end

function calc_err(df_group; obs = zeros(4))
    "Calculates model discrepancy error among
    ensemble members"
    #n = nrow(df_group)
    #data_matrix = Matrix(df_group)
    #errors = Vector{Float64}(undef, n)
    #idx = 1
    #Calculate average simulated output
    avg_output = mean(Matrix(df_group), dims=1)
    #Calculate model discrepancy
    mod_disc = model_error(vec(avg_output), vec(obs))
    #@inbounds for i in 1:n
    #    errors[idx] = model_error(Vector(view(data_matrix, i, :)), Vector(obs))
    #    idx += 1
    #end
    return mod_disc
end

function calc_var(df_group)
    "Calculates variance in errors among
    ensemble members"
    n = nrow(df_group)
    data_matrix = Matrix(df_group)
    errors = Vector{Float64}(undef, binomial(n, 2))
    idx = 1
    
    @inbounds for i in 1:n-1
        for j in i+1:n
            errors[idx] = model_error(Vector(view(data_matrix, i, :)), Vector(view(data_matrix, j, :)))
            idx += 1
        end
    end
    
    return var(errors)
end