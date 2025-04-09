### File holds functions needed for workflow analysis

function obj_err(f_x, z)
    return (f_x - z)^2    
end

function model_error(outputs::Vector, obs::Vector)
    return sum(obj_err.(outputs, obs))
end