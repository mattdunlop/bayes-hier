% Define the forward map, mapping the geometric field into the observation 
% space.

function l = ell(S,model)	

    % Identity model
    if strcmp(model,'id')
        l = model_id(S);
        
    % Groundwater flow model
    elseif strcmp(model,'gwf')
        l = model_gwf(S);        
        
    % EIT model
    % Note: requires EIDORS software
    elseif strcmp(model,'eit') 
        l = model_eit(S);
  
    else
        error('Invalid model specified');
    end
end