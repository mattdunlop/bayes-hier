function l = model_eit(S)
    % Load EIDORS if not already loaded
    if exist('show_fem','file') == 0
        run('/path/to/eidors/startup.m');
    end

    % Forward model parameters
    model_type = 'm2C';				% Choice of predefined 2D model (help mk_common_model)
    num_elec = 16;					% Number of electrodes
    num_stim = num_elec - 1;		% Number of lienarly independent stimulation patterns
    inj_pattern = [0,1];			% Injection pattern (help mk_stim_patterns)
    meas_pattern = '{mono}';		% Measurement pattern (help mk_stim_patterns)
    stim_current = 0.1;				% Drive current levels (Amps)
    contact_imp = 0.01;				% Contact impedences (Ohms)

    % Create the 2D model
    imdl= mk_common_model(model_type,num_elec);

    % Create homogeneous image
    img = mk_image(imdl);
 
    % Create stimulation patterns
    [stim,~] = mk_stim_patterns(num_elec,1,inj_pattern,meas_pattern,{'meas_current','balance_inj','balance_meas'},stim_current);

    % Set the contact impedences
    for j=1:num_elec
        img.fwd_model.electrode(j).z_contact = contact_imp;
        nodej = img.fwd_model.electrode(j).nodes;

        % Adjust electrode size (affected by mesh fidelity)
        img.fwd_model.electrode(j).nodes = [nodej,nodej+2,nodej+4]; 
    end
    
    img.fwd_model.stimulation = stim;
    img.fwd_solve.get_all_meas = 1;

    N = size(S,1);
    [X,Y] = meshgrid(-1+1/N:2/N:1-1/N,-1+1/N:2/N:1-1/N);
    c_field = @(x,y,z) interp2(X,Y,S,x,y,'spline');
    img.elem_data = elem_select(img.fwd_model,c_field);

    v = fwd_solve(img);	
    l = v.meas(1:num_stim*num_elec);
end