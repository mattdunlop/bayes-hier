function l = model_id(S)

    % Number of observation points
    J = 100;

    % Find appropriate entries of S and return them
    N = size(S,1);
    step = round(N/(sqrt(J)+1));

    obs_base = step:step:sqrt(J)*step;
    obs = obs_base + (step-1)*N;

    for p=2:sqrt(J);
        obs = [obs,obs_base + (step*p-1)*N];
    end

    l = S(obs)';

end