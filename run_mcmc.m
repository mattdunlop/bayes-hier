%% Define the parameters

% Model choice
model = 'id';                   % Choice of model: 'id', 'gwf' or 'eit'

% Sample resolution
N = 80;                         % Sampling performed on NxN square grid
NT = 320;						% Data generated on NTxNT square grid

% Prior parameters
prior.U.alpha = 4;				% (Conditional) field prior: N(0,(-Delta + tau^2)^(-alpha))
prior.tau.mean = 20;			% Lengthscale prior: N(mean,std^2)
prior.tau.std = 10;      

% Data parameters
data.gamma = 0.3;				% Noise standard deviation

% MCMC parameters
mcmc.beta = 0.05;				% pCN jump for U
mcmc.eps = 2;					% RWM jump for tau
mcmc.samples = 5000000;			% Number of samples to generate

% Output parameters
output.figures = 1;				% Display figures?
output.thinning = 100;			% How often to update output?

%% Create the data (or alternatively import from a file)

% Define the true geometric field
tauT = 15;
UT = gaussrnd(prior.U.alpha,tauT,NT);
UT_phys = make_lvl(idct2(reshape(UT,NT,NT)),tauT,prior.U.alpha);

% Create the data
data.clean = ell(UT_phys,model);
data.J = length(data.clean);
eta = normrnd(0,1,data.J,1);
data.noisy = data.clean + data.gamma*eta;


%% Define the potential
Phi = @(U,tau) norm((ell(make_lvl(idct2(reshape(U,N,N)),tau,prior.U.alpha),model) - data.noisy)./data.gamma)^2/2;


%% Run the MCMC

% Define the initial MCMC state
tau = 60;
U = gaussrnd(prior.U.alpha,tau,N);

phiU = Phi(U,tau);

% Track the acceptance rates for the proposals
mcmc.accepted.field = 0;
mcmc.accepted.tau = 0;

% Store the traces of tau and the lowest frequency modes
traceTau = zeros(mcmc.samples,1);
traceU = zeros(mcmc.samples,8);

% Perform the MCMC
for s=1:mcmc.samples
    traceTau(s) = tau;
	traceU(s,:) = [U(2),U(3),U(N+1),U(N+2),U(N+3),U(2*N+1),U(2*N+2),U(2*N+3)];
	
	% Provide output every output.thinning samples
	if mod(s,output.thinning) == 0
		fprintf('s = %i\n',s);
		fprintf('Acceptance rate for U:\t\t%f\n',mcmc.accepted.field/s);
		fprintf('Acceptance rate for tau:\t%f\n\n',mcmc.accepted.tau/s);
	
		if output.figures == 1
			display_figures;
		end
	end
	
    %% Update U | tau,y
	
	% Generate the pCN proposal
	UJump = gaussrnd(prior.U.alpha,tau,N);
	UProp = sqrt(1-mcmc.beta^2)*U + mcmc.beta*UJump;

	% Calculate the acceptance probability aProb
	phiUProp = Phi(UProp,tau);
	aProb = min(1,exp(phiU - phiUProp));

	% Accept proposal with probability aProb
	if unifrnd(0,1) < aProb
		U = UProp;
		phiU = phiUProp;
		mcmc.accepted.field = mcmc.accepted.field + 1;
    end
    
    %% Update tau | U,y

	% Generate the RWM proposal
    tauJump = normrnd(0,1);
    tauProp = tau + mcmc.eps*tauJump;
	
	% Calculate the acceptance probability aProb
	phiUProp = Phi(U,tauProp);
    
    [K1,K2] = meshgrid(0:N-1,0:N-1);
    eigL = reshape((pi^2*(K1.^2+K2.^2) + tau^2).^prior.U.alpha,N^2,1);
    eigLProp = reshape((pi^2*(K1.^2+K2.^2) + tauProp^2).^prior.U.alpha,N^2,1);
	summand = (eigL - eigLProp).*U.^2./N^2 + log(eigLProp./eigL);
    
    summand(1) = 0;
	propSum = sum(summand);
    aProb = min(1,exp(phiU - phiUProp + ((tau-prior.tau.mean)^2 - (tauProp - prior.tau.mean)^2)/(2*prior.tau.std^2) + propSum/2));
    
	% Accept proposal with probability aProb
    if unifrnd(0,1) < aProb
		tau = tauProp;
		phiU = phiUProp;
		mcmc.accepted.tau = mcmc.accepted.tau + 1;
	end	
	
end


































