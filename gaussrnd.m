% Return a sample of a Gaussian random field on [0,1]^2 with: 
%       mean 0
%       covariance operator C = (-Delta + tau^2)^(-alpha)
% where Delta is the Laplacian with zero Neumann boundary conditions.
% Returned is an N^2x1 vector of Fourier cosine coefficients, so that
% U = idct2(reshape(L,N,N)) is the sample in physical space.

function L = gaussrnd(alpha,tau,N)
	
	% Random variables in KL expansion
	xi = normrnd(0,1,N);
	
	% Define the (square root of) eigenvalues of the covariance operator
	[K1,K2] = meshgrid(0:N-1,0:N-1);
	coef = (pi^2*(K1.^2+K2.^2) + tau^2).^(-alpha/2);	
	
	% Construct the KL coefficients
	L = N*coef.*xi;
    L(1,1) = 0;
	
    L = reshape(L,N^2,1);
    
end