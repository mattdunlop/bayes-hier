% The construction map taking a level set field U to a piecewise constant
% geometric field V, with scaling of relative levels by the lengthscale
% tau.
% The example here is for three phases/two thresholding levels.

function V = make_lvl(U,tau,alpha)
	
	d = 2;								% Spatial dimension
	c = 0.1*tau^(d/2-alpha)*[1,-1];		% The thresholding levels, scaled by tau

	V = U;
    V(find(U>=c(1))) = 10;				% Value of phase 1
	V(find(c(1)>U & U>=c(2))) = 5;		% Value of phase 2
    V(find(c(2)>U)) = 1;				% Value of phase 3
end