% Some example code for providing figure output.
%
% If the software for the forward map (e.g. EIDORS) has its own plotting 
% tools, it may be preferable to use those.

% Define the output grids, and output the true image
if s == output.thinning
	[X,Y] = meshgrid(1/(2*N):1/N:1-1/(2*N),1/(2*N):1/N:1-1/(2*N));
	[XT,YT] = meshgrid(1/(2*NT):1/NT:1-1/(2*NT),1/(2*NT):1/NT:1-1/(2*NT));
	
	subplot(221);
	surf(XT,YT,UT_phys,'EdgeColor','None');view(2);axis square;
	axis off;
	title('True field');
end

% Output the current sample
subplot(222);
hh = surf(X,Y,make_lvl(idct2(reshape(U,N,N)),tau,prior.U.alpha),'EdgeColor','None');view(2);axis square;
axis off;
title('Current sample');

% Output the trace of tau
subplot(224);
plot(traceTau(1:s));
if s > 1
	xlim([1 s]);
end
axis square;
title('Inverse length scale');

pause(0.01);