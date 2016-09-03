MATLAB code to run dimension robust MCMC for hierarchical Bayesian inversion, as outlined 
in the paper `Hierarchical Bayesian Level Set Inversion' by Dunlop, Iglesias and Stuart. 
Three example forward models are provided: direct point observations, a groundwater flow 
model and an electrical impedance tomography model.

The following files are provided:

run_mcmc.m
	Perform the MCMC, given the parameters (grid resolution, number of samples, prior 
	smoothness, etc) defined at the start of the file.
	Optionally the output may be displayed as a figure. Saved in memory are the traces of 
	the length scale parameter and a few Fourier modes.

gaussrnd.m
	Generate a 2D sample from the Gaussian prior outlined in the paper, given a smoothness 
	parameter alpha, inverse length scale parameter tau, and grid size N.
	The sample is provided in Fourier space, reshaped into an N^2*1 vector.

make_lvl.m
	Take a square matrix representing a continuous function, and threshold at two levels, 
	returning a matrix representing a piecewise constant function.
	The values that the thresholded function takes are defined in this file.

ell.m
	Select the appropriate forward model, mapping the thresholded function to the output 
	measurements. The three models are contained in model_id.m, model_gwf.m and 
	model_eit.m.

model_id.m
	Perform direct observations of the piecewise constant field. Observations are 
	performed on a square grid of J points. J is defined in the file, and should be a
	perfect square.

model_gwf.m
	Solve using the steady state Darcy flow model, as outlined in the paper. Observations
	are performed on a square grid of J points. J is defined in the file, and should be a
	perfect square.

model_eit.m
	NOTE: Requires EIDORS software (http://eidors3d.sourceforge.net/). The path to the 
	software must be defined at the start of the file.
	
	EIT forward model, using the complete electrode model. In the example provided, a 
	a circular domain is used with 16 electrodes equally spaced around its boundary, and 
	the maximum 15 linear independent current stimulation patterns are applied.

display_figures.m
	Provide some visual output for the MCMC: the true field that generates the data, the 
	current MCMC sample, and the trace of the length scale parameter.
	Note that, for the EIT example, output of the full sample on a square domain is shown, 
	even though the simulation domain is circular. The code may be modified to use EIDORS 
	built in image output tools (as was done in the paper).
	
	