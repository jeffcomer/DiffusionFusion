DiffusionFusion
===============

This is a program for calculating position-dependent diffusivity functions from trajectory data. It was written by Jeff Comer <jeffcomer at gmail dot com>. The approach is described in:

Comer, Chipot and Gonzalez-Nilo (2013) J. Chem. Theory Comput., 9(2) 876-882, DOI: 10.1021/ct300867e


See the following for information on two-dimensional calculations and using reflecting boundary conditions, which are useful for an orientation angle:

Comer J, Schulten K, Chipot C (2014) Diffusive models of membrane permeation with explicit orientational freedom. J Chem Theory Comput 10(7):2710–2718. doi: 10.1021/ct500209j


See the following for information on use of the Smoluchowski and fractional Smoluchowski solver:

Chipot C, Comer J (2016) Subdiffusion in membrane permeation of small molecules. Sci. Rep. 6:35913. doi: 10.1038/srep35913


Installation
------------

You need the Gnu Scientific Library development package installed (for random numbers). Here is the Ubuntu package, which may or may not help you:
> sudo apt-get install gsl-bin libgsl0-dev

You may also have to install OpenMP support for g++ (multithreading). Here is the Ubuntu package, which may or may not help you:
> sudo apt-get install libgomp


Then try:
> make

On a Mac, this might work better:
> make mac

Just running diffusion fusion with no parameters will give you the usage information.

> ./diffusionFusion

See the Examples/ directory for more information.


Frequently Asked Questions
===========================

How do I change the lag time? Does it have to do with "timestep" or "hop"?
---------------------------------------------------------------------------

The lag time is not controlled by the options you suggested, but instead by the lag time in the input data. You can modify this by changing the "-stride" option to the "load" commands that take in the input data, or by striding inputs yourself using another program. 

"-hop" is merely an optimization that allows the Smoluchowski solver to calculate over a domain smaller than the entire domain. It is an integer that should be larger than sqrt(2\*largest_expected_diffusivity\*lag_time)/spatial_grid_size. You can always give a number as big as the total number of grid points or larger, but the calculation will be slower without the optimization.

"-timestep" is the time step for the Smoluchowski solver. It just needs to be small enough that the Crank-Nicolson scheme converges. The criterion for damping of oscillations is something like timestep/(spatial_grid_size^2)\*largest_expected_diffusivity < 1/2 ( https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method ). Note that smaller spatial grids require smaller timesteps.


My diffusivity profiles look very noisy for larger strides. Can I reformat the data so that every data point is used even at larger strides? For example, at stride 10, I could include distinct trajectories with frames (0, 10, 20, 30...) and frames (1, 11, 21, 31...).
-------------------------------------------------------------------------------------------------------------------

You could reformat to use every frame, but it would negligibly improve sampling since the trajectories would be highly correlated. For example, the trajectory from frames (0, 10, 20, 30...) and that from frames (1, 11, 21, 31...) will be very similar. It would be much better to obtain more data.

There are other ways to better use the available data. You could make your grid coarser, so that there is more data per diffusivity grid point. The coarseness of the grid is determined by the grid spacing of the input diffusivity grid). For example, you can make change the spacing in the file "coarse_oneHundred.dat":
> field diffuse cubic -step 2 -f coarse_oneHundred.dat -minVal 0.1

Another way would be to more strongly enforce the smoothness of the result by decreasing "-grad" in a smoothness prior. For example, you could reduce this value from 50 Å/ns to 10 Å/ns (the units of this gradient are units(diffusivity)/units(grid_spacing)):
> prior smooth diffuse -grad 10



