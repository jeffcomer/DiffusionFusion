DiffusionFusion
===============

This is a program for calculating position-dependent diffusivity functions from trajectory data. It was written by Jeff Comer <jeffcomer at gmail dot com>. The approach is described in:

Comer, Chipot and Gonzalez-Nilo (2013) J. Chem. Theory Comput., 9(2) 876-882, DOI: https://doi.org/10.1021/ct300867e


See the following for information on two-dimensional calculations and using reflecting boundary conditions, which are useful for an orientation angle:

Comer J, Schulten K, Chipot C (2014) Diffusive models of membrane permeation with explicit orientational freedom. J Chem Theory Comput 10(7):2710–2718. doi: https://doi.org/10.1021/ct500209j


See the following for information on use of the Smoluchowski and fractional Smoluchowski solver:

Chipot C, Comer J (2016) Subdiffusion in membrane permeation of small molecules. Sci. Rep. 6:35913. doi: https://doi.org/10.1038/srep35913


Apologies to Phyllis Nicolson ( https://en.wikipedia.org/wiki/Phyllis_Nicolson ) for misspelling her name in the source files. 


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
--------------------------

**How do I change the lag time? Does it have to do with "timestep" or "hop"?**

The lag time is not controlled by the options you suggested, but instead by the lag time in the input data. You can modify this by changing the "-stride" option to the "load" commands that take in the input data, or by striding inputs yourself using another program. 

"-hop" is merely an optimization that allows the Smoluchowski solver to calculate over a domain smaller than the entire domain. It is an integer that should be larger than sqrt(2\*largest_expected_diffusivity\*lag_time)/spatial_grid_size. You can always give a number as big as the total number of grid points or larger, but the calculation will be slower without the optimization.

"-timestep" is the time step for the Smoluchowski solver. It just needs to be small enough that the Crank-Nicolson scheme converges. The criterion for damping of oscillations is something like timestep/(spatial_grid_size^2)\*largest_expected_diffusivity < 1/2 ( https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method ). Note that smaller spatial grids require smaller timesteps.

&nbsp;

**My diffusivity profiles look very noisy for larger strides. Can I reformat the data so that every data point is used even at larger strides? For example, at stride 10, I could include distinct trajectories with frames (0, 10, 20, 30...) and frames (1, 11, 21, 31...).**

You could reformat to use every frame, but it would negligibly improve sampling since the trajectories would be highly correlated. For example, the trajectory from frames (0, 10, 20, 30...) and that from frames (1, 11, 21, 31...) will be very similar. It would be much better to obtain more data.

There are other ways to better use the available data. You could make your grid coarser, so that there is more data per diffusivity grid point. The coarseness of the grid is determined by the grid spacing of the input diffusivity grid). For example, you can make change the spacing in the file "coarse_oneHundred.dat":
> field diffuse cubic -step 2 -f coarse_oneHundred.dat -minVal 0.1
> 

Another way would be to more strongly enforce the smoothness of the result by decreasing "-grad" in a smoothness prior. For example, you could reduce this value from 50 Å/ns to 10 Å/ns (the units of this gradient are units(diffusivity)/units(grid_spacing)):
> prior smooth diffuse -grad 10


&nbsp;

**Why do you set the force to zero, " -force fz -col . " and when we should do this?**

Short answer: the Smoluchowski solvers ignore the biasing force information, so it doesn't matter what it is.

For the simplest method, "ccg" (described in the paper Comer, Chipot, Gonzalez https://doi.org/10.1021/ct300867e ), it is possible to have a different bias force for each event. This allows you to calculate the diffusivity before the ABF method has converged or for an arbitrary time- and position-dependent force (doesn't have to be ABF). The limitation of this method is that it assumes that nothing (neither the biasing force, nor the diffusivity, nor the intrinsic force) changes very much over the lag time. This does not work well for long lag times.
 
The Smoluchowski solvers work for arbitrarily long lag times (although if your lag time is longer than the typical time to diffuse across the entire domain, the probability distribution will be flat and no dynamical information can be obtained). These solvers do not support single biasing force values for each event and ignore the biasing force column (so it doesn't matter if you set it to zero). These solvers implement evolution of the diffusion equation over the whole domain, so biasing force values wouldn't provide enough information. They need a position-dependent biasing force function to be defined over the whole domain. I have written code that allows you to have different biasing force functions (as a function of position) for different parts of the trajectory (TrajSmolCrankBias.H), but this code is currently commented out because I never found it particularly useful. 

For most of my applications, I simply use the Smoluchowski solvers after the ABF method has mostly converged, so that bias force isn't changing much anymore. In this case, the PMF can be assumed to be completely flat (<F_bias + F_intrinsic> = 0) and we don't need to worry about a time-dependent bias or even the intrinsic force (the negative gradient of the potential of mean force).

&nbsp;

**What is the difference between -stride and -skip?**

The "-stride" option strides the original data. So if the time lag between samples in the input data is 1 ps and you use "-stride 10", the events used in the calculation have a time lag of 10 ps. The "-skip" option, on the other hand, is applied after the events have been calculated. So if your data has a time lag of 1 ps, using "-skip" does not change this. It simply keeps only one in every 10 events, which can be done to improve performance. Note that "-skip" and "-stride" can be used at the same time.

Note that, in my nomenclature, and event refers to an initial position, a spatial displacement, and a time lag. The reason for this definition is that the method is based on calculating the probability of the particle being displaced a certain amount after a time lag, given that the particle was at particular position initially. 

&nbsp;

**Why we need to write "mc alpha diffuse" (I believe the example only has "mc diffuse")?**

The "mc" command defines which fields the optimization is applied to. If you do only "mc diffuse", you will only optimize the diffusivity profile, while the alpha profile is kept constant (whatever you set it to initially). 

&nbsp;

**Why we need to define two "groups" with different strides?**

First of all, each Smoluchowski solver can only take data with a single lag time. If you want two different lag times, you have to put your data in two different "groups" and have two different solvers (called trajCost in the input file) that each work on a different group (and lag time).

The fractional Smoluchowski model is supposed to provide a consistent model that doesn't depend on lag time. Despite this, I noticed that I still saw different alpha/diffusivity values when I used a small lag time versus a longer lag time (although these differences were much smaller than differences seen with the non-fractional Smoluchowski model). So having two different groups with different lag times can allow you to obtain a consensus alpha/diffusivity over a range of lag times.

&nbsp;

**With two groups and two lag times, what is the effective lag time of the output diffusivity (if any) ?**

It should be some kind of consensus. Again, the hope of the fractional Smoluchowski model is to describe diffusion with two position-dependent functions D(x) and α(x), with no lag time dependence. However, this fractional Smoluchowski model may still be too simple. 


&nbsp;

**With two groups, how do we control the "coupling" of the alphas and diffuses (if at all)?** 

You should have one diffuse field and one alpha field that get optimized over both lag times. The costs just get added together, so if you want both lag times to have equal weight, they need to have equal numbers of events, which is not very natural (for the same trajectory, a larger lag time will have fewer events than a shorter lag time). The "-weight" option allows you to scale the trajCost to account for this.

&nbsp;

**How to restart a calculation from "preview" files?**

You can just restart the calculation using the already partially optimized fields (diffuse and alpha), instead of the uniform initial field (in my example, I give initial diffusivity values of 100 for all points).

&nbsp;

**What is the significance of the initial values of the "field diffuse", "field alpha" files?**

This is just the initial conditions that the optimization starts from. In an ideal world, the results would converge to the same diffuse(x) and alpha(x) functions regardless of the initial values; however, if they are not very reasonable to start, the optimization may take a very long time to converge or converge to something strange.




