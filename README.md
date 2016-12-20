DiffusionFusion
===============

This is a program for calculating position-dependent diffusivity functions from trajectory data. It was written by Jeff Comer <jeffcomer at gmail dot com>. The approach is described in:

Comer, Chipot and Gonzalez-Nilo (2013) J. Chem. Theory Comput., 9(2) 876-882, DOI: 10.1021/ct300867e


See the following for information on using reflecting boundary conditions, which are useful for an orientation angle:

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
