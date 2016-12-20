DiffusionFusion
===============

This is a program for calculating position-dependent diffusivity functions from trajectory data. It was written by Jeff Comer <jeffcomer at gmail dot com>. The approach is described in:

Comer, Chipot and Gonzalez-Nilo (2013) J. Chem. Theory Comput., 9(2) 876-882, DOI: 10.1021/ct300867e

Installation
------------

You need the Gnu Scientific Library development package installed (for random numbers). Here is the Ubuntu package, which may or may not help you:
$ sudo apt-get install gsl-bin libgsl0-dev

You may also have to install OpenMP support for g++ (multithreading). Here is the Ubuntu package, which may or may not help you:
$ sudo apt-get install libgomp


Then try:
make

On a Mac, this might work better:
make mac

Just running diffusion fusion with no parameters will give you the usage information.

> ./diffusionFusion

See the Examples/ directory for more information.
