This python code does MCMC exploration of source fluxes on an image ( here PACS, Herschel image)
The input files that are read are :
	1-Simulated Herschel science image ( This is the data in the Bayesian formalism )
	2-The RMS map of the simulated image ( in the Chi2 fitting, acts as the sigma for each pixel ) 
	3-relevant catalogues that the sources are read from, those catalogues provide the RA, Dec, and Redshift (spec_z or photo_z) and fluxes in specific bands
	4-The PSF that is used to construct the simulated image, or in the case of an observation, the PFS that is used for source flux extraction

The flow chart of the code is as follow:
1)Using pywcs , the code converts the RA and Dec to x & y coordinate on the image.
2)Each source with its pixel coordinate and its Flux is put in a matrix of size (NX,NY) where NX & NY are size of the image in each axis

	The construction of the simulated image has these two steps in addition:
		X)The matrix is convolved with the PSF
		X)Using the RMS map, a noise realization is generated and added to the image. The construction of the simulated image is done by this step.
3)sources with a given range of flux ( i.e. brighter than 1 mJy) is chosen in a given part of the image. (these sources are those that we want to guess their flux)
4) python implementation of MCMC affine invariant routine (Emcee) is coded up to explore the plausible flux values. 
5)The pixels that the specified sources fall into are set to zero each time during the MCMC chain.
6)In each step, MCMC guesses fluxes for all the sources asked to find their flux for, changes the value of pixels to the guesses flux
	convolves the image with the PSF, compute a chi2 with respect to the data ( simulated image or science image)
	and returns the ln(probabilty) or -chi2/2.

7) in order to confine the MCMC to search only reasonable flues, there is weak prior imposed on each source flux with a width of 
	3 dex. So that MCMC exploration is within that suggested window for each source.

8) Emcee defines walkers that act as chains in MCMC formalism, and asks for the number of steps that needs to be taken for a simulation.
	Each walker explores the n-dimensional space ( = number of sources to fit for ).
9)The result of the MCMC is saved to a file to be analysed later.
	Either a 2-d data (#sources , # steps) into a text or a data-cube ( #sources, #walkers, #steps) into a fits file.


	
