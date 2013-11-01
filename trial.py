import acor
import emcee
import hconvolve
from multiprocessing import *
import anfft
import emcee.utils as ut
import tables
import time
import pywcs
import sys
import pyfits 
from scipy.ndimage import interpolation
from itertools import *
from pygoods import *
import numpy as np

### reading the input data ( simulated image, PSF, RMS map, Catalogues etc.)
herschel=pyfits.getdata('/astro/ferguson1/safar/SAMS/simulated_images/simulated_pacs160_with_photo_z_zoom2_new_correct.fits')
RMS_pacs160=pyfits.getdata('/astro/ferguson1/safar/SAMS/GH_GOODS-South_PACS160_DR1/gh_goodss_dr1_160_err.fits')
header_pacs160=pyfits.getheader('/astro/candels9/user/ferguson/goods_s_rgb/images/gh_goodss_dr1_160_sci.fits')
kernel_pacs160=pyfits.getdata('/astro/ferguson1/safar/SAMS/simulated_images/psf_herschel.fits')
fir= tables.openFile('/astro/candels9/user/ferguson/mocksv3/FIR.hdf5')
sam_indices=np.loadtxt('/astro/ferguson1/safar/SAMS/simulated_images/sam_indices_new_try.dat')
cat_index=sam_indices[:,6].astype(int)
sam_index=sam_indices[:,0].astype(int)

cat=sextractor('/astro/ferguson1/safar/SAMS/simulated_images/photoz.dat')
photoz_cat=tables.openFile('/astro/ferguson1/safar/SAMS/simulated_images/photoz.hdf5')
specz_cat=tables.openFile('/astro/ferguson1/safar/SAMS/simulated_images/gs_all_tf_h_120919a_multi_specz.hdf5')
ra_cat_photo_z=photoz_cat.root.data.col('ra')
dec_cat_photo_z=photoz_cat.root.data.col('dec')
pacs160_sam=fir.root.data.col('pacs160')

###Reading the input data is done

NX = 701
NY = 601
ZOOM=2
kernel_pacs160_zoom2=interpolation.zoom(kernel_pacs160,ZOOM)
normfactor = sum(interpolation.zoom(kernel_pacs160_zoom2,1./ZOOM))

wcs_pacs160= pywcs.WCS(header_pacs160)
pacs160_image=np.zeros((NX*ZOOM,NY*ZOOM),dtype=np.float64)
xc_orig=NX/2.
yc_orig=NY/2.
xc_expanded = xc_orig*ZOOM
yc_expanded = yc_orig*ZOOM

y_cat_photo_z,x_cat_photo_z=wcs_pacs160.wcs_sky2pix(ra_cat_photo_z,dec_cat_photo_z,0)


####################selecting sources to fit for by their pixel coordinate###############################
#index_photo_z=((y_cat_photo_z<343)&(y_cat_photo_z>274)&(x_cat_photo_z<435)&(x_cat_photo_z>367)).nonzero()[0]
#index_photo_z=((y_cat_photo_z<320)&(y_cat_photo_z>290)&(x_cat_photo_z<410)&(x_cat_photo_z>380)).nonzero()[0]
#index_photo_z=((y_cat_photo_z<360)&(y_cat_photo_z>290)&(x_cat_photo_z<410)&(x_cat_photo_z>300)).nonzero()[0]
index_photo_z=((y_cat_photo_z<270)&(y_cat_photo_z>230)&(x_cat_photo_z<270)&(x_cat_photo_z>230)).nonzero()[0]

#defining the region to compute chi2 based on where the sources are on the image ± 10 pixels on each side.
x_min=x_cat_photo_z[index_photo_z].min().astype(int)-10
x_max=x_cat_photo_z[index_photo_z].max().astype(int)+10
y_min=y_cat_photo_z[index_photo_z].min().astype(int)-10
y_max=y_cat_photo_z[index_photo_z].max().astype(int)+10

#print x_cat_photo_z[index_photo_z]
#print y_cat_photo_z[index_photo_z]

#finding the corresponding galaxy in SAM MOCK catalogue that is used to construct to construct the simulated image.
corresponding_index=[]
for i in range(0,len(index_photo_z)):
	try:
		(cat_index==index_photo_z[i]).nonzero()[0][0]
		corresponding_index.append((cat_index==index_photo_z[i]).nonzero()[0][0])
	except IndexError:
		continue

cat_corresponding_index=cat_index[corresponding_index]
sam_corresponding_index=sam_index[corresponding_index]
#print "index_photo_z",index_photo_z
#print "cat_corresponding_index",cat_corresponding_index
#print "sam_corresponding_index",sam_corresponding_index

corresponding_fluxes=np.log10(pacs160_sam[sam_corresponding_index])
#print corresponding_fluxes

#select sources brighter than a given limit ( in log10(Jy))
desired_index=(corresponding_fluxes >-3.6).nonzero()[0]
print desired_index
desired=corresponding_fluxes[desired_index]
print desired
ndim=len(desired)
#sys.exit()

#populate the image with the sources that are used to construct the simulated image, later on the pixels that contain sources
#we are fitting for will be set to zero in MCMC chain but other sources are left at their true values.
for i in izip(sam_index,cat_index):
        if i[0]!=-99:
        	yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[i[1]],cat.dec[i[1]],0)
#	        yo,xo=wcs_pacs160.wcs_sky2pix(photoz_cat.root.data.col('ra')[i[1]],photoz_cat.root.data.col('dec')[i[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]+=pacs160_sam[i[0]]
                

#=================constructing the PSF Kernel==================
r1, c1 = pacs160_image.shape
r2, c2 = kernel_pacs160_zoom2.shape
if r2 % 2 == 0:
	r = r1 + r2/2
else:
	r = r1 + (r2 + 1) / 2
if c2 % 2 == 0:
	c = c1 + c2/2 + 1
else:
	c = c1 + (c2 + 1) / 2

kernel_p = hconvolve.padzero2d_k(kernel_pacs160_zoom2, r, c)
fftkernel=anfft.fftn(kernel_p)
#===========This will be used in the MCMC routine, it is calculated here to avoid doing it each time in the loop===========

chi2=np.zeros(len(index_photo_z))

#explorefluxes is the function that returns the ln(probability)
def explorefluxes(region_coeff):
        #weak prior of 3 dex wide on each source is imposed to limit the exploration of walkers with that plausible range.
        if (np.any(abs(region_coeff-desired)>1.5)):
                return -np.inf
        #at each step the pixels that contain the sources of interest is set to zero.
        for qq in izip(sam_corresponding_index[desired_index],cat_corresponding_index[desired_index]):
                yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[qq[1]],cat.dec[qq[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]=0

        # Pixel values are set to the sources fluxes guesses by the walkers at each step
        pp=0
        for qq in izip(sam_corresponding_index[desired_index],cat_corresponding_index[desired_index]):
                yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[qq[1]],cat.dec[qq[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]+=(10**region_coeff[pp])
                pp+=1
                
        #convolution is performed after sources are put in their relevent pixel
        image_p = hconvolve.padzero2d_i(pacs160_image, r, c) # taking the FFT of the image
        fftimage = anfft.fftn(image_p)*fftkernel#  mutiply by FFT of Kernel
        test_k_temp = anfft.ifftn(fftimage)[:r1,:c1].real # Inverse FFT the result and take the real part
        test_k = interpolation.zoom(test_k_temp,1./ZOOM)/normfactor  # zoom out by the factor PSF was zoomed in
        
        #compute the chi2 at each step
        test_chi2=0
        for i in range(x_min,x_max):
                for j in range(y_min,y_max):
                        test_chi2+=(test_k[i,j]-herschel[i,j])**2/(RMS_pacs160[i,j]**2)
        print  -1*(test_chi2/2)
        return -1*(test_chi2/2)

nwalkers = (ndim)*16
# Choose an initial set of positions for the walkers.

a=np.ones(ndim)*0.9
#p0= [desired+0.5 + (np.random.rand(ndim)*(-2)+1)*a for i in range(nwalkers)]
p0= [desired + (np.random.rand(ndim)*(-2)+1)*a for i in range(nwalkers)]

# Initialize the sampler with the chosen specs.

sampler = emcee.EnsembleSampler(nwalkers, ndim, explorefluxes,threads=16)
#Run 10 steps as a burn-in.
print time.ctime()
pos, prob, state = sampler.run_mcmc(p0, 10)
print time.ctime()
# Reset the chain to remove the burn-in samples.
sampler.reset()

# Starting from the final position in the burn-in chain, sample for 50 steps.
t1=time.ctime()
sampler.run_mcmc(pos, 50, rstate0=state)
t2=time.ctime()
print t1
print t2
print 'finished sampling'
print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))

#try:
#    import matplotlib.pyplot as pl
#except ImportError:
#    print("Try installing matplotlib to generate some sweet plots...")
#else:
#plt.hist(sampler.flatchain[:,0], 100)
#    pl.show()
#plt.savefig('hist.png')

#np.savetxt('flat_sim_herschel_witherror_psfcorrect.txt',sampler.flatchain)
np.savetxt('flat_coeff_3.6_tophat_np_inf_1.5_sources_startingattrueflux_patch1.txt',sampler.flatchain)
pyfits.writeto('data-cube-patch1.fits',sampler.chain)

def padzero2d_i(a, nr, nc):
   """
   a: array to be padded with zero
   nr: number of rows of the output
   nc: number of columns of the output
   pad zeros to the end of image
   """
   s = list(a.shape)

   index = [slice(None)] * len(s)
   index[0] = slice(0, s[0])
   index[1] = slice(0, s[1])
   s[0] = nr
   s[1] = nc
   z = np.zeros(s, a.dtype.char)
   z[index] = a
   return z
