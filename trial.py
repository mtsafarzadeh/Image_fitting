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
#herschel=pyfits.getdata('/astro/ferguson1/safar/SAMS/simulated_images/simulated_pacs160_with_photo_z_zoom2.fits')
herschel=pyfits.getdata('/astro/ferguson1/safar/SAMS/simulated_images/simulated_pacs160_with_photo_z_zoom2_new_correct.fits')
RMS_pacs160=pyfits.getdata('/astro/ferguson1/safar/SAMS/GH_GOODS-South_PACS160_DR1/gh_goodss_dr1_160_err.fits')
header_pacs160=pyfits.getheader('/astro/candels9/user/ferguson/goods_s_rgb/images/gh_goodss_dr1_160_sci.fits')
fir= tables.openFile('/astro/candels9/user/ferguson/mocksv3/FIR.hdf5')
sam_indices=np.loadtxt('/astro/ferguson1/safar/SAMS/simulated_images/sam_indices_new_try.dat')

kernel_pacs160=pyfits.getdata('/astro/ferguson1/safar/SAMS/simulated_images/psf_herschel.fits')

NX = 701
NY = 601
ZOOM=2
kernel_pacs160_zoom2=interpolation.zoom(kernel_pacs160,ZOOM)
normfactor = sum(interpolation.zoom(kernel_pacs160_zoom2,1./ZOOM))

wcs_pacs160= pywcs.WCS(header_pacs160)
#y_sam_pacs160,x_sam_pacs160=wcs_pacs160.wcs_sky2pix(ra,dec,0)
pacs160_image=np.zeros((NX*ZOOM,NY*ZOOM),dtype=np.float64)
xc_orig=NX/2.
yc_orig=NY/2.
xc_expanded = xc_orig*ZOOM
yc_expanded = yc_orig*ZOOM

cat_index=sam_indices[:,6].astype(int)
sam_index=sam_indices[:,0].astype(int)
photoz_cat=tables.openFile('/astro/ferguson1/safar/SAMS/simulated_images/photoz.hdf5')
specz_cat=tables.openFile('/astro/ferguson1/safar/SAMS/simulated_images/gs_all_tf_h_120919a_multi_specz.hdf5')
ra_cat_photo_z=photoz_cat.root.data.col('ra')
dec_cat_photo_z=photoz_cat.root.data.col('dec')
pacs160_sam=fir.root.data.col('pacs160')


y_cat_photo_z,x_cat_photo_z=wcs_pacs160.wcs_sky2pix(ra_cat_photo_z,dec_cat_photo_z,0)
#index_photo_z=((y_cat_photo_z<343)&(y_cat_photo_z>274)&(x_cat_photo_z<435)&(x_cat_photo_z>367)).nonzero()[0]
#index_photo_z=((y_cat_photo_z<320)&(y_cat_photo_z>290)&(x_cat_photo_z<410)&(x_cat_photo_z>380)).nonzero()[0]
#index_photo_z=((y_cat_photo_z<360)&(y_cat_photo_z>290)&(x_cat_photo_z<410)&(x_cat_photo_z>300)).nonzero()[0]
index_photo_z=((y_cat_photo_z<270)&(y_cat_photo_z>230)&(x_cat_photo_z<270)&(x_cat_photo_z>230)).nonzero()[0]


x_min=x_cat_photo_z[index_photo_z].min().astype(int)-10
x_max=x_cat_photo_z[index_photo_z].max().astype(int)+10
y_min=y_cat_photo_z[index_photo_z].min().astype(int)-10
y_max=y_cat_photo_z[index_photo_z].max().astype(int)+10

#print x_cat_photo_z[index_photo_z]
#print y_cat_photo_z[index_photo_z]
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
#desired_index=(corresponding_fluxes >-3.5).nonzero()[0]
desired_index=(corresponding_fluxes >-3.6).nonzero()[0]
print desired_index
desired=corresponding_fluxes[desired_index]
print desired
ndim=len(desired)
#sys.exit()

#the fluxes of 15 bright sources(>-3.2) from previous MCMC run:
#np.array([-3.06699883, -2.848357  , -2.44119659, -2.43732354, -3.05025237,
#       -1.79345365, -1.88282073, -3.02891038, -3.096002  , -2.05671181,
#       -2.37651056, -3.0015159 , -2.89873345, -3.16458786, -2.4151817 ])

#this is 15 sources fluxes from prvious MCMC run + the new 5 faintest are shifted 1 order of magnitude up
'''
desired=np.array([-3.06699883, -2.848357  ,-3.33368926+1, -2.44119659, -3.49915531+1,
       -2.43732354, -3.05025237,-1.79345365, -1.88282073, -3.02891038, -3.096002  , -2.05671181,
       -2.37651056, -3.0015159 ,-3.20800659+1, -2.89873345, -3.16458786, -3.22111799+1,-3.49912911+1,-2.4151817 ])
'''
faint=np.array([ 2,  4, 14, 17, 18])
bright=np.array([ 0,  1,  3,  5,  6,  7,  8,  9, 10, 11, 12, 13, 15, 16, 19])

#desired=np.array([-2.10314408, -2.01585577, -2.44257923, -2.44882814, -2.07945325,
#      -1.79581421, -1.88503517, -3.23616362, -2.83184716, -2.05116142,
#      -2.37277159, -2.11795198, -2.70044516, -2.17463291, -2.48023795])
#this is based on the result of 10 sources brighter than -3. and those 5 sources fainter are started at an order of magnitude higher.
#bright=np.array([ 2,  3,  5,  6,  7,  8,  9, 10, 12, 14])
#faint=np.array([ 0,  1,  4, 11, 13])

#after one run with faint sources started at wrong position, these are the final medians that are now used to re-run the same MCMC run.
#desired=np.array([-3.16486443, -2.79629887, -2.4238941 , -2.3909444 , -3.15671008,
#       -1.79138201, -1.8753516 , -3.21013146, -3.1073628 , -2.02274163,
#       -2.43814733, -2.9329209 , -2.69014608, -3.10745301, -2.61806779])


#pacs160_image=np.zeros(701*601).reshape(701,601)
cat=sextractor('/astro/ferguson1/safar/SAMS/simulated_images/photoz.dat')
for i in izip(sam_index,cat_index):
        if i[0]!=-99:
        	yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[i[1]],cat.dec[i[1]],0)
#	        yo,xo=wcs_pacs160.wcs_sky2pix(photoz_cat.root.data.col('ra')[i[1]],photoz_cat.root.data.col('dec')[i[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]+=pacs160_sam[i[0]]
#===================================
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
#===================================
#bins=np.logspace(-10,3,49)
chi2=np.zeros(len(index_photo_z))
def digthis(region_coeff):

#       for i in range(0,len(ra_sam)):
#               image[x_sam[i].round(),y_sam[i].round()]=flux_sam[i]+(10**region_coeff[0])# adding background level to the pixels       
        #print region_coeff
#               if (np.any(region_coeff)>3 or np.any(region_coeff<-10)):
        if (np.any(abs(region_coeff-desired)>1.5)):
                return -np.inf
#       if (np.any(abs(region_coeff[faint]-desired[faint])>1)):
#               return -np.inf

#        for qq in index_photo_z:
                #image[x_cat_photo_z[qq].round(),y_cat_photo_z[qq].round()]=(10**region_coeff[pp]) # because the zero's index is used for background.

        for qq in izip(sam_corresponding_index[desired_index],cat_corresponding_index[desired_index]):
                yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[qq[1]],cat.dec[qq[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]=0


        pp=0
        for qq in izip(sam_corresponding_index[desired_index],cat_corresponding_index[desired_index]):
                yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[qq[1]],cat.dec[qq[1]],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
                pacs160_image[xn[0],yn[0]]=pacs160_image[xn[0],yn[0]]+(10**region_coeff[pp])
                pp+=1

#               chi2[pp]=np.log(np.interp(region_coeff[pp+1],np.log10(bins),pdffile[qq,:]/pdffile[qq,:].max()) ) # this is -chi2/2 in fact that we retrun.
#               return len((chi2==-np.inf).nonzero()[0]*(-1e6)

        image_p = hconvolve.padzero2d_i(pacs160_image, r, c)

        fftimage = anfft.fftn(image_p)*fftkernel# * anfft.fftn(kernel_p)
        test_k_temp = anfft.ifftn(fftimage)[:r1,:c1].real
        test_k = interpolation.zoom(test_k_temp,1./ZOOM)/normfactor
        test_chi2=0
        for i in range(x_min,x_max):
                for j in range(y_min,y_max):
                        test_chi2+=(test_k[i,j]-herschel[i,j])**2/(RMS_pacs160[i,j]**2)
        print  -1*(test_chi2/2)
        return -1*(test_chi2/2)#+chi2.sum()

nwalkers = (ndim)*16
# Choose an initial set of positions for the walkers.
#to choose between 10^-8 and 1.
#p0= [np.random.rand(ndim)+[.6,-450] for i in xrange(nwalkers)]
#p0 = [5000*np.random.rand(ndim)+1000 for i in xrange(nwalkers)]

#p0= [np.random.rand(ndim)*-10  for i in xrange(nwalkers)]

#the median is the result of MCMC run on 7 most bright (log(flux) > -2.5)sources. these are used to go fainter. these are found by running MCMC at wrong place for all of them.

#p0=[endofpreviousrun+0.4*np.random.randn(ndim) for i in range(nwalkers)]
#p0= [desired +0.1*np.random.randn(ndim) for i in range(nwalkers)]

#bright=np.array([ 2,  3,  5,  6,  7,  8,  9, 10, 12, 14])
#faint=np.array([ 0,  1,  4, 11, 13])
a=np.ones(ndim)*0.9
#a[bright]=0.2
#p0= [desired+0.5 + (np.random.rand(ndim)*(-2)+1)*a for i in range(nwalkers)]
p0= [desired + (np.random.rand(ndim)*(-2)+1)*a for i in range(nwalkers)]
#p0= [desired +0.2*np.random.randn(ndim) for i in range(nwalkers)]

#p0=[np.array([-5.20417812 ,-7.89563203 ,-3.03564876 ,-4.70245622 ,-3.36751449 ,-7.78029145 ,-1.06330097 ,-5.55249664 ,-3.97559012 ,-4.87208369 ,-8.51483771 ,-6.24832324 ,-8.84646587 ,-4.35006721])for i in xrange(nwalkers)]
#[  205,-482,1119,-2722,70.50230457,1811,3488,37]
#pool = ut.MPIPool()
#if not pool.is_master():
#    pool.wait()
#    sys.exit(0)

# Initialize the sampler with the chosen specs.

sampler = emcee.EnsembleSampler(nwalkers, ndim, digthis,threads=16)
#Run 100 steps as a burn-in.
print time.ctime()
pos, prob, state = sampler.run_mcmc(p0, 10)
print time.ctime()
# Reset the chain to remove the burn-in samples.
sampler.reset()

# Starting from the final position in the burn-in chain, sample for 1000
# steps.
# I changed it to 100
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
