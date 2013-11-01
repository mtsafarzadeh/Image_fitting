import numpy as np
import pywcs
from pygoods import *
import pyfits
import hconvolve
import anfft
import pygoods
from itertools import *
from scipy.ndimage import interpolation
import tables
#cat=sextractor('/astro/ferguson1/safar/SAMS/catalog/gs_all_tf_h_120919a_multi_specz.cat')
NX = 701
NY = 601
ZOOM=2


cat=pygoods.sextractor('photoz.dat')
fir= tables.openFile('/astro/candels9/user/ferguson/mocksv3/FIR.hdf5')
pacs160_sam=fir.root.data.col('pacs160')
#pacs160_sam=np.loadtxt('/astro/ferguson1/safar/SAMS/newrun/FIR_pacs160.dat')
header_pacs160=pyfits.getheader('/astro/candels9/user/ferguson/goods_s_rgb/images/gh_goodss_dr1_160_sci.fits')
kernel_pacs160=pyfits.getdata('../test_mcmc/psf_herschel.fits')
RMS_pacs160=pyfits.getdata('/astro/ferguson1/safar/SAMS/GH_GOODS-South_PACS160_DR1/gh_goodss_dr1_160_err.fits')

kernel_pacs160_zoom2=interpolation.zoom(kernel_pacs160,ZOOM)

normfactor = sum(interpolation.zoom(kernel_pacs160_zoom2,1./ZOOM))
#header_irac_ch2=pyfits.getheader('/data/raid2/spitzer_archive/redux/irac/s12/v0p30/s12_ch2_v0.30_sci.fits')
#RMS_irac_ch2=pyfits.getdata('/data/raid2/spitzer_archive/redux/irac/s12/v0p30/s12_ch1_v0.30_rms.fits')

noise_pacs160=np.random.normal(0,RMS_pacs160)

wcs_pacs160= pywcs.WCS(header_pacs160)
#y_sam_pacs160,x_sam_pacs160=wcs_pacs160.wcs_sky2pix(ra,dec,0)
pacs160_image=np.zeros((NX*ZOOM,NY*ZOOM),dtype=np.float64)
xc_orig=NX/2.
yc_orig=NY/2.
xc_expanded = xc_orig*ZOOM
yc_expanded = yc_orig*ZOOM

#indices=np.loadtxt('Sam_indices_used_new.dat')
indices=np.loadtxt('sam_indices_new_try.dat')
cat_index=indices[:,6].astype(int)
sam_index=indices[:,0].astype(int)
foundmatch=indices[:,2].astype(int)
p=0
'''
for i,j,k in cat_index,sam_index,foundmatch:
        if foundmatch==0:
#        irac_ch1_image[x_sam_irac_ch1[i].round(),y_sam_irac_ch1[i].round()]=irac_ch1_sam[j]
                irac_ch1_image[wcs_irac_ch1.wcs_sky2pix(cat.ra[cat_index],cat.dec[cat_index],0)]=irac_ch1_sam[j]
'''
'''
for i in cat_index:
        if foundmatch[p]==0:
                y,x=wcs_pacs160.wcs_sky2pix(cat.ra[i],cat.dec[i],0)
                pacs160_image[x[0],y[0]]=pacs160_sam[sam_index[p]]
                gtk.append(pacs160_sam[sam_index[p]])
        p+=1
'''
for i in izip(sam_index,cat_index):
        if i[0]!=-99:
                yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[i[1]-1],cat.dec[i[1]-1],0)
                xn,yn = ZOOM*(xo-xc_orig)+xc_expanded,ZOOM*(yo-yc_orig)+yc_expanded
		pacs160_image[xn[0],yn[0]]=pacs160_sam[i[0]]

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
image_p = hconvolve.padzero2d_i(pacs160_image, r, c)
fftimage = anfft.fftn(image_p)*fftkernel# * anfft.fftn(kernel_p)
image_conv = anfft.ifftn(fftimage)[:r1,:c1].real

image_blk = interpolation.zoom(image_conv,1./ZOOM)/normfactor

final=image_blk+noise_pacs160

pyfits.writeto('/astro/ferguson1/safar/SAMS/simulated_images/simulated_pacs160_with_photo_z_zoom2_new_correct_ID.fits',final,header_pacs160)

