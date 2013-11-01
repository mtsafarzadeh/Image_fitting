import pywcs
import numpy as np
import pygoods
import tables
from itertools import *
import pyfits
indices=np.loadtxt('sam_indices_new_try.dat')
header_pacs160=pyfits.getheader('/astro/candels9/user/ferguson/goods_s_rgb/images/gh_goodss_dr1_160_sci.fits')
fir= tables.openFile('/astro/candels9/user/ferguson/mocksv3/FIR.hdf5')
pacs160_sam=fir.root.data.col('pacs160')
irac_ch1_sam=fir.root.data.col('irac_ch1')

wcs_pacs160= pywcs.WCS(header_pacs160)

cat_index=indices[:,6].astype(int)
sam_index=indices[:,0].astype(int)
foundmatch=indices[:,2].astype(int)

cat=pygoods.sextractor('photoz.dat')
sim_herschel_cataloge_2=open('sim_herschel_cataloge_new.dat','w+')
for (real,sim) in izip(cat_index,sam_index):
	if (sim!= -99):
		yo,xo=wcs_pacs160.wcs_sky2pix(cat.ra[real],cat.dec[real],0)
		print>>sim_herschel_cataloge_2, cat.id[real], '\t',sim,'\t',cat.ra[real],'\t',cat.dec[real],'\t',yo[0],'\t',xo[0],'\t', pacs160_sam[sim],'\t', irac_ch1_sam[sim] 
sim_herschel_cataloge_2.close()

