import tables
from pygoods import *
import numpy as np
import multiprocessing
import threading
def readhdf5():
    fir= tables.openFile('/astro/candels9/user/ferguson/mocksv3/FIR.hdf5')
    d=   tables.openFile('/astro/candels9/user/ferguson/mocksv3/galphotdust.hdf5')
    p=   tables.openFile('/astro/candels9/user/ferguson/mocksv3/galprop.hdf5')
    h=   tables.openFile('/astro/candels9/user/ferguson/mocksv3/halos.hdf5')  

    return fir,d,p,h

def find_a_model(z0,z1,h0,h1):
    """ Find a SAM galaxy within the range of redshift and h band magnitude of the hubble image """
    fir,d,p,h = readhdf5()
    z_sam= fir.root.data.col('z')
    h_sam=d.root.data.col('wfc3f160w')
    bc = tables.Expr('(z_sam>%.3f) & (z_sam<%.3f) & (h_sam>%.3f) & (h_sam<%.3f)' % (z0,z1,h0,h1)).eval()
    igal = np.nonzero(bc)
    try:
    	rr=np.random.randint(0,len(igal[0]))
    except ValueError:
	return -99,len(igal[0]) 
    return igal[0][rr], len(igal[0])

def get_variables():
	photoz_cat=tables.openFile('photoz.hdf5')
	specz_cat=tables.openFile('gs_all_tf_h_120919a_multi_specz.hdf5')
	photo_z=photoz_cat.root.data.col('photo_z')
	spec_z=photoz_cat.root.data.col('spec_z')
	h_band_hst_jy=specz_cat.root.data.col('wfc3_f160w_flux')
	h_band_hst_ab=23.9-2.5*np.log10(h_band_hst_jy)

	return photo_z,spec_z,h_band_hst_ab
def run_the_program(q):
	file=open('/astro/ferguson1/safar/SAMS/simulated_images/sam_indices_new_try.dat','a')
	chunck=q.get(block=False)
	for g in xrange(chunck,chunck+len(photo_z)/10):
		z0,z1=photo_z[g]-0.05,photo_z[g]+0.05
		h0,h1=h_band_hst_ab[g]-.1,h_band_hst_ab[g]+.1
		result=find_a_model(z0,z1,h0,h1)	
		print>>file, result[0],'\t',result[1],'\t',z0,'\t',z1,'\t',"%.5g" %h0,'\t',"%.5g" %h1,'\t',g
		print result[0],'\t',result[1],'\t',z0,'\t',z1,'\t',"%.5g" %h0,'\t',"%.5g" %h1,'\t',g
	return
if __name__=='__main__':
	nworkers=10
	photo_z,spec_z,h_band_hst_ab=get_variables()
	work_queue=multiprocessing.Queue()
	for i in xrange(0,len(photo_z),len(photo_z)/nworkers):
		work_queue.put(i)

	file=open('/astro/ferguson1/safar/SAMS/simulated_images/sam_indices_new_try.dat','w')
	processes=[multiprocessing.Process(target=run_the_program,args=(work_queue,)) for i in range(nworkers)]
	for p in processes:
		p.start()
	for p in processes:
		p.join()
	file.close()
'''

if __name__=='__main__':
	jobs=[]
	for i in range(10):
		p=threading.Thread(target=run_the_program)
		jobs.append(p)
		p.start()
'''
