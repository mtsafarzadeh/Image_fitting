from numpy import *
import kdmatch
from pygoods import *

def ab_to_microjy(ab):
    return 10.**((23.9-ab)/2.5)

def microjy_to_ab(mjy):
    return -2.5*log10(mjy)+23.9

def readcat():
    catdir = './'
    imdir = './images/'
    sim = sextractor(catdir+'sim_herschel_cataloge_new.dat')
    hanae = sextractor(catdir+'hanae-catalogue_SEX.cat')
    return hanae, sim


def matchandjoin(cat1,cat2,image="/astro/candels9/user/ferguson/goods_s_rgb/images/gh_goodss_dr1_160_sci.fits",prefix="m",sep=1.0,matchflag="matchflag"):
#image="/astro/candels1/www/data/goodss/mosaics/current/goods_s_all_combined_v1.0/60mas/gs_all_candels_wfc3_f160w_060mas_v1.0_drz.fits"
    """ Match and join two catalogs. The catalogs are matched in RA and DEC.
        They have to have columns named ra and dec. 
        Arguments:
           cat1 -- the reference catalog; data from matching rows from cat2 are added to cat 1.
           cat2 -- the other catalog
           image -- a reference image to help with coordinate transformations (optional)
           prefix -- a prefix to add to each column name that is appended to the original catalog
           sep -- maximum separation in arcsec
           matchflag -- the name to give the flag indicating a match in the new cat1
    """
    #sep,i = kdmatch.matchradec(cat1.ra,cat1.dec,cat2.ra,cat2.dec,fitsfile=image)
    d,i,x0,y0,x1,x2= kdmatch.matchradec(cat1.ra,cat1.dec,cat2.ra,cat2.dec,fitsfile=image)
    '''#near = sep*3600 < sep # 0.2 arcseconds
    near = d < 10 # 0.2 arcseconds
    #i11match = where(sep*3600 < sep)[0]
    i11match = i[:,0]
#    i12match = i[nonzero(near)]
    i12match= i[:,0]
    cat1.addcolumn(matchflag,zeros(len(cat1.dec),dtype=bool),'%1d',"")
    for k in cat2._colnames:
        dt = cat2.__dict__[k].dtype
        newcol = prefix+"_"+k
        if dt == "int32":
            cat1.addcolumn(newcol,zeros(len(cat1.dec),dtype=dt),"%8d")
        else:
            cat1.addcolumn(newcol,zeros(len(cat1.dec),dtype=dt),"%15.6f")
        cat1.__dict__[newcol][i12match] = cat2.__dict__[k][i11match] # Insert the values when a match is found
        cat1.__dict__[matchflag][i12match] = True
    return cat1, i11match, i12match,d,i,x0, y0,x1,x2
'''
    return d,i,x0, y0,x1,x2
