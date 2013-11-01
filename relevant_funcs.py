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
