
hanae, sim= readcat()
d,i,x0,y0,x1,y1=matchandjoin(sim,hanae,prefix='sim',matchflag='match')

predicted=np.zeros(len(d[:,0]))
for p in range(0,len(d[:,0])-1 ):
    ind=(d[p,:]<2).nonzero()[0]
#    ind=np.array([0]) # closest match 	
    indd=argsort(sim.irac_ch1[i[p,:][ind]])[-1]
    predicted[p]=sim.pacs160[i[p,:][indd]]
'''
plt.scatter(predicted,hanae.flux,1)
plt.ylim(0.0000001,.1)
plt.xlim(0.0000001,.1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('the closest match only by distance')
plt.ylabel('Hanae-catalogue')
'''
faint=(hanae.flux<0.005).nonzero()[0]
#plt.hist(predicted[faint]-hanae.flux[faint],100,histtype='step',normed=True,color='blue')
#plt.hist( (predicted[faint]-hanae.flux[faint])/hanae.fluxerr[faint],100,histtype='step',normed=True,color='blue')
bright=(hanae.flux>0.005).nonzero()[0]
#plt.hist(predicted[bright]-hanae.flux[bright],100,histtype='step',normed=True,color='green')
#plt.hist( (predicted[bright]-hanae.flux[bright])/hanae.fluxerr[bright],100,histtype='step',normed=True,color='green')

sigmas=np.loadtxt('sigmas_predictions.dat')
s16=sigmas[:,0]
s84=sigmas[:,1]
s50=sigmas[:,6]
s_true=sigmas[:,7]
plt.hist((10**s50-10**s_true)/((10**s84-10**s16)*2.),histtype='step',color='red',normed=True)

plt.xlabel('(MCMC_derived_flux-True_flux, Jy)/error')
#plt.title('Blue--faint(<0.005 Jy), green--bright(>0.005 Jy), red--our method')
plt.title(' red--our method')
