import numpy as np
import sys
import warnings
from input import get_imax
from output import *

def fit_average_errors(fit,g,N_degs,r):
    '''
    Takes in the average of all fits for each polynomial degree
    and the average of the extrapolated PCF:s, and computes the 
    mean absolute error and mean square error of fit against
    PCF for each of the polynomial degrees.
    '''    

    imax=get_imax(r,1.)

    fit_errors=np.mean(np.absolute((fit[:imax,:].T-g[:imax]).T),axis=0)
    fit_sqerrors=np.mean(((fit[:imax,:].T-g[:imax]).T)**2,axis=0)
        
    return fit_errors,fit_sqerrors

def cross_validation_error(fitset,gex,r,args):

    N_degs=fitset.shape[2]
    Npcf=fitset.shape[1]
    Nx=fitset.shape[0]

    gex=np.array(gex).T
    
    if(args.valset1[0]>-.1):
        fits=np.zeros((Nx,2,N_degs))
        fits[:,0,:]=np.mean(fitset[:,args.valset1,:],axis=1)
        fits[:,1,:]=np.mean(fitset[:,args.valset2,:],axis=1)
        g=np.zeros((Nx,2))
        g[:,0]=np.mean(gex[:,args.valset1],axis=1)
        g[:,1]=np.mean(gex[:,args.valset2],axis=1)

    else:
        fits=fitset
        g=gex

    imax=get_imax(r,1.)
    imax0=get_imax(r,-1)
    '''
    for i in range(imax0,imax):
        for j in range(g.shape[1]):
            x=g[i,j]
            y=np.mean(fits[i,:,0])
            #print("----------------")
            #print(x,y)
            #print(np.absolute(x-y))
            if(i>10 and False):
                sys.exit()
            if(np.absolute(x-y)>0.2*y):
                x1=max(0,i-5)
                x2=min(Nx,i+5)
                g[i,j]=np.mean(g[x1:x2,j])
                #print(np.absolute(g[i,j]-y))
    '''         

    
    cverror=[]
    for i in range(N_degs):
        mse=[]
        for itraining in range(fits.shape[1]):
            fit=np.exp(fits[imax0:imax,itraining,i])
            for ivalidation in range(fits.shape[1]):
                if(ivalidation==itraining):
                    continue
                error=np.absolute(fit-g[imax0:imax,ivalidation])**1
                for x in range(error.shape[0]):
                    error[x]*=4.*np.pi*r[i]**2
                              
                mse.append(np.mean(error))
        mse=sum(mse)/len(mse)
        cverror.append(mse)

    cverror2=[]
    for i in range(N_degs):
        mse=[]
        for itraining in range(fits.shape[1]):
            fit=np.exp(fits[imax0:imax,itraining,i])
            for ivalidation in range(fits.shape[1]):
                if(ivalidation==itraining):
                    continue
                error=np.absolute(fit-g[imax0:imax,ivalidation])**2
                for x in range(error.shape[0]):
                    error[x]*=4.*np.pi*r[i]**2
                    mse.append(np.mean(error))
                #mse.append(np.mean(np.absolute(fit-g[imax0:imax,ivalidation])**2))
        mse=sum(mse)/len(mse)
        cverror2.append(mse)
        
    return cverror,cverror2
        
def fit_statistics(args,fits,g,r,r_range):

    # The dimensions of our current experiment
    Nx=len(r)
    Npcf=len(g)
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)
    bond_distance=args.fit_range*args.lat_vec
    
    # Averages of the PCFs and fits
    g_average=np.mean(g,axis=0)
    try:
        fit_average=np.exp(np.mean(fits[:get_imax(r,args.explim),:,:],axis=1))
    except:
        print("SUCKABLIEVT")
        print(np.mean(fits,axis=1))
        
    fit_errors, fit_sqerrors = fit_average_errors(fit_average,g_average,N_degs,r)

    cverror,cverror2=cross_validation_error(fits,g,r,args)
    
    imax=get_imax(r,r_range)
    if args.plot == 1:
        plot_results(p_degs,N_degs,g_average,fit_average,r,imax)

    return fit_errors,fit_sqerrors,cverror,cverror2

def mean_and_error(g,ws):
    '''
    Input: (N_pcf,N_degs)-order matrix of PCF-zeroes.
    Output: 
      * Means over separate PCF's.
      * Mean average errors of PCF's.
      * Standard deviations of PCF's.

    
    '''
    gt=[]
    for i in range(g.shape[0]):
        for j in range(int(ws[i])):
            gt.append(g[i])

    gtemp=np.array(gt)
    
    gmean=np.mean(gtemp,axis=0)
    mean_error=np.mean(np.absolute(gmean-gtemp),axis=0)
    standard_deviation=np.std(gtemp,axis=0)/np.sqrt(g.shape[0])
        
    return gmean,mean_error,standard_deviation
