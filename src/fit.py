import numpy as np
import sys
from input import get_imax
from scipy.optimize import curve_fit

def polynomial_for_kimball3(x,a0,a2,a3):
    return a0-x+a2*x**2+a3*x**3

def polynomial_for_kimball5(x,a0,a2,a3,a4,a5):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5

def polynomial_for_kimball7(x,a0,a2,a3,a4,a5,a6,a7):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7

def polynomial_for_kimball9(x,a0,a2,a3,a4,a5,a6,a7,a8,a9):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9

def polynomial_for_kimball11(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11

def polynomial_for_kimball13(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13

def polynomial_for_kimball15(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15

def polynomial_for_kimball17(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15+a16*x**16+a17*x**17

def fit_func(r,x,y,polynomial,imax,args):
    if(args.fitscale>0):
        if(args.opt_method=='lm'):
            sigma=4*np.pi*x[:imax]
            popt,pcov=curve_fit(polynomial,x[:imax],y[:imax],method=args.opt_method,sigma=sigma**-1)
        else:
            sigma=4*np.pi*x[:imax]
            popt,pcov=curve_fit(polynomial,x[:imax],y[:imax],method=args.opt_method,jac='2-point',sigma=sigma**-1)
    else:
        popt,pcov=curve_fit(polynomial,x[:imax],y[:imax],method=args.opt_method)
    return polynomial(r,*popt),popt

def get_fit(r,rclean,g,deg,imax,args):

    kimball_polynomials={
        "3" : polynomial_for_kimball3,
        "5" : polynomial_for_kimball5,
        "7" : polynomial_for_kimball7,
        "9" : polynomial_for_kimball9,
        "11" : polynomial_for_kimball11,
        "13" : polynomial_for_kimball13,
        "15" : polynomial_for_kimball15,
	"17" : polynomial_for_kimball17
    }
    
    # Fits at bondlength
    fit=np.zeros((len(r),))

    fit,popt=fit_func(r,rclean,g,kimball_polynomials[str(deg)],imax,args)

    return fit,popt

def logofg(gs,r,Nx,Npcf,args,crossval):

    arrglog=[]
    r_extrapolated=[]

    for ipcf in range(Npcf):
        glog=np.zeros((Nx,))
        rcopy=r
        i=0
        ind=0
        g=gs[ipcf]
        N_removed=0
        Nxcopy=Nx
        while(i<Nxcopy):
            if(g[i]>0):
                glog[ind]=np.log(g[i])
                ind+=1
                # The following lines are an option for the if-condition section above
                #if(gv[i]>0 and gd[i]>0):# Data is good
                #    glog[ind]=2*np.log(gd[i])-np.log(gv[i])
                #    ind+=1
                
            else: # logarithm cannot be taken        
                if(args.verbosity>0):
                   N_removed+=1
                rcopy=np.delete(rcopy,i); Nxcopy-=1 
                glog=np.delete(glog,ind+1)
            i+=1
        if(args.verbosity>0 and not(crossval)):
            print("PCF {}: {} out of {} points removed because 2*gdmc[i]-gvmc[i] < 0.".format(ipcf,N_removed,Nx))
        arrglog.append(glog)
        r_extrapolated.append(rcopy)

    return arrglog,r_extrapolated
    
def do_fit(r,r_range,p_degs,g,args,crossval=False):

    Nx=len(r)
    Npcf=len(g)
    
    N_degs = len(p_degs)

    # Logarithm of the gs, with negative g-values removed
    glog,rex=logofg(g,r,Nx,Npcf,args,crossval)
    
    # The array to hold the fits
    logfits=np.zeros((Nx,Npcf,N_degs))
    
    # Optimal polynomial coefficients
    cpol=[]

    # Collect fits and fitting polynomial coefficients
    for i in range(Npcf):
        imax=get_imax(rex[i],r_range)
        cpol_deg=[]
        for d in range(N_degs):
            logfits[:,i,d],popt=get_fit(r,rex[i],glog[i],p_degs[d],imax,args)
            cpol_deg.append(popt)
            
        cpol.append(cpol_deg)

    try:
        fits=np.exp(logfits[:get_imax(r,args.explim),:,:])
    except:
        print("Error in exponentiation,")
        print(np.mean(logfits,axis=1))
        sys.exit("Try to limit the length scale of the exponentiated values.")
    
    return fits,logfits,glog,rex,cpol

