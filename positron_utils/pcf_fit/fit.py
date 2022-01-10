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

def polynomial_for_kimball19(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a19):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15+a16*x**16+a17*x**17+a19*x**19

def polynomial_for_kimball21(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a19,a21):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15+a16*x**16+a17*x**17+a19*x**19+a21*x**21

def polynomial_for_kimball23(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a19,a21,a23):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15+a16*x**16+a17*x**17+a19*x**19+a21*x**21+a23*x**23

def polynomial_for_kimball25(x,a0,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a19,a21,a23,a25):
    return a0-x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8+a9*x**9+a10*x**10+a11*x**11+a12*x**12+a13*x**13+a14*x**14+a15*x**15+a16*x**16+a17*x**17+a19*x**19+a21*x**21+a23*x**23+a25*x**25

def fit_func3(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball3,x[:imax],y[:imax],method=mthd)
    return polynomial_for_kimball3(r,*popt),popt

def fit_func5(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball5,x[:imax],y[:imax],method=mthd)
    return polynomial_for_kimball5(r,*popt),popt
        
def fit_func7(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball7,x[:imax],y[:imax],method=mthd)
    return polynomial_for_kimball7(r,*popt),popt

def fit_func9(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball9,x[:imax],y[:imax],method=mthd)
    return polynomial_for_kimball9(r,*popt),popt

def fit_func11(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball11,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball11(r,*popt),popt

def fit_func13(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball13,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball13(r,*popt),popt

def fit_func15(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball15,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball15(r,*popt),popt

def fit_func17(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball17,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball17(r,*popt),popt

def fit_func19(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball19,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball19(r,*popt),popt

def fit_func21(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball21,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball21(r,*popt),popt

def fit_func23(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball23,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball23(r,*popt),popt

def fit_func25(r,x,y,imax,mthd):
    popt,pcov=curve_fit(polynomial_for_kimball25,x[:imax],y[:imax],method=mthd,maxfev=20000)
    return polynomial_for_kimball25(r,*popt),popt

def get_fit(r,rclean,g,deg,imax,mthd):
    # Number of degrees to be tested
    
    # Fits at bondlength
    fit=np.zeros((len(r),))

    if(deg==3):
        fit,popt=fit_func3(r,rclean,g,imax,mthd)
    elif(deg==5):
        fit,popt=fit_func5(r,rclean,g,imax,mthd)
    elif(deg==7):
        fit,popt=fit_func7(r,rclean,g,imax,mthd)
    elif(deg==9):
        fit,popt=fit_func9(r,rclean,g,imax,mthd)
    elif(deg==11):
        fit,popt=fit_func11(r,rclean,g,imax,mthd)
    elif(deg==13):
        fit,popt=fit_func13(r,rclean,g,imax,mthd)
    elif(deg==15):
        fit,popt=fit_func15(r,rclean,g,imax,mthd)
    elif(deg==17):
        fit,popt=fit_func17(r,rclean,g,imax,mthd)
    elif(deg==19):
        fit,popt=fit_func19(r,rclean,g,imax,mthd)
    elif(deg==21):
        fit,popt=fit_func21(r,rclean,g,imax,mthd)
    elif(deg==23):
        fit,popt=fit_func23(r,rclean,g,imax,mthd)
    elif(deg==25):
        fit,popt=fit_func25(r,rclean,g,imax,mthd)
    else:
        sys.exit("This degree of polynomial, "+str(deg)+', is not implemented.')
        
    return fit,popt

def logofg(gs,r,Nx,Npcf,args):

    arrglog=[]
    r_extrapolated=[]

    for ipcf in range(Npcf):
        glog=np.zeros((Nx,))
        rcopy=r
        i=0
        ind=0
        g=gs[ipcf]
        while(i<Nx):
            if(g[i]>0):
                glog[ind]=np.log(g[i])
                ind+=1
                # The following lines are an option for the if-condition section above
                #if(gv[i]>0 and gd[i]>0):# Data is good
                #    glog[ind]=2*np.log(gd[i])-np.log(gv[i])
                #    ind+=1
                
            else: # logarithm cannot be taken        
                if(args.verbosity>0):
                    print("ind: "+str(ipcf)+" i: "+str(i)+", 2gd-gv: "+str(g[i]))
                rcopy=np.delete(rcopy,i)
                glog=np.delete(glog,ind+1)
            i+=1
        
        arrglog.append(glog)
        r_extrapolated.append(rcopy)

    return arrglog,r_extrapolated
    
def do_fit(r,r_range,g,args):

    Nx=len(r)
    Npcf=len(g)
    
    # List of polynomial orders to be tested
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)

    glog,rex=logofg(g,r,Nx,Npcf,args)
    
    # The array to hold the fits
    logfits=np.zeros((Nx,Npcf,N_degs))
    
    # Optimal polynomial coefficients
    cpol=[]

    # Collect fits and fitting polynomial coefficients
    for i in range(Npcf):
        imax=get_imax(rex[i],r_range)
        cpol_deg=[]
        for d in range(N_degs):
            logfits[:,i,d],popt=get_fit(r,rex[i],glog[i],p_degs[d],imax,args.opt_method)
            cpol_deg.append(popt)
        cpol.append(cpol_deg)

    try:
        fits=np.exp(logfits[:get_imax(r,args.explim),:,:])
    except:
        print("SUCKABLIEVT")
        print(np.mean(logfits,axis=1))
        sys.exit("Try to limit the scope of exponentiation.")
        
    return fits,logfits,glog,rex,cpol

