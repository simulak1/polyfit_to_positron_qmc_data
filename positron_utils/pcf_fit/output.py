import numpy as np
import matplotlib.pyplot as plt

def plot_crossval_results(MSE,pd,rf):

    Nd=pd.shape[0]
    Nr=rf.shape[0]

    for i in range(Nd):
        plt.plot(rf,MSE[:,i],'-*',label='Polynomial degree {}'.format(pd[i]))
    plt.xlabel('Fitting range')
    plt.legend()
    

def plot_results(args,fits,logfits,gex,glog,r,rex,popt):

    from input import get_imax
    
    def _plot_main(fignum,r,imax,g,f,Nd,p_degs,Npcf):
        plt.figure(1)
        plt.plot(r,g,'k-',linewidth=2,label='Extrapolated data')
        for i in range(N_degs):
            plt.plot(r[:imax],f[:imax,i],label=str(p_degs[i])+'-order pol.')
        plt.grid()
        plt.xlabel('Positron-electron distance (Bohr)')
        plt.ylabel('Pair correlation function')
        plt.title('Graphs of {} averaged fits and pcfs.'.format(Npcf))
        plt.legend()

    def _plot_twist(Nx,Ny,Npcf,Nd,rex,r,imax,glog,logfits):
        fig,ax=plt.subplots(Nx,Ny)
        ind=0
        for i in range(Nx):
            for j in range(Ny):
                if(ind<Npcf):
                    ax[i,j].plot(rex[ind],glog[ind],'k-')
                    for k in range(N_degs):
                        ax[i,j].plot(r[:imax],logfits[:imax,ind,k],label="Pol. deg. {}".format(p_degs[k]))
                    ind+=1
    
    def _plot_error(Nx,Ny,Npcf,Nd,rex,glog,popt):
        fig,ax=plt.subplots(Nx,Ny)
        fig2,ax2=plt.subplots(Nx,Ny)
        ind=0
        for i in range(Nx):
            for j in range(Ny):
                if(ind<Npcf):
                    for k in range(N_degs):
                        coeff=np.insert(popt[ind][k],1,-1)
                        fit=np.polyval(coeff[::-1],rex[ind])
                        error=fit-glog[ind]
                        #error2=error*4*np.pi*rex[ind]**1
                        error2=error*rex[ind]**1
                        #ax[i,j].plot(rex[ind],error,label="Fit {}, deg {}".format(ind,k))
                        ax[i,j].plot(rex[ind][:get_imax(rex[ind],2)],error[:get_imax(rex[ind],2)],label="Fit {}, deg {}".format(ind,k))
                        ax[i,j].legend()
                        ax2[i,j].plot(rex[ind][:get_imax(rex[ind],2)],error2[:get_imax(rex[ind],2)],label="Fit {}, deg {}".format(ind,k))
                        ax2[i,j].legend()
                    ind+=1
                    
    
    r_range=args.fit_range*args.lat_vec

    imax=get_imax(r,r_range)
    
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)

    Npcf=fits.shape[1]

    g_average=np.mean(gex,axis=0)

    fit_average=np.mean(fits,axis=1)

    print(args.plot)
    if(args.plot>0):
        _plot_main(1,r,imax,g_average,fit_average,N_degs,p_degs,Npcf)
    if(args.plot>1):
        _plot_twist(3,3,Npcf,N_degs,rex,r,imax,glog,logfits)
    if(args.plot>2):
        _plot_error(3,3,Npcf,N_degs,rex,glog,popt)
    
    return


def make_table(args,m,mt,fe,fsqe,e,std,stdt,gzeros,lifetimes):
    corepart=args.corepart

    print(" ")
    print("latex table:")
    print(" ")
    print("\\begin{table*}[!t]")
    print("\caption{ $g(0)$ and $\tau (ps)$, computed with PCF histogram reblocked by *** and fitted to range "+str(args.fit_range*args.lat_vec)+" au.}")
    string="X"*(m.shape[0]+1)
    print("\\begin{tabularx}{\\textwidth}{"+string+"}")
    print("\hline\hline")
    str1=" "
    str2="Fit error "
    str3="Fit squared error "
    str4="g(0) from PCF "
    str5="STD "
    str6="Lifetime "
    str7="STD "
    for i in range(m.shape[0]):
        deg=args.min_pol+i*2
        str1+="& {}-deg pol. ".format(deg)
        str2+="& {0:.4f} ".format(fe[i])
        str3+="& {0:.4f} ".format(fsqe[i])
        str4+="& {0:.4f} ".format(m[i])
        str5+="& {0:.4f} ".format(std[i])
        str6+="& {0:.4f} ".format(mt[i])
        str7+="& {0:.4f} ".format(stdt[i])

    print(str1+"\\\\")
    print("\hline")
    print(str4+"\\\\")
    print(str5+"\\\\")
    print("\hline")
    print(str6+"\\\\")
    print(str7+"\\\\")
    print("\hline")
    print(str2+"\\\\")
    print(str3+"\\\\")
    print("\hline\hline")
    print("\end{tabularx}")
    print("\end{table*}")
    print(" ")

def print_output(args,m,mt,fe,fsqe,cve,cve2,e,std,stdt,gzeros,lifetimes):
    
    corepart=args.corepart
    for i in range(m.shape[0]):
        deg=args.min_pol+i*2
        print("="*40)
        print("Polynomial degree: {}".format(deg))
        print("Mean squared error: {:.5f}, cross-validation MSE: {:.5f}".format(fsqe[i],cve2[i]))
        if(i==np.argmin(fe)):
            print("--> BEST")
        print("-"*40)
        print("PCF        : {0:.4f}".format(m[i]))
        if(args.verbosity>0):
            print("Mean error : {0:.4f}".format(e[i]))
        print("STD        : {0:.4f}".format(std[i]))
        print("-"*40)
        if(abs(corepart-1.0)>0.0000001):
            print("Lifetime   : {0:.4f} <-core-corrected".format(corepart*mt[i]))
        else:
            print("Lifetime   : {0:.4f}".format(mt[i]))
        if(args.verbosity>0):
            print("Mean error : {0:.4f}".format(corepart*e[i]))
        print("STD        : {0:.4f}".format(corepart*stdt[i]))
        if(args.verbosity>0):
            if(args.verbosity>1):
                print("Separate PCFs:")
                for j in range(gzeros.shape[0]):
                    print("- {0:.4f}".format(gzeros[j,i]))
            print("Separate lifetimes:")
            for j in range(lifetimes.shape[0]):
                print("- {0:.4f}".format(lifetimes[j,i]))
        print(" ")
