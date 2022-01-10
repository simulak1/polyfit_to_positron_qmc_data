import numpy as np
import matplotlib.pyplot as plt

def plot_results(p_degs,N_degs,g_average,fit_average,r,imax):

    plt.plot(r,g_average,'k-',linewidth=2,label='Extrapolated data')
    for i in range(N_degs):
        plt.plot(r[:imax],np.exp(fit_average[:imax,i]),label=str(p_degs[i])+'-order pol.')
    plt.grid()
    plt.xlabel('Positron-electron distance (Bohr) vai hartree?')
    plt.ylabel('Pair correlation function')
    plt.title('Pair correlation function of positron-electron pairs in diamond-phase silicon')
    plt.legend()
    plt.show()
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

def print_output(args,m,mt,fe,fsqe,e,std,stdt,gzeros,lifetimes):
    corepart=args.corepart
    for i in range(m.shape[0]):
        deg=args.min_pol+i*2
        print("="*40)
        print("Polynomial degree: {}".format(deg))
        print("Mean error: {0:.5f}, Mean squared error: {0:.5f}".format(fe[i],fsqe[i]))
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
