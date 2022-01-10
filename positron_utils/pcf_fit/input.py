import numpy as np

def get_imax(r,rmax):
    imax=r.shape[0]
    for i in range(len(r)):
        if(r[i]>rmax):
            imax=i
            return imax
    return imax

def extrapolate(args,g_vmc,g_dmc,N,metal,r):
    '''
    Takes in PCFs calculated with VMC and DMC, and
    cancels the first-order errors on the wavefunction
    by returning an array with elements g(x)=2*g_dmc(x)-g_vmc(x).
    Before extrapolation, the arrays are multiplied by N/(N-1)
    to assure an asymptotic behaviour g(x)=1 as x-> inf.
    '''

    if(metal==0):
        S=1
    else:
        S=N/(N-1)

    Npcf=g_vmc.shape[1]
        
    # THE LOGARITHM of extrapolated g
    arrglog=[]
    # Actual, extrapolated g
    arrg=[]
    # The r-arrays that are cleaned from offlier points
    r_extrapolated=[]
    # Extrapolate the QMC result
    for ipcf in range(Npcf):
        glog=np.zeros((len(g_vmc),))
        g=np.zeros((len(g_vmc),))
        gv=np.zeros((len(g_vmc),))
        gd=np.zeros((len(g_vmc),))
        rcopy=r
        for i in range(len(r)):
            gv[i]=S*g_vmc[i,ipcf]
            gd[i]=S*g_dmc[i,ipcf]
            g[i]=2*gd[i]-gv[i]
        i=0
        ind=0
        # Now go throught all the extrapolated values in g
        while(i<len(r)):
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
        arrg.append(g)
        arrglog.append(glog)
        r_extrapolated.append(rcopy)
        
    return arrg,arrglog,r_extrapolated


def get_input(args):

    r=np.load(args.rfile)
    g_vmc=np.load(args.vmcfile)
    g_dmc=np.load(args.dmcfile)
    if((len(r)!=len(g_vmc))or(len(r)!=len(g_dmc))):
        print("Error: input data is not consistent.")
        print(r.shape)
        print(g_vmc.shape)
        print(g_dmc.shape)
        sys.exit()
    if(args.omit_pcf[0]>-0.1):
        g_vmc=np.delete(g_vmc,args.omit_pcf,axis=1)
        g_dmc=np.delete(g_dmc,args.omit_pcf,axis=1)
    gex,glogex,r_ex=extrapolate(args,g_vmc,g_dmc,args.num_e,args.metal,r)

    return r,r_ex,gex,glogex
