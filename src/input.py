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
        
    # Actual, extrapolated g
    arrg=[]
    # Extrapolate the QMC result
    for ipcf in range(Npcf):
        g=np.zeros((len(g_vmc),))
        gv=np.zeros((len(g_vmc),))
        gd=np.zeros((len(g_vmc),))
        for i in range(len(r)):
            gv[i]=S*g_vmc[i,ipcf]
            gd[i]=S*g_dmc[i,ipcf]
            g[i]=2*gd[i]-gv[i]
        i=0
        ind=0
        # Now go throught all the extrapolated values in g
        arrg.append(g)
        
        
    return arrg


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
    gex=extrapolate(args,g_vmc,g_dmc,args.num_e,args.metal,r)

    if(args.weight_file=='nofile'):
        ws=np.ones((g_vmc.shape[1],))
        wtot=np.sum(ws)
    else:
        print("-------------------------------------------------------------------------")
        print("WARNING: You are using the weights-file "+args.weight_file+".")
        print("         the total weigth should be # of twists X # simulations/twist. ")
        print("-------------------------------------------------------------------------")
        if(args.wtot<0):
            sys.exit("You must give total weight when a weight file is present")
        wtot=args.wtot
        ws=np.loadtxt(args.weight_file)
        ws=wtot*ws

    
    return r,gex,ws
