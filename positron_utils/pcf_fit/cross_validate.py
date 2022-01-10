import numpy as np
import matplotlib.pyplot as plt
import sys
from fit import *
from input import get_input
from statistics import cross_validation_error

def printProgressBar(i,max,postText):
    n_bar =10 #size of progress bar
    j= i/max
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}")
    sys.stdout.flush()

def cross_validation(args,r,gex,ws,rfit):

    Nx=len(r)
    Npcf=len(gex)
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)

    MAE_matrix=np.zeros((rfit.shape[0],N_degs))
    MSE_matrix=np.zeros((rfit.shape[0],N_degs))
    
    irmax=get_imax(r,1.)
    ir1=get_imax(r,0.0)
    
    for ir in range(rfit.shape[0]):
        printProgressBar(ir,rfit.shape[0],"Validation progress")

        ri=rfit[ir]

        fits,logfits,glog,rex,opt_pol_coeff=do_fit(r,ri,p_degs,gex,args,True)
        
        mae,mse=cross_validation_error(fits,np.array(gex)[:,:get_imax(r,args.explim)],r,args)
        
        MAE_matrix[ir,:]=np.array(mae)
        MSE_matrix[ir,:]=np.array(mse)

                                       
    print(50*'-')
    print("\n Cross-validated mean-square errors on the fits in the selected ranges:")
    print(" ")

    str1="          | "
    for deg in p_degs:
        str1+="     {}      |".format(deg)

    print(str1)
    print(50*"-")

    for	ir in range(rfit.shape[0]):
        ri=rfit[ir]
        s=" {} Bohr | ".format(ri)
        for ideg in range(N_degs):
            if(np.abs(MSE_matrix[ir,ideg]-MSE_matrix.min())<0.0000000001):
                s+=" *{:.6f}* |".format(MSE_matrix[ir,ideg])
                irmin=ri
                idegmin=p_degs[ideg]
            else:
                s+="  {:.6f}  |".format(MSE_matrix[ir,ideg])
        print(s)
    print(50*"-")

    return irmin,idegmin,MSE_matrix
