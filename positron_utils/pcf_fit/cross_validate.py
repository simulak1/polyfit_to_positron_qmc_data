import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from fit import *
from input import get_input
from statistics import cross_validation_error

def get_args():

    args_parser = argparse.ArgumentParser()
    
    # Data files arguments
    args_parser.add_argument(
        '--vmcfile',
        help='Local path to vmc data',
        required=True,
        type=str
    )

    
    args_parser.add_argument(
        '--dmcfile',
        help='Local path to dmc data',
        required=True,
        type=str
    )

    args_parser.add_argument(
        '--rfile',
        help='Local path to dmc data',
        required=True,
        type=str
    )

    args_parser.add_argument(
        '--max-pol',
        help='Maximum polynomial order.',
        required=False,
        type=int,
        default=3
    )

    args_parser.add_argument(
        '--min-pol',
        help='Minimum polynomial order.',
        required=False,
        type=int,
        default=3
    )



    args_parser.add_argument(
        '--num-e',
        help='Number of electrons',
        required=False,
        type=int,
        default=3
    )

    args_parser.add_argument(
        '--fitscale',
        help='Scale the errors according to particle distance in the fitting.',
        required=False,
        type=int,
	default=0
    )
    
    args_parser.add_argument(
        '--fit-range-min',
        help='Minimum fit range',
        required=True,
        type=float
    )

    args_parser.add_argument(
	'--fit-range-max',
        help='Maximum fit range',
        required=True,
	type=float
    )

    args_parser.add_argument(
        '--fit-range-dx',
	help='Fit range interval',
        required=True,
	type=float
    )
    
    args_parser.add_argument(
        '--opt-method',
	help='lm,trf or dogbox',
	required=False,
        type=str,
        default='lm'
    )
    
    args_parser.add_argument(
        '--volume',
        help='Volume of the simulation cell.',
        required=True,
        type=float
    )

    args_parser.add_argument(
        '--metal',
        help='System is metallic (1) or not (0). Affects to the scaling of the extrapolated PCF. Default: 0.',
        required=False,
        type=int,
        default=0
    )
    
    args_parser.add_argument(
        '--verbosity',
        help='Controls the amount of info printed.',
        required=False,
        type=int,
        default=0
    )

    args_parser.add_argument(
	'--omit_pcf',
        help='skip pcfs.',
        required=False,
        nargs='+',
	type=int,
        default=[-1]
    )

    args_parser.add_argument(
        '--weight-file',
        help='give weights for the pcf histograms in the numpy array.',
	required=False,
        type=str,
        default='nofile'
    )
    
    args_parser.add_argument(
	'--valset1',
	help='Cross-validation set 1, by averaging.',
        required=False,
        nargs='+',
        type=int,
        default=[-1]
    )

    args_parser.add_argument(
	'--valset2',
        help='Cross-validation set 2, by averaging.',
        required=False,
        nargs='+',
        type=int,
        default=[-1]
    )

    args_parser.add_argument(
        '--explim',
        help='After fitting, we sometimes want to exponentiate. This is the upper limit for exp interval.',
        type=float,
        default=100.0
    )
    
    return args_parser.parse_args()

def printProgressBar(i,max,postText):
    n_bar =10 #size of progress bar
    j= i/max
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}")
    sys.stdout.flush()

def main():

    args=get_args()

    r,gex,ws=get_input(args)

    Nx=len(r)
    Npcf=len(gex)
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)

    rfit=np.arange(args.fit_range_min,args.fit_range_max,args.fit_range_dx)
    
    print("==== Cross-validation of PCF data fits ===")
    print(" ")
    print("Fit ranges to be tested: ")
    print(rfit)
    print("Polynomial orders to be tested:")
    print(p_degs)
    print(" ")
    
    MAE_matrix=np.zeros((rfit.shape[0],N_degs))
    MSE_matrix=np.zeros((rfit.shape[0],N_degs))
    
    irmax=get_imax(r,1.)
    ir1=get_imax(r,0.0)
    
    for ir in range(rfit.shape[0]):
        printProgressBar(ir,rfit.shape[0],"Validation progress")

        ri=rfit[ir]

        fits,logfits,glog,rex,opt_pol_coeff=do_fit(r,ri,gex,args)

        mae,mse=cross_validation_error(fits,gex,r,args)
        
        MAE_matrix[ir,:]=np.array(mae)
        MSE_matrix[ir,:]=np.array(mse)

                                       
    print("\n Number of PCF histograms: {} ".format(Npcf))
    print(" ")

    print("Cross-validated mean-average errors on the fits in the reange 0-1 au.:")
    print(" ")

    str1="     | "
    for deg in p_degs:
        str1+="     {}      |".format(deg)

    print(str1)
    print(50*"-")

    for ir in range(rfit.shape[0]):
        ri=rfit[ir]
        s=" {} | ".format(ri)
        for ideg in range(N_degs):
            if(np.abs(MAE_matrix[ir,ideg]-MAE_matrix.min())<0.0000000001):
                s+=" *{:.6f}* |".format(MAE_matrix[ir,ideg])
            else:
                s+="  {:.6f}  |".format(MAE_matrix[ir,ideg])
        print(s)
    print(50*"-")

    print("Cross-validated mean-square errors on the fits in the reange 0-1 au.:")
    print(" ")

    str1="     | "
    for deg in p_degs:
        str1+="     {}      |".format(deg)

    print(str1)
    print(50*"-")

    for	ir in range(rfit.shape[0]):
        ri=rfit[ir]
        s=" {} | ".format(ri)
        for ideg in range(N_degs):
            if(np.abs(MSE_matrix[ir,ideg]-MSE_matrix.min())<0.0000000001):
                s+=" *{:.6f}* |".format(MSE_matrix[ir,ideg])
            else:
                s+="  {:.6f}  |".format(MSE_matrix[ir,ideg])
        print(s)
    print(50*"-")
    
    
if __name__ == '__main__':
    main()

