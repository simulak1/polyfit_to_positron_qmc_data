import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from fit import *
from input import get_input

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
    
    return args_parser.parse_args()

def main():

    args=get_args()

    r,r_ex,gex,glogex=get_input(args)

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
    
    MSE_matrix=np.zeros((rfit.shape[0],N_degs))

    irmax=get_imax(r,2.0)
    
    for ir in range(rfit.shape[0]):
        ri=rfit[ir]

        fits,opt_pol_coeff=do_fit(r,r_ex,ri,glogex,args)
        
        for ideg in range(N_degs):
            print("Validating fit range {}, polynomial deg. {}".format(ri,p_degs[ideg]))
            mse=[]
            for itraining in range(Npcf):
                fit=np.exp(fits[:irmax,itraining,ideg])
                for ivalidation in range(Npcf):
                    if(ivalidation==itraining):
                        continue
                    mse.append(np.mean((fit-gex[ivalidation][:irmax])**2))

            MSE_matrix[ir,ideg]=sum(mse)/len(mse)

    print("Number of PCF histograms: {} ".format(Npcf))
    print(" ")
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

