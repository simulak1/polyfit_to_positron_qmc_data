import numpy as np
import sys
import warnings
import argparse
from fit import do_fit
from output import *
from input import get_input,get_imax
warnings.simplefilter("error")

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
        '--plot',
        help='Plot PCF(s)?',
        required=False,
        type=int,
        default=0
    )
    
    args_parser.add_argument(
        '--num-e',
        help='Number of electrons',
        required=False,
        type=int,
        default=3
    )
    
    args_parser.add_argument(
        '--lat-vec',
        help='Maximum polynomial order.',
        required=True,
        type=float
    )


    args_parser.add_argument(
        '--fit-range',
        help='Maximum polynomial order.',
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
        '--corepart',
        help='tau_valence/tau_total.',
        type=float,
        default=1.0
    )

    args_parser.add_argument(
	'--table',
        help='Printout table in latex.',
        required=False,
	type=int,
        default=[-1]
    )

    args_parser.add_argument(
        '--cross-val',
        help='Perform a cross-validation over the fits.',
        required=False,
        type=int,
        default=[-1]
    )

    return args_parser.parse_args()
                    
def fit_statistics(args,fits,g,r):

    # The dimensions of our current experiment
    Nx=len(r)
    Npcf=len(g)
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    N_degs = len(p_degs)
    bond_distance=args.fit_range*args.lat_vec
    
    # Averages of the PCFs and fits
    fit_average=np.zeros((Nx,N_degs))
    g_average=np.zeros((Nx,))
    for i in range(Npcf):
        g_average[:]=g_average[:]+g[i]
        for d in range(N_degs):
            fit_average[:,d]=fit_average[:,d]+fits[:,i,d]

    g_average=g_average/Npcf
    fit_average=fit_average/Npcf
            
    # Compute the fitting error
    fit_errors=np.zeros((N_degs,))
    fit_sqerrors=np.zeros((N_degs,))
    imax=get_imax(r,2)
    for i in range(N_degs):
        for j in range(imax):
            fit_errors[i]=fit_errors[i]+np.absolute(np.exp(fit_average[j,i])-g_average[j])
            fit_sqerrors[i]=fit_sqerrors[i]+(np.exp(fit_average[j,i])-g_average[j])**2
    fit_errors=fit_errors/imax
    fit_sqerrors=fit_sqerrors/imax

    
    if args.plot == 1:
        plot_results(p_degs,N_degs,g_average,fit_average,r,imax)

    return fit_errors,fit_sqerrors

def mean_and_error(g):
    '''
    Input: (N_pcf,N_degs)-order matrix of PCF-zeroes.
    Output: 
      * Means over separate PCF's.
      * Mean average errors of PCF's.
      * Standard deviations of PCF's.

    
    '''
    gmean=np.mean(g,axis=0)
    mean_error=np.mean(np.absolute(gmean-g),axis=0)
    standard_deviation=np.std(g,axis=0)/np.sqrt(g.shape[0])#np.sqrt(np.mean((gmean-g)**2,axis=0))/np.sqrt(g.shape[0])
        
    return gmean,mean_error,standard_deviation


def main():

    args=get_args()
    
    r,r_ex,gex,glogex=get_input(args)

    r_range=args.fit_range*args.lat_vec
    
    # Fitting
    fits,opt_pol_coeff=do_fit(r,r_ex,r_range,glogex,args)

    # Fitting statistics
    fe,fsqe=fit_statistics(args,fits,gex,r)
    
    # Get g(0) values and statistics
    gzeros=np.exp(fits[0,:,:])
    m,e,std=mean_and_error(gzeros)
    
    # Lifetime statistics
    lifetimes=1000.0*(100.617/2*args.num_e/args.volume*gzeros)**-1
    mt,et,stdt=mean_and_error(lifetimes)

    if(args.table==1):
        make_table(args,m,mt,fe,fsqe,e,std,stdt,gzeros,lifetimes)
    else:
        print_output(args,m,mt,fe,fsqe,e,std,stdt,gzeros,lifetimes)
        
    sys.exit('All done.')

if __name__ == '__main__':
    main()
        
