import numpy as np
import sys
import warnings
import argparse
from fit import do_fit
from output import *
from statistics import *
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
        '--fitscale',
        help='Scale the errors according to particle distance in the fitting.',
        required=False,
        type=int,
        default=-1
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

    args_parser.add_argument(
        '--weight-file',
        help='give weights for the pcf histograms in the numpy array.',
	required=False,
        type=str,
	default='nofile'
    )

    args_parser.add_argument(
        '--wtot',
	help='total weight.',
        type=float,
        default=-1.0
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
                    
def main():

    args=get_args()

    # PCF matrices and weights
    r,gex,ws=get_input(args)

    print("Number of PCFs: {}".format(len(gex)))

    # Fitting range
    r_range=args.fit_range*args.lat_vec
    
    # Fitting
    fits,logfits,glog,rex,opt_pol_coeff=do_fit(r,r_range,gex,args)

    # Fitting statistics and plotting. NOTE! plotting should be done in separate function
    fe,fsqe,cve,cve2=fit_statistics(args,fits,gex,r,r_range)    
    
    # Get g(0) values and statistics
    gzeros=fits[0,:,:]

    m,e,std=mean_and_error(gzeros,ws)
    
    # Lifetime statistics
    coeff1=100.617
    coeff2=100.93952105134674
    lifetimes=1000.0*(coeff1/2*args.num_e/args.volume*gzeros)**-1
    mt,et,stdt=mean_and_error(lifetimes,ws)

    if(args.table==1):
        make_table(args,m,mt,fe,fsqe,e,std,stdt,gzeros,lifetimes)
    else:
        print_output(args,m,mt,fe,fsqe,cve,cve2,e,std,stdt,gzeros,lifetimes)

    if args.plot == 1:
        plot_results(args,fits,logfits,gex,glog,r,rex,opt_pol_coeff)
        
    sys.exit('All done.')

if __name__ == '__main__':
    main()
        
