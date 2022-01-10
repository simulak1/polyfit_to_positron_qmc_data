import numpy as np
import sys
import warnings
import argparse
from fit import do_fit
from output import *
from statistics import *
from input import get_input,get_imax
from cross_validate import cross_validation
warnings.simplefilter("error")

def get_args():

    args_parser = argparse.ArgumentParser()

    # Data files arguments

    args_parser.add_argument(
	'--task',
        help='''
        1: Polynomial fit
        2: Cross-validation
        3: Both, 1 is performed based on output from 2.
        ''',
        required=False,
	type=int,
        default=1
    )

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
        '--fit-range-min',
        help='Minimum fit range', 
        required=False,
        type=float
    )                                                                                   

    args_parser.add_argument(
	'--fit-range-max',
	help='Maximum fit range',
        required=False,
        type=float
    )

    args_parser.add_argument(
        '--fit-range-dx',
        help='Fit range interval',
	required=False,
        type=float
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
        required=False,
        type=float,
        default=1
    )


    args_parser.add_argument(
        '--fit-range',
        help='Maximum polynomial order.',
        required=False,
        type=float,
        default=1
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

    print(" ")
    print("#################################################")
    print("# Polynomial fit of positron-electron PCF       #")
    print("#################################################")
    print(" ")
    print("System:")
    print("-------")
    print("Volume                 : {} ".format(args.volume))
    print("Number of electrons    : {}".format(args.num_e))
    print("Is metal (0=no, 1=yes)?: {} ".format(args.metal))
    if(args.weight_file != "nofile"):
        print("Twist averaged data, weight file : "+args.weight_file)
    print(" ")
    print("Loaded: ")
    print(args.vmcfile)
    print(args.dmcfile)
    print(args.rfile)
    print("Shape of PCF files: [{},{}]".format(len(gex),gex[0].shape))
    print(" ")
    p_degs=np.arange(args.min_pol,args.max_pol+1,2)
    if(args.task==1):
        print("Optimization: ")
        print("-------")
        r_range=args.fit_range*args.lat_vec
        print("Fitting range             : {}".format(r_range))
        print("Maxmimum polynomial order : {}".format(args.max_pol))
        
    elif(args.task==2):
        rfit=np.arange(args.fit_range_min,args.fit_range_max,args.fit_range_dx)

        print("==== Cross-validation of PCF data fits ===")
        print(" ")
        print("Fit ranges to be tested: ")
        print(rfit)
        print("Polynomial orders to be tested:")
        print(p_degs)
        print(" ")
    elif(args.task==3):
        rfit=np.arange(args.fit_range_min,args.fit_range_max,args.fit_range_dx)

        print("==== Cross-validation of PCF data fits ===")
        print(" ")
        print("Fit ranges to be tested: ")
        print(rfit)
        print("Polynomial orders to be tested:")
        print(p_degs)
        print(" ")
        print("Polynomial fitting is performed on top of the cross-validation")
    else:
        sys.exit("Unknown task, give option 1-3.")
        
    print("Optimization method       :"+args.opt_method)
    if(args.fitscale>0):
        print("Weighting of the residuals applied")

    if(args.task>1):
        r_range,minpol=cross_validation(args,r,gex,ws,rfit)
        p_degs=np.arange(minpol,minpol+1,2)

    if(args.task==1 or args.task==3):
        # Fitting
        fits,logfits,glog,rex,opt_pol_coeff=do_fit(r,r_range,p_degs,gex,args)

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

    if(args.plot==1):
        plt.show()
        
    sys.exit('All done.')

if __name__ == '__main__':
    main()
        
