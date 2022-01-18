import argparse
import numpy as np
from fit import fit_func,polynomial_for_kimball3

def get_args():

    args_parser = argparse.ArgumentParser()

    args_parser.add_argument(
        '--fitscale',
        type=float,
        default=-1.
    )

    args_parser.add_argument(
        '--opt-method',
        type=str,
        default='lm'
    )

    return args_parser.parse_args()

def test_fit_func():

    args=get_args()
    
    '''
    Test to fit 1+-x+x**3
    '''
    x=np.arange(-15,15,.1)
    y=1.-x+2.*x**3
    r=np.array([-2,1])
    fit,popt=fit_func(r,x,y,polynomial_for_kimball3,4,args)
    print(popt)
    assert abs(fit[0]+13.)<.1
