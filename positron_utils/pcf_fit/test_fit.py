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
    Test to fit 1+x**3
    x=-1,0,1,2
    y=0,1,2,9
    r=-2,1
    returned=-7,2 ??
    '''
    x=np.array([-1.,0.,1.,2.])
    y=np.array([0.,1.,2.,9.])
    r=np.array([-2,1])
    fit,popt=fit_func(r,x,y,polynomial_for_kimball3,4,args)
    assert abs(fit[0]+7.)<.1
