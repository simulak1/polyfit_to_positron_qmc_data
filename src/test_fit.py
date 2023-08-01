import argparse
import numpy as np
from fit import do_fit,fit_func,polynomial_for_kimball3,polynomial_for_kimball5,polynomial_for_kimball7
from statistics import mean_and_error,fit_statistics


def get_args1():

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

def get_args2():

    args_parser = argparse.ArgumentParser()

    args_parser.add_argument(
        '--fitscale',
        type=float,
        default=.4
    )

    args_parser.add_argument(
        '--lat-vec',
        type=float,
        default=.4
    )

    args_parser.add_argument(
        '--opt-method',
        type=str,
        default='lm'
    )

    args_parser.add_argument(
        '--explim',
        type=float,
	default=1
    )

    args_parser.add_argument(
	'--verbosity',
        type=int,
	default=0.
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
	'--max-pol',
        help='Maximum polynomial order.',
	required=False,
        type=int,
        default=7
    )



    args_parser.add_argument(
        '--min-pol',
        help='Minimum polynomial order.',
	required=False,
        type=int,
        default=3
    )


    return args_parser.parse_args()

def simple_function(x):
    return 1.-x+2.*x**3

def function_with_noise1(x,noise=.0):
    y=9-x-1.5*x**2+.35*x**3
    if(noise>0.0001):
        a=1./x*np.random.normal(.0,noise,x.shape[0])
        return y+a
    else:
        return y

def function_with_noise2(x,noise=.0):
    y=5-x-.3*x**2-.3*x**3+.04*x**4+.012*x**5
    if(noise>0.0001):
        a=1./x*np.random.normal(.0,noise,x.shape[0])
        return y+a
    else:
        return y

def function_with_noise3(x,noise=.0):
    a0=5
    a1=-1.
    a2=-4.
    a3=-3.4 #1.                                                                                      
    a4=2.9 #.00000001                                                                                
    a5=-0.35
    a6=-.054
    a7=.0082
    y=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7

    if(noise>0.0001):
        a=1./x*np.random.normal(.0,noise,x.shape[0])
        return y+a
    else:
        return y
    
def test_simple_function():

    args=get_args1()
    
    x=np.arange(-5,5,.1)
    r=np.array([-2,1])

    y=simple_function(x)
    fit,popt=fit_func(r,x,y,polynomial_for_kimball3,-1,args)

    assert abs(fit[0]+13.) < .01
    assert(abs(fit[1]-2.)) < .01
    assert abs(popt[0]-1)  < .01
    assert abs(popt[1])    < .01
    assert abs(popt[2]-2.) < .01

def test_function_with_noise1():
    
    x1=np.arange(.01,5.2,.02)
    x2=np.arange(.0000001,5.1,.02)
    
    y1=function_with_noise1(x1,1.)
    y2=function_with_noise1(x2)

    args=get_args2()
    fit1,popt1=fit_func(x2,x1,y1,polynomial_for_kimball3,-1,args)
    mse=np.mean((fit1-y2)**2)

    assert mse<0.05
    assert abs(popt1[0]-9)   <.5
    assert abs(popt1[1]+1.5) <.5
    assert abs(popt1[2]-.35) <.05

def test_function_with_noise2():

    x1=np.arange(.01,5.2,.02)
    x2=np.arange(.0000001,5.1,.02)

    y1=function_with_noise2(x1,1.)
    y2=function_with_noise2(x2)

    args=get_args2()
    fit1,popt1=fit_func(x2,x1,y1,polynomial_for_kimball5,-1,args)
    mse=np.mean((fit1-y2)**2)

    assert mse<0.08
    assert abs(popt1[0]-5)   <1.

def test_function_with_noise3():

    x1=np.arange(.01,5.2,.02)
    x2=np.arange(.0000001,5.1,.02)

    y1=function_with_noise3(x1,1.)
    y2=function_with_noise3(x2)

    args=get_args2()
    fit1,popt1=fit_func(x2,x1,y1,polynomial_for_kimball7,-1,args)
    mse=np.mean((fit1-y2)**2)

    assert mse<0.15
    assert abs(popt1[0]-5)   <2.0    

def test_do_fit_for_f1():

    x1=np.arange(.01,5.2,.02)

    y_clean=function_with_noise1(x1)

    args=get_args2()
    y=[]
    for i in range(16):
        y1=np.exp(function_with_noise1(x1))
        a=1./x1*np.random.normal(.0,1,x1.shape[0])
        y.append(y1+a)

    fits,logfits,glog,rex,cpol=do_fit(x1,5,[3,5,7],y,args)

    zeros=np.log(fits[0,:,:])
    ws=np.ones((16,))
    m,e,std=mean_and_error(zeros,ws)
    
    assert abs(m[0]-9) < .05
    assert abs(m[1]-9) < .05
    assert abs(m[2]-9) < .05


def test_do_fit_for_f2():

    x1=np.arange(.01,5.2,.02)

    y_clean=function_with_noise2(x1)

    args=get_args2()
    y=[]
    for i in range(16):
        y1=np.exp(function_with_noise2(x1))
        a=1./x1*np.random.normal(.0,1,x1.shape[0])
        y.append(y1+a)

    fits,logfits,glog,rex,cpol=do_fit(x1,1,[3,5,7],y,args)

    zeros=np.log(fits[0,:,:])
    ws=np.ones((16,))
    m,e,std=mean_and_error(zeros,ws)

    assert abs(m[0]-5) < .05
    assert abs(m[1]-5) < .5
    assert abs(m[2]-5) < .1

def test_do_fit_for_f3():

    x1=np.arange(.01,5.2,.02)

    y_clean=function_with_noise3(x1)

    args=get_args2()
    y=[]
    for i in range(16):
        y1=np.exp(function_with_noise3(x1))
        a=1./x1*np.random.normal(.0,1,x1.shape[0])
        y.append(y1+a)

    fits,logfits,glog,rex,cpol=do_fit(x1,1,[3,5,7],y,args)

    zeros=np.log(fits[0,:,:])
    ws=np.ones((16,))
    m,e,std=mean_and_error(zeros,ws)

    assert abs(m[0]-5) < 1.7
    assert abs(m[1]-5) < 3
    assert abs(m[2]-5) < 4
