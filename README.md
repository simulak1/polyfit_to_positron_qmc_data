# polyfit_to_positron_qmc_data

A program to fit N:th order polynomial to quantum Monte Carlo-based, translationally and rotationally averaged positron-electron correlation functions (PCFs).

## Folders 

* positron_utils/pcf_fit
Contains the source code of the program. To print program info, run `python positron_utils/pcf_fit/plyfit.py -h`
* examples
Example data of silicon and lithium with different quantum wave function approximations. Each folder within `examples` contains 4 files: radial grid r_1.npy, variational Monte Carlo PCF data gvmc_1.npy, diffusion Monte Carlo PCF data gdmc_1.npy and a weight file for the relative weights of the twists corresponding to single PCFs.
* scripts
Helpful scripts to run the main program. Information for running the scripts can be found from within the scripts analyse.sh and validate.sh, and from scripts/README
* extra
Some extra things about the program

## Quick start to analyse example data

### Create local environment

1. Clone the repository
2. Go to repository
#### Pip

3.1 `pip install virtualenv`
3.1 `python -m venv virtualenv; source virtualenv/bin/activate`
3.3 `pip install -r requirements.txt`

#### Anaconda

3. `conda env create -f environment.yml`

### Validate the fitting parameters, fit, and compute the lifetime estimates for bulk silicon

4. cd scripts
5. ./validate.sh ../examples/silicon_SJ 2 8 0.5 7
* *The best parameter combinations are either with polynomial order 3 and fitting range of 2 Bohr or order 5 and range 4.5* *
6. ./analyse.sh ../examples/silicon_SJ 8 0.270107093552E+03 4.4 5 0.000155019074496
* *I get lifetime of 238.8ps with the 5:th-order polynomial

## Program background and theory

### Background

This program has been the result of development project on implementing positron simulation capabilities to the quantum Monte Carlo simulation package CASINO (https://vallico.net/casinoqmc/). The program has been detached from the CASINO's positron simulation branch (not public, but under development), and the development history has been copied from the branch manually. Thus the history has the correct chronological order, but incorrect dates and times for individual commits (except from 10.1.2022 onwards the commit times are correct). This program has been intensively used in the soon to be published article 'Quantum Monte Carlo Simulation of Positron Lifetimes in Solids' (url will be copied here when the article is published).

### Theory and numerics

Positron lifetime calculation with quantum Monte Carlo requires a well-defined wave function with which to accumulate positron-electron PCF. In CASINO manual (https://casinoqmc.net/casino_manual_dir/casino_manual.pdf), the PCF is defined in equation (398). For positron lifetimes, a rotationally and translationally averaged PCFs (Eq. (407)) must be computed by collecting the positron-electron distances in bins at each Metropolis step used in statistic accumulation (the procedure is explained in the CASINO manual section 34.3.5). With the obtained PCF, we can estimate the positron annihilation rate in the sample by taking the value of PCF at zero, g(0), multiplying by the electron density in the simulation cell and by a constant prefactor defined by quantum electrodynamics to define the dependence of the annihilation rate to the cross-section of positron and electron at contact. Here the prefactor that is being used is ~100.6 The positron lifetime is the inverse of the annihilation rate. 

The accumulated PCF data is very noisy as it approaches zero. The strenght of the noise is dependent on 1/r, where r is the positron-electron distance. This behaviour leads to both strong noise near zero distance and divergences. Therefore the statistics at zero is bad-behaved and needs to be improved. To improve it, an N:th-order polynomial can be fitted against the logarithm of the PCF data. For a polynomial form a0+a1*x+a2*x^2+...+aN*x^N, we require a_1=-1 to satisfy the Kimball cusp conditions (the cusp conditions also require the polynomial to be fitted to log of the PCF). Then we must choose the fitting range and the order of the fitting polynomial. Then the value exp(a0) can be estimated as the value of g(0).

To improve statistics further, PCF data should be accumulated with multiple independent QMC runs. This results to the numpy arrays of the form (Npcf,Nx) found in the examples-directory, where Npcf is the number of independently accumulated PCF functions, and Nx is the number of radial points in the accumulated PCFs. The radial values are given in the separate r_1.npy-files. Then the fitting can be done separately to the individual PCFs, resulting to Npcf values exp(a0^i), and we can take the average of the obtained lifetime estimates. 

The use of independently accumulated PCFs also enable cross-validation of different fitting parameters. For example, if Nx=16, we can take 2 sets of averages of 8 individual PCF functions, fit polynomials with different parameters to the logarithm of the another set and compute mean-square-errors of the residuals between the obtained fits and the PCF values not used in the fit. This way we can obtain independent estimates of the model performances. 

When the program is being run to produce lifetime estimates, 3 kinds of errors are printed to output (with verbosity=0).

1. The errors of the lifetime and PCF values obtained with different degrees of polynomials. These errors are just the standard deviations of the independently obtained NPCF lifetime values and g(0):s. 
2. The mean-square-errors of the average fit (averaged over all of the separate fits) against the average PCF (averaged over the separate PCFs).
3. The mean-square-cross-validation-error. Given Npcf PCF functions, we can divide them to e.g. 2 sets of same number of different PCFs averaged, fit with some fitting range and polynomial order to both averaged PCFs, and then compute 2 values of mean-square-errors obtained by taking residuals of the fits against the averaged PCFs not used in the fit. 





