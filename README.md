# polyfit_to_positron_qmc_data

A program to fit N:th order polynomial to quantum Monte Carlo-based, translationally and rotationally averaged positron-electron correlation functions (PCFs).

## Folders 

* positron_utils/pcf_fit
Contains the source code of the program. To print program info, run 'python positron_utils/pcf_fit/plyfit.py -h' 
* examples
Example data of silicon and lithium with different quantum wave function approximations. Each folder within 'examples' contains 4 files: radial grid r_1.npy, variational Monte Carlo PCF data gvmc_1.npy, diffusion Monte Carlo PCF data gdmc_1.npy and a weight file for the relative weights of the twists corresponding to single PCFs.
* scripts
Helpful scripts to run the main program. Information for running the scripts can be found from within the scripts analyse.sh and validate.sh, and from scripts/README
* extra
Some extra things about the program

## Quick start to analyse example data

### Create local environment

1) Clone the repository
2) Go to repository
3) 'pip install virtualenv'
4) 'python -m venv virtualenv; source virtualenv/bin/activate'
5) pip install -r requirements.txt

### Validate the fitting parameters, fit, and compute the lifetime estimates for bulk silicon

6) cd scripts
7) ./validate.sh ../examples/silicon_SJ 2 8 0.5 7
* *The best parameter combinations are either with polynomial order 3 and fitting range of 2 Bohr or order 5 and range 4.5* *
8) ./analyse.sh ../examples/silicon_SJ 8 0.270107093552E+03 4.4 5
* *I get lifetime of 238.8ps with the 5:th-order polynomial

## Program background and theory

### Background

This program has been the result of development project on implementing positron simulation capabilities to the quantum Monte Carlo simulation package CASINO (https://vallico.net/casinoqmc/). The program has been detached from the CASINO's positron simulation branch (not public, but under development), and the development history has been copied from the branch manually. Thus the history has the correct chronological order, but incorrect dates and times for individual commits (except from 10.1.2022 onwards the commit times are correct). This program has been intensively used in the soon to be published article 'Quantum Monte Carlo Simulation of Positron Lifetimes in Solids' (url will be copied here when the article is published).

### Theory

Positron lifetime calculation with quantum Monte Carlo requires a well-defined wave function with which to accumulate positron-electron PCF. In CASINO manual (https://casinoqmc.net/casino_manual_dir/casino_manual.pdf), the PCF is defined in equation (398). For positron lifetimes, a rotationally and translationally averaged PCFs (Eq. (407)) must be computed by collecting the positron-electron distances in bins at each Metropolis step used in statistic accumulation (the procedure is explained in the CASINO manual section 34.3.5).§  
