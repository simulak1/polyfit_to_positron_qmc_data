 #!/bin/bash

# Write the path to the src-directory here
SRC=../src

if [ "$#" -ne 5 ]; then
    echo " "
    echo "Usage: "
    echo "./$0 <folder> <pcf_reblock> <fit_range-min> <fit_range-max> <fit_range-dx> <max-pol> , "
    echo "where"
    echo "<folder>:      found from local folder as system_size/opt_level,"
    echo "<fit_range-min>  : Minimum fit range"
    echo "<fit_range-max>  : Maximum fit range"
    echo "<fit_range-dx>   : Fit range interval"
    echo "<max-pol>        : maximum polynomial order"
    
    echo " "
    exit
fi

if [ ! -d ${1} ]; then
    echo "Error: PCF directory not found"
    exit 1
fi

r=$1/r_1.npy
gv=$1/gvmc_1.npy
gd=$1/gdmc_1.npy

if [ ! -f $r  ]; then
    echo "Error: $r not found."
    exit 2
fi

if [ ! -f $gv  ]; then
    echo "Error: $r not found."
    exit 3
fi

if [ ! -f $gd  ]; then
    echo "Error: $r not found."
    exit 4
fi


python $SRC/plyfit.py --task 2 --vmcfile $gv --dmcfile $gd --rfile $r --max-pol $5 --min-pol 3  --num-e 1 --fit-range-min $2 --fit-range-max $3 --fit-range-dx $4 --volume 2 --verbosity 0 --metal 0 --plot 1 --explim 1.0 --opt-method 'lm' --fitscale 1 --valset1 0 1 2 3 4 5 6 7 --valset2 8 9 10 11 12 13 14 15 


