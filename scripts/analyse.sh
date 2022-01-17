 #!/bin/bash

# Moodify this to point to the src-directory
#SRC=<path-to-source-dir>
SRC=../positron_utils/pcf_fit

if [ "$#" -ne 7 ]; then
    echo " "
    echo "Usage: "
    echo "./$0 <folder> <num-e> <volume> <fit_range> <max-pol> <corepart> <element>, "
    echo "where"
    echo "<folder>:      found from local folder as system_size/opt_level,"
    echo "<num-e>: Number of electrons"
    echo "<volume>: Volume"
    echo "<fit_range>  : specified currently as fraction of the lattice constant"
    echo "<max-pol>    : maximum polynomial order"
    echo "<corepart>    : core electron contribution modifier"
    echo "<element>    : element selected"
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

# These are the fractional contributions of the core electrons to the
# total annihilation rate in silicon, computed with different positron
# correlation functionals using the Atsup simulation package. 
si_corepart_ND=0.000162020075261
si_corepart_BN=0.000167294203076
si_corepart_B15=0.000155019074496
si_corepart_B95=0.000110971904772
si_corepart_KUR=0.00014124373212

# These are the coreparts for lithium
li_corepart_ND=0.00044730984928999994
li_corepart_BN=0.000455641003118
li_corepart_B15=0.000352885464408
li_corepart_B95=0.000214330966927
li_corepart_KUR=0.000334737553995

echo $6 $7

if [ $6 = "ECP"  ]; then
    corepart=1

elif [ $6 = "ND" ] && [ $7 = "Si"  ]; then
    corepart=$si_corepart_ND

elif [ $6 = "BN" ] && [ $7 = "Si"  ]; then
    corepart=$si_corepart_BN

elif [ $6 = "B15" ] && [ $7 = "Si"  ]; then
    corepart=$si_corepart_B15

elif [ $6 = "B95" ] && [ $7 = "Si"  ]; then
    corepart=$si_corepart_B95

elif [ $6 = "KUR" ] && [ $7 = "Si"  ]; then
    corepart=$si_corepart_KUR

elif [ $6 = "ND" ] && [ $7 = "Li"  ]; then
    corepart=$li_corepart_ND

elif [ $6 = "BN" ] && [ $7 = "Li"  ]; then
    corepart=$li_corepart_BN

elif [ $6 = "B15" ] && [ $7 = "Li"  ]; then
    corepart=$li_corepart_B15

elif [ $6 = "B95" ] && [ $7 = "Li"  ]; then
    corepart=$li_corepart_B95

elif [ $6 = "KUR" ] && [ $7 = "Li"  ]; then
    corepart=$li_corepart_KUR

else 
    echo "Core contribution option not found."
    exit 5
fi


# Silicon primitive cell volume 
V=0.270107093552E+03
# Number of electrons in silicon primitive cell with AREP pseudos
N=8
# Number of electrons in silicon primitive cell with ECP pseudos
N=24

# Lithium primitive cell volume
V=0.787919086440E+04
# Number of electrons in lithium primitive cell with AREP pseudos
N=54
# Number of electrons in all-electron simulation
N=162

python $SRC/plyfit.py --vmcfile $gv --dmcfile $gd --rfile $r --max-pol $5 --min-pol 3 --plot 1 --lat-vec 1 --fit-range $4 --num-e $2 --volume $3 --verbosity 0 --metal 0 --corepart $corepart  --table 0 --weight-file $1/weights.txt --wtot=128 --explim 4.5 --opt-method='lm' --fitscale 1

echo "Core correction done via" $6 core contribution=$corepart


