#!/bin/bash


#### Vapor names, assuming that the cluster set file is named according to them ####
vapors=("A" "N")

#### Include ions? ####
l_incl_ions=1

#### Always use constant vapor concentrations? ####
# Constant concentrations can also be set in ACDC when l_const_vapor=0, but compilation may be slower for complex systems
l_const_vapor=0

#### Fixed T and RH values - comment out when using varying input T ####
#temperature=280
#rh=20

#### Suffix of the cluster system files ####

suffix="_example"


# End of the input section
##################################################################################################

# Create the Perl option string

perl_opt=""
cluster_file="input_"
vapor_suffix="_"

if [ -n "$temperature" ]; then
    perl_opt+=" --temperature $temperature"
else
    perl_opt+=" --variable_temp"
fi

if [ -n "$rh" ]; then
    perl_opt+=" --rh $rh"
fi

for vapor in "${vapors[@]}"; do
    perl_opt+=" --cs_only 1$vapor,0"
    [ $l_const_vapor -eq 1 ] && perl_opt+=" --no_eq 1$vapor"
    
    cluster_file+="$vapor"
    vapor_suffix+="$vapor"
done

cluster_file+="narrow_neutral"

if [ $l_incl_ions -eq 1 ]; then
    perl_opt+=" --variable_ion_source"
    cluster_file+="_neg_pos"
    vapor_suffix+="_ions"
else
    vapor_suffix+="_noions"
fi

cluster_file+=".inp"


# Generate the equations

perl acdc_2020_04_28.pl --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e ./Perl_input/HS298.15K_example.txt --dip ./Perl_input/dip_pol_298.15K_example.txt `echo "$perl_opt --i ./Perl_input/$cluster_file --append $vapor_suffix$suffix"`
