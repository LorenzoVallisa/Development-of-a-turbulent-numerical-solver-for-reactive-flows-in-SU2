#!/bin/bash
#PBS -S /bin/bash
#
# This is a bash script, thus every row is a bash command
# unless the first character of the row is a "#" :
# in this case the row is a bash comment.

mpirun -n 8 /usr/local_fin/bin/SU2_CFD my_combustion_no_chem.cfg > output_file.txt 2>&1
mpirun -n 8 /usr/local_fin/bin/SU2_CFD my_combustion_first_chem.cfg > output_file_1.txt 2>&1
mpirun -n 8 /usr/local_fin/bin/SU2_CFD my_combustion_second_chem.cfg > output_file_2.txt 2>&1
mpirun -n 8 /usr/local_fin/bin/SU2_SOL my_combustion_second_chem.cfg
#---------------------------------------------------------------------#
