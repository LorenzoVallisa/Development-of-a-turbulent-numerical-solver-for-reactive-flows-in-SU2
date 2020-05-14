#!/bin/bash
#PBS -S /bin/bash
#
# This is a bash script, thus every row is a bash command
# unless the first character of the row is a "#" :
# in this case the row is a bash comment.

mpirun -n 4 /vagrant/SU2_Combustion_Turbulent/SU2_CFD/bin/SU2_CFD my_combustion_no_chem.cfg #> output1.txt 2>&1
mpirun -n 4 /vagrant/SU2_Combustion_Turbulent/SU2_CFD/bin/SU2_CFD my_combustion_first_chem_PaSR.cfg #> output2.txt 2>&1
mpirun -n 4 /vagrant/SU2_Combustion_Turbulent/SU2_CFD/bin/SU2_CFD my_combustion_second_chem_PaSR.cfg #> output3.txt 2>&1
mpirun -n 4 /vagrant/SU2_Combustion_Turbulent/SU2_SOL/bin/SU2_SOL my_combustion_second_chem_PaSR.cfg #> output_file_4.txt 2>&1
#---------------------------------------------------------------------#
