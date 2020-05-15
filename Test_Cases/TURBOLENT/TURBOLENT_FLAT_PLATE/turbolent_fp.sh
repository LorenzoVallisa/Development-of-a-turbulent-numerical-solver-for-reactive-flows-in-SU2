#!/bin/bash
#PBS -S /bin/bash
#
# This is a bash script, thus every row is a bash command
# unless the first character of the row is a "#" :
# in this case the row is a bash comment.

#mpirun -n 4 /vagrant/SU2_Turbulent_Combustion/SU2_CFD/bin/SU2_CFD my_turbulent_flatplate_air.cfg #> output_1.txt 2>&1
mpirun -n 4 /vagrant/SU2_Turbulent_Combustion/SU2_SOL/bin/SU2_SOL my_turbulent_flatplate_air.cfg #> output_2.txt 2>&1
#---------------------------------------------------------------------#
