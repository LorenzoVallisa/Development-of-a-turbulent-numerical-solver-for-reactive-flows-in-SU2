# PACS Project (8 CFU): Development of a turbulent numerical solver for reactive flows

In this project the capability to simulate multispecies turbulent reactive flows has been added to the SU2 Suite. The project infact comprises both a challenging modelling of turbulence phenomenon happening inside a jet reactor, in which chemical reactions are happening due to the combustion process, as well as a remarkable ability to understand and implement complex C++ structures within one of the biggest and most elaborated software written for Computational Fluid Dynamics. Under the physics point of view, the interaction between fuel, oxidizer and combustion products have been modeled using a multispecies approach: one mass conservation equation for each species, two equations for momentum (2D simulation) and one for energy conservation. Turbulent combustion has been modeled using a PaSR approach: indeed the source term of mass conservation equations accounts for the time needed by species to diffuse before reacting, species reacting in a shorter time than the one needed by turbulence to transport them will happen only on a small part of the computational cells, allowing thus the fuel not to fully react at once, but to be able to further expand in the domain. The most challenging part though was the one connected to the implementation of the code and the consequent stability workaround of the numerical scheme. Moreover, the size of SU2 and its polymorphic architecture make it a tough task, even for expert software developer, to add additional working features. After a first phase of understanding how the full system of drivers was interacting with iteration classes, it was in fact clear that a simple “ just add the code for the turbulent closure” was definitely not enough: the implementation of an entire new solver was compulsory if a communication between the reactive laminar solver and the turbulent one (in SU2, as well as in all other CFD software, RANS and mean flow field sets of equations are solved separately and sequentially). The main goal of second phase was therefore to identify all the patterns involved in a SU2 solver processing and find a most efficient way to introduce a new solver without compromising the full structure of the code. Once successfully implemented a reactive-rans solver, the following obstacle to overcome was to retrieve the turbulent methods that a single fluid Navier-Stokes solver class had included thanks to the polymorphic structure, but that were unfortunately hidden by the previous laminar implementation of a multispecies reactive Navier-Stokes solver. PaSR constants (the portion of the cell devoted to combustion) were introduced to account for turbulent combustion, and forsee a lower bound (which can be set within .cfg file), avoiding thus to assume too low and hence not physical values for combustion process to happen.

**Installation** \
\
In order to execute this version of the code it is required a C++ compiler with c++11 standard.
Please check that in the m4 folder the parmetis.m4, metis.m4, compiler.m4, codi.m4, cgns.m4, ax_tls.m4 files are present.
For installation please refer to the following stages:
  1. cd /path/to/SU2
  2. ./bootstrap
  3. ./configure --prefix=/path/to/install/SU2 CXXFLAGS="-O3" LIBS="-lstdc++fs" --enable-mpi --with-cxx=/path/to/mpicxx --with-cc=/path/to/mpicc
  4. make -j 8 install

 As it can be easily noticed it is required to know the location of the MPI compiler; in case it is not available also a simple build is possible substituting step 3 and 4 as it follows:\
  3. ./configure --prefix=/path/to/install/SU2 CXXFLAGS="-O3" LIBS="-lstdc++fs"\
  4. make\
  5. make install

The --prefix option defines the location that the executables will be installed (in a folder named bin/ within your chosen install location from –prefix) and so it requires user access. If the --prefix option is not specified, the code will be installed in /usr/local/bin, which may require admin access.

Make sure to note the SU2_RUN and SU2_HOME environment variables displayed at the conclusion of configure. It is recommended that you add the SU2_RUN and SU2_HOME variables to your ~/.bashrc file and update your PATH variable to include the install location ($SU2_RUN, specified by --prefix).

**Execution** \
\
In order to run the software is necessary a file in cfg format with options for a particular problem specified.
The syntax to run a simulation is

SU2_CFD your_config_file.cfg

or in case MPI is available

mpirun -n 8 SU2_CFD your_config_file.cfg

In this case the software does not return a file to be postprocessed but it is required to launch

mpirun -n 8 SU2_SOL your_config_file.cfg

in order to obtain it.

Please refer to the .sh file you can find within the Test Case folder for examples.

SU2 is capable of outputting solution files that can be visualized in a number of formats, including ParaView (.vtk) and Tecplot (.dat for ASCII, .plt for binary).
For a PDE analysis (direct solution), these files might look like the following:
  flow.dat or flow.vtk: full volume flow solution.

**IMPORTANT REMARK**\
\
In case you want to modify input parameters inside the .cfg file be careful NOT to set Mach number to values too low, moreover if using MPI
the value set in the .cfg will appear only on one node. We tried to solve that issue but due to a lack of time we were not able to fix it.
Sorry.
