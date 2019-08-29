# PACS Project (8 CFU): Development of a laminar numerical solver for reactive flows

In this project the capability to simulate multispecies flows (both reactive and not) has been added to the SU2 Suite; in order to reach this goal these main files that rely on physics have been added: numerics_reactive.hpp, solver_reactive.hpp and variable_reactive.hpp (in the include folder of SU2_CFD directory) and numerics_direct_reactive.cpp, solver_direct_reactive.cpp and variable_direct_reactive.cpp (in the src folder of SU2_CFD directory).
In the solver files we mainly adapted the routines that execute loop for computing fluxes, imposing boundary conditions, obtaining gradients through Green-Gauss or Least-Sqaures and solving linear system.
In the variable files we added two classes (CReactiveEulerVariable and CReactiveNSVariable) for storing the state in case of Euler or Navier-Stokes simulations respectively, while in the numerics files we added classes for computing convective, diffusive and source residuals.

Moreover the original files config_structure, driver_structure, integration_time, integration_structure and output_structure have been modified accordingly in order to make the software interact correctly with this new kind of problem and all additions are segnalated by a comment with the word NOTE.
The file option_structure.hpp has been enriched with the class COptionInlet_MassFrac in order to read mass fractions at inlet which have the difficulty that they are known only at run-time besides some other flags to indicate that we are solving a multispecies problem: even in this case these new variables are denoted by a comment with NOTE.  

Eventually in the directory Common we added two subfolders for the include part and two subfolder for the src part: in Framework we implemented the library that computes the physical and chemical properties and the factory through we call it besides some useful exceptions, while in Tools there is an implementation of cubic spline for computing transport and thermodynamic properties and the routines needed to read the reactions from a textfile.

**Installation** \
\
In order to execute this version of the code it is required a C++ compiler with c++11 standard. For installation please refer to the following stages:
  1. cd /path/to/SU2
  2. ./bootstrap
  3. ./configure --prefix=/path/to/install/SU2 CXXFLAGS="-O3" LIBS="-lstdc++fs" --enable-mpi --with-cxx=/path/to/mpicxx
  4. make -j 8 install

 As it can be easily noticed it is required to know the location of the MPI compiler; in case it is not available also a simple build is possible substituting step 3 and 4 as it follows:
  3. ./configure --prefix=/path/to/install/SU2 CXXFLAGS="-O3" LIBS="-lstdc++fs"
  4. make
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

SU2 is capable of outputting solution files that can be visualized in a number of formats, including ParaView (.vtk) and Tecplot (.dat for ASCII, .plt for binary).
For a PDE analysis (direct solution), these files might look like the following:
  flow.dat or flow.vtk: full volume flow solution.
