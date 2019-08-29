#include "../include/solver_reactive.hpp"
#include "../include/numerics_reactive.hpp"
#include "../../externals/Eigen/IterativeLinearSolvers"

#include "../../Common/include/Framework/factory.hpp"
#include "../../Common/include/Framework/not_implemented_exception.hpp"

#include <algorithm>
#include <limits>
#include <cmath>

namespace {
  /*!
   * \brief Compute area for the face specified by outward normal
   */
  su2double ComputeArea(su2double* Normal, const unsigned short nDim) {
    su2double Area = std::inner_product(Normal, Normal + nDim, Normal, 0.0);
    return std::sqrt(Area);
  }
} /*-- End of unnamed namespace ---*/

/*--- Set the pointer to the library to compute physical-chemical properties to null.
      NOTE: It will be set correctly in the constructor ---*/
CReactiveEulerSolver::LibraryPtr CReactiveEulerSolver::library = NULL;

//
//
/*--- Default constructor ---*/
//
//
CReactiveEulerSolver::CReactiveEulerSolver(): CSolver(), nSpecies(), nPrimVarLim(), space_centered(), implicit(), grid_movement(),
                                              least_squares(), second_order(), limiter(), Density_Inf(), Pressure_Inf(), Temperature_Inf(),
                                              T_INDEX_PRIM(), VX_INDEX_PRIM(), P_INDEX_PRIM(), RHO_INDEX_PRIM(), H_INDEX_PRIM(),
                                              A_INDEX_PRIM(), RHOS_INDEX_PRIM(), RHO_INDEX_SOL(), RHOVX_INDEX_SOL(), RHOE_INDEX_SOL(),
                                              RHOS_INDEX_SOL(), T_INDEX_GRAD(), VX_INDEX_GRAD(), P_INDEX_GRAD(),
                                              T_INDEX_LIM(), VX_INDEX_LIM(), P_INDEX_LIM() {
  /*--- Built-in types variables which is better to set since they are used frequently ---*/
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  IterLinSolver = 0;

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint = 0;
  nPointDomain = 0;
  nVar = 0;
  nPrimVar = 0;
  nPrimVarGrad = 0;
}

//
//
/*--- Class constructor ---*/
//
//
CReactiveEulerSolver::CReactiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh): CSolver() {
  /*--- Local variables ---*/
  unsigned long iPoint;
  unsigned short iVar, iDim;

  /*--- Built-in types variables which is better to set since they are used frequently ---*/
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  IterLinSolver = 0;

  /*--- Load library ---*/
  if(iMesh == MESH_0) {
    Common::Factory<Framework::PhysicalChemicalLibrary> fact(config->GetLibraryName(), config->GetConfigLibFile(),
                                                             config->GetLibraryPath());
    library = LibraryPtr(fact.GetLibraryPtr());
    library->Setup();
  }
  nSpecies = library->GetnSpecies();

  /*--- Set other variables ---*/
  bool restart = (config->GetRestart() || config->GetRestart_Flow());

  /*--- Read and store the restart metadata. ---*/
  if(!(!restart || iMesh != MESH_0))
    Read_SU2_Restart_Metadata(geometry, config);

  /*--- Continue setting other variables ---*/
  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho, rho*vx, rho*vy, rho*vz, rho*E, rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T, vx, vy, vz, P, rho, h, a, Y1,...YNs)^T ---*/
  nPrimVarLim = nDim + 2; /*--- Gradient Primitive variables (T, vx, vy, vz, P)^T ---*/
  nPrimVarGrad = nPrimVarLim; /*--- Primitive variables to limit (T, vx, vy, vz, P)^T ---*/

  /*--- Read freestream mass fractions. ---*/
  MassFrac_Inf = RealVec(config->GetMassFrac_FreeStream(), config->GetMassFrac_FreeStream() + nSpecies);

  /*--- Check right order of species in configuration file ----*/
  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
  if(iMesh == MESH_0 && rank == MASTER_NODE)
    Check_FreeStream_Species_Order(config);

  /*--- Perform the non-dimensionalization for the flow equations using the specified reference values. ---*/
  US_System = (config->GetSystemMeasurements() == US);
  SetNondimensionalization(geometry, config, iMesh);

	/*--- Allocate a CVariable array foreach node of the mesh ---*/
	node = new CVariable*[nPoint];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new su2double[nVar];
  Residual_RMS = new su2double[nVar];
  Residual_Max = new su2double[nVar];
  Residual_i = new su2double[nVar];
  Residual_j = new su2double[nVar];
  Res_Conv = new su2double[nVar];
  Res_Sour = new su2double[nVar];
  Vector_i = new su2double[nDim];
  Vector_j = new su2double[nDim];
  Vector   = new su2double[nDim];

  /*--- Set residual vectors to zero for safety ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    Residual[iVar] = 0.0;
    Residual_RMS[iVar] = 0.0;
    Residual_Max[iVar] = 0.0;
    Residual_i[iVar] = 0.0;
    Residual_j[iVar] = 0.0;
    Res_Conv[iVar] = 0.0;
    Res_Sour[iVar] = 0.0;
  }

  /*--- Define some structures for locating max residuals ---*/
  Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
  for(iVar = 0; iVar < nVar; ++iVar) {
    Point_Max[iVar] = 0;
    Point_Max_Coord[iVar] = new su2double[nDim];
    std::fill(Point_Max_Coord[iVar], Point_Max_Coord[iVar] + nDim, 0.0);
  }

  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];
  Solution_i = new su2double[nVar];
  Solution_j = new su2double[nVar];

  for(iVar = 0; iVar < nVar; ++iVar) {
    Solution[iVar] = 0.0;
    Solution_i[iVar] = 0.0;
    Solution_j[iVar] = 0.0;
  }

	/*--- Allocate vectors related to the primitive ---*/
  PrimVar_i.resize(nPrimVarGrad);
  PrimVar_j.resize(nPrimVarGrad);
  PrimVar_Vertex.resize(nPrimVarGrad);

  Prim_i.resize(nPrimVarLim);
  Prim_j.resize(nPrimVarLim);
  Primitive.resize(nPrimVarLim);

  /*--- Allocate auxiliary vectors ---*/
  Ys_i.resize(nSpecies);
  Ys_j.resize(nSpecies);
  Ys.resize(nSpecies);

  /*--- Initialize the solution & residual CVectors ---*/
 	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
 	if(config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

	/*--- Allocate Jacobians for implicit time-stepping ---*/
	if(config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    implicit = true;
		Jacobian_i = new su2double* [nVar];
		Jacobian_j = new su2double* [nVar];
		for(iVar = 0; iVar < nVar; ++iVar) {
			Jacobian_i[iVar] = new su2double[nVar];
			Jacobian_j[iVar] = new su2double[nVar];
		}

    if(rank == MASTER_NODE)
      std::cout << "Initialize Jacobian structure (Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  }
  else
    implicit = false;

  /*--- Allocate matrices for gradient computation by least squares ---*/
	if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    least_squares = true;

    Smatrix = new su2double* [nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      Smatrix[iDim] = new su2double[nDim];

    Cvector = new su2double* [nPrimVarGrad];
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      Cvector[iVar] = new su2double[nDim];
  }
  else
    least_squares = false;

  /*--- Prepare the arrays to store the value of the characteristic primitive variables at the boundaries ---*/
  nMarker = config->GetnMarker_All();
  CharacPrimVar = new su2double** [nMarker];
  unsigned short iMarker;
  unsigned long iVertex;
  for(iMarker = 0; iMarker < nMarker; ++iMarker) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double[nPrimVar];
      for(iVar = 0; iVar < nPrimVar; ++iVar) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex.reserve(nMarker);
  for(iMarker = 0; iMarker < nMarker; ++iMarker)
    nVertex.push_back(geometry->nVertex[iMarker]);

  /*--- Set other variables ---*/
  grid_movement = config->GetGrid_Movement();
  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  second_order = (config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  if(second_order) {
    Primitive_i.resize(nPrimVar);
    Primitive_j.resize(nPrimVar);
    if(implicit) {
      Secondary_i.resize(nVar);
      Secondary_j.resize(nVar);
    }
  }
  unsigned long ExtIter = config->GetExtIter();
  limiter = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter());

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(), config->GetVelocity_FreeStreamND() + nDim);

  /*--- Initialize the far-field state ----*/
  node_infty = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, nSpecies,
                                          nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);

  /*--- Set local variables to access indexes in proper way for all vectors and arrays. ---*/
  T_INDEX_PRIM    = CReactiveEulerVariable::GetT_INDEX_PRIM();
  VX_INDEX_PRIM   = CReactiveEulerVariable::GetVX_INDEX_PRIM();
  P_INDEX_PRIM    = CReactiveEulerVariable::GetP_INDEX_PRIM();
  RHO_INDEX_PRIM  = CReactiveEulerVariable::GetRHO_INDEX_PRIM();
  H_INDEX_PRIM    = CReactiveEulerVariable::GetH_INDEX_PRIM();
  A_INDEX_PRIM    = CReactiveEulerVariable::GetA_INDEX_PRIM();
  RHOS_INDEX_PRIM = CReactiveEulerVariable::GetRHOS_INDEX_PRIM();

  RHO_INDEX_SOL   = CReactiveEulerVariable::GetRHO_INDEX_SOL();
  RHOVX_INDEX_SOL = CReactiveEulerVariable::GetRHOVX_INDEX_SOL();
  RHOE_INDEX_SOL  = CReactiveEulerVariable::GetRHOE_INDEX_SOL();
  RHOS_INDEX_SOL  = CReactiveEulerVariable::GetRHOS_INDEX_SOL();

  T_INDEX_GRAD    = CReactiveEulerVariable::GetT_INDEX_GRAD();
  VX_INDEX_GRAD   = CReactiveEulerVariable::GetVX_INDEX_GRAD();
  P_INDEX_GRAD    = CReactiveEulerVariable::GetP_INDEX_GRAD();

  T_INDEX_LIM     = CReactiveEulerVariable::GetT_INDEX_LIM();
  VX_INDEX_LIM    = CReactiveEulerVariable::GetVX_INDEX_LIM();
  P_INDEX_LIM     = CReactiveEulerVariable::GetP_INDEX_LIM();

  if(!restart || iMesh != MESH_0) {
    /*--- Initialize the solution to the far-field state everywhere. ---*/
    for(iPoint = 0; iPoint < nPoint; ++iPoint)
      node[iPoint] = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
                                                nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);
  }
  else
    Load_Restart(geometry, config);

  /*--- Allocate vectors for conserved variable limits ---*/
  Lower_Limit.resize(nVar, 0.0);
  Upper_Limit.resize(nVar, 1.0/EPS);

  std::fill(Lower_Limit.begin() + RHOVX_INDEX_SOL, Lower_Limit.begin() + (RHOVX_INDEX_SOL + nDim), -1.0/EPS);
  Lower_Limit[RHOE_INDEX_SOL] = -1.0/EPS;

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
}

//
//
/*--- Class destructor ---*/
//
//
CReactiveEulerSolver::~CReactiveEulerSolver() {
  /*--- Delete free-stream values node ---*/
  if(node_infty != NULL)
    delete node_infty;

  /*--- Delete Characteristic primitive variables ---*/
  if(CharacPrimVar != NULL) {
    for(unsigned short iMarker = 0; iMarker < nMarker; ++iMarker) {
      for(unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex)
        delete[] CharacPrimVar[iMarker][iVertex];
      delete[] CharacPrimVar[iMarker];
    }
    delete[] CharacPrimVar;
  }
}

//
//
/*--- Check coherence in species order declaration. ---*/
//
//
void CReactiveEulerSolver::Check_FreeStream_Species_Order(CConfig* config) {
  /*--- Check coherence in number of species from mass fractions freestream ---*/
  SU2_Assert(nSpecies == config->GetnSpecies(), "Wrong number of species in the configuration file");

  /*--- Check the coherence between the order declared in configuration file and the one declared in the library in case of multispecies ---*/
  auto Config_Order = config->GetSpeciesOrder();

  /*--- Check coherence between the number of species detected from library and from the declaration order first ---*/
  SU2_Assert(nSpecies == config->GetnSpecies(), "Wrong number of species in the configuration file");
  std::string name_species;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    name_species = Config_Order[iSpecies];
    unsigned short idx_library = library->GetIndexSpecies(name_species);
    SU2_Assert(idx_library == iSpecies, "The species " + name_species + " is declared in different position in the library");
  }
}

//
//
/*--- Check whether there is any non physical point in the initialization. ---*/
//
//
void CReactiveEulerSolver::Check_FreeStream_Solution(CConfig* config) {
  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  bool nonPhys_infty = node_infty->SetPrimVar(config);
  SU2_Assert(nonPhys_infty == false, "Assigned freestream values are not compatible with physics");

  unsigned long counter_local = 0, counter_global;
  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Check current solution ---*/
    bool nonPhys = node[iPoint]->SetPrimVar(config);
    /*--- If current solution is not physical use free-stream values ---*/
    if(nonPhys) {
      su2double rho, hs_Tot;
      su2double sqvel = std::inner_product(Velocity_Inf.cbegin(), Velocity_Inf.cend(), Velocity_Inf.cbegin(), 0.0);

      /*--- Compute density from supplied quantities ---*/
      rho = library->ComputeDensity(Temperature_Inf, Pressure_Inf, MassFrac_Inf)*config->GetGas_Constant_Ref();
      Solution[RHO_INDEX_SOL] = rho;
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Solution[RHOS_INDEX_SOL + iSpecies] = MassFrac_Inf[iSpecies]*rho;

      /*--- Compute momentum from supplied primitive quanitites ---*/
      for(unsigned short iDim = 0; iDim < nDim; ++iDim)
        Solution[RHOVX_INDEX_SOL + iDim] = rho*Velocity_Inf[iDim];

      /*---- Compute energy from temperature ---*/
      su2double dim_temp = Temperature_Inf*config->GetTemperature_Ref();
      if(US_System)
        dim_temp *= 5.0/9.0;
      hs_Tot = library->ComputeEnthalpy(dim_temp, MassFrac_Inf)/config->GetEnergy_Ref();
      if(US_System)
        hs_Tot *= 3.28084*3.28084;
      Solution[RHOE_INDEX_SOL] = rho*(hs_Tot + 0.5*sqvel) - Pressure_Inf;

      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);

      counter_local++;
    }
  } /*--- End of for-loop over iPoint ---*/

  /*--- Warning message about non-physical points ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
    #else
      counter_global = counter_local;
    #endif
    if(rank == MASTER_NODE && counter_global != 0)
      std::cout << "Warning. The original solution contains "<< counter_global << " points that are not physical." << std::endl;
  }
} /*--- End of Check_FreeStream_Solution() function ---*/

//
//
/*--- Read basic infos in case of restart simulation. ---*/
//
//
void CReactiveEulerSolver::Read_SU2_Restart_Metadata(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  std::string filename = config->GetSolution_FlowFileName();
  int Unst_RestartIter;

  /*--- Multizone problems require the number of the zone to be appended. ---*/
  if(nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Modify file name for a dual-time unsteady restart ---*/
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if(dual_time) {
    if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    else
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 2;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Modify file name for a simple unsteady restart ---*/
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  if(time_stepping) {
    Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Open the restart file, throw an error if this fails. ---*/
  int rank = MASTER_NODE;
  std::ifstream restart_file;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
  restart_file.open(filename.data(), std::ios::in);
  if(restart_file.fail()) {
    if(rank == MASTER_NODE)
      std::cout<<"There is no flow restart file!! "<<filename.data()<<"."<< endl;
    std::exit(EXIT_FAILURE);
  }

  unsigned long iPoint_Global;
  std::string text_line;

  /*--- The first line is the header (General description) ---*/
  std::getline(restart_file, text_line);

  /*--- Space for the solution ---*/
  for(iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global)
    std::getline(restart_file, text_line);

  /*--- Space for extra info (if any) ---*/
  std::string::size_type position;
  while(std::getline(restart_file, text_line)) {
    /*--- Angle of attack ---*/
    position = text_line.find ("AOA=",0);
    if(position != std::string::npos) {
      text_line.erase(0,4);
      su2double AoA_ = std::stod(text_line);
      if(config->GetDiscard_InFiles() == false) {
        if(config->GetAoA() != AoA_ &&  rank == MASTER_NODE) {
          std::cout.precision(6);
          std::cout << std::fixed <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << std::endl;
          std::cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " <<
                       AoA_ + config->GetAoA_Offset() << " deg." << std::endl;
        }
        config->SetAoA(AoA_ + config->GetAoA_Offset());
      }
      else {
        if(config->GetAoA() != AoA_ &&  rank == MASTER_NODE)
          std::cout <<"WARNING: Discarding the AoA in the solution file." << std::endl;
      }
    } /*--- End AoA ---*/

    /*--- Sideslip angle ---*/
    position = text_line.find ("SIDESLIP_ANGLE=",0);
    if(position != std::string::npos) {
      text_line.erase(0,15);
      su2double AoS_ = std::stod(text_line);
      if(config->GetDiscard_InFiles() == false) {
        if(config->GetAoS() != AoS_ &&  rank == MASTER_NODE) {
          std::cout.precision(6);
          std::cout << std::fixed <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << std::endl;
          std::cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " <<
                       AoS_ + config->GetAoS_Offset() << " deg." << std::endl;
        }
        config->SetAoS(AoS_ + config->GetAoS_Offset());
      }
      else {
        if(config->GetAoS() != AoS_ &&  rank == MASTER_NODE)
          std::cout <<"WARNING: Discarding the AoS in the solution file." << std::endl;
      }
    } /*--- End sideslip angle ---*/

    /*--- BCThrust angle ---*/
    position = text_line.find("INITIAL_BCTHRUST=",0);
    if(position != std::string::npos) {
      text_line.erase(0,17);
      su2double BCThrust_ = std::stod(text_line);
      if(config->GetDiscard_InFiles() == false) {
        if(config->GetInitial_BCThrust() != BCThrust_ &&  rank == MASTER_NODE)
          std::cout <<"WARNING: ACDC will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << std::endl;
        config->SetInitial_BCThrust(BCThrust_);
      }
      else {
        if(config->GetInitial_BCThrust() != BCThrust_ &&  rank == MASTER_NODE)
          std::cout <<"WARNING: Discarding the BC Thrust in the solution file." << std::endl;
      }
    } /*--- End of BCThrust angle ---*/

    /*--- External iteration ---*/
    position = text_line.find("EXT_ITER=",0);
    if(position != std::string::npos) {
      text_line.erase(0,9);
      unsigned long ExtIter_ = std::stoi(text_line);
      if(!config->GetContinuous_Adjoint() && !config->GetDiscrete_Adjoint())
        config->SetExtIter_OffSet(ExtIter_);
    }

  } /*--- End of extra infos reading ---*/

  /*--- Close the restart file... NOTE: we will open this file again... ---*/
  restart_file.close();
}

//
//
/*--- Read the file to restart. ---*/
//
//
void CReactiveEulerSolver::Load_Restart(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  std::string filename = config->GetSolution_FlowFileName();
  int Unst_RestartIter;

  /*--- Multizone problems require the number of the zone to be appended. ---*/
  if(nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Modify file name for a dual-time unsteady restart ---*/
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if(dual_time) {
    if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    else
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 2;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Modify file name for a simple unsteady restart ---*/
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  if (time_stepping) {
    Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Open the restart file, throw an error if this fails. ---*/
  int rank = MASTER_NODE;
  std::ifstream restart_file;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
  restart_file.open(filename.data(), std::ios::in);
  if(restart_file.fail()) {
    if(rank == MASTER_NODE)
      std::cout<<"There is no flow restart file!! "<<filename.data()<<"."<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*--- In case this is a parallel simulation, we need to perform the
        Global2Local index transformation first. ---*/
  std::map<unsigned long, unsigned long> Global2Local;
  std::map<unsigned long, unsigned long>::const_iterator MI;

  /*--- Now fill array with the transform values only for local points ---*/
  unsigned long iPoint;
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

  /*--- Read all lines in the restart file ---*/
  unsigned long iPoint_Local, iPoint_Global;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  std::string text_line;

  /*--- The first line is the header ---*/
  std::getline(restart_file, text_line);
  unsigned long index;
  su2double dull_val;
  for(iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global) {
    std::getline(restart_file, text_line);
    std::istringstream point_line(text_line);
    if(iPoint_Global >= geometry->GetGlobal_nPointDomain()) {
      sbuf_NotMatching = 1;
      break;
    }

    /*--- Retrieve local index. If this node from the restart file lives
          on the current processor, we will load and instantiate the variables. ---*/
    MI = Global2Local.find(iPoint_Global);
    if(MI != Global2Local.end()) {
      iPoint_Local = Global2Local[iPoint_Global];
      if(nDim == 2)
        point_line >> index >> dull_val >> dull_val >> Solution[RHO_INDEX_SOL] >> Solution[RHOVX_INDEX_SOL] >>
                      Solution[RHOVX_INDEX_SOL + 1] >> Solution[RHOE_INDEX_SOL];
      else if(nDim == 3)
        point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[RHO_INDEX_SOL] >> Solution[RHOVX_INDEX_SOL] >>
                      Solution[RHOVX_INDEX_SOL + 1] >> Solution[RHOVX_INDEX_SOL + 2] >> Solution[RHOE_INDEX_SOL];
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        point_line >> Solution[RHOS_INDEX_SOL + iSpecies];

      node[iPoint_Local] = new CReactiveEulerVariable(Solution, nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);
      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/
  if(iPoint_Global_Local < nPointDomain)
    sbuf_NotMatching = 1;

  #ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
  #else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if(rbuf_NotMatching != 0) {
    if(rank == MASTER_NODE) {
      std::cout<<std::endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }

    #ifndef HAVE_MPI
      std::exit(EXIT_FAILURE);
    #else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    #endif
  }

  /*--- Instantiate the variable class with an arbitrary solution
        at any halo/periodic nodes. The initial solution can be arbitrary,
        because a send/recv is performed immediately in the solver. ---*/
  for(iPoint = nPointDomain; iPoint < nPoint; ++iPoint)
    node[iPoint] = new CReactiveEulerVariable(Solution, nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);

  /*--- Close the restart file ---*/
  restart_file.close();
}

//
//
/*--- Set the initial solution in case of restart ---*/
//
//
void CReactiveEulerSolver::SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig *config,
                                               unsigned long ExtIter) {
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, Solution[nVar];

  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- If restart solution, then interpolate the flow solution to
        all the multigrid levels, this is important with the dual time strategy ---*/
  if(restart && ExtIter == 0) {
    for(iMesh = 1; iMesh <= config->GetnMGLevels(); ++iMesh) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); ++iPoint) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for(iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); ++iChildren) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          auto Solution_Fine = solver_container[iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution();
          for(iVar = 0; iVar < nVar; ++iVar) {
            Solution[iVar] = Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->SetSolution(Solution);
      }
      solver_container[iMesh][FLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
  } /*--- End of restart ---*/

  /*--- The value of the solution for the first iteration of the dual time ---*/
  if(dual_time && (ExtIter == 0 || (restart && static_cast<long>(ExtIter) == config->GetUnst_RestartIter()))) {
    /*--- Push back the initial condition to previous solution containers
          for a 1st-order restart or when simply intitializing to freestream. ---*/
    for(iMesh = 0; iMesh <= config->GetnMGLevels(); ++iMesh) {
      for(iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); ++iPoint) {
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }

    if(restart && static_cast<long>(ExtIter) == config->GetUnst_RestartIter() && config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {
      /*--- Load an additional restart file for a 2nd-order restart ---*/
      solver_container[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver_container, config, SU2_TYPE::Int(config->GetUnst_RestartIter() - 1));

      /*--- Push back this new solution to time level N. ---*/
      for(iMesh = 0; iMesh <= config->GetnMGLevels(); ++iMesh) {
        for(iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); ++iPoint)
          solver_container[iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
      }
    }
  }
}

//
//
/*--- Set the fluid solver nondimensionalization. ---*/
void CReactiveEulerSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {
   su2double Temperature_FreeStream = 0.0, ModVel_FreeStream = 0.0, Energy_FreeStream = 0.0,
             Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, SoundSpeed_FreeStream = 0.0;

   su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
             Temperature_Ref = 0.0, Time_Ref = 0.0, Gas_Constant_Ref = 0.0, Energy_Ref = 0.0;

   su2double Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0, ModVel_FreeStreamND = 0.0,
             Temperature_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0, Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

   su2double Velocity_FreeStreamND[nDim];

   unsigned short iDim;

   /*--- Check if the simulation is unsteady ---*/
   bool unsteady = config->GetUnsteady_Simulation();

   int rank = MASTER_NODE;
   #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   #endif

 	 /*--- Compressible non dimensionalization ---*/
   Pressure_FreeStream = config->GetPressure_FreeStream();
   Temperature_FreeStream = config->GetTemperature_FreeStream();
   Density_FreeStream = library->ComputeDensity(Temperature_FreeStream, Pressure_FreeStream, MassFrac_Inf);
   config->SetDensity_FreeStream(Density_FreeStream);
   SoundSpeed_FreeStream = library->ComputeFrozenSoundSpeed(Temperature_FreeStream, MassFrac_Inf);

   /*--- Compute the modulus of the free stream velocity ---*/
   ModVel_FreeStream = std::inner_product(config->GetVelocity_FreeStream(), config->GetVelocity_FreeStream() + nDim,
                                          config->GetVelocity_FreeStream(), 0.0);
   ModVel_FreeStream = std::sqrt(ModVel_FreeStream);
   config->SetModVel_FreeStream(ModVel_FreeStream);

   /*--- Compute the free stream energy ---*/
   Energy_FreeStream = library->ComputeEnergy(Temperature_FreeStream, MassFrac_Inf) + 0.5*ModVel_FreeStream*ModVel_FreeStream;

   /*--- Compute non dimensional quantities: Notice that the grid is in meters. ---*/
   if(config->GetRef_NonDim() == DIMENSIONAL) {
     Pressure_Ref = 1.0;
     Density_Ref  = 1.0;
     Temperature_Ref = 1.0;
     Length_Ref = 1.0;
   }
   else if(config->GetRef_NonDim() == REFERENCE){
     Pressure_Ref     = config->GetPressure_Ref();
     Density_Ref      = config->GetDensity_Ref();
     Temperature_Ref  = config->GetTemperature_Ref();
     Length_Ref       = config->GetLength_Ref();
   }
   else
     throw std::out_of_range("Invalid option for adimensionalitazion for multispecies flows.");

   /*--- Set the computed reference quantitites ---*/
   config->SetPressure_Ref(Pressure_Ref);
   config->SetDensity_Ref(Density_Ref);
   config->SetTemperature_Ref(Temperature_Ref);
   config->SetLength_Ref(Length_Ref);

   /*--- Compute the reference quantities from the supplied ones for all other variables ---*/
   Velocity_Ref = std::sqrt(Pressure_Ref/Density_Ref);
   config->SetVelocity_Ref(Velocity_Ref);

   Time_Ref = Length_Ref/Velocity_Ref;
   config->SetTime_Ref(Time_Ref);

   Energy_Ref = Velocity_Ref*Velocity_Ref;
   config->SetEnergy_Ref(Energy_Ref);

   Gas_Constant_Ref = Energy_Ref/Temperature_Ref;
   config->SetGas_Constant_Ref(Gas_Constant_Ref);

   /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/
   Pressure_FreeStreamND = Pressure_FreeStream/Pressure_Ref;
   config->SetPressure_FreeStreamND(Pressure_FreeStreamND);

   Density_FreeStreamND  = Density_FreeStream/Density_Ref;
   config->SetDensity_FreeStreamND(Density_FreeStreamND);

   for(iDim = 0; iDim < nDim; ++iDim) {
     Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref;
     config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
   }

   Temperature_FreeStreamND = Temperature_FreeStream/Temperature_Ref;
   config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

   Energy_FreeStreamND = Energy_FreeStream/Energy_Ref;
   config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

   ModVel_FreeStreamND = std::inner_product(Velocity_FreeStreamND, Velocity_FreeStreamND + nDim, Velocity_FreeStreamND, 0.0);
   ModVel_FreeStreamND = std::sqrt(ModVel_FreeStreamND);
   config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

   Total_UnstTimeND = config->GetTotal_UnstTime()/Time_Ref;
   config->SetTotal_UnstTimeND(Total_UnstTimeND);

   Delta_UnstTimeND = config->GetDelta_UnstTime()/Time_Ref;
   config->SetDelta_UnstTimeND(Delta_UnstTimeND);

   /*--- Write output to the console ifthis is the master node and first domain ---*/
   if(config->GetConsole_Output_Verb() == VERB_HIGH && rank == MASTER_NODE && iMesh == MESH_0) {
     std::cout.precision(6);

     /*--- Print out reference values. ---*/
     std::cout <<"-- Reference values:"<< std::endl;

     bool SI_Measurement = (config->GetSystemMeasurements() == SI);
     bool US_Measurament = US_System;

     std::cout << "Reference specific gas constant: " << Gas_Constant_Ref;
     if(SI_Measurement)
       std::cout << " N.m/kg.K." << std::endl;
     else if(US_Measurament)
       std::cout << " lbf.ft/slug.R." << std::endl;

     std::cout << "Reference pressure: " << Pressure_Ref;
     if(SI_Measurement)
       std::cout << " Pa." << std::endl;
     else if(US_Measurament)
       std::cout << " psf." << std::endl;

     std::cout << "Reference temperature: " << Temperature_Ref;
     if(SI_Measurement)
       std::cout << " K." << std::endl;
     else if(US_Measurament)
       std::cout << " R." << std::endl;

     std::cout << "Reference density: " << Density_Ref;
     if(SI_Measurement)
       std::cout << " kg/m^3." << std::endl;
     else if(US_Measurament)
       std::cout << " slug/ft^3." << std::endl;

     std::cout << "Reference velocity: " << Velocity_Ref;
     if(SI_Measurement)
       std::cout << " m/s." << std::endl;
     else if(US_Measurament)
       std::cout << " ft/s." << std::endl;

     std::cout << "Reference energy per unit mass: " << Energy_Ref;
     if(SI_Measurement)
       std::cout << " m^2/s^2." << std::endl;
     else if(US_Measurament)
       std::cout << " ft^2/s^2." << std::endl;

     if(unsteady)
       std::cout << "Reference time: " << Time_Ref <<" s." << std::endl;

     /*--- Print out resulting non-dim values here. ---*/
     std::cout << "-- Resulting non-dimensional state:" << std::endl;

     config->SetMach(ModVel_FreeStream/SoundSpeed_FreeStream);
     std::cout << "Mach number (non-dim): " << config->GetMach() << std::endl;

     std::cout << "Free-stream temperature (non-dim): " << Temperature_FreeStreamND << std::endl;

     std::cout << "Free-stream pressure (non-dim): " << Pressure_FreeStreamND << std::endl;

     std::cout << "Free-stream density (non-dim): " << Density_FreeStreamND << std::endl;

     std::cout << "Free-stream velocity (non-dim): (";
     for(iDim = 0; iDim < nDim; ++iDim)
      iDim < (nDim -1) ? std::cout << Velocity_FreeStreamND[iDim] << ", " :
                         std::cout << Velocity_FreeStreamND[iDim] << "). "<<std::endl;

     std::cout << "Magnitude (non-dim): " << ModVel_FreeStreamND << std::endl;

     std::cout << "Free-stream total energy per unit mass (non-dim): " << Energy_FreeStreamND << std::endl;

     if(unsteady) {
       std::cout << "Total time (non-dim): " << Total_UnstTimeND << std::endl;
       std::cout << "Time step (non-dim): " << Delta_UnstTimeND << std::endl;
     }
     std::cout << std::endl;
   }
 }

//
//
/*--- Set primitive variables according to the current state in the preprocessing. ---*/
//
//
unsigned long CReactiveEulerSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {
  unsigned long iPoint, ErrorCounter = 0;
  bool NonPhys = false;

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Initialize the non-physical points vector ---*/
    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies + nDim + 5, (T, vx, vy, vz, P, rho, h, a, Y1,....YNs) ---*/
    NonPhys = node[iPoint]->SetPrimVar(config);
    bool ignition = config->GetIgnition();
    /*--- Check if ignition is desired ---*/
    if(ignition) {
      if(config->GetExtIter() < config->GetIgnitionIter()) {
        unsigned short Fuel_Index = config->GetFuelIndex();
        unsigned short Oxidizer_Index = config->GetOxidizerIndex();
        if(node[iPoint]->GetMassFraction(Fuel_Index) > 0.4 && node[iPoint]->GetMassFraction(Oxidizer_Index) > 0.2 &&
           node[iPoint]->GetTemperature() < config->GetIgnitionTemperature())
           bool tmp = node[iPoint]->SetTemperature(config->GetIgnitionTemperature());
      }
    }

    if(NonPhys) {
      node[iPoint]->SetNon_Physical(true);
      ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/
    if(!Output)
      LinSysRes.SetBlock_Zero(iPoint);
  }

  return ErrorCounter;
}

//
//
/*--- Preprocessing stage to set primitive variables and eventually gradient of primitive variables. ---*/
//
//
void CReactiveEulerSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                                         unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  /*--- Set the primitive variables ---*/
  unsigned long ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Upwind second order reconstruction ---*/
  if((second_order && !space_centered) && iMesh == MESH_0 && !Output) {
    /*--- Gradient computation ---*/
    if(least_squares)
      SetPrimitive_Gradient_LS(geometry, config);
    else if(config->GetKind_Gradient_Method() == GREEN_GAUSS)
      SetPrimitive_Gradient_GG(geometry, config);

    /*--- Limiter computation ---*/
    if(limiter && iMesh == MESH_0 && !Output)
      SetPrimitive_Limiter(geometry, config);
  }

  /*--- Artificial dissipation ---*/
  if(space_centered && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");

  /*--- Initialize the Jacobian matrices ---*/
  if(implicit)
    Jacobian.SetValZero();

  /*--- Error message ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      SU2_MPI::Allreduce(MPI_IN_PLACE, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    #endif
    if(iMesh == MESH_0)
      config->SetNonphysical_Points(ErrorCounter);
  }
}

//
//
/*--- Set gradient of primitive variables according to Green-Gauss. ---*/
//
//
void CReactiveEulerSolver::SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;

  su2double PrimVar_Average, Partial_Res;
  su2double Volume;
  su2double Normal[nDim];

  /*--- Gradient of the following primitive variables: [T, vx, vy, vz, p]^T ---*/
  /*--- Set Gradient_Primitive to zero ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
  	node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
  	iPoint = geometry->edge[iEdge]->GetNode(0);
  	jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
    PrimVar_i[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_j[T_INDEX_GRAD] = node[jPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_i[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
    PrimVar_j[P_INDEX_GRAD] = node[jPoint]->GetPrimitive(P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim) {
      PrimVar_i[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
      PrimVar_j[VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
    }

  	geometry->edge[iEdge]->GetNormal(Normal);

    /*--- Compute face averaged value as discrete counter-part of Green-Gauss theorem ---*/
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      PrimVar_Average = 0.5*(PrimVar_i[iVar] + PrimVar_j[iVar]);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Partial_Res = PrimVar_Average*Normal[iDim];
        if(geometry->node[iPoint]->GetDomain())
  			  node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if(geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  	  }
    }
  } /*--- End loop interior edges ---*/

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
  	for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
  		iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  		if(geometry->node[iPoint]->GetDomain()) {
        /*--- Get primitives from CVariable. NOTE: On the boundary we cannot average so we use this value directly ---*/
        PrimVar_Vertex[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
        PrimVar_Vertex[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
        for(iDim = 0; iDim < nDim; ++iDim)
          PrimVar_Vertex[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);

        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
        for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
          for(iDim = 0; iDim < nDim; ++iDim) {
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
  				  node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
        	}
        }
  		}
  	} /*--- End of iVertex for loop ---*/
  } /*--- End of iMarker for loop ---*/

  /*--- Update gradient value ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();
    SU2_Assert(Volume > EPS,"The measure of the volume is not consistent(smaller or equal to zero)");
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  	  for(iDim = 0; iDim < nDim; ++iDim)
  		  node[iPoint]->SetGradient_Primitive(iVar, iDim, node[iPoint]->GetGradient_Primitive(iVar,iDim)/Volume);
    }
  }

  Set_MPI_Primitive_Gradient(geometry, config);
}

//
//
/*--- Set gradient primitive variables according to weighted least squares. ---*/
//
//
void CReactiveEulerSolver::SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iDim, jDim, iNeigh;
	unsigned long iPoint, jPoint;

  su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, detR2, z11, z12, z13, z22, z23, z33;
  su2double weight;
  su2double Coord_ij[nDim];

  bool singular;

  /*--- Gradient of the following primitive variables: [T,u,v,w,p]^T ---*/
  /*--- Loop over points of the grid ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		/*--- Set the value of singular ---*/
		singular = false;

    /*--- Get coordinates ---*/
		auto Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get primitives from CVariable ---*/
    PrimVar_i[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_i[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim)
      PrimVar_i[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);

    /*--- Inizialization of variables ---*/
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      std::fill(Cvector[iVar], Cvector[iVar] + nDim, 0.0);

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			auto Coord_j = geometry->node[jPoint]->GetCoord();

      /*--- Get primitives from CVariable ---*/
      PrimVar_j[T_INDEX_GRAD] = node[jPoint]->GetPrimitive(T_INDEX_PRIM);
      PrimVar_j[P_INDEX_GRAD] = node[jPoint]->GetPrimitive(P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        PrimVar_j[VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);

      /*--- Compute components of vector that links node i to node j ---*/
		  std::transform(Coord_j, Coord_j + nDim, Coord_i, Coord_ij, std::minus<su2double>());

      weight = std::inner_product(Coord_ij, Coord_ij + nDim, Coord_ij, 0.0);

			/*--- Sumations for entries of upper triangular matrix R ---*/
      if(weight > EPS) {
        r11 += Coord_ij[0]*Coord_ij[0]/weight;
        r12 += Coord_ij[0]*Coord_ij[1]/weight;
        r22 += Coord_ij[1]*Coord_ij[1]/weight;
        if(nDim == 3) {
          r13 += Coord_ij[0]*Coord_ij[2]/weight;
          r23_a += Coord_ij[1]*Coord_ij[2]/weight;
          r23_b += Coord_ij[0]*Coord_ij[2]/weight;
          r33 += Coord_ij[2]*Coord_ij[2]/weight;
        }
        /*--- Entries of c:= transpose(A)*b ---*/
        for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
          for(iDim = 0; iDim < nDim; ++iDim)
            Cvector[iVar][iDim] += Coord_ij[iDim]*(PrimVar_j[iVar] - PrimVar_i[iVar])/weight;
      }
    } /*--- End of iNeigh for loop---*/

		/*--- Entries of upper triangular matrix R ---*/
    if(r11 > EPS)
      r11 = std::sqrt(r11);
    else
      r11 = 0.0;

    if(std::abs(r11) > EPS)
      r12 = r12/r11;
    else
      r12 = 0.0;

    if(r22 - r12*r12 > EPS)
      r22 = std::sqrt(r22 - r12*r12);
    else
      r22 = 0.0;

    if(nDim == 3) {
      if(std::abs(r11) > EPS)
        r13 = r13/r11;
      else
        r13 = 0.0;

      if(std::abs(r22) > EPS && std::abs(r11*r22) > EPS)
        r23 = r23_a/r22 - r23_b*r12/(r11*r22);
      else
        r23 = 0.0;
      if(r33 - r23*r23 - r13*r13 > EPS)
        r33 = std::sqrt(r33 - r23*r23 - r13*r13);
      else r33 = 0.0;
    }

    /*--- Compute determinant ---*/
    if(nDim == 2)
      detR2 = (r11*r22)*(r11*r22);
    else
      detR2 = (r11*r22*r33)*(r11*r22*r33);

    /*--- Detect singular matrices ---*/
    if(std::abs(detR2) < EPS) {
      detR2 = 1.0;
      singular = true;
    }

		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    if(singular) {
      for(iDim = 0; iDim < nDim; ++iDim)
        std::fill(Smatrix[iDim], Smatrix[iDim] + nDim, 0.0);
    }
    else {
      if(nDim == 2) {
        Smatrix[0][0] = (r12*r12 + r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23 - r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11 + z12*z12 + z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22 + z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22 + z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }

    /*--- Computation of the gradient: C*S ---*/
    su2double result;
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
		  for(iDim = 0; iDim < nDim; ++iDim) {
        result = 0.0;
        for(jDim = 0; jDim < nDim; ++jDim)
          result += Cvector[iVar][jDim]*Smatrix[iDim][jDim];
        node[iPoint]->SetGradient_Primitive(iVar, iDim, result);
      }
    }

	} /*--- End of iPoint for loop ---*/

	Set_MPI_Primitive_Gradient(geometry, config);
}


//
//
/*--- Set limiter of primitive variables ---*/
//
//
void CReactiveEulerSolver::SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;

  su2double dave, LimK, eps1, eps2, dm, dp, du, y, limiter_val;

  dave = config->GetRefElemLength();
  LimK = config->GetLimiterCoeff();

  if(!limiter) {
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        node[iPoint]->SetLimiter_Primitive(iVar, 1.0);
      }
    }
  }
  else {
    /*--- Initialize solution max, solution min and limiter in entire domain ---*/
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
      }
    }

    /*--- Establish bounts for Spekreijse monotonicity by finding max/min values of neighbor variables ---*/
    for(iEdge = 0; iEdge < geometry-> GetnEdge(); ++iEdge) {
      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry-> edge[iEdge]->GetNode(0);
      jPoint = geometry-> edge[iEdge]->GetNode(1);

      /*--- Get primitive variables ---*/
      Prim_i[T_INDEX_LIM] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
      Prim_i[P_INDEX_LIM] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
      Prim_j[T_INDEX_LIM] = node[jPoint]->GetPrimitive(T_INDEX_PRIM);
      Prim_j[P_INDEX_LIM] = node[jPoint]->GetPrimitive(P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Prim_i[VX_INDEX_LIM + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
        Prim_j[VX_INDEX_LIM + iDim] = node[jPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
      }

      /*--- Compute the max and min values for nodes i and j ---*/
      for(iVar = 0; iVar < nPrimVarLim; ++iVar){
        du = (Prim_j[iVar] - Prim_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, std::min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, std::max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, std::min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, std::max(node[jPoint]->GetSolution_Max(iVar), -du));
      }
    }
  } /*-- End of prepprocessing phase for limiting procedure ---*/

  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  if(config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      auto Coord_i = geometry->node[iPoint]->GetCoord();
      auto Coord_j = geometry->node[jPoint]->GetCoord();

      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        /*--- Compute the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_j[iDim] - Coord_i[iDim])*Gradient_i[iVar][iDim];

        if(dm < EPS)
          limiter_val = 2.0;
        else {
          if(dm > EPS)
            dp = node[iPoint]->GetSolution_Max(iVar);
          else
            dp = node[iPoint]->GetSolution_Min(iVar);
          limiter_val = dp/dm;
        }

        if(limiter_val < node[iPoint]->GetLimiter_Primitive(iVar))
          node[iPoint]->SetLimiter_Primitive(iVar, limiter_val);

        /*--- Compute the interface right gradient, delta+ (dp) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if(dm < EPS)
          limiter_val = 2.0;
        else {
          if(dm > EPS)
            dp = node[jPoint]->GetSolution_Max(iVar);
          else
            dp = node[jPoint]->GetSolution_Min(iVar);
          limiter_val = dp/dm;
        }

        if(limiter < node[jPoint]->GetLimiter_Primitive(iVar))
          node[jPoint]->SetLimiter_Primitive(iVar, limiter_val);

      } /*--- End of limiter computation for a single edge ---*/

    } /*--- End of iEdge for loop ---*/

    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter_val = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter_val);
      }
    }
  } /*--- End of Barth-Jespersen limiter ---*/

  /*--- Venkatakrishnan limiter ---*/
  if(config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) {
    /*--- Allocate memory for the max and min primitive value --*/
    su2double  LocalMinPrimitive[nPrimVarLim], GlobalMinPrimitive[nPrimVarLim],
               LocalMaxPrimitive[nPrimVarLim], GlobalMaxPrimitive[nPrimVarLim];

    /*--- Compute the max value and min value of the solution ---*/
    std::fill(LocalMinPrimitive, LocalMinPrimitive + nPrimVarLim, std::numeric_limits<su2double>::max());
    std::fill(LocalMaxPrimitive, LocalMaxPrimitive + nPrimVarLim, std::numeric_limits<su2double>::min());

    /*--- Loop over the internal domain ---*/
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      /*--- Get the involved primitive variables ---*/
      Primitive[T_INDEX_LIM] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
      Primitive[P_INDEX_LIM] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        Primitive[VX_INDEX_LIM + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);

      std::transform(Primitive.cbegin(), Primitive.cend(), LocalMinPrimitive, LocalMinPrimitive, std::less<su2double>());
      std::transform(Primitive.cbegin(), Primitive.cend(), LocalMaxPrimitive, LocalMaxPrimitive, std::greater<su2double>());
    }

    #ifdef HAVE_MPI
      SU2_MPI::Allreduce(LocalMinPrimitive, GlobalMinPrimitive, nPrimVarLim, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(LocalMaxPrimitive, GlobalMaxPrimitive, nPrimVarLim, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #else
      std::copy(LocalMinPrimitive, LocalMinPrimitive + nPrimVarLim, GlobalMinPrimitive);
      std::copy(LocalMaxPrimitive, LocalMaxPrimitive + nPrimVarLim, GlobalMaxPrimitive);
    #endif

    /*--- Loop over the interior edges ---*/
    for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      auto Coord_i = geometry->node[iPoint]->GetCoord();
      auto Coord_j = geometry->node[jPoint]->GetCoord();

      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        eps1 = LimK*dave;
        eps2 = eps1*eps1*eps1;

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        if(dm > 0.0)
          dp = node[iPoint]->GetSolution_Max(iVar);
        else
          dp = node[iPoint]->GetSolution_Min(iVar);

        limiter_val = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if(limiter_val < node[iPoint]->GetLimiter_Primitive(iVar))
          node[iPoint]->SetLimiter_Primitive(iVar, limiter_val);

        /*-- Repeat for point j on the edge ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if(dm > 0.0)
          dp = node[jPoint]->GetSolution_Max(iVar);
        else
          dp = node[jPoint]->GetSolution_Min(iVar);

        limiter_val = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if(limiter_val < node[jPoint]->GetLimiter_Primitive(iVar))
          node[jPoint]->SetLimiter_Primitive(iVar, limiter_val);

      } /*--- End of iVar for loop ---*/

    }
  } /*--- End of Venkatakrishnan limiter ---*/

  /*--- Limiter MPI ---*/
  Set_MPI_Primitive_Limiter(geometry, config);
}

//
//
/*--- Call MPI to set solution for parallel simulations. ---*/
//
//
void CReactiveEulerSolver::Set_MPI_Solution(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  #ifdef HAVE_MPI
	 int send_to, receive_from;
   MPI_Status status;
  #endif

	for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
		if(config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE && config->GetMarker_All_SendRecv(iMarker) > 0) {
			MarkerS = iMarker;
      MarkerR = iMarker + 1;

      #ifdef HAVE_MPI
        send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
			  receive_from = std::abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      #endif

		  nVertexS = geometry->nVertex[MarkerS];
      nBufferS_Vector = nVertexS*nVar;
      nVertexR = geometry->nVertex[MarkerR];

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Send_U.resize(nBufferS_Vector);

      /*--- Copy the solution that should be sended ---*/
      for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for(iVar = 0; iVar < nVar; ++iVar)
          Buffer_Send_U[iVar*nVertexS + iVertex] = node[iPoint]->GetSolution(iVar);
      }

      #ifdef HAVE_MPI
        unsigned long nBufferR_Vector = nVertexR*nVar;
        Buffer_Receive_U.resize(nBufferR_Vector);
		    /*--- Send/Receive information using Sendrecv ---*/
        SU2_MPI::Sendrecv(Buffer_Send_U.data(), nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                          Buffer_Receive_U.data(), nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      #else
        Buffer_Receive_U = Buffer_Send_U;
      #endif

      /*--- Do the coordinate transformation ---*/
      for(iVertex = 0; iVertex < nVertexR; ++iVertex) {
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        auto angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];        phi    = angles[1];          psi    = angles[2];
        cosTheta = std::cos(theta);  cosPhi = std::cos(phi);      cosPsi = std::cos(psi);
        sinTheta = std::sin(theta);  sinPhi = std::sin(phi);      sinPsi = std::sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
        ordering is rotation about the x-axis, y-axis,then z-axis. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;

        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;

        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for(iVar = 0; iVar < nVar; ++iVar)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if(nDim == 2) {
          Solution[RHOVX_INDEX_SOL] = rotMatrix[0][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                      rotMatrix[0][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[1][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL+1)*nVertexR + iVertex];
        }
        else {
          Solution[RHOVX_INDEX_SOL] = rotMatrix[0][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                      rotMatrix[0][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                      rotMatrix[0][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[1][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                          rotMatrix[1][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 2] = rotMatrix[2][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[2][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                          rotMatrix[2][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        node[iPoint]->SetSolution(Solution);
      } /*--- End of iVertex for loop ---*/
    }
	} /*--- End of iMarker for loop ---*/
}

//
//
/*--- Call MPI to set solution for parallel simulations. ---*/
//
//
void CReactiveEulerSolver::Set_MPI_Solution_Old(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  #ifdef HAVE_MPI
	 int send_to, receive_from;
   MPI_Status status;
  #endif

	for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
		if(config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE && config->GetMarker_All_SendRecv(iMarker) > 0) {
			MarkerS = iMarker;
      MarkerR = iMarker + 1;

      #ifdef HAVE_MPI
        send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
			  receive_from = std::abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      #endif

		  nVertexS = geometry->nVertex[MarkerS];
      nBufferS_Vector = nVertexS*nVar;
      nVertexR = geometry->nVertex[MarkerR];

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Send_U.resize(nBufferS_Vector);

      /*--- Copy the solution that should be sended ---*/
      for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for(iVar = 0; iVar < nVar; ++iVar)
          Buffer_Send_U[iVar*nVertexS + iVertex] = node[iPoint]->GetSolution_Old(iVar);
      }

      #ifdef HAVE_MPI
        unsigned long nBufferR_Vector = nVertexR*nVar;
        Buffer_Receive_U.resize(nBufferR_Vector);
		    /*--- Send/Receive information using Sendrecv ---*/
        SU2_MPI::Sendrecv(Buffer_Send_U.data(), nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                          Buffer_Receive_U.data(), nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      #else
        Buffer_Receive_U = Buffer_Send_U;
      #endif

      /*--- Do the coordinate transformation ---*/
      for(iVertex = 0; iVertex < nVertexR; ++iVertex) {
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        auto angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];        phi    = angles[1];          psi    = angles[2];
        cosTheta = std::cos(theta);  cosPhi = std::cos(phi);      cosPsi = std::cos(psi);
        sinTheta = std::sin(theta);  sinPhi = std::sin(phi);      sinPsi = std::sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
        ordering is rotation about the x-axis, y-axis,then z-axis. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;

        rotMatrix[0][1] = cosPhi*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;

        rotMatrix[0][2] = -sinPhi;
        rotMatrix[1][2] = sinTheta*cosPhi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for(iVar = 0; iVar < nVar; ++iVar)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Rotate the momentum components. ---*/
        if(nDim == 2) {
          Solution[RHOVX_INDEX_SOL] = rotMatrix[0][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                      rotMatrix[0][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[1][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL+1)*nVertexR + iVertex];
        }
        else {
          Solution[RHOVX_INDEX_SOL] = rotMatrix[0][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                      rotMatrix[0][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                      rotMatrix[0][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[1][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                          rotMatrix[1][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];

          Solution[RHOVX_INDEX_SOL + 2] = rotMatrix[2][0]*Buffer_Receive_U[RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                          rotMatrix[2][1]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 1)*nVertexR + iVertex] +
                                          rotMatrix[2][2]*Buffer_Receive_U[(RHOVX_INDEX_SOL + 2)*nVertexR + iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        node[iPoint]->SetSolution_Old(Solution);
      } /*--- End of iVertex for loop ---*/
    }
	} /*--- End of iMarker for loop ---*/
}

//
//
/*--- Call MPI to set limiter for parallel simulations ---*/
//
//
void CReactiveEulerSolver::Set_MPI_Primitive_Limiter(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint;
  unsigned long nVertexS, nVertexR, nBufferS_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  /*--- Declare vector to store the limiter ---*/
  su2double Limiter[nPrimVarLim];

  #ifdef HAVE_MPI
    int send_to, receive_from;
    MPI_Status status;
  #endif

  /*--- Loop over all send/receive boundaries ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE && config->GetMarker_All_SendRecv(iMarker) > 0) {
      MarkerS = iMarker;
      MarkerR = iMarker + 1;

      #ifdef HAVE_MPI
        send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
        receive_from = std::abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      #endif

			  nVertexS = geometry->nVertex[MarkerS];
        nBufferS_Vector = nVertexS*nPrimVarLim;
        nVertexR = geometry->nVertex[MarkerR];

        /*--- Allocate Receive and send buffers  ---*/
        Buffer_Send_U.resize(nBufferS_Vector);

        /*--- Copy the solution old that should be sended ---*/
        for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
          iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
          for(iVar = 0; iVar < nPrimVarLim; ++iVar)
            Buffer_Send_U[iVar*nVertexS + iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
        }

        #ifdef HAVE_MPI
          unsigned long nBufferR_Vector = nVertexR*nPrimVarLim;
          Buffer_Receive_U.resize(nBufferR_Vector);
          /*--- Send/Receive information using Sendrecv ---*/
          SU2_MPI::Sendrecv(Buffer_Send_U.data(), nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                            Buffer_Receive_U.data(), nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
        #else
          Buffer_Receive_U = Buffer_Send_U;
        #endif

        /*--- Do the coordinate transformation ---*/
        for(iVertex = 0; iVertex < nVertexR; ++iVertex) {
          /*--- Find point and its type of transformation ---*/
          iPoint          = geometry->vertex[MarkerR][iVertex]->GetNode();
          iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

          /*--- Retrieve the supplied periodic information. ---*/
          auto angles = config->GetPeriodicRotation(iPeriodic_Index);

          /*--- Store angles separately for clarity. ---*/
          theta    = angles[0];        phi    = angles[1];      psi    = angles[2];
          cosTheta = std::cos(theta);  cosPhi = std::cos(phi);  cosPsi = std::cos(psi);
          sinTheta = std::sin(theta);  sinPhi = std::sin(phi);  sinPsi = std::sin(psi);

          /*--- Compute the rotation matrix. Note that the implicit
          ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;

          rotMatrix[0][1] = cosPhi*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;

          rotMatrix[0][2] = -sinPhi;
          rotMatrix[1][2] = sinTheta*cosPhi;
          rotMatrix[2][2] = cosTheta*cosPhi;

          /*--- Copy conserved variables before performing transformation. ---*/
          for(iVar = 0; iVar < nPrimVarLim; ++iVar)
            Limiter[iVar] = Buffer_Receive_U[iVar*nVertexR + iVertex];

          /*--- Rotate the momentum components. ---*/
          if(nDim == 2) {
            Limiter[VX_INDEX_LIM] =  rotMatrix[0][0]*Buffer_Receive_U[VX_INDEX_GRAD*nVertexR + iVertex] +
                                     rotMatrix[0][1]*Buffer_Receive_U[(VX_INDEX_GRAD + 1)*nVertexR + iVertex];

            Limiter[VX_INDEX_LIM + 1] = rotMatrix[1][0]*Buffer_Receive_U[VX_INDEX_GRAD*nVertexR + iVertex] +
                                        rotMatrix[1][1]*Buffer_Receive_U[(VX_INDEX_GRAD + 1)*nVertexR + iVertex];
          }

          else {
            Limiter[VX_INDEX_LIM] =  rotMatrix[0][0]*Buffer_Receive_U[VX_INDEX_GRAD*nVertexR + iVertex] +
                                     rotMatrix[0][1]*Buffer_Receive_U[(VX_INDEX_GRAD + 1)*nVertexR + iVertex] +
                                     rotMatrix[0][2]*Buffer_Receive_U[(VX_INDEX_GRAD + 2)*nVertexR + iVertex];

            Limiter[VX_INDEX_LIM + 1] = rotMatrix[1][0]*Buffer_Receive_U[VX_INDEX_GRAD*nVertexR + iVertex] +
                                        rotMatrix[1][1]*Buffer_Receive_U[(VX_INDEX_GRAD + 1)*nVertexR + iVertex] +
                                        rotMatrix[1][2]*Buffer_Receive_U[(VX_INDEX_GRAD + 2)*nVertexR + iVertex];

            Limiter[VX_INDEX_LIM + 2] =  rotMatrix[2][0]*Buffer_Receive_U[VX_INDEX_GRAD*nVertexR + iVertex] +
                                         rotMatrix[2][1]*Buffer_Receive_U[(VX_INDEX_GRAD + 1)*nVertexR + iVertex] +
                                         rotMatrix[2][2]*Buffer_Receive_U[(VX_INDEX_GRAD + 2)*nVertexR + iVertex];
        }

        /*--- Copy transformed limited variables back into buffer. ---*/
        for(iVar = 0; iVar < nPrimVarLim; ++iVar)
          node[iPoint]->SetLimiter_Primitive(iVar, Limiter[iVar]);

      } /*--- End if iVertex for loop ---*/
    }
  } /*--- End if iMarker for loop ---*/
}

//
//
/*--- Call MPI to set gradient of primitive variables for parallel simulations. ---*/
//
//
void CReactiveEulerSolver::Set_MPI_Primitive_Gradient(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  /*--- Declare matrix to store the matrix ---*/
  su2double Gradient[nPrimVarGrad][nDim];

  #ifdef HAVE_MPI
    int send_to, receive_from;
    MPI_Status status;
  #endif

	for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
		if(config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE && config->GetMarker_All_SendRecv(iMarker) > 0) {
			MarkerS = iMarker;
      MarkerR = iMarker + 1;

      #ifdef HAVE_MPI
        send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
        receive_from = std::abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      #endif

			  nVertexS = geometry->nVertex[MarkerS];
			  nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;
        nVertexR = geometry->nVertex[MarkerR];

        /*--- Allocate Receive and send buffers  ---*/
        Buffer_Send_U.resize(nBufferS_Vector);

        /*--- Copy the solution old that should be sended ---*/
        for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
          iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
          for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
            for(iDim = 0; iDim < nDim; ++iDim)
              Buffer_Send_U[iDim*nPrimVarGrad*nVertexS + iVar*nVertexS + iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
        }

        #ifdef HAVE_MPI
          unsigned long nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
          Buffer_Receive_U.resize(nBufferR_Vector);
          /*--- Send/Receive information using Sendrecv ---*/
          SU2_MPI::Sendrecv(Buffer_Send_U.data(), nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                            Buffer_Receive_U.data(), nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
        #else
          Buffer_Receive_U = Buffer_Send_U;
        #endif

        /*--- Do the coordinate transformation ---*/
        for(iVertex = 0; iVertex < nVertexR; ++iVertex) {
          /*--- Find point and its type of transformation ---*/
          iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
          iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

          /*--- Retrieve the supplied periodic information. ---*/
          auto angles = config->GetPeriodicRotation(iPeriodic_Index);

          /*--- Store angles separately for clarity. ---*/
          theta    = angles[0];        phi    = angles[1];      psi    = angles[2];
          cosTheta = std::cos(theta);  cosPhi = std::cos(phi);  cosPsi = std::cos(psi);
          sinTheta = std::sin(theta);  sinPhi = std::sin(phi);  sinPsi = std::sin(psi);

          /*--- Compute the rotation matrix. Note that the implicit
          ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;

          rotMatrix[0][1] = cosPhi*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;

          rotMatrix[0][2] = -sinPhi;
          rotMatrix[1][2] = sinTheta*cosPhi;
          rotMatrix[2][2] = cosTheta*cosPhi;

          /*--- Copy conserved variables before performing transformation. ---*/
          for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
            for(iDim = 0; iDim < nDim; ++iDim)
              Gradient[iVar][iDim] = Buffer_Receive_U[iDim*nPrimVarGrad*nVertexR + iVar*nVertexR+iVertex];

          /*--- Need to rotate the gradients for all conserved variables. ---*/
          for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
            if(nDim == 2) {
              Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_U[0*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[0][1]*Buffer_Receive_U[1*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex];

              Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_U[0*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[1][1]*Buffer_Receive_U[1*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex];
            }
            else {
              Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_U[0*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[0][1]*Buffer_Receive_U[1*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[0][2]*Buffer_Receive_U[2*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex];

              Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_U[0*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[1][1]*Buffer_Receive_U[1*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[1][2]*Buffer_Receive_U[2*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex];

              Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_U[0*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[2][1]*Buffer_Receive_U[1*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex] +
                                  rotMatrix[2][2]*Buffer_Receive_U[2*nPrimVarGrad*nVertexR + iVar*nVertexR + iVertex];
            }
         }

         /*--- Store the received information ---*/
         for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
           for(iDim = 0; iDim < nDim; ++iDim)
             node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);

      } /*--- End of iVertex for loop ---*/
    }
	} /*--- End of iMarker for loop ---*/
}

//
//
/*--- Set time step according to maximum inviscid eigenvalue. ---*/
//
//
void CReactiveEulerSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned long Iteration) {
  /*--- Local variables ---*/
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker;

  su2double Area, Volume, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda;
  su2double Local_Delta_Time, Global_Delta_Time = 1.0e6;
  su2double Normal[nDim];

  bool time_steping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  Min_Delta_Time = 1.0e6;
  Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero and compute sound speed ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    node[iPoint]->SetMax_Lambda_Inv(0.0);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    geometry->edge[iEdge]->GetNormal(Normal);
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5*(node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5*(node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());;

    /*--- Adjustment for grid movement ---*/
    if(grid_movement) {
      auto GridVel_i = geometry->node[iPoint]->GetGridVel();
      auto GridVel_j = geometry->node[jPoint]->GetGridVel();
      su2double ProjVel_i = std::inner_product(GridVel_i, GridVel_i + nDim, Normal, 0.0);
      su2double ProjVel_j = std::inner_product(GridVel_j, GridVel_j + nDim, Normal, 0.0);
      Mean_ProjVel -= 0.5*(ProjVel_i + ProjVel_j);
    }

    /*--- Inviscid contribution ---*/
    Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;

    if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);
  } /*--- End loop interior edges ---*/

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed();

        /*--- Adjustment for grid movement ---*/
        if(grid_movement) {
          auto GridVel = geometry->node[iPoint]->GetGridVel();
          su2double ProjVel = std::inner_product(GridVel, GridVel + nDim, Normal, 0.0);
          Mean_ProjVel -= ProjVel;
        }

        /*--- Inviscid contribution ---*/
        Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
        if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
    }
  } /*--- End loop boundary edges ---*/

  /*--- Each element uses their own speed for a steady state simulation ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();

    if(Volume > EPS) {
      Local_Delta_Time = config->GetCFL(iMesh)*Volume/node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
      if(Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      su2double rbuf_time, sbuf_time;
      sbuf_time = Min_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Min_Delta_Time = rbuf_time;

      sbuf_time = Max_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Max_Delta_Time = rbuf_time;
    #endif
  }

  /*--- Check if there is any element with only one neighbor...
        a CV that is inside another CV ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		if(geometry->node[iPoint]->GetnPoint() == 1)
			node[iPoint]->SetDelta_Time(Min_Delta_Time);
	}

  if(time_steping) {
    #ifdef HAVE_MPI
      su2double rbuf_time, sbuf_time;
      sbuf_time = Global_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Global_Delta_Time = rbuf_time;
    #endif
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      /*--- Sets the regular CFL equal to the unsteady CFL ---*/
      config->SetCFL(iMesh, config->GetUnst_CFL());

      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
            it computes the time step based on the unsteady CFL ---*/
      if(config->GetCFL(iMesh) < EPS)
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      else
        node[iPoint]->SetDelta_Time(Global_Delta_Time);
    }
  } /*--- End of time stepping checking ---*/

  if(dual_time) {
    su2double Global_Delta_UnstTimeND;
    if(Iteration == 0 && config->GetUnst_CFL() != 0.0 && iMesh == MESH_0) {
      Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

      #ifdef HAVE_MPI
        su2double rbuf_time, sbuf_time;
        sbuf_time = Global_Delta_UnstTimeND;
        SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        Global_Delta_UnstTimeND = rbuf_time;
      #endif

      config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
    }

    /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      if(!implicit) {
        Local_Delta_Time = std::min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  } /*--- End of dual time checking ---*/
}

//
//
/*--- Set residual for dual time strategies. ---*/
//
//
void CReactiveEulerSolver::SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                                unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  if(config->GetExtIter() == 0)
    SU2_Assert(Residual != NULL, "The array for residual has not been allocated");

  /*--- Local variables ---*/
  unsigned short iVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  su2double Volume_nM1, Volume_nP1, TimeStep, Residual_GCL;
  su2double Normal[nDim];

  /*--- Store the physical time step ---*/
  TimeStep = config->GetDelta_UnstTimeND();

  /*--- Compute the dual time-stepping source term for static meshes ---*/
  if(!grid_movement) {
    /*--- Loop over all nodes (excluding halos) ---*/
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
            we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
            previous solutions that are stored in memory. ---*/
      auto U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      auto U_time_n   = node[iPoint]->GetSolution_time_n();
      auto U_time_nP1 = node[iPoint]->GetSolution();

      /*--- CV volume at time n+1. As we are on a static mesh, the volume
            of the CV will remained fixed for all time steps. ---*/
      Volume_nP1 = geometry->node[iPoint]->GetVolume();

      /*--- Compute the dual time-stepping source term based on the chosen
            time discretization scheme (1st- or 2nd-order).---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*Volume_nP1 / TimeStep;
        if(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (3.0*U_time_nP1[iVar] - 4.0*U_time_n[iVar] + U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }

      /*--- Store the residual and compute the Jacobian contribution due  to the dual time source term. ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      if(implicit) {
        for(iVar = 0; iVar < nVar; ++iVar) {
          std::fill(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, 0.0);
          if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

  else {
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
          dynamically deforming), the Geometric Conservation Law (GCL) should be
          satisfied in conjunction with the ALE formulation of the governing
          equations. The GCL prevents accuracy issues caused by grid motion, i.e.
          a uniform free-stream should be preserved through a moving grid. First,
          we will loop over the edges and boundaries to compute the GCL component
          of the dual time source term that depends on grid velocities. ---*/
    for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
      /*--- Get indices for nodes i & j plus the face normal ---*/
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      geometry->edge[iEdge]->GetNormal(Normal);

      /*--- Grid velocities stored at nodes i & j ---*/
      auto GridVel_i = geometry->node[iPoint]->GetGridVel();
      auto GridVel_j = geometry->node[jPoint]->GetGridVel();

      /*--- Compute the GCL term by averaging the grid velocities at the
            edge mid-point and dotting with the face normal. ---*/
      Residual_GCL = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim)
        Residual_GCL += 0.5*(GridVel_i[iDim] + GridVel_j[iDim])*Normal[iDim];

      /*--- Compute the GCL component of the source term for node i ---*/
      auto U_time_n = node[iPoint]->GetSolution_time_n();
      for(iVar = 0; iVar < nVar; ++iVar)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Compute the GCL component of the source term for node j ---*/
      U_time_n = node[jPoint]->GetSolution_time_n();
      for(iVar = 0; iVar < nVar; ++iVar)
        Residual[iVar] = U_time_n[iVar]*Residual_GCL;
      LinSysRes.SubtractBlock(jPoint, Residual);
    } /*---End loop interior edges ---*/

    /*---   Loop over the boundary edges ---*/
    for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
      if(config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
        for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          /*--- Get the index for node i plus the boundary face normal ---*/
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

          /*--- Grid velocities stored at boundary node i ---*/
          auto GridVel = geometry->node[iPoint]->GetGridVel();

          /*--- Compute the GCL term by dotting the grid velocity with the face normal.
                The normal is negated to match the boundary convention. ---*/
          Residual_GCL = -std::inner_product(GridVel, GridVel + nDim, Normal, 0.0);

          /*--- Compute the GCL component of the source term for node i ---*/
          auto U_time_n = node[iPoint]->GetSolution_time_n();
          for(iVar = 0; iVar < nVar; ++iVar)
            Residual[iVar] = U_time_n[iVar]*Residual_GCL;
          LinSysRes.AddBlock(iPoint, Residual);
        }
      }
    } /*--- End loop boundary edges ---*/

    /*--- Loop over all nodes (excluding halos) to compute the remainder  of the dual time-stepping source term. ---*/
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
            we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
            previous solutions that are stored in memory. ---*/
      auto U_time_nM1 = node[iPoint]->GetSolution_time_n1();
      auto U_time_n   = node[iPoint]->GetSolution_time_n();
      auto U_time_nP1 = node[iPoint]->GetSolution();

      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
            grids, the volumes will change. On rigidly transforming grids, the
            volumes will remain constant. ---*/
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();

      /*--- Compute the dual time-stepping source residual. Due to the
            introduction of the GCL term above, the remainder of the source residual
            due to the time discretization has a new form.---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(Volume_nP1/TimeStep);
        if(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Residual[iVar] = (U_time_nP1[iVar] - U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep)) +
                           (U_time_nM1[iVar] - U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }

      /*--- Store the residual and compute the Jacobian contribution due
            to the dual time source term. ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      if(implicit) {
        for(iVar = 0; iVar < nVar; ++iVar) {
          std::fill(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, 0.0);
          if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if(config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    } /*--- End loop over all nodes ---*/
  } /*--- End of dynamic mesh part ---*/
}


//
//
/*--- Iteration of implicit Euler method. ---*/
//
//
void CReactiveEulerSolver::ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol = 0;

  su2double Volume, Delta;

  /*--- Set maximum residual to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Read the residual ---*/
    auto Local_Res_TruncError = node[iPoint]->GetResTruncError();
    /*--- Read the volume ---*/
    Volume = geometry->node[iPoint]->GetVolume();

    if(node[iPoint]->GetDelta_Time() > EPS) {
      Delta = Volume/node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for(iVar = 0; iVar < nVar; ++iVar) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        Local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + Local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, std::abs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/
  for(iPoint = nPointDomain; iPoint < nPoint; ++iPoint) {
    for(iVar = 0; iVar < nVar; ++iVar) {
			total_index = iPoint*nVar + iVar;
			LinSysRes[total_index] = 0.0;
			LinSysSol[total_index] = 0.0;
		}
  }

  /*--- Solve or smooth the linear system ---*/
  CSysSolve system;
  IterLinSol = system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/
  SetIterLinSolver(IterLinSol);

  /*--- Update solution (system written in terms of increments) ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		for(iVar = 0; iVar < nVar; ++iVar)
			node[iPoint]->AddClippedSolution(iVar, config->GetRelaxation_Factor_Flow()*LinSysSol[iPoint*nVar+iVar],
                                       Lower_Limit[iVar], Upper_Limit[iVar]);
  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*--- Iteration of explicit Euler method. ---*/
//
//
void CReactiveEulerSolver::ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar;
  unsigned long iPoint;

  su2double Res, Volume, Delta = 0.0;

  /*--- Set residual to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();
    if(Volume > EPS)
      Delta = node[iPoint]->GetDelta_Time()/Volume;

    auto Local_Res_TruncError = node[iPoint]->GetResTruncError();
    auto Local_Residual = LinSysRes.GetBlock(iPoint);

    for(iVar = 0; iVar < nVar; ++iVar) {
      Res = Local_Residual[iVar] + Local_Res_TruncError[iVar];
      node[iPoint]->AddClippedSolution(iVar, -Res*Delta, Lower_Limit[iVar], Upper_Limit[iVar]);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, std::abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*--- Iteration of explicit RK method. ---*/
//
//
void CReactiveEulerSolver::ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep) {
  /*--- Local variables ---*/
  unsigned short iVar;
  unsigned long iPoint;

  su2double Res, Volume, Delta = 0.0;

  auto RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

  /*--- Set residual to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();
    if(Volume > EPS)
      Delta = node[iPoint]->GetDelta_Time()/Volume;

    auto Local_Res_TruncError = node[iPoint]->GetResTruncError();
    auto Local_Residual = LinSysRes.GetBlock(iPoint);

    for(iVar = 0; iVar < nVar; ++iVar) {
      Res = Local_Residual[iVar] + Local_Res_TruncError[iVar];
      node[iPoint]->AddClippedSolution(iVar, -Res*Delta*RK_AlphaCoeff, Lower_Limit[iVar], Upper_Limit[iVar]);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, std::abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*--- Set freestream solution all over the domain. ---*/
//
//
void CReactiveEulerSolver::SetFreeStream_Solution(CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint;
  unsigned short iDim, iSpecies;

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Set density ---*/
    node[iPoint]->SetSolution(RHO_INDEX_SOL, Density_Inf);
    /*--- Set momentum ---*/
    for(iDim = 0; iDim < nDim; ++iDim)
      node[iPoint]->SetSolution(RHOVX_INDEX_SOL + iDim, Density_Inf*Velocity_Inf[iDim]);
    /*--- Set energy ---*/
    node[iPoint]->SetSolution(RHOE_INDEX_SOL, Density_Inf*config->GetEnergy_FreeStreamND());
    /*--- Set partial densities ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      node[iPoint]->SetSolution(RHOS_INDEX_SOL + iSpecies, Density_Inf*MassFrac_Inf[iSpecies]);
  }
}

//
//
/*--- Centered residual convective term. ---*/
//
//
void CReactiveEulerSolver::Centered_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                             CConfig* config, unsigned short iMesh, unsigned short iRKStep) {

  throw Common::NotImplemented("Function not implemented: Centered Residual\n");
}

//
//
/*--- Compute Upwind residual convective term. ---*/
//
//
void CReactiveEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                           CConfig* config, unsigned short iMesh) {
  /*--- Local variables ---*/
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar, iSpecies;

  su2double Project_Grad_i, Project_Grad_j;

  /*--- Loop over all the edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Retrive points in the edge and normal vector ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Get primitive variables ---*/
    auto V_i = node[iPoint]->GetPrimitive();
    auto V_j = node[jPoint]->GetPrimitive();

    if(second_order) {
      if(config->GetExtIter() == 0) {
        SU2_Assert(Vector_i != NULL, "The array to compute the vector from i to j has not been allocated");
        SU2_Assert(Vector_j != NULL, "The array to compute the vector from j to i has not been allocated");
      }
      for(iDim = 0; iDim < nDim; ++iDim) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = -Vector_i[iDim];
      }
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();

      /*--- Auxiliary vector to store reconstruction ---*/
      su2double Prim_Recon_i[nPrimVarLim];
      su2double Prim_Recon_j[nPrimVarLim];

      Prim_Recon_i[T_INDEX_LIM] = V_i[T_INDEX_PRIM];
      Prim_Recon_j[T_INDEX_LIM] = V_j[T_INDEX_PRIM];
      Prim_Recon_i[P_INDEX_LIM] = V_i[P_INDEX_PRIM];
      Prim_Recon_j[P_INDEX_LIM] = V_j[P_INDEX_PRIM];
      for(iDim = 0; iDim < nDim; ++iDim) {
        Prim_Recon_i[VX_INDEX_LIM + iDim] = V_i[VX_INDEX_PRIM + iDim];
        Prim_Recon_j[VX_INDEX_LIM + iDim] = V_j[VX_INDEX_PRIM + iDim];
      }

      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        Project_Grad_i = Project_Grad_j = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }

        if(limiter) {
          auto Limiter_i = node[iPoint]->GetLimiter_Primitive();
          auto Limiter_j = node[jPoint]->GetLimiter_Primitive();
          Prim_Recon_i[iVar] += Limiter_i[iVar]*Project_Grad_i;
          Prim_Recon_j[iVar] += Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Prim_Recon_i[iVar] += Project_Grad_i;
          Prim_Recon_j[iVar] += Project_Grad_j;
        }
      }

      /*--- Check for non physical solution after reconstruction ---*/
      bool non_phys_i = false;
      if(Prim_Recon_i[T_INDEX_LIM] > EPS)
        Primitive_i[T_INDEX_PRIM] = Prim_Recon_i[T_INDEX_LIM];
      else
        non_phys_i = true;
      if(!non_phys_i) {
        if(Prim_Recon_i[P_INDEX_LIM] > EPS)
          Primitive_i[P_INDEX_PRIM] = Prim_Recon_i[P_INDEX_LIM];
        else
          non_phys_i = true;
      }

      bool non_phys_j = false;
      if(Prim_Recon_j[T_INDEX_LIM] > EPS)
        Primitive_j[T_INDEX_PRIM] = Prim_Recon_j[T_INDEX_LIM];
      else
        non_phys_j = true;
      if(!non_phys_j) {
        if(Prim_Recon_i[P_INDEX_LIM] > EPS)
          Primitive_j[P_INDEX_PRIM] = Prim_Recon_j[P_INDEX_LIM];
        else
          non_phys_j = true;
      }

      /*--- Compute other primitive variables accordingly to the reconstruction ---*/
      su2double dim_temp_i, Gamma_i, sq_vel_i;
      if(!non_phys_i) {
        /*--- Copy reconstructed velocities ---*/
        std::copy(Prim_Recon_i + VX_INDEX_LIM, Prim_Recon_i + (VX_INDEX_LIM + nDim), Primitive_i.begin() + VX_INDEX_PRIM);

        /*--- Build a thermodynamic consistent state from the reconstructed variables ---*/
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Primitive_i[RHOS_INDEX_PRIM + iSpecies] = V_i[RHOS_INDEX_PRIM + iSpecies];
          Ys_i[iSpecies] = V_i[RHOS_INDEX_PRIM + iSpecies];
        }
        su2double rho_recon_i = library->ComputeDensity(Primitive_i[T_INDEX_PRIM], Primitive_i[P_INDEX_PRIM], Ys_i);
        rho_recon_i *= config->GetGas_Constant_Ref();
        Primitive_i[RHO_INDEX_PRIM] = rho_recon_i;
        dim_temp_i = Primitive_i[T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp_i *= 5.0/9.0;
        Primitive_i[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp_i, Ys_i)/config->GetEnergy_Ref();
        if(US_System)
          Primitive_i[H_INDEX_PRIM] *= 3.28084*3.28084;
        sq_vel_i = std::inner_product(Primitive_i.cbegin() + VX_INDEX_PRIM, Primitive_i.cbegin() + (VX_INDEX_PRIM + nDim),
                                      Primitive_i.cbegin() + VX_INDEX_PRIM, 0.0);
        Primitive_i[H_INDEX_PRIM] += 0.5*sq_vel_i;
        Gamma_i = library->ComputeFrozenGamma(dim_temp_i, Ys_i);
        Primitive_i[A_INDEX_PRIM] = std::sqrt(Gamma_i*Primitive_i[P_INDEX_PRIM]/rho_recon_i);
      }
      else
        std::copy(V_i, V_i + nPrimVar, Primitive_i.data());

      /*--- Compute other primitive variables accordingly to the reconstruction ---*/
      su2double dim_temp_j, Gamma_j, sq_vel_j;
      if(!non_phys_j) {
        /*--- Copy reconstructed velocities ---*/
        std::copy(Prim_Recon_j + VX_INDEX_LIM, Prim_Recon_j + (VX_INDEX_LIM + nDim), Primitive_j.begin() + VX_INDEX_PRIM);

        /*--- Build a thermodynamic consistent state from the recosntructed variables ---*/
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Primitive_j[RHOS_INDEX_PRIM + iSpecies] = V_j[RHOS_INDEX_PRIM + iSpecies];
          Ys_j[iSpecies] = V_j[RHOS_INDEX_PRIM + iSpecies];
        }
        su2double rho_recon_j = library->ComputeDensity(Primitive_j[T_INDEX_PRIM], Primitive_j[P_INDEX_PRIM], Ys_j);
        rho_recon_j *= config->GetGas_Constant_Ref();
        Primitive_j[RHO_INDEX_PRIM] = rho_recon_j;
        dim_temp_j = Primitive_j[T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp_j *= 5.0/9.0;
        Primitive_j[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp_j, Ys_j)/config->GetEnergy_Ref();
        if(US_System)
          Primitive_j[H_INDEX_PRIM] *= 3.28084*3.28084;
        sq_vel_j = std::inner_product(Primitive_j.cbegin() + VX_INDEX_PRIM, Primitive_j.cbegin() + (VX_INDEX_PRIM + nDim),
                                      Primitive_j.cbegin() + VX_INDEX_PRIM, 0.0);
        Primitive_j[H_INDEX_PRIM] += 0.5*sq_vel_j;
        Gamma_j = library->ComputeFrozenGamma(dim_temp_j, Ys_j);
        Primitive_j[A_INDEX_PRIM] = std::sqrt(Gamma_j*Primitive_j[P_INDEX_PRIM]/rho_recon_j);
      }
      else
        std::copy(V_j, V_j + nPrimVar, Primitive_j.data());

      numerics->SetPrimitive(Primitive_i.data(), Primitive_j.data());

      if(implicit) {
        /*--- Set pressure derivatives if everything is ok... ---*/
        if(!non_phys_i) {
          /*--- Derivatives with respect to density ---*/
          Secondary_i[RHO_INDEX_SOL] = (Gamma_i - 1.0)*0.5*sq_vel_i;

          /*--- Derivatives with respect to momentum ---*/
          for(iDim = 0; iDim < nDim; ++iDim)
            Secondary_i[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma_i)*Primitive_i[VX_INDEX_PRIM + iDim];

          /*--- Derivatives with respect to energy ---*/
          Secondary_i[RHOE_INDEX_SOL] = Gamma_i - 1.0;

          /*--- Derivatives with respect to partial densities ---*/
          for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            Secondary_i[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp_i, Gamma_i, iSpecies)/config->GetEnergy_Ref();
            if(US_System)
              Secondary_i[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
          }
        }
        /*--- ...otherwise use the original values. ---*/
        else
          std::copy(node[iPoint]->GetdPdU(), node[iPoint]->GetdPdU() + nVar, Secondary_i.data());

        if(!non_phys_j) {
          /*--- Derivatives with respect to density ---*/
          Secondary_j[RHO_INDEX_SOL] = (Gamma_j - 1.0)*0.5*sq_vel_j;
          /*--- Derivatives with respect to momentum ---*/
          for(iDim = 0; iDim < nDim; ++iDim)
            Secondary_j[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma_j)*Primitive_j[VX_INDEX_PRIM + iDim];

          /*--- Derivatives with respect to energy ---*/
          Secondary_j[RHOE_INDEX_SOL] = Gamma_j - 1.0;

          /*--- Derivatives with respect to partial densities ---*/
          for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Secondary_j[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp_j, Gamma_j, iSpecies)/config->GetEnergy_Ref();
          if(US_System)
            Secondary_j[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
          }
        }
        else
          std::copy(node[jPoint]->GetdPdU(), node[jPoint]->GetdPdU() + nVar, Secondary_j.data());

        numerics->SetSecondary(Secondary_i.data(), Secondary_j.data());
      }
    } /*--- End second order reconstruction ---*/

    else {
      /*--- Set primitive variables without reconstruction ---*/
      numerics->SetPrimitive(V_i, V_j);
      if(implicit)
        numerics->SetSecondary(node[iPoint]->GetdPdU(), node[jPoint]->GetdPdU());
    }

    /*--- Set grid velocity ---*/
    if(grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());

    /*--- Compute the residual ---*/
    numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = !std::none_of(Res_Conv, Res_Conv + nVar, [](su2double elem){return std::isnan(elem);});
    if(implicit) {
      if(!err) {
        for(iVar = 0; iVar < nVar; ++iVar) {
          err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
          err = err  || !std::none_of(Jacobian_j[iVar], Jacobian_j[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
          if(err)
            break;
        }
      }
    }

    /*--- Update residual value ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);

      /*--- Set implicit Jacobians ---*/
      if(implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
      }
    }
    else
      throw std::runtime_error("NaN found in the upwind residual");
  } /*--- End loop over edges ---*/

  /*--- Warning message about non-physical reconstructions ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  #else
    counter_global = counter_local;
  #endif
    if(iMesh == MESH_0)
      config->SetNonphysical_Reconstr(counter_global);
  }
}

//
//
/*--- Compute Residual source term. ---*/
//
//
void CReactiveEulerSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CNumerics* second_numerics,
                                           CConfig* config, unsigned short iMesh) {
  /*--- Local variables ---*/
  unsigned long iPoint, counter_local = 0, counter_global = 0;

  /*--- Set residual to zero for safety ---*/
  std::fill(Res_Sour, Res_Sour + nVar, 0.0);

  /*--- Loop over all points ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Load the primitive variables and temperature derivatives ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[iPoint]->GetPrimitive());
    if(implicit)
      numerics->SetSecondary(node[iPoint]->GetdTdU(), node[iPoint]->GetdTdU());

    /*--- Load the volume of the dual mesh cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Compute the source residual ---*/
    numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = !std::none_of(Res_Sour, Res_Sour + nVar, [](su2double elem){return std::isnan(elem);});
    if(implicit) {
      if(!err) {
        for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
          err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
          if(err)
            break;
        }
      }
    }

    /*--- Add the source residual to the total ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Res_Sour);

      /*--- Add the implicit Jacobian contribution ---*/
      if(implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
    else
      throw std::runtime_error("NaN found in the source residual");
  } /*--- End loop over points ---*/

  /*--- Warning message about non-physical reconstructions ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  #else
    counter_global = counter_local;
  #endif
    if(iMesh == MESH_0)
      config->SetNonphysical_Reconstr(counter_global);
  }
}

//
//
/*--- Solid wall boundary condition for Euler imposed in a weak form (i.e. via residual). ---*/
//
//
void CReactiveEulerSolver::BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                         CConfig* config, unsigned short val_marker) {
  if(config->GetExtIter() == 0)
    SU2_Assert(Residual != NULL,"The array to save residual for boundary conditions has not been allocated");

  unsigned long iVertex, iPoint;
  unsigned short iDim;

  su2double Area, Pressure;
  su2double Normal[nDim], UnitNormal[nDim];

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if(geometry->node[iPoint]->GetDomain()) {
      /*--- Set residual to zero for safety ---*/
      std::fill(Residual, Residual + nVar, 0.0);

			/*--- Normal vector for this vertex (negative for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

			/*--- Compute parameters from the geometry ---*/
      Area = ::ComputeArea(Normal, nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			/*--- Retrieve the pressure on the vertex ---*/
      Pressure = node[iPoint]->GetPressure();

			/*--- Compute the residual ---*/
			for(iDim = 0; iDim < nDim; ++iDim)
				Residual[RHOVX_INDEX_SOL + iDim] = Pressure*UnitNormal[iDim]*Area;

      /*--- Correction to energy due to mowing wall: extra term due to pressure contribution ---*/
      su2double ProjGridVel = 0.0;
      if(grid_movement) {
        auto GridVel = geometry->node[iPoint]->GetGridVel();
        ProjGridVel = std::inner_product(GridVel, GridVel + nDim, UnitNormal, 0.0);
        Residual[RHOE_INDEX_SOL] = Pressure*ProjGridVel*Area;
      }

      /*--- Add value to the residual ---*/
  		LinSysRes.AddBlock(iPoint, Residual);

      if(implicit) {
        unsigned short iVar, jVar;

        /*--- Set Jacobian to zero ---*/
        for(iVar = 0; iVar < nVar; ++iVar)
          for(jVar = 0; jVar < nVar; ++jVar)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Retrieve pressure derivatives ---*/
        auto dPdU = node[iPoint]->GetdPdU();

        /*--- Contribution of momentum equation ---*/
        for(iDim = 0; iDim < nDim; ++iDim)
          for(iVar = 0; iVar < nVar; ++iVar)
            Jacobian_i[RHOVX_INDEX_SOL + iDim][iVar] = dPdU[iVar]*UnitNormal[iDim]*Area;

        /*--- Extra contribution to energy in case of grid_movement ---*/
        if(grid_movement) {
          for(iVar = 0; iVar < nVar; ++iVar)
            Jacobian_i[RHOE_INDEX_SOL][iVar] = dPdU[iVar]*ProjGridVel*Area;
        }

        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }

    } /*--- End of if Geometry() ---*/

  } /*--- End of iVertex for loop ---*/
}

//
//
/*--- Impose supersonic inlet boundary condition. ---*/
//
//
void CReactiveEulerSolver::BC_Supersonic_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                               CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  /*--- Local variables ---*/
  unsigned long iVertex, iPoint, Point_Normal;

  su2double V_aux[nPrimVar];
  su2double Normal[nDim];

  std::string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool viscous = config->GetViscous();

  su2double Density, Pressure, Temperature, Enthalpy, SoundSpeed;

  /*--- Supersonic inlet flow: there are no outgoing characteristics,
        so all flow variables can be imposed at the inlet.
        First, retrieve the specified values for the primitive variables. ---*/
  Ys = RealVec(config->GetInlet_MassFrac(Marker_Tag), config->GetInlet_MassFrac(Marker_Tag) + nSpecies);
  if(config->GetExtIter() == 0) {
    SU2_Assert(config->GetnSpecies_Inlet() == nSpecies, "Wrong number of species detected at inlet");
    SU2_Assert(std::abs(std::accumulate(Ys.cbegin(), Ys.cend(), 0.0) - 1.0) < 1.0e-10,
               "The mass fractions don't sum up to 1 in inlet boundary");
  }
  Temperature = config->GetInlet_Temperature(Marker_Tag);
	Pressure    = config->GetInlet_Pressure(Marker_Tag);
	auto Velocity = config->GetInlet_Velocity(Marker_Tag);

	/*--- Density at the inlet from the gas law ---*/
	Density = library->ComputeDensity(Temperature, Pressure, Ys);

  /*--- Compute the energy from the specified state ---*/
  su2double dim_temp = Temperature;
  if(US_System)
    dim_temp *= 5.0/9.0;
  Enthalpy = library->ComputeEnthalpy(dim_temp, Ys);
  SoundSpeed = library->ComputeFrozenSoundSpeed(dim_temp, Ys);
  if(US_System) {
    Enthalpy *= 3.28084*3.28084;
    SoundSpeed *= 3.28084;
  }

  /*--- Non-dimensionalize the inputs if necessary. ---*/
	Temperature /= config->GetTemperature_Ref();
	Pressure /= config->GetPressure_Ref();
	Density /= config->GetDensity_Ref();
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    Velocity[iDim] /= config->GetVelocity_Ref();
  Enthalpy /= config->GetEnergy_Ref();
  SoundSpeed /= config->GetVelocity_Ref();

  /*--- Primitive variables using the derived quantities ---*/
  V_aux[T_INDEX_PRIM] = Temperature;
  std::copy(Velocity, Velocity + nDim, V_aux + VX_INDEX_PRIM);
  V_aux[P_INDEX_PRIM] = Pressure;
  V_aux[RHO_INDEX_PRIM] = Density;
  su2double Velocity2 = std::inner_product(Velocity, Velocity + nDim, Velocity, 0.0);
  V_aux[H_INDEX_PRIM] = Enthalpy + 0.5*Velocity2;
  V_aux[A_INDEX_PRIM] = SoundSpeed;
  std::copy(Ys.cbegin(), Ys.cend(), V_aux + RHOS_INDEX_PRIM);

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if(geometry->node[iPoint]->GetDomain()) {
       /*--- Normal vector for this vertex (negative for outward convention) ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      std::transform(Normal, Normal + nDim, Normal, std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      /*--- Pull primitve variable from actual point ---*/
      auto V_domain = node[iPoint]->GetPrimitive();

      /*--- Get characteristic variables on the boundary ---*/
      auto V_inlet = GetCharacPrimVar(val_marker, iVertex);
      std::copy(V_aux, V_aux + nPrimVar, V_inlet);

      /*--- Set primitive variables for residual ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      if(implicit) {
        su2double Secondary[nVar];

        /*--- Derivatves with respect to density ---*/
        su2double Gamma = library->ComputeFrozenGamma_FromSoundSpeed(V_inlet[T_INDEX_PRIM], V_inlet[A_INDEX_PRIM], Ys);
        Gamma *= config->GetGas_Constant_Ref();
        Secondary[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*Velocity2;

        /*--- Derivatives with respect to momentum ---*/
        for(unsigned short iDim = 0; iDim < nDim; ++iDim)
          Secondary[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V_inlet[VX_INDEX_PRIM + iDim];

        /*--- Derivatives with respect to energy ---*/
        Secondary[RHOE_INDEX_SOL] = Gamma - 1.0;

        /*--- Derivatives with respect to partial densities ---*/
        for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Secondary[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp, Gamma, iSpecies)/config->GetEnergy_Ref();
          if(US_System)
            Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
        }

        conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), Secondary);
      }

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      bool err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      /*--- Apply the residual to the convective part ---*/
      if(!err) {
        LinSysRes.AddBlock(iPoint, Residual);

        if(implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the convective residual of supersonic inòet boundary condition");

      /*--- Viscous contribution ---*/
  		if(viscous) {
        /*--- Index of the closest interior node ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Set the normal vector and the coordinates ---*/
  			visc_numerics->SetNormal(Normal);
  			visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

  			/*--- Set primitive variables and gradient ---*/
  			visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

        /*--- Set temperature derivatives if implicit ---*/
        if(implicit) {
          su2double Secondary[nVar];

          su2double dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          su2double Cv = library->ComputeCV(dim_temp, Ys)/config->GetGas_Constant_Ref();
          if(US_System)
            Cv *= 3.28084*3.28084*5.0/9.0;
          su2double rhoCv = V_inlet[RHO_INDEX_PRIM]*Cv;
          su2double sq_vel = std::inner_product(V_inlet + VX_INDEX_PRIM, V_inlet + (VX_INDEX_PRIM + nDim), V_inlet + VX_INDEX_PRIM, 0.0);

          /*--- Set temperature derivatives ---*/
          Secondary[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
          for(unsigned short iDim = 0; iDim < nDim; ++iDim)
            Secondary[RHOVX_INDEX_SOL + iDim] = -V_inlet[VX_INDEX_PRIM + iDim]/rhoCv;
          Secondary[RHOE_INDEX_SOL] = 1.0/rhoCv;
          for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            Secondary[RHOS_INDEX_SOL + iSpecies] = -library->ComputePartialEnergy(dim_temp, iSpecies)/(config->GetEnergy_Ref()*rhoCv);
            if(US_System)
              Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
          }
          visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), Secondary);
        }

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[iPoint]->GetDiffusionCoeff());

        /*--- Compute the residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        /*--- Check for NaNs before applying the residual to the linear system ---*/
        err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
        if(implicit) {
          if(!err) {
            for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
              err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
              if(err)
                break;
            }
          }
        }

        /*--- Apply the residual to the viscous part ---*/
        if(!err) {
          LinSysRes.SubtractBlock(iPoint, Residual);

  			  /*--- Jacobian contribution for implicit integration ---*/
  			  if(implicit)
            Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        else
          throw std::runtime_error("NaN found in the viscous residual of supersonic inlet boundary condition.");
      }
    } /*--- End of if GetDomain() ---*/
  } /*--- End of iVertex for loop ---*/
}

//
//
/*!
  *\brief Inlet boundary conditions (Total conditions or mass flow according to specific simulation)
  */
//
//
void CReactiveEulerSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  throw Common::NotImplemented("Far field boundary condition not implemented yet for multispecies flows.");
}

//
//
/*--- Impose subsonic inlet boundary condition (Total conditions or mass flow according to specific simulation). ---*/
//
//
void CReactiveEulerSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  /*--- Local variables ---*/
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Riemann, SoundSpeed, Vn, Vel_Mag, Gamma, Gamma_Minus_One, dim_temp, alpha, Area;
  su2double Normal[nDim], UnitNormal[nDim];

  unsigned short Kind_Inlet = config->GetKind_Inlet();
  std::string Marker_Tag    = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();

  Ys = RealVec(config->GetInlet_MassFrac(Marker_Tag), config->GetInlet_MassFrac(Marker_Tag) + nSpecies);
  if(config->GetExtIter() == 0) {
    SU2_Assert(config->GetnMarker_Inlet() == config->GetnMarker_Inlet_MassFrac(), "Wrong number of inlet markers detected");
    SU2_Assert(config->GetnSpecies_Inlet() == nSpecies, "Wrong number of species detected at inlet");
    SU2_Assert(std::abs(std::accumulate(Ys.cbegin(), Ys.cend(), 0.0) - 1.0) < 1.0e-10,
               "The mass fractions don't sum up to 1 in inlet boundary");
  }
  auto Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      Area = ::ComputeArea(Normal,nDim);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Normal[iDim] = -Normal[iDim];
        UnitNormal[iDim] = Normal[iDim]/Area;
      }
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at this boundary node ---*/
      auto V_domain = node[iPoint]->GetPrimitive();

      /*--- Get characteristic variables on the boundary ---*/
      auto V_inlet = GetCharacPrimVar(val_marker, iVertex);

      /*--- Build the fictitious intlet state based on characteristics ---*/
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
            therefore we can specify all but one state variable at the inlet.
            The outgoing Riemann invariant provides the final piece of info. ---*/
      switch (Kind_Inlet) {
        /*--- Total properties have been specified at the inlet. ---*/
        case TOTAL_CONDITIONS: {
          /*--- Retrieve the specified mass flow for the inlet. ---*/
          su2double Ttot = config->GetInlet_Ttotal(Marker_Tag)/config->GetTemperature_Ref();
          su2double Ptot = config->GetInlet_Ptotal(Marker_Tag)/config->GetPressure_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/
          su2double Velocity2 = Vn = 0.0;
          for(iDim = 0; iDim < nDim; iDim++) {
            Velocity2 += V_domain[VX_INDEX_PRIM + iDim]*V_domain[VX_INDEX_PRIM + iDim];
            Vn += V_domain[VX_INDEX_PRIM + iDim]*UnitNormal[iDim];
          }
          SoundSpeed = V_domain[A_INDEX_PRIM];
          Gamma = node[iPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0;
          dim_temp = Ttot*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          su2double Gamma_Tot = library->ComputeFrozenGamma(dim_temp, Ys);

          /*--- Use harmonic average between gamma at node and the total gamma ---*/
          Gamma = 2.0/(1.0/Gamma + 1.0/Gamma_Tot);
          Gamma_Minus_One = Gamma - 1.0;
          Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

          /*--- Compute total enthalpy from total temperature ---*/
          su2double Tot_Enthalpy = library->ComputeEnthalpy(dim_temp, Ys);
          V_inlet[H_INDEX_PRIM] = Tot_Enthalpy;

          /*--- Compute inner product between unit normal and flow direction and save current gas constant ---*/
          alpha = std::inner_product(UnitNormal, UnitNormal + nDim, Flow_Dir, 0.0);
          su2double Rgas = library->ComputeRgas(Ys)/config->GetGas_Constant_Ref();

          /*--- Auxiliary function to impose adiabaticity ---*/
          auto f = std::function<su2double(su2double)>([&](su2double T){
            su2double hb = library->ComputeEnthalpy(T, Ys);
            su2double cb = std::sqrt(Gamma*Rgas*T);
            su2double Vb = (Riemann - 2.0*cb/Gamma_Minus_One)/alpha;
            return hb + 0.5*Vb*Vb;
          });

          /*--- Set parameters for secant method to find temperature ---*/
          bool NRconvg, Bconvg;
          su2double NRtol = 1.0e-9;    // Tolerance for the Secant method
          su2double Btol = 1.0e-6;    // Tolerance for the Bisection method
          unsigned short maxNIter = 15;        // Maximum Secant method iterations
          unsigned short maxBIter = 100;        // Maximum Bisection method iterations
          unsigned short iIter;

          su2double Told = Ttot + 1.0;
          su2double Tcurr = Ttot;
          su2double Tnew;
          NRconvg = false;

          /*--- Execute a secant root-finding method to find the outlet temperature (TRAPEZOIDAL) ---*/
          for(iIter = 0; iIter < maxNIter; ++iIter) {
            su2double tmp = f(Tcurr);
            su2double F = tmp - Tot_Enthalpy;
            su2double dF = tmp - f(Told);
            Tnew = Tcurr - F*(Tcurr - Told)/dF;

            /*--- Check for convergence ---*/
            if(std::abs(Tnew - Tcurr) < NRtol) {
              NRconvg = true;
              break;
            }
            else {
              Told  = Tcurr;
              Tcurr = Tnew;
            }
          }

          if(NRconvg)
            V_inlet[T_INDEX_PRIM] = Tcurr;
          else {
            /*--- Execute the bisection root-finding method ---*/
            Bconvg = false;
            su2double Ta = 300.0/config->GetTemperature_Ref();
            su2double Tb = Ttot;
            for(iIter = 0; iIter < maxBIter; ++iIter) {
              Tcurr = (Ta + Tb)/2.0;
              su2double F = f(Tcurr) - Tot_Enthalpy;

              if(std::abs(F) < Btol) {
                V_inlet[T_INDEX_PRIM] = Tcurr;
                Bconvg = true;
                break;
              }
              else {
                if(F > 0.0)
                  Ta = Tcurr;
                else
                  Tb = Tcurr;
              }
            }

            /*--- If absolutely no convergence, then something is going really wrong ---*/
            if(!Bconvg)
              throw std::runtime_error("Convergence not achieved for bisection method in inlet boundary condition");
          }

          /*--- Set inlet state ---*/
          su2double rho_tot = Ptot/(Rgas*Ttot);
          V_inlet[RHO_INDEX_PRIM] = rho_tot*std::pow(V_inlet[T_INDEX_PRIM]/Ttot, 1.0/Gamma_Minus_One);
          V_inlet[P_INDEX_PRIM] = V_inlet[RHO_INDEX_PRIM]*Rgas*V_inlet[T_INDEX_PRIM];
          V_inlet[A_INDEX_PRIM] = std::sqrt(V_inlet[T_INDEX_PRIM]*Gamma*Rgas);
          Vel_Mag = std::abs((Riemann - 2.0*V_inlet[A_INDEX_PRIM]/Gamma_Minus_One)/alpha);
          for(iDim = 0; iDim < nDim; ++iDim)
            V_inlet[VX_INDEX_PRIM + iDim] = Vel_Mag*Flow_Dir[iDim];
          std::copy(Ys.cbegin(), Ys.cend(), V_inlet + RHOS_INDEX_PRIM);

          break;
        }

        /*--- Mass flow has been specified at the inlet. ---*/
        case MASS_FLOW: {
          /*--- Local variables ---*/
          su2double Density, Pressure;

          /*--- Retrieve the specified mass flow for the inlet. ---*/
          Density = config->GetInlet_Ttotal(Marker_Tag);
          Vel_Mag = config->GetInlet_Ptotal(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/
          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/
          Pressure = V_domain[P_INDEX_PRIM];
          SoundSpeed = V_domain[A_INDEX_PRIM];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
                from the domain interior. ---*/
          Gamma = node[iPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0;
          Gamma_Minus_One = Gamma - 1.0;
          Vn = std::inner_product(V_domain + VX_INDEX_PRIM, V_domain + (VX_INDEX_PRIM + nDim), UnitNormal, 0.0);
          Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

          /*--- Speed of sound squared for fictitious inlet state ---*/
          alpha = std::inner_product(UnitNormal, UnitNormal + nDim, Flow_Dir, 0.0);
          SoundSpeed = Riemann - Vel_Mag*alpha;
          SoundSpeed = std::max(0.0, 0.5*Gamma_Minus_One*SoundSpeed);

          /*--- Pressure for the fictitious inlet state ---*/
          Pressure = SoundSpeed*SoundSpeed*Density/Gamma;

          /*--- Set primitive variables using the derived quantities ---*/
          V_inlet[T_INDEX_PRIM] = library->ComputeTemperature(Pressure, Density, Ys)*config->GetGas_Constant_Ref();
          for(iDim = 0; iDim < nDim; ++iDim)
            V_inlet[VX_INDEX_PRIM + iDim] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[P_INDEX_PRIM] = Pressure;
          V_inlet[RHO_INDEX_PRIM] = Density;
          dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          V_inlet[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
          if(US_System)
            V_inlet[H_INDEX_PRIM] *= 3.28084*3.28084;
          V_inlet[H_INDEX_PRIM] += 0.5*Vel_Mag*Vel_Mag;
          V_inlet[A_INDEX_PRIM] = SoundSpeed;
          std::copy(Ys.cbegin(), Ys.cend(), V_inlet + RHOS_INDEX_PRIM);

          break;
        }

        /*--- Mass flow has been specified at the inlet. ---*/
        case TEMPERATURE_IMPOSE: {
          /*--- Local variables ---*/
          su2double Temperature;

          /*--- Retrieve the specified mass flow for the inlet. ---*/
          Temperature = config->GetInlet_Ttotal(Marker_Tag);
          Vel_Mag = config->GetInlet_Ptotal(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/
          Temperature /= config->GetTemperature_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Set primitive variables using the derived quantities ---*/
          V_inlet[T_INDEX_PRIM] = Temperature;
          for(iDim = 0; iDim < nDim; ++iDim)
            V_inlet[VX_INDEX_PRIM + iDim] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[P_INDEX_PRIM] = V_domain[P_INDEX_PRIM];
          V_inlet[RHO_INDEX_PRIM] = library->ComputeDensity(V_inlet[T_INDEX_PRIM], V_inlet[P_INDEX_PRIM], Ys)*config->GetGas_Constant_Ref();
          dim_temp = Temperature*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          V_inlet[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
          if(US_System)
            V_inlet[H_INDEX_PRIM] *= 3.28084*3.28084;
          V_inlet[H_INDEX_PRIM] += 0.5*Vel_Mag*Vel_Mag;
          V_inlet[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp, Ys)/config->GetVelocity_Ref();
          if(US_System)
            V_inlet[A_INDEX_PRIM] *= 3.28084;
          std::copy(Ys.cbegin(), Ys.cend(), V_inlet + RHOS_INDEX_PRIM);

          break;
        }

      }

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);

      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      if(implicit) {
        su2double Secondary[nVar];

        /*--- Derivatives with respect to density ---*/
        su2double Velocity2 = Vel_Mag*Vel_Mag;
        Secondary[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*Velocity2;

        /*--- Derivatives with respect to momentum ---*/
        for(iDim = 0; iDim < nDim; ++iDim)
          Secondary[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V_inlet[VX_INDEX_PRIM + iDim];

        /*--- Derivatives with respect to energy ---*/
        Secondary[RHOE_INDEX_SOL] = Gamma - 1.0;

        /*--- Derivatives with respect to partial densities ---*/
        dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Secondary[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp, Gamma, iSpecies)/config->GetEnergy_Ref();
          if(US_System)
            Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
        }
        conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), Secondary);
      }

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      bool err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      if(!err) {
        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the convective residual of subsonic inlet boundary condition.");

      /*--- Viscous contribution ---*/
      if(viscous) {
        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables, and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

        /*--- Set temperature derivatives if implicit ---*/
        if(implicit) {
          su2double Secondary[nVar];

          su2double dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          su2double Cv = library->ComputeCV(dim_temp, Ys)/config->GetGas_Constant_Ref();
          if(US_System)
            Cv *= 3.28084*3.28084*5.0/9.0;
          su2double rhoCv = V_inlet[RHO_INDEX_PRIM]*Cv;
          su2double sq_vel = std::inner_product(V_inlet + VX_INDEX_PRIM, V_inlet + (VX_INDEX_PRIM + nDim), V_inlet + VX_INDEX_PRIM, 0.0);

          /*--- Set temperature derivatives ---*/
          Secondary[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
          for(unsigned short iDim = 0; iDim < nDim; ++iDim)
            Secondary[RHOVX_INDEX_SOL + iDim] = -V_inlet[VX_INDEX_PRIM + iDim]/rhoCv;
          Secondary[RHOE_INDEX_SOL] = 1.0/rhoCv;
          for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            Secondary[RHOS_INDEX_SOL + iSpecies] = -library->ComputePartialEnergy(dim_temp, iSpecies)/(config->GetEnergy_Ref()*rhoCv);
            if(US_System)
              Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
          }
          visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), Secondary);
        }

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[iPoint]->GetDiffusionCoeff());

        /*--- Compute residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        /*--- Check for NaNs before applying the residual to the linear system ---*/
        err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
        if(implicit) {
          if(!err) {
            for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
              err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
              if(err)
                break;
            }
          }
        }

        if(!err) {
          /*--- Update residual ---*/
          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/
          if(implicit)
            Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        else
          throw std::runtime_error("NaN found in the viscous residual of subsonic inlet boundary condition.");
      }

    } /*--- End of if GetDomain() ---*/

  } /*--- End of iVertex for loop ---*/
}

//
//
/*--- Impose supersonic outlet boundary condition. ---*/
//
//
void CReactiveEulerSolver::BC_Supersonic_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                                CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  /*--- Local variables ---*/
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Normal[nDim];

  std::string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool viscous = config->GetViscous();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      std::transform(Normal, Normal + nDim, Normal, std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/
      auto V_domain = node[iPoint]->GetPrimitive();

      /*--- Get characteristic variables on the boundary ---*/
      auto V_outlet = GetCharacPrimVar(val_marker, iVertex);

      /*--- Supersonic exit flow: there are no incoming characteristics,
            so no boundary condition is necessary. Set outlet state to current
            state so that upwinding handles the direction of propagation. ---*/
      std::copy(V_domain, V_domain + nPrimVar, V_outlet);

      /*--- Set primitives and grid movement in the convective numerics class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      if(implicit)
        conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), node[iPoint]->GetdPdU());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      bool err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      if(!err) {
        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the convective residual of supersonic outlet boundary condition.");

      /*--- Viscous contribution ---*/
      if(viscous) {
        /*--- Index of the closest interior node ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

        /*--- Set temperature derivatives if implicit ---*/
        if(implicit)
          visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), node[iPoint]->GetdTdU());

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[iPoint]->GetDiffusionCoeff());

        /*--- Compute viscous residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        /*--- Check for NaNs before applying the residual to the linear system ---*/
        err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
        if(implicit) {
          if(!err) {
            for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
              err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
              if(err)
                break;
            }
          }
        }

        if(!err) {
          /*--- Update residual ---*/
          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/
          if(implicit)
            Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
      }

    } /*--- End of if GetDomain() ---*/

  } /*--- End of iVertex for loop ---*/
}


//
//
/*--- Impose outlet boundary condition (supersonic or subosnic according to specific simulation). ---*/
//
//
void CReactiveEulerSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                                     CConfig* config, unsigned short val_marker) {
  /*--- Local variables ---*/
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Pressure, Velocity2, Density, SoundSpeed, Gamma, Mach_Exit, Area;
  su2double Velocity[nDim], Normal[nDim], UnitNormal[nDim];

  std::string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool viscous = config->GetViscous();
  bool gravity = config->GetGravityForce();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      Area = ::ComputeArea(Normal,nDim);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Normal[iDim] = -Normal[iDim];
        UnitNormal[iDim] = Normal[iDim]/Area;
      }
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/
      auto V_domain = node[iPoint]->GetPrimitive();

      /*--- Get characteristic variables on the boundary ---*/
      auto V_outlet = GetCharacPrimVar(val_marker, iVertex);

      /*--- Check whether the flow is supersonic at the exit. The type
            of boundary update depends on this. ---*/
      Density = V_domain[RHO_INDEX_PRIM];
      Velocity2 = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim) {
        Velocity[iDim] = V_domain[VX_INDEX_PRIM + iDim];
        Velocity2 += Velocity[iDim]*Velocity[iDim];
      }
      Pressure = V_domain[P_INDEX_PRIM];
      Gamma = node[iPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0;
      SoundSpeed = std::sqrt(Gamma*Pressure/Density);
      Mach_Exit = std::sqrt(Velocity2)/SoundSpeed;

      if(Mach_Exit >= 1.0) {
        /*--- Supersonic exit flow: there are no incoming characteristics,
              so no boundary condition is necessary. Set outlet state to current
              state so that upwinding handles the direction of propagation. ---*/
        std::copy(V_domain, V_domain + nPrimVar, V_outlet);
      }
      else {
        /*--- Subsonic exit flow: there is one incoming characteristic,
              therefore one variable can be specified (back pressure) and is used
              to update the conservative variables. Compute the entropy and the
              acoustic Riemann variable. These invariants, as well as the
              tangential velocity components, are extrapolated. ---*/
        su2double Vn, Vn_Exit, Entropy, Riemann, Gamma_Minus_One, Temperature, dim_temp;

        /*--- Compute entropy and Riemann invariant ---*/
        Entropy = Pressure*std::pow(1.0/Density, Gamma);
        Vn = std::inner_product(Velocity, Velocity + nDim, UnitNormal, 0.0);
        Gamma_Minus_One = Gamma - 1.0;
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Build the fictitious intlet state based on characteristics ---*/
        /*--- Retrieve the specified back pressure for this outlet. ---*/
        su2double P_Exit = config->GetOutlet_Pressure(Marker_Tag);
        if(gravity)
          P_Exit -= geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;

        /*--- Non-dim. the inputs if necessary. ---*/
        P_Exit /= config->GetPressure_Ref();

        /*--- Compute the new fictitious state ---*/
        V_outlet[P_INDEX_PRIM] = P_Exit;
        V_outlet[RHO_INDEX_PRIM] = std::pow(P_Exit/Entropy, 1.0/Gamma);
        SoundSpeed = std::sqrt(Gamma*P_Exit/V_outlet[RHO_INDEX_PRIM]);
        Vn_Exit = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        Velocity2  = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim) {
          Velocity[iDim] += (Vn_Exit - Vn)*UnitNormal[iDim];
          Velocity2 += Velocity[iDim]*Velocity[iDim];
          V_outlet[VX_INDEX_PRIM + iDim] = Velocity[iDim];
        }
        std::copy(V_domain + RHOS_INDEX_PRIM, V_domain + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
        V_outlet[T_INDEX_PRIM] = library->ComputeTemperature(P_Exit, V_outlet[RHO_INDEX_PRIM], Ys)*config->GetGas_Constant_Ref();
        Temperature = V_outlet[T_INDEX_PRIM];
        dim_temp = Temperature*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        V_outlet[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
        if(US_System)
          V_outlet[H_INDEX_PRIM] *= 3.28084*3.28084;
        V_outlet[H_INDEX_PRIM] += 0.5*Velocity2;
        V_outlet[A_INDEX_PRIM] = SoundSpeed;
        std::copy(Ys.cbegin(), Ys.cend(), V_outlet + RHOS_INDEX_PRIM);
      } /*--- End subsonic outlet preliminary computations ---*/

      /*--- Set primitives and grid movement in the convective numerics class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      if(implicit) {
        if(Mach_Exit >= 1.0)
          conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), node[iPoint]->GetdPdU());
        else {
          su2double Secondary[nVar];

          /*--- Derivatves with respect to density ---*/
          Secondary[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*Velocity2;

          /*--- Derivatives with respect to momentum ---*/
          for(iDim = 0; iDim < nDim; ++iDim)
            Secondary[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V_outlet[VX_INDEX_PRIM + iDim];

          /*--- Derivatives with respect to energy ---*/
          Secondary[RHOE_INDEX_SOL] = Gamma - 1.0;

          /*--- Derivatives with respect to partial densities ---*/
          su2double dim_temp = V_outlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
          if(US_System)
            dim_temp *= 5.0/9.0;
          for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            Secondary[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp, Gamma, iSpecies)/config->GetEnergy_Ref();
            if(US_System)
              Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
          }
          conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), Secondary);
        }
      }

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      bool err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      if(!err) {
        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the convective residual of outlet boundary condition.");

      /*--- Viscous contribution ---*/
      if(viscous) {
        /*--- Index of the closest interior node ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Set the normal vector and the coordinates ---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables and gradient ---*/
        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

        /*--- Set temperature derivatives if implicit ---*/
        if(implicit) {
          if(Mach_Exit >= 1.0)
            visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), node[iPoint]->GetdTdU());
          else {
            su2double Secondary[nVar];

            su2double dim_temp = V_outlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
            if(US_System)
              dim_temp *= 5.0/9.0;
            su2double Cv = library->ComputeCV(dim_temp, Ys)/config->GetGas_Constant_Ref();
            if(US_System)
              Cv *= 3.28084*3.28084*5.0/9.0;
            su2double rhoCv = V_outlet[RHO_INDEX_PRIM]*Cv;
            su2double sq_vel = std::inner_product(V_outlet + VX_INDEX_PRIM, V_outlet + (VX_INDEX_PRIM + nDim),
                                                  V_outlet + VX_INDEX_PRIM, 0.0);

            /*--- Set temperature derivatives ---*/
            Secondary[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
            for(unsigned short iDim = 0; iDim < nDim; ++iDim)
              Secondary[RHOVX_INDEX_SOL + iDim] = -V_outlet[VX_INDEX_PRIM + iDim]/rhoCv;
            Secondary[RHOE_INDEX_SOL] = 1.0/rhoCv;
            for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
              Secondary[RHOS_INDEX_SOL + iSpecies] = -library->ComputePartialEnergy(dim_temp, iSpecies)/(config->GetEnergy_Ref()*rhoCv);
              if(US_System)
                Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
            }
            visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), Secondary);
          }
        } /*--- End of implicit part computations ---*/

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[iPoint]->GetDiffusionCoeff());

        /*--- Compute residual ---*/
        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

        /*--- Check for NaNs before applying the residual to the linear system ---*/
        err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
        if(implicit) {
          if(!err) {
            for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
              err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
              if(err)
                break;
            }
          }
        }

        if(!err) {
          /*--- Update residual ---*/
          LinSysRes.SubtractBlock(iPoint, Residual);

          /*--- Jacobian contribution for implicit integration ---*/
          if(implicit)
            Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }
        else
          throw std::runtime_error("NaN found in the viscous residual of outlet boundary condition.");
      }

    } /*--- End of if GetDomain() ---*/

  } /*--- End of iVertex for loop ---*/
}


//
//
/*--- Class CReactiveNSSolver constructor ---*/
//
//
CReactiveNSSolver::CReactiveNSSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh): CReactiveEulerSolver() {
  /*--- Local variables ---*/
  unsigned long iPoint;
  unsigned short iVar, iDim;

  /*--- Initialize built-in types often used ---*/
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  IterLinSolver = 0;

  /*--- Load library ---*/
  if(iMesh == MESH_0) {
    library = LibraryPtr(new Framework::ReactingModelLibrary(config->GetConfigLibFile(), config->GetLibraryPath()));
    library->Setup();
  }
  nSpecies = library->GetnSpecies();

  /*--- Set other variables ---*/
  bool restart = (config->GetRestart() || config->GetRestart_Flow());

  /*--- Read and store the restart metadata. ---*/
  if(!(!restart || iMesh != MESH_0))
    Read_SU2_Restart_Metadata(geometry, config);

  /*--- Continue setting other variables ---*/
  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho, rho*vx, rho*vy, rho*vz, rho*E, rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T, vx, vy, vz, P, rho, h, a, Y1,...YNs)^T ---*/
  nPrimVarLim = nDim + 2;  /*--- Primitive variables to limit (T, vx, vy, vz, P)^T ---*/
  nPrimVarGrad = nPrimVarLim + nSpecies; /*--- Gradient Primitive variables (T, vx, vy, vz, P, X1....XNs)^T ---*/

  /*--- Read freestream mass fractions. ---*/
  MassFrac_Inf = RealVec(config->GetMassFrac_FreeStream(), config->GetMassFrac_FreeStream() + nSpecies);

  /*--- Check right order of species in configuration file ----*/
  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
  if(iMesh == MESH_0 && rank == MASTER_NODE)
    Check_FreeStream_Species_Order(config);

  /*--- Perform the non-dimensionalization for the flow equations using the specified reference values. ---*/
  US_System = (config->GetSystemMeasurements() == US);
  SetNondimensionalization(geometry, config, iMesh);

	/*--- Allocate a CVariable array foreach node of the mesh ---*/
	node = new CVariable*[nPoint];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new su2double[nVar];
  Residual_RMS = new su2double[nVar];
  Residual_Max = new su2double[nVar];
  Residual_i = new su2double[nVar];
  Residual_j = new su2double[nVar];
  Res_Conv = new su2double[nVar];
  Res_Sour = new su2double[nVar];
  Res_Visc = new su2double[nVar];
  Vector_i = new su2double[nDim];
  Vector_j = new su2double[nDim];
  Vector   = new su2double[nDim];

  /*--- Set residual vectors to zero for safety ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    Residual[iVar] = 0.0;
    Residual_RMS[iVar] = 0.0;
    Residual_Max[iVar] = 0.0;
    Residual_i[iVar] = 0.0;
    Residual_j[iVar] = 0.0;
    Res_Conv[iVar] = 0.0;
    Res_Sour[iVar] = 0.0;
  }

  /*--- Define some structures for locating max residuals ---*/
  Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
  for(iVar = 0; iVar < nVar; ++iVar) {
    Point_Max[iVar] = 0;
    Point_Max_Coord[iVar] = new su2double[nDim];
    std::fill(Point_Max_Coord[iVar], Point_Max_Coord[iVar] + nDim, 0.0);
  }

  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];
  Solution_i = new su2double[nVar];
  Solution_j = new su2double[nVar];

  for(iVar = 0; iVar < nVar; ++iVar) {
    Solution[iVar] = 0.0;
    Solution_i[iVar] = 0.0;
    Solution_j[iVar] = 0.0;
  }

  /*--- Allocate vectors related to the primitive variables ---*/
  PrimVar_i.resize(nPrimVarGrad);
  PrimVar_j.resize(nPrimVarGrad);
  PrimVar_Vertex.resize(nPrimVarGrad);

  Prim_i.resize(nPrimVarLim);
  Prim_j.resize(nPrimVarLim);
  Primitive.resize(nPrimVarLim);

  /*--- Allocate auxiliary vectors ---*/
  Ys_i.resize(nSpecies);
  Ys_j.resize(nSpecies);
  Ys.resize(nSpecies);

  /*--- Initialize the solution & residual CVectors ---*/
 	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
 	if(config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

	/*--- Allocate Jacobians for implicit time-stepping ---*/
	if(config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    implicit = true;
		Jacobian_i = new su2double* [nVar];
		Jacobian_j = new su2double* [nVar];
		for(iVar = 0; iVar < nVar; ++iVar) {
			Jacobian_i[iVar] = new su2double[nVar];
			Jacobian_j[iVar] = new su2double[nVar];
		}

    if(rank == MASTER_NODE)
      std::cout << "Initialize Jacobian structure (Navier-Stokes). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  }
  else
    implicit = false;

  /*--- Allocate arrays for gradient computation by least squares ---*/
	if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    least_squares = true;

    Smatrix = new su2double* [nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      Smatrix[iDim] = new su2double [nDim];

    Cvector = new su2double* [nPrimVarGrad];
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      Cvector[iVar] = new su2double[nDim];
  }
	else
    least_squares = false;

  /*--- Store the value of the characteristic primitive variables at the boundaries ---*/
  nMarker = config->GetnMarker_All();
  CharacPrimVar = new su2double** [nMarker];
  unsigned short iMarker;
  unsigned long iVertex;
  for(iMarker = 0; iMarker < nMarker; ++iMarker) {
    CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
    for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
      CharacPrimVar[iMarker][iVertex] = new su2double[nPrimVar];
      for(iVar = 0; iVar < nPrimVar; ++iVar) {
        CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
      }
    }
  }

  /*--- Store the number of vertices on each marker for deallocation later ---*/
  nVertex.reserve(nMarker);
  for(iMarker = 0; iMarker < nMarker; ++iMarker)
    nVertex.push_back(geometry->nVertex[iMarker]);

  /*--- Set other variables ---*/
  grid_movement = config->GetGrid_Movement();
  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  second_order = (config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  if(second_order) {
    Primitive_i.resize(nPrimVar);
    Primitive_j.resize(nPrimVar);
    if(implicit) {
      Secondary_i.resize(nVar);
      Secondary_j.resize(nVar);
    }
  }
  unsigned long ExtIter = config->GetExtIter();
  limiter = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter());

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(), config->GetVelocity_FreeStreamND() + nDim);

  Viscosity_Inf    = config->GetViscosity_FreeStreamND();

  /*--- Initialize the far field state ---*/
  node_infty = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, Viscosity_Inf,
                                       nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);

  /*--- Set local variables to access indices in proper arrays ---*/
  T_INDEX_PRIM    = CReactiveNSVariable::GetT_INDEX_PRIM();
  VX_INDEX_PRIM   = CReactiveNSVariable::GetVX_INDEX_PRIM();
  P_INDEX_PRIM    = CReactiveNSVariable::GetP_INDEX_PRIM();
  RHO_INDEX_PRIM  = CReactiveNSVariable::GetRHO_INDEX_PRIM();
  H_INDEX_PRIM    = CReactiveNSVariable::GetH_INDEX_PRIM();
  A_INDEX_PRIM    = CReactiveNSVariable::GetA_INDEX_PRIM();
  RHOS_INDEX_PRIM = CReactiveNSVariable::GetRHOS_INDEX_PRIM();

  RHO_INDEX_SOL   = CReactiveNSVariable::GetRHO_INDEX_SOL();
  RHOVX_INDEX_SOL = CReactiveNSVariable::GetRHOVX_INDEX_SOL();
  RHOE_INDEX_SOL  = CReactiveNSVariable::GetRHOE_INDEX_SOL();
  RHOS_INDEX_SOL  = CReactiveNSVariable::GetRHOS_INDEX_SOL();

  T_INDEX_GRAD    = CReactiveNSVariable::GetT_INDEX_GRAD();
  VX_INDEX_GRAD   = CReactiveNSVariable::GetVX_INDEX_GRAD();
  P_INDEX_GRAD    = CReactiveNSVariable::GetP_INDEX_GRAD();
  RHOS_INDEX_GRAD = CReactiveNSVariable::GetRHOS_INDEX_GRAD();

  T_INDEX_LIM     = CReactiveNSVariable::GetT_INDEX_LIM();
  VX_INDEX_LIM    = CReactiveNSVariable::GetVX_INDEX_LIM();
  P_INDEX_LIM     = CReactiveNSVariable::GetP_INDEX_LIM();

  if(!restart || iMesh != MESH_0) {
    /*--- Initialize the solution to the far-field state everywhere. ---*/
    for(iPoint = 0; iPoint < nPoint; ++iPoint)
      node[iPoint] = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, Viscosity_Inf,
                                             nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);
  }
  else
    Load_Restart(geometry, config);

  /*--- Allocate vectors for conserved variable limits ---*/
  Lower_Limit.resize(nVar, 0.0);
  Upper_Limit.resize(nVar, 1.0/EPS);

  std::fill(Lower_Limit.begin() + RHOVX_INDEX_SOL, Lower_Limit.begin() + (RHOVX_INDEX_SOL + nDim), -1.0/EPS);
  Lower_Limit[RHOE_INDEX_SOL] = -1.0/EPS;

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
}

//
//
/*--- Read the file to restart. ---*/
//
//
void CReactiveNSSolver::Load_Restart(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  std::string filename = config->GetSolution_FlowFileName();
  int Unst_RestartIter;

  /*--- Multizone problems require the number of the zone to be appended. ---*/
  if(nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Modify file name for a dual-time unsteady restart ---*/
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if(dual_time) {
    if(config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    else
      Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 2;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Modify file name for a simple unsteady restart ---*/
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  if (time_stepping) {
    Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter()) - 1;
    filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
  }

  /*--- Open the restart file, throw an error if this fails. ---*/
  int rank = MASTER_NODE;
  std::ifstream restart_file;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
  restart_file.open(filename.data(), std::ios::in);
  if(restart_file.fail()) {
    if(rank == MASTER_NODE)
      std::cout<<"There is no flow restart file!! "<<filename.data()<<"."<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*--- In case this is a parallel simulation, we need to perform the
        Global2Local index transformation first. ---*/
  std::map<unsigned long, unsigned long> Global2Local;
  std::map<unsigned long, unsigned long>::const_iterator MI;

  /*--- Now fill array with the transform values only for local points ---*/
  unsigned long iPoint;
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

  /*--- Read all lines in the restart file ---*/
  unsigned long iPoint_Local, iPoint_Global;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  std::string text_line;

  /*--- The first line is the header ---*/
  std::getline(restart_file, text_line);
  unsigned long index;
  su2double dull_val;
  for(iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); ++iPoint_Global) {
    std::getline(restart_file, text_line);
    std::istringstream point_line(text_line);
    if(iPoint_Global >= geometry->GetGlobal_nPointDomain()) {
      sbuf_NotMatching = 1;
      break;
    }

    /*--- Retrieve local index. If this node from the restart file lives
          on the current processor, we will load and instantiate the vars. ---*/
    MI = Global2Local.find(iPoint_Global);
    if(MI != Global2Local.end()) {
      iPoint_Local = Global2Local[iPoint_Global];
      if(nDim == 2)
        point_line >> index >> dull_val >> dull_val >> Solution[RHO_INDEX_SOL] >> Solution[RHOVX_INDEX_SOL] >>
                      Solution[RHOVX_INDEX_SOL + 1] >> Solution[RHOE_INDEX_SOL];
      else if(nDim == 3)
        point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[RHO_INDEX_SOL] >> Solution[RHOVX_INDEX_SOL] >>
                      Solution[RHOVX_INDEX_SOL + 1] >> Solution[RHOVX_INDEX_SOL + 2] >> Solution[RHOE_INDEX_SOL];
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        point_line >> Solution[RHOS_INDEX_SOL + iSpecies];

      node[iPoint_Local] = new CReactiveNSVariable(Solution, nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);
      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/
  if(iPoint_Global_Local < nPointDomain)
    sbuf_NotMatching = 1;

  #ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
  #else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
  #endif
  if(rbuf_NotMatching != 0) {
    if(rank == MASTER_NODE) {
      std::cout<<std::endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }

    #ifndef HAVE_MPI
      std::exit(EXIT_FAILURE);
    #else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    #endif
  }

  /*--- Instantiate the variable class with an arbitrary solution
        at any halo/periodic nodes. The initial solution can be arbitrary,
        because a send/recv is performed immediately in the solver. ---*/
  for(iPoint = nPointDomain; iPoint < nPoint; ++iPoint)
    node[iPoint] = new CReactiveNSVariable(Solution, nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, library, config);

  /*--- Close the restart file ---*/
  restart_file.close();
}

//
//
/*--- Set the fluid solver nondimensionalization. ---*/
//
//
void CReactiveNSSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {
  /*--- Adimensionalitazion for the invisicid part does not change and so we use the Euler version. ---*/
  CReactiveEulerSolver::SetNondimensionalization(geometry, config, iMesh);

  /*--- Get the rank of an eventual parallel simulation ---*/
  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  /*--- Local variables for the adimensionalitazion of viscous contribution ---*/
  su2double Viscosity_FreeStream = 0.0;
  su2double Viscosity_Ref = 0.0, Conductivity_Ref = 0.0;
  su2double Viscosity_FreeStreamND = 0.0;

  Viscosity_Ref = config->GetDensity_Ref()*config->GetVelocity_Ref()*config->GetLength_Ref();
  config->SetViscosity_Ref(Viscosity_Ref);

  Conductivity_Ref = Viscosity_Ref*config->GetGas_Constant_Ref();
  config->SetConductivity_Ref(Conductivity_Ref);

  su2double dim_temp = config->GetTemperature_FreeStream();

  bool SI_Measurement = (config->GetSystemMeasurements() == SI);
  bool US_Measurament = US_System;

  if(US_Measurament)
    dim_temp *= 5.0/9.0;
  Viscosity_FreeStream = library->ComputeEta(dim_temp, MassFrac_Inf);
  if(US_Measurament)
    Viscosity_FreeStream *= 0.02088553108;
  config->SetViscosity_FreeStream(Viscosity_FreeStream);

  Viscosity_FreeStreamND = Viscosity_FreeStream/Viscosity_Ref;
  config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  /*--- Write output to the console if this is the master node and first domain ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH && rank == MASTER_NODE && iMesh == MESH_0) {
    std::cout.precision(6);

    std::cout<< "Reference viscosity: " << config->GetViscosity_Ref();
    if(SI_Measurement)
      std::cout << " N.s/m^2."<< std::endl;
    else if(US_Measurament)
      std::cout<< " lbf.s/ft^2."<< std::endl;

    std::cout << "Reference conductivity: " << config->GetConductivity_Ref();
    if(SI_Measurement)
      std::cout << " W/m.K."<< std::endl;
    else if(US_Measurament)
      std::cout << " lbf/s.R."<< std::endl;

    config->SetReynolds(config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetLength_Reynolds()/Viscosity_FreeStream);
    std::cout<< "Reynolds number (non-dim): " << config->GetReynolds() <<std::endl;

    std::cout<< "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << std::endl;

    std::cout<<std::endl;
  }
}

//
//
/*--- Preprocessing stage to set primitive variables and gradient of primitive variables. ---*/
//
//
void CReactiveNSSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                                      unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  /*--- Set the primitive variables ---*/
  unsigned long ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Gradient computation ---*/
  if(least_squares)
    SetPrimitive_Gradient_LS(geometry, config);
  else if(config->GetKind_Gradient_Method() == GREEN_GAUSS)
    SetPrimitive_Gradient_GG(geometry, config);

  /*--- Limiter computation ---*/
  if(limiter && iMesh == MESH_0 && !Output)
    SetPrimitive_Limiter(geometry, config);

  /*--- Artificial dissipation ---*/
  if(space_centered && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");

  /*--- Initialize the Jacobian matrices ---*/
  if(implicit)
    Jacobian.SetValZero();

  /*--- Error message ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      SU2_MPI::Allreduce(MPI_IN_PLACE, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    #endif
    if(iMesh == MESH_0)
      config->SetNonphysical_Points(ErrorCounter);
  }
}

//
//
/*--- Set gradient of primitive variables according to Green-Gauss. ---*/
//
//
void CReactiveNSSolver::SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iSpecies, iVar, iMarker;

  su2double PrimVar_Average, Partial_Res;
  su2double Volume;
  su2double Normal[nDim];

  /*--- Set Gradient_Primitive to zero ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
  	node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
  	iPoint = geometry->edge[iEdge]->GetNode(0);
  	jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
    PrimVar_i[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_j[T_INDEX_GRAD] = node[jPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_i[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
    PrimVar_j[P_INDEX_GRAD] = node[jPoint]->GetPrimitive(P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim) {
      PrimVar_i[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
      PrimVar_j[VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
    }
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      PrimVar_i[RHOS_INDEX_GRAD + iSpecies] = node[iPoint]->GetPrimitive(RHOS_INDEX_PRIM + iSpecies);
      PrimVar_j[RHOS_INDEX_GRAD + iSpecies] = node[iPoint]->GetPrimitive(RHOS_INDEX_PRIM + iSpecies);
    }

    /*--- Compute molar fractions for iPoint and jPoint ---*/
    Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.cbegin() + RHOS_INDEX_GRAD, PrimVar_i.cbegin() + (RHOS_INDEX_GRAD + nSpecies)));
    Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.cbegin() + RHOS_INDEX_GRAD, PrimVar_j.cbegin() + (RHOS_INDEX_GRAD + nSpecies)));
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      PrimVar_i[RHOS_INDEX_GRAD + iSpecies] = Xs_i[iSpecies];
      PrimVar_j[RHOS_INDEX_GRAD + iSpecies] = Xs_j[iSpecies];
    }

  	geometry->edge[iEdge]->GetNormal(Normal);

    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      /*--- Compute face-averaged value as discrete counter-part of Green-Gauss theorem. ---*/
      PrimVar_Average = 0.5*(PrimVar_i[iVar] + PrimVar_j[iVar]);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Partial_Res = PrimVar_Average*Normal[iDim];
  		  if(geometry->node[iPoint]->GetDomain())
  			  node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
        if(geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  	  }
    }
  } /*--- End loop interior edges ---*/

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
  	for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
  		iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  		if(geometry->node[iPoint]->GetDomain()) {
        /*--- Get primitives from CVariable ---*/
        PrimVar_Vertex[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
        PrimVar_Vertex[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
        for(iDim = 0; iDim < nDim; ++iDim)
          PrimVar_Vertex[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          PrimVar_Vertex[RHOS_INDEX_GRAD + iSpecies] = node[iPoint]->GetPrimitive(RHOS_INDEX_PRIM + iSpecies);

        /*--- Compute molar fractions ---*/
        Xs = library->GetMolarFromMass(RealVec(PrimVar_Vertex.cbegin() + RHOS_INDEX_GRAD,
                                               PrimVar_Vertex.cbegin() + (RHOS_INDEX_GRAD + nSpecies)));
        std::copy(Xs.cbegin(), Xs.cend(), PrimVar_Vertex.begin() + RHOS_INDEX_GRAD);

        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
        for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
          for(iDim = 0; iDim < nDim; ++iDim) {
            /*--- NOTE: On the boundary we cannot compute face-averaged value and so we use directly the stored value. ---*/
            Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
  				  node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
        	}
        }
  		}
  	} /*--- End of iVertex for loop ---*/
  } /*--- End of iMarker for loop ---*/

  /*--- Update gradient value ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();
    SU2_Assert(Volume > EPS,"The measure of the volume is not consistent(smaller or equal to zero)");
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  	  for(iDim = 0; iDim < nDim; ++iDim)
  		  node[iPoint]->SetGradient_Primitive(iVar, iDim, node[iPoint]->GetGradient_Primitive(iVar,iDim)/Volume);
    }
  }

  Set_MPI_Primitive_Gradient(geometry, config);
}

//
//
/*--- Set gradient of primitive variables according to least squares. ---*/
//
//
void CReactiveNSSolver::SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned short iVar, iDim, jDim, iSpecies, iNeigh;
	unsigned long iPoint, jPoint;

  su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, detR2, z11, z12, z13, z22, z23, z33;
  su2double weight;
  su2double Coord_ij[nDim];

  bool singular;

  /*--- Gradient of the following primitive variables: [T, vx, vy, vz, p, X_1,...X_Ns]^T ---*/
  /*--- Loop over points of the grid ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		/*--- Set the value of singular ---*/
		singular = false;

    /*--- Get coordinates ---*/
		auto Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get primitives from CVariable ---*/
    PrimVar_i[T_INDEX_GRAD] = node[iPoint]->GetPrimitive(T_INDEX_PRIM);
    PrimVar_i[P_INDEX_GRAD] = node[iPoint]->GetPrimitive(P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim)
      PrimVar_i[VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      PrimVar_i[RHOS_INDEX_GRAD + iSpecies] = node[iPoint]->GetPrimitive(RHOS_INDEX_PRIM + iSpecies);

    /*--- Compute molar fractions from mass fractions for iPoint ---*/
    Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.cbegin() + RHOS_INDEX_GRAD, PrimVar_i.cbegin() + (RHOS_INDEX_GRAD + nSpecies)));
    std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.begin() + RHOS_INDEX_GRAD);

    /*--- Inizialization of variables ---*/
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      std::fill(Cvector[iVar], Cvector[iVar] + nDim, 0.0);

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      /*--- Get coordinates ---*/
			auto Coord_j = geometry->node[jPoint]->GetCoord();

      /*--- Get primitives from CVariable ---*/
      PrimVar_j[T_INDEX_GRAD] = node[jPoint]->GetPrimitive(T_INDEX_PRIM);
      PrimVar_j[P_INDEX_GRAD] = node[jPoint]->GetPrimitive(P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        PrimVar_j[VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(VX_INDEX_PRIM + iDim);
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        PrimVar_j[RHOS_INDEX_GRAD + iSpecies] = node[jPoint]->GetPrimitive(RHOS_INDEX_PRIM + iSpecies);

      /*--- Compute molar fractions from mass fractions for jPoint ---*/
      Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.cbegin() + RHOS_INDEX_GRAD, PrimVar_j.cbegin() + (RHOS_INDEX_GRAD + nSpecies)));
      std::copy(Xs_j.cbegin(), Xs_j.cend(), PrimVar_j.begin() + RHOS_INDEX_GRAD);

      /*--- Compute components of the vector that links i to j ---*/
      std::transform(Coord_j, Coord_j + nDim, Coord_i, Coord_ij, std::minus<su2double>());

      weight = std::inner_product(Coord_ij, Coord_ij + nDim, Coord_ij, 0.0);

			/*--- Sumations for entries of upper triangular matrix R ---*/
      if(weight > EPS){
        r11 += Coord_ij[0]*Coord_ij[0]/weight;
        r12 += Coord_ij[0]*Coord_ij[1]/weight;
        r22 += Coord_ij[1]*Coord_ij[1]/weight;
        if(nDim == 3) {
          r13 += Coord_ij[0]*Coord_ij[2]/weight;
          r23_a += Coord_ij[1]*Coord_ij[2]/weight;
          r23_b += Coord_ij[0]*Coord_ij[2]/weight;
          r33 += Coord_ij[2]*Coord_ij[2]/weight;
        }
        /*--- Entries of c:= transpose(A)*b ---*/
        for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
          for(iDim = 0; iDim < nDim; ++iDim)
            Cvector[iVar][iDim] += Coord_ij[iDim]*(PrimVar_j[iVar] - PrimVar_i[iVar])/weight;
      }
    } /*--- End iNeigh for loop ---*/

		/*--- Entries of upper triangular matrix R ---*/
    if(r11 > EPS)
      r11 = std::sqrt(r11);
    else
      r11 = 0.0;

    if(std::abs(r11) > EPS)
      r12 = r12/r11;
    else
      r12 = 0.0;

    if(r22 - r12*r12 > EPS)
      r22 = std::sqrt(r22 - r12*r12);
    else
      r22 = 0.0;

    if(nDim == 3) {
      if(std::abs(r11) > EPS)
        r13 = r13/r11;
      else
        r13 = 0.0;

      if(std::abs(r22) > EPS && std::abs(r11*r22) > EPS)
        r23 = r23_a/r22 - r23_b*r12/(r11*r22);
      else
        r23 = 0.0;
      if(r33 - r23*r23 - r13*r13 > EPS)
        r33 = std::sqrt(r33 - r23*r23 - r13*r13);
      else r33 = 0.0;
    }

    /*--- Compute determinant ---*/
    if(nDim == 2)
      detR2 = (r11*r22)*(r11*r22);
    else
      detR2 = (r11*r22*r33)*(r11*r22*r33);

    /*--- Detect singular matrices ---*/
    if(std::abs(detR2) < EPS) {
      detR2 = 1.0;
      singular = true;
    }

		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    if(singular) {
      for(iDim = 0; iDim < nDim; ++iDim)
        std::fill(Smatrix[iDim], Smatrix[iDim] + nDim, 0.0);
    }
    else {
      if(nDim == 2) {
        Smatrix[0][0] = (r12*r12 + r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11 + z12*z12 + z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22 + z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22 + z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }

	  /*--- Computation of the gradient: C*S ---*/
    su2double result;
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
		  for(iDim = 0; iDim < nDim; ++iDim) {
        result = 0.0;
        for(jDim = 0; jDim < nDim; ++jDim)
          result += Cvector[iVar][jDim]*Smatrix[iDim][jDim];
        node[iPoint]->SetGradient_Primitive(iVar, iDim, result);
      }
    }

	} /*--- End of iPoint for loop ---*/

  Set_MPI_Primitive_Gradient(geometry,config);
}

//
//
/*--- Set time step using maximum inviscid and viscous eigenvalues. ---*/
//
//
void CReactiveNSSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                     unsigned short iMesh, unsigned long Iteration) {
  /*--- Local variables ---*/
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker, iSpecies;

  su2double Area, Volume;
  su2double Local_Delta_Time, Global_Delta_Time = 1.0e6;
  su2double Local_Delta_Time_Visc;
  su2double Mean_SoundSpeed, Mean_ProjVel;
  su2double Mean_LaminarVisc, Mean_ThermalCond, Mean_Density, Mean_CV;
  su2double Lambda, Lambda_1, Lambda_2;
  su2double Normal[nDim];

  const su2double K_v = 0.25;

  bool time_steping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  Min_Delta_Time = 1.0e6;
  Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid and viscous eigenvalue to zero before computing it ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    node[iPoint]->SetMax_Lambda_Inv(0.0);
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    geometry->edge[iEdge]->GetNormal(Normal);
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5*(node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5*(node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());

    /*--- Adjustment for grid movement ---*/
    if(grid_movement) {
      auto GridVel_i = geometry->node[iPoint]->GetGridVel();
      auto GridVel_j = geometry->node[jPoint]->GetGridVel();
      su2double ProjVel_i = std::inner_product(GridVel_i, GridVel_i + nDim, Normal, 0.0);
      su2double ProjVel_j = std::inner_product(GridVel_j, GridVel_j + nDim, Normal, 0.0);
      Mean_ProjVel -= 0.5*(ProjVel_i + ProjVel_j);
    }

    Mean_Density = 0.5*(node[iPoint]->GetDensity() + node[jPoint]->GetDensity());
    Mean_ThermalCond = 0.5*(node[iPoint]->GetThermalConductivity() + node[jPoint]->GetThermalConductivity());
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_CV = 0.5*(node[iPoint]->GetSpecificHeatCp()/(node[iPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0) +
                   node[jPoint]->GetSpecificHeatCp()/(node[jPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0));

    /*--- Inviscid contribution ---*/
    Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
    if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);

    /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
  	Lambda_1 = 4.0/3.0*Mean_LaminarVisc;
  	Lambda_2 = Mean_ThermalCond/Mean_CV;
  	Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
  	if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);
  	if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Visc(Lambda);
  } /*--- End loop interior edges ---*/

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed();

        /*--- Adjustment for grid movement ---*/
        if(grid_movement) {
          auto GridVel = geometry->node[iPoint]->GetGridVel();
          su2double ProjVel = std::inner_product(GridVel, GridVel + nDim, Normal, 0.0);
          Mean_ProjVel -= ProjVel;
        }

        Mean_Density = node[iPoint]->GetDensity();
        Mean_ThermalCond = node[iPoint]->GetThermalConductivity();
        Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
        Mean_CV = node[iPoint]->GetSpecificHeatCp()/(node[iPoint]->GetdPdU()[RHOE_INDEX_SOL] + 1.0);

        /*--- Inviscid contribution ---*/
        Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
        if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Inv(Lambda);

        /*--- Viscous contribution ---*/
      	Lambda_1 = 4.0/3.0*Mean_LaminarVisc;
      	Lambda_2 = Mean_ThermalCond/Mean_CV;
      	Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      	if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Visc(Lambda);
      }
    }
  } /*--- End loop boundary edges ---*/

  /*--- Each element uses their own speed for a steady state simulation ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();

    if(Volume > EPS) {
      Local_Delta_Time = config->GetCFL(iMesh)*Volume/node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Volume*Volume/node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = std::min(Local_Delta_Time, Local_Delta_Time_Visc);
      Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
      if(Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      su2double rbuf_time, sbuf_time;
      sbuf_time = Min_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Min_Delta_Time = rbuf_time;

      sbuf_time = Max_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Max_Delta_Time = rbuf_time;
    #endif
  }

  /*--- Check if there is any element with only one neighbor...
        a CV that is inside another CV ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		if(geometry->node[iPoint]->GetnPoint() == 1)
			node[iPoint]->SetDelta_Time(Min_Delta_Time);
	}

  if(time_steping) {
    #ifdef HAVE_MPI
      su2double rbuf_time, sbuf_time;
      sbuf_time = Global_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Global_Delta_Time = rbuf_time;
    #endif
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  } /*--- End of time stepping checking ---*/

  if(dual_time) {
    su2double Global_Delta_UnstTimeND;
    if(Iteration == 0 && config->GetUnst_CFL() != 0.0 && iMesh == MESH_0) {
      Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

      #ifdef HAVE_MPI
        su2double rbuf_time, sbuf_time;
        sbuf_time = Global_Delta_UnstTimeND;
        SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
        SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
        Global_Delta_UnstTimeND = rbuf_time;
      #endif

      config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
    }

    /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
    for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
      if(!implicit) {
        Local_Delta_Time = std::min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
    }
  } /*--- End of dual-time check ---*/
}

//
//
/*--- Compute viscous residual at each interface. ---*/
//
//
void CReactiveNSSolver::Viscous_Residual(CGeometry* geometry, CSolver** solution_container, CNumerics* numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  /*--- Local variables ---*/
	unsigned long iPoint, jPoint, iEdge;

  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Set Primitive variables and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());

    /*--- Set temperature derivatives if implicit ---*/
    if(implicit)
      numerics->SetSecondary(node[iPoint]->GetdTdU(), node[jPoint]->GetdTdU());

    /*--- Set primitve variables limited ---*/
    if(limiter)
      numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());

    /*--- Set laminar viscosity at node i and j ---*/
    numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());

    /*--- Set thermal conductivity at node i and j ---*/
    numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[jPoint]->GetThermalConductivity());

    /*--- Set species binary diffusion coefficients at node i and j ---*/
    numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[jPoint]->GetDiffusionCoeff());

    /*--- Compute the residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = !std::none_of(Res_Visc, Res_Visc + nVar, [](su2double elem){return std::isnan(elem);});
    if(implicit) {
      if(!err) {
        for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
          err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
          err = err  || !std::none_of(Jacobian_j[iVar], Jacobian_j[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
          if(err)
            break;
        }
      }
    }

    /*--- Update the residual and Jacobian in case of implicit simulation ---*/
    if(!err) {
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      LinSysRes.AddBlock(jPoint, Res_Visc);
      if(implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
      }
    }
    else
      throw std::runtime_error("NaN found in the viscous residual");
  } /*--- End for loop over interior edges ---*/
}

//
//
/*--- Impose an isothermal wall boundary condition. ---*/
//
//
void CReactiveNSSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                           CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0)
    SU2_Assert(Vector != NULL,"The array to store velocity for boundary conditions has not been allocated");

  /*--- Local variables ---*/
  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double ktr;
  su2double Tj, dTdn, Twall;
  su2double dij, Area;
  su2double Normal[nDim], UnitNormal[nDim];

  /*--- Identify the boundary ---*/
	auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

	/*--- Retrieve the specified wall temperature ---*/
	Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

	/*--- Loop over boundary points to calculate energy flux ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if(geometry->node[iPoint]->GetDomain()) {
			/*--- Compute dual-grid area and boundary normal ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			Area = ::ComputeArea(Normal, nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			/*--- Compute closest normal neighbor ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between the wall and the normal neighbour ---*/
      auto Coord_i = geometry->node[iPoint]->GetCoord();
      auto Coord_j = geometry->node[Point_Normal]->GetCoord();
      dij = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = std::sqrt(dij);

      /*--- Store the corrected velocity at the wall which will be zero (u = 0)
            unless there is grid motion (u = u_wall) ---*/
      if(grid_movement) {
        auto GridVel = geometry->node[iPoint]->GetGridVel();
        std::copy(GridVel, GridVel + nDim, Vector);
      }
      else
        std::fill(Vector, Vector + nDim, 0.0);

      /*--- Initialize residual to zero ---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
			node[iPoint]->SetVelocity_Old(Vector);
			for(iDim = 0; iDim < nDim; ++iDim) {
        LinSysRes.SetBlock_Zero(iPoint, RHOVX_INDEX_SOL + iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(RHOVX_INDEX_SOL + iDim);
      }

      /*--- Extract the interior temperature and the thermal conductivity for the boundary node ---*/
      Tj   = node[Point_Normal]->GetTemperature();
      ktr  = node[iPoint]->GetThermalConductivity();

      /*--- Compute normal gradient with finite difference approximation ---*/
      dTdn = (Twall - Tj)/dij;

      /*--- Apply to the linear system ---*/
      Res_Visc[RHOE_INDEX_SOL] = ktr*dTdn*Area;

      if(implicit) {
        for(iVar = 0; iVar < nVar; ++iVar)
          for(jVar = 0; jVar < nVar; ++jVar)
            Jacobian_i[iVar][jVar] = 0.0;

        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for(iVar = RHOVX_INDEX_SOL; iVar < RHOVX_INDEX_SOL + nDim; ++iVar) {
          unsigned long total_index = iPoint*nVar + iVar;
          Jacobian.DeleteValsRowi(total_index);
        }

        auto dTdU = node[Point_Normal]->GetdTdU();
        Jacobian_i[RHOE_INDEX_SOL][RHO_INDEX_SOL] = -ktr*dTdU[RHO_INDEX_SOL]/dij*Area;
        Jacobian_i[RHOE_INDEX_SOL][RHOE_INDEX_SOL] = -ktr*dTdU[RHOE_INDEX_SOL]/dij*Area;
        for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          Jacobian_i[RHOE_INDEX_SOL][RHOS_INDEX_SOL + iSpecies] = -ktr*dTdU[RHOS_INDEX_SOL + iSpecies]/dij*Area;

        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }

      /*--- If the wall is moving, there are additional residual contributions
            due to pressure (p v_{wall}\cdot n) and shear stress ((tau*v_{wall}) \cdot n). ---*/
      if(grid_movement) {
        /*--- Get the projected grid velocity at the current boundary node ---*/
        su2double ProjGridVel = std::inner_product(Vector, Vector + nDim, UnitNormal, 0.0);

        /*--- Retrieve other primitive quantities and viscosities ---*/
        su2double Pressure = node[iPoint]->GetPressure();
        su2double mu = node[iPoint]->GetLaminarViscosity();

        unsigned short jDim;
        su2double Grad_Vel[nDim][nDim];
        for(iDim = 0; iDim < nDim; ++iDim)
          for(jDim = 0; jDim < nDim; ++jDim)
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(VX_INDEX_GRAD + iDim, jDim);

        /*--- Divergence of the velocity ---*/
        su2double div_vel = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/
        su2double tau[nDim][nDim];
        for(iDim = 0; iDim < nDim; iDim++) {
          /*--- Compute tau---*/
          for(jDim = 0; jDim < nDim; jDim++)
            tau[iDim][jDim] = mu*(Grad_Vel[iDim][jDim] + Grad_Vel[jDim][iDim]);
          tau[iDim][iDim] -= mu*TWO3*div_vel;
        }

        /*--- Dot product of the stress tensor with the grid velocity ---*/
        su2double tau_vel[nDim];
        std::fill(tau_vel, tau_vel + nDim, 0.0);
        for(iDim = 0; iDim < nDim; ++iDim) {
          for(jDim = 0; jDim < nDim; ++jDim)
            tau_vel[iDim] += tau[iDim][jDim]*Vector[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        Res_Conv[RHOE_INDEX_SOL] = Pressure*ProjGridVel*Area;
        for(iDim = 0; iDim < nDim; iDim++)
          Res_Visc[RHOE_INDEX_SOL] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/
        if(implicit) {
          for(iVar = 0; iVar < nVar; ++iVar)
            for(jVar = 0; jVar < nVar; ++jVar)
              Jacobian_i[iVar][jVar] = 0.0;

          /*--- Convective residual Jacobian related to pressure term ---*/
          auto dPdU = node[iPoint]->GetdPdU();
          for(iVar = 0; iVar < nVar; ++iVar)
            Jacobian_i[RHOE_INDEX_SOL][iVar] = dPdU[iVar]*ProjGridVel*Area;

          /*--- Add the block to the Global Jacobian structure ---*/
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Now the Jacobian contribution related to the shear stress ---*/
          for(iVar = 0; iVar < nVar; ++iVar)
            for(jVar = 0; jVar < nVar; ++jVar)
              Jacobian_i[iVar][jVar] = 0.0;

          su2double factor = mu*Area/(node[iPoint]->GetDensity()*dij);
          su2double theta = std::inner_product(UnitNormal, UnitNormal + nDim, UnitNormal, 0.0);

          if(nDim == 2) {
            su2double thetax = theta + UnitNormal[0]*UnitNormal[0]/3.0;
            su2double thetay = theta + UnitNormal[1]*UnitNormal[1]/3.0;

            su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;

            su2double pix = Vector[0]*thetax + Vector[1]*etaz;
            su2double piy = Vector[0]*etaz   + Vector[1]*thetay;

            Jacobian_i[RHOE_INDEX_SOL][RHO_INDEX_SOL] = -factor*(-pix*Vector[0] + piy*Vector[1]);
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = -factor*pix;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = -factor*piy;
          }
          else {
            su2double thetax = theta + UnitNormal[0]*UnitNormal[0]/3.0;
            su2double thetay = theta + UnitNormal[1]*UnitNormal[1]/3.0;
            su2double thetaz = theta + UnitNormal[2]*UnitNormal[2]/3.0;

            su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            su2double etax = UnitNormal[1]*UnitNormal[2]/3.0;
            su2double etay = UnitNormal[0]*UnitNormal[2]/3.0;

            su2double pix = Vector[0]*thetax + Vector[1]*etaz   + Vector[2]*etay;
            su2double piy = Vector[0]*etaz   + Vector[1]*thetay + Vector[2]*etax;
            su2double piz = Vector[0]*etay   + Vector[1]*etax   + Vector[2]*thetaz;

            Jacobian_i[RHOE_INDEX_SOL][RHO_INDEX_SOL] = -factor*(-pix*Vector[0] + piy*Vector[1] + piz*Vector[2]);
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = -factor*pix;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = -factor*piy;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 2] = -factor*piz;
          }

          /*--- Subtract the block from the Global Jacobian structure ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

        }
      }

      /*--- Apply residuals to the linear system ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if(implicit) {
        for(iVar = RHOVX_INDEX_SOL; iVar < RHOVX_INDEX_SOL + nDim; ++iVar) {
          unsigned long total_index = iPoint*nVar + iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    } /*--- End of GetDomain() check ---*/
  } /*--- End of iVertex for loop ---*/
}

//
//
/*--- Impose a heat flux wall boundary condition. ---*/
//
//
void CReactiveNSSolver::BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                         CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0)
    SU2_Assert(Vector != NULL,"The array to store velocity for boundary conditions has not been allocated");

  /*--- Local variables ---*/
  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Area;
  su2double Normal[nDim], UnitNormal[nDim];

  /*--- Identify the boundary ---*/
	auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

	/*--- Retrieve the specified wall temperature ---*/
  su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

	/*--- Loop over boundary points to calculate energy flux ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if(geometry->node[iPoint]->GetDomain()) {
			/*--- Compute dual-grid area and boundary normal ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			Area = ::ComputeArea(Normal, nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			 /*--- Store the corrected velocity at the wall which will be zero (u = 0)
            unless there is grid motion (u = u_wall) ---*/
      if(grid_movement) {
        auto GridVel = geometry->node[iPoint]->GetGridVel();
        std::copy(GridVel, GridVel + nDim, Vector);
      }
      else
        std::fill(Vector, Vector + nDim, 0.0);

      /*--- Initialize residual to zero ---*/
      for(iVar = 0; iVar < nVar; ++iVar) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
			node[iPoint]->SetVelocity_Old(Vector);
			for(iDim = 0; iDim < nDim; ++iDim) {
        LinSysRes.SetBlock_Zero(iPoint, RHOVX_INDEX_SOL + iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(RHOVX_INDEX_SOL + iDim);
      }

      /*--- Apply to the linear system ---*/
      Res_Visc[RHOE_INDEX_SOL] = Wall_HeatFlux*Area;

      /*--- If the wall is moving, there are additional residual contributions
            due to pressure (p v_{wall}\cdot n) and shear stress ((tau*v_{wall}) \cdot n). ---*/
      if(grid_movement) {
        /*--- Get the projected grid velocity at the current boundary node ---*/
        su2double ProjGridVel = std::inner_product(Vector, Vector + nDim, UnitNormal, 0.0);

        /*--- Retrieve other primitive quantities and viscosities ---*/
        su2double Pressure = node[iPoint]->GetPressure();
        su2double mu = node[iPoint]->GetLaminarViscosity();

        unsigned short jDim;
        su2double Grad_Vel[nDim][nDim];
        for(iDim = 0; iDim < nDim; ++iDim)
          for(jDim = 0; jDim < nDim; ++jDim)
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(VX_INDEX_GRAD + iDim, jDim);

        /*--- Divergence of the velocity ---*/
        su2double div_vel = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/
        su2double tau[nDim][nDim];
        for(iDim = 0; iDim < nDim; iDim++) {
          /*--- Compute tau---*/
          for(jDim = 0; jDim < nDim; jDim++)
            tau[iDim][jDim] = mu*(Grad_Vel[iDim][jDim] + Grad_Vel[jDim][iDim]);
          tau[iDim][iDim] -= mu*TWO3*div_vel;
        }

        /*--- Dot product of the stress tensor with the grid velocity ---*/
        su2double tau_vel[nDim];
        std::fill(tau_vel, tau_vel + nDim, 0.0);
        for(iDim = 0; iDim < nDim; ++iDim) {
          for(jDim = 0; jDim < nDim; ++jDim)
            tau_vel[iDim] += tau[iDim][jDim]*Vector[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        Res_Conv[RHOE_INDEX_SOL] = Pressure*ProjGridVel*Area;
        for(iDim = 0; iDim < nDim; iDim++)
          Res_Visc[RHOE_INDEX_SOL] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/
        if(implicit) {
          for(iVar = 0; iVar < nVar; ++iVar)
            for(jVar = 0; jVar < nVar; ++jVar)
              Jacobian_i[iVar][jVar] = 0.0;

          /*--- Convective residual Jacobian related to pressure term ---*/
          auto dPdU = node[iPoint]->GetdPdU();
          for(iVar = 0; iVar < nVar; ++iVar)
            Jacobian_i[RHOE_INDEX_SOL][iVar] = dPdU[iVar]*ProjGridVel*Area;

          /*--- Add the block to the Global Jacobian structure ---*/
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

          /*--- Now the Jacobian contribution related to the shear stress ---*/
          for(iVar = 0; iVar < nVar; ++iVar)
            for(jVar = 0; jVar < nVar; ++jVar)
              Jacobian_i[iVar][jVar] = 0.0;

          /*--- Compute closest normal neighbor ---*/
          Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

          /*--- Compute distance between wall & normal neighbor ---*/
          auto Coord_i = geometry->node[iPoint]->GetCoord();
          auto Coord_j = geometry->node[Point_Normal]->GetCoord();
          su2double dij = 0.0;
          for(iDim = 0; iDim < nDim; ++iDim)
            dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
          dij = std::sqrt(dij);

          su2double factor = mu*Area/(node[iPoint]->GetDensity()*dij);
          su2double theta = std::inner_product(UnitNormal, UnitNormal + nDim, UnitNormal, 0.0);

          if(nDim == 2) {
            su2double thetax = theta + UnitNormal[0]*UnitNormal[0]/3.0;
            su2double thetay = theta + UnitNormal[1]*UnitNormal[1]/3.0;

            su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;

            su2double pix = Vector[0]*thetax + Vector[1]*etaz;
            su2double piy = Vector[0]*etaz   + Vector[1]*thetay;

            Jacobian_i[RHOE_INDEX_SOL][RHO_INDEX_SOL] = -factor*(-pix*Vector[0] + piy*Vector[1]);
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = -factor*pix;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = -factor*piy;
          }
          else {
            su2double thetax = theta + UnitNormal[0]*UnitNormal[0]/3.0;
            su2double thetay = theta + UnitNormal[1]*UnitNormal[1]/3.0;
            su2double thetaz = theta + UnitNormal[2]*UnitNormal[2]/3.0;

            su2double etaz = UnitNormal[0]*UnitNormal[1]/3.0;
            su2double etax = UnitNormal[1]*UnitNormal[2]/3.0;
            su2double etay = UnitNormal[0]*UnitNormal[2]/3.0;

            su2double pix = Vector[0]*thetax + Vector[1]*etaz   + Vector[2]*etay;
            su2double piy = Vector[0]*etaz   + Vector[1]*thetay + Vector[2]*etax;
            su2double piz = Vector[0]*etay   + Vector[1]*etax   + Vector[2]*thetaz;

            Jacobian_i[RHOE_INDEX_SOL][RHO_INDEX_SOL] = -factor*(-pix*Vector[0] + piy*Vector[1] + piz*Vector[2]);
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = -factor*pix;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = -factor*piy;
            Jacobian_i[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 2] = -factor*piz;
          }

          /*--- Subtract the block from the Global Jacobian structure ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

        }
      }

      /*--- Apply residuals to the linear system ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if(implicit) {
        for(iVar = RHOVX_INDEX_SOL; iVar < RHOVX_INDEX_SOL + nDim; ++iVar) {
          unsigned long total_index = iPoint*nVar + iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    } /*--- End of GetDomain() check ---*/
  } /*--- End of iVertex for loop ---*/
}

//
//
/*--- Solution Of Stefan-Maxwell equations ---*/
//
//
CReactiveNSSolver::Vec CReactiveNSSolver::Solve_SM(const su2double val_density, const su2double val_alpha, const RealMatrix& val_Dij,
                                                   const RealVec& val_xs, const Vec& val_grad_xs, const RealVec& val_ys) {
  /*--- Set tolerance ---*/
  const su2double toll = 1e-11;

  /*--- Rename for convenience ---*/
  su2double rho = val_density;
  su2double alpha = val_alpha;
  RealMatrix Gamma_tilde(nSpecies, nSpecies);

  /*--- Compute original matrix of Stefan-Maxwell equations ---*/
  auto Gamma = library->GetGamma(rho, val_xs, val_ys, val_Dij);

  /*--- Add artificial diffusion part ---*/
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
      Gamma_tilde(iSpecies,jSpecies) = Gamma(iSpecies,jSpecies) + alpha*val_ys[iSpecies];

  Eigen::BiCGSTAB<RealMatrix> bicg(Gamma_tilde);
  bicg.setTolerance(toll);
  return bicg.solve(-val_grad_xs);
}

//
//
/*--- Regression boundary condition. NOTE: TO BE CLIPPED... ---*/
//
//
void CReactiveNSSolver::BC_Engine_Inflow(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                         CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0) {
    SU2_Assert(Res_Conv != NULL,"The array of convective residual for boundary conditions has not been allocated");
    SU2_Assert(Res_Visc != NULL,"The array of viscous residual for boundary conditions has not been allocated");
  }

  if(grid_movement)
    throw Common::NotImplemented("Regression BC not available in case of grid movement");

  /*--- Local variables ---*/
  unsigned short iDim, iSpecies, jSpecies;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Area;
  su2double Normal[nDim], UnitNormal[nDim], Coord_ij[nDim];
  su2double pg, rho_g, Tg, Tg_old, Tg_tmp, Vg, Vg_old;
  RealVec Ys_g(nSpecies), Xs_g, Ys_int(nSpecies), Xs_int;
  su2double Ys_g_old[nSpecies], Ys_g_tmp[nSpecies];
  RealMatrix fprime_Ys(nSpecies,nSpecies);
  Vec Jd, Jd_pert, f_Ys(nSpecies), deltaf_Ys, Grad_Xs_norm(nSpecies);

  /*--- Identify the boundary ---*/
	auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  SU2_Assert(config->GetnMarker_EngineInflow() == 1, "Only one fuel inflow is allowed");

  /*--- Read surface and fuel data ---*/
  SU2_Assert(config->GetnSpecies_Inflow() == nSpecies, "Wrong number of species detected at inflow");
  Ys = RealVec(config->GetInflow_MassFrac(Marker_Tag), config->GetInflow_MassFrac(Marker_Tag) + nSpecies);
  SU2_Assert(std::abs(std::accumulate(Ys.cbegin(), Ys.cend(), 0.0) - 1.0) < 1.0e-10,
             "The mass fractions don't sum up to 1 at inlet boundary");
  auto Flow_Dir = config->GetVelocityDir_Inflow();
  const su2double rho_s = config->GetDensity_Fuel()/config->GetDensity_Ref();
  const su2double Cp_s = config->GetSpecificHeat_Fuel()/config->GetEnergy_Ref()*config->GetTemperature_Ref();
  const su2double h_pf = config->GetEnthalpy_Fuel()/config->GetEnergy_Ref();
  const su2double kappa_s = config->GetConductivity_Fuel()/config->GetConductivity_Ref();
  const su2double T0 = config->GetTemperature_Fuel()/config->GetTemperature_Ref();

  /*--- Set fuel regression rate data ---*/
  if(config->GetExtIter() == 0)
    library->ReadDataFuel(config->GetFuelData_File());

  /*--- Auxiliary function for energy balance at wall ---*/
  auto f_wall = std::function<su2double(su2double)>([&](double T){
    su2double dim_temp = T*config->GetTemperature_Ref();
    if(US_System)
      dim_temp *= 5.0/9.0;
    su2double rb = library->ComputeRegressionRate(dim_temp)/config->GetVelocity_Ref();
    if(US_System)
      rb *= 3.28084;
    auto hg_s = library->ComputePartialEnthalpy(dim_temp);
    for(auto& elem: hg_s)
      elem /= config->GetEnergy_Ref();
    if(US_System) {
      for(auto& elem: hg_s)
        elem *= 3.28084*3.28084;
    }
    su2double hg = std::inner_product(hg_s.cbegin(), hg_s.cend(), Ys_g.cbegin(), 0.0);
    su2double hgf = std::inner_product(hg_s.cbegin(), hg_s.cend(), Ys.cbegin(), 0.0);
    su2double kappa_g = library->ComputeLambda(dim_temp, Ys_g)/config->GetConductivity_Ref();
    if(US_System)
      kappa_g /= 8.006796;
    su2double Grad_Tw = -rho_s*Cp_s*rb*(T - T0)/kappa_s;
    su2double partial_res = 0.0;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      partial_res += (Ys[iSpecies] - Ys_g[iSpecies])*hg_s[iSpecies];
    return rb*rho_s*(hg - hgf + h_pf + partial_res) + Grad_Tw*(kappa_s - kappa_g);
  });

  /*--- Limits for regula-falsi method ---*/
  su2double Tmin = config->GetTemperatureMin();
  su2double Tmax = config->GetTemperatureMax();

  /*--- Loop over boundary points to compute flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if(geometry->node[iPoint]->GetDomain()) {
			/*--- Compute dual-grid area and boundary normal ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
			Area = ::ComputeArea(Normal, nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			/*--- Compute closest normal neighbor ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      auto Coord_i = geometry->node[iPoint]->GetCoord();
      auto Coord_j = geometry->node[Point_Normal]->GetCoord();

      for(iDim = 0; iDim < nDim; ++iDim)
        Coord_ij[iDim] = std::abs(Coord_j[iDim] - Coord_i[iDim]);

      /*--- Get characteristic variables at inlet ---*/
      auto V_inlet = GetCharacPrimVar(val_marker, iVertex);

      /*--- Extrapolate pressure from the internal node ---*/
      pg = node[Point_Normal]->GetPressure();
      su2double dim_press = pg*config->GetPressure_Ref()/101325.0;
      if(US_System)
        dim_press *= 47.8803;

      /*--- Set iterations for complete non-linear system ---*/
      unsigned short maxIter = 10;
      unsigned short iIter;
      const su2double toll_sys = 1.0e-3;
      bool conv_sys = false;

      /*--- Set iterations for mass fractions and temperature system ---*/
      unsigned short maxSubIter = 15;
      unsigned short iSubIter;
      const su2double toll_subsys = 1.0e-6;
      bool conv_subsys = false;
      Tg = node[Point_Normal]->GetTemperature();
      Vg = 0.0;
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Ys_g[iSpecies] = Ys_int[iSpecies] = node[Point_Normal]->GetMassFraction(iSpecies);
      Xs_int = library->GetMolarFromMass(Ys_int);

      /*--- Start non-linear system solution procedure ---*/
      for(iIter = 0; iIter < maxIter; ++iIter) {
        /*--- Save previous iteration values ---*/
        Tg_old = Tg;
        Vg_old = Vg;
        std::copy(Ys_g.cbegin(), Ys_g.cend(), Ys_g_old);

        su2double dim_temp = Tg*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        auto val_Dij = library->GetDij_SM(dim_temp, dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1.0e4);
        su2double alpha_Dij = 1.0/(node[Point_Normal]->GetDensity()*val_Dij.maxCoeff());

        /*--- Solve iteratively mass fractions ---*/
        su2double omega_bar = rho_s*library->ComputeRegressionRate(dim_temp)/config->GetVelocity_Ref();
        if(US_System)
          omega_bar *= 3.28084;
        for(iSubIter = 0; iSubIter < maxSubIter; ++iSubIter) {
          /*--- Save previous subiteration values ---*/
          std::copy(Ys_g.cbegin(), Ys_g.cend(), Ys_g_tmp);

          /*--- NOTE: Compute molar fractions and its gradient ----*/
          Xs_g = library->GetMolarFromMass(Ys_g);

          /*--- Solve S-M equations ---*/
          su2double rho_g = library->ComputeDensity(Tg, pg, Ys_g)*config->GetGas_Constant_Ref();
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            Grad_Xs_norm[iSpecies] = 0.0;
            for(iDim = 0; iDim < nDim; ++iDim)
              Grad_Xs_norm[iSpecies] += (Xs_int[iSpecies] - Xs_g[iSpecies])/Coord_ij[iDim]*UnitNormal[iDim];
          }
          Jd = Solve_SM(rho_g, alpha_Dij, val_Dij, Xs_g, Grad_Xs_norm, Ys_g);

          const su2double eps = 1.0e-7; /*--- Perturbation for numerical Jacobian ---*/
          const su2double toll = 1.0e-9; /*--- Tolerance to solve iteratively the linear system to find update step ---*/
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            /*--- Set current value of flux ---*/
            f_Ys[iSpecies] = (Ys_g[iSpecies] - Ys[iSpecies])*omega_bar + Jd[iSpecies];

            /*--- Perturb mass fraction to compute perturbed diffusion fluxes ---*/
            if(std::abs(Ys_g[iSpecies]) <= 1.0)
              Ys_g[iSpecies] += eps;
            else
              Ys_g[iSpecies] += eps*Ys_g[iSpecies];
            Xs_g = library->GetMolarFromMass(Ys_g);
            for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
              Grad_Xs_norm[jSpecies] = 0.0;
              for(iDim = 0; iDim < nDim; ++iDim)
                Grad_Xs_norm[jSpecies] += (Xs_int[jSpecies] - Xs_g[jSpecies])/Coord_ij[iDim]*UnitNormal[iDim];
            }
            Jd_pert = Solve_SM(rho_g, alpha_Dij, val_Dij, Xs_g, Grad_Xs_norm, Ys_g);

            /*--- Populate jacobian ---*/
            for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
              fprime_Ys(jSpecies,iSpecies) = (Jd_pert[jSpecies] - Jd[jSpecies])/(Ys_g[iSpecies] - Ys_g_tmp[iSpecies]);
            fprime_Ys(iSpecies,iSpecies) += omega_bar;

            /*--- Restore value of mass fraction to pass to subsequent derivatives ---*/
            Ys_g[iSpecies] = Ys_g_tmp[iSpecies];
          }

          /*--- Find update step and apply to mass fractions ---*/
          Eigen::BiCGSTAB<RealMatrix> bicg(fprime_Ys);
          bicg.setTolerance(toll);
          deltaf_Ys = bicg.solve(f_Ys);
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
            Ys_g[iSpecies] += -deltaf_Ys[iSpecies];

          /*--- Check convergecne ---*/
          su2double err_Ys = std::numeric_limits<su2double>::min();
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
            err_Ys = std::max(err_Ys, std::abs(Ys_g[iSpecies] - Ys_g_tmp[iSpecies])/Ys_g_tmp[iSpecies]);
          if(err_Ys < toll_subsys) {
            conv_subsys = true;
            break;
          }
        }
        if(!conv_subsys)
          throw std::runtime_error("Convergence not achieved for mass fractions system at inflow");

        /*--- Solve iteratively wall temperature ---*/
        conv_subsys = false;
        su2double Ta = Tmin;
        su2double Tb = Tmax;
        su2double fa = f_wall(Ta);
        su2double fb = f_wall(Tb);
        for(iSubIter = 0; iSubIter < maxSubIter; ++iSubIter) {
          /*--- Save previous subiteration value ---*/
          Tg_tmp = Tg;

          /*--- Regula falsi scheme update ----*/
          Tg = std::min(std::max((Ta*fb - Tb*fa)/(fb - fa), Tmax), Tmin);

          /*--- Clip temperature in order to avoid too large variations ---*/
          Tg = std::min(std::max(Tg,Tg_tmp - 20.0),Tg_tmp + 20.0);

          /*--- Check convergence ---*/
          if(std::abs(Tg - Tg_tmp)/Tg_tmp < toll_subsys) {
            conv_subsys = true;
            break;
          }

          /*--- Update limits ---*/
          if(f_wall(Tg)*fa < 0) {
            Tb = Tg;
            fb = f_wall(Tb);
          }
          else {
            Ta = Tg;
            fa = f_wall(Ta);
          }
        }
        if(!conv_subsys)
          throw std::runtime_error("Convergence not achieved for wall temperature at inflow");

        /*--- Compute wall density ---*/
        rho_g = library->ComputeDensity(Tg, pg, Ys_g)*config->GetGas_Constant_Ref();

        /*--- Compute wall velocity ---*/
        su2double alpha = std::inner_product(UnitNormal, UnitNormal + nDim, Flow_Dir, 0.0);
        dim_temp = Tg*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        su2double rb = library->ComputeRegressionRate(dim_temp)/config->GetVelocity_Ref();
        if(US_System)
          rb *= 3.28084;
        Vg = (rho_s - rho_g)*rb/(rho_g*alpha);

        /*--- Check global convergecne ---*/
        su2double err_Ys = std::numeric_limits<su2double>::min();
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          err_Ys = std::max(err_Ys, std::abs(Ys_g[iSpecies] - Ys_g_old[iSpecies])/Ys_g_old[iSpecies]);
        su2double err = std::max(std::abs(Vg_old - Vg)/Vg_old, std::max(std::abs(Tg_old - Tg)/Tg_old, err_Ys));
        if(err < toll_sys) {
          conv_sys = true;
          break;
        }
      } /*--- End full system cycle ---*/
      if(!conv_sys)
        throw std::runtime_error("Convergence not achieved for regression boundary condition system");

      /*--- Set inlet state ---*/
      V_inlet[T_INDEX_PRIM] = Tg;
      for(iDim = 0; iDim < nDim; ++iDim)
        V_inlet[VX_INDEX_PRIM + iDim] = Vg*Flow_Dir[iDim];
      V_inlet[P_INDEX_PRIM] = pg;
      V_inlet[RHO_INDEX_PRIM] = rho_g;
      su2double dim_temp = Tg*config->GetTemperature_Ref();
      if(US_System)
        dim_temp *= 5.0/9.0;
      V_inlet[H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp, Ys_g)/config->GetEnergy_Ref();
      if(US_System)
        V_inlet[H_INDEX_PRIM] *= 3.28084*3.28084;
      V_inlet[H_INDEX_PRIM] += 0.5*Vg*Vg;
      V_inlet[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp, Ys_g)/config->GetVelocity_Ref();
      if(US_System)
        V_inlet[A_INDEX_PRIM] *= 3.28084;
      std::copy(Ys.cbegin(), Ys.cend(), V_inlet + RHOS_INDEX_PRIM);

      /*--- Set normal for convective flux ---*/
      conv_numerics->SetNormal(Normal);

      /*--- Set primitive for convective flux ---*/
      conv_numerics->SetPrimitive(node[iPoint]->GetPrimitive(), V_inlet);

      if(implicit) {
        su2double Secondary[nVar];
        su2double Gamma = library->ComputeFrozenGamma_FromSoundSpeed(V_inlet[T_INDEX_PRIM], V_inlet[A_INDEX_PRIM], Ys);
        Gamma *= config->GetGas_Constant_Ref();

        /*--- Derivatives with respect to density ---*/
        Secondary[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*Vg*Vg;

        /*--- Derivatives with respect to momentum ---*/
        for(iDim = 0; iDim < nDim; ++iDim)
          Secondary[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V_inlet[VX_INDEX_PRIM + iDim];

        /*--- Derivatives with respect to energy ---*/
        Secondary[RHOE_INDEX_SOL] = Gamma - 1.0;

        /*--- Derivatives with respect to partial densities ---*/
        su2double dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Secondary[RHOS_INDEX_SOL + iSpecies] = library->ComputedP_dYs(dim_temp, Gamma, iSpecies)/config->GetEnergy_Ref();
          if(US_System)
            Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
        }
        conv_numerics->SetSecondary(node[iPoint]->GetdPdU(), Secondary);
      }

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      bool err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      if(!err) {
        /*--- Update residual value ---*/
        LinSysRes.AddBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the convective residual of regression boundary condition.");

      /*--- Viscous contribution ---*/
      /*--- Set the normal vector and the coordinates ---*/
      visc_numerics->SetNormal(Normal);
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

      /*--- Primitive variables, and gradient ---*/
      visc_numerics->SetPrimitive(node[iPoint]->GetPrimitive(), V_inlet);
      visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

      /*--- Set temperature derivatives if implicit ---*/
      if(implicit) {
        su2double Secondary[nVar];

        su2double dim_temp = V_inlet[T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        su2double Cv = library->ComputeCV(dim_temp, Ys)/config->GetGas_Constant_Ref();
        if(US_System)
          Cv *= 3.28084*3.28084*5.0/9.0;
        su2double rhoCv = V_inlet[RHO_INDEX_PRIM]*Cv;
        su2double sq_vel = std::inner_product(V_inlet + VX_INDEX_PRIM, V_inlet + (VX_INDEX_PRIM + nDim), V_inlet + VX_INDEX_PRIM, 0.0);

        /*--- Set temperature derivatives ---*/
        Secondary[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
        for(unsigned short iDim = 0; iDim < nDim; ++iDim)
          Secondary[RHOVX_INDEX_SOL + iDim] = -V_inlet[VX_INDEX_PRIM + iDim]/rhoCv;
        Secondary[RHOE_INDEX_SOL] = 1.0/rhoCv;
        for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          Secondary[RHOS_INDEX_SOL + iSpecies] = -library->ComputePartialEnergy(dim_temp, iSpecies)/(config->GetEnergy_Ref()*rhoCv);
          if(US_System)
            Secondary[RHOS_INDEX_SOL + iSpecies] *= 3.28084*3.28084;
        }
        visc_numerics->SetSecondary(node[iPoint]->GetdTdU(), Secondary);
      }

      /*--- Laminar viscosity ---*/
      visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

      /*--- Thermal conductivity ---*/
      visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

      /*--- Species binary coefficients ---*/
      visc_numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[iPoint]->GetDiffusionCoeff());

      /*--- Compute residual ---*/
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Check for NaNs before applying the residual to the linear system ---*/
      err = !std::none_of(Residual, Residual + nVar, [](su2double elem){return std::isnan(elem);});
      if(implicit) {
        if(!err) {
          for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
            err = !std::none_of(Jacobian_i[iVar], Jacobian_i[iVar] + nVar, [](su2double elem){return std::isnan(elem);});
            if(err)
              break;
          }
        }
      }

      if(!err) {
        /*--- Update residual ---*/
        LinSysRes.SubtractBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
      else
        throw std::runtime_error("NaN found in the viscous residual of regression boundary condition");
    }
  } /*--- End loop boundary points ---*/
}
