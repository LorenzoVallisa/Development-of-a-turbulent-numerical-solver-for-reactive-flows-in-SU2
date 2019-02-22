#include "../include/solver_reactive.hpp"
#include "../include/numerics_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/not_implemented_exception.hpp"

#include <algorithm>
#include <limits>
#include <cmath>
#include <iterator>

namespace {
  /*!
   * \brief Compute area for the current normal
   */
  su2double ComputeArea(su2double* Normal, const unsigned short nDim) {
    su2double Area = std::inner_product(Normal, Normal + nDim, Normal, 0.0);
    Area = std::sqrt(Area);
    return Area;
  }
} /*-- End of unnamed namespace ---*/

/*--- Get the library to compute physical-chemical properties ---*/
CReactiveEulerSolver::LibraryPtr CReactiveEulerSolver::library = CReactiveEulerVariable::GetLibrary();

//
//
/*!
  *\brief Default constructor
  */
//
//
CReactiveEulerSolver::CReactiveEulerSolver(): CSolver(), nSpecies(), nPrimVarLim(), space_centered(), implicit(), grid_movement(),
                                              least_squares(), second_order(), limiter(), Density_Inf(), Pressure_Inf(), Temperature_Inf() {
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
/*!
  *\brief Class constructor
  */
//
//
CReactiveEulerSolver::CReactiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh):
                      CSolver(), nSpecies(library->GetNSpecies()) {
  unsigned long iPoint;
  unsigned short iVar,iDim;

  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  IterLinSolver = 0;

  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  std::string file_name = config->GetSolution_FlowFileName();

  if(!(!restart || iMesh != MESH_0)) {
    /*--- Modify file name for a time stepping unsteady restart ---*/
    if(time_stepping) {
      auto Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter());
      file_name = config->GetUnsteady_FileName(file_name, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/
    //Read_SU2_Restart_Metadata(geometry, config, false, file_name);
  }

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho,rho*vx,rho*vy,rho*vz,rho*E,rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T,vx,vy,vz,P,rho,h,a,rho1,...rhoNs)^T ---*/
  nPrimVarLim = nDim + 2; /*--- Gradient Primitive variables (T,vx,vy,vz,P)^T ---*/
  nPrimVarGrad = nPrimVarLim; /*--- Primitive variables to limit (T,vx,vy,vz,P)^T ---*/

  /*--- Perform the non-dimensionalization for the flow equations using the specified reference values. ---*/
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

  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];
  Solution_i = new su2double[nVar];
  Solution_j = new su2double[nVar];

	/*--- Allocate vectors related to the primitive ---*/
  PrimVar_i.resize(nPrimVarGrad);
  PrimVar_j.resize(nPrimVarGrad);
  PrimVar_Vertex.resize(nPrimVarGrad);

  Prim_i.resize(nPrimVarLim);
  Prim_j.resize(nPrimVarLim);
  Primitive.resize(nPrimVarLim);

  /*--- Allocate vectors for conserved variable limits ---*/
  Lower_Limit.resize(nPrimVar, 0.0);
  Upper_Limit.resize(nPrimVar, 1.0/EPS);

  std::fill(Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL,
            Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL + nDim, -1.0/EPS);
  Lower_Limit[CReactiveEulerVariable::RHOE_INDEX_SOL] = -1.0/EPS;

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
		Jacobian_i = new su2double*[nVar];
		Jacobian_j = new su2double*[nVar];
		for(iVar = 0; iVar < nVar; ++iVar) {
			Jacobian_i[iVar] = new su2double[nVar];
			Jacobian_j[iVar] = new su2double[nVar];
		}
  }
  else
    implicit = false;

  /*--- Allocate matrices for gradient computation by least squares ---*/
	if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    least_squares = true;

    Smatrix = new su2double*[nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      Smatrix[iDim] = new su2double[nDim];

    Cvector = new su2double*[nPrimVarGrad];
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      Cvector[iVar] = new su2double[nDim];
  }
  else
    least_squares = false;

  grid_movement = config->GetGrid_Movement();
  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  second_order = (config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  unsigned long ExtIter = config->GetExtIter();
  limiter = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter());

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(),config->GetVelocity_FreeStreamND() + nDim);
  MassFrac_Inf     = RealVec(config->GetMassFrac_FreeStream(),config->GetMassFrac_FreeStream() + nSpecies);

  node_infty = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, nSpecies,
                                          nPrimVar, nPrimVarGrad, nPrimVarLim, config);

  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  auto direct_diff = config->GetDirectDiff();
  switch(direct_diff) {
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
    break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
    break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
    break;
    default:
    break;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  for(iPoint = 0; iPoint < nPoint; ++iPoint)
    node[iPoint] = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
                                              nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, config);

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
}


//
//
/*!
 * \brief Check whether there is any non physical point in the initialization
 */
//
//
void CReactiveEulerSolver::Check_FreeStream_Solution(CConfig* config) {
  int rank;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  bool check_infty = node_infty->SetPrimVar(config);
  SU2_Assert(check_infty == true, "Assigned freestream values are not compatible with physics");

  unsigned long counter_local = 0, counter_global;
  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    bool nonPhys = node[iPoint]->SetPrimVar(config);
    if(nonPhys) {
      su2double rho, hs_Tot;
      su2double sqvel = std::inner_product(Velocity_Inf.cbegin(),Velocity_Inf.cend(),Velocity_Inf.cbegin(),0.0);

      /*--- Compute density from supplied quantities ---*/
      //rho = library->ComputeDensity(Pressure_Inf,Temperature_Inf,MassFrac_Inf);
      rho *= config->GetGas_Constant_Ref();
      Solution[CReactiveEulerVariable::RHO_INDEX_SOL] = rho;
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Solution[CReactiveEulerVariable::RHOS_INDEX_SOL + iSpecies] = MassFrac_Inf[iSpecies]*rho;

      /*--- Compute momentum from supplied primitive quanitites ---*/
      for(unsigned short iDim = 0; iDim < nDim; ++iDim)
        Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = rho*Velocity_Inf[iDim];

      /*---- Compute energy from temperature ---*/
      su2double dim_temp = Temperature_Inf*config->GetTemperature_Ref();
      bool US_System = (config->GetSystemMeasurements() == US);
      if(US_System)
        dim_temp *= 5.0/9.0;
      //hs_Tot = library->ComputeEnthalpy(Temperature_Inf,MassFrac_Inf)/config->GetEnergy_Ref()
      if(US_System)
        hs_Tot *= 3.28084*3.28084;
      Solution[CReactiveEulerVariable::RHOE_INDEX_SOL] = rho*(hs_Tot + 0.5*sqvel) - Pressure_Inf;

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
}


//
//
/*!
 * \brief Set the fluid solver nondimensionalization.
 */
void CReactiveEulerSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {
   su2double Temperature_FreeStream = 0.0, ModVel_FreeStream = 0.0, Energy_FreeStream = 0.0,
             Density_FreeStream = 0.0, Pressure_FreeStream = 0.0, SoundSpeed_FreeStream = 0.0;

   su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
             Temperature_Ref = 0.0, Time_Ref = 0.0, Gas_Constant_Ref = 0.0, Energy_Ref = 0.0;

   su2double Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,ModVel_FreeStreamND = 0.0,
             Temperature_FreeStreamND = 0.0, Energy_FreeStreamND = 0.0, Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

   su2double Velocity_FreeStreamND[nDim];

   unsigned short iDim;

   /*--- Check if the simulation is unsteady ---*/
   bool unsteady = config->GetUnsteady_Simulation();

   int rank;
   #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   #endif

 	 /*--- Compressible non dimensionalization ---*/
   Pressure_FreeStream = config->GetPressure_FreeStream();
   Temperature_FreeStream = config->GetTemperature_FreeStream();
   Density_FreeStream = library->ComputeDensity(Pressure_FreeStream,Temperature_FreeStream,MassFrac_Inf);
   config->SetDensity_FreeStream(Density_FreeStream);
   //SoundSpeed_FreeStream = library->ComputeFrozenSoundSpeed(Temperature_FreeStream,MassFrac_Inf);

   /*--- Compute the modulus of the free stream velocity ---*/
   ModVel_FreeStream = std::inner_product(config->GetVelocity_FreeStream(), config->GetVelocity_FreeStream() + nDim,
                                          config->GetVelocity_FreeStream(), 0.0);
   ModVel_FreeStream = std::sqrt(ModVel_FreeStream);
   config->SetModVel_FreeStream(ModVel_FreeStream);

   /*--- Compute the free stream energy ---*/
   //Energy_FreeStream = library->ComputeEnergy(Temperature_FreeStream,MassFrac_Inf) + 0.5*ModVel_FreeStream*ModVel_FreeStream;

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
     throw std::out_of_range("Invalid option for adimensionalitazion");

   config->SetPressure_Ref(Pressure_Ref);
   config->SetDensity_Ref(Density_Ref);
   config->SetTemperature_Ref(Temperature_Ref);
   config->SetLength_Ref(Length_Ref);

   Velocity_Ref = std::sqrt(Pressure_Ref/Density_Ref);
   config->SetVelocity_Ref(Velocity_Ref);

   Time_Ref = Length_Ref/Velocity_Ref;
   config->SetTime_Ref(Time_Ref);

   //Energy_Ref = library->ComputeEnergy(Temperature_Ref,MassFrac_Inf) + 0.5*Velocity_Ref*Velocity_Ref;
   config->SetEnergy_Ref(Energy_Ref);

   Gas_Constant_Ref = Velocity_Ref*Velocity_Ref/Temperature_Ref;
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

   ModVel_FreeStreamND = std::inner_product(Velocity_FreeStreamND,Velocity_FreeStreamND + nDim,Velocity_FreeStreamND,0.0);
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
     bool US_Measurament = (config->GetSystemMeasurements() == US);

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
 /*!
   *\brief Setting primitive variables in the preprocessing
   */
 //
 //
unsigned long CReactiveEulerSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {
  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Initialize the non-physical points vector ---*/
    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies + nDim + 5, (T, vx, vy, vz, P, rho, h, a, Y1,....YNs) ---*/
    RightSol = node[iPoint]->SetPrimVar(config);

    if(!RightSol) {
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
/*!
  *\brief Variables preprocessing
  */
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
    throw Common::NotImplemented("Implicit scheme not implemented. Switch to explicit");

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
/*!
 *\brief Set gradient primitive variables according to Green-Gauss
 */
//
//
void CReactiveEulerSolver::SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, jDim, iVar, iMarker;

  su2double PrimVar_Average, Partial_Res;
  su2double Volume;
  su2double Normal[nDim];

  /*--- Gradient of the following primitive variables: [T,u,v,w,p]^T ---*/
  /*--- Set Gradient_Primitive to zero ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
  	node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
  	iPoint = geometry->edge[iEdge]->GetNode(0);
  	jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
    PrimVar_i[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
    PrimVar_j[CReactiveEulerVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
    PrimVar_i[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
    PrimVar_j[CReactiveEulerVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim) {
      PrimVar_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
      PrimVar_j[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
    }

  	Normal = geometry->edge[iEdge]->GetNormal();

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
        /*--- Get primitives from CVariable ---*/
        PrimVar_Vertex[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
        PrimVar_Vertex[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
        for(iDim = 0; iDim < nDim; ++iDim)
          PrimVar_Vertex[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] =
          node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
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
/*!
 *\brief Set gradient primitive variables according to least squares
 */
//
//
void CReactiveEulerSolver::SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) {
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
    PrimVar_i[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
    PrimVar_i[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim)
      PrimVar_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

    /*--- Inizialization of variables ---*/
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      std::fill(Cvector[iVar], Cvector[iVar] + nDim, 0.0);

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		//AD::StartPreacc();
    //AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    //AD::SetPreaccIn(Coord_i, nDim);

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			auto Coord_j = geometry->node[jPoint]->GetCoord();

      /*--- Get primitives from CVariable ---*/
      PrimVar_j[CReactiveEulerVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      PrimVar_j[CReactiveEulerVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        PrimVar_j[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

			//AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
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
            Cvector[iVar][iDim] += Coord_ij[iDim]*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
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
        node[iPoint]->SetGradient_Primitive(iVar,iDim,result);
      }
    }

	//	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
  //  AD::EndPreacc();
  } /*--- End of iPoint for loop ---*/

	Set_MPI_Primitive_Gradient(geometry, config);
}


//
//
/*!
 *\brief Set limiter of primitive variables
 */
//
//
void CReactiveEulerSolver::SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) {
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;

  su2double dave, LimK, eps1, eps2, dm, dp, du, y, limiter;

  dave = config->GetRefElemLength();
  LimK = config->GetLimiterCoeff();

  if(!limiter) {
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint){
      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        node[iPoint]->SetLimiter_Primitive(iVar,1.0);
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
      Prim_i[CReactiveEulerVariable::T_INDEX_LIM] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Prim_i[CReactiveEulerVariable::P_INDEX_LIM] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      Prim_j[CReactiveEulerVariable::T_INDEX_LIM] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Prim_j[CReactiveEulerVariable::P_INDEX_LIM] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Prim_i[CReactiveEulerVariable::VX_INDEX_LIM + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
        Prim_j[CReactiveEulerVariable::VX_INDEX_LIM + iDim] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
      }

      /*--- Compute the max and min values for nodes i and j ---*/
      for(iVar = 0; iVar < nPrimVarLim; ++iVar){
        du = (Prim_j[iVar] - Prim_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
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

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_j[iDim] - Coord_i[iDim])*Gradient_i[iVar][iDim];

        if(dm == 0.0)
          limiter = 2.0;
        else {
          if(dm > 0.0)
            dp = node[iPoint]->GetSolution_Max(iVar);
          else
            dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if(limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          //AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if(dm == 0.0)
          limiter = 2.0;
        else {
          if(dm > 0.0)
            dp = node[jPoint]->GetSolution_Max(iVar);
          else
            dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if(limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          //AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
      } /*--- End of limiter computation for a single edge ---*/
      //AD::EndPreacc();
    } /*--- End of iEdge for loop ---*/

    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
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
      Primitive[CReactiveEulerVariable::T_INDEX_LIM] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Primitive[CReactiveEulerVariable::P_INDEX_LIM] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        Primitive[CReactiveEulerVariable::VX_INDEX_LIM + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

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

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(eps2);

      for(iVar = 0; iVar < nPrimVarLim; ++iVar) {
        eps1 = LimK*dave;
        eps2 = eps1*eps1*eps1;

        //AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        if(dm > 0.0)
          dp = node[iPoint]->GetSolution_Max(iVar);
        else
          dp = node[iPoint]->GetSolution_Min(iVar);

        limiter = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if(limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          //AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }

        /*-- Repeat for point j on the edge ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if(dm > 0.0)
          dp = node[jPoint]->GetSolution_Max(iVar);
        else
          dp = node[jPoint]->GetSolution_Min(iVar);

        limiter = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

        if(limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
        //AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }
      } /*--- End of iVar for loop ---*/
      //AD::EndPreacc();
    }
  }

  /*--- Limiter MPI ---*/
  Set_MPI_Primitive_Limiter(geometry, config);
}

//
//
/*!
 *\brief Call MPI to set solution for parallel simulations
 */
//
//
void CReactiveEulerSolver::Set_MPI_Solution(CGeometry* geometry, CConfig* config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  #ifdef HAVE_MPI
	 int send_to, receive_from;
   SU2_MPI::Status status;
  #endif

	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		if(config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE && config->GetMarker_All_SendRecv(iMarker) > 0) {
			MarkerS = iMarker;
      MarkerR = iMarker+1;

      #ifdef HAVE_MPI
        send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;
			  receive_from = std::abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;
      #endif

		  nVertexS = geometry->nVertex[MarkerS];
      nBufferS_Vector = nVertexS*nVar;
      nVertexR = geometry->nVertex[MarkerR];
		  nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Send_U.resize(nBufferS_Vector);

      /*--- Copy the solution that should be sended ---*/
      for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for(iVar = 0; iVar < nVar; ++iVar)
          Buffer_Send_U[iVar*nVertexS + iVertex] = node[iPoint]->GetSolution(iVar);
      }

      #ifdef HAVE_MPI
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
          Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL] = rotMatrix[0][0]*
                                                              Buffer_Receive_U[CReactiveEulerVariable::RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                                              rotMatrix[0][1]*
                                                              Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 1)*nVertexR +
                                                                                iVertex];

          Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*
                                                                  Buffer_Receive_U[CReactiveEulerVariable::RHOVX_INDEX_SOL*nVertexR +
                                                                                  iVertex] +
                                                                  rotMatrix[1][1]*
                                                                  Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL+1)*nVertexR +
                                                                                    iVertex];
        }
        else {
          Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL] = rotMatrix[0][0]*
                                                              Buffer_Receive_U[CReactiveEulerVariable::RHOVX_INDEX_SOL*nVertexR + iVertex] +
                                                              rotMatrix[0][1]*
                                                              Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 1)*nVertexR +
                                                                                iVertex] +
                                                              rotMatrix[0][2]*
                                                              Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 2)*nVertexR +
                                                                                iVertex];

         Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + 1] = rotMatrix[1][0]*
                                                                 Buffer_Receive_U[CReactiveEulerVariable::RHOVX_INDEX_SOL*nVertexR +
                                                                                  iVertex] +
                                                                 rotMatrix[1][1]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 1)*nVertexR +
                                                                                   iVertex] +
                                                                 rotMatrix[1][2]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 2)*nVertexR +
                                                                                   iVertex];

         Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + 2] = rotMatrix[2][0]*
                                                                 Buffer_Receive_U[CReactiveEulerVariable::RHOVX_INDEX_SOL*nVertexR +
                                                                                  iVertex] +
                                                                 rotMatrix[2][1]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 1)*nVertexR +
                                                                                   iVertex] +
                                                                 rotMatrix[2][2]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::RHOVX_INDEX_SOL + 2)*nVertexR +
                                                                                   iVertex];
        }

        /*--- Copy transformed conserved variables back into buffer. ---*/
        node[iPoint]->SetSolution(Solution);
      } /*--- End of iVertex for loop ---*/
    }
	} /*--- End of iMarker for loop ---*/
}

//
//
/*!
 *\brief Call MPI to set limiter for parallel simulations
 */
//
//
void CReactiveEulerSolver::Set_MPI_Primitive_Limiter(CGeometry* geometry, CConfig* config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint;
  unsigned long nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};

  /*--- Declare vector to store the limiter ---*/
  su2double Limiter[nPrimVarLim];

  #ifdef HAVE_MPI
    int send_to, receive_from;
    SU2_MPI::Status status;
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
        nBufferR_Vector = nVertexR*nPrimVarLim;

        /*--- Allocate Receive and send buffers  ---*/
        Buffer_Send_U.resize(nBufferS_Vector);

        /*--- Copy the solution old that should be sended ---*/
        for(iVertex = 0; iVertex < nVertexS; ++iVertex) {
          iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
          for(iVar = 0; iVar < nPrimVarLim; ++iVar)
            Buffer_Send_U[iVar*nVertexS + iVertex] = node[iPoint]->GetLimiter_Primitive(iVar);
        }

        #ifdef HAVE_MPI
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
            Limiter[CReactiveEulerVariable::VX_INDEX_LIM] =  rotMatrix[0][0]*
                                                             Buffer_Receive_U[CReactiveEulerVariable::VX_INDEX_GRAD*nVertexR + iVertex] +
                                                             rotMatrix[0][1]*
                                                             Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 1)*nVertexR +
                                                                                   iVertex];

            Limiter[CReactiveEulerVariable::VX_INDEX_LIM + 1] = rotMatrix[1][0]*
                                                                Buffer_Receive_U[CReactiveEulerVariable::VX_INDEX_GRAD*nVertexR +
                                                                                      iVertex] +
                                                                rotMatrix[1][1]*
                                                                Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 1)*nVertexR +
                                                                                       iVertex];
          }

          else {
            Limiter[CReactiveEulerVariable::VX_INDEX_LIM] =  rotMatrix[0][0]*
                                                             Buffer_Receive_U[CReactiveEulerVariable::VX_INDEX_GRAD*nVertexR + iVertex] +
                                                             rotMatrix[0][1]*
                                                             Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 1)*nVertexR +
                                                                                   iVertex] +
                                                             rotMatrix[0][2]*
                                                             Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 2)*nVertexR +
                                                                                   iVertex];

            Limiter[CReactiveEulerVariable::VX_INDEX_LIM + 1] = rotMatrix[1][0]*
                                                                Buffer_Receive_U[CReactiveEulerVariable::VX_INDEX_GRAD*nVertexR +
                                                                                     iVertex] +
                                                                rotMatrix[1][1]*
                                                                Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 1)*nVertexR +
                                                                                      iVertex] +
                                                                rotMatrix[1][2]*
                                                                Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 2)*nVertexR +
                                                                                      iVertex];

            Limiter[CReactiveEulerVariable::VX_INDEX_LIM + 2] =  rotMatrix[2][0]*
                                                                 Buffer_Receive_U[CReactiveEulerVariable::VX_INDEX_GRAD*nVertexR +
                                                                                      iVertex] +
                                                                 rotMatrix[2][1]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 1)*nVertexR +
                                                                                       iVertex] +
                                                                 rotMatrix[2][2]*
                                                                 Buffer_Receive_U[(CReactiveEulerVariable::VX_INDEX_GRAD + 2)*nVertexR +
                                                                                       iVertex];
        }

        /*--- Copy transformed limited variables back into buffer. ---*/
        for(iVar = 0; iVar < nPrimVarLim; ++iVar)
          node[iPoint]->SetLimiter_Primitive(iVar,Limiter[iVar]);

      } /*--- End if iVertex for loop ---*/
    }
  } /*--- End if iMarker for loop ---*/
}

//
//
/*!
 *\brief Call MPI to set gradient of primitive variables for parallel simulations
 */
//
//
void CReactiveEulerSolver::Set_MPI_Primitive_Gradient(CGeometry* geometry, CConfig* config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;

  su2double theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi;
  su2double rotMatrix[3][3] = {};
  su2double Gradient[nPrimVarGrad][nDim];

  #ifdef HAVE_MPI
    int send_to, receive_from;
    SU2_MPI::Status status;
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
        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;

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
/*!
 *\brief Setting time step
 */
//
//
void CReactiveEulerSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned long Iteration) {
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker;

  su2double Area, Volume, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda;
  su2double Local_Delta_Time, Global_Delta_Time = 1E6;
  su2double Normal[nDim], GridVel[nDim], GridVel_i[nDim], GridVel_j[nDim];

  bool time_steping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  Min_Delta_Time = 1.E6;
  Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero and compute sound speed ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    node[iPoint]->SetMax_Lambda_Inv(0.0);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5*(node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5*(node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());

    /*--- Adjustment for grid movement ---*/
    if(grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      su2double ProjVel_i = std::inner_product(GridVel_i,GridVel_i + nDim, Normal, 0.0);
      su2double ProjVel_j = std::inner_product(GridVel_j,GridVel_j + nDim, Normal, 0.0);
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
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed();

        /*--- Adjustment for grid movement ---*/
        if(grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          su2double ProjVel = std::inner_product(GridVel,GridVel + nDim, Normal, 0.0);
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
      Local_Delta_Time = config->GetCFL(iMesh)*Volume / node[iPoint]->GetMax_Lambda_Inv();
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
      config->SetCFL(iMesh,config->GetUnst_CFL());

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
/*!
 *\brief Set residual for dual time strategies
 */
//
//
void CReactiveEulerSolver::SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                                unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) {
  SU2_Assert(Residual != NULL, "The array for residual has not been allocated");

  /*--- Local variables ---*/
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;

  su2double Volume_nM1, Volume_nP1, TimeStep, Residual_GCL;
  su2double Normal[nDim], GridVel[nDim], GridVel_i[nDim], GridVel_j[nDim];

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
      if(implicit)
        throw Common::NotImplemented("Implicit computation for dual time not implemented");
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
      Normal = geometry->edge[iEdge]->GetNormal();

      /*--- Grid velocities stored at nodes i & j ---*/
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();

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
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          /*--- Grid velocities stored at boundary node i ---*/
          GridVel = geometry->node[iPoint]->GetGridVel();

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
      if(implicit)
        throw Common::NotImplemented("Implicit computation for upwind residual not implemented");
    } /*--- End loop over all nodes ---*/
  } /*--- End of dynamic mesh part ---*/
}


//
//
/*!
 *\brief Iteration of implicit Euler method
 */
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
      Delta = Volume / node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      std::fill(Local_Res_TruncError,Local_Res_TruncError + nVar,0.0);
      for(iVar = 0; iVar < nVar; ++iVar) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
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
			node[iPoint]->AddClippedSolution(iVar, LinSysSol[iPoint*nVar+iVar], Lower_Limit[iVar], Upper_Limit[iVar]);
  }

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*!
 *\brief Iteration of explicit Euler method
 */
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
      Delta = node[iPoint]->GetDelta_Time() / Volume;

    auto Res_TruncError = node[iPoint]->GetResTruncError();
    auto Residual = LinSysRes.GetBlock(iPoint);

    for(iVar = 0; iVar < nVar; ++iVar) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddClippedSolution(iVar, -Res*Delta, Lower_Limit[iVar], Upper_Limit[iVar]);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, std::abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*!
 *\brief Iteration of explicit RK method
 */
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
      Delta = node[iPoint]->GetDelta_Time() / Volume;

    auto Res_TruncError = node[iPoint]->GetResTruncError();
    auto Residual = LinSysRes.GetBlock(iPoint);

    for(iVar = 0; iVar < nVar; ++iVar) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddClippedSolution(iVar, -Res*Delta*RK_AlphaCoeff, Lower_Limit[iVar], Upper_Limit[iVar]);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, std::abs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
}

//
//
/*!
 *\brief Centered residual convective term
 */
//
//
void CReactiveEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                             CConfig *config, unsigned short iMesh, unsigned short iRKStep) {

  throw Common::NotImplemented("Function not implemented: Centered Residual\n");
}

//
//
/*!
 *\brief Upwind residual convective term
 */
//
//
void CReactiveEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                           CConfig* config, unsigned short iMesh) {
  unsigned long iEdge, iPoint, jPoint, counter_local = 0, counter_global = 0;
  unsigned short iDim, iVar;

  su2double Project_Grad_i, Project_Grad_j;

  /*--- Initialize the upwind convective residual to zero ---*/
  SU2_Assert(Res_Conv != NULL,"The array for source residual has not been allocated");
  std::fill(Res_Conv, Res_Conv + nVar, 0.0);

  /*--- Loop over all the edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Get primitive variables ---*/
    auto Primitive_i = node[iPoint]->GetPrimitive();
    auto Primitive_j = node[jPoint]->GetPrimitive();

    if(second_order) {
      SU2_Assert(Vector_i != NULL, "The array to compute the vector from i to j has not been allocated");
      SU2_Assert(Vector_j != NULL, "The array to compute the vector from j to i has not been allocated");
      for(iDim = 0; iDim < nDim; ++iDim) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = - Vector_i[iDim];
      }
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[iPoint]->GetGradient_Primitive();

      /*--- Auxiliary vector to store reconstruction ---*/
      su2double Prim_Recon_i[nPrimVarLim];
      su2double Prim_Recon_j[nPrimVarLim];

      Prim_Recon_i[CReactiveEulerVariable::T_INDEX_LIM] = Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM];
      Prim_Recon_j[CReactiveEulerVariable::T_INDEX_LIM] = Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM];
      Prim_Recon_i[CReactiveEulerVariable::P_INDEX_LIM] = Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM];
      Prim_Recon_j[CReactiveEulerVariable::P_INDEX_LIM] = Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM];
      std::copy(Primitive_i + CReactiveEulerVariable::VX_INDEX_PRIM, Primitive_i + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                Prim_Recon_i + CReactiveEulerVariable::VX_INDEX_LIM);
      std::copy(Primitive_j + CReactiveEulerVariable::VX_INDEX_PRIM, Primitive_j + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                Prim_Recon_j + CReactiveEulerVariable::VX_INDEX_LIM);

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
      if(Prim_Recon_i[CReactiveEulerVariable::T_INDEX_LIM] > EPS)
        Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM] = Prim_Recon_i[CReactiveEulerVariable::T_INDEX_LIM];
      else
        non_phys_i = true;
      if(!non_phys_i) {
        if(Prim_Recon_i[CReactiveEulerVariable::P_INDEX_LIM] > EPS)
          Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM] = Prim_Recon_i[CReactiveEulerVariable::P_INDEX_LIM];
        else
          non_phys_i = true;
      }

      bool non_phys_j = false;
      if(Prim_Recon_j[CReactiveEulerVariable::T_INDEX_LIM] > EPS)
        Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM] = Prim_Recon_j[CReactiveEulerVariable::T_INDEX_LIM];
      else
        non_phys_j = false;
      if(!non_phys_j) {
        if(Prim_Recon_i[CReactiveEulerVariable::P_INDEX_LIM] > EPS)
          Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM] = Prim_Recon_j[CReactiveEulerVariable::P_INDEX_LIM];
        else
          non_phys_j = false;
      }

      /*--- Compute other primitive variables accordingly to the reconstruction ---*/
      if(!non_phys_i) {
        su2double rho_recon_i;
        std::copy(Primitive_i + CReactiveEulerVariable::RHOS_INDEX_PRIM,
                  Primitive_i + (CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies), Ys_i.begin());
        //rho_recon_i = library->ComputeDensity(Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM],
        //                                      Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM],Ys_i);
        rho_recon_i *= config->GetGas_Constant_Ref();
        Primitive_i[CReactiveEulerVariable::RHO_INDEX_PRIM] = rho_recon_i;
        su2double dim_temp_i;
        bool US_System = (config->GetSystemMeasurements() == US);
        dim_temp_i = Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp_i *= 5.0/9.0;
        //Primitive_i[CReactiveEulerVariable::H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp_i,Ys_i);
        Primitive_i[CReactiveEulerVariable::H_INDEX_PRIM] /= config->GetEnergy_Ref();
        if(US_System)
          Primitive_i[CReactiveEulerVariable::H_INDEX_PRIM] *= 3.28084*3.28084;
        su2double sq_vel_i = std::inner_product(Primitive_i + CReactiveEulerVariable::VX_INDEX_PRIM,
                                                Primitive_i + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                                                Primitive_i + CReactiveEulerVariable::VX_INDEX_PRIM, 0.0);
        Primitive_i[CReactiveEulerVariable::H_INDEX_PRIM] += 0.5*sq_vel_i;
        //Gamma_i = library->ComputeFrozenGamma(dim_temp_i,Ys_i);
        //Primitive_i[CReactiveEulerVariable::A_INDEX_PRIM] = std::sqrt(Gamma_i*Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM]/rho_recon_i);
      }

      if(!non_phys_j) {
        su2double rho_recon_j;
        std::copy(Primitive_j + CReactiveEulerVariable::RHOS_INDEX_PRIM,
                  Primitive_j + (CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies), Ys_j.begin());
        //rho_recon_j = library->ComputeDensity(Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM],
        //                                      Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM],Ys_j);
        rho_recon_j *= config->GetGas_Constant_Ref();
        Primitive_j[CReactiveEulerVariable::RHO_INDEX_PRIM] = rho_recon_j;
        su2double dim_temp_j;
        bool US_System = (config->GetSystemMeasurements() == US);
        dim_temp_j = Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp_j *= 5.0/9.0;
        //Primitive_j[CReactiveEulerVariable::H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp_j,Ys_j);
        Primitive_j[CReactiveEulerVariable::H_INDEX_PRIM] /= config->GetEnergy_Ref();
        if(US_System)
          Primitive_j[CReactiveEulerVariable::H_INDEX_PRIM] *= 3.28084*3.28084;
        su2double sq_vel_j = std::inner_product(Primitive_j + CReactiveEulerVariable::VX_INDEX_PRIM,
                                                Primitive_j + CReactiveEulerVariable::VX_INDEX_PRIM + nDim,
                                                Primitive_j + CReactiveEulerVariable::VX_INDEX_PRIM, 0.0);
        Primitive_j[CReactiveEulerVariable::H_INDEX_PRIM] += 0.5*sq_vel_j;
        //Gamma_j = library->ComputeFrozenGamma(dim_temp_j,Ys_j);
        //Primitive_j[CReactiveEulerVariable::A_INDEX_PRIM] = std::sqrt(Gamma_j*Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM]/rho_recon_j);
      }

      numerics->SetPrimitive(Primitive_i, Primitive_j);
    } /*--- End second order reconstruction ---*/

    else {
      /*--- Set primitive variables without reconstruction ---*/
      numerics->SetPrimitive(Primitive_i, Primitive_j);
    }

    /*--- Compute the residual ---*/
    try {
      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      dynamic_cast<CUpwReactiveAUSM*>(numerics)->SetExplicit();
      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err;
    err = !std::none_of(Res_Conv,Res_Conv + nVar,[](su2double elem){return std::isnan(elem);});
    if(implicit)
      throw Common::NotImplemented("Implicit computation for upwind residual not implemented");

    /*--- Update residual value ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);

      /*--- Set implicit Jacobians ---*/
      if(implicit)
        throw Common::NotImplemented("Implicit computation for upwind residual not implemented");
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
/*!
 *\brief Residual source term
 */
//
//
void CReactiveEulerSolver::Source_Residual(CGeometry* geometry, CSolver** solver_container,
                                           CNumerics* numerics, CNumerics* second_numerics,
                                           CConfig* config, unsigned short iMesh) {
  /*--- Local variables ---*/
  unsigned long iPoint, counter_local = 0, counter_global=0;

  /*--- Initialize the source residual to zero ---*/
  SU2_Assert(Res_Sour != NULL,"The array for source residual has not been allocated");
  std::fill(Res_Sour, Res_Sour + nVar, 0.0);

  /*--- Loop over all points ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Load the primitive variables ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(),node[iPoint]->GetPrimitive());

    /*--- Load the volume of the dual mesh cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Compute the source residual ---*/
    try {
      numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      dynamic_cast<CSourceReactive*>(numerics)->SetExplicit();
      numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err;
    err = !std::none_of(Res_Sour,Res_Sour + nVar,[](su2double elem){return std::isnan(elem);});

    /*--- Add the source residual to the total ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Add the implicit Jacobian contribution ---*/
      if(implicit)
        throw Common::NotImplemented("Implicit computation for source residual not implemented");
    }
    else
      throw std::runtime_error("NaN found in the source residual");
  }

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
/*!
  *\brief Set freestream solution all over the domain
  */
//
//
void CReactiveEulerSolver::SetFreeStream_Solution(CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint;
  unsigned short iDim, iSpecies;

  for(iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Set density ---*/
    node[iPoint]->SetSolution(CReactiveEulerVariable::RHO_INDEX_SOL, Density_Inf);
    /*--- Set momentum ---*/
    for(iDim = 0; iDim < nDim; ++iDim)
      node[iPoint]->SetSolution(CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim, Density_Inf*Velocity_Inf[iDim]);
    /*--- Set energy ---*/
    node[iPoint]->SetSolution(CReactiveEulerVariable::RHOE_INDEX_SOL, Density_Inf*config->GetEnergy_FreeStreamND());
    /*--- Set partial densities ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      node[iPoint]->SetSolution(CReactiveEulerVariable::RHOS_INDEX_SOL, Density_Inf*MassFrac_Inf[iSpecies]);
  }
}


//
//
/*!
  *\brief Solid wall boundary conditions for Euler
  */
//
//
void CReactiveEulerSolver::BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                         CConfig* config,unsigned short val_marker) {
  SU2_Assert(Residual != NULL,"The array to save residual for boundary conditions has not been allocated");

  unsigned long iVertex,iPoint;
  unsigned short iDim;

  su2double Area,Pressure;
  su2double Normal[nDim], UnitNormal[nDim], GridVel[nDim];

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    std::fill(Residual,Residual + nVar,0.0);
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if(geometry->node[iPoint]->GetDomain()) {
			/*--- Normal vector for this vertex (negative for outward convention) ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

			/*--- Compute parameters from the geometry ---*/
      Area = ::ComputeArea(Normal,nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			/*--- Retrieve the pressure on the vertex ---*/
      Pressure = node[iPoint]->GetPressure();

			/*--- Compute the residual ---*/
			for(iDim = 0; iDim < nDim; ++iDim)
				Residual[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = Pressure*UnitNormal[iDim]*Area;

      /*--- Correction to energy due to mowing wall: extra term due to pressure contribution ---*/
      if(grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        su2double ProjGridVel = std::inner_product(GridVel, GridVel + nDim, UnitNormal, 0.0);
        Residual[CReactiveEulerVariable::RHOE_INDEX_SOL] += Pressure*ProjGridVel*Area;
      }

      /*--- Add value to the residual ---*/
  		LinSysRes.AddBlock(iPoint, Residual);

      if(implicit)
        throw Common::NotImplemented("Implicit computation for solid wall boundary conditions for Euler not implemented");
    }
  } /*--- End of iVertex for loop ---*/
}

//
//
/*!
  *\brief Supersonic inlet boundary conditions
  */
//
//
void CReactiveEulerSolver::BC_Supersonic_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                               CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
	unsigned long iVertex, iPoint, Point_Normal;

  su2double V_inlet[nPrimVar], V_domain[nPrimVar];
  su2double Normal[nDim];

  std::string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool viscous = config->GetViscous();

  su2double Density, Pressure, Temperature, Enthalpy, SoundSpeed;

  /*--- Supersonic inlet flow: there are no outgoing characteristics,
        so all flow variables can be imposed at the inlet.
        First, retrieve the specified values for the primitive variables. ---*/
	Temperature = config->GetInlet_Temperature(Marker_Tag);
	Pressure    = config->GetInlet_Pressure(Marker_Tag);
	auto Velocity = config->GetInlet_Velocity(Marker_Tag);
  Ys = config->GetInlet_MassFrac(Marker_Tag);

	/*--- Density at the inlet from the gas law ---*/
	//Density = library->ComputeDensity(Temperature,Pressure);

  /*--- Compute the energy from the specified state ---*/
  bool US_System = (config->GetSystemMeasurements() == US);
  su2double dim_temp = Temperature;
  if(US_System)
    dim_temp *= 5.0/9.0;
  //Enthalpy = library->ComputeEnthalpy(dim_temp,MassFrac);
  //SoundSpeed = library->ComputeFrozenSoundSpeed(dim_temp,MassFrac);
  //if(US_System) {
  //  Enthalpy *= 3.28084*3.28084;
  //  SoundSpeed *= 3.28084;
  //}

  /*--- Non-dim. the inputs if necessary. ---*/
	Temperature /= config->GetTemperature_Ref();
	Pressure /= config->GetPressure_Ref();
	Density /= config->GetDensity_Ref();
  Velocity /= config->GetVelocity_Ref();
  //Enthalpy /= config->GetEnergy_Ref();
  //SoundSpeed /= config->GetVelocity_Ref();

  /*--- Primitive variables using the derived quantities ---*/
  V_inlet[CReactiveEulerVariable::T_INDEX_PRIM] = Temperature;
  std::copy(Velocity, Velocity + nDim, V_inlet + CReactiveEulerVariable::VX_INDEX_PRIM);
  V_inlet[CReactiveEulerVariable::P_INDEX_PRIM] = Pressure;
  V_inlet[CReactiveEulerVariable::RHO_INDEX_PRIM] = Density;
  su2double Velocity2 = std::inner_product(Velocity, Velocity + nDim, Velocity, 0.0);
  V_inlet[CReactiveEulerVariable::H_INDEX_PRIM] = Enthalpy + 0.5*Velocity2;
  V_inlet[CReactiveEulerVariable::A_INDEX_PRIM] = SoundSpeed;
  std::copy(Ys.cbegin(), Ys.cend(), V_inlet + CReactiveEulerVariable::RHOS_INDEX_PRIM);

  /*--- Loop over all the vertices on this boundary (val_marker) ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if(geometry->node[iPoint]->GetDomain()) {
       /*--- Normal vector for this vertex (negative for outward convention) ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      std::transform(Normal,Normal + nDim,Normal,std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      /*--- Pull primitve variable from actual point ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Pull primitve variable from actual point ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      try {
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }
      catch(const std::exception& e) {
        std::cout<<e.what()<<std::endl;
        dynamic_cast<CUpwReactiveAUSM*>(conv_numerics)->SetExplicit();
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }
      LinSysRes.AddBlock(iPoint, Residual);

      if(implicit)
        throw Common::NotImplemented("Implicit computation for supersonic inlet BC not implemented");

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

          /*--- Laminar viscosity ---*/
          visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

          /*--- Thermal conductivity ---*/
          visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

          /*--- Species binary coefficients ---*/
          auto avggrad_numerics = dynamic_cast<CAvgGradReactive_Flow*>(visc_numerics);
          SU2_Assert(avggrad_numerics != NULL, "The cast to compute the viscous flux has not been successfull");
          auto ns_variable = dynamic_cast<CReactiveNSVariable*>(node[iPoint]);
          SU2_Assert(ns_variable != NULL, "The cast compute the binary diffusion coefficients has not been successfull");
          avggrad_numerics->SetBinaryDiffCoeff(ns_variable->GetBinaryDiffusionCoeff(), ns_variable->GetBinaryDiffusionCoeff());
          try {
            visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
          }
          catch(const std::exception& e) {
            std::cout<<e.what()<<std::endl;
            avggrad_numerics->SetExplicit();
            visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
          }
          LinSysRes.SubtractBlock(iPoint, Residual);

  				/*--- Jacobian contribution for implicit integration ---*/
  				if(implicit)
  					throw Common::NotImplemented("Implicit computation for supersonic inlet BC not implemented");
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
void CReactiveEulerSolver::BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned short iDim;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double P_Total, T_Total, Velocity2, H_Total, Temperature, Riemann, Pressure, Density,
            Energy, Mach2, SoundSpeed2, SoundSpeed_Total2, Vel_Mag, Gamma, Gamma_Minus_One,
            alpha, aa, bb, cc, dd, Area;
  su2double Normal[nDim], UnitNormal[nDim], Velocity[nDim], Flow_Dir[nDim];
  su2double V_inlet[nPrimVar], V_domain[nPrimVar];

  unsigned short Kind_Inlet = config->GetKind_Inlet();
  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);
  bool viscous              = config->GetViscous();
  bool gravity              = config->GetGravityForce();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = ::ComputeArea(Normal,nDim)
      for(iDim = 0; iDim < nDim; iDim++) {
        Normal[iDim] = -Normal[iDim];
        UnitNormal[iDim] = Normal[iDim]/Area;
      }
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Build the fictitious intlet state based on characteristics ---*/
      /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
            therefore we can specify all but one state variable at the inlet.
            The outgoing Riemann invariant provides the final piece of info. ---*/
      switch (Kind_Inlet) {
        /*--- Total properties have been specified at the inlet. ---*/
        case TOTAL_CONDITIONS:
          /*--- Retrieve the specified total conditions for this inlet. ---*/
          P_Total = config->GetInlet_Ptotal(Marker_Tag);
          if(gravity)
            P_Total -= geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;

          T_Total  = config->GetInlet_Ttotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/
          P_Total /= config->GetPressure_Ref();
          T_Total /= config->GetTemperature_Ref();

          /*--- Store primitives and set some variables for clarity. ---*/
          Density = V_domain[CReactiveEulerVariable::RHO_INDEX_PRIM];
          Velocity2 = 0.0;
          for(iDim = 0; iDim < nDim; iDim++) {
            Velocity[iDim] = V_domain[CReactiveEulerVariable::VX_INDEX_PRIM + iDim];
            Velocity2 += Velocity[iDim]*Velocity[iDim];
          }
          Pressure    = V_domain[nDim+1];
          Energy      = V_domain[CReactiveEulerVariable::H_INDEX_PRIM] - Pressure/Density;
          //H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
          SoundSpeed2 = Gamma*Pressure/Density;

          /*--- Compute the acoustic Riemann invariant that is extrapolated
                from the domain interior. ---*/

          Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
          for(iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Total speed of sound ---*/

          SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

          /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/

          alpha = 0.0;
          for(iDim = 0; iDim < nDim; iDim++)
            alpha += UnitNormal[iDim]*Flow_Dir[iDim];

          /*--- Coefficients in the quadratic equation for the velocity ---*/

          aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
          bb = -1.0*Gamma_Minus_One*alpha*Riemann;
          cc =  0.5*Gamma_Minus_One*Riemann*Riemann
              -2.0*SoundSpeed_Total2/Gamma_Minus_One;

          /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/

          dd = bb*bb - 4.0*aa*cc;
          dd = sqrt(max(0.0, dd));
          Vel_Mag   = (-bb + dd)/(2.0*aa);
          Vel_Mag   = max(0.0, Vel_Mag);
          Velocity2 = Vel_Mag*Vel_Mag;

          /*--- Compute speed of sound from total speed of sound eqn. ---*/

          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/

          Mach2 = Velocity2/SoundSpeed2;
          Mach2 = min(1.0, Mach2);
          Velocity2   = Mach2*SoundSpeed2;
          Vel_Mag     = sqrt(Velocity2);
          SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

          /*--- Compute new velocity vector at the inlet ---*/

          for(iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

          /*--- Static temperature from the speed of sound relation ---*/

          Temperature = SoundSpeed2/(Gamma*Gas_Constant);

          /*--- Static pressure using isentropic relation at a point ---*/

          Pressure = P_Total*pow((Temperature/T_Total), Gamma/Gamma_Minus_One);

          /*--- Density at the inlet from the gas law ---*/

          Density = Pressure/(Gas_Constant*Temperature);

          /*--- Using pressure, density, & velocity, compute the energy ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Velocity2;
          if(tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Temperature;
          for(iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Velocity[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;

          /*--- Mass flow has been specified at the inlet. ---*/

        case MASS_FLOW:

          /*--- Retrieve the specified mass flow for the inlet. ---*/

          Density  = config->GetInlet_Ttotal(Marker_Tag);
          Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
          Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

          /*--- Non-dim. the inputs if necessary. ---*/

          Density /= config->GetDensity_Ref();
          Vel_Mag /= config->GetVelocity_Ref();

          /*--- Get primitives from current inlet state. ---*/

          for(iDim = 0; iDim < nDim; iDim++)
            Velocity[iDim] = node[iPoint]->GetVelocity(iDim);
          Pressure    = node[iPoint]->GetPressure();
          SoundSpeed2 = Gamma*Pressure/V_domain[nDim+2];

          /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/

          Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
          for(iDim = 0; iDim < nDim; iDim++)
            Riemann += Velocity[iDim]*UnitNormal[iDim];

          /*--- Speed of sound squared for fictitious inlet state ---*/

          SoundSpeed2 = Riemann;
          for(iDim = 0; iDim < nDim; iDim++)
            SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitNormal[iDim];

          SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
          SoundSpeed2 = SoundSpeed2*SoundSpeed2;

          /*--- Pressure for the fictitious inlet state ---*/

          Pressure = SoundSpeed2*Density/Gamma;

          /*--- Energy for the fictitious inlet state ---*/

          Energy = Pressure/(Density*Gamma_Minus_One) + 0.5*Vel_Mag*Vel_Mag;
          if(tkeNeeded) Energy += GetTke_Inf();

          /*--- Primitive variables, using the derived quantities ---*/

          V_inlet[0] = Pressure / ( Gas_Constant * Density);
          for(iDim = 0; iDim < nDim; iDim++)
            V_inlet[iDim+1] = Vel_Mag*Flow_Dir[iDim];
          V_inlet[nDim+1] = Pressure;
          V_inlet[nDim+2] = Density;
          V_inlet[nDim+3] = Energy + Pressure/Density;

          break;
      }

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/

      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if(implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Roe Turkel preconditioning, set the value of beta ---*/

      if(config->GetKind_Upwind() == TURKEL)
        node[iPoint]->SetPreconditioner_Beta(conv_numerics->GetPrecond_Beta());

      /*--- Viscous contribution ---*/

      if(viscous) {

        /*--- Set laminar and eddy viscosity at the infinity ---*/

        V_inlet[nDim+5] = node[iPoint]->GetLaminarViscosity();
        V_inlet[nDim+6] = node[iPoint]->GetEddyViscosity();

        /*--- Set the normal vector and the coordinates ---*/

        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());

        /*--- Primitive variables, and gradient ---*/

        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[iPoint]->GetGradient_Primitive());

        /*--- Turbulent kinetic energy ---*/

        if(config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));

        /*--- Compute and update residual ---*/

        visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        LinSysRes.SubtractBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/

        if(implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

      }

    }
  }

  /*--- Free locally allocated memory ---*/

  delete [] Normal;

}

//
//
/*!
  *\brief Outlet boundary conditions (supersonic or subosnic according to specific simulation)
  */
//
//
void CReactiveEulerSolver::BC_Supersonic_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                                CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Normal[nDim];
  su2double V_outlet[nPrimVar], V_domain[nPrimVar];

  std::string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  bool viscous = config->GetViscous();
  bool gravity = config->GetGravityForce();

  /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if(geometry->node[iPoint]->GetDomain()) {
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      std::transform(Normal,Normal + nDim,Normal,std::negate<su2double>());
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Supersonic exit flow: there are no incoming characteristics,
            so no boundary condition is necessary. Set outlet state to current
            state so that upwinding handles the direction of propagation. ---*/
      std::copy(V_domain, V_domain + nPrimVar, V_outlet);

      /*--- Set primitives and grid movement in the convective numerics class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      try {
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }
      catch(const std::exception& e) {
        std::cout<<e.what()<<std::endl;
        dynamic_cast<CUpwReactiveAUSM*>(conv_numerics)->SetExplicit();
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if(implicit)
        throw Common::NotImplemented("Implicit computation for outlet BC not implemented");

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

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        auto avggrad_numerics = dynamic_cast<CAvgGradReactive_Flow*>(visc_numerics);
        SU2_Assert(avggrad_numerics != NULL, "The cast to compute the viscous flux has not been successfull");
        auto ns_variable = dynamic_cast<CReactiveNSVariable*>(node[iPoint]);
        SU2_Assert(ns_variable != NULL, "The cast compute the binary diffusion coefficients has not been successfull");
        avggrad_numerics->SetBinaryDiffCoeff(ns_variable->GetBinaryDiffusionCoeff(), ns_variable->GetBinaryDiffusionCoeff());

        /*--- Compute and update residual ---*/
        try {
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        }
        catch(const std::exception& e) {
          std::cout<<e.what()<<std::endl;
          avggrad_numerics->SetExplicit();
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        }
        LinSysRes.SubtractBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          throw Common::NotImplemented("Implicit computation for supersonic inlet BC not implemented");
      }
    } /*--- End of if GetDomain() ---*/
  } /*--- End of iVertex for loop ---*/
}


//
//
/*!
  *\brief Outlet boundary conditions (supersonic or subosnic according to specific simulation)
  */
//
//
void CReactiveEulerSolver::BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                                     CConfig* config, unsigned short val_marker) {
  unsigned short iVar, iDim, jDim, iSpecies;
  unsigned long iVertex, iPoint, Point_Normal;

  su2double Pressure, P_Exit, Velocity2, Entropy, Density, Energy, dim_temp,
            Riemann, Vn, SoundSpeed, Gamma, Gamma_Minus_One, Mach_Exit, Vn_Exit, Area;
  su2double V_outlet[nPrimVar], V_domain[nPrimVar];
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
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = ::ComputeArea(Normal,nDim);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Normal[iDim] = -Normal[iDim];
        UnitNormal[iDim] = Normal[iDim]/Area;
      }
      conv_numerics->SetNormal(Normal);

      /*--- Current solution at this boundary node ---*/
      V_domain = node[iPoint]->GetPrimitive();

      /*--- Build the fictitious intlet state based on characteristics ---*/
      /*--- Retrieve the specified back pressure for this outlet. ---*/
      P_Exit = config->GetOutlet_Pressure(Marker_Tag);
      if(gravity)
        P_Exit -= geometry->node[iPoint]->GetCoord(nDim-1)*STANDART_GRAVITY;

      /*--- Non-dim. the inputs if necessary. ---*/
      P_Exit /= config->GetPressure_Ref();

      /*--- Check whether the flow is supersonic at the exit. The type
            of boundary update depends on this. ---*/
      Density = V_domain[CReactiveEulerVariable::RHO_INDEX_PRIM];
      std::copy(V_domain + CReactiveEulerVariable::VX_INDEX_PRIM,
                V_domain + CReactiveEulerVariable::VX_INDEX_PRIM + nDim, Velocity);
      Velocity2 = std::inner_product(Velocity, Velocity + nDim, Velocity, 0.0);
      Vn = std::inner_product(UnitNormal, UnitNormal + nDim, Velocity, 0.0);
      Pressure = V_domain[CReactiveEulerVariable::P_INDEX_PRIM];
      std::copy(V_domain + CReactiveEulerVariable::RHOS_INDEX_PRIM,
                V_domain + CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies, Ys.begin());
      bool US_System = (config->GetSystemMeasurements() == US);
      dim_temp = V_domain[CReactiveEulerVariable::T_INDEX_PRIM]*config->GetTemperature_Ref();
      if(US_System)
        dim_temp *= 5.0/9.0;
      //Gamma = library->ComputeFrozenGamma(dim_temp,Ys);
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
        Entropy = Pressure*std::pow(1.0/Density, Gamma);
        Gamma_Minus_One = Gamma - 1.0;
        Riemann = Vn + 2.0*SoundSpeed/Gamma_Minus_One;

        /*--- Compute the new fictious state at the outlet ---*/
        Pressure = P_Exit;

        /*--- Secant method to determine the temperature ---*/
        double a;
        //const double a = Pressure*config->GetGas_Constant_Ref()/library->ComputeRgas(Ys);
        const double b = Pressure/Entropy;
        bool NRconvg, Bconvg;
        su2double NRtol = 1.0E-6;    // Tolerance for the Secant method
        su2double Btol = 1.0E-4;    // Tolerance for the Bisection method
        unsigned short maxNIter = 5;        // Maximum Secant method iterations
        unsigned short maxBIter = 32;        // Maximum Bisection method iterations
        unsigned short iIter;

        su2double T = V_domain[CReactiveEulerVariable::T_INDEX_PRIM];
        su2double Told = T + 1.0;
        su2double Tnew;
        NRconvg = false;
        su2double gamma;
        su2double dim_temp_old, gamma_old;

        /*--- Execute a secant root-finding method to find the outlet temperature ---*/
        for(iIter = 0; iIter < maxNIter; ++iIter) {
          dim_temp = T*config->GetTemperature_Ref();
          dim_temp_old = Told*config->GetTemperature_Ref();
          if(US_System) {
            dim_temp *= 5.0/9.0;
            dim_temp_old *= 5.0/9.0;
          }
          //gamma = library->ComputeFrozenGamma(dim_temp,Ys);
          //gamma_old = library->ComputeFrozenGamma(dim_temp_old,Ys);
          su2double tmp = std::pow(a/T,gamma);
          su2double f = tmp - b;
          su2double df = tmp - std::pow(a/Told,gamma_old);
          Tnew = T - f*(T-Told)/df;

          /*--- Check for convergence ---*/
          if(std::abs(Tnew - T) < NRtol) {
            NRconvg = true;
            break;
          }
          else {
            Told = T;
            T = Tnew;
          }
        }

        if(NRconvg)
          V_outlet[CReactiveEulerVariable::T_INDEX_PRIM] = T;
        else {
          /*--- Execute the bisection root-finding method ---*/
          Bconvg = false;
          su2double Ta = 50;
          su2double Tb = 8.0e4;
          for(iIter = 0; iIter < maxBIter; ++iIter) {
            T = (Ta + Tb)/2.0;
            dim_temp = T*config->GetTemperature_Ref();;
            if(US_System)
              dim_temp *= 5.0/9.0;
            //gamma = library->ComputeFrozenGamma(dim_temp,Ys);

            su2double f = std::pow(a/T,gamma) - b;

            if(std::abs(f) < Btol) {
              V_outlet[CReactiveEulerVariable::T_INDEX_PRIM] = T;
              Bconvg = true;
              break;
            }
            else {
              if(f > 0)
                Ta = T;
              else
                Tb = T;
            }

          /*--- If absolutely no convergence, then something is going really wrong ---*/
          if(!Bconvg)
            throw std::runtime_error("Convergence not achieved for bisection method");
          }
        }

        /*--- Complete the fictious state at the outlet ---*/
        Density  = a/V_outlet[CReactiveEulerVariable::T_INDEX_PRIM];
        dim_temp = V_outlet[CReactiveEulerVariable::T_INDEX_PRIM]*config->GetTemperature_Ref();
        if(US_System)
          dim_temp *= 5.0/9.0;
        //Gamma = library->Frozen_GammaSoundSpeed(dim_temp,Ys);
        Gamma_Minus_One = Gamma - 1.0;
        SoundSpeed = std::sqrt(Gamma*Pressure/Density);
        Vn_Exit = Riemann - 2.0*SoundSpeed/Gamma_Minus_One;
        su2double diff_v = Vn_Exit - Vn;
        for(iDim = 0; iDim < nDim; ++iDim)
          Velocity[iDim] += diff_v*UnitNormal[iDim];
        Velocity2 = std::inner_product(Velocity, Velocity + nDim, Velocity, 0.0);

        /*--- Conservative variables, using the derived quantities ---*/
        std::copy(Velocity, Velocity + nDim, V_outlet + CReactiveEulerVariable::VX_INDEX_PRIM);
        V_outlet[CReactiveEulerVariable::P_INDEX_PRIM] = Pressure;
        V_outlet[CReactiveEulerVariable::RHO_INDEX_PRIM] = Density;
        //V_outlet[CReactiveEulerVariable::H_INDEX_PRIM] = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
        if(US_System)
          V_outlet[CReactiveEulerVariable::H_INDEX_PRIM] *= 3.28084*3.28084;
        V_outlet[CReactiveEulerVariable::H_INDEX_PRIM] += 0.5*Velocity2;
        V_outlet[CReactiveEulerVariable::A_INDEX_PRIM] = SoundSpeed;
        std::copy(Ys.cbegin(), Ys.cend(), V_outlet + CReactiveEulerVariable::RHOS_INDEX_PRIM);
      } /*--- End subsonic outlet preliminary computations ---*/

      /*--- Set primitives and grid movement in the convective numerics class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      if(grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      try {
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }
      catch(const std::exception& e) {
        std::cout<<e.what()<<std::endl;
        dynamic_cast<CUpwReactiveAUSM*>(conv_numerics)->SetExplicit();
        conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      }

      /*--- Update residual value ---*/
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      if(implicit)
        throw Common::NotImplemented("Implicit computation for outlet BC not implemented");

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

        /*--- Laminar viscosity ---*/
        visc_numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[iPoint]->GetLaminarViscosity());

        /*--- Thermal conductivity ---*/
        visc_numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[iPoint]->GetThermalConductivity());

        /*--- Species binary coefficients ---*/
        auto avggrad_numerics = dynamic_cast<CAvgGradReactive_Flow*>(visc_numerics);
        SU2_Assert(avggrad_numerics != NULL, "The cast to compute the viscous flux has not been successfull");
        auto ns_variable = dynamic_cast<CReactiveNSVariable*>(node[iPoint]);
        SU2_Assert(ns_variable != NULL, "The cast compute the binary diffusion coefficients has not been successfull");
        avggrad_numerics->SetBinaryDiffCoeff(ns_variable->GetBinaryDiffusionCoeff(), ns_variable->GetBinaryDiffusionCoeff());

        /*--- Compute and update residual ---*/
        try {
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        }
        catch(const std::exception& e) {
          std::cout<<e.what()<<std::endl;
          avggrad_numerics->SetExplicit();
          visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
        }
        LinSysRes.SubtractBlock(iPoint, Residual);

        /*--- Jacobian contribution for implicit integration ---*/
        if(implicit)
          throw Common::NotImplemented("Implicit computation for supersonic inlet BC not implemented");
      }
    } /*--- End of if GetDomain() ---*/
  } /*--- End of iVertex for loop ---*/
}


//
//
/*!
  *\brief Class CReactiveNSSolver constructor
  */
//
//
CReactiveNSSolver::CReactiveNSSolver(CGeometry* geometry, CConfig* config,unsigned short iMesh): CReactiveEulerSolver() {
  unsigned long iPoint;
  unsigned short iVar,iDim;

  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  IterLinSolver = 0;

  nSpecies = library->GetNSpecies();

  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  std::string file_name = config->GetSolution_FlowFileName();

  if(!(!restart || iMesh != MESH_0)) {

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if(time_stepping) {
      auto Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter());
      file_name = config->GetUnsteady_FileName(file_name, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/
    //Read_SU2_Restart_Metadata(geometry, config, false, file_name);
  }

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho,rho*vx,rho*vy,rho*vz,rho*E,rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T,vx,vy,vz,P,rho,h,a,rho1,...rhoNs)^T ---*/
  nPrimVarLim = nDim + 2;  /*--- Primitive variables to limit (T,vx,vy,vz,P)^T ---*/
  nPrimVarGrad = nPrimVarLim + nSpecies; /*--- Gradient Primitive variables (T,vx,vy,vz,P,X1....XNs)^T ---*/
  //nPrimVarGrad = nPrimVarLim + nSpecies + 1; /*--- Gradient Primitive variables (T,vx,vy,vz,P,rho,Y1....YNs)^T ---*/

  /*--- Perform the non-dimensionalization for the flow equations using the specified reference values. ---*/
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

  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];
  Solution_i = new su2double[nVar];
  Solution_j = new su2double[nVar];

  /*--- Allocate vectors related to the primitive ---*/
  PrimVar_i.resize(nPrimVarGrad);
  PrimVar_j.resize(nPrimVarGrad);
  PrimVar_Vertex.resize(nPrimVarGrad);

  Prim_i.resize(nPrimVarLim);
  Prim_j.resize(nPrimVarLim);
  Primitive.resize(nPrimVarLim);

  /*--- Allocate vectors for conserved variable limits ---*/
  Lower_Limit.resize(nPrimVar, 0.0);
  Upper_Limit.resize(nPrimVar, 1.0/EPS);

  std::fill(Lower_Limit.begin() + CReactiveNSVariable::RHOVX_INDEX_SOL,
            Lower_Limit.begin() + CReactiveNSVariable::RHOVX_INDEX_SOL + nDim, -1.0/EPS);
  Lower_Limit[CReactiveNSVariable::RHOE_INDEX_SOL] = -1.0/EPS;

  /*--- Allocate auxiliary vectors ---*/
  Ys_i.resize(nSpecies);
  Ys_j.resize(nSpecies);
  Ys.resize(nSpecies);

  /*--- Initialize the solution & residual CVectors ---*/
 	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure forstoring extra information ---*/
 	if(config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

	/*--- Allocate Jacobians forimplicit time-stepping ---*/
	if(config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    implicit = true;
		Jacobian_i = new su2double* [nVar];
		Jacobian_j = new su2double* [nVar];
		for(iVar = 0; iVar < nVar; ++iVar) {
			Jacobian_i[iVar] = new su2double [nVar];
			Jacobian_j[iVar] = new su2double [nVar];
		}
  }
  else
    implicit = false;

  /*--- Allocate arrays for gradient computation by least squares ---*/
	if(config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    least_squares = true;

    Smatrix = new su2double*[nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      Smatrix[iDim] = new su2double[nDim];

    Cvector = new su2double*[nPrimVarGrad];
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      Cvector[iVar] = new su2double[nDim];
  }
	else
    least_squares = false;

  grid_movement = config->GetGrid_Movement();
  space_centered = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  second_order = (config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  unsigned long ExtIter = config->GetExtIter();
  limiter = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter());

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(),config->GetVelocity_FreeStreamND() + nDim);
  MassFrac_Inf     = RealVec(config->GetMassFrac_FreeStream(),config->GetMassFrac_FreeStream() + nSpecies);

  Viscosity_Inf    = config->GetViscosity_FreeStreamND();

  node_infty = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
                                       nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, config);

  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  auto direct_diff = config->GetDirectDiff();
  switch(direct_diff) {
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
    break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
    break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
    break;
    default:
    break;
  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  for(iPoint = 0; iPoint < nPoint; ++iPoint)
    node[iPoint] = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
                                           nDim, nVar, nSpecies, nPrimVar, nPrimVarGrad, nPrimVarLim, config);

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	Set_MPI_Solution(geometry, config);
}

//
//
/*!
  *\brief Variables preprocessing
  */
//
//
void CReactiveNSSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
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
    else
      throw std::out_of_range("Unknown option in the computation of the gradient");

    /*--- Limiter computation ---*/
    if(limiter && iMesh == MESH_0 && !Output)
      SetPrimitive_Limiter(geometry, config);
  }

  /*--- Artificial dissipation ---*/
  if(space_centered && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");

  /*--- Initialize the Jacobian matrices ---*/
  if(implicit)
    throw Common::NotImplemented("Implicit scheme not implemented. Switch to explicit");

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
/*!
 *\brief Setting time step
 */
//
//
void CReactiveNSSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                     unsigned short iMesh, unsigned long Iteration) {
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker, iSpecies;

  su2double Area, Volume;
  su2double Local_Delta_Time, Global_Delta_Time = 1E6;
  su2double Local_Delta_Time_Visc;
  su2double Mean_SoundSpeed, Mean_ProjVel;
  su2double Mean_LaminarVisc, Mean_ThermalCond, Mean_Density, Mean_CV;
  su2double Lambda, Lambda_1, Lambda_2;
  su2double Normal[nDim], GridVel[iDim], GridVel_i[nDim], GridVel_j[nDim]

  const su2double K_v = 0.5;
  const su2double FOUR_OVER_THREE = 4.0/3.0;

  bool time_steping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  Min_Delta_Time = 1.E6;
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

    Normal = geometry->edge[iEdge]->GetNormal();
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5*(node[iPoint]->GetProjVel(Normal) + node[jPoint]->GetProjVel(Normal));
    Mean_SoundSpeed = 0.5*(node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());

    /*--- Adjustment for grid movement ---*/
    if(grid_movement) {
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      su2double ProjVel_i = std::inner_product(GridVel_i,GridVel_i + nDim, Normal, 0.0);
      su2double ProjVel_j = std::inner_product(GridVel_j,GridVel_j + nDim, Normal, 0.0);
      Mean_ProjVel -= 0.5*(ProjVel_i + ProjVel_j);
    }

    Mean_Density = 0.5*(node[iPoint]->GetDensity() + node[jPoint]->GetDensity());
    Mean_ThermalCond = 0.5*(node[iPoint]->GetThermalConductivity() + node[jPoint]->GetThermalConductivity());
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Ys_i[iSpecies] = node[iPoint]->GetMassFraction(iSpecies);
      Ys_j[iSpecies] = node[jPoint]->GetMassFraction(iSpecies);
    }
    /*
    su2double dim_temp_i = node[iPoint]->GetTemperature()*config->GetTemperature_Ref();
    su2double dim_temp_j = node[jPoint]->GetTemperature()*config->GetTemperature_Ref();
    bool US_System = (config->GetSystemMeasurements() == US);
    if(US_System) {
      dim_temp_i *= 5.0/9.0;
      dit_temp_j *= 5.0/9.0;
    }
    Mean_CV = 0.5*(library->ComputeCV(dim_temp_i,Ys_i) +  library->ComputeCV(dim_temp_j,Ys_j));
    Mean_CV /= config->GetEnergy_Ref()*config->GetTemperature_Ref():
    if(US_System)
      Mean_CV *= 3.28084*3.28084*9.0/5.0;
    */

    /*--- Inviscid contribution ---*/
    Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
    if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);

    /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
  	Lambda_1 = FOUR_OVER_THREE*Mean_LaminarVisc;
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
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal);
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed();

        /*--- Adjustment for grid movement ---*/
        if(grid_movement) {
          GridVel = geometry->node[iPoint]->GetGridVel();
          su2double ProjVel = std::inner_product(GridVel,GridVel + nDim, Normal, 0.0);
          Mean_ProjVel -= ProjVel;
        }

        Mean_Density = node[iPoint]->GetDensity();
        Mean_ThermalCond = node[iPoint]->GetThermalConductivity();
        Mean_LaminarVisc = node[iPoint]->GetLaminarViscosity();
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          Ys_i[iSpecies] = node[iPoint]->GetMassFraction(iSpecies);
        /*
        dim_temp_i = node[iPoint]->GetTemperature()*config->GetTemperature_Ref();
        bool US_System = (config->GetSystemMeasurements() == US);
        if(US_System)
          dim_temp_i *= 5.0/9.0;
        Mean_CV = library->ComputeCV(dim_temp_i,Ys_i);
        Mean_CV /= config->GetEnergy_Ref()*config->GetTemperature_Ref():
        if(US_System)
          Mean_CV *= 3.28084*3.28084*9.0/5.0;
        */

        /*--- Inviscid contribution ---*/
        Lambda = (std::abs(Mean_ProjVel) + Mean_SoundSpeed)*Area;
        if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Inv(Lambda);

        /*--- Viscous contribution ---*/
      	Lambda_1 = FOUR_OVER_THREE*Mean_LaminarVisc;
      	Lambda_2 = Mean_ThermalCond/Mean_CV;
      	Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      	if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Visc(Lambda);
      	if(geometry->node[jPoint]->GetDomain())
          node[jPoint]->AddMax_Lambda_Visc(Lambda);
      }
    }
  } /*--- End loop boundary edges ---*/

  /*--- Each element uses their own speed for a steady state simulation ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();

    if(Volume > EPS) {
      Local_Delta_Time = config->GetCFL(iMesh)*Volume / node[iPoint]->GetMax_Lambda_Inv();
      Local_Delta_Time_Visc = config->GetCFL(iMesh)*K_v*Volume*Volume/ node[iPoint]->GetMax_Lambda_Visc();
      Local_Delta_Time = std::min(Local_Delta_Time,Local_Delta_Time_Visc);
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
      config->SetCFL(iMesh,config->GetUnst_CFL());

      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
            it computes the time step based on the unsteady CFL ---*/
      if(config->GetCFL(iMesh) < EPS)
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      else
        node[iPoint]->SetDelta_Time(Global_Delta_Time);
    }
  } /*--- End of tiem_stepping check ---*/

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
/*!
 * \brief Set the fluid solver nondimensionalization.
 */
//
//
void CReactiveNSSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {
  /*--- Adimensionalitazion for the invisicid part ---*/
  CReactiveEulerSolver::SetNondimensionalization(geometry,config,iMesh);

  /*--- Get the rank of an eventual parallel simulation ---*/
  int rank;
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
  bool US_Measurament = (config->GetSystemMeasurements() == US);

  if(US_Measurament)
    dim_temp *= 5.0/9.0;
  //Viscosity_FreeStream = library->ComputeViscosity(dim_temp,config->GetMassFrac_FreeStream());
  if(US_Measurament)
    Viscosity_FreeStream *= 0.02088553108;
  config->SetViscosity_FreeStream(Viscosity_FreeStream);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;
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

    config->SetReynolds(config->GetDensity_FreeStream()*config->GetModVel_FreeStream()*config->GetLength_Ref()/Viscosity_FreeStream);
    std::cout<< "Reynolds number (non-dim): " << config->GetReynolds() <<std::endl;

    std::cout<< "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << std::endl;

    std::cout<<std::endl;
  }
}

//
//
/*!
 *\brief Viscous residual
 */
//
//
void CReactiveNSSolver::Viscous_Residual(CGeometry* geometry, CSolver** solution_container, CNumerics* numerics,
                                         CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iPoint, jPoint, iEdge;
  unsigned short iDim, jDim, iSpecies;

  /*--- Initialize the viscous residual to zero ---*/
  SU2_Assert(Res_Visc != NULL,"The array for viscous residual has not been allocated");
  std::fill(Res_Visc, Res_Visc + nVar, 0.0);

  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Set Primitive variables and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());

    /*--- Set primitve variables limited ---*/
    numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());

    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[jPoint]->GetThermalConductivity());

    /*--- Species binary coefficients ---*/
    auto avggrad_numerics = dynamic_cast<CAvgGradReactive_Flow*>(numerics);
    SU2_Assert(avggrad_numerics != NULL, "The cast to compute the viscous flux has not been successfull");
    avggrad_numerics->SetBinaryDiffCoeff(dynamic_cast<CReactiveNSVariable*>(node[iPoint])->GetBinaryDiffusionCoeff(),
                                         dynamic_cast<CReactiveNSVariable*>(node[jPoint])->GetBinaryDiffusionCoeff());

    /*--- Compute the residual ---*/
    try {
      numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      avggrad_numerics->SetExplicit();
      numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = false;
    err = !std::none_of(Res_Visc,Res_Visc + nVar,[](su2double elem){return std::isnan(elem);});

    /*--- Update the residual and Jacobian ---*/
    if(!err) {
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      LinSysRes.AddBlock(jPoint, Res_Visc);
      if(implicit)
        throw Common::NotImplemented("Implicit computation for viscous residual not implemented");
    }
  } /*--- End for loop over interior edges ---*/
}

//
//
/*!
 *\brief Set gradient primitive variables according to least squares
 */
//
//
void CReactiveNSSolver::SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) {
  /*--- Local variables ---*/
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, jDim, iSpecies, iVar, iMarker;

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
    PrimVar_i[CReactiveNSVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::T_INDEX_PRIM);
    PrimVar_j[CReactiveNSVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveNSVariable::T_INDEX_PRIM);
    PrimVar_i[CReactiveNSVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::P_INDEX_PRIM);
    PrimVar_j[CReactiveNSVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveNSVariable::P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim) {
      PrimVar_i[CReactiveNSVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveNSVariable::VX_INDEX_PRIM + iDim);
      PrimVar_j[CReactiveNSVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveNSVariable::VX_INDEX_PRIM + iDim);
    }
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      PrimVar_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies] =
      node[iPoint]->GetPrimitive(CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies);

      PrimVar_j[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies] =
      node[iPoint]->GetPrimitive(CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies);
    }

    /*--- Compute molar fractions for iPoint and jPoint ---*/
    //Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD,
    //                                          PrimVar_i.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD + nSpecies));
    //Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD,
    //                                          PrimVar_j.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD + nSpecies));
    std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.begin() + CReactiveNSVariable::RHOS_INDEX_GRAD);
    std::copy(Xs_j.cbegin(), Xs_j.cend(), PrimVar_j.begin() + CReactiveNSVariable::RHOS_INDEX_GRAD);

  	Normal = geometry->edge[iEdge]->GetNormal();

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
        /*--- Get primitives from CVariable ---*/
        PrimVar_Vertex[CReactiveNSVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::T_INDEX_PRIM);
        PrimVar_Vertex[CReactiveNSVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::P_INDEX_PRIM);
        for(iDim = 0; iDim < nDim; ++iDim)
          PrimVar_Vertex[CReactiveNSVariable::VX_INDEX_GRAD + iDim] =
          node[iPoint]->GetPrimitive(CReactiveNSVariable::VX_INDEX_PRIM + iDim);
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          PrimVar_Vertex[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies] =
          node[iPoint]->GetPrimitive(CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies);

        /*--- Compute molar fractions ---*/
        //Xs = library->GetMolarFromMass(RealVec(PrimVar_Vertex.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD,
        //                                       PrimVar_Vertex.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD + nSpecies));
        std::copy(Xs.cbegin(), Xs.cend(), PrimVar_Vertex.begin() + CReactiveNSVariable::RHOS_INDEX_GRAD);

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
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
    SU2_Assert(Volume > EPS,"The measure of the volume is not consistent(smaller or equalt to zero)");
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  	  for(iDim = 0; iDim < nDim; ++iDim)
  		  node[iPoint]->SetGradient_Primitive(iVar, iDim, node[iPoint]->GetGradient_Primitive(iVar,iDim)/Volume);
    }
  }

  Set_MPI_Primitive_Gradient(geometry, config);
}

//
//
/*!
 *\brief Set gradient primitive variables according to least squares
 */
//
//
void CReactiveNSSolver::SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) {
  unsigned short iVar, iDim, jDim, iSpecies, iNeigh;
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
    PrimVar_i[CReactiveNSVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::T_INDEX_PRIM);
    PrimVar_i[CReactiveNSVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveNSVariable::P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim)
      PrimVar_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveNSVariable::VX_INDEX_PRIM + iDim);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      PrimVar_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies] =
      node[iPoint]->GetPrimitive(CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies);

    /*--- Compute molar fractions for iPoint ---*/
    //Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD,
    //                                         PrimVar_i.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD + nSpecies));
    std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.begin() + CReactiveNSVariable::RHOS_INDEX_GRAD);

    /*--- Inizialization of variables ---*/
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      std::fill(Cvector[iVar], Cvector[iVar] + nDim, 0.0);

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		//AD::StartPreacc();
    //AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    //AD::SetPreaccIn(Coord_i, nDim);

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      /*--- Get coordinates ---*/
			auto Coord_j = geometry->node[jPoint]->GetCoord();

      /*--- Get primitives from CVariable ---*/
      PrimVar_j[CReactiveNSVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveNSVariable::T_INDEX_PRIM);
      PrimVar_j[CReactiveNSVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveNSVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        PrimVar_j[CReactiveNSVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveNSVariable::VX_INDEX_PRIM + iDim);
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        PrimVar_j[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies] =
        node[jPoint]->GetPrimitive(CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies);

      /*--- Compute molar fractions for jPoint ---*/
      //Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD,
      //                                         PrimVar_j.cbegin() + CReactiveNSVariable::RHOS_INDEX_GRAD + nSpecies));
      std::copy(Xs_j.cbegin(), Xs_i.cend(), PrimVar_j.begin() + CReactiveNSVariable::RHOS_INDEX_GRAD);

      //AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

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
      if(r33-r23*r23-r13*r13 > EPS)
        r33 = std::sqrt(r33-r23*r23-r13*r13);
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
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }

	  /*--- Computation of the gradient: C*S ---*/
    su2double result;
    for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      result = 0.0;
		  for(iDim = 0; iDim < nDim; ++iDim) {
        for(jDim = 0; jDim < nDim; ++jDim)
          result += Cvector[iVar][jDim]*Smatrix[iDim][jDim];
        node[iPoint]->SetGradient_Primitive(iVar,iDim,result);
      }
    }

    //	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    //  AD::EndPreacc();
	} /*--- End of iPoint for loop ---*/

  Set_MPI_Primitive_Gradient(geometry,config);
}

//
//
/*!
 *\brief Isothermal wall boundary condition
 */
//
//
void CReactiveNSSolver::BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                           CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) {
  SU2_Assert(Res_Conv != NULL,"The array of convective residual for boundary conditions has not been allocated");
  SU2_Assert(Res_Visc != NULL,"The array of viscous residual for boundary conditions has not been allocated");
  SU2_Assert(Vector != NULL,"The array to store velocity for boundary conditions has not been allocated");

  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr;
  su2double Ti, Tj, dTdn, Twall;
  su2double dij, Area, C = 5.0;
  su2double Normal[nDim], UnitNormal[nDim], GridVel[nDim];

	/*--- Identify the boundary ---*/
	auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

	/*--- Retrieve the specified wall temperature ---*/
	Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

	/*--- Loop over boundary points to calculate energy flux ---*/
	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; ++iVertex) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
		if(geometry->node[iPoint]->GetDomain()) {

			/*--- Compute dual-grid area and boundary normal ---*/
			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
			Area = ::ComputeArea(Normal,nDim);
      for(iDim = 0; iDim < nDim; ++iDim)
        UnitNormal[iDim] = -Normal[iDim]/Area;

			/*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      auto Coord_i = geometry->node[iPoint]->GetCoord();
      auto Coord_j = geometry->node[jPoint]->GetCoord();
      dij = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = std::sqrt(dij);

      /*--- Store the corrected velocity at the wall which will be zero (u = 0)
            unless there is grid motion (u = u_wall) ---*/
      if(grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        std::copy(GridVel, GridVel + nDim, Vector);
      }
      else
        std::fill(Vector, Vector + nDim, 0.0);

      /*--- Initialize residual to zero ---*/
      std::fill(Res_Visc, Res_Visc + nVar, 0.0);

      /*--- Set the residual, truncation error and velocity value on the boundary ---*/
			node[iPoint]->SetVelocity_Old(Vector);
			for(iDim = 0; iDim < nDim; ++iDim) {
        LinSysRes.SetBlock_Zero(iPoint, CReactiveNSVariable::RHOVX_INDEX_SOL + iDim);
        node[iPoint]->SetVal_ResTruncError_Zero(CReactiveNSVariable::RHOVX_INDEX_SOL + iDim);
      }

      /*--- Extract the interior temperature and the thermal conductivity for the boundary node ---*/
      Tj   = node[jPoint]->GetTemperature();
      ktr  = node[iPoint]->GetThermalConductivity();

      /*--- Compute normal gradient ---*/
      dTdn = (Twall - Tj)/dij;

      /*--- Apply to the linear system ---*/
      Res_Visc[CReactiveNSVariable::RHOE_INDEX_SOL] = ktr*dTdn*Area;
      //1)dTdn = (Ti - Tj)/dij;
      //Res_Visc[CReactiveNSVariable::RHOE_INDEX_SOL] = (ktr*dTdn + C*(Twall-Ti)/(2^config->GetExtIter()*dij))*Area;
      //2)dTdn = (Ti - Tj)/dij;
      //Ti = node[iPoint]->GetTemperature();
      //Res_Visc[CReactiveNSVariable::RHOE_INDEX_SOL] = (ktr*dTdn + ktr*C*(Twall-Ti)/(2^config->GetExtIter()*dij))*Area;
      //3)dTdn = (Twall - Tj)/dij;
      //Res_Visc[CReactiveNSVariable::RHOE_INDEX_SOL] = ktr*dTdn*Area;
      //bool nonPhys_temp = node[iPoint]->SetTemperature(Twall);
      //SU2_Assert(nonPhys_temp == false,"The wall temperature to impose is not a valid temperature");

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if(implicit)
        throw Common::NotImplemented("Implicit computations for isothermal wall not implemented");

      /*--- If the wall is moving, there are additional residual contributions
            due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/
      if(grid_movement) {
        /*--- Get the grid velocity at the current boundary node ---*/
        GridVel = geometry->node[iPoint]->GetGridVel();
        su2double ProjGridVel = std::inner_product(GridVel, GridVel + nDim, UnitNormal, 0.0);

        /*--- Retrieve other primitive quantities and viscosities ---*/
        su2double Density  = node[iPoint]->GetDensity();
        su2double Pressure = node[iPoint]->GetPressure();
        su2double laminar_viscosity = node[iPoint]->GetLaminarViscosity();

        unsigned short jDim;
        su2double Grad_Vel[nDim][nDim] = {};
        for(iDim = 0; iDim < nDim; ++iDim)
          for(jDim = 0; jDim < nDim; ++jDim)
            Grad_Vel[iDim][jDim] = node[iPoint]->GetGradient_Primitive(CReactiveNSVariable::VX_INDEX_GRAD + iDim, jDim);

        /*--- Divergence of the velocity ---*/
        su2double div_vel = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          div_vel += Grad_Vel[iDim][iDim];

        /*--- Compute the viscous stress tensor ---*/
        su2double tau[nDim][nDim] = {};
        su2double delta[nDim][nDim] = {};

        for(iDim = 0; iDim < nDim; iDim++) {
          /*--- Fill delta's row of zeros for safety ---*/
          std::fill(delta[iDim], delta[iDim] + nDim, 0.0);
          /*--- Compute dleta and tau---*/
          delta[iDim][iDim] = 1.0;
          for(jDim = 0; jDim < nDim; jDim++)
            tau[iDim][jDim] = laminar_viscosity*(Grad_Vel[iDim][jDim] + Grad_Vel[jDim][iDim]) -
                                                 TWO3*laminar_viscosity*div_vel*delta[iDim][jDim];
        }

        /*--- Dot product of the stress tensor with the grid velocity ---*/
        su2double tau_vel[nDim] = {};
        std::fill(tau_vel, tau_vel + nDim, 0.0);
        for(iDim = 0; iDim < nDim; ++iDim) {
          for(jDim = 0; jDim < nDim; ++jDim)
            tau_vel[iDim] += tau[iDim][jDim]*GridVel[jDim];
        }

        /*--- Compute the convective and viscous residuals (energy eqn.) ---*/
        std::fill(Res_Conv, Res_Conv + nVar, 0.0);
        Res_Conv[CReactiveNSVariable::RHOE_INDEX_SOL] = Pressure*ProjGridVel*Area;
        for(iDim = 0; iDim < nDim; iDim++)
          Res_Visc[CReactiveNSVariable::RHOE_INDEX_SOL] += tau_vel[iDim]*UnitNormal[iDim]*Area;

        /*--- Implicit Jacobian contributions due to moving walls ---*/
        if(implicit)
          throw Common::NotImplemented("Implicit computations for isothermal wall not implemented");
      }
    }
  } /*--- End of iVertex for loop ---*/
}
