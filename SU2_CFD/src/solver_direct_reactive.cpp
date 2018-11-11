#include "../include/solver_reactive.hpp"
#include "../include/numerics_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/move_pointer.hpp"
#include "../../Common/include/default_initialization.hpp"

#include <algorithm>
#include <limits>
#include <cmath>
#include <iterator>

namespace {
  using SmartArr = CReactiveEulerSolver::SmartArr;
  /*!
   * \brief Compute area forthe current normal
   */
  su2double ComputeArea(const SmartArr& Normal,const unsigned short nDim) {

    su2double Area = std::inner_product(Normal.get(),Normal.get() + nDim, Normal.get(), 0.0);
    Area = std::sqrt(Area);
    return Area;

  }
} /*-- End of unnamed namespace ---*/

//
//
/*!
  *\brief Default constructor
  */
//
//
CReactiveEulerSolver::CReactiveEulerSolver():CSolver(),nSpecies(),space_centered(),implicit(),least_squares(),second_order(),limiter(),
                                             Density_Inf(),Pressure_Inf(),Temperature_Inf() {
  std::tie(nSecondaryVar,nSecondaryVarGrad,nVarGrad,IterLinSolver) = Common::repeat<4,decltype(nVarGrad)>(decltype(nVarGrad)());

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  std::tie(nPoint,nPointDomain,nVar,nPrimVar,nPrimVarGrad) = Common::repeat<5,decltype(nVar)>(decltype(nVar)());
}

//
//
/*!
  *\brief Class constructor
  */
//
//
CReactiveEulerSolver::CReactiveEulerSolver(std::shared_ptr<CGeometry> geometry, std::shared_ptr<CConfig> config,unsigned short iMesh):
                      CSolver(),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {
  unsigned long iPoint;
  unsigned short iVar,iDim;

  std::tie(nSecondaryVar,nSecondaryVarGrad,nVarGrad,IterLinSolver) = Common::repeat<4,decltype(nVarGrad)>(decltype(nVarGrad)());

  bool restart = config->GetRestart() || config->GetRestart_Flow();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  std::string file_name = config->GetSolution_FlowFileName();

  if(!(!restart || iMesh != MESH_0)) {

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if (time_stepping) {
      auto Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter());
      file_name = config->GetUnsteady_FileName(file_name, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/
    //Read_SU2_Restart_Metadata(geometry.get(), config.get(), false, file_name);
  }

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho,rho*vx,rho*vy,rho*vz,rho*E,rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T,vx,vy,vz,P,rho,h,a,rhoCv,rho1,...rhoNs)^T ---*/
  nPrimVarGrad = nDim + 2; /*--- Gradient Primitive variables (T,vx,vy,vz,P,rho)^T ---*/

  /*--- Perform the non-dimensionalization forthe flow equations using the specified reference values. ---*/
  SetNondimensionalization(geometry.get(), config.get(), iMesh);

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

	/*--- Define some structure for locating max residuals ---*/
  /*
  Point_Max = new unsigned long[nVar];

  Point_Max_Coord = new su2double*[nVar];
	for(iVar = 0; iVar < nVar; ++iVar)
		Point_Max_Coord[iVar] = new su2double[nDim];
  */

  /*--- Allocate vectors related to the solution ---*/
  Sol_i.resize(nVar);
  Sol_j.resize(nVar);
  Primitive_i.resize(nPrimVar);
  Primitive_j.resize(nPrimVar);

  /*--- Allocate arrays forconserved variable limits ---*/
  Lower_Limit.resize(nPrimVar);
  Upper_Limit.resize(nPrimVar);

  std::fill(Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL,Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL+nDim, -1E16);
  std::fill(Upper_Limit.begin(),Upper_Limit.end(),1E16);

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
		for(iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new su2double [nVar];
			Jacobian_j[iVar] = new su2double [nVar];
		}
  }
  else
    implicit = false;

  /*--- Allocate arrays for gradient computation by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      least_squares = true;
		  /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		  Smatrix = new su2double* [nDim];
      for(iDim = 0; iDim < nDim; ++iDim)
        Smatrix[iDim] = new su2double[nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
		  Cvector = new su2double* [nPrimVarGrad];
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
        Cvector[iVar] = new su2double[nDim];
	}
  else
    least_squares = false;

  space_centered = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  second_order = config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER;
  unsigned long ExtIter = config->GetExtIter();
  limiter = config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter();

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(),config->GetVelocity_FreeStreamND() + nDim);
  MassFrac_Inf     = RealVec(config->GetMassFrac_FreeStream(),config->GetMassFrac_FreeStream() + nSpecies);

  node_infty = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

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
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	//Set_MPI_Solution(geometry.get(), config.get());

}
//
//

//
//
/*!
 * \brief Check if there is any non physical point in the initialization
 */
//
//
void CReactiveEulerSolver::Check_FreeStream_Solution(std::shared_ptr<CConfig> config) {

  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  RealVec Solution(nVar);

  bool check_infty = node_infty->SetPrimVar(config.get());

  unsigned long counter_local = 0, counter_global;
  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    bool nonPhys = node[iPoint]->SetPrimVar(config.get());

    if(nonPhys) {

      su2double rho,Cp_i, hTot;
      su2double sqvel = std::inner_product(Velocity_Inf.cbegin(),Velocity_Inf.cend(),Velocity_Inf.cbegin(),0.0);

      /*--- Compute density from supplied quantities ---*/
      //rho = library->ComputeDensity(MassFrac_Inf,Pressure_Inf,Temperature_Inf);
      Solution[CReactiveEulerVariable::RHO_INDEX_SOL] = rho;
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Solution[CReactiveEulerVariable::RHOS_INDEX_SOL+iSpecies] = MassFrac_Inf[iSpecies]*rho;

      /*--- Calculate momentum energy from supplied primitive quanitites ---*/
      for(unsigned short iDim = 0; iDim < nDim; ++iDim)
        Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = rho*Velocity_Inf[iDim];
      Solution[CReactiveEulerVariable::RHOE_INDEX_SOL] = hTot + 0.5*sqvel;

      node[iPoint]->SetSolution(Solution.data());
      node[iPoint]->SetSolution_Old(Solution.data());

      counter_local++;
    }
  }

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

//
//
/*!
 * \brief Set the fluid solver nondimensionalization.
 */
 void CReactiveEulerSolver::SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) {

   su2double Temperature_FreeStream = 0.0, ModVel_FreeStream = 0.0, Energy_FreeStream = 0.0,
             Density_FreeStream = 0.0, Pressure_FreeStream = 0.0,Mach2Vel_FreeStream = 0.0,
             Gamma = 0.0, Gamma_Minus_One = 0.0;

   su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
             Temperature_Ref = 0.0, Time_Ref = 0.0, Gas_Constant_Ref = 0.0, Energy_Ref = 0.0;

   su2double Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,ModVel_FreeStreamND = 0.0,
             Temperature_FreeStreamND = 0.0, Gas_Constant = 0.0, Gas_ConstantND = 0.0,
             Energy_FreeStreamND = 0.0, Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;

   RealVec  Velocity_FreeStreamND(nDim);

   unsigned short iDim;

   /*--- Local variables ---*/
   bool unsteady = config->GetUnsteady_Simulation();

   int rank = MASTER_NODE;
   #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   #endif

 	 /*--- Compressible non dimensionalization ---*/
   Pressure_FreeStream = config->GetPressure_FreeStream();
   Temperature_FreeStream = config->GetTemperature_FreeStream();
   Gas_Constant = library->GetRgas();
   Density_FreeStream = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
   config->SetDensity_FreeStream(Density_FreeStream);
   //library->Gamma_FrozenSoundSpeed(Temperature_FreeStream,Pressure_FreeStream,Density_FreeStream,Gamma,Mach2Vel_FreeStream);
   Gamma_Minus_One = Gamma - 1.0;

   /*--- Compute the Free Stream velocity, using the Mach number ---*/
   /*
   if(nDim == 2) {
     config->GetVelocity_FreeStream()[0] = std::cos(Alpha)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[1] = std::sin(Alpha)*Mach*Mach2Vel_FreeStream;
   }
   if(nDim == 3) {
     config->GetVelocity_FreeStream()[0] = std::cos(Alpha)*std::cos(Beta)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[1] = std::sin(Beta)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[2] = std::sin(Alpha)*std::cos(Beta)*Mach*Mach2Vel_FreeStream;
   }
   */

   /*--- Compute the modulus of the free stream velocity ---*/
   ModVel_FreeStream = std::inner_product(config->GetVelocity_FreeStream(),config->GetVelocity_FreeStream() + nDim,
                                          config->GetVelocity_FreeStream(),0.0);
   ModVel_FreeStream = std::sqrt(ModVel_FreeStream);
   config->SetModVel_FreeStream(ModVel_FreeStream);

   Energy_FreeStream = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+ 0.5*ModVel_FreeStream*ModVel_FreeStream;

   /*--- Compute non dimensional quantities: Notice that the grid is in meters. ---*/
   if(config->GetRef_NonDim() == DIMENSIONAL) {
     Pressure_Ref      = 1.0;
     Density_Ref       = 1.0;
     Temperature_Ref   = 1.0;
     Length_Ref        = 1.0;
   }

   else {
     Pressure_Ref     = config->GetPressure_Ref();
     Density_Ref      = config->GetDensity_Ref();
     Temperature_Ref  = config->GetTemperature_Ref();
     Length_Ref       = config->GetLength_Ref();
   }
   /*
   else if(config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
     Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
     Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
     Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
   }
   else if(config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
     Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
     Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
     Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
   }
   else if(config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
     Pressure_Ref      = Mach*Mach*Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/(Gamma*(M_inf)^2)
     Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
     Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
   }
   */
   config->SetPressure_Ref(Pressure_Ref);
   config->SetDensity_Ref(Density_Ref);
   config->SetTemperature_Ref(Temperature_Ref);
   config->SetLength_Ref(Length_Ref);

   Velocity_Ref = std::sqrt(Pressure_Ref/Density_Ref);
   config->SetVelocity_Ref(Velocity_Ref);

   Time_Ref = Length_Ref/Velocity_Ref;
   config->SetTime_Ref(Time_Ref);

   Energy_Ref = Pressure_Ref/(Density_Ref*Gamma_Minus_One)+ 0.5*Velocity_Ref*Velocity_Ref;
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

   Gas_ConstantND = Gas_Constant/Gas_Constant_Ref;
   config->SetGas_ConstantND(Gas_ConstantND);

   ModVel_FreeStreamND = std::inner_product(Velocity_FreeStreamND.cbegin(),Velocity_FreeStreamND.cend(),Velocity_FreeStreamND.cbegin(),0.0);
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

     bool SI_Measurement = config->GetSystemMeasurements() == SI;
     bool US_Measuremanet = config->GetSystemMeasurements() == US;
     if(!SI_Measurement && !US_Measuremanet)
      throw std::out_of_range("Unknown option in the unit measure for adimensionalitazion");

     std::cout << "Reference specific gas constant: " << Gas_Constant_Ref;
     if(SI_Measurement)
       std::cout << " N.m/kg.K." << std::endl;
     else if(US_Measuremanet)
       std::cout << " lbf.ft/slug.R." << std::endl;

     std::cout << "Reference pressure: " << Pressure_Ref;
     if(SI_Measurement)
       std::cout << " Pa." << std::endl;
     else if(US_Measuremanet)
       std::cout << " psf." << std::endl;

     std::cout << "Reference temperature: " << Temperature_Ref;
     if(SI_Measurement)
       std::cout << " K." << std::endl;
     else if(US_Measuremanet)
       std::cout << " R." << std::endl;

     std::cout << "Reference density: " << Density_Ref;
     if(SI_Measurement)
       std::cout << " kg/m^3." << std::endl;
     else if(US_Measuremanet)
       std::cout << " slug/ft^3." << std::endl;

     std::cout << "Reference velocity: " << Velocity_Ref;
     if(SI_Measurement)
       std::cout << " m/s." << std::endl;
     else if(US_Measuremanet)
       std::cout << " ft/s." << std::endl;

     std::cout << "Reference energy per unit mass: " << Energy_Ref;
     if(SI_Measurement)
       std::cout << " m^2/s^2." << std::endl;
     else if(US_Measuremanet)
       std::cout << " ft^2/s^2." << std::endl;

     if(unsteady)
       std::cout << "Reference time: " << Time_Ref <<" s." << std::endl;

     /*--- Print out resulting non-dim values here. ---*/
     std::cout << "-- Resulting non-dimensional state:" << std::endl;
     std::cout << "Mach number (non-dim): " << config->GetMach() << std::endl;

     std::cout << "Specific gas constant (non-dim): " << Gas_ConstantND << std::endl;

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
  *\brief Variables preprocessing
  */
//
//
void CReactiveEulerSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                         unsigned short iMesh, unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem, bool Output) {

  bool interface = config->GetnMarker_InterfaceBound() != 0;

  /*--- Compute Interface MPI ---*/
  if(interface)
    Set_MPI_Interface(geometry, config);

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
/*!
 *\brief Set gradient primitive variables according to Green-Gauss
 */
//
//
void CReactiveEulerSolver::SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) {
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, iVar, iMarker;
  su2double PrimVar_Average, Partial_Gradient, Partial_Res;
  RealVec PrimVar_Vertex(nPrimVarGrad), /*--- Gradient of the following primitive variables: [T,u,v,w,p]^T ---*/
          PrimVar_i(nPrimVarGrad),
          PrimVar_j(nPrimVarGrad);

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
      PrimVar_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
    }

  	auto Normal = Common::wrap_in_unique(geometry-> edge[iEdge]->GetNormal());
  	for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  		PrimVar_Average =  0.5 * (PrimVar_i[iVar] + PrimVar_j[iVar]);
  		for(iDim = 0; iDim < nDim; ++iDim) {
  			Partial_Res = PrimVar_Average*Normal[iDim];
  			if(geometry->node[iPoint]->GetDomain())
  				node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
  			if(geometry->node[jPoint]->GetDomain())
  				node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  		}
  	}
  }

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
  	for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
  		iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  		if(geometry->node[iPoint]->GetDomain()) {
        /*--- Get primitives from CVariable ---*/
        PrimVar_Vertex[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
        PrimVar_Vertex[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
        for(iDim = 0; iDim < nDim; ++iDim)
          PrimVar_Vertex[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

        auto Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
  			for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
  				for(iDim = 0; iDim < nDim; ++iDim) {
  					Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
  					node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  				}
  		  }
  		}
  	}

  	/*--- Update gradient value ---*/
  	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
  		for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  			for(iDim = 0; iDim < nDim; ++iDim) {
  				Partial_Gradient = node[iPoint]->GetGradient_Primitive(iVar,iDim) / (geometry->node[iPoint]->GetVolume());
  				node[iPoint]->SetGradient_Primitive(iVar, iDim, Partial_Gradient);
  			}
  		}
  	}

  	//Set_MPI_Primitive_Gradient(geometry, config);

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
	su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, product, detR2, z11, z12, z13, z22, z23, z33;
  su2double weight;
  bool singular;
  RealVec PrimVar_i(nPrimVarGrad), /*--- Gradient of the following primitive variables: [T,u,v,w,p]^T ---*/
          PrimVar_j(nPrimVarGrad);

	/*--- Loop over points of the grid ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		/*--- Set the value of singular ---*/
		singular = false;

    /*--- Get coordinates ---*/
		auto Coord_i = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());

    /*--- Get primitives from CVariable ---*/
    PrimVar_i[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
    PrimVar_i[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
    for(iDim = 0; iDim < nDim; ++iDim)
      PrimVar_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

		/*--- Inizialization of variables ---*/
		for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
      for(iDim = 0; iDim < nDim; ++iDim)
        Cvector[iVar][iDim] = 0.0;

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		//AD::StartPreacc();
    //AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    //AD::SetPreaccIn(Coord_i, nDim);

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			auto Coord_j = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      /*--- Get primitives from CVariable ---*/
      PrimVar_j[CReactiveEulerVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      PrimVar_j[CReactiveEulerVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        PrimVar_j[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

			//AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);
      RealVec Coord_ij(nDim);
      std::transform(Coord_j.get(),Coord_j.get() + nDim,Coord_i.get(),Coord_ij.begin(),std::minus<su2double>());

      weight = std::inner_product(Coord_ij.cbegin(),Coord_ij.cend(),Coord_ij.cbegin(),0.0);

			/*--- Sumations for entries of upper triangular matrix R ---*/
      if(std::abs(weight) > EPS){
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
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
      }
    } /*--- End of iNeigh for loop---*/

		/*--- Entries of upper triangular matrix R ---*/
    if(r11 >= 0.0)
      r11 = std::sqrt(r11);
    else
      r11 = 0.0;

    if(std::abs(r11) > EPS)
      r12 = r12/r11;
    else
      r12 = 0.0;

    if(r22-r12*r12 >= 0.0)
      r22 = std::sqrt(r22-r12*r12);
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
      if(r33-r23*r23-r13*r13 >= 0.0)
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
        for(jDim = 0; jDim < nDim; ++jDim)
          Smatrix[iDim][jDim] = 0.0;
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

		/*--- Computation of the gradient: S*c ---*/
		for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
			for(iDim = 0; iDim < nDim; ++iDim) {
				product = 0.0;
				for(jDim = 0; jDim < nDim; ++jDim)
					product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
				node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
			}
		}

	//	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
  //  AD::EndPreacc();
} /*--- End iPoint for loop ---*/

	//Set_MPI_Primitive_Gradient(geometry, config);
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
  su2double  dave, LimK, eps1, eps2, dm, dp, du, y, limiter;
  RealVec Prim_i(nPrimVarGrad),
          Prim_j(nPrimVarGrad),
          Primitive(nPrimVarGrad);
  //  RealMatrix Gradient_i, Gradient_j;

  dave = config->GetRefElemLength();
  LimK = config->GetLimiterCoeff();

  if(!limiter) {
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint){
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        node[iPoint]->SetLimiter_Primitive(iVar,1.0);
      }
    }
  }

  else {
    /*--- Initialize solution max, solution min and limiter in entire domain ---*/
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
      }
    }

    /*--- Establish bounts for Spekreijse monotonicity by finding max/min values
          of neighbor variables ---*/
    for(iEdge = 0; iEdge < geometry-> GetnEdge(); ++iEdge) {
      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry-> edge[iEdge]->GetNode(0);
      jPoint = geometry-> edge[iEdge]->GetNode(1);

      /*--- Get primitive variables ---*/
      Prim_i[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Prim_i[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      Prim_j[CReactiveEulerVariable::T_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Prim_j[CReactiveEulerVariable::P_INDEX_GRAD] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim) {
        Prim_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
        Prim_j[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[jPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);
      }

      /*--- Compute the max and min values fornodes i & j ---*/
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar){
        du = (Primitive_j[iVar]-Primitive_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
      }
    }
  }

  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  if(config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      auto Coord_i = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());
      auto Coord_j = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

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

      }

      //AD::EndPreacc();

    }

    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }

  }

  /*--- Venkatakrishnan limiter ---*/
  if(config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) {

    /*--- Allocate memory for the max and min primitive value --*/
    RealVec LocalMinPrimitive(nPrimVarGrad), GlobalMinPrimitive(nPrimVarGrad),
            LocalMaxPrimitive(nPrimVarGrad), GlobalMaxPrimitive(nPrimVarGrad);

    /*--- Compute the max value and min value of the solution ---*/
    std::fill(LocalMinPrimitive.begin(),LocalMinPrimitive.end(),std::numeric_limits<su2double>::max());
    std::fill(LocalMaxPrimitive.begin(),LocalMaxPrimitive.end(),std::numeric_limits<su2double>::min());

    /*--- Loop over the internal domain ---*/
    for(iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      /*--- Get the involved primitive variables ---*/
      Primitive[CReactiveEulerVariable::T_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::T_INDEX_PRIM);
      Primitive[CReactiveEulerVariable::P_INDEX_GRAD] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::P_INDEX_PRIM);
      for(iDim = 0; iDim < nDim; ++iDim)
        Primitive[CReactiveEulerVariable::VX_INDEX_GRAD + iDim] = node[iPoint]->GetPrimitive(CReactiveEulerVariable::VX_INDEX_PRIM + iDim);

      std::transform(Primitive.cbegin(), Primitive.cend(), LocalMinPrimitive.begin(),
                      std::back_insert_iterator<RealVec>(LocalMinPrimitive),std::less<su2double>());
      std::transform(Primitive.cbegin(), Primitive.cend(), LocalMaxPrimitive.begin(),
                      std::back_insert_iterator<RealVec>(LocalMaxPrimitive),std::greater<su2double>());
    }

    #ifdef HAVE_MPI
      SU2_MPI::Allreduce(LocalMinPrimitive.data(), GlobalMinPrimitive.data(), nPrimVarGrad, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(LocalMaxPrimitive.data(), GlobalMaxPrimitive.data(), nPrimVarGrad, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #else
    GlobalMinPrimitive = LocalMinPrimitive;
    GlobalMaxPrimitive = LocalMaxPrimitive;
    #endif

    /*--- Loop over the interior edges ---*/
    for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      auto Coord_i = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());
      auto Coord_j = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(eps2);

      for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
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
    }
    //AD::EndPreacc();
    }
  }

  /*--- Limiter MPI ---*/
  //Set_MPI_Primitive_Limiter(geometry, config);
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

  su2double Area, Volume, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda;
  su2double Local_Delta_Time, Global_Delta_Time = 1E6;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker;

  bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time = config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  Min_Delta_Time = 1.E6;
  Max_Delta_Time = 0.0;

  /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint)
    node[iPoint]->SetMax_Lambda_Inv(0.0);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Point identification, Normal vector and area ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    auto Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());

    /*--- Inviscid contribution ---*/
    Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed*Area;

    if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);

  }

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {

        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        auto Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = node[iPoint]->GetProjVel(Normal.get());
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

        /*--- Inviscid contribution ---*/
        Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
        if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
    }
  }

  /*--- Each element uses their own speed for a steady state simulation ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();

    if(std::abs(Volume) > EPS) {
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
	for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		if (geometry->node[iPoint]->GetnPoint() == 1)
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
      if(std::abs(config->GetCFL(iMesh)) < EPS)
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      else
        node[iPoint]->SetDelta_Time(Global_Delta_Time);

    }
  }

  if(dual_time)
    throw Common::NotImplemented("Dual time strategies are not implemented");
}

//
//
/*!
 *\brief Iteration of implicit Euler method
 */
//
//

void CReactiveEulerSolver::ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol = 0;

  /*--- Set maximum residual to zero ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Read the residual ---*/
    auto Local_Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    /*--- Read the volume ---*/
    su2double Volume = geometry->node[iPoint]->GetVolume();

    if(std::abs(node[iPoint]->GetDelta_Time()) > EPS) {
      su2double Delta = Volume / node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      std::fill(Local_Res_TruncError.get(),Local_Res_TruncError.get() + nVar,0.0);
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
    for (iVar = 0; iVar < nVar; iVar++) {
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
			node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
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

  su2double Res;
  unsigned short iVar;
  unsigned long iPoint;

  for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    su2double Volume = geometry->node[iPoint]->GetVolume();
    su2double Delta = node[iPoint]->GetDelta_Time() / Volume;

    auto Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    auto Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for(iVar = 0; iVar < nVar; ++iVar) {
        Res = Residual[iVar] + Res_TruncError[iVar];
        node[iPoint]->AddSolution(iVar, -Res*Delta);
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

void CReactiveEulerSolver::ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container,
                                        CConfig* config, unsigned short iRKStep) {
  su2double Res;
  unsigned short iVar;
  unsigned long iPoint;

  auto RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

  for(iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    su2double Volume = geometry->node[iPoint]->GetVolume();
    su2double Delta = node[iPoint]->GetDelta_Time() / Volume;

    auto Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    auto Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for(iVar = 0; iVar < nVar; ++iVar) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
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
  unsigned long iEdge, iPoint, jPoint,counter_local=0, counter_global=0;
  unsigned short iDim, iVar;

  su2double Project_Grad_i, Project_Grad_j;

  /*--- Initialize the upwind convective residual to zero ---*/
  SU2_Assert(Res_Sour != NULL,"The array for source residual has not been allocated");
  std::fill(Residual,Residual + nVar,0.0);

  /*--- Loop over all the edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Get primitive variables ---*/
    Primitive_i = RealVec(node[iPoint]->GetPrimitive(),node[iPoint]->GetPrimitive() + nPrimVar);
    Primitive_j = RealVec(node[jPoint]->GetPrimitive(),node[jPoint]->GetPrimitive() + nPrimVar);

    if(second_order) {
      SU2_Assert(Vector_i != NULL, "The array to compute the vector from i to j has not been allocated");
      SU2_Assert(Vector_j != NULL, "The array to compute the vector from j to i has not been allocated");
      for(iDim = 0; iDim < nDim; ++iDim) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = - Vector_i[iDim];
      }
      //auto Gradient_i = static_cast<CReactiveEulerVariable*>(node[iPoint])->GetGradient_Primitive();
      //auto Gradient_j = static_cast<CReactiveEulerVariable*>(node[jPoint])->GetGradient_Primitive();
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[iPoint]->GetGradient_Primitive();

      /*--- Auxiliary vector to store reconstruction ---*/
      RealVec Prim_Recon_i(nPrimVarGrad);
      RealVec Prim_Recon_j(nPrimVarGrad);

      Prim_Recon_i[CReactiveEulerVariable::T_INDEX_GRAD] = Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM];
      Prim_Recon_j[CReactiveEulerVariable::T_INDEX_GRAD] = Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM];
      Prim_Recon_i[CReactiveEulerVariable::P_INDEX_GRAD] = Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM];
      Prim_Recon_j[CReactiveEulerVariable::P_INDEX_GRAD] = Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM];
      std::copy(Primitive_i.cbegin() + CReactiveEulerVariable::VX_INDEX_PRIM, Primitive_i.cbegin() + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                Prim_Recon_i.begin() + CReactiveEulerVariable::VX_INDEX_GRAD);
      std::copy(Primitive_j.cbegin() + CReactiveEulerVariable::VX_INDEX_PRIM, Primitive_j.cbegin() + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                Prim_Recon_j.begin() + CReactiveEulerVariable::VX_INDEX_GRAD);

      for(iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        /*
        Project_Grad_i = std::inner_product(Vector_i,Vector_i + nDim,Gradient_i.cbegin() + iVar*nDim,0.0);
        Project_Grad_j = std::inner_product(Vector_j,Vector_j + nDim,Gradient_j.cbegin() + iVar*nDim,0.0);
        */
        Project_Grad_i = Project_Grad_j = 0.0;
        for(iDim = 0; iDim < nDim; ++iDim) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }

        if(limiter) {
          auto Limiter_i = Common::wrap_in_unique(node[iPoint]->GetLimiter_Primitive());
          auto Limiter_j = Common::wrap_in_unique(node[jPoint]->GetLimiter_Primitive());
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
      if(Prim_Recon_i[CReactiveEulerVariable::T_INDEX_GRAD] > EPS)
        Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM] = Prim_Recon_i[CReactiveEulerVariable::T_INDEX_GRAD];
      else
        non_phys_i = true;
      if(Prim_Recon_i[CReactiveEulerVariable::P_INDEX_GRAD] > EPS)
        Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM] = Prim_Recon_i[CReactiveEulerVariable::P_INDEX_GRAD];
      else
        non_phys_i = true;

      bool non_phys_j = false;
      if(Prim_Recon_j[CReactiveEulerVariable::T_INDEX_GRAD] > EPS)
        Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM] = Prim_Recon_j[CReactiveEulerVariable::T_INDEX_GRAD];
      else
        non_phys_j = false;
      if(Prim_Recon_i[CReactiveEulerVariable::P_INDEX_GRAD] > EPS)
        Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM] = Prim_Recon_j[CReactiveEulerVariable::P_INDEX_GRAD];
      else
        non_phys_j = false;

      /*--- Compute other primitive variables accordingly to the reconstruction ---*/
      if(!non_phys_i) {
        su2double rho_recon_i;
        //rho_recon_i = library->ComputeDensity(Yi,Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM],Primitive_i[CReactiveEulerVariable::P_INDEX_PRIM]);
        std::for_each(Primitive_i.begin() + CReactiveEulerVariable::RHOS_INDEX_PRIM,
                      Primitive_i.begin() + CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies,
                      [=](su2double elem){elem *= rho_recon_i/Primitive_i[CReactiveEulerVariable::RHO_INDEX_PRIM];});
        Primitive_i[CReactiveEulerVariable::RHO_INDEX_PRIM] = rho_recon_i;
        //Primitive_i[CReactiveEulerVariable::H_INDEX_PRIM] = library->ComputeEnthalpy(Primitive_i[CReactiveEulerVariable::T_INDEX_PRIM]);
        //library->Frozen_GammaSoundSpeed(.....................);
      }

      if(!non_phys_j) {
        su2double rho_recon_j;
        //rho_recon_j = library->ComputeDensity(Yj,Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM],Primitive_j[CReactiveEulerVariable::P_INDEX_PRIM]);
        std::for_each(Primitive_j.begin() + CReactiveEulerVariable::RHOS_INDEX_PRIM,
                      Primitive_j.begin() + CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies,
                      [=](su2double elem){elem *= rho_recon_j/Primitive_j[CReactiveEulerVariable::RHO_INDEX_PRIM];});
        Primitive_j[CReactiveEulerVariable::RHO_INDEX_PRIM] = rho_recon_j;
        //Primitive_j[CReactiveEulerVariable::H_INDEX_PRIM] = library->ComputeEnthalpy(Primitive_j[CReactiveEulerVariable::T_INDEX_PRIM]);
        //library->Frozen_GammaSoundSpeed(.....................);
      }

      numerics->SetPrimitive(Primitive_i.data(), Primitive_j.data());
    } /*--- End second order reconstruction ---*/
    else {
      /*--- Set primitive variables without reconstruction ---*/
      numerics->SetPrimitive(Primitive_i.data(), Primitive_j.data());
    }

    /*--- Compute the residual ---*/
    try {
      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      static_cast<CUpwReactiveAUSM*>(numerics)->SetExplicit();
      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err;
    err = !std::none_of(Res_Conv,Res_Conv + nVar,[](su2double elem){return std::isnan(elem);});
    if (implicit)
      throw Common::NotImplemented("Implicit computation for upwind residual not implemented");
      /*
      for(iVar = 0; iVar < nVar; ++iVar)
        for(jVar = 0; jVar < nVar; ++jVar)
          if(std::isnan(Jacobian_i[iVar][jVar]) || std::isnan(Jacobian_j[iVar][jVar]))
            err = true;
            */

    /*--- Update residual value ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);

      /*--- Set implicit Jacobians ---*/
      if(implicit) {
        throw Common::NotImplemented("Implicit computation for upwind residual not implemented");
        /*
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
        */
      }
    }
    else
      throw std::runtime_error("NaN found in the upwind residual");
  }
  /*--- Warning message about non-physical reconstructions ---*/
  if(config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  #else
    counter_global = counter_local;
  #endif
    if (iMesh == MESH_0)
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
  unsigned long iPoint,counter_local=0, counter_global=0;

  /*--- Initialize the source residual to zero ---*/
  SU2_Assert(Res_Sour != NULL,"The array for source residual has not been allocated");
  std::fill(Res_Sour,Res_Sour + nVar,0.0);

  /*--- Loop over all points ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Load the conservative variables ---*/
    numerics->SetConservative(node[iPoint]->GetSolution(),node[iPoint]->GetSolution());
    /*--- Load the volume of the dual mesh cell ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    /*--- Compute the source residual ---*/
    try {
      numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      static_cast<CSourceReactive*>(numerics)->SetExplicit();
      numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err;
    err = !std::none_of(Res_Sour,Res_Sour + nVar,[](su2double elem){return std::isnan(elem);});
    if(implicit)
      throw Common::NotImplemented("Implicit computation for source residual not implemented");
      /*
      for(iVar = 0; iVar < nVar; ++iVar)
        for(jVar = 0; jVar < nVar; ++jVar)
          if(std::isnan(Jacobian_i[iVar][jVar]))
            err = true;
            */

    /*--- Add the source residual to the total ---*/
    if(!err) {
      LinSysRes.AddBlock(iPoint, Residual);
      /*--- Add the implicit Jacobian contribution ---*/
      if(implicit)
        throw Common::NotImplemented("Implicit computation for source residual not implemented");
        //Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
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
    if (iMesh == MESH_0)
      config->SetNonphysical_Reconstr(counter_global);
  }
}

//
//
/*!
  *\brief Class constructor
  */
//
//
CReactiveNSSolver::CReactiveNSSolver(std::shared_ptr<CGeometry> geometry, std::shared_ptr<CConfig> config,unsigned short iMesh):
                   CReactiveEulerSolver() {
  unsigned long iPoint;
  unsigned short iVar,iDim;

  std::tie(nSecondaryVar,nSecondaryVarGrad,nVarGrad,IterLinSolver) = Common::repeat<4,decltype(nVarGrad)>(decltype(nVarGrad)());

  library = LibraryPtr(new Framework::ReactingModelLibrary(config->GetLibraryName()));
  nSpecies = library->GetNSpecies();

  bool restart = config->GetRestart() || config->GetRestart_Flow();
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  std::string file_name = config->GetSolution_FlowFileName();

  if(!(!restart || iMesh != MESH_0)) {

    /*--- Modify file name for a time stepping unsteady restart ---*/
    if (time_stepping) {
      auto Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter());
      file_name = config->GetUnsteady_FileName(file_name, Unst_RestartIter);
    }

    /*--- Read and store the restart metadata. ---*/
    //Read_SU2_Restart_Metadata(geometry.get(), config.get(), false, file_name);
  }

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho,rho*vx,rho*vy,rho*vz,rho*E,rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T,vx,vy,vz,P,rho,h,a,rhoCv,rho1,...rhoNs)^T ---*/
  nPrimVarGrad = nSpecies + 2; /*--- Gradient Primitive variables (T,vx,vy,vz,P,rho)^T ---*/
  nPrimVarAvgGrad = nSpecies + nDim + 1; /*--- Gradient variables for average gradient computation (T,vx,vy,vz,Y1....YNs)^T ---*/

  /*--- Perform the non-dimensionalization forthe flow equations using the specified reference values. ---*/
  SetNondimensionalization(geometry.get(), config.get(), iMesh);

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

	/*--- Define some structure for locating max residuals ---*/
  /*
  Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
	for(iVar = 0; iVar < nVar; ++iVar)
		Point_Max_Coord[iVar] = new su2double[nDim];
	*/

  /*--- Allocate vectors related to the solution ---*/
  Sol_i.resize(nVar);
  Sol_j.resize(nVar);
  Primitive_i.resize(nPrimVar);
  Primitive_j.resize(nPrimVar);

  /*--- Allocate arrays forconserved variable limits ---*/
  Lower_Limit.resize(nPrimVar);
  Upper_Limit.resize(nPrimVar);

  std::fill(Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL,Lower_Limit.begin() + CReactiveEulerVariable::RHOVX_INDEX_SOL+nDim, -1E16);
  std::fill(Upper_Limit.begin(),Upper_Limit.end(),1E16);

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
		for(iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new su2double [nVar];
			Jacobian_j[iVar] = new su2double [nVar];
		}
  }
  else
    implicit = false;

  /*--- Allocate arrays for gradient computation by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      least_squares = true;
		  /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		  Smatrix = new su2double* [nDim];
      for(iDim = 0; iDim < nDim; ++iDim)
        Smatrix[iDim] = new su2double[nDim];
      /*--- c vector := transpose(WA)*(Wb) ---*/
		  Cvector = new su2double* [nPrimVarGrad];
      for(iVar = 0; iVar < nPrimVarGrad; ++iVar)
        Cvector[iVar] = new su2double[nDim];
	}
  else
    least_squares = false;

  space_centered = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  second_order = config->GetSpatialOrder_Flow() == SECOND_ORDER || config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER;
  unsigned long ExtIter = config->GetExtIter();
  limiter = config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER && ExtIter <= config->GetLimiterIter();

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();

  Velocity_Inf     = RealVec(config->GetVelocity_FreeStreamND(),config->GetVelocity_FreeStreamND() + nDim);
  MassFrac_Inf     = RealVec(config->GetMassFrac_FreeStream(),config->GetMassFrac_FreeStream() + nSpecies);

  Viscosity_Inf    = config->GetViscosity_FreeStreamND();

  node_infty = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

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
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CReactiveNSVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf, nDim, nVar, config);

  /*--- Use a function to check that the initial solution is physical ---*/
  Check_FreeStream_Solution(config);

  /*--- MPI solution ---*/
	//Set_MPI_Solution(geometry.get(), config.get());

}
//
//

//
//
/*!
  *\brief Variables preprocessing
  */
//
//
void CReactiveNSSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                      unsigned short iMesh, unsigned short iRKStep,
                                      unsigned short RunTime_EqSystem, bool Output) {

  bool interface = config->GetnMarker_InterfaceBound() != 0;
  bool nearfield = config->GetnMarker_NearFieldBound() != 0;

  /*--- Compute Interface MPI ---*/
  if(interface)
    Set_MPI_Interface(geometry, config);

  /*--- Compute NearField MPI ---*/
  if(nearfield)
    Set_MPI_Nearfield(geometry, config);

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
/*!
 *\brief Setting time step
 */
//
//

void CReactiveNSSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                     unsigned short iMesh, unsigned long Iteration) {

  su2double Area, Volume;
  su2double Local_Delta_Time, Global_Delta_Time = 1E6;
  su2double Local_Delta_Time_Visc;
  su2double Mean_SoundSpeed, Mean_ProjVel;
  su2double Mean_LaminarVisc, Mean_ThermalCond, Mean_Density,Mean_CV;
  su2double Lambda, Lambda_1, Lambda_2, K_v = 0.5;
  unsigned long iEdge, iVertex, iPoint, jPoint;
  unsigned short iMarker;

  bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  bool dual_time = config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

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

    auto Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/
    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());
    Mean_Density = 0.5 * (node[iPoint]->GetDensity() + node[jPoint]->GetDensity());
    Mean_ThermalCond = 0.5 * (node[iPoint]->GetThermalConductivity() + node[jPoint]->GetThermalConductivity());
    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    //Mean_CV = 0.5*(library->ComputeCV(node[iPoint]->GetTemperature()) + library->ComputeCV(node[iPoint]->GetTemperature()));

    /*--- Inviscid contribution ---*/
    Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed*Area;
    if(geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if(geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Visc(Lambda);

    /*--- Determine the viscous spectral radius and apply it to the control volume ---*/
  	Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
  	Lambda_2 = Mean_ThermalCond/Mean_CV;
  	Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
  	if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);
  	if (geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Visc(Lambda);

  }

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {

        /*--- Point identification, Normal vector and area ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        auto Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/
        Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
        Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed());
        Mean_Density = 0.5 * (node[iPoint]->GetDensity() + node[jPoint]->GetDensity());
        Mean_ThermalCond = 0.5 * (node[iPoint]->GetThermalConductivity() + node[jPoint]->GetThermalConductivity());
        Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
        //Mean_CV = 0.5*(library->ComputeCV(node[iPoint]->GetTemperature()) + library->ComputeCV(node[iPoint]->GetTemperature()));

        /*--- Inviscid contribution ---*/
        Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
        if(geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Visc(Lambda);

        /*--- Viscous contribution ---*/
      	Lambda_1 = (4.0/3.0)*(Mean_LaminarVisc);
      	Lambda_2 = Mean_ThermalCond/Mean_CV;
      	Lambda = (Lambda_1 + Lambda_2)*Area*Area/Mean_Density;
      	if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Visc(Lambda);
      	if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->AddMax_Lambda_Visc(Lambda);
      }
    }
  }

  /*--- Each element uses their own speed for a steady state simulation ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    Volume = geometry->node[iPoint]->GetVolume();

    if(std::abs(Volume) > EPS) {
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
	for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		if (geometry->node[iPoint]->GetnPoint() == 1)
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
      if(std::abs(config->GetCFL(iMesh)) < EPS)
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      else
        node[iPoint]->SetDelta_Time(Global_Delta_Time);

    }
  }

  if(dual_time)
    throw Common::NotImplemented("Dual time strategies are not implemented");
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
  static_cast<CReactiveEulerSolver*>(this)->SetNondimensionalization(geometry,config,iMesh);

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

  //Viscosity_FreeStream = library->ComputeViscosity(config->GetTemperature_FreeStream(),config->GetMassFrac_FreeStream());
  config->SetViscosity_FreeStream(Viscosity_FreeStream);

  Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;
  config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

  /*--- Constant viscosity model ---*/
  //config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

  /* constant thermal conductivity model */
  //config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

  /*--- Write output to the console if this is the master node and first domain ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH && rank == MASTER_NODE && iMesh == MESH_0) {
    std::cout.precision(6);

    std::cout<< "Reference viscosity: " << config->GetViscosity_Ref();
    if(config->GetSystemMeasurements() == SI)
      std::cout << " N.s/m^2."<< std::endl;
    else if(config->GetSystemMeasurements() == US)
      std::cout<< " lbf.s/ft^2."<< std::endl;
    else
      throw std::out_of_range("Unknown option in the unit measure");

    std::cout << "Reference conductivity: " << config->GetConductivity_Ref();
    if(config->GetSystemMeasurements() == SI)
      std::cout << " W/m^2.K."<< std::endl;
    else if (config->GetSystemMeasurements() == US)
      std::cout << " lbf/ft.s.R."<< std::endl;
    else
      throw std::out_of_range("Unknown option in the unit measure");

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

  /*--- Initialize the viscous residual to zero ---*/
  SU2_Assert(Res_Visc != NULL,"The array for viscous residual has not been allocated");
  std::fill(Res_Visc,Res_Visc + nVar,0.0);

  for (iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),geometry->node[jPoint]->GetCoord() );
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Set Primitive variables and gradient ---*/
    numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
    static_cast<CAvgGradReactive_Flow*>(numerics)->SetGradient_AvgPrimitive(static_cast<CReactiveNSVariable*>(node[iPoint])->GetAvgGradient(),
                                                                            static_cast<CReactiveNSVariable*>(node[jPoint])->GetAvgGradient());

    /*--- Species diffusion coefficients ---*/
    numerics->SetDiffusionCoeff(node[iPoint]->GetDiffusionCoeff(), node[jPoint]->GetDiffusionCoeff());

    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(node[iPoint]->GetLaminarViscosity(), node[jPoint]->GetLaminarViscosity());

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(node[iPoint]->GetThermalConductivity(), node[jPoint]->GetThermalConductivity());

    /*--- Compute the residual ---*/
    try {
      numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    }
    catch(const std::exception& e) {
      std::cout<<e.what()<<std::endl;
      static_cast<CAvgGradReactive_Flow*>(numerics)->SetExplicit();
      numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);
    }

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    bool err = false;
    err = !std::none_of(Res_Visc,Res_Visc + nVar,[](su2double elem){return std::isnan(elem);});
    if(implicit)
      throw Common::NotImplemented("Implicit computation for viscous residual not implemented");
      /*
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
              (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
            err = true;
      */

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      LinSysRes.AddBlock(jPoint, Res_Visc);
      if(implicit)
        throw Common::NotImplemented("Implicit computation for viscous residual not implemented");
        /*
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
        */
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
  /*--- The primitive for reconstruction procedure are the same and so we exploit the Euler version ---*/
  static_cast<CReactiveEulerSolver*>(this)->SetPrimitive_Gradient_GG(geometry,config);

  /*--- Now we foucs on the gradient of primitives required in average gradient to compute the residual ---*/
  unsigned long iPoint, jPoint, iEdge, iVertex;
  unsigned short iDim, jDim, iMarker, iSpecies;
  su2double PrimVar_Average, Partial_Gradient, Partial_Res;

  /*--- Set Gradient Primitive for average to zero for each node and extract the already computed for temperature and velocity ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient_Zero();
    for(iDim = 0; iDim < nDim; ++iDim) {
      static_cast<CReactiveNSVariable*>(node[iPoint])->
      SetAvgGradient(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim,node[iPoint]->GetGradient_Primitive(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim));
      for(jDim = 0; jDim < nDim; ++jDim)
        static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::VX_INDEX_AVGGRAD + iDim,jDim,
                                                         node[iPoint]->GetGradient_Primitive(CReactiveNSVariable::VX_INDEX_AVGGRAD+iDim,jDim));
    }
  }

  RealVec PrimAvgVar_i(nSpecies);
  RealVec PrimAvgVar_j(nSpecies);

  /*--- Loop interior edges ---*/
  for(iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Pull primitives from CVariable ---*/
    std::copy(node[iPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM,
              node[iPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies,
              PrimAvgVar_i.begin());
    std::copy(node[jPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM,
              node[jPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies,
              PrimAvgVar_j.begin());

    /*--- Apply Green-Gauss algorithm ---*/
    auto Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
    // Mixture density
    PrimVar_Average = 0.5 * (node[iPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM) +
                             node[jPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM));
    for(iDim = 0; iDim < nDim; ++iDim) {
      Partial_Res = PrimVar_Average*Normal[iDim];
      if(geometry->node[iPoint]->GetDomain())
        static_cast<CReactiveNSVariable*>(node[iPoint])->AddAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD, iDim, Partial_Res);
      if(geometry->node[jPoint]->GetDomain())
        static_cast<CReactiveNSVariable*>(node[jPoint])->SubtractAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD, iDim, Partial_Res);
    }
    // Partial densities
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    	PrimVar_Average =  0.5 * (PrimAvgVar_i[iSpecies] + PrimAvgVar_j[iSpecies]);
    	for(iDim = 0; iDim < nDim; ++iDim) {
    		Partial_Res = PrimVar_Average*Normal[iDim];
    		if(geometry->node[iPoint]->GetDomain())
    			static_cast<CReactiveNSVariable*>(node[iPoint])->AddAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies, iDim, Partial_Res);
    		if(geometry->node[jPoint]->GetDomain())
    			static_cast<CReactiveNSVariable*>(node[jPoint])->SubtractAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies, iDim, Partial_Res);
    		}
    	}
    }

  /*--- Loop boundary edges ---*/
  for(iMarker = 0; iMarker < geometry->GetnMarker(); ++iMarker) {
    for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); ++iVertex) {
    	iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    	if(geometry->node[iPoint]->GetDomain()) {
        auto Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
        // Mixture density
        PrimVar_Average = 0.5 * (node[iPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM) +
                                 node[jPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM));
        for(iDim = 0; iDim < nDim; ++iDim) {
          Partial_Res = PrimVar_Average*Normal[iDim];
          if(geometry->node[iPoint]->GetDomain())
            static_cast<CReactiveNSVariable*>(node[iPoint])->AddAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD, iDim, Partial_Res);
          if(geometry->node[jPoint]->GetDomain())
            static_cast<CReactiveNSVariable*>(node[jPoint])->SubtractAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD, iDim, Partial_Res);
        }
        // Partial densities
        for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
          PrimVar_Average =  0.5 * (PrimAvgVar_i[iSpecies] + PrimAvgVar_j[iSpecies]);
          for(iDim = 0; iDim < nDim; ++iDim) {
            Partial_Res = PrimVar_Average*Normal[iDim];
            if(geometry->node[iPoint]->GetDomain())
              static_cast<CReactiveNSVariable*>(node[iPoint])->
              AddAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies, iDim, Partial_Res);
            if(geometry->node[jPoint]->GetDomain())
              static_cast<CReactiveNSVariable*>(node[jPoint])->
              SubtractAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies, iDim, Partial_Res);
          }
        }
    	} /*--- End of GetDomain() ---*/
    } /*--- End of nVertex() ---*/
  } /*--- End of nMarker() ---*/

  /*--- Update gradient value ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    su2double Volume = geometry->node[iPoint]->GetVolume();
  	for(iDim = 0; iDim < nDim; ++iDim) {
  		Partial_Gradient = static_cast<CReactiveNSVariable*>(node[iPoint])->GetAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD,iDim) / Volume;
  		static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD, iDim, Partial_Gradient);
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        Partial_Gradient = static_cast<CReactiveNSVariable*>(node[iPoint])->
                           GetAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies,iDim) / Volume;
        static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies, iDim, Partial_Gradient);
      }
    }
  }
} /*---End of function ---*/

//
//
/*!
 *\brief Set gradient primitive variables according to Green-Gauss
 */
//
//
void CReactiveNSSolver::SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) {
  /*--- The primitive for reconstruction procedure are the same and so we exploit the Euler version ---*/
  static_cast<CReactiveEulerSolver*>(this)->SetPrimitive_Gradient_LS(geometry,config);

  /*--- Now we foucs on the gradient of primitives required in average gradient to compute the residual ---*/
  unsigned short iVar, iDim, jDim, iSpecies, iNeigh;
	unsigned long iPoint, jPoint;
	su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, detR2, z11, z12, z13, z22, z23, z33;
  su2double weight;

  /*--- Set Gradient Primitive for average to zero for each node and extract the already computed for temperature and velocity ---*/
  for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient_Zero();
    for(iDim = 0; iDim < nDim; ++iDim) {
      static_cast<CReactiveNSVariable*>(node[iPoint])->
      SetAvgGradient(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim,node[iPoint]->GetGradient_Primitive(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim));
      for(jDim = 0; jDim < nDim; ++jDim)
        static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::VX_INDEX_AVGGRAD + iDim,jDim,
                                                         node[iPoint]->GetGradient_Primitive(CReactiveNSVariable::VX_INDEX_AVGGRAD+iDim,jDim));
    }
  }

  bool singular;
  RealVec PrimAvgVar_i(nSpecies + 1), /*--- Gradient of the following primitive variables: [T,u,v,w,p]^T ---*/
          PrimAvgVar_j(nSpecies + 1);
  RealMatrix C_Mat(nSpecies+1,nDim);
  RealMatrix S_Mat(nDim,nDim);

	/*--- Loop over points of the grid ---*/
	for(iPoint = 0; iPoint < nPointDomain; ++iPoint) {
		/*--- Set the value of singular ---*/
		singular = false;

    /*--- Get coordinates ---*/
		auto Coord_i = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());

    /*--- Get primitives from CVariable ---*/
    PrimAvgVar_i[0] = node[iPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM);
    std::copy(node[iPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM,
              node[iPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies,
              PrimAvgVar_i.begin() + 1);

    /*--- Inizialization of variables ---*/
		std::fill(C_Mat.begin(),C_Mat.end(),0.0);

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		//AD::StartPreacc();
    //AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    //AD::SetPreaccIn(Coord_i, nDim);

		for(iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      /*--- Get coordinates ---*/
			auto Coord_j = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      /*--- Get primitives from CVariable ---*/
      PrimAvgVar_j[0] = node[jPoint]->GetPrimitive(CReactiveNSVariable::RHO_INDEX_PRIM);
      std::copy(node[jPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM,
                node[jPoint]->GetPrimitive() + CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies,
                PrimAvgVar_j.begin() + 1);

			//AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

      RealVec Coord_ij(nDim);
      std::transform(Coord_j.get(),Coord_j.get() + nDim,Coord_i.get(),Coord_ij.begin(),std::minus<su2double>());

      weight = std::inner_product(Coord_ij.cbegin(),Coord_ij.cend(),Coord_ij.cbegin(),0.0);

			/*--- Sumations for entries of upper triangular matrix R ---*/
      if(std::abs(weight) > EPS){
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
        for(iVar = 0; iVar < nSpecies + 1; ++iVar)
          for(iDim = 0; iDim < nDim; ++iDim)
            C_Mat(iVar,iDim) += Coord_ij[iDim]*(PrimAvgVar_j[iVar] - PrimAvgVar_i[iVar])/weight;
      }
    } /*--- End iNeigh for loop ---*/

		/*--- Entries of upper triangular matrix R ---*/
    if(r11 >= 0.0)
      r11 = std::sqrt(r11);
    else
      r11 = 0.0;

    if(std::abs(r11) > EPS)
      r12 = r12/r11;
    else
      r12 = 0.0;

    if(r22 - r12*r12 >= 0.0)
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
      if(r33-r23*r23-r13*r13 >= 0.0)
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
    if(singular)
      std::fill(S_Mat.begin(),S_Mat.end(),0.0);
    else {
      if(nDim == 2) {
        S_Mat(0,0) = (r12*r12+r22*r22)/detR2;
        S_Mat(0,1) = -r11*r12/detR2;
        S_Mat(1,0) = S_Mat(0,1);
        S_Mat(1,1) = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        S_Mat(0,0) = (z11*z11+z12*z12+z13*z13)/detR2;
        S_Mat(0,1) = (z12*z22+z13*z23)/detR2;
        S_Mat(0,2) = (z13*z33)/detR2;
        S_Mat(1,0) = S_Mat(0,1);
        S_Mat(1,1) = (z22*z22+z23*z23)/detR2;
        S_Mat(1,2) = (z23*z33)/detR2;
        S_Mat(2,0) = S_Mat(0,2);
        S_Mat(2,1) = S_Mat(1,2);
        S_Mat(2,2) = (z33*z33)/detR2;
      }
    }

	  /*--- Computation of the gradient: S*c ---*/
    auto result = C_Mat*S_Mat;
    for(iDim = 0; iDim < nDim; ++iDim) {
      static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::RHO_INDEX_AVGGRAD,iDim,result(0,iDim));
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        static_cast<CReactiveNSVariable*>(node[iPoint])->SetAvgGradient(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies,
                                                                        iDim,result(iSpecies + 1,iDim));
    }
	  //	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
    //  AD::EndPreacc();
	} /*--- End of iPoint for loop ---*/
} /*--- End of the function ---*/
