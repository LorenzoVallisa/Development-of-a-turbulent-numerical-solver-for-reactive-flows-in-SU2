#include "../include/solver_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/move_pointer.hpp"

namespace {
  using SmartArr = CReactiveEulerSolver::SmartArr;
  /*!
   * \brief Compute area for the current normal
   */
  su2double ComputeArea(const SmartArr& Normal,const unsigned short nDim) {
    su2double Area = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
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
CReactiveEulerSolver::CReactiveEulerSolver():CSolver(),nSpecies(),nMarker(),space_centered(),implicit(),least_squares(),Alpha(),Beta(),
                                             Gamma(),Gamma_Minus_One(),Mach_Inf(),Density_Inf(),Pressure_Inf(),Temperature_Inf() {

  IterLinSolver = 0;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;
  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;
  nMarker = 0;
  nPoint  = 0;
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
CReactiveEulerSolver::CReactiveEulerSolver(std::shared_ptr<CGeometry> geometry, std::shared_ptr<CConfig> config,unsigned short iMesh):
                      CSolver(),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim, iMarker;

  IterLinSolver = 0;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nVarGrad = 0;

  Max_Delta_Time = 0.0;
  Min_Delta_Time = 1.E6;

  nMarker = config->GetnMarker_All();
  nPoint  = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nDim = geometry->GetnDim();

  nVar = nSpecies + nDim + 2; /*--- Conserved variables (rho,rho*vx,rho*vy,rho*vz,rho*E,rho1,...rhoNs)^T ---*/
  nPrimVar = nSpecies + nDim + 5; /*--- Primitive variables (T,vx,vy,vz,P,rho,h,a,rho1,...rhoNs)^T ---*/
  nPrimVarGrad = nDim + 2; /*--- Gradient Primitive variables (T,vx,vy,vz,P,rho)^T ---*/

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/
  SetNondimensionalization(geometry.get(), config.get(), iMesh);

	/*--- Store the number of vertices on each marker for deallocation ---*/
	nVertex = std::unique_ptr<unsigned long[]>(new unsigned long[nMarker]);
	for (iMarker = 0; iMarker < nMarker; ++iMarker)
		nVertex[iMarker] = geometry->nVertex[iMarker];

	/*--- Allocate a CVariable array for each node of the mesh ---*/
	node = new CVariable*[nPoint];

	/*--- Define some auxiliary vectors related to the residual ---*/
	Residual = new su2double[nVar];
  std::fill(Residual,Residual + nVar,0.0);

  Residual_RMS = new su2double[nVar];
  std::fill(Residual_RMS,Residual_RMS + nVar,0.0);

  Residual_Max = new su2double[nVar];
  std::fill(Residual_Max,Residual_Max + nVar,0.0);

  Residual_i = new su2double[nVar];
  std::fill(Residual_i,Residual_i + nVar,0.0);

  Residual_j = new su2double[nVar];
  std::fill(Residual_j,Residual_j + nVar,0.0);

  Res_Conv = new su2double[nVar];
  std::fill(Res_Conv,Res_Conv + nVar,0.0);

  //Res_Visc = new su2double[nVar];
  //std::fill(Res_Visc,Res_Visc + nVar,0.0);

  Res_Sour = new su2double[nVar];
  std::fill(Res_Sour,Res_Sour + nVar,0.0);

	/*--- Define some structure for locating max residuals ---*/
	Point_Max = new unsigned long[nVar];
  std::fill(Point_Max,Point_Max + nVar,0.0);

  Point_Max_Coord = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; ++iVar){
		Point_Max_Coord[iVar] = new su2double[nDim];
		for (iDim = 0; iDim < nDim; ++iDim)
      Point_Max_Coord[iVar][iDim] = 0.0;
	}

  /*--- Allocate vectors related to the solution ---*/
  Sol_i = SmartArr(new su2double[nVar]);
  Sol_j = SmartArr(new su2double[nVar]);
  Primitive_i = SmartArr(new su2double[nPrimVar]);
  Primitive_j = SmartArr(new su2double[nPrimVar]);

  /*--- Allocate arrays for conserved variable limits ---*/
  Lower_Limit.resize(nVar);
  Upper_Limit.resize(nVar);

  Upper_Limit[CReactiveEulerVariable::RHO_INDEX_SOL] = 1E16;
  for(iVar = CReactiveEulerVariable::RHOVX_INDEX_SOL; iVar < CReactiveEulerVariable::RHOVX_INDEX_SOL + nDim; ++iVar) {
    Lower_Limit[iVar] = -1E16;
    Upper_Limit[iVar] = 1E16;
  }
  std::fill(Upper_Limit.begin() + CReactiveEulerVariable::RHOE_INDEX_SOL,Upper_Limit.end(),1E16);

  /*--- Initialize the solution & residual CVectors ---*/
 	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Create the structure for storing extra information ---*/
 	if (config->GetExtraOutput()) {
    nOutputVariables = nVar;
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
  }

	/*--- Allocate Jacobians for implicit time-stepping ---*/
	if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {
    implicit = true;
		Jacobian_i = new su2double* [nVar];
		Jacobian_j = new su2double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_i[iVar] = new su2double [nVar];
			Jacobian_j[iVar] = new su2double [nVar];
		}
  }
  else
    implicit = false;

  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    space_centered = true;
  else
    space_centered = false;

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
    least_squares = true;
  else
    least_squares = false;

  Gamma = 0.0;
  Gamma_Minus_One = 0.0;

  Density_Inf      = config->GetDensity_FreeStreamND();
  Pressure_Inf     = config->GetPressure_FreeStreamND();
	Temperature_Inf  = config->GetTemperature_FreeStreamND();
  Mach_Inf         = config->GetMach();

  Velocity_Inf     = Common::wrap_in_unique(config->GetVelocity_FreeStreamND());
  MassFrac_Inf     = Common::wrap_in_unique(config->GetMassFrac_FreeStream());

  Alpha  = config->GetAoA()*PI_NUMBER/180.0;
  Beta   = config->GetAoS()*PI_NUMBER/180.0;

  //node_infty = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf,Velocity_Inf, Temperature_Inf,
  //                                        nDim, nVar,nPrimVar, nPrimVarGrad, config);

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  //for (iPoint = 0; iPoint < nPoint; iPoint++)
    //node[iPoint] = new CReactiveEulerVariable(Pressure_Inf, MassFrac_Inf, Velocity_Inf, Temperature_Inf,
    //                                          nDim, nVar, nPrimVar, nPrimVarGrad,config);

  /*--- Use a function to check that the initial solution is physical ---*/

  Check_Initial_Solution(config);

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
void CReactiveEulerSolver::Check_Initial_Solution(std::shared_ptr<CConfig> config) {

  int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  SmartArr Solution(new su2double[nVar]);

  //bool check_infty = node_infty->SetPrimVar(config.get());

  unsigned long counter_local = 0, counter_global;
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    bool nonPhys = node[iPoint]->SetPrimVar(config.get());

    if (nonPhys) {

      unsigned short iSpecies,iDim;
      su2double Ru, T, rhoCvtr, sqvel, rhoE, denom, conc;
      su2double rho, rhos, Ef, Sound_Speed;
      SmartArr Ms, hf, xi, Tref;
      SmartArr Mvec_Inf(new su2double[nDim]);
      if (nDim == 2) {
        Mvec_Inf[0] = std::cos(Alpha)*Mach_Inf;
        Mvec_Inf[1] = std::sin(Alpha)*Mach_Inf;
      }
      if (nDim == 3) {
        Mvec_Inf[0] = std::cos(Alpha)*std::cos(Beta)*Mach_Inf;
        Mvec_Inf[1] = std::sin(Beta)*Mach_Inf;
        Mvec_Inf[2] = std::sin(Alpha)*std::cos(Beta)*Mach_Inf;
      }

      /*--- Load variables from the config class --*/
      Ms = Common::wrap_in_unique(config->GetMolar_Mass());  // Species molar mass
      hf = Common::wrap_in_unique(config->GetEnthalpy_Formation()); // Formation enthalpy [J/kg]
      xi = Common::wrap_in_unique(config->GetRotationModes());  // Rotation modes
      Tref = Common::wrap_in_unique(config->GetRefTemperature());  // Thermodynamic reference temperature for each species [K]

      /*--- Rename & initialize for convenience ---*/
      Ru      = (*library).R_ungas;         // Universal gas constant [J/(kmol*K)]
      T       = config->GetTemperature_FreeStream();           // Translational-rotational temperature [K]
      sqvel   = 0.0;                            // Velocity^2 [m2/s2]
      rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
      rhoCvtr = 0.0;
      denom   = 0.0;
      conc    = 0.0;

      /*--- Calculate mixture density from supplied primitive quantities ---*/
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        denom += MassFrac_Inf[iSpecies] * (Ru/Ms[iSpecies]) * T;
      rho = Pressure_Inf / denom;

      /*--- Calculate sound speed and extract velocities ---*/
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        conc += MassFrac_Inf[iSpecies]*rho/Ms[iSpecies];
        rhoCvtr += rho*MassFrac_Inf[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
      }
      Sound_Speed = std::sqrt((1.0 + Ru/rhoCvtr*conc) * Pressure_Inf/rho);
      for (iDim = 0; iDim < nDim; ++iDim)
        sqvel += Mvec_Inf[iDim]*Sound_Speed*Mvec_Inf[iDim]*Sound_Speed;

      /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
      for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        // Species density
        rhos = MassFrac_Inf[iSpecies]*rho;

        // Species formation energy
        Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

        // Mixture total energy
        rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies]) + Ef + 0.5*sqvel);

        // Set in the solution
        Solution[CReactiveEulerVariable::RHOS_INDEX_SOL+iSpecies] = rhos;
      }

      Solution[CReactiveEulerVariable::RHO_INDEX_SOL] = rho;
      for (iDim = 0; iDim < nDim; ++iDim)
        Solution[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = rho*Mvec_Inf[iDim]*Sound_Speed;
      Solution[CReactiveEulerVariable::RHOE_INDEX_SOL] = rhoE;

      node[iPoint]->SetSolution(Solution.get());
      node[iPoint]->SetSolution_Old(Solution.get());

      counter_local++;
    }
  }

  /*--- Warning message about non-physical points ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
  #ifdef HAVE_MPI
    SU2_MPI::Reduce(&counter_local, &counter_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  #else
    counter_global = counter_local;
  #endif
    if (rank == MASTER_NODE and counter_global != 0)
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
             ModVel_FreeStreamND = 0.0, Density_FreeStream = 0.0, Pressure_FreeStream = 0.0,
             Mach2Vel_FreeStream = 0.0;
             //Viscosity_FreeStream;
   su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Velocity_Ref = 0.0,
             Temperature_Ref = 0.0, Time_Ref = 0.0, Gas_Constant_Ref = 0.0, Energy_Ref = 0.0;
             //Viscosity_Ref = 0.0,  //Conductivity_Ref = 0.0;
   su2double Pressure_FreeStreamND = 0.0, Density_FreeStreamND = 0.0,
             Temperature_FreeStreamND = 0.0, Gas_Constant = 0.0, Gas_ConstantND = 0.0,
             Energy_FreeStreamND = 0.0, Total_UnstTimeND = 0.0, Delta_UnstTimeND = 0.0;
             //Viscosity_FreeStreamND;

   SmartArr  Velocity_FreeStreamND(new su2double[nDim]);

   unsigned short iDim;

   /*--- Local variables ---*/

   bool unsteady   = config->GetUnsteady_Simulation();
   su2double Mach  = config->GetMach();
   //bool viscous    = config->GetViscous();

   int rank = MASTER_NODE;
   #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   #endif

 	 /*--- Compressible non dimensionalization ---*/

   /*--- Compute the Free Stream velocity, using the Mach number ---*/

   Pressure_FreeStream     = config->GetPressure_FreeStream();
   Temperature_FreeStream  = config->GetTemperature_FreeStream();
   Gas_Constant  				  = library->GetRgas();
   Gamma = library->Gamma();
   Gamma_Minus_One = Gamma - 1.0;
 	 Mach2Vel_FreeStream = std::sqrt(Gamma*Gas_Constant*Temperature_FreeStream);

   /*--- Compute the Free Stream velocity, using the Mach number ---*/

   if (nDim == 2) {
     config->GetVelocity_FreeStream()[0] = std::cos(Alpha)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[1] = std::sin(Alpha)*Mach*Mach2Vel_FreeStream;
   }
   if (nDim == 3) {
     config->GetVelocity_FreeStream()[0] = std::cos(Alpha)*std::cos(Beta)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[1] = std::sin(Beta)*Mach*Mach2Vel_FreeStream;
     config->GetVelocity_FreeStream()[2] = std::sin(Alpha)*std::cos(Beta)*Mach*Mach2Vel_FreeStream;
   }


   /*--- Compute the modulus of the free stream velocity ---*/

   ModVel_FreeStream = 0.0;
   for (iDim = 0; iDim < nDim; ++iDim)
     ModVel_FreeStream += config->GetVelocity_FreeStream()[iDim]*config->GetVelocity_FreeStream()[iDim];
   ModVel_FreeStream = std::sqrt(ModVel_FreeStream);
   config->SetModVel_FreeStream(ModVel_FreeStream);


 	/*--- Compute the density ---*/

  Density_FreeStream = Pressure_FreeStream/(Gas_Constant*Temperature_FreeStream);
  config->SetDensity_FreeStream(Density_FreeStream);


   /*--- Viscous initialization ---*/

   //if (viscous) {

     /*--- The dimensional viscosity is needed to determine the free-stream conditions.
           To accomplish this, simply set the non-dimensional coefficients to the
           dimensional ones. This will be overruled later.---*/
     /*
     config->SetMu_RefND(config->GetMu_Ref());
     config->SetMu_Temperature_RefND(config->GetMu_Temperature_Ref());
     config->SetMu_SND(config->GetMu_S());

     config->SetMu_ConstantND(config->GetMu_Constant());

    }



   }
   else {
   */
   /*--- For inviscid flow, energy is calculated from the specified
        FreeStream quantities using the proper gas law. ---*/
   Energy_FreeStream = Pressure_FreeStream/(Density_FreeStream*Gamma_Minus_One)+ 0.5*ModVel_FreeStream*ModVel_FreeStream ;

  // }

   /*--- Compute non dimensional quantities: Notice that the grid is in meters. ---*/

   if (config->GetRef_NonDim() == DIMENSIONAL) {
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
   else if (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE) {
     Pressure_Ref      = Pressure_FreeStream;     // Pressure_FreeStream = 1.0
     Density_Ref       = Density_FreeStream;      // Density_FreeStream = 1.0
     Temperature_Ref   = Temperature_FreeStream;  // Temperature_FreeStream = 1.0
   }
   else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH) {
     Pressure_Ref      = Gamma*Pressure_FreeStream; // Pressure_FreeStream = 1.0/Gamma
     Density_Ref       = Density_FreeStream;        // Density_FreeStream = 1.0
     Temperature_Ref   = Temperature_FreeStream;    // Temp_FreeStream = 1.0
   }
   else if (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE) {
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

   //Viscosity_Ref = Density_Ref*Velocity_Ref*Length_Ref;
   //config->SetViscosity_Ref(Viscosity_Ref);

   //Conductivity_Ref = Viscosity_Ref*Gas_Constant_Ref;
   //config->SetConductivity_Ref(Conductivity_Ref);

   /*--- Divide by reference values, to compute the non-dimensional free-stream values ---*/

   Pressure_FreeStreamND = Pressure_FreeStream/Pressure_Ref;
   config->SetPressure_FreeStreamND(Pressure_FreeStreamND);

   Density_FreeStreamND  = Density_FreeStream/Density_Ref;
   config->SetDensity_FreeStreamND(Density_FreeStreamND);

   for (iDim = 0; iDim < nDim; ++iDim) {
     Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref;
     config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
   }

   Temperature_FreeStreamND = Temperature_FreeStream/Temperature_Ref;
   config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

   Energy_FreeStreamND = Temperature_FreeStream/Energy_Ref;
   config->SetEnergy_FreeStreamND(Energy_FreeStreamND);

   Gas_ConstantND = Gas_Constant/Gas_Constant_Ref;
   config->SetGas_ConstantND(Gas_ConstantND);

   ModVel_FreeStreamND = 0.0;
   for (iDim = 0; iDim < nDim; ++iDim)
    ModVel_FreeStreamND += Velocity_FreeStreamND[iDim]*Velocity_FreeStreamND[iDim];
   ModVel_FreeStreamND  = std::sqrt(ModVel_FreeStreamND);
   config->SetModVel_FreeStreamND(ModVel_FreeStreamND);

   //Viscosity_FreeStreamND = Viscosity_FreeStream / Viscosity_Ref;   config->SetViscosity_FreeStreamND(Viscosity_FreeStreamND);

   /*--- Initialize the dimensionless Fluid Model that will be used to solve the dimensionless problem ---*/

   //if (viscous) {

     /*--- Constant viscosity model ---*/
     //config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);

     /* constant thermal conductivity model */
     //config->SetKt_ConstantND(config->GetKt_Constant()/Conductivity_Ref);

   //}

   Total_UnstTimeND = config->GetTotal_UnstTime()/Time_Ref;
   config->SetTotal_UnstTimeND(Total_UnstTimeND);

   Delta_UnstTimeND = config->GetDelta_UnstTime()/Time_Ref;
   config->SetDelta_UnstTimeND(Delta_UnstTimeND);

   /*--- Write output to the console if this is the master node and first domain ---*/

   if (config->GetConsole_Output_Verb() == VERB_HIGH and rank == MASTER_NODE and iMesh == MESH_0) {

     std::cout.precision(6);

     /*--- Print out reference values. ---*/

     std::cout <<"-- Reference values:"<< std::endl;

     bool SI_Measurement = config->GetSystemMeasurements() == SI;
     bool US_Measuremanet = config->GetSystemMeasurements() == US;

     std::cout << "Reference specific gas constant: " << Gas_Constant_Ref;
     if (SI_Measurement)
       std::cout << " N.m/kg.K." << std::endl;
     else if (US_Measuremanet)
       std::cout << " lbf.ft/slug.R." << std::endl;

     std::cout << "Reference pressure: " << Pressure_Ref;
     if (SI_Measurement)
       std::cout << " Pa." << std::endl;
     else if (US_Measuremanet)
       std::cout << " psf." << std::endl;

     std::cout << "Reference temperature: " << Temperature_Ref;
     if (SI_Measurement)
       std::cout << " K." << std::endl;
     else if (US_Measuremanet)
       std::cout << " R." << std::endl;

     std::cout << "Reference density: " << Density_Ref;
     if (SI_Measurement)
       std::cout << " kg/m^3." << std::endl;
     else if (US_Measuremanet)
       std::cout << " slug/ft^3." << std::endl;

     std::cout << "Reference velocity: " << Velocity_Ref;
     if (SI_Measurement)
       std::cout << " m/s." << std::endl;
     else if (US_Measuremanet)
       std::cout << " ft/s." << std::endl;

     std::cout << "Reference energy per unit mass: " << Energy_Ref;
     if (SI_Measurement)
       std::cout << " m^2/s^2." << std::endl;
     else if (US_Measuremanet)
       std::cout << " ft^2/s^2." << std::endl;

    /*
     if (viscous) {
       std::cout << "Reference viscosity: " << config->GetViscosity_Ref();
       if (config->GetSystemMeasurements() == SI)
         std::cout << " N.s/m^2." << std::endl;
       else if (config->GetSystemMeasurements() == US)
         std::cout << " lbf.s/ft^2." << std::endl;

       std::cout << "Reference conductivity: " << config->GetConductivity_Ref();
       if (config->GetSystemMeasurements() == SI)
         std::cout << " W/m^2.K." << std::endl;
       else if (config->GetSystemMeasurements() == US)
         std::cout << " lbf/ft.s.R." << std::endl;
     }
     */

     if (unsteady)
       std::cout << "Reference time: " << Time_Ref <<" s." << std::endl;

     /*--- Print out resulting non-dim values here. ---*/

     std::cout << "-- Resulting non-dimensional state:" << std::endl;
     std::cout << "Mach number (non-dim): " << config->GetMach() << std::endl;

     //if (viscous)
      //std::cout << "Reynolds number (non-dim): " << config->GetReynolds() <<std::endl;

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

     //if (viscous)
      // std::cout << "Free-stream viscosity (non-dim): " << config->GetViscosity_FreeStreamND() << std::endl;

     if (unsteady) {
       std::cout << "Total time (non-dim): " << Total_UnstTimeND << std::endl;
       std::cout << "Time step (non-dim): " << Delta_UnstTimeND << std::endl;
     }

     std::cout << std::endl;

   }

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
void CReactiveEulerSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                         unsigned short iMesh, unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem, bool Output) {

  unsigned long ErrorCounter = 0;

  unsigned long ExtIter = config->GetExtIter();
  bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  //bool low_fidelity     = (config->GetLowFidelitySim() && (iMesh == MESH_1));
  bool second_order     = (config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER);
  bool limiter          = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER and ExtIter <= config->GetLimiterIter());
                          //and !low_fidelity;
  bool center           = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;

  /*--- Set the primitive variables ---*/
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

  /*--- Upwind second order reconstruction ---*/

  if ((second_order and !center) and iMesh == MESH_0 and !Output) {

    /*--- Gradient computation ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
      SetPrimitive_Gradient_GG(geometry, config);
    }
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
      SetPrimitive_Gradient_LS(geometry, config);
    }

    /*--- Limiter computation ---*/

    if (limiter and iMesh == MESH_0 and !Output)
      SetPrimitive_Limiter(geometry, config);
  }

  /*--- Artificial dissipation ---*/

  if (center and !Output)
      throw Common::NotImplemented("Centered convective scheme not implemented\n");

  /*--- Initialize the Jacobian matrices ---*/

  if (implicit)
    Jacobian.SetValZero();

  /*--- Error message ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
    #ifdef HAVE_MPI
      unsigned long MyErrorCounter = ErrorCounter;
      ErrorCounter = 0;
      SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    #endif
    if (iMesh == MESH_0)
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
  unsigned short iDim, iSpecies, iVar, iMarker;
  su2double PrimVar_Average, Partial_Gradient, Partial_Res;
  su2double rho_i,rho_j;
  RealVec PrimVar_Vertex(nPrimVarGrad); /*---Primitive variables: [T,u,v,w,Y1, ..., YNs]^T ---*/
  RealVec PrimVar_i(nPrimVarGrad);
  RealVec PrimVar_j(nPrimVarGrad);
  SmartArr Normal;

  /*--- Set Gradient_Primitive to zero ---*/
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint)
  	node[iPoint]->SetGradient_PrimitiveZero(nPrimVarGrad);

  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {
  	iPoint = geometry->edge[iEdge]->GetNode(0);
  	jPoint = geometry->edge[iEdge]->GetNode(1);

  /*--- Pull primitives from CVariable ---*/
    for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  		PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);
  		PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);
  	}

  	Normal = Common::wrap_in_unique(geometry-> edge[iEdge]->GetNormal());
  	for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  		PrimVar_Average =  0.5 * (PrimVar_i[iVar] + PrimVar_j[iVar]);
  		for (iDim = 0; iDim < nDim; ++iDim) {
  			Partial_Res = PrimVar_Average*Normal[iDim];
  			if (geometry->node[iPoint]->GetDomain())
  				node[iPoint]->AddGradient_Primitive(iVar, iDim, Partial_Res);
  			if (geometry->node[jPoint]->GetDomain())
  				node[jPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  		}
  	}
  }

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
  	for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
  		iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
  		if (geometry->node[iPoint]->GetDomain()) {

        /*--- Get primitives from CVariable ---*/
  			for (iVar = 0; iVar < nPrimVarGrad; ++iVar)
  				PrimVar_Vertex[iVar] = node[iPoint]->GetPrimitive(iVar);

        /*--- Modify species density to mass concentration ---*/
        rho_i = node[iPoint]->GetPrimitive(CReactiveEulerVariable::RHO_INDEX_PRIM);
        for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          PrimVar_Vertex[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] /= rho_i;

  			Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
  			for (iVar = 0; iVar < nPrimVarGrad; ++iVar)
  				for (iDim = 0; iDim < nDim; ++iDim) {
  					Partial_Res = PrimVar_Vertex[iVar]*Normal[iDim];
  					node[iPoint]->SubtractGradient_Primitive(iVar, iDim, Partial_Res);
  				}
  		  }
  		}
  	}

  	/*--- Update gradient value ---*/
  	for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
  		for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
  			for (iDim = 0; iDim < nDim; ++iDim) {
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
  su2double rho_i, rho_j, weight;
  bool singular;
  RealVec PrimVar_i(nPrimVarGrad), PrimVar_j(nPrimVarGrad); /*---Primitive variables: [T,u,v,w,Y1, ..., YNs]^T ---*/
  SmartArr Coord_i, Coord_j;

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

		/*--- Set the value of singulare ---*/
		singular = false;

    /*--- Get coordinates ---*/
		Coord_i = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());

    /*--- Get primitives from CVariable ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			PrimVar_i[iVar] = node[iPoint]->GetPrimitive(iVar);

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nPrimVarGrad; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				Cvector[iVar][iDim] = 0.0;

		r11 = 0.0; r12   = 0.0; r13   = 0.0; r22 = 0.0;
		r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		//AD::StartPreacc();
    //AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
    //AD::SetPreaccIn(Coord_i, nDim);

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); ++iNeigh) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

			for (iVar = 0; iVar < nPrimVarGrad; iVar++)
				PrimVar_j[iVar] = node[jPoint]->GetPrimitive(iVar);

			//AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

			weight = 0.0;
			for (iDim = 0; iDim < nDim; ++iDim)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

			/*--- Sumations for entries of upper triangular matrix R ---*/
      if (weight != 0.0){
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }

        /*--- Entries of c:= transpose(A)*b ---*/
        for (iVar = 0; iVar < nPrimVarGrad; ++iVar)
          for (iDim = 0; iDim < nDim; ++iDim)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim]) *
                                   (PrimVar_j[iVar]-PrimVar_i[iVar])/weight;
      }
    }

		/*--- Entries of upper triangular matrix R ---*/
    if (r11 >= 0.0)
      r11 = sqrt(r11);
    else
      r11 = 0.0;

    if (r11 != 0.0)
      r12 = r12/r11;
    else
      r12 = 0.0;

    if (r22-r12*r12 >= 0.0)
      r22 = sqrt(r22-r12*r12);
    else
      r22 = 0.0;

    if (nDim == 3) {
      if (r11 != 0.0)
        r13 = r13/r11;
      else
        r13 = 0.0;

      if (r22 != 0.0 and r11*r22 != 0.0)
        r23 = r23_a/r22 - r23_b*r12/(r11*r22);
      else
        r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0)
        r33 = sqrt(r33-r23*r23-r13*r13);
      else r33 = 0.0;
    }

    /*--- Compute determinant ---*/
    if (nDim == 2)
      detR2 = (r11*r22)*(r11*r22);
    else
      detR2 = (r11*r22*r33)*(r11*r22*r33);

    /*--- Detect singular matrices ---*/
    if (std::abs(detR2) <= EPS) {
      detR2 = 1.0;
      singular = true;
    }

		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    if (singular) {
      for (iDim = 0; iDim < nDim; ++iDim)
        for (jDim = 0; jDim < nDim; ++jDim)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
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
		for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
			for (iDim = 0; iDim < nDim; ++iDim) {
				product = 0.0;
				for (jDim = 0; jDim < nDim; ++jDim)
					product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
				node[iPoint]->SetGradient_Primitive(iVar, iDim, product);
			}
		}

	//	AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
  //  AD::EndPreacc();
	}

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
  SmartArr   Primitive, Primitive_i, Primitive_j;
  SmartArr   Coord_i, Coord_j;
//  RealMatrix Gradient_i, Gradient_j;

  dave = config->GetRefElemLength();
  LimK = config->GetLimiterCoeff();

  /*if (config->GetKind_SlopeLimit_Flow() == NO_LIMITER) {

    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint){
      for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        node[iPoint]->SetLimiter_Primitive(iVar,1.0);
      }
    }
  }
  */
  //else {

    /*--- Initialize solution max, solution min and limiter in entire domain ---*/
    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter_Primitive(iVar, 2.0);
      }
    }

    /*--- Establish bounts for Spekreijse monotonicity by finding max/min values
          of neighbor variables ---*/
    for (iEdge = 0; iEdge < geometry-> GetnEdge(); ++iEdge) {

      /*--- Point identification, Normal vector and area ---*/
      iPoint = geometry-> edge[iEdge]->GetNode(0);
      jPoint = geometry-> edge[iEdge]->GetNode(1);

      /*--- Get primitive variables ---*/
      Primitive_i = Common::wrap_in_unique(node[iPoint]->GetPrimitive());
      Primitive_j = Common::wrap_in_unique(node[jPoint]->GetPrimitive());

      /*--- Compute the max and min values for nodes i & j ---*/
      for (iVar = 0; iVar < nPrimVarGrad; ++iVar){
        du = (Primitive_j[iVar]-Primitive_i[iVar]);
        node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
        node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
        node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
        node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
      }
    }
  //}

  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());
      Coord_j    = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

        //AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        //AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

        if (dm == 0.0)
          limiter = 2.0;
        else {
          if(dm > 0.0)
            dp = node[iPoint]->GetSolution_Max(iVar);
          else
            dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
          node[iPoint]->SetLimiter_Primitive(iVar, limiter);
          //AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
        }

        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

        if (dm == 0.0)
          limiter = 2.0;
        else {
          if (dm > 0.0)
            dp = node[jPoint]->GetSolution_Max(iVar);
          else
            dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }

        if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
          node[jPoint]->SetLimiter_Primitive(iVar, limiter);
          //AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
        }

      }

      //AD::EndPreacc();

    }

    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        y =  node[iPoint]->GetLimiter_Primitive(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
      }
    }

  }

  /*--- Venkatakrishnan limiter ---*/
  if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN) {

    /*--- Allocate memory for the max and min primitive value --*/
    RealVec LocalMinPrimitive(nPrimVarGrad), GlobalMinPrimitive(nPrimVarGrad),
            LocalMaxPrimitive(nPrimVarGrad), GlobalMaxPrimitive(nPrimVarGrad);

    /*--- Compute the max value and min value of the solution ---*/
    Primitive = Common::wrap_in_unique(node[0]->GetPrimitive());
    for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      LocalMinPrimitive[iVar] = Primitive[iVar];
      LocalMaxPrimitive[iVar] = Primitive[iVar];
    }

    for (iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {

      /*--- Get the primitive variables ---*/
      Primitive = Common::wrap_in_unique(node[iPoint]->GetPrimitive());

      for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
        LocalMinPrimitive[iVar] = std::min(LocalMinPrimitive[iVar], Primitive[iVar]);
        LocalMaxPrimitive[iVar] = std::max(LocalMaxPrimitive[iVar], Primitive[iVar]);
      }

    }

    #ifdef HAVE_MPI
      SU2_MPI::Allreduce(LocalMinPrimitive, GlobalMinPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(LocalMaxPrimitive, GlobalMaxPrimitive, nPrimVarGrad, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #else
    for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
      GlobalMinPrimitive[iVar] = LocalMinPrimitive[iVar];
      GlobalMaxPrimitive[iVar] = LocalMaxPrimitive[iVar];
    }
    #endif

    for (iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {

      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      auto Gradient_i = node[iPoint]->GetGradient_Primitive();
      auto Gradient_j = node[jPoint]->GetGradient_Primitive();
      Coord_i    = Common::wrap_in_unique(geometry->node[iPoint]->GetCoord());
      Coord_j    = Common::wrap_in_unique(geometry->node[jPoint]->GetCoord());

      //AD::StartPreacc();
      //AD::SetPreaccIn(Gradient_i, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Gradient_j, nPrimVarGrad, nDim);
      //AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
      //AD::SetPreaccIn(eps2);

      eps1 = LimK*dave;
      eps2 = eps1*eps1*eps1;

      //AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
      //AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
      //AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
      //AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

      /*--- Calculate the interface left gradient, delta- (dm) ---*/
      dm = 0.0;
      for (iDim = 0; iDim < nDim; ++iDim)
        dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];

      /*--- Calculate the interface right gradient, delta+ (dp) ---*/
      if (dm > 0.0)
        dp = node[iPoint]->GetSolution_Max(iVar);
      else
        dp = node[iPoint]->GetSolution_Min(iVar);

      limiter = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

      if (limiter < node[iPoint]->GetLimiter_Primitive(iVar)) {
        node[iPoint]->SetLimiter_Primitive(iVar, limiter);
        //AD::SetPreaccOut(node[iPoint]->GetLimiter_Primitive()[iVar]);
      }

      /*-- Repeat for point j on the edge ---*/
      dm = 0.0;
      for (iDim = 0; iDim < nDim; ++iDim)
        dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];

      if (dm > 0.0)
        dp = node[jPoint]->GetSolution_Max(iVar);
      else
        dp = node[jPoint]->GetSolution_Min(iVar);

      limiter = (dp*dp + 2.0*dp*dm + eps2)/(dp*dp + dp*dm + 2.0*dm*dm + eps2);

      if (limiter < node[jPoint]->GetLimiter_Primitive(iVar)) {
        node[jPoint]->SetLimiter_Primitive(iVar, limiter);
        //AD::SetPreaccOut(node[jPoint]->GetLimiter_Primitive()[iVar]);
      }

    }

    //AD::EndPreacc();

    }

  /*--- Limiter MPI ---*/
  //Set_MPI_Primitive_Limiter(geometry, config);

}

//
//
/*!
 *\brief Set primitive variables
 */
//
//

unsigned long CReactiveEulerSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies+nDim+5, (T,vx,vy,vz,P,rho,h,c,rho1,...rhoNs) ---*/

    RightSol = node[iPoint]->SetPrimVar();

    if (!RightSol) {
      node[iPoint]->SetNon_Physical(true);
      ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output)
      LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;

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

   su2double Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda;
   su2double Local_Delta_Time, Global_Delta_Time = 1E6;
   SmartArr Normal;
   unsigned long iEdge, iVertex, iPoint, jPoint;
   unsigned short iDim, iMarker;

   bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;

   Min_Delta_Time = 1.E6;
   Max_Delta_Time = 0.0;

   /*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

   for (iPoint = 0; iPoint < nPointDomain; iPoint++)
     node[iPoint]->SetMax_Lambda_Inv(0.0);


  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/

    Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
    Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;

    /*--- Inviscid contribution ---*/
    Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;

    if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);
    if (geometry->node[jPoint]->GetDomain())
      node[jPoint]->AddMax_Lambda_Inv(Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- Point identification, Normal vector and area ---*/

        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
        Area = ::ComputeArea(Normal,nDim);

        /*--- Mean Values ---*/

        Mean_ProjVel = node[iPoint]->GetProjVel(Normal.get());
        Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

        /*--- Inviscid contribution ---*/
        Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddMax_Lambda_Inv(Lambda);
      }
    }
  }

  /*--- Each element uses their own speed for a steady state simulation ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();

    if (Vol != 0.0) {
      Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
      Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
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

  if(time_steping) {
    #ifdef HAVE_MPI
      su2double rbuf_time, sbuf_time;
      sbuf_time = Global_Delta_Time;
      SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      Global_Delta_Time = rbuf_time;
    #endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Sets the regular CFL equal to the unsteady CFL ---*/
      config->SetCFL(iMesh,config->GetUnst_CFL());

      /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
            it computes the time step based on the unsteady CFL ---*/
      if (config->GetCFL(iMesh) == 0.0)
        node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
      else
        node[iPoint]->SetDelta_Time(Global_Delta_Time);

    }

  }

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
	SmartArr local_Res_TruncError;

  /*--- Set maximum residual to zero ---*/
    for (iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/
  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    /*--- Read the residual ---*/
    local_Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    /*--- Read the volume ---*/
    auto Vol = geometry->node[iPoint]->GetVolume();

    if (node[iPoint]->GetDelta_Time() != 0.0) {
      auto Delta = Vol / node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; ++iVar) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; ++iVar) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, std::abs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = nPointDomain; iPoint < nPoint; ++iPoint) {
    for (iVar = 0; iVar < nVar; ++iVar) {
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
	for (iPoint = 0; iPoint < nPointDomain; ++iPoint)
		for (iVar = 0; iVar < nVar; ++iVar)
			node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);


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
  SmartArr Residual, Res_TruncError;
  unsigned short iVar;
  unsigned long iPoint;

  for (iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    auto Vol = geometry->node[iPoint]->GetVolume();
    auto Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for (iVar = 0; iVar < nVar; ++iVar) {
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
  SmartArr Residual, Res_TruncError;
  unsigned short iVar;
  unsigned long iPoint;

  auto RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);

  for (iVar = 0; iVar < nVar; ++iVar) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Update the solution ---*/

  for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    auto Vol = geometry->node[iPoint]->GetVolume();
    auto Delta = node[iPoint]->GetDelta_Time() / Vol;

    Res_TruncError = Common::wrap_in_unique(node[iPoint]->GetResTruncError());
    Residual = Common::wrap_in_unique(LinSysRes.GetBlock(iPoint));

    for (iVar = 0; iVar < nVar; ++iVar) {
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

    unsigned long iEdge, iPoint, jPoint;
    unsigned short iDim, iVar;

    SmartArr Limiter_i,Limiter_j;
    su2double Project_Grad_i, Project_Grad_j;

    bool implicit     = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
    //bool low_fidelity = (config->GetLowFidelitySim() && (iMesh == MESH_1));
    bool second_order = ((config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)
                        and iMesh == MESH_0);
    bool limiter      = config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER;

    /*--- Loop over all the edges ---*/

    for (iEdge = 0; iEdge < geometry->GetnEdge(); ++iEdge) {

      /*--- Points in edge and normal vectors ---*/

      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

      /*--- Get primitive variables ---*/

      Primitive_i = Common::wrap_in_unique(node[iPoint]->GetPrimitive());
      Primitive_j = Common::wrap_in_unique(node[jPoint]->GetPrimitive());

      if (second_order) {

        for (iDim = 0; iDim < nDim; ++iDim) {
          Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
          Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
        }

      //  node[iPoint]->GetGradient_Primitive();
      //  node[jPoint]->GetGradient_Primitive();

        if (limiter) {
          Limiter_i = Common::wrap_in_unique(node[iPoint]->GetLimiter_Primitive());
          Limiter_j = Common::wrap_in_unique(node[jPoint]->GetLimiter_Primitive());
        }

        for (iVar = 0; iVar < nPrimVarGrad; ++iVar) {
          Project_Grad_i = 0.0;
          Project_Grad_j = 0.0;
          for (iDim = 0; iDim < nDim; ++iDim) {
            Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
            Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
          }
          if (limiter) {
            Primitive_i[iVar] += Limiter_i[iVar]*Project_Grad_i;
            Primitive_j[iVar] += Limiter_j[iVar]*Project_Grad_j;
          }
          else {
            Primitive_i[iVar] += Project_Grad_i;
            Primitive_j[iVar] += Project_Grad_j;
          }
        }
        numerics->SetPrimitive(Primitive_i.get(), Primitive_j.get());
      }
      else {

        /*--- Set primitive variables without reconstruction ---*/

        numerics->SetPrimitive(Primitive_i.get(), Primitive_j.get());

      }
      /*--- Set conserative variables ---*/

      numerics->SetConservative(node[iPoint]->GetSolution(),node[jPoint]->GetSolution());

      /*--- Compute the residual ---*/

      numerics->ComputeResidual(Res_Conv, Jacobian_i, Jacobian_j, config);

      /*--- Update residual value ---*/

      LinSysRes.AddBlock(iPoint, Res_Conv);
      LinSysRes.SubtractBlock(jPoint, Res_Conv);

      /*--- Set implicit Jacobians ---*/

      if (implicit) {
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
      }
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

unsigned long iPoint;
unsigned short iVar;

bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

/*--- Initialize the source residual to zero ---*/
std::fill(Residual,Residual + nVar,0.0);

/*--- Loop over all points ---*/
for (iPoint = 0; iPoint < nPointDomain; ++iPoint) {

  /*--- Load the conservative variables ---*/
  numerics->SetConservative(node[iPoint]->GetSolution(),node[iPoint]->GetSolution());

  /*--- Load the volume of the dual mesh cell ---*/
  numerics->SetVolume(geometry->node[iPoint]->GetVolume());

  /*--- Compute the rotating frame source residual ---*/
  numerics->ComputeResidual(Residual, Jacobian_i, config);

  /*--- Add the source residual to the total ---*/
  LinSysRes.AddBlock(iPoint, Residual);

  /*--- Add the implicit Jacobian contribution ---*/
  if (implicit)
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

}

}

//
//
/*--- Euler_wall_BC ---*/
//
//

void CReactiveEulerSolver::BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* solver,
                                         CConfig* config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
	unsigned short iDim, iVar;
	su2double Pressure;
  SmartArr Normal;

	bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			Normal = Common::wrap_in_unique(geometry->vertex[val_marker][iVertex]->GetNormal());
      su2double Area = ::ComputeArea(Normal,nDim);

      RealVec UnitaryNormal(3);
			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

				//Pressure = library->node[iPoint]->GetPressure(iSpecies);

				//Residual[0] = 0.0;
				for (iDim = 0; iDim < nDim; iDim++)
					Residual[iDim+1] =  Pressure*UnitaryNormal[iDim]*Area;

        //Residual[nDim+1] = 0.0;
				//Residual[nDim+2] = 0.0;


			/*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/

			if (implicit) {
				su2double a2;
				su2double phi;
				su2double Energy_el=0.0;

				//a2 = library->config->GetSpecies_Gamma(iSpecies)-1.0;
				//phi = a2*(library->0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

				for (iVar = 0; iVar < nVar; iVar++) {
					Jacobian_i[0][iVar] = 0.0;
					Jacobian_i[nDim+1][iVar] = 0.0;
					Jacobian_i[nDim+2][iVar] = 0.0;
        }

				for (iDim = 0; iDim < nDim; iDim++) {
					Jacobian_i[iDim+1][0] = -phi*Normal[iDim];

          for (unsigned short jDim = 0; jDim < nDim; jDim++)
					//	Jacobian_i[iDim+1][jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];

            Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
						Jacobian_i[iDim+1][nDim+2] = a2*Normal[iDim];
				}

			Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}
		}
	}
}

//
//
/*--- Far-field_BC ---*/
//
//

void CReactiveEulerSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) {

  unsigned long iVertex, iPoint;
	unsigned short iVar, iDim;
	RealVec U_domain(nVar), U_infty(nVar);

	bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

	/*--- Bucle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Interpolated solution at the wall ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				U_domain[iVar] = node[iPoint]->GetSolution(iVar);

			/*--- Solution at the infinity ---*/
			/*U_infty[0] = GetDensity_Inf(iSpecies);

			for (iDim = 0; iDim < nDim; iDim++)
				U_infty[iDim+1] = GetDensity_Velocity_Inf(iDim, iSpecies);

        U_infty[nDim+1] = GetDensity_Energy_Inf(iSpecies);*/

      /*--- Set the normal vector ---*/
			geometry->vertex[val_marker][iVertex]->GetNormal(Vector);
			for (iDim = 0; iDim < nDim; iDim++)
        Vector[iDim] = -Vector[iDim];

      conv_numerics->SetNormal(Vector);

			/*--- Compute the residual using an upwind scheme ---*/
			conv_numerics->SetConservative(U_domain.data(), U_infty.data());
			conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit)
				Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
		}
	}
}

//
//
/*--- Symmetry Plane BC ---*/
//
//

void CReactiveEulerSolver::BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container,
                                        CNumerics* conv_numerics, CNumerics* visc_numerics,
                                        CConfig* config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
	unsigned short iDim, iVar;
	su2double Pressure, Energy_el;
  SmartArr Normal;

	bool implicit = config->GetKind_TimeIntScheme() == EULER_IMPLICIT;

	/*--- Buckle over all the vertices ---*/
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- If the node belong to the domain ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/*--- Compute the projected residual ---*/
			Normal = Common::wrap_in_unique(geometry->vertex[val_marker][iVertex]->GetNormal());
      su2double Area = ::ComputeArea(Normal,nDim);

      RealVec UnitaryNormal(3);

			for (iDim = 0; iDim < nDim; iDim++)
				UnitaryNormal[iDim] = -Normal[iDim]/Area;

			//Residual[0] = 0.0;
			//Residual[nDim+1] = 0.0;
			//Residual[nDim+2] = 0.0;

      //Pressure = node[iPoint]->GetPressure(iSpecies);
			for (iDim = 0; iDim < nDim; iDim++)
				Residual[iDim+1] = Pressure*UnitaryNormal[iDim]*Area;


      /*--- Add value to the residual ---*/
			LinSysRes.AddBlock(iPoint, Residual);

			/*--- In case we are doing a implicit computation ---*/
			if (implicit) {
				su2double a2 = 0.0;
				su2double phi = 0.0;
				Energy_el = 0.0;

        //a2 = config->GetSpecies_Gamma(iSpecies)-1.0;
				//phi = a2*(0.5*node[iPoint]->GetVelocity2(iSpecies) - config->GetEnthalpy_Formation(iSpecies) - Energy_el);

				for (iVar =  0; iVar <  nVar; iVar++) {
					Jacobian_i[0][iVar] = 0.0;
					Jacobian_i[nDim+1][iVar] = 0.0;
					Jacobian_i[nDim+2][iVar] = 0.0;
				}
				for (iDim = 0; iDim < nDim; iDim++) {
					Jacobian_i[iDim+1][0] = -phi*Normal[iDim];
					//for (unsigned short jDim = 0; jDim < nDim; jDim++)
					//	Jacobian_i[iDim+1][loc + jDim+1] = a2*node[iPoint]->GetVelocity(jDim,iSpecies)*Normal[iDim];

          Jacobian_i[iDim+1][nDim+1] = -a2*Normal[iDim];
					Jacobian_i[iDim+1][nDim+2] = a2*Normal[iDim];
				}

				Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
			}

		}
	}
}
















































































































































































































































































































//
//
/*--- Preprocessing variables for NS ---*/
//
//

void CReactiveNSSolver::Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                      unsigned short iMesh, unsigned short iRKStep,
                                      unsigned short RunTime_EqSystem, bool Output) {

unsigned long ErrorCounter = 0;

unsigned long ExtIter = config->GetExtIter();
bool implicit         = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
bool second_order     = ((config->GetSpatialOrder_Flow() == SECOND_ORDER or config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER)
                         and iMesh == MESH_0);
bool limiter          = (config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER and ExtIter <= config->GetLimiterIter());
bool center           = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
bool limiter_visc     = config->GetViscous_Limiter_Flow();

/*--- Set the primitive variables ---*/

ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);


/*--- Gradient computation ---*/

if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
}
if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
}


/*--- Limiter computation ---*/

if ((limiter || limiter_visc) && (iMesh == MESH_0) && !Output)
  SetPrimitive_Limiter(geometry, config);


/*--- Artificial dissipation ---*/

if (center && !Output)
    throw Common::NotImplemented("Centered convective scheme not implemented\n");


/*--- Initialize the Jacobian matrices ---*/

if (implicit && !config->GetDiscrete_Adjoint())
  Jacobian.SetValZero();

}

//
//
/*--- Setting primitive variables NS ---*/
//
//

unsigned long CReactiveNSSolver::SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) {

  unsigned long iPoint, ErrorCounter = 0;
  bool RightSol = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);

    /*--- Compressible flow, primitive variables nSpecies+nDim+8, (T, vx, vy, vz, P, rho, h, c, lamMu, ThCond, Cp,rhos) ---*/

    RightSol = static_cast<CReactiveNSVariable*>(node[iPoint])->SetPrimVar();

    if (!RightSol) {
       node[iPoint]->SetNon_Physical(true);
       ErrorCounter++;
    }

    /*--- Initialize the convective, source and viscous residual vector ---*/

    if (!Output) LinSysRes.SetBlock_Zero(iPoint);

  }

  return ErrorCounter;

}

//
//
/*--- Setting time step NS ---*/
//
//

void CReactiveNSSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                                        unsigned short iMesh, unsigned long Iteration) {

su2double Area, Vol, Mean_SoundSpeed = 0.0, Mean_ProjVel = 0.0, Lambda, Local_Delta_Time,
          Lambda_1,Global_Delta_Time = 1E6,Mean_LaminarVisc = 0.0,Mean_Density = 0.0;
SmartArr Normal;
unsigned long iEdge, iVertex, iPoint, jPoint;
unsigned short iDim, iMarker;

bool implicit = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
bool time_steping = config->GetUnsteady_Simulation() == TIME_STEPPING;

Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
  node[iPoint]->SetMax_Lambda_Inv(0.0);
  node[iPoint]->SetMax_Lambda_Visc(0.0);
}

/*--- Loop interior edges ---*/
for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  /*--- Point identification, Normal vector and area ---*/

  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);

  Normal = Common::wrap_in_unique(geometry->edge[iEdge]->GetNormal());
  Area = ::ComputeArea(Normal,nDim);

  /*--- Mean Values ---*/

  Mean_ProjVel = 0.5 * (node[iPoint]->GetProjVel(Normal.get()) + node[jPoint]->GetProjVel(Normal.get()));
  Mean_SoundSpeed = 0.5 * (node[iPoint]->GetSoundSpeed() + node[jPoint]->GetSoundSpeed()) * Area;


  /*--- Inviscid contribution ---*/
  Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;

  if (geometry->node[iPoint]->GetDomain())
    node[iPoint]->AddMax_Lambda_Inv(Lambda);
  if (geometry->node[jPoint]->GetDomain())
    node[jPoint]->AddMax_Lambda_Inv(Lambda);

  /*--- Viscous contribution ---*/

  Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
  Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));

  Lambda_1 = (4.0/3.0)*Mean_LaminarVisc;
  Lambda = Lambda_1*Area*Area/Mean_Density;

  if (geometry->node[iPoint]->GetDomain())
    node[iPoint]->AddMax_Lambda_Visc(Lambda);
  if (geometry->node[jPoint]->GetDomain())
    node[jPoint]->AddMax_Lambda_Visc(Lambda);


}

/*--- Loop boundary edges ---*/

for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
  if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
  for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

    /*--- Point identification, Normal vector and area ---*/

    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    Normal = Common::wrap_in_unique(geometry->vertex[iMarker][iVertex]->GetNormal());
    Area = ::ComputeArea(Normal,nDim);

    /*--- Mean Values ---*/

    Mean_ProjVel = node[iPoint]->GetProjVel(Normal.get());
    Mean_SoundSpeed = node[iPoint]->GetSoundSpeed() * Area;

    /*--- Inviscid contribution ---*/
    Lambda = std::abs(Mean_ProjVel) + Mean_SoundSpeed;
    if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Inv(Lambda);

    /*--- Viscous contribution ---*/

    Mean_LaminarVisc = 0.5*(node[iPoint]->GetLaminarViscosity() + node[jPoint]->GetLaminarViscosity());
    Mean_Density     = 0.5*(node[iPoint]->GetSolution(0) + node[jPoint]->GetSolution(0));

    Lambda_1 = (4.0/3.0)*Mean_LaminarVisc;
    Lambda = Lambda_1 *Area*Area/Mean_Density;

    if (geometry->node[iPoint]->GetDomain())
      node[iPoint]->AddMax_Lambda_Visc(Lambda);

  }
}

/*--- Each element uses their own speed, steady state simulation ---*/

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

  Vol = geometry->node[iPoint]->GetVolume();

  if (Vol != 0.0) {
    Local_Delta_Time = config->GetCFL(iMesh)*Vol / node[iPoint]->GetMax_Lambda_Inv();
    Global_Delta_Time = std::min(Global_Delta_Time, Local_Delta_Time);
    Min_Delta_Time = std::min(Min_Delta_Time, Local_Delta_Time);
    Max_Delta_Time = std::max(Max_Delta_Time, Local_Delta_Time);
    if (Local_Delta_Time > config->GetMax_DeltaTime())
      Local_Delta_Time = config->GetMax_DeltaTime();
      node[iPoint]->SetDelta_Time(Local_Delta_Time);
  }
  else {
    node[iPoint]->SetDelta_Time(0.0);
  }
}

for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    if(time_steping) {
    /*--- Sets the regular CFL equal to the unsteady CFL ---*/
    config->SetCFL(iMesh,config->GetUnst_CFL());

    /*--- If the unsteady CFL is set to zero, it uses the defined unsteady time step, otherwise
    it computes the time step based on the unsteady CFL ---*/
    if (config->GetCFL(iMesh) == 0.0)
      node[iPoint]->SetDelta_Time(config->GetDelta_UnstTime());
    else
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
    }

  }
}

//
//
/*--- Viscous Residual ---*/
//
//

void CReactiveNSSolver::Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) {

unsigned long iPoint, jPoint, iEdge;

bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  /*--- Points, coordinates and normal vector in edge ---*/

  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);
  numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
  numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

  /*--- Primitive and secondary variables ---*/

  numerics->SetPrimitive(node[iPoint]->GetPrimitive(), node[jPoint]->GetPrimitive());
  numerics->SetSecondary(node[iPoint]->GetSecondary(), node[jPoint]->GetSecondary());

  /*--- Gradient and limiters ---*/

  numerics->SetPrimVarGradient(node[iPoint]->GetGradient_Primitive(), node[jPoint]->GetGradient_Primitive());
  numerics->SetPrimVarLimiter(node[iPoint]->GetLimiter_Primitive(), node[jPoint]->GetLimiter_Primitive());

  /*--- Compute and update residual ---*/

  numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

  LinSysRes.SubtractBlock(iPoint, Res_Visc);
  LinSysRes.AddBlock(jPoint, Res_Visc);

  /*--- Implicit part ---*/

  if (implicit) {
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
  }

}

}
