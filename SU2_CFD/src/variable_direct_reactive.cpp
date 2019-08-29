#include "../include/variable_reactive.hpp"

/*--- Declaration of static variables ---*/
unsigned short CReactiveEulerVariable::P_INDEX_PRIM = CReactiveEulerVariable::VX_INDEX_PRIM + CReactiveEulerVariable::nDim;
unsigned short CReactiveEulerVariable::RHO_INDEX_PRIM = CReactiveEulerVariable::P_INDEX_PRIM + 1;
unsigned short CReactiveEulerVariable::H_INDEX_PRIM = CReactiveEulerVariable::RHO_INDEX_PRIM + 1;
unsigned short CReactiveEulerVariable::A_INDEX_PRIM = CReactiveEulerVariable::H_INDEX_PRIM + 1;
unsigned short CReactiveEulerVariable::RHOS_INDEX_PRIM = CReactiveEulerVariable::A_INDEX_PRIM + 1;

unsigned short CReactiveEulerVariable::RHOE_INDEX_SOL = CReactiveEulerVariable::RHOVX_INDEX_SOL + CReactiveEulerVariable::nDim;
unsigned short CReactiveEulerVariable::RHOS_INDEX_SOL = CReactiveEulerVariable::RHOE_INDEX_SOL + 1;

unsigned short CReactiveEulerVariable::P_INDEX_GRAD = CReactiveEulerVariable::VX_INDEX_GRAD + CReactiveEulerVariable::nDim;

unsigned short CReactiveEulerVariable::P_INDEX_LIM = CReactiveEulerVariable::VX_INDEX_LIM + CReactiveEulerVariable::nDim;

unsigned short CReactiveNSVariable::RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

//
//
/*--- Class default constructor ---*/
//
//
CReactiveEulerVariable::CReactiveEulerVariable(): CVariable(), nSpecies(), nPrimVarLim(), US_System(false), Cp() {
  /*--- Set to zero some recurrent variables of elementary type ---*/
  nPrimVar = 0;
  nPrimVarGrad = 0;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
}

//
//
/*--- Class constructor to initialize dimensions of the problem ---*/
//
//
CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                               unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                               LibraryPtr lib_ptr, CConfig* config): CVariable(val_nDim, val_nvar, config),
                                                                                     library(lib_ptr), Cp() {
  /*--- Update indexes ---*/
  P_INDEX_PRIM = CReactiveEulerVariable::VX_INDEX_PRIM + CReactiveEulerVariable::nDim;
  RHO_INDEX_PRIM = P_INDEX_PRIM + 1;
  H_INDEX_PRIM = RHO_INDEX_PRIM + 1;
  A_INDEX_PRIM = H_INDEX_PRIM + 1;
  RHOS_INDEX_PRIM = A_INDEX_PRIM + 1;

  RHOE_INDEX_SOL = RHOVX_INDEX_SOL + CReactiveEulerVariable::nDim;
  RHOS_INDEX_SOL = RHOE_INDEX_SOL + 1;

  P_INDEX_GRAD = CReactiveEulerVariable::VX_INDEX_GRAD + CReactiveEulerVariable::nDim;

  P_INDEX_LIM = CReactiveEulerVariable::VX_INDEX_LIM + CReactiveEulerVariable::nDim;

  /*--- Set variables ---*/
  nSpecies = val_nSpecies;
  nPrimVar = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nPrimVarLim = val_nprimvarlim;

  US_System = (config->GetSystemMeasurements() == US);

  /*--- Allocate array related to solution ---*/
  Solution_Max = new su2double[nPrimVarLim];
  Solution_Min = new su2double[nPrimVarLim];

  for(unsigned short iVar = 0; iVar < nPrimVarLim; ++iVar) {
    Solution_Min[iVar] = 0.0;
    Solution_Max[iVar] = 0.0;
  }

  /*--- Allocate residual structures ---*/
  Res_TruncError = new su2double[nVar];
  std::fill(Res_TruncError, Res_TruncError + nVar, 0.0);

  /*--- Only for residual smoothing (multigrid) ---*/
  unsigned short nMGSmooth = 0;
  for(unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); ++iMesh)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if(nMGSmooth > 0) {
    Residual_Sum = new su2double[nVar];
    Residual_Old = new su2double[nVar];
  }

  /*--- Allocate and initializate solution for dual time strategy ---*/
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if (dual_time) {
    Solution_time_n = new su2double[nVar];
    Solution_time_n1 = new su2double[nVar];
  }

  /*--- Allocate primitive vectors ---*/
  Primitive.resize(nPrimVar);
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for(unsigned short iVar = 0; iVar < nPrimVarGrad; ++iVar) {
    Gradient_Primitive[iVar] = new su2double[nDim];
    std::fill(Gradient_Primitive[iVar], Gradient_Primitive[iVar] + nDim, 0.0);
  }
  Limiter_Primitive.resize(nPrimVarLim, 0.0);

  dTdU.resize(nVar);
  dPdU.resize(nVar);

  Ys.resize(nSpecies);
}

//
//
/*--- Class overloaded constructor (pressuure, temprature,, mass fractions and velocity initialization values) ---*/
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                               const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                               CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad,
                                                                      val_nprimvarlim, lib_ptr, config) {
  /*--- Check Solution and Solution_Old arrays allocation ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  /*--- Rename and initialize for convenience ---*/
  su2double T = val_temperature;   // Translational-rotational temperature
  su2double P = val_pressure;      // Pressure

  su2double rho, rhoE;

  /*--- Compute mixture density ---*/
  rho = library->ComputeDensity(T, P, val_massfrac)*config->GetGas_Constant_Ref();

  /*--- Compute energy (RHOE) from supplied primitive quanitites ---*/
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  su2double sqvel = std::inner_product(val_velocity.cbegin(), val_velocity.cend(), val_velocity.cbegin(), 0.0);
  su2double e_tot = library->ComputeEnergy(dim_temp, val_massfrac)/config->GetEnergy_Ref();
  if(US_System)
    e_tot *= 3.28084*3.28084;
  rhoE = rho*(e_tot + 0.5*sqvel);

  /*--- Initialize Solution and Solution_Old vectors ---*/
  /*--- Initialize mixture density and partial density ---*/
  Solution[RHO_INDEX_SOL] = Solution_Old[RHO_INDEX_SOL] = rho;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Solution[RHOS_INDEX_SOL + iSpecies] = Solution_Old[RHOS_INDEX_SOL + iSpecies] = rho*val_massfrac[iSpecies];

  /*--- Initialize momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    Solution[RHOVX_INDEX_SOL + iDim] = Solution_Old[RHOVX_INDEX_SOL + iDim] = rho*val_velocity[iDim];

  /*--- Initialize energy contribution ---*/
  Solution[RHOE_INDEX_SOL] = Solution_Old[RHOE_INDEX_SOL] = rhoE;

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if(dual_time) {
    for(unsigned short iVar = 0; iVar < nVar; ++iVar) {
      Solution_time_n[iVar] = Solution[iVar];
      Solution_time_n1[iVar] = Solution[iVar];
    }
  }

  /*--- Assign primitive variables ---*/
  Primitive.at(T_INDEX_PRIM) = T;
  Primitive.at(P_INDEX_PRIM) = P;
}

//
//
/*--- Class overloaded constructor (initialization vector with a complete state). ---*/
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                               CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad,
                                                                      val_nprimvarlim, lib_ptr, config) {
  /*--- Initialize Solution and Solution_Old array ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  unsigned short iVar;
  for(iVar = 0; iVar < nVar; ++iVar) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if (dual_time) {
    for(iVar = 0; iVar < nVar; ++iVar) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }

  /*--- Initialize T to the free stream for the secant method ---*/
  Primitive.at(T_INDEX_PRIM) = config->GetTemperature_FreeStream()/config->GetTemperature_Ref();
}

//
//
/*--- Class overloaded constructor (initialization array with a complete state). ---*/
//
//
CReactiveEulerVariable::CReactiveEulerVariable(su2double* val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                               CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad,
                                                                      val_nprimvarlim, lib_ptr, config) {
  /*--- Initialize Solution and Solution_Old arrays ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");
  SU2_Assert(val_solution != NULL,"The array to initialize Solution and Solution_Old has not been allocated");

  unsigned short iVar;
  for(iVar = 0; iVar < nVar; ++iVar) {
    Solution[iVar] = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  if (dual_time) {
    for(iVar = 0; iVar < nVar; ++iVar) {
      Solution_time_n[iVar] = val_solution[iVar];
      Solution_time_n1[iVar] = val_solution[iVar];
    }
  }

  /*--- Initialize T to the free stream for the secant method ---*/
  Primitive.at(T_INDEX_PRIM) = config->GetTemperature_FreeStream()/config->GetTemperature_Ref();
}


//
//
/*--- Class destructor ---*/
//
//
CReactiveEulerVariable::~CReactiveEulerVariable() {
  if(Gradient_Primitive != NULL) {
    for(unsigned short iVar = 0; iVar < nPrimVarGrad; ++iVar)
      delete[] Gradient_Primitive[iVar];
    delete[] Gradient_Primitive;
  }
}

//
//
/*--- Set primitive variables ---*/
//
//
bool CReactiveEulerVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables ---*/
  bool nonPhys = Cons2PrimVar(config, Solution, Primitive.data());

  /*--- Check for non physical solutions. NOTE: The first global iteration is taken into account by the solver class ---*/
  if(nonPhys && config->GetExtIter() > 0) {
    std::copy(Solution_Old, Solution_Old + nVar, Solution);
    bool nonPhys_old = Cons2PrimVar(config, Solution, Primitive.data());
    SU2_Assert(nonPhys_old == false, "Neither the old solution is feasible to set primitive variables: problem unsolvable");
  }

  /*--- Set specific heat at constant pressure ---*/
  su2double dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  su2double dim_a = Primitive[A_INDEX_PRIM]*config->GetVelocity_Ref();
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_a /= 3.28084;
  }
  Cp = library->ComputeCP_FromSoundSpeed(dim_temp, dim_a, Ys)/config->GetGas_Constant_Ref();
  if(US_System)
    Cp *= 3.28084*3.28084*5.0/9.0;

  /*--- Set temperature and pressure derivatives ---*/
  CalcdTdU(Primitive.data(), config, dTdU.data());
  CalcdPdU(Primitive.data(), config, dPdU.data());

  return nonPhys;
}

//
//
/*--- Pass from conserved to primitive variables ---*/
//
//
bool CReactiveEulerVariable::Cons2PrimVar(CConfig* config, su2double* U, su2double* V) {
  if(config->GetExtIter() == 0) {
    SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
    SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");
  }

  bool NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter;
  su2double rho, rhoE;
  su2double sqvel;
  su2double f, df, NRtol, Btol;
  su2double T, Told, Tnew, hs, hs_old, Tmin, Tmax;

  /*--- Conserved and primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoE, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, Y1, ..., YNs]^T

  /*--- Set boolean for non physical state---*/
  nonPhys = false;

  /*--- Assign species mass fraction and mixture density ---*/
  // NOTE: If any species densities are < 0, these values are re-assigned
  //       in the conserved vector to ensure positive density
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < 0.0) {
      U[RHOS_INDEX_SOL + iSpecies] = 1.0e-30;
      nonPhys = true;
    }
  }

  if(U[RHO_INDEX_SOL] < EPS) {
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL] = EPS;
    nonPhys = true;
  }
  else
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL];

  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL + iSpecies]/U[RHO_INDEX_SOL];

  /*--- Checking sum of mass fraction ---*/
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  nonPhys = nonPhys || (std::abs(std::accumulate(Ys.cbegin(), Ys.cend(), 0.0) - 1.0) > 0.1);

  /*--- Rename for convenience ---*/
  rho = U[RHO_INDEX_SOL];    // Density [Kg/m3]
  rhoE = U[RHOE_INDEX_SOL];   // Density*total energy per unit of mass [J/m3]

  /*--- Assign mixture velocity and compute squared velocity ---*/
  sqvel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim) {
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
    sqvel += V[VX_INDEX_PRIM + iDim]*V[VX_INDEX_PRIM + iDim];
  }

  /*--- Set temperature clipping values ---*/
  Tmin   = config->GetTemperatureMin()/config->GetTemperature_Ref();
  Tmax   = config->GetTemperatureMax()/config->GetTemperature_Ref();

  /*--- Set temperature secant algorithm paramters ---*/
  NRtol    = 1.0e-6;    // Tolerance for the Secant method
  Btol     = 1.0e-4;    // Tolerance for the Bisection method
  maxNIter = 7;        // Maximum Secant method iterations
  maxBIter = 32;        // Maximum Bisection method iterations
  NRconvg  = false;

  /*--- Translational-Rotational Temperature ---*/
  const su2double Rgas = library->ComputeRgas(Ys)/config->GetGas_Constant_Ref();
  const su2double C1 = (-rhoE + 0.5*rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  /*--- Pick initial state and start algorithm ---*/
  const su2double old_temp = V[T_INDEX_PRIM];
  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    try {
      /*--- Execute a secant root-finding method to find T ---*/
      su2double dim_temp, dim_temp_old;
      dim_temp = T*config->GetTemperature_Ref();
      dim_temp_old = Told*config->GetTemperature_Ref();
      if(US_System) {
        dim_temp *= 5.0/9.0;
        dim_temp_old *= 5.0/9.0;
      }
      hs_old = library->ComputeEnthalpy(dim_temp_old, Ys)/config->GetEnergy_Ref();
      hs = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
      if(US_System) {
        hs_old *= 3.28084*3.28084;
        hs *= 3.28084*3.28084;
      }

      f = T - C1 - C2*hs;
      df = T - Told + C2*(hs_old-hs);
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
    catch(const std::out_of_range& e) {
      /*--- Print error message ---*/
      std::cout<<e.what()<<" Trying with bisection method."<<std::endl;

      /*--- Execute bisection method ---*/
      su2double Ta = Tmin;
      su2double Tb = Tmax;
      for(iIter = 0; iIter < 10000; ++iIter) {
        T = (Ta + Tb)/2.0;
        su2double dim_temp = T*config->GetTemperature_Ref();;
        if(US_System)
          dim_temp *= 5.0/9.0;
        hs = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
        if(US_System)
          hs *= 3.28084*3.28084;
        f = T - C1 - C2*hs;

        if(std::abs(f) < Btol) {
          NRconvg = true;
          break;
        }
        else {
          if(f > 0)
            Ta = T;
          else
            Tb = T;
       }
     }
     if(NRconvg)
        break;
     else
        throw std::runtime_error("Convergence not achieved for bisection method after catching out of range");
    }
  }

  /*--- If the secant method has converged, assign the value of T.
        Otherwise, if no exception has been caught, execute a bisection root-finding method ---*/
  if(NRconvg)
    V[T_INDEX_PRIM] = T;
  else {
    Bconvg = false;
    su2double Ta = Tmin;
    su2double Tb = Tmax;
    for(iIter = 0; iIter < maxBIter; ++iIter) {
      T = (Ta + Tb)/2.0;
      su2double dim_temp = T*config->GetTemperature_Ref();;
      if(US_System)
        dim_temp *= 5.0/9.0;
      hs = library->ComputeEnthalpy(dim_temp, Ys)/config->GetEnergy_Ref();
      if(US_System)
        hs *= 3.28084*3.28084;
      f = T - C1 - C2*hs;

      if(std::abs(f) < Btol) {
        V[T_INDEX_PRIM] = T;
        Bconvg = true;
        break;
      }
      else {
        if(f > 0)
          Ta = T;
        else
          Tb = T;
      }
    }

    /*--- If absolutely no convergence, then something is going really wrong ---*/
    if(!Bconvg)
      throw std::runtime_error("Convergence not achieved for bisection method");
  }

  /*--- Avoid too large variation in temperature ---*/
  if(config->GetExtIter() > 0 && config->GetClipping_Temp())
     V[T_INDEX_PRIM] = std::min(std::max(V[T_INDEX_PRIM],0.95*old_temp),1.05*old_temp);

  /*--- Check if the solution found is inside the limits, otherwise set non physical point ---*/
  if(V[T_INDEX_PRIM] < Tmin) {
    V[T_INDEX_PRIM] = Tmin;
    nonPhys = true;
  }
  else if(V[T_INDEX_PRIM] > Tmax) {
    V[T_INDEX_PRIM] = Tmax;
    nonPhys = true;
  }

  /*--- Rename for convenience ---*/
  T = V[T_INDEX_PRIM];

  /*--- Pressure ---*/
  V[P_INDEX_PRIM] = rho*Rgas*T;
  if(V[P_INDEX_PRIM] < EPS) {
    V[P_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Sound speed ---*/
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  V[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp, Ys, V[P_INDEX_PRIM], rho);
  if(V[A_INDEX_PRIM] < EPS) {
    V[A_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Enthalpy ---*/
  V[H_INDEX_PRIM] = (U[RHOE_INDEX_SOL] + V[P_INDEX_PRIM])/rho;

  return nonPhys;
}

//
//
/*--- Compute temperature derivatives ---*/
//
//
void CReactiveEulerVariable::CalcdTdU(su2double* V, CConfig* config, su2double* dTdU) {
  /*--- Check memory allocation---*/
  if(config->GetExtIter() == 0)
    SU2_Assert(dTdU != NULL,"The array to store the derivatives of temperature w.r.t. conserved variables has not been allocated");

  /*--- Rename for convenience ---*/
  su2double rho = V[RHO_INDEX_PRIM];
  su2double T = V[T_INDEX_PRIM];

  /*--- Compute useful quantities ---*/
  su2double dim_cp = Cp*config->GetGas_Constant_Ref();
  if(US_System)
    dim_cp /= 3.28084*3.28084*5.0/9.0;
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM+ nSpecies), Ys.begin());
  su2double Cv = library->ComputeCV_FromCP(dim_cp, Ys)/config->GetGas_Constant_Ref();
  if(US_System)
    Cv *= 3.28084*3.28084*5.0/9.0;
  su2double rhoCv = rho*Cv;
  su2double sq_vel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;

  dTdYs = library->ComputePartialEnergy(dim_temp);
  for(auto& elem: dTdYs)
    elem /= config->GetEnergy_Ref();
  if(US_System) {
    for(auto& elem: dTdYs)
      elem *= 3.28084*3.28084;
  }

  /*--- Set temperature derivatives ---*/
  dTdU[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    dTdU[RHOVX_INDEX_SOL + iDim] = -V[VX_INDEX_PRIM + iDim]/rhoCv;
  dTdU[RHOE_INDEX_SOL] = 1.0/rhoCv;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    dTdU[RHOS_INDEX_SOL + iSpecies] = -dTdYs[iSpecies]/rhoCv;
}

//
//
/*--- Compute pressure derivatives ---*/
//
//
void CReactiveEulerVariable::CalcdPdU(su2double* V, CConfig* config, su2double* dPdU) {
  /*--- Check memory allocation---*/
  if(config->GetExtIter() == 0)
    SU2_Assert(dPdU != NULL,"The array to store the derivatives of pressure w.r.t. conserved variables has not been allocated");

  /*--- Useful quantities ---*/
  su2double dim_cp = Cp*config->GetGas_Constant_Ref();
  if(US_System)
    dim_cp /= 3.28084*3.28084*5.0/9.0;
  su2double Gamma = library->ComputeFrozenGamma_FromCP(dim_cp, Ys);
  su2double sq_vel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Rename for convenience ---*/
  su2double T = V[T_INDEX_PRIM];

  /*--- Set pressure derivatives ---*/
  dPdU[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*sq_vel;
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    dPdU[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V[VX_INDEX_PRIM + iDim];
  dPdU[RHOE_INDEX_SOL] = Gamma - 1.0;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    dPdU[RHOS_INDEX_SOL + iSpecies] = library->GetRiGas(iSpecies)/config->GetGas_Constant_Ref()*T - (Gamma - 1.0)*dTdYs[iSpecies];
}


//
//
/*--- Pass from primitive to conserved variables ---*/
//
//
void CReactiveEulerVariable::Prim2ConsVar(CConfig* config, su2double* V, su2double* U) {
  if(config->GetExtIter() == 0) {
    SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
    SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");
  }

  /*--- Set mixture density and species density ---*/
  U[RHO_INDEX_SOL] = V[RHO_INDEX_PRIM];
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    U[RHOS_INDEX_SOL + iSpecies] = V[RHO_INDEX_PRIM]*V[RHOS_INDEX_PRIM + iSpecies];

  /*--- Set momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM + iDim];

  /*--- Set energies ---*/
  U[RHOE_INDEX_SOL] = V[RHO_INDEX_PRIM]*V[H_INDEX_PRIM] - V[P_INDEX_PRIM];
}

//
//
/*--- Set mixture density ---*/
//
//
inline bool CReactiveEulerVariable::SetDensity(void) {
  /*--- Check memory allocation for safety ---*/
  SU2_Assert(Solution != NULL,"The array of solution variables has not been allocated");

  Primitive.at(RHO_INDEX_PRIM) = Solution[RHO_INDEX_SOL];

  if(Primitive.at(RHO_INDEX_PRIM) < EPS)
    return true;

  return false;
}

//
//
/*--- Set pressure (NOTE: it requires SetDensity() call) ---*/
//
//
bool CReactiveEulerVariable::SetPressure(CConfig* config) {
  /*--- Compute mixture pressure ---*/
  std::copy(Primitive.cbegin() + RHOS_INDEX_PRIM, Primitive.cbegin() + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  su2double Pressure = library->ComputePressure(Primitive.at(T_INDEX_PRIM), Primitive.at(RHO_INDEX_PRIM), Ys);
  Pressure /= config->GetGas_Constant_Ref();

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if(Pressure < EPS)
    return true;

  return false;
}

//
//
/*--- Set sound speed (NOTE: it requires SetTemperature() call) ---*/
//
//
bool CReactiveEulerVariable::SetSoundSpeed(CConfig* config) {
  /*--- Compute useful quantities ---*/
  su2double dim_cp = Cp*config->GetGas_Constant_Ref();
  su2double dim_temp = Primitive.at(T_INDEX_PRIM)*config->GetTemperature_Ref();
  if(US_System) {
    dim_cp /= 3.28084*3.28084*5.0/9.0;
    dim_temp *= 5.0/9.0;
  }
  std::copy(Primitive.cbegin() + RHOS_INDEX_PRIM, Primitive.cbegin() + (RHOS_INDEX_PRIM+ nSpecies), Ys.begin());

  /*--- Compute frozen sound spped ---*/
  su2double Sound_Speed = library->ComputeFrozenSoundSpeed_FromCP(dim_temp, dim_cp, Ys);
  Sound_Speed /= config->GetVelocity_Ref();
  if(US_System)
    Sound_Speed *= 3.28084;

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(A_INDEX_PRIM) = Sound_Speed;
  if(Sound_Speed < EPS)
    return true;

  return false;
}

//
//
/*--- Set the velocity vector from old solution ---*/
//
//
inline void CReactiveEulerVariable::SetVelocity_Old(su2double* val_velocity) {
  /*--- Check memory allocation for safety ---*/
  SU2_Assert(val_velocity != NULL,"The array of velocity val_velocity has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array of solution variables at previous step has not been allocated");

  for(unsigned short iDim=0; iDim < nDim; ++iDim)
    Solution_Old[RHOVX_INDEX_SOL+iDim] = val_velocity[iDim]*Primitive.at(RHO_INDEX_PRIM);
}

//
//
/*--- Class constructor to initialize dimensions of the problem. ---*/
//
//
CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                         unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                         LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar,
                                                                val_nprimvargrad, val_nprimvarlim, lib_ptr, config),
                                                                Laminar_Viscosity(), Thermal_Conductivity() {
  /*--- Update index ---*/
  RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

  /*--- Resize matrix ---*/
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}

//
//
/*--- Class overloaded constructor (initialization values for pressure, temperature, mass fractions, velocity, viscosity, conductivity). ---*/
//
//
CReactiveNSVariable::CReactiveNSVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                         const su2double val_temperature, const su2double val_viscosity, unsigned short val_nDim,
                                         unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                                         unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_pressure, val_massfrac, val_velocity, val_temperature, val_nDim,
                                                                val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad, val_nprimvarlim,
                                                                lib_ptr, config), Laminar_Viscosity(val_viscosity), Thermal_Conductivity() {
  /*--- Update index ---*/
  RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

  /*--- Resize matrix ---*/
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}

//
//
/*--- Class overloaded constructor (initialization vector with a complete state). ---*/
//
//
CReactiveNSVariable::CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_solution, val_nDim, val_nvar, val_nSpecies, val_nprimvar,
                                                                val_nprimvargrad, val_nprimvarlim, lib_ptr, config),
                                                                Laminar_Viscosity(), Thermal_Conductivity() {
  /*--- Update index ---*/
  RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

  /*--- Resize matrix ---*/
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}

//
//
/*--- Class overloaded constructor (initialization array with a complete state). ---*/
//
//
CReactiveNSVariable::CReactiveNSVariable(su2double* val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_solution, val_nDim, val_nvar, val_nSpecies, val_nprimvar,
                                                                val_nprimvargrad, val_nprimvarlim, lib_ptr, config),
                                                                Laminar_Viscosity(), Thermal_Conductivity() {
  /*--- Update index ---*/
  RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

  /*--- Resize matrix ---*/
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}


//
//
/*--- Set primitive variables ---*/
//
//
bool CReactiveNSVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables using Euler version since primitives are the same ---*/
  bool nonPhys = CReactiveEulerVariable::SetPrimVar(config);

  su2double dim_temp, dim_press;
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  /*--- Compute transport properties --- */
  Ys = GetMassFractions();
  Laminar_Viscosity = library->ComputeEta(dim_temp, Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  Thermal_Conductivity = library->ComputeLambda(dim_temp, Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  /*--- Compute binary diffusion coefficents. NOTE: The empirical formula employed in the library should return it in cm2/s ---*/
  Diffusion_Coeffs = library->GetDij_SM(dim_temp, dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1.0e4);
  if(US_System)
    Diffusion_Coeffs *= 3.28084*3.28084;

  return nonPhys;
}
