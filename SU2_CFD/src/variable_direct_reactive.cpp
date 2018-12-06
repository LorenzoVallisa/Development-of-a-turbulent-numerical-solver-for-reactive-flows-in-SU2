#include "../include/variable_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/move_pointer.hpp"
#include "../../Common/include/default_initialization.hpp"

//
//
/*!
 *\brief Class default constructor
 */
//
//
unsigned short CReactiveEulerVariable::nSpecies = 0;
CReactiveEulerVariable::LibraryPtr CReactiveEulerVariable::library = NULL;

CReactiveEulerVariable::CReactiveEulerVariable():CVariable() {
  std::tie(nPrimVar,nPrimVarGrad,nSecondaryVar,nSecondaryVarGrad,nPrimVarLim) = Common::repeat<5,decltype(nPrimVar)>(decltype(nPrimVar)());
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim,unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                                               unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config):
                        CVariable(val_nDim,val_nvar,config) {
  nSpecies = val_nSpecies;
  nPrimVar = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nPrimVarLim = val_nprimvarlim;

  library = LibraryPtr(new Framework::ReactingModelLibrary(config->GetConfigLibFile()));

  Primitive.resize(nPrimVar);
  Gradient_Primitive.resize(nPrimVarGrad,nDim);
  Limiter_Primitive.resize(nPrimVarLim);

  Limiter = new su2double [nVar];
  Solution_Max = new su2double[nPrimVarLim];
  Solution_Min = new su2double[nPrimVarLim];

}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                               const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, CConfig* config):
                                               CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                      val_nprimvarlim,config) {

  /*--- Rename and initialize for convenience ---*/
  su2double T = val_temperature;   // Translational-rotational temperature [K]
  su2double P = val_pressure;

  RealVec Yi = val_massfrac;
  su2double rho,rhoE;

  /*--- Compute mixture density ---*/
  //rho = library->ComputeDensity(T,P,Yi);
  rho *= config->GetGas_Constant_Ref();

  su2double dim_temp = T*config->GetTemperature_Ref();
  bool US_System = config->GetSystemMeasurements() ==US;
  if(US_System)
    dim_temp *= 5.0/9.0;
  //su2double Sound_Speed = library->ComputeFrozenSoundSpeed(dim_temp,Yi,P,rho);

  /*--- Calculate energy (RHOE) from supplied primitive quanitites ---*/
  /*
  su2double sqvel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    sqvel += val_mach[iDim]*Sound_Speed * val_mach[iDim]*Sound_speed;
  su2double e_tot = library->ComputeEnergy(dim_temp,Yi)/config->GetEnergy_Ref();
  rhoE = rho*(0.5*sqvel + e_tot);
  */

  /*--- Initialize Solution and Solution_Old vectors ---*/

  /*--- Initialize mixture density and partial density ---*/
  Solution[RHO_INDEX_SOL] = Solution_Old[RHO_INDEX_SOL] = rho;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Solution[RHOS_INDEX_SOL + iSpecies] = Solution_Old[RHOS_INDEX_SOL + iSpecies] = rho*val_massfrac[iSpecies];

  /*--- Initialize momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    Solution[RHOVX_INDEX_SOL + iDim] = Solution_Old[RHOVX_INDEX_SOL + iDim] = rho*val_velocity[iDim];

  /*--- Initialize energy contribution ---*/
  Solution[RHOE_INDEX_SOL] = Solution_Old[nSpecies+nDim] = rhoE;

  /*--- Assign primitive variables ---*/
  Primitive[T_INDEX_PRIM] = T;
  Primitive[P_INDEX_PRIM] = P;

}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, CConfig* config):
                                               CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                      val_nprimvarlim,config) {
  /*--- Initialize Solution and Solution_Old vectors ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  std::copy(val_solution.cbegin(),val_solution.cend(),Solution);
  std::copy(val_solution.cbegin(),val_solution.cend(),Solution_Old);

  /*--- Initialize T and P to the free stream for Secant method ---*/
  Primitive[T_INDEX_PRIM] = config->GetTemperature_FreeStream();
  Primitive[P_INDEX_PRIM] = config->GetPressure_FreeStream();

}

//
//
/*!
 *\brief Set primitive variables
 */
//
//
bool CReactiveEulerVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables ---*/
  bool nonPhys = Cons2PrimVar(config, Solution, Primitive.data());
  if(nonPhys) {
    std::copy(Solution_Old,Solution_Old + nVar,Solution);
    bool check_old = Cons2PrimVar(config,Solution,Primitive.data());
    SU2_Assert(check_old == true, "Neither the old solution is feasible to set primitive variables: problem unsolvable");
  }
  return nonPhys;

}

//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//
bool CReactiveEulerVariable::Cons2PrimVar(CConfig* config, su2double* U, su2double* V) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  bool NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter;
  su2double rho, rhoE;
  su2double sqvel;
  su2double f, df, NRtol, Btol;
  su2double T,Told,Tnew,hs,hs_old,Tmin,Tmax;

  /*--- Conserved and primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoE, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 8E4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0E-6;    // Tolerance for the Secant method
  Btol     = 1.0E-4;    // Tolerance for the Bisection method
  maxNIter = 5;        // Maximum Secant method iterations
  maxBIter = 32;        // Maximum Bisection method iterations

  /*--- Rename variables forconvenience ---*/
  rhoE  = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species and mixture density ---*/
  // Note: If any species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < EPS) {
      V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL + iSpecies] = EPS;
      nonPhys = true;
    }
    else
      V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL + iSpecies];
  }
  if(U[RHO_INDEX_SOL] < EPS) {
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL] = EPS;
    nonPhys = true;
  }
  else
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL];

  // Rename for convenience
  rho = U[RHO_INDEX_SOL];
  /*--- Assign mixture velocity ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
  sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Translational-Rotational Temperature ---*/
  RealVec Yi(U + RHOS_INDEX_SOL,U + (RHOS_INDEX_SOL + nSpecies));
  std::for_each(Yi.begin(),Yi.end(),[rho](su2double& elem){elem /= rho;});
  const su2double Rgas = library->ComputeRgas(Yi)/config->GetGas_Constant_Ref();
  const su2double C1 = (-rhoE + rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  bool US_System = config->GetSystemMeasurements() == US;
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a secant root-finding method to find T ---*/
    su2double dim_temp, dim_temp_old;
    dim_temp = T*config->GetTemperature_Ref();
    dim_temp_old = Told*config->GetTemperature_Ref();
    if(US_System) {
      dim_temp *= 5.0/9.0;
      dim_temp_old *= 5.0/9.0;
    }
    //hs_old = library->ComputeEnthalpy(dim_temp_old,Yi)/config->GetEnergy_Ref();
    //hs = library->ComputeEnthalpy(dim_temp,Yi)/config->GetEnergy_Ref();
    if(US_System) {
      hs_old *= 3.28084*3.28084;
      hs *= 3.28084*3.28084;
    }

    f = T - C1 - C2*hs;
    df = T - Told + C2*(hs_old-hs);
    Tnew = T - f*(T-Told)/df;

    // Check for convergence
    if(std::abs(Tnew - T) < NRtol) {
      NRconvg = true;
      break;
    }
    else {
      Told = T;
      T = Tnew;
    }
  }

  // If the secant method has converged, assign the value of T.
  // Otherwise, execute a bisection root-finding method
  if(NRconvg)
    V[T_INDEX_PRIM] = T;
  else {
    // Execute the bisection root-finding method
    Bconvg = false;
    su2double Ta = Tmin;
    su2double Tb = Tmax;
    for(iIter = 0; iIter < maxBIter; ++iIter) {
      T = (Ta + Tb)/2.0;
      su2double dim_temp = T*config->GetTemperature_Ref();;
      if(US_System)
        dim_temp *= 5.0/9.0;
      //hs = library->ComputeEnthalpy(dim_temp,Yi)/config->GetEnergy_Ref();
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

    // If absolutely no convergence, then something is going really wrong
    if(!Bconvg)
      throw std::runtime_error("Convergence not achieved for bisection method");
    }
  }

  if(V[T_INDEX_PRIM] < Tmin) {
    V[T_INDEX_PRIM] = Tmin;
    nonPhys = true;
  }
  else if(V[T_INDEX_PRIM] > Tmax) {
    V[T_INDEX_PRIM] = Tmax;
    nonPhys = true;
  }

  /*--- Pressure ---*/
  //V[P_INDEX_PRIM] = library->ComputePressure(V[RHO_INDEX_PRIM],V[T_INDEX_PRIM],Yi)/config->GetGas_Constant_Ref();
  if(V[P_INDEX_PRIM] < EPS) {
    V[P_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Sound speed ---*/
  su2double dim_temp = V[T_INDEX_PRIM]*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //V[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp,Yi,V[P_INDEX_PRIM],V[RHO_INDEX_PRIM]);
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
/*!
 *\brief Pass from primitive to conserved variables
 */
//
//
void CReactiveEulerVariable::Prim2ConsVar(CConfig* config, su2double* V, su2double* U) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  unsigned short iDim, iSpecies;

  /*--- Set mixture density and species density ---*/
  U[RHO_INDEX_SOL] = V[RHO_INDEX_PRIM];
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies), U + RHOS_INDEX_SOL);

  /*--- Set momentum ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM + iDim];

  /*--- Set energies ---*/
  U[RHOE_INDEX_SOL] = V[RHO_INDEX_PRIM]*V[H_INDEX_PRIM] - V[P_INDEX_PRIM];

}

//
//
/*!
 *\brief Get Gradient Primitive Variable
 */
//
//
su2double** CReactiveEulerVariable::GetGradient_Primitive(void) {
  std::vector<std::vector<double>> tmp_matrix(Gradient_Primitive.nbRows(),std::vector<double>(Gradient_Primitive.nbCols()));
  for(std::size_t i = 0; i < Gradient_Primitive.nbRows(); ++i)
    for(std::size_t j = 0; j < Gradient_Primitive.nbCols(); ++j)
      tmp_matrix[i][j] = Gradient_Primitive(i,j);

  std::vector<su2double*> tmp;
  for(auto&& elem:tmp_matrix)
    tmp.push_back(elem.data());

  return tmp.data();
}


//
//
/*!
 *\brief Set density
 */
//
//
bool CReactiveEulerVariable::SetDensity(void) {
  SU2_Assert(Solution != NULL,"The array of solution variables has not been allocated");
  Primitive.at(RHO_INDEX_PRIM) = Solution[RHO_INDEX_SOL];

  if(Primitive.at(RHO_INDEX_PRIM) < EPS)
    return true;
  return false;
}

//
//
/*!
 *\brief Set pressure (requires SetDensity() call)
 *///
//
bool CReactiveEulerVariable::SetPressure(CConfig* config) {
  /*--- Compute this gas mixture property from library ---*/
  su2double Pressure;
  //su2double Pressure = library->ComputePressure(Primitive.at(T_INDEX_PRIM),Primtive.at(RHO_INDEX_PRIM), GetMassFractions());
  //Pressure /= config->GetGas_Constant_Ref();

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if(Pressure < EPS)
    return true;
  return false;
}

//
//
/*!
 *\brief Set sound speed (requires SetDensity() call)
 *///
//
bool CReactiveEulerVariable::SetSoundSpeed(CConfig* config) {
  su2double dim_temp = Primitive.at(T_INDEX_PRIM)*config->GetTemperature_Ref();
  if(config->GetSystemMeasurements() == US)
    dim_temp *= 5.0/9.0;
  su2double Sound_Speed;
  //su2double Sound_Speed = library->ComputeFrozenSoundSpeed(dim_temp, GetMassFractions(), Primitive.at(P_INDEX_PRIM), Primitive.at(RHO_INDEX_PRIM));

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(A_INDEX_PRIM) = Sound_Speed;
  if(Sound_Speed < EPS)
    return true;
  return false;
}

//
//
/*!
 *\brief Set the velocity vector from old solution
 */
//
//
inline void CReactiveEulerVariable::SetVelocity_Old(su2double* val_velocity) {
  SU2_Assert(val_velocity != NULL,"The array of velocity val_velocity has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array of solution variables at previous step has not been allocated");

  for(unsigned short iDim=0; iDim < nDim; ++iDim)
    Solution_Old[RHOVX_INDEX_SOL+iDim] = val_velocity[iDim]*Primitive.at(RHO_INDEX_PRIM);
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                                         unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config):
                                         CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,val_nprimvarlim,config) {
  Laminar_Viscosity = 0.0;
  Thermal_Conductivity = 0.0;
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                         const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies,unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, CConfig* config):
                                         CReactiveEulerVariable(val_pressure,val_massfrac,val_velocity,val_temperature,val_nDim,
                                                                val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,val_nprimvarlim,config) {
  su2double dim_temp, dim_press;
  bool US_System = config->GetSystemMeasurements() == US;
  dim_temp = val_temperature*config->GetTemperature_Ref();
  dim_press = val_pressure*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }
  //Laminar_Viscosity = library->GetLambda(dim_temp,val_massfrac)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,val_massfrac)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);
}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                         unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                         CConfig* config):
                                         CReactiveEulerVariable(val_solution,val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                val_nprimvarlim,config) {
  //RealVec YS(val_solution.begin() + CReactiveEulerVariable::RHOS_INDEX_PRIM,val_solution.end());
  //Ys /= val_solution[CReactiveEulerVariable::RHO_INDEX_SOL];
  su2double dim_temp, dim_press;
  bool US_System = config->GetSystemMeasurements() == US;
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);
}

//
//
/*!
 *\brief Set primitive variables
 */
//
//
bool CReactiveNSVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables using Euler version since primitives are the same ---*/
  bool nonPhys = CReactiveEulerVariable::SetPrimVar(config);

  /*--- Compute transport properties --- */
  su2double dim_temp, dim_press;
  bool US_System = config->GetSystemMeasurements() == US;
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  return nonPhys;
}
