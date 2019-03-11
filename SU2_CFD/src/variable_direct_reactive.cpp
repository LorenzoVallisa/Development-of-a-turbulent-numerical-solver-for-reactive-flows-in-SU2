#include "../include/variable_reactive.hpp"
//#include "../../Common/include/reacting_model_library.hpp"

CReactiveEulerVariable::RealVec CReactiveEulerVariable::Ri = {};

//
//
/*!
 *\brief Class default constructor
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(): CVariable(), nSpecies(), nPrimVarLim(), US_System(false), Cp() {
  nPrimVar = 0;
  nPrimVarGrad = 0;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                               unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                               LibraryPtr lib_ptr, CConfig* config): CVariable(val_nDim,val_nvar,config), Cp(),
                                                                                     library(lib_ptr) {
  nSpecies = val_nSpecies;
  nPrimVar = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nPrimVarLim = val_nprimvarlim;

  US_System = (config->GetSystemMeasurements() == US);

  Primitive.resize(nPrimVar);
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for(unsigned short iVar = 0; iVar < nPrimVarGrad; ++iVar)
    Gradient_Primitive[iVar] = new su2double[nDim];
  Limiter_Primitive.resize(nPrimVarLim);

  dTdU.resize(nVar);
  dPdU.resize(nVar);

  Ys.resize(nSpecies);
  Ri = library->GetRigas();
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Ri[iSpecies] /= config->GetGas_Constant_Ref();
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
                                               unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                               CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad,
                                                                      val_nprimvarlim, lib_ptr, config) {

  /*--- Rename and initialize for convenience ---*/
  su2double T = val_temperature;   // Translational-rotational temperature [K]
  su2double P = val_pressure;

  su2double rho,rhoE;

  /*--- Compute mixture density ---*/
  su2double Rgas = std::inner_product(Ri.cbegin(), Ri.cend(), val_massfrac.cbegin(), 0.0);
  rho = P/(Rgas*T);

  /*--- Compute energy (RHOE) from supplied primitive quanitites ---*/
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  /*
  su2double sqvel = std::inner_product(val_velocity.cbegin(), val_velocity.cend(), val_velocity.cbegin(), 0.0);
  su2double e_tot = library->ComputeEnergy(dim_temp,val_massfrac)/config->GetEnergy_Ref();
  if(US_System)
    e_tot *= 3.28084*3.28084;
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
  Solution[RHOE_INDEX_SOL] = Solution_Old[RHOE_INDEX_SOL] = rhoE;

  /*--- Assign primitive variables ---*/
  Primitive.at(T_INDEX_PRIM) = T;
  Primitive.at(P_INDEX_PRIM) = P;
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
                                               unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                               CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad,
                                                                      val_nprimvarlim, lib_ptr, config) {
  /*--- Initialize Solution and Solution_Old vectors ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  std::copy(val_solution.cbegin(),val_solution.cend(),Solution);
  std::copy(val_solution.cbegin(),val_solution.cend(),Solution_Old);

  /*--- Initialize T to the free stream for the secant method ---*/
  Primitive.at(T_INDEX_PRIM) = config->GetTemperature_FreeStream();
}

//
//
/*!
 *\brief Class destructor
 */
//
//
CReactiveEulerVariable::~CReactiveEulerVariable() {
  if(Gradient_Primitive != NULL) {
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      delete[] Gradient_Primitive[iVar];
    delete[] Gradient_Primitive;
  }
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
    bool check_old = Cons2PrimVar(config, Solution, Primitive.data());
    SU2_Assert(check_old == true, "Neither the old solution is feasible to set primitive variables: problem unsolvable");
  }

  /*--- Set specific heat at constant pressure ---*/
  su2double Rgas = std::inner_product(Ri.cbegin(), Ri.cend(), Ys.cbegin(), 0.0);
  su2double T = Primitive[T_INDEX_PRIM];
  su2double a = Primitive[A_INDEX_PRIM];
  Cp = a*a*Rgas/(a*a - Rgas*T);

  /*--- Set temperature and pressure derivatives ---*/
  CalcdTdU(Primitive.data(), config, dTdU.data());
  CalcdPdU(Primitive.data(), config, dPdU.data());

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
  // V: [T, u, v, w, P, rho, h, a, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 6.0e4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0e-6;    // Tolerance for the Secant method
  Btol     = 1.0e-4;    // Tolerance for the Bisection method
  maxNIter = 7;        // Maximum Secant method iterations
  maxBIter = 32;        // Maximum Bisection method iterations

  /*--- Rename variables forconvenience ---*/
  rhoE = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species mass fraction and mixture density ---*/
  // NOTE: If any species densities are < 0, these values are re-assigned
  //       in the conserved vector to ensure positive density
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < EPS) {
      U[RHOS_INDEX_SOL + iSpecies] = EPS;
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
    V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL]/U[RHO_INDEX_SOL];

  /*--- Checking sum of mass fraction ---*/
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  nonPhys = nonPhys || (std::abs(std::accumulate(Ys.cbegin(),Ys.cend(),0.0) - 1.0) > EPS);

  /*--- Rename for convenience ---*/
  rho = V[RHO_INDEX_PRIM];

  /*--- Assign mixture velocity ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
  sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Translational-Rotational Temperature ---*/
  const su2double Rgas = std::inner_product(Ri.cbegin(), Ri.cend(), Ys.cbegin(), 0.0);
  const su2double C1 = (-rhoE + rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a secant root-finding method to find T ---*/
    su2double dim_temp, dim_temp_old;
    dim_temp = T*config->GetTemperature_Ref();
    dim_temp_old = Told*config->GetTemperature_Ref();
    if(US_System) {
      dim_temp *= 5.0/9.0;
      dim_temp_old *= 5.0/9.0;
    }
    //hs_old = library->ComputeEnthalpy(dim_temp_old,Ys)/config->GetEnergy_Ref();
    //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
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

  /*--- If the secant method has converged, assign the value of T.
        Otherwise execute a bisection root-finding method ---*/
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
      //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
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
  //V[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp, Ys, V[P_INDEX_PRIM], rho);
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
 *\brief Compute temperature derivatives
 */
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
  su2double rhoCv = rho*(Cp - std::inner_product(Ri.cbegin(), Ri.cend(), V + RHOS_INDEX_PRIM, 0.0));
  su2double sq_vel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //Int_Energies = library->ComputePartialEnergy(temp);
  for(auto& elem: Int_Energies)
    elem /= config->GetEnergy_Ref();
  if(US_System) {
    for(auto& elem: Int_Energies)
      elem *= 3.28084*3.28084;
  }

  /*--- Set temperature derivatives ---*/
  dTdU[RHO_INDEX_SOL] = 0.5*sq_vel/rhoCv;
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    dTdU[RHOVX_INDEX_SOL + iDim] = -V[VX_INDEX_PRIM + iDim]/rhoCv;
  dTdU[RHOE_INDEX_SOL] = 1.0/rhoCv;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    dTdU[RHOS_INDEX_SOL + iSpecies] = -Int_Energies[iSpecies]/rhoCv;
}

//
//
/*!
 *\brief Compute pressure derivatives
 */
//
//
void CReactiveEulerVariable::CalcdPdU(su2double* V, CConfig* config, su2double* dPdU) {
  /*--- Check memory allocation---*/
  if(config->GetExtIter() == 0)
    SU2_Assert(dPdU != NULL,"The array to store the derivatives of pressure w.r.t. conserved variables has not been allocated");

  /*--- Useful quantities ---*/
  su2double Gamma = Cp/(Cp - std::inner_product(Ri.cbegin(), Ri.cend(), V + RHOS_INDEX_PRIM, 0.0));
  su2double sq_vel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Rename for convenience ---*/
  su2double T = V[T_INDEX_PRIM];
  su2double rho = V[RHO_INDEX_PRIM];

  /*--- Set pressure derivatives ---*/
  dPdU[RHO_INDEX_SOL] = (Gamma - 1.0)*0.5*sq_vel;
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    dPdU[RHOVX_INDEX_SOL + iDim] = (1.0 - Gamma)*V[VX_INDEX_PRIM + iDim];
  dPdU[RHOE_INDEX_SOL] = Gamma - 1.0;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    dPdU[RHOS_INDEX_SOL + iSpecies] = Ri[iSpecies]*T - (Gamma - 1.0)*Int_Energies[iSpecies];
}


//
//
/*!
 *\brief Pass from primitive to conserved variables
 */
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
    U[RHOS_INDEX_SOL + iSpecies] = V[RHO_INDEX_PRIM]*V[RHOS_INDEX_PRIM];

  /*--- Set momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM + iDim];

  /*--- Set energies ---*/
  //su2double sigma = std::accumulate(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies));
  //U[RHOE_INDEX_SOL] = V[RHO_INDEX_PRIM]*V[H_INDEX_PRIM] - V[P_INDEX_PRIM]*sigma;
  U[RHOE_INDEX_SOL] = V[RHO_INDEX_PRIM]*V[H_INDEX_PRIM] - V[P_INDEX_PRIM];
}

//
//
/*!
 *\brief Set density
 */
//
//
inline bool CReactiveEulerVariable::SetDensity(void) {
  /*--- Check memory location for safety ---*/
  SU2_Assert(Solution != NULL,"The array of solution variables has not been allocated");

  Primitive.at(RHO_INDEX_PRIM) = Solution[RHO_INDEX_SOL];

  if(Primitive.at(RHO_INDEX_PRIM) < EPS)
    return true;

  return false;
}

//
//
/*!
 *\brief Set pressure (NOTE: it requires SetDensity() call)
 */
//
//
bool CReactiveEulerVariable::SetPressure(CConfig* config) {
  /*--- Compute mixture pressure ---*/
  su2double Rgas = std::inner_product(Ri.cbegin(), Ri.cend(), Primitive.cbegin() + RHOS_INDEX_PRIM, 0.0);
  su2double Pressure = Primitive.at(RHO_INDEX_PRIM)*Rgas*Primitive.at(T_INDEX_PRIM);

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if(Pressure < EPS)
    return true;

  return false;
}

//
//
/*!
 *\brief Set sound speed (NOTE: it requires SetDensity() call)
 */
//
//
bool CReactiveEulerVariable::SetSoundSpeed(CConfig* config) {
  /*--- Compute useful quantities ---*/
  su2double Rgas = std::inner_product(Ri.cbegin(), Ri.cend(), Primitive.cbegin() + RHOS_INDEX_PRIM, 0.0);
  su2double Gamma = Cp/(Cp - Rgas);

  /*--- Compute frozen sound spped ---*/
  su2double Sound_Speed = std::sqrt(Gamma*Rgas*Primitive.at(T_INDEX_PRIM));

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
  /*--- Check memory location for safety ---*/
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
CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                         unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                         LibraryPtr lib_ptr, CConfig* config): CReactiveEulerVariable(val_nDim, val_nvar, val_nSpecies,
                                                                               val_nprimvar, val_nprimvargrad, val_nprimvarlim, lib_ptr,
                                                                               config), Laminar_Viscosity(), Thermal_Conductivity() {
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
                                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_pressure, val_massfrac, val_velocity, val_temperature, val_nDim,
                                                                val_nvar, val_nSpecies, val_nprimvar, val_nprimvargrad, val_nprimvarlim,
                                                                lib_ptr, config), Laminar_Viscosity(), Thermal_Conductivity() {
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config):
                                         CReactiveEulerVariable(val_solution, val_nDim, val_nvar, val_nSpecies, val_nprimvar,
                                                                val_nprimvargrad, val_nprimvarlim, lib_ptr, config),
                                                                Laminar_Viscosity(), Thermal_Conductivity() {
  Diffusion_Coeffs.resize(nSpecies,nSpecies);
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

  su2double dim_temp, dim_press;
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  /*--- Compute transport properties --- */
  Ys = GetMassFractions();
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  return nonPhys;
}
