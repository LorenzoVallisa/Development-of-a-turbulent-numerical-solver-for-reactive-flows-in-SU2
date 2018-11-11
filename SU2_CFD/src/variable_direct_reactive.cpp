#include "../include/variable_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/move_pointer.hpp"
#include "../../Common/include/default_initialization.hpp"
#include "../../Common/include/su2_assert.hpp"

//
//
/*!
 *\brief Class default constructor
 */
//
//

CReactiveEulerVariable::CReactiveEulerVariable():CVariable(),nSpecies() {
  std::tie(nPrimVar,nPrimVarGrad,nSecondaryVar,nSecondaryVarGrad) = Common::repeat<4,decltype(nPrimVar)>(decltype(nPrimVar)());
}

//
//
/*!
 *\brief Class constructor
 */
//
//

CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim,unsigned short val_nvar, std::shared_ptr<CConfig> config):
                        CVariable(val_nDim,val_nvar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
                        nSpecies(library->GetNSpecies()) {
  nPrimVar = nSpecies + nDim + 5;
  nPrimVarGrad = nDim + 2;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;

  Primitive.resize(nPrimVar);
  Gradient_Primitive.resize(nPrimVarGrad,nDim);
  Limiter_Primitive.resize(nPrimVarGrad);

  Limiter = new su2double [nVar];
  Solution_Max = new su2double[nPrimVarGrad];
  Solution_Min = new su2double[nPrimVarGrad];

}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(su2double val_pressure,RealVec& val_massfrac,RealVec& val_velocity,
                                               su2double val_temperature,unsigned short val_nDim,unsigned short val_nvar,
                                               std::shared_ptr<CConfig> config): CReactiveEulerVariable(val_nDim,val_nvar,config) {

  /*--- Load variables from the config class --*/
  /*
  xi        = config->GetRotationModes();      // Rotational modes of energy storage
  Ms        = config->GetMolar_Mass();         // Species molar mass
  Tref      = config->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = config->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  */

  /*--- Rename and initialize for convenience ---*/
  su2double T = val_temperature;   // Translational-rotational temperature [K]
  su2double P = val_pressure;

  RealVec Yi = val_massfrac;
  //library->SetRiGas(Yi);
  //su2double Rgas = library->GetRgas();
  su2double rho,rhoE;

  /*--- Compute mixture density ---*/
  //rho = library->ComputeDensity(T,P);

  su2double Gamma = 0.0,Sound_Speed = 0.0;
  library->Gamma_FrozenSoundSpeed(T,P,rho,Gamma,Sound_Speed);

  /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
  /*
  su2double sqvel = 0.0;
  for (iDim = 0; iDim < nDim; ++iDim)
    sqvel += val_mach[iDim]*Sound_Speed * val_mach[iDim]*Sound_speed;

  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    // Species density
    rhos = val_massfrac[iSpecies]*rho;

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);

  }

	for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];
  }
  */

  /*--- Initialize Solution and Solution_Old vectors ---*/

  /*--- Initialize mixture density and partial density ---*/
  Solution[RHO_INDEX_SOL] = Solution_Old[RHO_INDEX_SOL] = rho;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Solution[RHOS_INDEX_SOL + iSpecies] = Solution_Old[RHOS_INDEX_SOL + iSpecies] = rho*val_massfrac[iSpecies];

  /*--- Initialize momentum ---*/
  for (unsigned short iDim = 0; iDim < nDim; ++iDim)
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
CReactiveEulerVariable::CReactiveEulerVariable(RealVec& val_solution,unsigned short val_nDim,unsigned short val_nvar,
                                               std::shared_ptr<CConfig> config): CReactiveEulerVariable(val_nDim,val_nvar,config) {

  /*--- Initialize Solution and Solution_Old vectors ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  std::copy(val_solution.cbegin(),val_solution.cend(),Solution);
  std::copy(val_solution.cbegin(),val_solution.cend(),Solution_Old);

  /*--- Initialize T and P to the free stream for Newton-Raphson method ---*/
  Primitive[T_INDEX_PRIM] = config->GetTemperature_FreeStream();
  Primitive[P_INDEX_PRIM] = config->GetPressure_FreeStream();

}



//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//

bool CReactiveEulerVariable::SetPrimVar(CConfig* config) {

  /*--- Convert conserved to primitive variables ---*/
  bool nonPhys = Cons2PrimVar(Solution, Primitive.data());
  if(nonPhys)
    std::copy(Solution_Old,Solution_Old + nVar,Solution);

  return nonPhys;

}

//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//
bool CReactiveEulerVariable::Cons2PrimVar(su2double* U, su2double* V) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  bool NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter;
  su2double rho, rhoE;
  su2double sqvel;
  su2double f, df, NRtol, Btol;
  su2double T,Told,Tnew,hs,hs_old,Tmin,Tmax;
  RealVec hs_species(nSpecies),hs_old_species(nSpecies);

  /*--- Conserved and primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoE, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 8E4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0E-6;    // Tolerance forthe Newton-Raphson method
  Btol     = 1.0E-4;    // Tolerance forthe Bisection method
  maxNIter = 5;        // Maximum Newton-Raphson iterations
  maxBIter = 32;        // Maximum Bisection method iterations

  /*--- Rename variables forconvenience ---*/
  rhoE  = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species and mixture density ---*/
  // Note: ifany species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  V[RHO_INDEX_PRIM] = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < 0.0) {
      V[RHOS_INDEX_PRIM + iSpecies] = 1E-20;
      U[RHOS_INDEX_SOL + iSpecies] = 1E-20;
      nonPhys = true;
    }
    else
      V[RHOS_INDEX_PRIM+iSpecies] = U[RHOS_INDEX_SOL + iSpecies];

    V[RHO_INDEX_PRIM]  += U[RHOS_INDEX_SOL + iSpecies];
  }
  U[RHO_INDEX_SOL] = V[RHO_INDEX_SOL];

  // Rename forconvenience
  rho = V[RHO_INDEX_PRIM];
  /*--- Assign mixture velocity ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
  sqvel = std::inner_product(V + VX_INDEX_PRIM,V + (VX_INDEX_PRIM + nDim),V + VX_INDEX_PRIM,0.0);

  /*--- Translational-Rotational Temperature ---*/
  RealVec Yi(U + RHOS_INDEX_SOL,U + (RHOS_INDEX_SOL + nSpecies));
  std::for_each(Yi.begin(),Yi.end(),[=](su2double elem){elem /= rho;});
  //library->SetRgas(Yi);
  const su2double Rgas = library->GetRgas();
  const su2double C1 = (-rhoE + rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a Newton-Raphson root-finding method to find T ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)  {
      hs_old_species[iSpecies] = CalcHs(Told,iSpecies);
      hs_species[iSpecies] = CalcHs(T,iSpecies);
    }
    hs_old = std::inner_product(Yi.cbegin(),Yi.cend(),hs_old_species.cbegin(),0.0);
    hs = std::inner_product(Yi.cbegin(),Yi.cend(),hs_species.cbegin(),0.0);

    f = T - C1 - C2*hs;
    df = T - Told + C2*(hs_old-hs);
    Tnew = T - f*(T-Told)/df;

    // Check forconvergence
    if(std::abs(Tnew - T) < NRtol) {
      NRconvg = true;
      break;
    }
    else {
      Told = T;
      T = Tnew;
    }
  }

  // ifthe Newton-Raphson method has converged, assign the value of T.
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
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        hs_species[iSpecies] = CalcHs(T,iSpecies);
      hs = std::inner_product(Yi.cbegin(),Yi.cend(),hs_species.cbegin(),0.0);
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

    // ifabsolutely no convergence, then something is going really wrong
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
  V[P_INDEX_PRIM] = 0.0;
  //library->ComputePressure(V[RHO_INDEX_PRIM],V[T_INDEX_PRIM],V[P_INDEX_PRIM]);
  if(V[P_INDEX_PRIM] < 0.0) {
    V[P_INDEX_PRIM] = 1E-20;
    nonPhys = true;
  }

  /*--- Partial derivatives of pressure and temperature ---*/
  //CalcdPdU(V, config, val_dPdU);
  //CalcdTdU(V, config, val_dTdU);

  /*--- Sound speed ---*/
  su2double Gamma = 0.0;
  su2double Sound_Speed = 0.0;
  library->Gamma_FrozenSoundSpeed(V[T_INDEX_PRIM],V[P_INDEX_PRIM],V[RHO_INDEX_PRIM],Gamma,Sound_Speed);
  V[A_INDEX_PRIM] = Sound_Speed;

  if(Sound_Speed < 0.0) {
    nonPhys = true;
    V[A_INDEX_PRIM] = EPS;
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
  su2double rhoE;
  //su2double sqvel,T,rhos,Ef;

  /*--- Rename and initialize for convenience ---*/
  //T = V[T_INDEX_PRIM]; // Translational-rotational temperature [K]
  //sqvel = 0.0;
  rhoE = 0.0;

  //sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM+nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Set mixture density and species density ---*/
  U[RHO_INDEX_SOL] = V[RHO_INDEX_PRIM];
  std::copy(V + RHOS_INDEX_SOL,V + (RHOS_INDEX_SOL + nSpecies),U + RHOS_INDEX_SOL);

  /*--- Set momentum ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM+iDim];

  /*--- Set the total energy ---*/
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    //rhos = U[RHOS_INDEX_SOL + iSpecies];

    // Species formation energy
    //Ef = hf[iSpecies] - Ri[iSpecies]*Tref[iSpecies];
    //library->ComputeCp(T);
    // Mixture total energy
    //rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ri[iSpecies] * (T-Tref[iSpecies])
    //                + Ef + 0.5*sqvel);

  }

  /*--- Set energies ---*/
  U[RHOE_INDEX_SOL] = rhoE;

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
  std::copy(Solution + RHOS_INDEX_SOL,Solution + RHOS_INDEX_SOL + nSpecies,Primitive.begin() + RHOS_INDEX_PRIM);
  su2double Density = std::accumulate(Solution + RHOS_INDEX_SOL,Solution + RHOS_INDEX_SOL + nSpecies,0.0);
  Primitive.at(RHO_INDEX_PRIM) = Density;

  if(std::abs(Density - 1.0) < EPS)
    return true;
  return false;

}

//
//
/*!
 *\brief Set pressure
 *///
//

bool CReactiveEulerVariable::SetPressure(CConfig* config) {
  su2double Pressure;
  /*--- Compute this gas mixture property from library ---*/
  //library->GetPressure(Primitive.at(T_INDEX_PRIM),Primitive.at(RHO_INDEX_PRIM),Pressure);

  /*--- Store computed values and check fora physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if(Pressure < 0.0)
    return true;
  return false;
}

//
//
/*!
 *\brief Set sound speed
 *///
//

bool CReactiveEulerVariable::SetSoundSpeed(void) {
  su2double Gamma = 0.0;
  su2double Sound_Speed = 0.0;
  library->Gamma_FrozenSoundSpeed(Primitive.at(T_INDEX_PRIM),Primitive.at(P_INDEX_PRIM),Primitive.at(RHO_INDEX_PRIM),Gamma,Sound_Speed);
  Primitive.at(A_INDEX_PRIM) = Sound_Speed;
  if(Sound_Speed < 0.0)
    return false;
  return true;

}

//
//
/*--- Compute specific enthalpy fordesired species (not including kinetic energy contribution)---*/
//
//

su2double CReactiveEulerVariable::CalcHs(su2double val_T, unsigned short iSpecies) {
  SU2_Assert(iSpecies < nSpecies,"The number of species you want to compute the static enthalpy exceeds the number of species in the mixture");

  double hs;
  //library->ComputeEnthalpy(......);

  return hs;

}

//
//
/*!
 *\brief Compute norm projected velocity
 */
//
//

inline su2double CReactiveEulerVariable::GetProjVel(su2double* val_vector) {
  SU2_Assert(val_vector != NULL,"The array of velocity val_vector has not been allocated");
  su2double ProjVel = 0.0;

  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
		ProjVel += Solution[RHOVX_INDEX_SOL+iDim]*val_vector[iDim]/Solution[RHO_INDEX_SOL];

	return ProjVel;
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

  for(unsigned short iDim=0; iDim<nDim; ++iDim)
    Solution_Old[RHOVX_INDEX_SOL+iDim] = val_velocity[iDim]*Primitive.at(RHO_INDEX_PRIM);
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim,unsigned short val_nvar, std::shared_ptr<CConfig> config):
                     CReactiveEulerVariable(val_nDim,val_nvar,config),Laminar_Viscosity(),Thermal_Conductivity() {
  nSecondaryVarGrad = nSpecies + nDim + 1;
  Diffusion_Coeffs.resize(nSpecies);
}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(su2double val_density, RealVec& val_massfrac, RealVec& val_velocity,
                                         su2double val_temperature,unsigned short val_nDim, unsigned short val_nvar,
                                         std::shared_ptr<CConfig> config): CReactiveEulerVariable(val_density,val_massfrac,val_velocity,
                                                                           val_temperature,val_nDim,val_nvar,config) {
  nPrimVarAvgGrad = nSpecies + nDim + 2;
  AvgGradient_Primitive.resize(nPrimVarAvgGrad);

  Diffusion_Coeffs.resize(nSpecies);
  Laminar_Viscosity = 0.0;
  Thermal_Conductivity = 0.0;
}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                         std::shared_ptr<CConfig> config): CReactiveEulerVariable(val_solution,val_nDim,val_nvar,config) {
  nPrimVarAvgGrad = nSpecies + nDim + 2;
  AvgGradient_Primitive.resize(nPrimVarAvgGrad);

  Diffusion_Coeffs.resize(nSpecies);
  Laminar_Viscosity = 0.0;
  Thermal_Conductivity = 0.0;
}
