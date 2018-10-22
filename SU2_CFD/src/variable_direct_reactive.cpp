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

CReactiveEulerVariable::CReactiveEulerVariable():CVariable(),nSpecies() {
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

CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim,unsigned short val_nvar, std::shared_ptr<CConfig> config):
                        CVariable(val_nDim,val_nvar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
                        nSpecies(library->GetNSpecies()) {
  nPrimVar = nSpecies + nDim + 5;
  nPrimVarGrad = nDim + 2;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;

  Primitive.resize(nPrimVar);
  Gradient_Primitive.resize(nPrimVarGrad,RealVec(nDim));
  Limiter_Primitive.resize(nPrimVarGrad);

  Solution_Min = new su2double[nPrimVarGrad];
  Solution_Max = new su2double[nPrimVarGrad];

  /*
  xi   = Common::wrap_in_unique(config->GetRotationModes());
  Ms   = Common::wrap_in_unique(config->GetMolar_Mass());
  hf   = Common::wrap_in_unique(config->GetEnthalpy_Formation());
  Tref = Common::wrap_in_unique(config->GetRefTemperature());

  Ri.reserve(nSpecies);
  const su2double Ru = (*library).R_ungas;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Ri.push_back(Ru/Ms[iSpecies]);

  library->SetRiGas(Ri);
  */

}

//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//

bool CReactiveEulerVariable::SetPrimVar(CConfig* config) {

  bool nonPhys, bkup;
  unsigned short iVar;

  /*--- Convert conserved to primitive variables ---*/
  nonPhys = Cons2PrimVar(config, Solution, Primitive.data(), dPdU.data(), dTdU.data());
  if (nonPhys) {
    for (iVar = 0; iVar < nVar; ++iVar)
      Solution[iVar] = Solution_Old[iVar];
    bkup = Cons2PrimVar(config, Solution, Primitive.data(), dPdU.data(), dTdU.data());
  }

  //SetVelocity2();

  return nonPhys;

  //	unsigned short iVar, iSpecies;
  //  bool check_dens, check_press, check_sos, check_temp, RightVol;
  //
  //  /*--- Initialize booleans that check for physical solutions ---*/
  //  check_dens  = false;
  //  check_press = false;
  //  check_sos   = false;
  //  check_temp  = false;
  //  RightVol    = true;
  //
  //  /*--- Calculate primitive variables ---*/
  //  // Solution:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  //  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  //  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  //  SetDensity();                             // Compute species & mixture density
  //	SetVelocity2();                           // Compute the square of the velocity (req. mixture density).
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    check_dens = ((Solution[iSpecies] < 0.0) || check_dens);  // Check the density
  //  check_temp  = SetTemperature(config);     // Compute temperatures (T & Tve) (req. mixture density).
  // 	check_press = SetPressure(config);        // Requires T & Tve computation.
  //  CalcdPdU(Primitive, config, dPdU);        // Requires density, pressure, rhoCvtr, & rhoCvve.
  //  CalcdTdU(Primitive, config, dTdU);
  //  CalcdTvedU(Primitive, config, dTvedU);
  //  check_sos   = SetSoundSpeed(config);      // Requires density, pressure, rhoCvtr, & rhoCvve.
  //
  //  /*--- Check that the solution has a physical meaning ---*/
  //  if (check_dens || check_press || check_sos || check_temp) {
  //
  //    /*--- Copy the old solution ---*/
  //    for (iVar = 0; iVar < nVar; iVar++)
  //      Solution[iVar] = Solution_Old[iVar];
  //
  //    /*--- Recompute the primitive variables ---*/
  //    SetDensity();                           // Compute mixture density
  //    SetVelocity2();                         // Compute square of the velocity (req. mixture density).
  //    check_temp  = SetTemperature(config);   // Compute temperatures (T & Tve)
  //    check_press = SetPressure(config);      // Requires T & Tve computation.
  //    CalcdPdU(Primitive, config, dPdU);                       // Requires density, pressure, rhoCvtr, & rhoCvve.
  //    CalcdTdU(Primitive, config, dTdU);
  //    CalcdTvedU(Primitive, config, dTvedU);
  //    check_sos   = SetSoundSpeed(config);    // Requires density & pressure computation.
  //
  //    RightVol = false;
  //  }
  //
  //  SetEnthalpy();                            // Requires density & pressure computation.
  //
  //
  //
  //  return RightVol;

}

//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//
bool CReactiveEulerVariable::Cons2PrimVar(CConfig* config, su2double* U, su2double* V, su2double* dPdU, su2double* dTdU) {

  bool errT, NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter, maxBTry, nBTry;
  su2double rho, rhoE, rhoE_f, rhoE_ref;
  su2double sqvel, rhoCvtr;
  su2double f, df, NRtol, Btol;
  su2double T,Told,Tnew,hs,hs_old,Tmin,Tmax;
  su2double radical2;

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoe, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, rhoCvtr, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  errT    = false;
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 8E4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0E-6;    // Tolerance for the Newton-Raphson method
  Btol     = 1.0E-4;    // Tolerance for the Bisection method
  maxNIter = 5;        // Maximum Newton-Raphson iterations
  maxBIter = 32;        // Maximum Bisection method iterations
  maxBTry  = 4;       // Maximum Bisection method attempts
  nBTry    = 0;       // Number Bisection method attempts done

  /*--- Rename variables for convenience ---*/
  rhoE   = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species & mixture density ---*/
  // Note: if any species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  V[RHO_INDEX_PRIM] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if (U[RHOS_INDEX_SOL + iSpecies] < 0.0) {
      V[RHOS_INDEX_PRIM + iSpecies] = 1E-20;
      U[RHOS_INDEX_SOL + iSpecies] = 1E-20;
      nonPhys = true;
    }
    else
      V[RHOS_INDEX_PRIM+iSpecies] = U[RHOS_INDEX_SOL + iSpecies];

    V[RHO_INDEX_PRIM]  += U[RHOS_INDEX_SOL + iSpecies];
  }
  U[RHO_INDEX_SOL] = V[RHO_INDEX_SOL];

  // Rename for convenience
  rho = V[RHO_INDEX_PRIM];
  /*--- Assign mixture velocity ---*/
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
    sqvel += V[VX_INDEX_PRIM + iDim]*V[VX_INDEX_PRIM + iDim];
  }

  /*--- Translational-Rotational Temperature ---*/
  RealVec Yi(U + RHOS_INDEX_SOL,U + (RHOS_INDEX_SOL + nSpecies));
  std::for_each(Yi.begin(),Yi.end(),[=](su2double& elem){elem /= rho;});
  library->SetRgas(Yi,Ri);
  const su2double Rgas = library->GetRgas();
  const su2double C1 = (-U[RHOE_INDEX_SOL] + rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a Newton-Raphson root-finding method to find T ---*/
    hs = 0.0;
    hs_old = 0.0;
    for(iSpecies = 0; iSpecies<nSpecies; ++iSpecies)  {
      hs_old += Yi[iSpecies]*CalcHs(Told,iSpecies);
      hs += Yi[iSpecies]*CalcHs(T,iSpecies);
    }

    f = T - C1 - C2*hs;
    df = T - Told + C2*(hs_old-hs);
    Tnew = T - f*(T-Told)/df;

    // Check for convergence
    if (std::abs(Tnew - T) < NRtol) {
      NRconvg = true;
      break;
    }
    else {
      Told = T;
      T = Tnew;
    }
  }

  // If the Newton-Raphson method has converged, assign the value of T.
  // Otherwise, execute a bisection root-finding method

  if (NRconvg)
    V[T_INDEX_PRIM] = T;
  else {

    // Execute the bisection root-finding method
    Bconvg = false;
    su2double Ta = Tmin;
    su2double Tb = Tmax;
    while(!Bconvg and nBTry < maxBTry) {
      for (iIter = 0; iIter < maxBIter; ++iIter) {
        T = (Ta + Tb)/2.0;
        hs = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
          hs += Yi[iSpecies] * CalcHs(T,iSpecies);
        f = T - C1 - C2*hs;

        if (std::abs(f) < Btol) {
          V[T_INDEX_PRIM] = T;
          Bconvg = true;
          break;
        }
        else {
          if (f > 0)
            Ta = T;
          else
            Tb = T;
        }
      }

      // If absolutely no convergence, then try increasing number of iterations

      if (!Bconvg) {
        maxBIter *= 2;
        nBTry++;
      }
    }

  }

  // Determine properties of the mixture at the current state
  /*
  rhoE_f   = 0.0;
  rhoE_ref = 0.0;
  rhoCvtr  = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; ++iSpecies++) {
    rhoCvtr  += U[RHOS_INDEX_SOL + iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ri[iSpecies];
    rhoE_ref += U[RHOS_INDEX_SOL + iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ri[iSpecies] * Tref[iSpecies];
    rhoE_f   += U[RHOS_INDEX_SOL + iSpecies] * (hf[iSpecies] - Ri[iSpecies]*Tref[iSpecies]);
  }

  // Calculate translational-rotational temperature
  V[T_INDEX_PRIM] = (rhoE - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;
  V[RHOCV_INDEX_PRIM] = rhoCvtr;

  // Determine if the temperature lies within the acceptable range
  if (V[T_INDEX_PRIM] < Tmin) {
    V[T_INDEX_PRIM] = Tmin;
    nonPhys = true;
    errT    = true;
  }
  else if (V[T_INDEX_PRIM] > Tmax){
    V[T_INDEX] = Tmax;
    nonPhys = true;
    errT    = true;
  }
  */

  /*--- Vibrational-Electronic Temperature ---*/

  // Check for non-physical solutions
  /*
  rhoEve_min = 0.0;
  rhoEve_max = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhoEve_min += U[iSpecies]*CalcEve(config, Tvemin, iSpecies);
    rhoEve_max += U[iSpecies]*CalcEve(config, Tvemax, iSpecies);
  }
  if (rhoEve < rhoEve_min) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemin;
  } else if (rhoEve > rhoEve_max) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemax;
  } else {
  */
    /*--- Execute a Newton-Raphson root-finding method to find Tve ---*/
    // Initialize to the translational-rotational temperature
    //Tve   = V[T_INDEX];

    // Execute the root-finding method
    //NRconvg = false;
//    for (iIter = 0; iIter < maxNIter; iIter++) {
//      rhoEve_t = 0.0;
//      rhoCvve  = 0.0;
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        val_eves[iSpecies]  = CalcEve(config, Tve, iSpecies);
//        val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
//        rhoEve_t += U[iSpecies]*val_eves[iSpecies];
//        rhoCvve  += U[iSpecies]*val_Cvves[iSpecies];
//      }
//
//      // Find the root
//      f  = U[nSpecies+nDim+1] - rhoEve_t;
//      df = -rhoCvve;
//      Tve2 = Tve - (f/df)*scale;
//
//      // Check for nonphysical steps
//      if ((Tve2 < Tvemin) || (Tve2 > Tvemax))
//        break;
////      if (Tve2 < Tvemin)
////        Tve2 = Tvemin;
////      else if (Tve2 > Tvemax)
////        Tve2 = Tvemax;
//
//      // Check for convergence
//      if (fabs(Tve2 - Tve) < NRtol) {
//        NRconvg = true;
//        break;
//      } else {
//        Tve = Tve2;
//      }
//    }

    // If the Newton-Raphson method has converged, assign the value of Tve.
    // Otherwise, execute a bisection root-finding method
    /*
    if (NRconvg)
      V[TVE_INDEX] = Tve;
    else {

      // Assign the bounds
      Tve_o = Tvemin;
      Tve2  = Tvemax;

      // Execute the root-finding method
      Bconvg = false;
      for (iIter = 0; iIter < maxBIter; iIter++) {
        Tve      = (Tve_o+Tve2)/2.0;
        rhoEve_t = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies] = CalcEve(config, Tve, iSpecies);
          rhoEve_t          += U[iSpecies] * val_eves[iSpecies];
        }

        if (fabs(rhoEve_t - U[nSpecies+nDim+1]) < Btol) {
          V[TVE_INDEX] = Tve;
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            val_Cvves[iSpecies] = CalcCvve(Tve, config, iSpecies);
          Bconvg = true;
          break;
        } else {
          if (rhoEve_t > rhoEve) Tve2 = Tve;
          else                  Tve_o = Tve;
        }
      }

      // If absolutely no convergence, then assign to the TR temperature
      if (!Bconvg) {
        V[TVE_INDEX] = V[T_INDEX];
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies]  = CalcEve(config, V[TVE_INDEX], iSpecies);
          val_Cvves[iSpecies] = CalcCvve(V[TVE_INDEX], config, iSpecies);
        }
      }
    }
  }
  */
  /*--- If there are clipped temperatures, correct the energy terms ---*/
//  if (errT) {
//    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
//                       - rhoE_ref + 0.5*rho*sqvel;
//  }
//  if (errTve) {
//    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
//                       - rhoE_ref + 0.5*rho*sqvel;
//    U[nSpecies+nDim+1] = rhoCvve*V[TVE_INDEX];
//  }


  /*--- Pressure ---*/
  V[P_INDEX_PRIM] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    V[P_INDEX_PRIM] += U[RHOS_INDEX_SOL + iSpecies]*Ri[iSpecies]*V[T_INDEX_PRIM];
  if (V[P_INDEX_PRIM] < 0.0) {
    V[P_INDEX_PRIM] = 1E-20;
    nonPhys = true;
  }

  /*--- Partial derivatives of pressure and temperature ---*/
  //CalcdPdU(  V, val_eves, config, val_dPdU  );
  //CalcdTdU(  V, config, val_dTdU  );

  /*--- Sound speed ---*/
  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    radical2 += V[RHOS_INDEX_PRIM+iSpecies]/rho * val_dPdU[iSpecies];
  for (iDim = 0; iDim < nDim; ++iDim)
    radical2 += V[VX_INDEX_PRIM+iDim]*val_dPdU[nSpecies+iDim];
  radical2 += (U[RHOE_INDEX_SOL]+V[P_INDEX_PRIM])/rho * val_dPdU[nSpecies+nDim];
  V[A_INDEX_PRIM] = std::sqrt(radical2);

  if (radical2 < 0.0) {
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

  unsigned short iDim, iSpecies;
  su2double T, sqvel, rhoE, Ef, rhos;

  /*--- Rename & initialize for convenience ---*/
  T       = V[T_INDEX_PRIM]; // Translational-rotational temperature [K]
  sqvel   = 0.0;                            // Velocity^2 [m2/s2]
  rhoE    = 0.0;                            // Mixture total energy per unit of mass [J/kg]

  sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM+nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Set mixture density and species density ---*/
  U[RHO_INDEX_SOL] = V[RHO_INDEX_PRIM];
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    U[RHOS_INDEX_SOL + iSpecies] = V[RHOS_INDEX_PRIM + iSpecies];

  /*--- Set momentum ---*/
  for (iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM+iDim];

  /*--- Set the total energy ---*/
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    rhos = U[RHOS_INDEX_SOL + iSpecies];

    // Species formation energy
    Ef = hf[iSpecies] - Ri[iSpecies]*Tref[iSpecies];

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ri[iSpecies] * (T-Tref[iSpecies])
                    + Ef + 0.5*sqvel);

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
  std::vector<su2double*> tmp;
  for(auto&& elem:Gradient_Primitive)
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

  unsigned short iSpecies;
  su2double Density = 0.0;

  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Primitive.at(RHOS_INDEX_PRIM+iSpecies) = Solution[RHOS_INDEX_SOL+iSpecies];
    Density += Solution[RHOS_INDEX_SOL+iSpecies];
  }
  Primitive.at(RHO_INDEX_PRIM) = Density;

  return false;
}

//
//
/*!
 *\brief Set pressure
 *///
//

void CReactiveEulerVariable::SetPressure(void) {

  unsigned short iSpecies;
  su2double Pressure = 0.0;
  RealVec Ms(nSpecies);

  /*--- Read gas mixture properties from library ---*/
  library->GetMolarMasses(Ms);

  /*--- Solve for mixture pressure using ideal gas law & Dalton's law ---*/
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Pressure += Solution[RHOS_INDEX_SOL+iSpecies] * Ri[iSpecies] * Primitive.at(T_INDEX_PRIM);

  /*--- Store computed values and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if (Pressure < 0.0)
    throw std::runtime_error("Non Physical solution");
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
/*--- Compute specific enthalpy for desired species (not including kinetic energy contribution)---*/
//
//

su2double CReactiveEulerVariable::CalcHs(su2double val_T, unsigned short iSpecies) {

  assert(iSpecies < nSpecies);

  su2double T, hs;

  /*--- Rename for convenience ---*/
  T = val_T;

  hs = hf[iSpecies] + (1.0 + 3.0/2.0+xi[iSpecies]/2.0)*Ri[iSpecies]*(T-Tref[iSpecies]);

  return hs;

}

//
//
/*!
 *\brief Compute squared velocity
 */
//
//

/*void CReactiveEulerVariable::SetVelocity2(void) {
  unsigned short iDim;

  Velocity2 = 0.0;
  su2double density = Solution[RHO_INDEX_SOL];
  for (iDim = 0; iDim < nDim; ++iDim) {
    Primitive[VX_INDEX_PRIM+iDim] = Solution[RHOVX_INDEX_SOL + iDim]/density;
    Velocity2 +=  Primitive[VX_INDEX_PRIM+iDim]*Primitive[VX_INDEX_PRIM+iDim];
  }
}
*/

//
//
/*!
 *\brief Compute norm projected velocity
 */
//
//

su2double CReactiveEulerVariable::GetProjVel(su2double* val_vector) {

  su2double ProjVel = 0.0;
  //su2double density = 0.0;
	unsigned short iDim, iSpecies;

  //for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //  density += Solution[RHOS_INDEX_SOL+iSpecies];
	for (iDim = 0; iDim < nDim; ++iDim)
		ProjVel += Solution[RHOVX_INDEX_SOL+iDim]*val_vector[iDim]/Solution[RHO_INDEX_SOL];

	return ProjVel;
}






























































































//
//
/*!
 *\brief Class constructor
 *///
//

CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim,unsigned short val_nvar, std::unique_ptr<CConfig>& config):
                     CReactiveEulerVariable(val_nDim,val_nvar,config) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStream();

}

/*void CReactiveNSVariable::SetMax_Lambda_Visc(const su2double val_max_lambda,const unsigned short iSpecies) {
  Max_Lambda_Visc_MultiSpecies[iSpecies-1] = val_max_lambda;
}

su2double CReactiveNSVariable::GetMax_Lambda_Visc(const unsigned short iSpecies) {
  return Max_Lambda_Visc_MultiSpecies[iSpecies-1];
}*/
