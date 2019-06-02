#include "../include/numerics_reactive.hpp"
#include "../../externals/Eigen/IterativeLinearSolvers"

#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/su2_assert.hpp"

#include <algorithm>
#include <iterator>

//
//
/*--- Constructor of the class CUpwReactiveAUSM ---*/
//
//
CUpwReactiveAUSM::CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr):
                  CNumerics(val_nDim, val_nVar, config), library(lib_ptr), nSpecies(library->GetnSpecies()) {
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set local variables to access indices in proper arrays ---*/
  T_INDEX_PRIM    = CReactiveEulerVariable::GetT_INDEX_PRIM();
  VX_INDEX_PRIM   = CReactiveEulerVariable::GetVX_INDEX_PRIM();
  P_INDEX_PRIM    = CReactiveEulerVariable::GetP_INDEX_PRIM();
  RHO_INDEX_PRIM  = CReactiveEulerVariable::GetRHO_INDEX_PRIM();
  H_INDEX_PRIM    = CReactiveEulerVariable::GetH_INDEX_PRIM();
  A_INDEX_PRIM    = CReactiveEulerVariable::GetT_INDEX_PRIM();
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

  /*--- Resize local vectors ---*/
  Phi_i.resize(nVar);
  Phi_j.resize(nVar);
}

//
//
/*--- Compute residual convective term ---*/
//
//
void CUpwReactiveAUSM::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) {
  //AD::StartPreacc();
  //AD::SetPreaccIn(V_i, nSpecies + nDim + 5);
  //AD::SetPreaccIn(V_j, nSpecies + nDim + 5);
  //AD::SetPreaccIn(Normal, nDim);
  if(config->GetExtIter() == 0) {
    SU2_Assert(val_residual != NULL,"The array of residual for convective flux has not been allocated");
    SU2_Assert(V_i != NULL,"The array of primitive variables at node i has not been allocated");
    SU2_Assert(V_j != NULL,"The array of primitive variables at node j has not been allocated");
  }

  unsigned short iDim, iVar, iSpecies; // Indexes for iterations.
	su2double ProjVelocity_i, ProjVelocity_j, // Projected velocities at node i and at node j.
            ProjGridVel_i = 0.0, ProjGridVel_j = 0.0; // Grid velocities at node i and at node j.

  /*--- Face area (norm of the normal vector) and unit normal vector ---*/
  Area = std::inner_product(Normal, Normal + nDim, Normal, 0.0);
  Area = std::sqrt(Area);

  for(iDim = 0; iDim < nDim; ++iDim)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Point i: extract energy, pressure, sound speed and enthalpy ---*/
  Density_i = V_i[RHO_INDEX_PRIM];
  Pressure_i = V_i[P_INDEX_PRIM];
  Enthalpy_i = V_i[H_INDEX_PRIM];
  SoundSpeed_i = V_i[A_INDEX_PRIM];

  /*--- Point j: extract energy, pressure, sound speed and ethalpy  ---*/
  Density_j = V_j[RHO_INDEX_PRIM];
  Pressure_j = V_j[P_INDEX_PRIM];
  Enthalpy_j = V_j[H_INDEX_PRIM];
  SoundSpeed_j = V_j[A_INDEX_PRIM];

  /*--- Projected velocities ---*/
  ProjVelocity_i = std::inner_product(V_i + VX_INDEX_PRIM, V_i + (VX_INDEX_PRIM + nDim), UnitNormal, 0.0);
  ProjVelocity_j = std::inner_product(V_j + VX_INDEX_PRIM, V_j + (VX_INDEX_PRIM + nDim), UnitNormal, 0.0);

  /*--- Grid movement correction ---*/
  bool grid_movement = config->GetGrid_Movement();
  if(grid_movement) {
    ProjGridVel_i = std::inner_product(GridVel_i, GridVel_i + nDim, UnitNormal, 0.0);
    ProjGridVel_j = std::inner_product(GridVel_j, GridVel_j + nDim, UnitNormal, 0.0);
    ProjVelocity_i -= ProjGridVel_i;
    ProjVelocity_j -= ProjGridVel_j;
  }

  /*--- Compute user defined Mach number ---*/
  su2double mRef, mF;
  mRef = config->GetMach();

  /*--- Compute Mean sound speed ---*/
  su2double MeanSoundSpeed = 0.5*(SoundSpeed_i + SoundSpeed_j);

  /*--- Compute Normal L/R Mach numbers ---*/
  su2double mL  = ProjVelocity_i/MeanSoundSpeed;
  su2double mR  = ProjVelocity_j/MeanSoundSpeed;

  /*--- Compute mean local Mach number and reference Mach number ---*/
  mF  = std::sqrt(0.5*(mL*mL + mR*mR));
  mRef = std::min(1.0, std::max(mF,mRef));

  /*--- Set constants ---*/
  const su2double fa = mRef*(2.0 - mRef);
  const su2double alpha = 3.0/16.0*(5.0*fa*fa - 4.0);
  const su2double beta = 0.125;

  su2double mLP, mRM, pLP, pRM;

  /*--- Compute adjusted polynomial function AUSM +-Up ---*/
  if(std::abs(mL) < 1.0) {
    mLP = 0.25*(mL + 1.0)*(mL + 1.0) + beta*(mL*mL - 1.0)*(mL*mL - 1.0);
    pLP = 0.25*(mL + 1.0)*(mL + 1.0)*(2.0 - mL) + alpha*mL*(mL*mL - 1.0);
  }
  else {
    mLP = 0.5*(mL + std::abs(mL));
    pLP = 0.5*(1.0 + std::abs(mL)/mL);
  }

  if(std::abs(mR) < 1.0) {
    mRM = -0.25*(mR - 1.0)*(mR - 1.0) - beta*(mR*mR - 1.0)*(mR*mR - 1.0);
    pRM = 0.25*(mR - 1.0)*(mR - 1.0)*(2.0 + mR) - alpha*mR*(mR*mR - 1.0)*(mR*mR - 1.0);
  }
  else {
    mRM = 0.5*(mR - std::abs(mR));
    pRM = 0.5*(1.0 - std::abs(mR)/mR);
  }

  /*--- Build average state ---*/
  const su2double kP = 0.25;
  const su2double sigma = 1.0;

  su2double m12 = mLP + mRM;
  m12 -= kP/fa*std::max(1.0 - sigma*mF*mF, 0.0)*(Pressure_j - Pressure_i)/(0.5*(Density_i + Density_j)*MeanSoundSpeed*MeanSoundSpeed);
  su2double mLF = 0.5*(m12 + std::abs(m12));
  su2double mRF = 0.5*(m12 - std::abs(m12));
  m12 = MeanSoundSpeed*(mLF*Density_i + mRF*Density_j);

  /*--- Compute the state at node i and at node j ---*/
  Phi_i[RHO_INDEX_SOL] = 1.0;
  Phi_j[RHO_INDEX_SOL] = 1.0;
  for(iDim = 0; iDim < nDim; ++iDim) {
    Phi_i[RHOVX_INDEX_SOL + iDim] = V_i[VX_INDEX_PRIM + iDim];
    Phi_j[RHOVX_INDEX_SOL + iDim] = V_j[VX_INDEX_PRIM + iDim];
  }
  Phi_i[RHOE_INDEX_SOL] = Enthalpy_i;
  Phi_j[RHOE_INDEX_SOL] = Enthalpy_j;
  if(grid_movement) {
    if(std::abs(ProjVelocity_i) > EPS)
      Phi_i[RHOE_INDEX_SOL] += ProjGridVel_i*Pressure_i/(Density_i*ProjVelocity_i);
    if(std::abs(ProjVelocity_j) > EPS)
      Phi_j[RHOE_INDEX_SOL] += ProjGridVel_j*Pressure_j/(Density_j*ProjVelocity_j);
  }
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Phi_i[RHOS_INDEX_SOL + iSpecies] = V_i[RHOS_INDEX_PRIM + iSpecies];
    Phi_j[RHOS_INDEX_SOL + iSpecies] = V_j[RHOS_INDEX_PRIM + iSpecies];
  }

  /*--- Compute the pure convective part of the numerical flux ---*/
  for(iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] = 0.5*(m12*(Phi_i[iVar] + Phi_j[iVar]) + std::abs(m12)*(Phi_i[iVar] - Phi_j[iVar]))*Area;

  /*--- Add to the numerical flux the pressure contribution ---*/
  const su2double Ku = 0.75;
  su2double pLF = pLP*Pressure_i + pRM*Pressure_j;
  pLF -= Ku*pLP*pRM*(Density_i + Density_j)*fa*MeanSoundSpeed*(ProjVelocity_j - ProjVelocity_i);
  for(iDim = 0; iDim < nDim; ++iDim) {
    val_residual[RHOVX_INDEX_SOL + iDim] += pLF*UnitNormal[iDim]*Area;
  }

  if(implicit) {
    /*--- Check memory allocation ---*/
    if(config->GetExtIter() == 0) {
      SU2_Assert(val_Jacobian_i != NULL,"The matrix for convective term Jacobian at node i has not been allocated");
      SU2_Assert(val_Jacobian_j != NULL,"The matrix for convective term Jacobian at node j has not been allocated");
      for(iVar = 0; iVar < nVar; ++iVar) {
        SU2_Assert(val_Jacobian_i[iVar] != NULL,
                   std::string("The row " + std::to_string(iVar) + " of convective term Jacobian at node i has not been allocated"));
        SU2_Assert(val_Jacobian_j[iVar] != NULL,
                   std::string("The row " + std::to_string(iVar) + " of convective term Jacobian at node i has not been allocated"));
      }
    }

    unsigned short jVar;

    /*--- Set Jacobian to zero ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      for(jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
      }
    }

    /*--- Set derivatives of Mach number ---*/
    su2double Mach_Left_Der[nVar], Mach_Right_Der[nVar];

    for(iVar = 0; iVar < nVar; ++iVar) {
      Mach_Left_Der[iVar] = 0.0;
      Mach_Right_Der[iVar] = 0.0;
    }
    Mach_Left_Der[RHO_INDEX_SOL] = -mL/Density_i;
    Mach_Right_Der[RHO_INDEX_SOL] = -mR/Density_j;
    for(iDim = 0; iDim < nDim; ++iDim) {
      Mach_Left_Der[RHOVX_INDEX_SOL + iDim] = UnitNormal[iDim]/(Density_i*MeanSoundSpeed);
      Mach_Right_Der[RHOVX_INDEX_SOL + iDim] = UnitNormal[iDim]/(Density_j*MeanSoundSpeed);
    }

    /*---Set Polynomials Mach derivatives ---*/
    su2double MachPol_Left_Der[nVar], MachPol_Right_Der[nVar];

    if(std::abs(mL) < 1.0) {
      for(iVar = 0; iVar < nVar; ++iVar)
        MachPol_Left_Der[iVar] = Mach_Left_Der[iVar]*(0.5*(mL + 1.0) + 4.0*beta*mL*(mL*mL - 1.0));
    }
    else {
      for(iVar = 0; iVar < nVar; ++iVar)
        MachPol_Left_Der[iVar] = Mach_Left_Der[iVar]*(0.5*(1.0 + std::abs(mL)/mL));
    }

    if(std::abs(mR) < 1.0) {
      for(iVar = 0; iVar < nVar; ++iVar)
        MachPol_Right_Der[iVar] = Mach_Right_Der[iVar]*(0.5*(1.0 - mR) + 4.0*beta*mR*(1.0 - mR*mR));
    }
    else {
      for(iVar = 0; iVar < nVar; ++iVar)
        MachPol_Right_Der[iVar] = Mach_Right_Der[iVar]*(0.5*(1.0 - std::abs(mR)/mR));
    }

    /*---Set scaling factor derivatives ---*/
    su2double Scaling_Left_Der[nVar], Scaling_Right_Der[nVar];

    for(iVar = 0; iVar < nVar; ++iVar) {
      Scaling_Left_Der[iVar] = 0.0;
      Scaling_Right_Der[iVar] = 0.0;
    }

    if(mF == mRef) {
      for(iVar = 0; iVar < nVar; ++iVar) {
        Scaling_Left_Der[iVar]  = Mach_Left_Der[iVar]*mL*(1.0 - mF)/mF;
        Scaling_Right_Der[iVar] = Mach_Right_Der[iVar]*mR*(1.0 - mF)/mF;
      }
    }

    /*--- Set convective extra term derivatives ---*/
    su2double MachExtra_Left_Der[nVar], MachExtra_Right_Der[nVar];
    su2double MeanDensity = 0.5*(Density_i + Density_j);
    su2double factor = std::max(1.0 - sigma*mF*mF, 0.0);
    for(iVar = 0; iVar < nVar; ++iVar) {
      MachExtra_Left_Der[iVar]  = kP/(MeanSoundSpeed*MeanSoundSpeed*fa*fa*MeanDensity*MeanDensity)*
                                  ((sigma*factor*mL*Mach_Left_Der[iVar]*(Pressure_i - Pressure_j)*fa*MeanDensity) -
                                   (factor*S_i[iVar]*fa*MeanDensity) + (factor*(Pressure_j - Pressure_i)*MeanDensity*Scaling_Left_Der[iVar]));
      MachExtra_Right_Der[iVar] = kP/(MeanSoundSpeed*MeanSoundSpeed*fa*fa*MeanDensity*MeanDensity)*
                                 ((sigma*factor*mR*Mach_Right_Der[iVar]*(Pressure_i - Pressure_j)*fa*MeanDensity) +
                                  (factor*S_j[iVar]*fa*MeanDensity) + (factor*(Pressure_j - Pressure_i)*MeanDensity*Scaling_Right_Der[iVar]));
    }
    MachExtra_Left_Der[RHO_INDEX_SOL]  += kP/(MeanSoundSpeed*MeanSoundSpeed*fa*fa*MeanDensity*MeanDensity)*
                                          0.5*factor*(Pressure_j - Pressure_i)*fa;
    MachExtra_Right_Der[RHO_INDEX_SOL] += kP/(MeanSoundSpeed*MeanSoundSpeed*fa*fa*MeanDensity*MeanDensity)*
                                          0.5*factor*(Pressure_j - Pressure_i)*fa;

    /*---Set mass fluxes derivatives ---*/
    su2double MassPlus_Left_Der[nVar],  MassPlus_Right_Der[nVar],
              MassMinus_Left_Der[nVar], MassMinus_Right_Der[nVar];

    su2double sign_m12 = 0.0;
    if(m12 != 0.0)
      sign_m12 = std::abs(m12)/m12;
    for(iVar = 0; iVar < nVar; ++iVar) {
      MassPlus_Left_Der[iVar]    = 0.5*(MachPol_Left_Der[iVar]  - MachExtra_Left_Der[iVar])*(1.0 + sign_m12);
      MassMinus_Left_Der[iVar]   = 0.5*(MachPol_Left_Der[iVar]  - MachExtra_Left_Der[iVar])*(1.0 - sign_m12);
      MassPlus_Right_Der[iVar]   = 0.5*(MachPol_Right_Der[iVar] - MachExtra_Right_Der[iVar])*(1.0 + sign_m12);
      MassMinus_Right_Der[iVar]  = 0.5*(MachPol_Right_Der[iVar] - MachExtra_Right_Der[iVar])*(1.0 - sign_m12);
    }

    /*--- Set Jacobian for the convective part---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      for(jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] += MeanSoundSpeed*((MassPlus_Left_Der[jVar]*Density_i*Phi_i[iVar]) +
                                                      (MassMinus_Left_Der[jVar]*Density_j*Phi_j[iVar]));
        val_Jacobian_j[iVar][jVar] += MeanSoundSpeed*((MassPlus_Right_Der[jVar]*Density_i*Phi_i[iVar]) +
                                                      (MassMinus_Right_Der[jVar]*Density_j*Phi_j[iVar]));
      }
    }

    /*
    std::cout<<"first addition"<<std::endl;
    std::cout<<"val Jacobian i:"<<std::endl;
    for(iVar = 0; iVar < nVar; ++iVar) {
      for(jVar = 0; jVar< nVar; ++jVar)
        std::cout<<val_Jacobian_i[iVar][jVar]<<"    ";
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"val Jacobian j:"<<std::endl;
    for(iVar = 0; iVar < nVar; ++iVar) {
      for(jVar = 0; jVar< nVar; ++jVar)
        std::cout<<val_Jacobian_j[iVar][jVar]<<"    ";
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
    */

    /*--- Add pressure contribution to energy term ---*/
    su2double m12_P = 0.5*(m12 + std::abs(m12));
    su2double m12_M = 0.5*(m12 - std::abs(m12));
    for(iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[RHOE_INDEX_SOL][iVar] += MeanSoundSpeed*m12_P*S_i[iVar];
      val_Jacobian_j[RHOE_INDEX_SOL][iVar] += MeanSoundSpeed*m12_M*S_j[iVar];
    }

    /*--- Add contribution to the diagonal term ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[iVar][iVar] += MeanSoundSpeed*m12_P;
      val_Jacobian_j[iVar][iVar] += MeanSoundSpeed*m12_M;
    }

    /*--- Set polynomial pressure derivatives ---*/
    su2double PressPol_Left_Der[nVar], PressPol_Right_Der[nVar];

    if(std::abs(mL) < 1.0) {
      for(iVar = 0; iVar < nVar; ++iVar)
        PressPol_Left_Der[iVar] = 0.25*(mL + 1.0)*(3.0*(1.0 - mL) + 4.0*alpha*(5.0*mL*mL - 1.0)*(mL - 1.0))*Mach_Left_Der[iVar];
    }
    if(std::abs(mR) < 1.0) {
      for(iVar = 0; iVar < nVar; ++iVar)
        PressPol_Right_Der[iVar] = 0.25*(mL - 1.0)*(3.0*(1.0 + mL) + 4.0*alpha*(1.0 - 5.0*mL*mL)*(mL + 1.0))*Mach_Right_Der[iVar];
    }

    /*--- Set pressure extra term derivatives ---*/
    su2double PressExtra_Left_Der[nVar], PressExtra_Right_Der[nVar];
    for(iVar = 0; iVar < nVar; ++iVar) {
      PressExtra_Left_Der[iVar]  = Ku*pRM*MeanSoundSpeed*
                                   ((PressPol_Left_Der[iVar]*(Density_i + Density_j)*fa*(ProjVelocity_j - ProjVelocity_i)) +
                                    (pLP*(Density_i + Density_j)*(ProjVelocity_j - ProjVelocity_i)*Scaling_Left_Der[iVar]));
      PressExtra_Right_Der[iVar] = Ku*pLP*MeanSoundSpeed*
                                   ((PressPol_Right_Der[iVar]*(Density_i + Density_j)*fa*(ProjVelocity_j - ProjVelocity_i)) +
                                    (pRM*(Density_i + Density_j)*(ProjVelocity_j - ProjVelocity_i)*Scaling_Right_Der[iVar]));
    }
    PressExtra_Left_Der[RHO_INDEX_SOL]  += Ku*pRM*MeanSoundSpeed*pLP*fa*
                                           ((ProjVelocity_j - ProjVelocity_i) - (Density_i + Density_j)*ProjVelocity_i/Density_i);
    PressExtra_Right_Der[RHO_INDEX_SOL] += Ku*pLP*MeanSoundSpeed*pRM*fa*
                                           ((ProjVelocity_j - ProjVelocity_i) + (Density_i + Density_j)*ProjVelocity_j/Density_j);
    for(iDim = 0; iDim < nDim; ++iDim) {
      PressExtra_Left_Der[RHOVX_INDEX_SOL + iDim]  -= Ku*pRM*MeanSoundSpeed*pLP*(Density_i + Density_j)*fa*UnitNormal[iDim]/Density_i;
      PressExtra_Right_Der[RHOVX_INDEX_SOL + iDim] += Ku*pLP*MeanSoundSpeed*pRM*(Density_i + Density_j)*fa*UnitNormal[iDim]/Density_j;
    }

    /*--- Set whole pressure derivatives ---*/
    su2double Pressure_Left_Der[nVar], Pressure_Right_Der[nVar];
    for(iVar = 0; iVar < nVar; ++iVar) {
      Pressure_Left_Der[iVar]  = pLP*S_i[iVar] + Pressure_i*PressPol_Left_Der[iVar] - PressExtra_Left_Der[iVar];
      Pressure_Right_Der[iVar] = pRM*S_j[iVar] + Pressure_j*PressPol_Right_Der[iVar] - PressExtra_Right_Der[iVar];
    }

    /*--- Add global pressure contribution to the Jacobian ---*/
    for(iDim = 0; iDim < nDim; ++iDim) {
      for(jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[RHOVX_INDEX_SOL + iDim][jVar] += UnitNormal[iDim]*Pressure_Left_Der[jVar];
        val_Jacobian_j[RHOVX_INDEX_SOL + iDim][jVar] += UnitNormal[iDim]*Pressure_Right_Der[jVar];
      }
    }

    /*--- Integrate over dual-face area ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      for(jVar = 0; jVar < nVar; ++jVar) {
        val_Jacobian_i[iVar][jVar] *= Area;
        val_Jacobian_j[iVar][jVar] *= Area;
      }
    }

  }

  //AD::SetPreaccOut(val_residual, nVar);
  //AD::EndPreacc();
}

//
//
/*--- Constructor of the class CAvgGradReactive_Boundary ---*/
//
//
CAvgGradReactive_Boundary::CAvgGradReactive_Boundary(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr):
                           CNumerics(val_nDim, val_nVar, config), library(lib_ptr), nSpecies(library->GetnSpecies()) {
  /*--- Set local variables ---*/
  Laminar_Viscosity_i = Laminar_Viscosity_j = 0.0;
  Thermal_Conductivity_i = Thermal_Conductivity_j = 0.0;

  nPrimVar = nSpecies + nDim + 5;
  nPrimVarAvgGrad = nSpecies + nDim + 1;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set local variables to access indices in proper arrays ---*/
  T_INDEX_PRIM    = CReactiveNSVariable::GetT_INDEX_PRIM();
  VX_INDEX_PRIM   = CReactiveNSVariable::GetVX_INDEX_PRIM();
  P_INDEX_PRIM    = CReactiveNSVariable::GetP_INDEX_PRIM();
  RHO_INDEX_PRIM  = CReactiveNSVariable::GetRHO_INDEX_PRIM();
  H_INDEX_PRIM    = CReactiveNSVariable::GetH_INDEX_PRIM();
  A_INDEX_PRIM    = CReactiveNSVariable::GetT_INDEX_PRIM();
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

  T_INDEX_AVGGRAD    = 0;
  VX_INDEX_AVGGRAD   = 1;
  RHOS_INDEX_AVGGRAD = VX_INDEX_AVGGRAD + CReactiveNSVariable::GetnDim();

  PrimVar_i.resize(nPrimVar);
  PrimVar_j.resize(nPrimVar);

  Xs.resize(nSpecies);

  Mean_GradPrimVar.resize(nPrimVarAvgGrad, nDim);

  Gamma_tilde.resize(nSpecies,nSpecies);

  Ds_i.resize(nSpecies);
  Ds_j.resize(nSpecies);
  Ys_i.resize(nSpecies);
  Ys_j.resize(nSpecies);
}

//
//
/*--- Solution of Stefan-Maxwell equation using artificial diffusion modified matrix ---*/
//
//
void CAvgGradReactive_Boundary::Solve_SM(const su2double val_density, const su2double val_alpha, const RealMatrix& val_Dij,
                                         const RealVec& val_xs, const Vec& val_grad_xs, const RealVec& val_ys) {
  const su2double toll = 1e-11;

  /*--- Rename for convenience ---*/
  su2double rho = val_density;
  su2double alpha = val_alpha;

  /*--- Compute original matrix of Stefan-Maxwell equations ---*/
  Gamma = library->GetGamma(rho, val_xs, val_ys, val_Dij);

  /*--- Add artificial diffusion part ---*/
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
      Gamma_tilde(iSpecies,jSpecies) = Gamma(iSpecies,jSpecies) + alpha*val_ys[iSpecies];

  Eigen::BiCGSTAB<RealMatrix> bicg(Gamma_tilde);
  bicg.setTolerance(toll);
  Jd = bicg.solve(val_grad_xs);
}

//
//
/*--- Compute projection of viscous fluxes using Ramshaw self-consistent modification. ---*/
//
//
void CAvgGradReactive_Boundary::GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, su2double* val_normal,
                                                   const su2double val_viscosity, const su2double val_therm_conductivity,
                                                   const RealVec& val_diffusion_coeff, CConfig* config) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0) {
    SU2_Assert(Proj_Flux_Tensor != NULL, "The array for the projected viscous flux has not been allocated");
    SU2_Assert(Flux_Tensor != NULL, "The matrix for the viscous flux tensor has not been allocated");
    SU2_Assert(tau != NULL,"The matrix for the stress tensor has not been allocated");
  }

  /*--- Local variables ---*/
 	unsigned short iSpecies, iVar, iDim, jDim;

  su2double mu, ktr, div_vel;
  su2double rho, T;

  /*--- Initialize to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    Proj_Flux_Tensor[iVar] = 0.0;
    std::fill(Flux_Tensor[iVar], Flux_Tensor[iVar] + nDim, 0.0);
  }

  /*--- Rename for convenience ---*/
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  rho = val_primvar[RHO_INDEX_PRIM];
  T   = val_primvar[T_INDEX_PRIM];

  /*--- Compute partial enthalpies ---*/
  bool US_System = (config->GetSystemMeasurements() == US);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  hs = library->ComputePartialEnthalpy(dim_temp);
  for(auto& elem: hs)
    elem /= config->GetEnergy_Ref();
  if(US_System) {
    for(auto& elem: hs)
      elem *= 3.28084*3.28084;
  }

  /*--- Compute the velocity divergence ---*/
  div_vel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    div_vel += val_grad_primvar(VX_INDEX_AVGGRAD + iDim,iDim);

  /*--- Pre-compute mixture quantities ---*/
  su2double Normalization_Vec[nDim];
  for(iDim = 0; iDim < nDim; ++iDim) {
    Normalization_Vec[iDim] = 0.0;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Normalization_Vec[iDim] += val_diffusion_coeff[iSpecies]*val_grad_primvar(RHOS_INDEX_AVGGRAD + iSpecies,iDim);
  }

  /*--- Compute the viscous stress tensor ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    std::fill(tau[iDim], tau[iDim] + nDim, 0.0);

  for(iDim = 0; iDim < nDim; ++iDim) {
    for(jDim = 0; jDim < nDim; ++jDim)
      tau[iDim][jDim] += mu*(val_grad_primvar(VX_INDEX_AVGGRAD + jDim,iDim) + val_grad_primvar(VX_INDEX_AVGGRAD + iDim,jDim));
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*
  std::cout<<"val tau"<<std::endl;
  for(iDim = 0; iDim < nDim; ++iDim) {
    for(jDim = 0; jDim < nDim; ++jDim)
      std::cout<<tau[iDim][jDim]<<"     ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  */

  /*
  std::cout<<"val mu:   "<<mu<<std::endl;
  std::cout<<"val conductivity:   "<<ktr<<std::endl;
  std::cout<<"val gradprimvar"<<std::endl;
  for(iVar = 0; iVar < nPrimVarAvgGrad; ++iVar) {
    for(iDim = 0; iDim < nDim; ++iDim)
      std::cout<<val_grad_primvar(iVar,iDim)<<"       ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  */

  /*--- Populate entries in the viscous flux vector ---*/
  /*--- Density contribution ---*/
  for(iDim = 0; iDim < nDim; ++iDim) {
    /*--- Density contribution ---*/
    Flux_Tensor[RHO_INDEX_SOL][iDim] = 0.0;

    /*--- Shear stress related terms ---*/
    Flux_Tensor[RHOE_INDEX_SOL][iDim] = 0.0;
    for(jDim = 0; jDim < nDim; ++jDim) {
      Flux_Tensor[RHOVX_INDEX_SOL + jDim][iDim] = tau[iDim][jDim];
      Flux_Tensor[RHOE_INDEX_SOL][iDim] += tau[iDim][jDim]*val_primvar[VX_INDEX_PRIM + jDim];
    }

    /*--- Species diffusion velocity ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Flux_Tensor[RHOS_INDEX_SOL + iSpecies][iDim] = rho*
                                                     (val_diffusion_coeff[iSpecies]*val_grad_primvar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) -
                                                      val_primvar[RHOS_INDEX_PRIM + iSpecies]*Normalization_Vec[iDim]);

      /*--- Heat flux due to species diffusion term ---*/
      Flux_Tensor[RHOE_INDEX_SOL][iDim] += Flux_Tensor[RHOS_INDEX_SOL + iSpecies][iDim]*hs[iSpecies];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[RHOE_INDEX_SOL][iDim] += ktr*val_grad_primvar(T_INDEX_AVGGRAD,iDim);
  }

  for(iVar = 0; iVar < nVar; ++iVar) {
    for(iDim = 0; iDim < nDim; ++iDim)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
  }
}

//
//
/*--- Compute projection of viscous fluxes solving Stefan-Maxwell equations ---*/
//
//
void CAvgGradReactive_Boundary::GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, su2double* val_normal,
                                                   const su2double val_viscosity, const su2double val_thermal_conductivity,
                                                   const RealMatrix& val_Dij, CConfig* config) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0) {
    SU2_Assert(Proj_Flux_Tensor != NULL, "The array for the projected viscous flux has not been allocated");
    SU2_Assert(Flux_Tensor != NULL, "The matrix for the viscous flux tensor has not been allocated");
    SU2_Assert(tau != NULL,"The matrix for the stress tensor has not been allocated");
  }

  /*--- Local variables ---*/
 	unsigned short iSpecies, iVar, iDim, jDim;

  su2double mu, ktr, div_vel;
  su2double rho, T;

  /*--- Initialize to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    Proj_Flux_Tensor[iVar] = 0.0;
    std::fill(Flux_Tensor[iVar], Flux_Tensor[iVar] + nDim, 0.0);
  }

  /*--- Rename for convenience ---*/
  rho = val_primvar[RHO_INDEX_PRIM];
  mu  = val_viscosity;
  ktr = val_thermal_conductivity;
  T   = val_primvar[T_INDEX_PRIM];

  /*--- Compute partial enthalpies ---*/
  bool US_System = (config->GetSystemMeasurements() == US);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  hs = library->ComputePartialEnthalpy(dim_temp);
  for(auto& elem: hs)
    elem /= config->GetEnergy_Ref();
  if(US_System) {
    for(auto& elem: hs)
      elem *= 3.28084*3.28084;
  }

  /*-- NOTE: We need a non-standard primitive vector with molar fractions instead of partial densities ---*/
  /*--- Extract molar fractions, their gradient and mass fractions ---*/
  std::copy(val_primvar.data() + RHOS_INDEX_PRIM, val_primvar.data() + (RHOS_INDEX_PRIM + nSpecies), Xs.begin());
  Grad_Xs = val_grad_primvar.block(RHOS_INDEX_AVGGRAD, 0, nSpecies, nDim);
  Ys = library->GetMassFromMolar(Xs);

  /*--- Compute the velocity divergence ---*/
  div_vel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    div_vel += val_grad_primvar(VX_INDEX_AVGGRAD + iDim,iDim);

  /*--- Compute the viscous stress tensor ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    std::fill(tau[iDim], tau[iDim] + nDim, 0.0);

  for(iDim = 0; iDim < nDim; ++iDim) {
    for(jDim = 0; jDim < nDim; ++jDim)
      tau[iDim][jDim] += mu*(val_grad_primvar(VX_INDEX_AVGGRAD + jDim,iDim) + val_grad_primvar(VX_INDEX_AVGGRAD + iDim,jDim));
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux tensor ---*/
  const su2double alpha = 1.0/val_Dij.maxCoeff();
  const su2double sigma = std::accumulate(Ys.cbegin(), Ys.cend(), 0.0);
  const su2double sigma_alpha = sigma*alpha;
  for(iDim = 0; iDim < nDim; ++iDim) {
    /*--- Density contribution ---*/
    auto Grad_Xs_iDim = Grad_Xs.col(iDim);
    su2double Grad_sigma_iDim = Grad_Xs_iDim.sum();

    Flux_Tensor[RHO_INDEX_SOL][iDim] = rho*Grad_sigma_iDim/sigma_alpha;

    /*--- Shear stress related terms ---*/
    Flux_Tensor[RHOE_INDEX_SOL][iDim] = 0.0;
    for(jDim = 0; jDim < nDim; ++jDim) {
      Flux_Tensor[RHOVX_INDEX_SOL + jDim][iDim] = tau[iDim][jDim];
      Flux_Tensor[RHOE_INDEX_SOL][iDim] += tau[iDim][jDim]*val_primvar[VX_INDEX_PRIM + jDim];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[RHOE_INDEX_SOL][iDim] += ktr*val_grad_primvar(T_INDEX_AVGGRAD,iDim);

    /*--- Heat flux due to species diffusion term ---*/
    Solve_SM(rho, alpha, val_Dij, Xs, Grad_Xs_iDim, Ys);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Flux_Tensor[RHOS_INDEX_SOL + iSpecies][iDim] = Jd[iSpecies];
      Flux_Tensor[RHOE_INDEX_SOL][iDim] += Flux_Tensor[RHOS_INDEX_SOL + iSpecies][iDim]*hs[iSpecies];
    }
  }

  for(iVar = 0; iVar < nVar; ++iVar) {
    for(iDim = 0; iDim < nDim; ++iDim)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
  }
}

//
//
/*--- Compute approximate Jacobian of projected viscous flux. ---*/
//
//
void CAvgGradReactive_Boundary::GetViscousProjJacs(const Vec& val_Mean_PrimVar, const su2double val_laminar_viscosity,
                                                   const su2double val_thermal_conductivity, const Vec& val_diffusion_coeff,
                                                   const su2double val_dist_ij, const su2double val_dS, su2double* val_normal,
                                                   su2double* val_Proj_Visc_Flux, su2double** val_Proj_Jac_Tensor_i,
                                                   su2double** val_Proj_Jac_Tensor_j, CConfig* config) {
  /*--- Indexes for iteration ---*/
  unsigned short iDim, iVar, jVar, kVar, iSpecies, jSpecies;

  /*--- Local variables ---*/
  su2double theta = 0.0;
  su2double mu, ktr, dij, rho, rho_i, rho_j,T, dim_temp;
  su2double dFdVi[nVar][nVar], dFdVj[nVar][nVar], dVdUi[nVar][nVar], dVdUj[nVar][nVar];

  for(iDim = 0; iDim < nDim; ++iDim)
    theta += val_normal[iDim]*val_normal[iDim];

  /*--- Set Jacobian matrixes to zero ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    for(jVar = 0; jVar < nVar; ++jVar) {
      val_Proj_Jac_Tensor_i[iVar][jVar] = 0.0;
      val_Proj_Jac_Tensor_j[iVar][jVar] = 0.0;
      dFdVi[iVar][jVar] = 0.0;
      dFdVj[iVar][jVar] = 0.0;
    }
  }

  /*--- Rename for convenience ---*/
  dij = val_dist_ij;
  T = val_Mean_PrimVar[T_INDEX_PRIM];
  rho = val_Mean_PrimVar[RHO_INDEX_PRIM];
  rho_i = V_i[RHO_INDEX_PRIM];
  rho_j = V_j[RHO_INDEX_PRIM];
  mu = val_laminar_viscosity;
  ktr = val_thermal_conductivity;
  Ds = val_diffusion_coeff;

  /*--- Save mass fractions, partial enthalpies and specific heats at constant pressure---*/
  bool US_System = (config->GetSystemMeasurements() == US);
  dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  hs = library->ComputePartialEnthalpy(dim_temp);
  Cps = library->ComputeCps(dim_temp);
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Ys_i[iSpecies] = V_i[RHOS_INDEX_PRIM + iSpecies];
    Ys_j[iSpecies] = V_j[RHOS_INDEX_PRIM + iSpecies];
    hs[iSpecies]  /= config->GetEnergy_Ref();
    Cps[iSpecies] /= config->GetGas_Constant_Ref();
  }
  if(US_System) {
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      hs[iSpecies] *= 3.28084*3.28084;
      Cps[iSpecies] *= 3.28084*3.28084*5.0/9.0;
    }
  }

  /*--- Calculate useful diffusion parameters ---*/
  //Summation term of the diffusion fluxes
  su2double sumY_i = 0.0;
  su2double sumY_j = 0.0;
  su2double sum_Xi_Yj = 0.0;
  su2double sum_Xj_Yi = 0.0;
  su2double sigma_i = 0.0;
  su2double sigma_j = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    sumY_i += Ds[iSpecies]*Ys_i[iSpecies];
    sumY_j += Ds[iSpecies]*Ys_j[iSpecies];
    sum_Xi_Yj += Ds[iSpecies]*Xs_i[iSpecies]*Ys_j[iSpecies]/Xs_j[iSpecies];
    sum_Xj_Yi += Ds[iSpecies]*Xs_j[iSpecies]*Ys_i[iSpecies]/Xs_i[iSpecies];
    sigma_i += Ys_i[iSpecies];
    sigma_j += Ys_j[iSpecies];
  }

  /*--- Diffusion Jacobian matrices ---*/
  su2double dJdr_j[nSpecies][nSpecies], dJdr_i[nSpecies][nSpecies];
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    for(jSpecies  = 0; jSpecies < nSpecies; ++jSpecies) {
      /*--- Set to zero ---*/
      dJdr_i[iSpecies][jSpecies] = 0.0;
      dJdr_j[iSpecies][jSpecies] = 0.0;

      /*--- Compute first term ---*/
      dJdr_i[iSpecies][jSpecies] += -0.125*Ds[iSpecies]*theta*(1.0 + (rho_j/rho_i))*
                                     (Ys_i[iSpecies] + Ys_j[iSpecies]*Xs_i[iSpecies]*sigma_j/(sigma_i*Xs_j[iSpecies]));
      dJdr_j[iSpecies][jSpecies] +=  0.125*Ds[iSpecies]*theta*(1.0 + (rho_i/rho_j))*
                                     (Ys_j[iSpecies] + Ys_i[iSpecies]*Xs_j[iSpecies]*sigma_i/(sigma_j*Xs_i[iSpecies]));

      /*--- Compute second term ---*/
      dJdr_i[iSpecies][jSpecies] += 1.0/16*theta*(1.0 + (rho_j/rho_i))*(Ys_i[iSpecies] + Ys_j[iSpecies])*
                                    ((sumY_i + sigma_j*sum_Xi_Yj/sigma_i) +
                                     Ds[jSpecies]*(sigma_i + sigma_j*Ys_j[jSpecies]*Xs_i[jSpecies]/(Xs_j[iSpecies]*Ys_i[iSpecies])));
      dJdr_j[iSpecies][jSpecies] += -1.0/16*theta*(1.0 + (rho_i/rho_j))*(Ys_i[iSpecies] + Ys_j[iSpecies])*
                                    ((sumY_j + sigma_i*sum_Xj_Yi/sigma_j) +
                                     Ds[jSpecies]*(sigma_j + sigma_i*Ys_i[jSpecies]*Xs_j[jSpecies]/(Xs_i[iSpecies]*Ys_j[iSpecies])));
    }

    /*--- Add extra part first term ---*/
    dJdr_i[iSpecies][iSpecies] += -0.125*Ds[iSpecies]*theta*(1.0 + (rho_j/rho_i))*
                                  (sigma_i + sigma_j*Ys_j[iSpecies]*Xs_i[iSpecies]/(Ys_i[iSpecies]*Xs_j[iSpecies]));
    dJdr_j[iSpecies][iSpecies] += 0.125*Ds[iSpecies]*theta*(1.0 + (rho_i/rho_j))*
                                  (sigma_j + sigma_i*Ys_i[iSpecies]*Xs_j[iSpecies]/(Ys_j[iSpecies]*Xs_i[iSpecies]));
  }

  /*--- Calculate transformation matrix ---*/
  //Mixture density;
  dVdUi[RHO_INDEX_SOL][RHO_INDEX_SOL] = 1.0;
  dVdUj[RHO_INDEX_SOL][RHO_INDEX_SOL] = 1.0;
  // Partial densities
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    dVdUi[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + iSpecies] = 1.0;
    dVdUj[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + iSpecies] = 1.0;
  }
  for(iDim = 0; iDim < nDim; ++iDim) {
    dVdUi[RHOVX_INDEX_SOL + iDim][RHO_INDEX_SOL] = -V_i[VX_INDEX_PRIM + iDim]/V_i[RHO_INDEX_PRIM];
    dVdUi[RHOVX_INDEX_SOL + iDim][RHOVX_INDEX_SOL + iDim] = 1.0/V_i[RHO_INDEX_PRIM];
    dVdUj[RHOVX_INDEX_SOL + iDim][RHO_INDEX_SOL] = -V_j[VX_INDEX_PRIM + iDim]/V_j[RHO_INDEX_PRIM];
    dVdUj[RHOVX_INDEX_SOL + iDim][RHOVX_INDEX_SOL + iDim] = 1.0/V_j[RHO_INDEX_PRIM];
  }
  for(iVar = 0; iVar < nVar; ++iVar) {
    dVdUi[RHOE_INDEX_SOL][iVar] = S_i[iVar];
    dVdUj[RHOE_INDEX_SOL][iVar] = S_j[iVar];
  }

  /*--- Compute Jacobian with respect to primitives ---*/
  if(nDim == 2) {
    su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
    su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;

    su2double etaz = val_normal[0]*val_normal[1]/3.0;

    su2double pix = val_Mean_PrimVar[VX_INDEX_PRIM]*thetax + val_Mean_PrimVar[VX_INDEX_PRIM + 1]*etaz;
    su2double piy = val_Mean_PrimVar[VX_INDEX_PRIM]*etaz   + val_Mean_PrimVar[VX_INDEX_PRIM + 1]*thetay;

    /*--- Populate primitive Jacobian ---*/
    // X-momentum
    dFdVj[RHOVX_INDEX_SOL][RHO_INDEX_SOL] = -mu*pix/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL][RHOVX_INDEX_SOL] = mu*thetax/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL][RHOVX_INDEX_SOL + 1] = mu*etaz/(rho*dij)*val_dS;

    // Y-momentum
    dFdVj[RHOVX_INDEX_SOL + 1][RHO_INDEX_SOL] = -mu*piy/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL + 1][RHOVX_INDEX_SOL] = mu*etaz/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL + 1][RHOVX_INDEX_SOL + 1] = mu*thetay/(rho*dij)*val_dS;

    // Energy
    dFdVj[RHOE_INDEX_SOL][RHO_INDEX_SOL] = (-mu/rho*(pix*val_Mean_PrimVar[VX_INDEX_PRIM] + piy*val_Mean_PrimVar[VX_INDEX_PRIM + 1]) +
                                             ktr*theta)/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = (mu/rho*pix + ktr*theta*val_Mean_PrimVar[VX_INDEX_PRIM])/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = (mu/rho*piy + ktr*theta*val_Mean_PrimVar[VX_INDEX_PRIM + 1])/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOE_INDEX_SOL] = ktr*theta/dij*val_dS;

  } /*--- End onf nDim = 2 ---*/
  else {
    su2double thetax = theta + val_normal[0]*val_normal[0]/3.0;
    su2double thetay = theta + val_normal[1]*val_normal[1]/3.0;
    su2double thetaz = theta + val_normal[2]*val_normal[2]/3.0;

    su2double etax = val_normal[1]*val_normal[2]/3.0;
    su2double etay = val_normal[0]*val_normal[2]/3.0;
    su2double etaz = val_normal[0]*val_normal[1]/3.0;

    su2double pix = val_Mean_PrimVar[VX_INDEX_PRIM]*thetax + val_Mean_PrimVar[VX_INDEX_PRIM + 1]*etaz   +
                    val_Mean_PrimVar[VX_INDEX_PRIM + 2]*etay;
    su2double piy = val_Mean_PrimVar[VX_INDEX_PRIM]*etaz   + val_Mean_PrimVar[VX_INDEX_PRIM + 1]*thetay +
                    val_Mean_PrimVar[VX_INDEX_PRIM + 2]*etax;
    su2double piz = val_Mean_PrimVar[VX_INDEX_PRIM]*etay   + val_Mean_PrimVar[VX_INDEX_PRIM + 1]*etax   +
                    val_Mean_PrimVar[VX_INDEX_PRIM + 2]*thetaz;

    /*--- Populate primitive Jacobian ---*/
    // X-momentum
    dFdVj[RHOVX_INDEX_SOL][RHO_INDEX_SOL] = -mu*pix/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL][RHOVX_INDEX_SOL] = mu*thetax/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL][RHOVX_INDEX_SOL + 1] = mu*etaz/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL][RHOVX_INDEX_SOL + 2] = mu*etay/(rho*dij)*val_dS;

    // Y-momentum
    dFdVj[RHOVX_INDEX_SOL + 1][RHO_INDEX_SOL] = -mu*piy/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL + 1][RHOVX_INDEX_SOL] = mu*etaz/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL + 1][RHOVX_INDEX_SOL + 1] = mu*thetay/(rho*dij)*val_dS;
    dFdVj[RHOVX_INDEX_SOL + 1][RHOVX_INDEX_SOL + 2] = mu*etax/(rho*dij)*val_dS;

    // Energy
    dFdVj[RHOE_INDEX_SOL][RHO_INDEX_SOL] = (-mu/rho*(pix*val_Mean_PrimVar[VX_INDEX_PRIM] + piy*val_Mean_PrimVar[VX_INDEX_PRIM + 1] +
                                                     piz*val_Mean_PrimVar[VX_INDEX_PRIM + 2]) +  ktr*theta)/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL] = (mu/rho*pix + ktr*theta*val_Mean_PrimVar[VX_INDEX_PRIM])/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = (mu/rho*piy + ktr*theta*val_Mean_PrimVar[VX_INDEX_PRIM + 1])/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + 1] = (mu/rho*piz + ktr*theta*val_Mean_PrimVar[VX_INDEX_PRIM + 2])/dij*val_dS;
    dFdVj[RHOE_INDEX_SOL][RHOE_INDEX_SOL] = ktr*theta/dij*val_dS;
  } /*--- End of nDim = 3 ---*/

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

  // Common terms
  for(iDim = 0; iDim < nDim; ++iDim) {
    dFdVi[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + iDim] += 0.5*val_Proj_Visc_Flux[RHOVX_INDEX_SOL + iDim];
    dFdVj[RHOE_INDEX_SOL][RHOVX_INDEX_SOL + iDim] += 0.5*val_Proj_Visc_Flux[RHOVX_INDEX_SOL + iDim];
  }
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    dFdVi[RHOE_INDEX_SOL][RHOE_INDEX_SOL] += 0.5*val_Proj_Visc_Flux[RHOS_INDEX_SOL + iSpecies]*Cps[iSpecies];
    dFdVj[RHOE_INDEX_SOL][RHOE_INDEX_SOL] += 0.5*val_Proj_Visc_Flux[RHOS_INDEX_SOL + iSpecies]*Cps[iSpecies];
  }

  // Unique terms
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
      dFdVj[RHOS_INDEX_SOL + iSpecies][RHO_INDEX_SOL] += (1.0 + rho_j/rho_i)*theta*
                                                         (0.25*Ds[iSpecies]*
                                                         (sigma_i*Ys_i[iSpecies] + sigma_j*Ys_j[iSpecies]*Xs_i[iSpecies]/Xs_j[iSpecies])
                                                          - 0.125*(Ys_i[iSpecies] + Ys_j[iSpecies])*(sigma_i*sumY_i + sigma_j*sum_Xi_Yj)
                                                         )*val_dS/dij;
      dFdVi[RHOS_INDEX_SOL + iSpecies][RHO_INDEX_SOL] += -(1.0 + rho_i/rho_j)*theta*
                                                         (0.25*Ds[iSpecies]*
                                                         (sigma_j*Ys_j[iSpecies] + sigma_i*Ys_i[iSpecies]*Xs_j[iSpecies]/Xs_i[iSpecies])
                                                          - 0.125*(Ys_i[iSpecies] + Ys_j[iSpecies])*(sigma_j*sumY_j + sigma_i*sum_Xj_Yi)
                                                         )*val_dS/dij;
      dFdVj[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + jSpecies] += -dJdr_j[iSpecies][jSpecies]*val_dS/dij;
      dFdVi[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + jSpecies] += -dJdr_i[iSpecies][jSpecies]*val_dS/dij;
      dFdVj[RHO_INDEX_SOL][RHOS_INDEX_SOL + jSpecies] += dFdVj[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + jSpecies];
      dFdVi[RHO_INDEX_SOL][RHOS_INDEX_SOL + jSpecies] += dFdVi[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + jSpecies];
      dFdVj[RHOE_INDEX_SOL][RHOS_INDEX_SOL + iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS/dij;
      dFdVi[RHOE_INDEX_SOL][RHOS_INDEX_SOL + iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS/dij;
    }
  }

  for(iVar = 0; iVar < nVar; ++iVar) {
    for(jVar = 0; jVar < nVar; ++jVar) {
      for(kVar = 0; kVar < nVar; ++kVar) {
        val_Proj_Jac_Tensor_i[iVar][jVar] += dFdVi[iVar][kVar]*dVdUi[kVar][jVar];
        val_Proj_Jac_Tensor_j[iVar][jVar] += dFdVj[iVar][kVar]*dVdUj[kVar][jVar];
      }
    }
  }

} /*--- End of function ---*/

//
//
/*--- Compute residual for the viscous part for boundary conditions----*/
//
//
void CAvgGradReactive_Boundary::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                                                su2double** val_Jacobian_j, CConfig* config) {
  /*--- Indexes for iterations ---*/
  unsigned short iVar, iDim, jDim, iSpecies;

  /*--- Local variables ---*/
  su2double Mean_Laminar_Viscosity, Mean_Thermal_Conductivity;

  /*--- Mean transport coefficients ---*/
  Mean_Laminar_Viscosity = 2.0/(1.0/Laminar_Viscosity_i + 1.0/Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 2.0/(1.0/Thermal_Conductivity_i + 1.0/Thermal_Conductivity_j);
  Dij_i = Eigen::Map<RealMatrix>(Diffusion_Coeff_i, nSpecies, nSpecies);
  Dij_j = Eigen::Map<RealMatrix>(Diffusion_Coeff_j, nSpecies, nSpecies);
  Mean_Dij = 2.0/(Dij_i.cwiseInverse() + Dij_j.cwiseInverse()).array();

  /*--- Set to NULL for memory deallocation (otherwise double free because of base class) ---*/
  Diffusion_Coeff_i = Diffusion_Coeff_j = NULL;

  /*--- Copy primitive varaibles ---*/
  for(iVar = 0; iVar < nPrimVar; ++iVar) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
  }

  /*-- Use Molar fractions instead of mass fractions ---*/
  Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.data() + RHOS_INDEX_PRIM, PrimVar_i.data() + (RHOS_INDEX_PRIM + nSpecies)));
  Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.data() + RHOS_INDEX_PRIM, PrimVar_j.data() + (RHOS_INDEX_PRIM + nSpecies)));
  std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.data() + RHOS_INDEX_PRIM);
  std::copy(Xs_j.cbegin(), Xs_j.cend(), PrimVar_j.data() + RHOS_INDEX_PRIM);

  /*--- Compute the mean ---*/
  Mean_PrimVar = 0.5*(PrimVar_i + PrimVar_j);

  /*--- Mean gradient approximation ---*/
  for(iDim = 0; iDim < nDim; ++iDim) {
    /*--- Temperature ---*/
    Mean_GradPrimVar(T_INDEX_AVGGRAD,iDim) = 0.5*(PrimVar_Grad_i[T_INDEX_GRAD][iDim] + PrimVar_Grad_j[T_INDEX_GRAD][iDim]);

    /*--- Velocities ---*/
    for(jDim = 0; jDim < nDim; ++jDim)
      Mean_GradPrimVar(VX_INDEX_AVGGRAD + jDim,iDim) = 0.5*(PrimVar_Grad_i[VX_INDEX_GRAD + jDim][iDim] +
                                                            PrimVar_Grad_j[VX_INDEX_GRAD + jDim][iDim]);

    /*--- Molar Fractions ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Mean_GradPrimVar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) = 0.5*(PrimVar_Grad_i[RHOS_INDEX_GRAD + iSpecies][iDim] +
                                                                  PrimVar_Grad_j[RHOS_INDEX_GRAD + iSpecies][iDim]);
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Thermal_Conductivity, Mean_Dij, config);

	/*--- Update viscous residual with the species projected flux ---*/
  std::copy(Proj_Flux_Tensor, Proj_Flux_Tensor + nVar, val_residual);

  /*--- Implicit part ---*/
  if(implicit) {
    su2double denom_i, denom_j;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      denom_i = denom_j = 0.0;
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        if(jSpecies !!= iSpecies) {
          denom_i += Xs_i[jSpecies]/Dij_i(iSpecies, jSpecies);
          denom_j += Xs_j[jSpecies]/Dij_j(iSpecies, jSpecies);
        }
        Ds_i[iSpecies] = (1.0 - Xs_i[iSpecies])/denom_i;
        Ds_j[iSpecies] = (1.0 - Xs_j[iSpecies])/denom_j;
      }
    }

    /*--- NOTE: In case of single species present a NaN appears for that species so we check and we correct ---*/
    auto it_nan_i = std::find_if(Ds_i.data(), Ds_i.data() + Ds_i.size(), [](su2double elem){return std::isnan(elem);});
    if(it_nan_i != (Ds_i.data() + Ds_i.size()))
      *it_nan_i = 0.0;
    auto it_nan_j = std::find_if(Ds_j.data(), Ds_j.data() + Ds_j.size(), [](su2double elem){return std::isnan(elem);});
    if(it_nan_j != (Ds_j.data() + Ds_j.size()))
      *it_nan_j = 0.0;

    /*--- Compute the average of diffusion coefficients ---*/
    Ds = 0.5*(Ds_i + Ds_j);

    /*--- Compute distance between i and j ---*/
    su2double dist_ij_2 = 0.0;
    for(iDim = 0; iDim < nDim; iDim++)
      dist_ij_2 += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

    /*--- Compute area and unit normal ---*/
    su2double Area = std::inner_product(Normal, Normal + nDim, Normal, 0.0);
    Area = std::sqrt(Area);

    su2double UnitNormal[nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      UnitNormal[iDim] = Normal[iDim]/Area;

    GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Thermal_Conductivity, Ds,
                       std::sqrt(dist_ij_2), Area, UnitNormal, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j, config);
  }
}

//
//
/*--- Constructor of the class CAvgGradReactive_Flow ---*/
//
//
CAvgGradReactive_Flow::CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr):
                       CAvgGradReactive_Boundary(val_nDim, val_nVar, config, lib_ptr) {
  /*--- Resize local vectors ---*/
  Edge_Vector.resize(nDim);
  Diff_PrimVar.resize(nPrimVarAvgGrad);

  limiter = config->GetViscous_Limiter_Flow();
}

//
//
/*--- Compute residual for the viscous part at interface between two nodes. ----*/
//
//
void CAvgGradReactive_Flow::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) {
  /*--- Check memory allocation ---*/
  if(config->GetExtIter() == 0) {
    SU2_Assert(V_i != NULL, "The array for the primitive variables at node i has not been allocated");
    SU2_Assert(V_j != NULL, "The array for the primitive variables at node j has not been allocated");

    SU2_Assert(Mean_GradPrimVar.rows() == nPrimVarAvgGrad, "The number of rows in the mean gradient is not correct");
    SU2_Assert(Mean_GradPrimVar.cols() == nDim, "The number of columns in the mean gradient is not correct");
  }

  /*--- Indexes for iterations ---*/
  unsigned short iVar, iDim, jDim, iSpecies;

  /*--- Local variables ---*/
  su2double Mean_Laminar_Viscosity, Mean_Thermal_Conductivity;

  /*--- Mean transport coefficients ---*/
  Mean_Laminar_Viscosity = 2.0/(1.0/Laminar_Viscosity_i + 1.0/Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 2.0/(1.0/Thermal_Conductivity_i + 1.0/Thermal_Conductivity_j);
  Dij_i = Eigen::Map<RealMatrix>(Diffusion_Coeff_i, nSpecies, nSpecies);
  Dij_j = Eigen::Map<RealMatrix>(Diffusion_Coeff_j, nSpecies, nSpecies);
  Mean_Dij = 2.0/(Dij_i.cwiseInverse() + Dij_j.cwiseInverse()).array();

  /*--- Set to NULL for memory deallocation (otherwise double free because of base class) ---*/
  Diffusion_Coeff_i = Diffusion_Coeff_j = NULL;

  /*--- Copy primitive varaibles ---*/
  for(iVar = 0; iVar < nPrimVar; ++iVar) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
  }

  /*-- Use molar fractions instead of mass fractions ---*/
  Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.data() + RHOS_INDEX_PRIM, PrimVar_i.data() + (RHOS_INDEX_PRIM + nSpecies)));
  Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.data() + RHOS_INDEX_PRIM, PrimVar_j.data() + (RHOS_INDEX_PRIM + nSpecies)));
  std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.data() + RHOS_INDEX_PRIM);
  std::copy(Xs_j.cbegin(), Xs_j.cend(), PrimVar_j.data() + RHOS_INDEX_PRIM);

  /*--- Compute the mean ---*/
  Mean_PrimVar = 0.5*(PrimVar_i + PrimVar_j);

  /*--- Compute the vector from i to j ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    Edge_Vector[iDim] = Coord_j[iDim] - Coord_i[iDim];

  /*--- Mean gradient approximation ---*/
  if(!limiter) {
    for(iDim = 0; iDim < nDim; ++iDim) {
      /*--- Temperature ---*/
      Mean_GradPrimVar(T_INDEX_AVGGRAD,iDim) = 0.5*(PrimVar_Grad_i[T_INDEX_GRAD][iDim] + PrimVar_Grad_j[T_INDEX_GRAD][iDim]);

      /*--- Velocities ---*/
      for(jDim = 0; jDim < nDim; ++jDim)
        Mean_GradPrimVar(VX_INDEX_AVGGRAD + jDim,iDim) = 0.5*(PrimVar_Grad_i[VX_INDEX_GRAD + jDim][iDim] +
                                                              PrimVar_Grad_j[VX_INDEX_GRAD + jDim][iDim]);

      /*--- Molar Fractions ---*/
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Mean_GradPrimVar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) = 0.5*(PrimVar_Grad_i[RHOS_INDEX_GRAD + iSpecies][iDim] +
                                                                    PrimVar_Grad_j[RHOS_INDEX_GRAD + iSpecies][iDim]);
    }
  }
  else {
    /*--- Check memory allocation ---*/
    if(config->GetExtIter() == 0) {
      SU2_Assert(PrimVar_Lim_i != NULL, "The array for the primitive variables at node i has not been allocated");
      SU2_Assert(PrimVar_Lim_j != NULL, "The array for the primitive variables at node j has not been allocated");
    }

    for(iDim = 0; iDim < nDim; ++iDim) {
      /*--- Temperature ---*/
      Mean_GradPrimVar(T_INDEX_AVGGRAD,iDim) = 0.5*(PrimVar_Grad_i[T_INDEX_GRAD][iDim]*PrimVar_Lim_i[T_INDEX_LIM] +
                                                    PrimVar_Grad_j[T_INDEX_GRAD][iDim]*PrimVar_Lim_j[T_INDEX_LIM]);
      /*--- Velocities ---*/
      for(jDim = 0; jDim < nDim; ++jDim)
        Mean_GradPrimVar(VX_INDEX_AVGGRAD + jDim,iDim) = 0.5*(PrimVar_Grad_i[VX_INDEX_GRAD + jDim][iDim]*PrimVar_Lim_i[VX_INDEX_LIM + jDim] +
                                                              PrimVar_Grad_j[VX_INDEX_GRAD + jDim][iDim]*PrimVar_Lim_j[VX_INDEX_LIM + jDim]);
      /*--- Molar Fractions ---*/
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Mean_GradPrimVar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) = 0.5*(PrimVar_Grad_i[RHOS_INDEX_GRAD + iSpecies][iDim] +
                                                                    PrimVar_Grad_j[RHOS_INDEX_GRAD + iSpecies][iDim]);
    }
  }

  Proj_Mean_GradPrimVar_Edge = Mean_GradPrimVar*Edge_Vector;
  su2double dist_ij_2 = std::inner_product(Edge_Vector.data(), Edge_Vector.data() + Edge_Vector.size(), Edge_Vector.data(), 0.0);
  if(dist_ij_2 > EPS) {
    Diff_PrimVar[T_INDEX_AVGGRAD] = PrimVar_j[T_INDEX_PRIM] - PrimVar_i[T_INDEX_PRIM];
    /*--- Difference of velocities ---*/
    for(iDim = 0; iDim < nDim; ++iDim)
      Diff_PrimVar[VX_INDEX_AVGGRAD + iDim] = PrimVar_j[VX_INDEX_PRIM + iDim] - PrimVar_i[VX_INDEX_PRIM + iDim];
    /*--- Difference of mole fractions ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Diff_PrimVar[RHOS_INDEX_AVGGRAD + iSpecies] = V_j[RHOS_INDEX_PRIM + iSpecies] - V_i[RHOS_INDEX_PRIM + iSpecies];

    Mean_GradPrimVar -= (Proj_Mean_GradPrimVar_Edge - Diff_PrimVar)*Edge_Vector.transpose()/dist_ij_2;
  }
  else
    throw std::runtime_error("Error: You are trying to compute flux between a node and itself");

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Laminar_Viscosity, Mean_Thermal_Conductivity, Mean_Dij, config);

	/*--- Update viscous residual with the species projected flux ---*/
  std::copy(Proj_Flux_Tensor, Proj_Flux_Tensor + nVar, val_residual);

	/*--- Implicit part ---*/
	if(implicit) {
    su2double denom_i, denom_j;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      denom_i = denom_j = 0.0;
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        if(jSpecies != iSpecies) {
          denom_i += Xs_i[jSpecies]/Dij_i(iSpecies, jSpecies);
          denom_j += Xs_j[jSpecies]/Dij_j(iSpecies, jSpecies);
        }
        Ds_i[iSpecies] = (1.0 - Xs_i[iSpecies])/denom_i;
        Ds_j[iSpecies] = (1.0 - Xs_j[iSpecies])/denom_j;
      }
    }

    /*--- NOTE: In case of single species present a NaN appears for that species so we check and we correct ---*/
    auto it_nan_i = std::find_if(Ds_i.data(), Ds_i.data() + Ds_i.size(), [](su2double elem){return std::isnan(elem);});
    if(it_nan_i != (Ds_i.data() + Ds_i.size()))
      *it_nan_i = 0.0;
    auto it_nan_j = std::find_if(Ds_j.data(), Ds_j.data() + Ds_j.size(), [](su2double elem){return std::isnan(elem);});
    if(it_nan_j != (Ds_j.data() + Ds_j.size()))
      *it_nan_j = 0.0;

    /*--- Compute the average of diffusion coefficients ---*/
    Ds = 0.5*(Ds_i + Ds_j);

    /*--- Compute area and unit normal ---*/
    su2double Area = std::inner_product(Normal, Normal + nDim, Normal, 0.0);
    Area = std::sqrt(Area);

    su2double UnitNormal[nDim];
    for(iDim = 0; iDim < nDim; ++iDim)
      UnitNormal[iDim] = Normal[iDim]/Area;

    GetViscousProjJacs(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Thermal_Conductivity, Ds,
                       std::sqrt(dist_ij_2), Area, UnitNormal, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j, config);
  }
}

//
//
/*--- Constructor for the class CSourceChemistry ---*/
//
//
CSourceReactive::CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, CConfig* config, LibraryPtr lib_ptr):
                 CNumerics(val_nDim, val_nVar, config), library(lib_ptr), nSpecies(library->GetnSpecies()) {
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set local variables to access indices in proper arrays ---*/
  T_INDEX_PRIM    = CReactiveEulerVariable::GetT_INDEX_PRIM();
  VX_INDEX_PRIM   = CReactiveEulerVariable::GetVX_INDEX_PRIM();
  P_INDEX_PRIM    = CReactiveEulerVariable::GetP_INDEX_PRIM();
  RHO_INDEX_PRIM  = CReactiveEulerVariable::GetRHO_INDEX_PRIM();
  H_INDEX_PRIM    = CReactiveEulerVariable::GetH_INDEX_PRIM();
  A_INDEX_PRIM    = CReactiveEulerVariable::GetT_INDEX_PRIM();
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

  /*--- Resize local vectors ---*/
  Ys.resize(nSpecies);
}

//
//
/*--- Compute residual for chemistry source term ---*/
//
//
void CSourceReactive::ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) {
  /*--- Memory allocation check ---*/
  if(config->GetExtIter() == 0) {
    SU2_Assert(val_residual != NULL,"The array for residuals has not been allocated");
    SU2_Assert(V_i != NULL,"The array of primitive variables has not been allocated");
  }

  /*--- Local varaibles ---*/
  unsigned short iSpecies;

  su2double rho, temp, dim_temp, dim_rho;

  /*--- Nonequilibrium chemistry source term from library ---*/
  std::copy(V_i + RHOS_INDEX_PRIM, V_i + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  rho = V_i[RHO_INDEX_PRIM];
  temp = V_i[T_INDEX_PRIM];
  dim_temp = temp*config->GetTemperature_Ref();
  dim_rho = rho*config->GetDensity_Ref();
  bool US_System = (config->GetSystemMeasurements() == US);
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_rho *= 3.28084*3.28084*3.28084/0.0685218;
  }

  omega = library->GetMassProductionTerm(dim_temp, dim_rho, Ys);

  /*--- Assign to the residual. NOTE: We need to invert the sign since it is a residual that will be ADDED to the total one ---*/
  std::copy(omega.cbegin(), omega.cend(), val_residual + RHOS_INDEX_SOL);
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    val_residual[RHOS_INDEX_SOL + iSpecies] *= -Volume/(config->GetDensity_Ref()/config->GetTime_Ref());
  if(US_System) {
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      val_residual[RHOS_INDEX_SOL + iSpecies] /= 3.28084*3.28084*3.28084/0.0685218;
  }

  /*--- Set to zero the source residual for the other equations ---*/
  std::fill(val_residual, val_residual + RHOS_INDEX_SOL, 0.0);

  /*--- Implicit computation ---*/
  if(implicit) {
    /*--- Indexes for iteration ---*/
    unsigned short iDim, iVar, jSpecies;

    /*--- Check memory allocation ---*/
    if(config->GetExtIter() == 0) {
      SU2_Assert(S_i != NULL,"The array of primitive variables derivatives has not been allocated");
      SU2_Assert(val_Jacobian_i != NULL,"The matrix for source term jacobian has not been allocated");
      for(iVar = 0; iVar < nVar; ++iVar)
        SU2_Assert(val_Jacobian_i[iVar] != NULL,
                   std::string("The row " + std::to_string(iVar) + " of source chemistry jacobian has not been allocated"));
    }

    /*--- No source from density, momentum and total energy ---*/
    for(iVar = 0; iVar < nVar; ++iVar) {
      val_Jacobian_i[RHO_INDEX_SOL][iVar] = 0.0;
      for(iDim = 0; iDim < nDim; ++iDim)
        val_Jacobian_i[RHOVX_INDEX_SOL + iDim][iVar] = 0.0;
      val_Jacobian_i[RHOE_INDEX_SOL][iVar] = 0.0;
    }

    /*--- Jacobian from partial densities equations.
          NOTE: We need to invert the sign because we want the derivatives of a residual that will be ADDED to the total one ---*/
    source_jac = library->GetSourceJacobian(dim_temp, dim_rho);
    su2double fixed;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      fixed = source_jac(iSpecies,0)*config->GetTime_Ref()*config->GetTemperature_Ref()/config->GetDensity_Ref();
      if(US_System)
        fixed *= (5.0/9.0)*(0.0685218/(3.28084*3.28084*3.28084));
      val_Jacobian_i[RHOS_INDEX_SOL + iSpecies][RHO_INDEX_SOL] = -fixed*S_i[RHO_INDEX_SOL]*Volume;
      for(iDim = 0; iDim < nDim; ++iDim)
        val_Jacobian_i[RHOS_INDEX_SOL + iSpecies][RHOVX_INDEX_SOL + iDim] = -fixed*S_i[RHOVX_INDEX_SOL + iDim]*Volume;
      val_Jacobian_i[RHOS_INDEX_SOL + iSpecies][RHOE_INDEX_SOL] = -fixed*S_i[RHOE_INDEX_SOL]*Volume;
      for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
        val_Jacobian_i[RHOS_INDEX_SOL + iSpecies][RHOS_INDEX_SOL + jSpecies] =
        -fixed*S_i[RHOS_INDEX_SOL + jSpecies]*Volume - source_jac(iSpecies,jSpecies + 1)*config->GetTime_Ref()*Volume;
    }
  } /*--- End of implicit computation ---*/
}
