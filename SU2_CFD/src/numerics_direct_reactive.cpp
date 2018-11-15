#include "../include/numerics_reactive.hpp"
#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/move_pointer.hpp"
#include "../../Common/include/su2_assert.hpp"

#include <algorithm>

namespace {

  /*!
   * \brief Compute unit normal for the current cell interface
   */
  void Compute_Outward_UnitNormal(unsigned short nDim,su2double* Normal,su2double*& UnitNormal) {
    SU2_Assert(Normal != NULL,"Thw array for Normal has not been allocated");
    SU2_Assert(UnitNormal != NULL,"Thw array for Unit Normal has not been allocated");
    /*--- Face area (norm or the normal vector) ---*/
    su2double Area = std::inner_product(Normal,Normal + nDim, Normal, 0.0);
    Area = std::sqrt(Area);
    /*-- Unit Normal ---*/
    for (unsigned short iDim = 0; iDim < nDim; ++iDim)
      UnitNormal[iDim] = Normal[iDim]/Area;
    //UnitNormal = Normal/Area;
  }

} /*-- End of unnamed namespace ---*/

//
//
/*!
 * \brief Constructor of the class CUpwReactiveAUSM
 */
CUpwReactiveAUSM::CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, std::shared_ptr<CConfig> config):
    CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
    nSpecies(library->GetNSpecies()) {

  implicit = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;

}

//
//
/*!
 * \brief Compute residual convective term
 */
//
//

void CUpwReactiveAUSM::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                                       su2double** val_Jacobian_j, CConfig* config) {


  //AD::StartPreacc();
  //AD::SetPreaccIn(V_i, nSpecies + nDim+2);
  //AD::SetPreaccIn(V_j, nSpecies + nDim+2);
  //AD::SetPreaccIn(Normal, nDim);

  SU2_Assert(val_residual != NULL,"The array of residual for convective flux has not been allocated");
  //SU2_Assert(U_i != NULL,"The array of conserved variables at node i has not been allocated");
  //SU2_Assert(U_j != NULL,"The array of conserved variables at node j has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node i has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node j has not been allocated");

  unsigned short iDim, iVar, iSpecies; /*!< \brief Indexes for iterations. */
  su2double sq_vel_i,sq_vel_j,  /*!< \brief squared velocity. */
            //Energy_i, Energy_j, /*!< \brief Energy at node i and at node j. */
            Temperature_i, Temperature_j, /*!< \brief Temperature at node i and at node j. */
	          ProjVelocity_i, ProjVelocity_j; /*!< \brief Projected velocities at node i and at node j. */

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  /*--- Point i: compute energy,pressure,sound speed and enthalpy  ---*/
	sq_vel_i = std::inner_product(V_i + CReactiveEulerVariable::VX_INDEX_PRIM,V_i + CReactiveEulerVariable::VX_INDEX_PRIM + nDim,
                                V_i + CReactiveEulerVariable::VX_INDEX_PRIM,0.0);

  Density_i = V_i[CReactiveEulerVariable::RHO_INDEX_SOL];
  Pressure_i = V_i[CReactiveEulerVariable::P_INDEX_PRIM];
  //Energy_i = U_i[CReactiveEulerVariable::RHOE_INDEX_SOL]/Density_i;
  //Enthalpy_i = Energy_i + Pressure_i/Density_i + 0.5*sq_vel_i;
  Enthalpy_i = V_i[CReactiveEulerVariable::H_INDEX_PRIM];
  Temperature_i = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  library->Gamma_FrozenSoundSpeed(Temperature_i,Pressure_i,Density_i,Gamma,SoundSpeed_i);

  /*--- Point j: compute squared velocity,energy,pressure,sound speed and ethalpy  ---*/
  sq_vel_j = std::inner_product(V_j + CReactiveEulerVariable::VX_INDEX_PRIM,V_j + CReactiveEulerVariable::VX_INDEX_PRIM + nDim,
                                V_j + CReactiveEulerVariable::VX_INDEX_PRIM,0.0);

  Density_j = V_j[CReactiveEulerVariable::RHO_INDEX_SOL];
  Pressure_j = V_j[CReactiveEulerVariable::P_INDEX_PRIM];
  //Energy_j = U_j[CReactiveEulerVariable::RHOE_INDEX_SOL]/Density_j;
  //Enthalpy_j = Energy_j + Pressure_j/Density_j + 0.5*sq_vel_j;
  Enthalpy_j = V_j[CReactiveEulerVariable::H_INDEX_PRIM];
  Temperature_i = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  library->Gamma_FrozenSoundSpeed(Temperature_j,Pressure_j,Density_j,Gamma,SoundSpeed_j);

  /*--- Projected velocities ---*/
  ProjVelocity_i = std::inner_product(V_i + CReactiveEulerVariable::VX_INDEX_PRIM,V_i + CReactiveEulerVariable::VX_INDEX_PRIM + nDim,
                                      UnitNormal,0.0);
  ProjVelocity_j = std::inner_product(V_j + CReactiveEulerVariable::VX_INDEX_PRIM,V_j + CReactiveEulerVariable::VX_INDEX_PRIM + nDim,
                                      UnitNormal,0.0);

  /*--- Calculate L/R Mach numbers ---*/
  su2double mL = std::sqrt(sq_vel_i)/SoundSpeed_i;
  su2double mR = std::sqrt(sq_vel_j)/SoundSpeed_j;
  /*--- Calculate Mean sound speed ---*/
  su2double MeanSoundSpeed = 0.5*(SoundSpeed_i+SoundSpeed_j);
  /*--- Calculate adjusted polynomial function AUSM +-Up ---*/
  su2double mRef,mF;
  if(mL>=1.0 || mR >= 1.0)
    mRef = 1.0;
  else {
    mF = 0.5*(mL+mR);
    mRef = std::min(1.0,std::max(mF*mF,1e-4));
  }

  mF = std::sqrt((sq_vel_i + sq_vel_j)/(2.0*MeanSoundSpeed*MeanSoundSpeed));
  mRef = std::min(1.0,std::max(mF,mRef));
  const su2double fa = mRef*(2.0-mRef);
  const su2double alpha = 3.0/16.0*(5.0*fa*fa-4.0);

  /*--- Calculate Normal L/R Mach numbers ---*/
  mL  = ProjVelocity_i/MeanSoundSpeed;
  mR  = ProjVelocity_j/MeanSoundSpeed;

  su2double mLP,mRM,pLP,pRM;

  if (std::abs(mL) < 1.0) {
    mLP = 0.25*(mL+1.0)*(mL+1.0) + 0.125*(mL*mL-1.0)*(mL*mL-1.0);
    pLP = 0.25*(mL+1.0)*(mL+1.0)*(2.0-mL) + alpha*mL*(mL*mL-1.0);
  }
  else {
    mLP = 0.5*(mL+abs(mL));
    pLP = 0.5*(1.0+std::abs(mL)/mL);
  }

  if (std::abs(mR) < 1.0) {
    mRM = -0.25*(mR-1.0)*(mR-1.0) - 0.125*(mR*mR-1.0)*(mR*mR-1.0);
    pRM = 0.25*(mR-1.0)*(mR-1.0)*(mR+2.0) - alpha*mR*(mR*mR-1.0)*(mR*mR-1.0);
  }
  else {
    mRM = 0.5*(mR-abs(mR));
    pRM = 0.5*(1.0-std::abs(mR)/mR);
  }

  /*--- Build average state ---*/
  const su2double kP = 0.25;
  const su2double sigma = 1.0;

  mF = mLP + mRM;
  su2double mP = kP/fa*std::max(1.0-sigma*mF*mF,0.0)*(Pressure_i-Pressure_j)/(Pressure_i+Pressure_j);
  mF += mP;
  su2double mLF = 0.5*(mF+std::abs(mF));
  su2double mRF = 0.5*(mF-std::abs(mF));

  //su2double pF = pLP + pRM;
  su2double pLF = pLP*Pressure_i + pRM*Pressure_j;
  mF = MeanSoundSpeed*(mLF*Density_i + mRF*Density_j);
  //auto Phi = abs(mF);

  RealVec Ys_i,Ys_j;
  Ys_i.resize(nVar);
  Ys_j.resize(nVar);
  Ys_i[CReactiveEulerVariable::RHO_INDEX_SOL] = Density_i;
  Ys_j[CReactiveEulerVariable::RHO_INDEX_SOL] = Density_j;
  for(iDim = 0; iDim < nDim; ++iDim) {
    Ys_i[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = Density_i*V_i[CReactiveEulerVariable::VX_INDEX_PRIM + iDim];
    Ys_j[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = Density_j*V_j[CReactiveEulerVariable::VX_INDEX_PRIM + iDim];
  }
  Ys_i[CReactiveEulerVariable::RHOE_INDEX_SOL] = Density_i*Enthalpy_i;
  Ys_j[CReactiveEulerVariable::RHOE_INDEX_SOL] = Density_j*Enthalpy_j;
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Ys_i[CReactiveEulerVariable::RHOS_INDEX_SOL + iSpecies] = V_i[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies];
    Ys_j[CReactiveEulerVariable::RHOS_INDEX_SOL + iSpecies] = V_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies];
  }

  /*--- Calculate the numerical flux ---*/

  for (iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] = 0.5*(mF*(Ys_i[iVar]+Ys_j[iVar]) - std::abs(mF)*(Ys_i[iVar]-Ys_j[iVar]))*Area;

  const su2double Ku = 1.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    val_residual[CReactiveEulerVariable::VX_INDEX_PRIM+iDim] +=
    (pLF*UnitNormal[iDim] - Ku*pLP*pRM*(Density_i+Density_j)*mF*(ProjVelocity_i)*UnitNormal[iDim])*Area;
  }

  /*
  RealVec residual(val_residual, val_residual + nVar);
  residual = 0.5*(mF*(Ys_i + Ys_j)) - std::abs(mF)*(Ys_i - Ys_j)*Area;
  const su2double Ku = 1.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    residual[CReactiveEulerVariable::VX_INDEX_PRIM+iDim] +=
    (pLF*UnitNormal[iDim] - Ku*pLP*pRM*(Density_i+Density_j)*mF*(ProjVelocity_i)*UnitNormal[iDim])*Area;
  }

  val_residual = residual.data();
  */

  if(implicit) {
    SU2_Assert(val_Jacobian_i != NULL,"The jacobian at node i has not been allocated");
    SU2_Assert(val_Jacobian_j != NULL,"The jacobian at node j has not been allocated");

    /*--- Initialize the Jacobians ---*/
    for (iVar = 0; iVar < nVar; ++iVar) {
        std::fill(val_Jacobian_i[iVar],val_Jacobian_i[iVar] + nDim,0.0);
        std::fill(val_Jacobian_j[iVar],val_Jacobian_j[iVar] + nDim,0.0);
    }

    throw Common::NotImplemented("Implicit method for convective flux still not implemented. Setting explicit");

  }

  //AD::SetPreaccOut(val_residual, nVar);
  //AD::EndPreacc();
}

//
//
/*!
 * \brief Constructor of the class CAvgGradReactive_Flow
 */
//
//

CAvgGradReactive_Flow::CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, std::shared_ptr<CConfig> config):
  CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  nPrimVar = nSpecies + nDim + 5;
  nPrimVarAvgGrad = nSpecies + nDim + 1;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Mean_PrimVar.resize(nPrimVar);
  PrimVar_i.resize(nPrimVar);
  PrimVar_j.resize(nPrimVar);
  //Mean_GradPrimVar.resize(nPrimVarAvgGrad,RealVec(nDim));
  //GradPrimVar_i.resize(nPrimVarAvgGrad + 1,RealVec(nDim));
  //GradPrimVar_i.resize(nPrimVarAvgGrad + 1,RealVec(nDim));

  GradPrimVar_i.resize(nPrimVarAvgGrad + 1,nDim);
  GradPrimVar_j.resize(nPrimVarAvgGrad + 1,nDim);
  Mean_GradPrimVar.resize(nPrimVarAvgGrad,nDim);

}

//
//
/*!
 * \brief Compute projection of viscous fluxes
 */
//
//
void CAvgGradReactive_Flow::GetViscousProjFlux(RealVec& val_primvar, RealMatrix& val_grad_primvar,
                                               SmartArr val_normal, RealVec& val_diffusion_coeff,
                                               su2double val_viscosity,su2double val_therm_conductivity) {
  SU2_Assert(Proj_Flux_Tensor != NULL, "The array for the projected viscous flux has not been allocated");
  SU2_Assert(Flux_Tensor != NULL, "The matrix for the viscous flux tensor has not been allocated");
  SU2_Assert(tau != NULL,"The matrix for the stress tensor has not been allocated");
  // Requires a non-standard primitive vector with mass fractions instead of partial densities
 	unsigned short iSpecies, iVar, iDim, jDim;
  RealVec Ds, V;
  RealMatrix GV;
  su2double mu, ktr, div_vel;
  su2double rho, T, P;

  /*--- Initialize ---*/
  std::fill(Proj_Flux_Tensor,Proj_Flux_Tensor + nVar,0.0);
  for (iVar = 0; iVar < nVar; ++iVar)
    std::fill(Flux_Tensor[iVar],Flux_Tensor[iVar] + nDim, 0.0);

  /*--- Rename for convenience ---*/
  V   = val_primvar;
  GV  = val_grad_primvar;
  Ds  = val_diffusion_coeff;
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  rho = val_primvar[CReactiveNSVariable::RHO_INDEX_PRIM];
  //T   = val_primvar[CReactiveNSVariable::T_INDEX_PRIM];
  //P   = val_primvar[CReactiveNSVariable::P_INDEX_PRIM];
  RealVec hs(nSpecies);
  //library->GetSpeciesTotEnthalpies(T,P,hs);

  /*--- Calculate the velocity divergence ---*/
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; ++iDim)
    div_vel += GV(CReactiveNSVariable::VX_INDEX_GRAD + iDim,iDim);

  /*--- Pre-compute mixture quantities ---*/
  RealVec Vec(nDim);
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Vec[iDim] += rho*Ds[iSpecies]*GV(CReactiveNSVariable::RHOS_INDEX_AVGGRAD+iSpecies,iDim);
  }

  /*--- Compute the viscous stress tensor ---*/
  for (iDim = 0; iDim < nDim; ++iDim)
    std::fill(tau[iDim],tau[iDim] + nDim,0.0);


  for (iDim = 0 ; iDim < nDim; ++iDim) {
    for (jDim = 0 ; jDim < nDim; ++jDim)
      tau[iDim][jDim] += mu * (GV(CReactiveNSVariable::VX_INDEX_GRAD+jDim,iDim) + GV(CReactiveNSVariable::VX_INDEX_GRAD+iDim,jDim));
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux vector ---*/

  /*--- Density contribution ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {

    /*--- Density contribution ---*/
    Flux_Tensor[CReactiveNSVariable::RHO_INDEX_SOL][iDim] = 0.0;

    /*--- Shear stress related terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] = 0.0;
    for (jDim = 0; jDim < nDim; ++jDim) {
      Flux_Tensor[CReactiveNSVariable::RHOVX_INDEX_SOL + jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += tau[iDim][jDim]*V[CReactiveNSVariable::VX_INDEX_PRIM+jDim];
    }

    /*--- Species diffusion velocity ---*/
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim] = rho*Ds[iSpecies]*GV(CReactiveNSVariable::RHOS_INDEX_AVGGRAD+iSpecies,iDim)
                                                                          - V[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies]*Vec[iDim];

    /*--- Diffusion terms ---*/
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim] * hs[iSpecies];

    /*--- Heat transfer terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim]   += ktr*GV(CReactiveNSVariable::T_INDEX_GRAD,iDim);
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    for (iDim = 0; iDim < nDim; ++iDim)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
  }

}

//
//
/*!
 * \brief Compute projection of viscous fluxes
 */
//
//
void CAvgGradReactive_Flow::GetViscousProjJacs(RealVec& val_Mean_PrimVar, RealMatrix& val_Mean_GradPrimVar,
                                               RealVec& val_diffusion_coeff, su2double val_laminar_viscosity,su2double val_thermal_conductivity,
                                               su2double val_dist_ij, SmartArr val_normal, su2double val_dS,
                                               su2double* val_Proj_Visc_Flux,su2double** val_Proj_Jac_Tensor_i,su2double** val_Proj_Jac_Tensor_j,
                                               CConfig* config) {
  SU2_Assert(val_Proj_Jac_Tensor_i != NULL,"The matrix for projected jacobian at node i has not been allocated");
  SU2_Assert(val_Proj_Jac_Tensor_j != NULL,"The matrix for projected jacobian at node j has not been allocated");

  throw Common::NotImplemented("Calcutation of Jacobians has not been implemented");
}

//
//
/*!
 * \brief Compute residual for viscous term
 */
//
//

void CAvgGradReactive_Flow::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                                            su2double** val_Jacobian_j,CConfig* config) {
  SU2_Assert(V_i != NULL,"The array for the primitive variables at node i has not been allocated");
  SU2_Assert(V_j != NULL,"The array for the primitive variables at node j has not been allocated");

  SU2_Assert(GradPrimVar_i.nbRows() == nPrimVarAvgGrad + 1,"The number of rows in the gradient of varaible i is not correct");
  SU2_Assert(GradPrimVar_i.nbCols() == nDim,"The number of columns in the gradient of varaible i is not correct");
  SU2_Assert(GradPrimVar_j.nbRows() == nPrimVarAvgGrad + 1,"The number of rows in the gradient of varaible j is not correct");
  SU2_Assert(GradPrimVar_j.nbCols() == nDim,"The number of columns in the gradient of varaible j is not correct");
  SU2_Assert(Mean_GradPrimVar.nbRows() == nPrimVarAvgGrad,"The number of rows in the mean gradient is not correct");
  SU2_Assert(Mean_GradPrimVar.nbCols() == nDim,"The number of columns in the mean gradient is not correct");

  unsigned short iDim,jDim,iSpecies; /*!< \brief Indexes for iterations. */

  RealVec Diffusion_Coefficient_Species_i(nSpecies), /*!< \brief Diffusion coefficients species at node i. */
          Diffusion_Coefficient_Species_j(nSpecies),  /*!< \brief Diffusion coefficients species at node j. */
          Mean_Diffusion_Coefficient(nSpecies); /*!< \brief Mean value of diffusion coefficient. */

  su2double Mean_Laminar_Viscosity; /*!< \brief Mean value of laminar viscosity. */
  su2double Mean_Thermal_Conductivity;  /*!< \brief Mean value of thermal conductivity. */

  //Diffusion_Coefficient_Species_i = library->GetRhoUdiff(rho,kappa);
  //Diffusion_Coefficient_Species_i = library->GetRhoUdiff(rho,kappa);

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  //for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //  Mean_Diffusion_Coefficient[iSpecies] = 1.0/(1.0/Diffusion_Coefficient_Species_i[iSpecies] + 1.0/Diffusion_Coefficient_Species_j[iSpecies]);

  /*--- Mean transport coefficients ---*/
  Mean_Diffusion_Coefficient = 1.0/(1.0/Diffusion_Coefficient_Species_i + 1.0/Diffusion_Coefficient_Species_j);
  Mean_Laminar_Viscosity = 1.0/(1.0/Laminar_Viscosity_i + 1.0/Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 1.0/(1.0/Thermal_Conductivity_i + 1.0/Thermal_Conductivity_j);

  /*--- Set the primitive variables and compute the mean ---*/
  std::copy(V_i,V_i + nPrimVar,PrimVar_i.begin());
  std::copy(V_j,V_j + nPrimVar,PrimVar_j.begin());
  Mean_PrimVar = 0.5*(PrimVar_i + PrimVar_j);
  //for (iVar = 0; iVar < nVar; iVar++)
  //  Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);


  /*--- Mean gradient approximation ---*/
  // Mass fraction
  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    PrimVar_i[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies] = V_i[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies] /
                                                                  V_i[CReactiveNSVariable::RHO_INDEX_PRIM];
    PrimVar_j[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies] = V_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] /
                                                                  V_j[CReactiveNSVariable::RHO_INDEX_PRIM];
    Mean_PrimVar[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies] = 0.5*(PrimVar_i[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies] +
                                                                       PrimVar_j[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies]);
    for (iDim = 0; iDim < nDim; ++iDim) {
      Mean_GradPrimVar(CReactiveNSVariable::RHOS_INDEX_AVGGRAD,iDim) = 0.5*(1.0/V_i[CReactiveNSVariable::RHO_INDEX_PRIM]*
                                                                    (GradPrimVar_i(CReactiveNSVariable::RHOS_INDEX_AVGGRAD + iSpecies,iDim) -
                                                                     PrimVar_i[CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies]*
                                                                     GradPrimVar_i(CReactiveNSVariable::RHO_INDEX_AVGGRAD,iDim)) +
                                                                     1.0/V_j[CReactiveNSVariable::RHO_INDEX_PRIM]*
                                                                     (GradPrimVar_j(CReactiveNSVariable::RHO_INDEX_AVGGRAD,iDim) -
                                                                      PrimVar_j[CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies]*
                                                                      GradPrimVar_j(CReactiveNSVariable::RHO_INDEX_AVGGRAD,iDim)));
    }
  }
  //Temperature
  for(iDim = 0; iDim < nDim; ++iDim)
    Mean_GradPrimVar(CReactiveNSVariable::T_INDEX_GRAD,iDim) = 0.5*(GradPrimVar_i(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim) +
                                                                    GradPrimVar_j(CReactiveNSVariable::T_INDEX_AVGGRAD,iDim));
  //Velocities
  for(iDim = 0; iDim < nDim; ++iDim) {
    for(jDim = 0; jDim < nDim; ++jDim)
    Mean_GradPrimVar(CReactiveNSVariable::VX_INDEX_GRAD + iDim,jDim) = 0.5*(GradPrimVar_i(CReactiveNSVariable::VX_INDEX_AVGGRAD + iDim,jDim) +
                                                                            GradPrimVar_j(CReactiveNSVariable::VX_INDEX_AVGGRAD + iDim,jDim));
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Common::wrap_in_unique(Normal),
                     Mean_Diffusion_Coefficient, Mean_Laminar_Viscosity,Mean_Thermal_Conductivity);


	/*--- Update viscous residual with the species projected flux ---*/
  std::copy(Proj_Flux_Tensor,Proj_Flux_Tensor + nVar,val_residual);

	/*--- Implicit part ---*/
	if(implicit) {

    throw Common::NotImplemented("Implicit method still not implemented. Setting explicit");

    RealVec Edge_Vector(nDim);

    /*--- Compute vector going from iPoint to jPoint ---*/
    for (iDim = 0; iDim < nDim; iDim++)
    	Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    su2double dist_ij = std::inner_product(Edge_Vector.cbegin(),Edge_Vector.cend(),Edge_Vector.cbegin(),0.0);
    dist_ij = std::sqrt(dist_ij);

    GetViscousProjJacs(Mean_PrimVar,Mean_GradPrimVar,Mean_Diffusion_Coefficient,Mean_Laminar_Viscosity,Mean_Thermal_Conductivity,
                       dist_ij, Common::wrap_in_unique(UnitNormal), Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j,config);
  }

}

//
//
/*!
 * \brief Constructor for the class CSourceChemistry
 */
//
//

CSourceReactive::CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, std::shared_ptr<CConfig> config):
  CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
  nSpecies(library->GetNSpecies()) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

}

//
//
/*!
 * \brief Compute residual for chemistru term
 */
//
//

void CSourceReactive::ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) {

  SU2_Assert(val_residual != NULL,"The array for residuals has not been allocated");
  /*--- Initialize residual array ---*/
  std::fill(val_residual, val_residual + nVar, 0.0);

  /*--- Nonequilibrium chemistry from library ---*/
  RealVec Ys(U_i + CReactiveEulerVariable::RHOS_INDEX_SOL,U_i + (CReactiveEulerVariable::RHOS_INDEX_SOL + nSpecies));
  //auto res = library->GetMassProductionTerm(temp,rho,Ys);
  //std::copy(res.cbegin(),res.cend(),val_residual);
  std::for_each(val_residual,val_residual + nDim,[=](su2double elem){return elem *= -Volume;});
  //L'UNICA COSA A CUI STARE ATTENTI SONO I SEGNI PERCHE' QUESTO E' UN RESIDUO QUINDI MI SA CHE DEVO INVERTIRLI

  /*--- Implicit computation ---*/
  if (implicit) {
    SU2_Assert(val_Jacobian_i != NULL,"The matrix for Jacobian at node i has not been allocated");
    /*--- Initialize Jacobian array for implicit computation ---*/
    for (unsigned short iVar = 0; iVar < nVar; ++iVar)
      std::fill(val_Jacobian_i[iVar],val_Jacobian_i[iVar] + nVar,0.0);

    throw Common::NotImplemented("Implicit method for source chemistry residual not yet implemented. Setting explicit");
  }

}
