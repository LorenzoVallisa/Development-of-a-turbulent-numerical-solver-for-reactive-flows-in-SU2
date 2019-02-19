#include "../include/numerics_reactive.hpp"
#include <Eigen/IterativeLinearSolvers>

#include "../../Common/include/reacting_model_library.hpp"
#include "../../Common/include/not_implemented_exception.hpp"
#include "../../Common/include/move_pointer.hpp"
#include "../../Common/include/su2_assert.hpp"
#include "../../Common/include/default_initialization.hpp"

#include <algorithm>
#include <iterator>

namespace {
  /*!
   * \brief Compute unit normal for the current cell interface
   */
  void Compute_Outward_UnitNormal(unsigned short nDim, su2double* Normal, su2double* UnitNormal) {
    SU2_Assert(Normal != NULL,"The array for Normal has not been allocated");
    SU2_Assert(UnitNormal != NULL,"The array for Unit Normal has not been allocated");

    /*--- Face area (norm of the normal vector) ---*/
    su2double Area = std::inner_product(Normal,Normal + nDim, Normal, 0.0);
    Area = std::sqrt(Area);

    for(unsigned short iDim = 0; iDim < nDim; ++iDim)
      UnitNormal[iDim] = Normal[iDim]/Area;
  }
} /*-- End of unnamed namespace ---*/

//
//
/*!
 * \brief Constructor of the class CUpwReactiveAUSM
 */
CUpwReactiveAUSM::CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
                  CNumerics(val_nDim,val_nVar,config),library(CReactiveEulerVariable::GetLibrary()),nSpecies(library->GetNSpecies()) {
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Phi_i.resize(nVar);
  Phi_j.resize(nVar);
}

//
//
/*!
 * \brief Compute residual convective term
 */
//
//
void CUpwReactiveAUSM::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) {
  //AD::StartPreacc();
  //AD::SetPreaccIn(V_i, nSpecies + nDim + 5);
  //AD::SetPreaccIn(V_j, nSpecies + nDim + 5);
  //AD::SetPreaccIn(Normal, nDim);

  SU2_Assert(val_residual != NULL,"The array of residual for convective flux has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node i has not been allocated");
  SU2_Assert(V_j != NULL,"The array of primitive variables at node j has not been allocated");

  unsigned short iDim, iVar, iSpecies; /*!< \brief Indexes for iterations. */
  su2double sq_vel_i, sq_vel_j,  /*!< \brief squared velocity. */
            Temperature_i, Temperature_j, /*!< \brief Temperature at node i and at node j. */
            Sound_Speed_i, Sound_Speed_j, /*!< \brief Sound speed at node i and at node j. */
	          ProjVelocity_i, ProjVelocity_j; /*!< \brief Projected velocities at node i and at node j. */

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  /*--- Point i: compute energy,pressure,sound speed and enthalpy  ---*/
	sq_vel_i = std::inner_product(V_i + CReactiveEulerVariable::VX_INDEX_PRIM, V_i + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                                V_i + CReactiveEulerVariable::VX_INDEX_PRIM,0.0);
  Density_i = V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];
  Pressure_i = V_i[CReactiveEulerVariable::P_INDEX_PRIM];
  Enthalpy_i = V_i[CReactiveEulerVariable::H_INDEX_PRIM];
  Temperature_i = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  Sound_Speed_i = V_i[CReactiveEulerVariable::A_INDEX_PRIM];

  /*--- Point j: compute squared velocity,energy,pressure,sound speed and ethalpy  ---*/
  sq_vel_j = std::inner_product(V_j + CReactiveEulerVariable::VX_INDEX_PRIM, V_j + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                                V_j + CReactiveEulerVariable::VX_INDEX_PRIM,0.0);
  Density_j = V_j[CReactiveEulerVariable::RHO_INDEX_PRIM];
  Pressure_j = V_j[CReactiveEulerVariable::P_INDEX_PRIM];
  Enthalpy_j = V_j[CReactiveEulerVariable::H_INDEX_PRIM];
  Temperature_j = V_j[CReactiveEulerVariable::T_INDEX_PRIM];
  Sound_Speed_j = V_j[CReactiveEulerVariable::A_INDEX_PRIM];

  /*--- Projected velocities ---*/
  ProjVelocity_i = std::inner_product(V_i + CReactiveEulerVariable::VX_INDEX_PRIM, V_i + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                                      UnitNormal,0.0);
  ProjVelocity_j = std::inner_product(V_j + CReactiveEulerVariable::VX_INDEX_PRIM, V_j + (CReactiveEulerVariable::VX_INDEX_PRIM + nDim),
                                      UnitNormal,0.0);

  /*
  if(grid_movement) {
    ProjVelocity_i -= std::inner_product(GridVel_i,GridVel_i + nDim,UnitNormal,0.0);
    ProjVelocity_j -= std::inner_product(GridVel_j,GridVel_j + nDim,UnitNormal,0.0);
  }
  */

  /*--- Calculate L/R Mach numbers ---*/
  su2double mL = std::sqrt(sq_vel_i)/SoundSpeed_i;
  su2double mR = std::sqrt(sq_vel_j)/SoundSpeed_j;

  /*--- Calculate Mean sound speed ---*/
  su2double MeanSoundSpeed = 0.5*(SoundSpeed_i + SoundSpeed_j);

  /*--- Calculate adjusted polynomial function AUSM +-Up ---*/
  su2double mRef,mF;
  if(mL >= 1.0 || mR >= 1.0)
    mRef = 1.0;
  else
    mRef = std::min(1.0,std::max(1e-4,(0.25*(mL+mR)*(mL+mR))));

  /*--- Calculate Normal L/R Mach numbers ---*/
  mL  = ProjVelocity_i/MeanSoundSpeed;
  mR  = ProjVelocity_j/MeanSoundSpeed;
  mF = std::sqrt(0.5*(mL*mL + mR*mR));
  mRef = std::min(1.0,std::max(mF,mRef));
  const su2double fa = mRef*(2.0 - mRef);
  const su2double alpha = 3.0/16.0*(5.0*fa*fa - 4.0);
  const double beta = 0.125;

  su2double mLP,mRM,pLP,pRM;

  if(std::abs(mL) < 1.0) {
    mLP = 0.25*(mL + 1.0)*(mL + 1.0) + beta*(mL*mL - 1.0)*(mL*mL - 1.0);
    pLP = 0.25*(mL + 1.0)*(mL + 1.0)*(2.0 - mL) + alpha*mL*(mL*mL - 1.0);
  }
  else {
    mLP = 0.5*(mL + abs(mL));
    pLP = 0.5*(1.0 + std::abs(mL)/mL);
  }

  if(std::abs(mR) < 1.0) {
    mRM = -0.25*(mR - 1.0)*(mR - 1.0) - beta*(mR*mR - 1.0)*(mR*mR - 1.0);
    pRM = 0.25*(mR - 1.0)*(mR - 1.0)*(mR + 2.0) - alpha*mR*(mR*mR - 1.0)*(mR*mR - 1.0);
  }
  else {
    mRM = 0.5*(mR - abs(mR));
    pRM = 0.5*(1.0 - std::abs(mR)/mR);
  }

  /*--- Build average state ---*/
  const su2double kP = 0.25;
  const su2double sigma = 1.0;

  su2double m12 = mLP + mRM;
  su2double mP = kP/fa*std::max(1.0 - sigma*mF*mF, 0.0)*(Pressure_i - Pressure_j)/(Pressure_i + Pressure_j);
  m12 += mP;
  su2double mLF = 0.5*(m12 + std::abs(m12));
  su2double mRF = 0.5*(m12 - std::abs(m12));

  su2double pLF = pLP*Pressure_i + pRM*Pressure_j;
  m12 = MeanSoundSpeed*(mLF*Density_i + mRF*Density_j);

  /*--- Compute the state at node i and at node j ---*/
  Phi_i[CReactiveEulerVariable::RHO_INDEX_SOL] = 1.0;
  Phi_j[CReactiveEulerVariable::RHO_INDEX_SOL] = 1.0;
  for(iDim = 0; iDim < nDim; ++iDim) {
    Phi_i[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = V_i[CReactiveEulerVariable::VX_INDEX_PRIM + iDim];
    Phi_j[CReactiveEulerVariable::RHOVX_INDEX_SOL + iDim] = V_j[CReactiveEulerVariable::VX_INDEX_PRIM + iDim];
  }
  Phi_i[CReactiveEulerVariable::RHOE_INDEX_SOL] = Enthalpy_i;
  Phi_j[CReactiveEulerVariable::RHOE_INDEX_SOL] = Enthalpy_j;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Phi_i[CReactiveEulerVariable::RHOS_INDEX_SOL + iSpecies] = V_i[CReactiveEulerVariable::RHOS_INDEX_PRIM + iSpecies];
    Phi_j[CReactiveEulerVariable::RHOS_INDEX_SOL + iSpecies] = V_j[CReactiveEulerVariable::RHOS_INDEX_PRIM + iSpecies];
  }

  /*--- Calculate the numerical flux ---*/
  for(iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] = 0.5*(m12*(Phi_i[iVar] + Phi_j[iVar]) - std::abs(m12)*(Phi_j[iVar] - Phi_i[iVar]))*Area;

  const su2double Ku = 0.75;
  for(iDim = 0; iDim < nDim; ++iDim) {
    val_residual[CReactiveEulerVariable::VX_INDEX_PRIM + iDim] +=
    (pLF*UnitNormal[iDim] - Ku*pLP*pRM*(Density_i + Density_j)*m12*ProjVelocity_i*UnitNormal[iDim])*Area;
  }

  if(implicit)
    throw Common::NotImplemented("Implicit method for convective flux still not implemented. Setting explicit");

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
CAvgGradReactive_Flow::CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
  CNumerics(val_nDim,val_nVar,config),library(CReactiveNSVariable::GetLibrary()),nSpecies(library->GetNSpecies()) {

  Laminar_Viscosity_i = Laminar_Viscosity_j = 0.0;
  Thermal_Conductivity_i = Thermal_Conductivity_j = 0.0;

  nPrimVar = nSpecies + nDim + 5;
  nPrimVarAvgGrad = nSpecies + nDim + 1;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  limiter = config->GetViscous_Limiter_Flow();

  Edge_Vector.resize(nDim);

  PrimVar_i.resize(nPrimVar);
  PrimVar_j.resize(nPrimVar);
  Diff_PrimVar.resize(nPrimVarAvgGrad);

  Xs.resize(nSpecies);

  Mean_GradPrimVar.resize(nPrimVarAvgGrad, nDim);

  Gamma.resize(nSpecies,nSpecies);
  Gamma_tilde.resize(nSpecies,nSpecies);
}

//
//
/*!
 * \brief Solution of Stefan-Maxwell equation using artificial diffusion modified matrix
 */
//
//
CAvgGradReactive_Flow::Vec CAvgGradReactive_Flow::Solve_SM(const su2double val_density, const su2double val_alpha, const RealMatrix& val_Dij,
                                                           const RealVec& val_xs, const Vec& val_grad_xs, const RealVec& val_ys) {
  su2double toll = 1e-11;

  /*--- Rename for convenience ---*/
  su2double rho = val_density;
  su2double alpha = val_alpha;

  /*--- Compute original matrix of Stefan-Maxwell equations ---*/
  //Gamma = library->GetGamma(rho,val_xs,val_ys,val_Dij)

  /*--- Add artificial diffusion part ---*/
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
      Gamma_tilde(iSpecies,jSpecies) = Gamma(iSpecies,jSpecies) + alpha*val_ys[iSpecies];

  Eigen::BiCGSTAB<RealMatrix> bicg(Gamma_tilde);
  bicg.setTolerance(toll);
  Vec Jd = bicg.solve(val_grad_xs);

  return Jd;
}

//
//
/*!
 * \brief Compute projection of viscous fluxes using Ramshaw self-consistent modification
 */
//
//
void CAvgGradReactive_Flow::GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, SmartArr val_normal,
                                               const su2double val_viscosity, const su2double val_therm_conductivity,
                                               const RealVec& val_diffusion_coeff, CConfig* config) {
  SU2_Assert(Proj_Flux_Tensor != NULL, "The array for the projected viscous flux has not been allocated");
  SU2_Assert(Flux_Tensor != NULL, "The matrix for the viscous flux tensor has not been allocated");
  SU2_Assert(tau != NULL,"The matrix for the stress tensor has not been allocated");

  /*--- We need a non-standard primitive vector with mass fractions instead of partial densities ---*/
 	unsigned short iSpecies, iVar, iDim, jDim;
  su2double mu, ktr, div_vel;
  su2double rho, T;

  /*--- Initialize ---*/
  std::fill(Proj_Flux_Tensor,Proj_Flux_Tensor + nVar,0.0);
  for(iVar = 0; iVar < nVar; ++iVar)
    std::fill(Flux_Tensor[iVar],Flux_Tensor[iVar] + nDim, 0.0);

  /*--- Rename for convenience ---*/
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  rho = val_primvar[CReactiveNSVariable::RHO_INDEX_PRIM];
  T   = val_primvar[CReactiveNSVariable::T_INDEX_PRIM];

  /*--- Compute partial enthalpies ---*/
  bool US_System = (config->GetSystemMeasurements() == SI);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //hs = library->ComputePartialEnthalpy(dim_temp);
  //std::transform(hs.begin(),hs.end(),hs.begin(),[config](su2double elem){return elem*config->GetEnergy_Ref();});
  //if(US_System)
  //  std::transform(hs.begin(),hs.end(),hs.begin(),[config](su2double elem){return elem*3.28084*3.28084;});

  /*--- Compute the velocity divergence ---*/
  div_vel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    div_vel += val_grad_primvar(VX_INDEX_AVGGRAD + iDim,iDim);

  /*--- Pre-compute mixture quantities ---*/
  RealVec Normalization_Vec(nDim);
  for(iDim = 0; iDim < nDim; ++iDim) {
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Normalization_Vec[iDim] += val_diffusion_coeff[iSpecies]*val_grad_primvar(RHOS_INDEX_AVGGRAD + iSpecies,iDim);
  }

  /*--- Compute the viscous stress tensor ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    std::fill(tau[iDim],tau[iDim] + nDim,0.0);

  for(iDim = 0; iDim < nDim; ++iDim) {
    for(jDim = 0; jDim < nDim; ++jDim)
      tau[iDim][jDim] += mu*(val_grad_primvar(VX_INDEX_AVGGRAD + jDim,iDim) + val_grad_primvar(VX_INDEX_AVGGRAD + iDim,jDim));
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux vector ---*/
  /*--- Density contribution ---*/
  for(iDim = 0; iDim < nDim; ++iDim) {
    /*--- Density contribution ---*/
    Flux_Tensor[CReactiveNSVariable::RHO_INDEX_SOL][iDim] = 0.0;

    /*--- Shear stress related terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] = 0.0;
    for(jDim = 0; jDim < nDim; ++jDim) {
      Flux_Tensor[CReactiveNSVariable::RHOVX_INDEX_SOL + jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += tau[iDim][jDim]*val_primvar[CReactiveNSVariable::VX_INDEX_PRIM+jDim];
    }

    /*--- Species diffusion velocity ---*/
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim] = rho*(val_diffusion_coeff[iSpecies]*
                                                                               val_grad_primvar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) -
                                                                               val_primvar[CReactiveNSVariable::RHOS_INDEX_PRIM+iSpecies]*
                                                                               Normalization_Vec[iDim]);
      /*--- Heat flux due to species diffusion term ---*/
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] +=
      Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim]*hs[iSpecies];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += ktr*val_grad_primvar(T_INDEX_AVGGRAD,iDim);
  }

  for(iVar = 0; iVar < nVar; ++iVar) {
    for(iDim = 0; iDim < nDim; ++iDim)
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
void CAvgGradReactive_Flow::GetViscousProjFlux(const Vec& val_primvar, const RealMatrix& val_grad_primvar, SmartArr val_normal,
                                               const double val_viscosity, const double val_thermal_conductivity,
                                               const RealMatrix& val_Dij, CConfig* config) {
  SU2_Assert(Proj_Flux_Tensor != NULL, "The array for the projected viscous flux has not been allocated");
  SU2_Assert(Flux_Tensor != NULL, "The matrix for the viscous flux tensor has not been allocated");
  SU2_Assert(tau != NULL,"The matrix for the stress tensor has not been allocated");

  /*--- We need a non-standard primitive vector with molar fractions instead of partial densities ---*/
 	unsigned short iSpecies, iVar, iDim, jDim;
  su2double mu, ktr, div_vel;
  su2double rho, T;

  /*--- Initialize ---*/
  std::fill(Proj_Flux_Tensor,Proj_Flux_Tensor + nVar,0.0);
  for(iVar = 0; iVar < nVar; ++iVar)
    std::fill(Flux_Tensor[iVar],Flux_Tensor[iVar] + nDim, 0.0);

  /*--- Rename for convenience ---*/
  rho = val_primvar[CReactiveNSVariable::RHO_INDEX_PRIM];
  mu  = val_viscosity;
  ktr = val_thermal_conductivity;
  T = val_primvar[CReactiveNSVariable::T_INDEX_PRIM];

  /*--- Compute partial enthalpies ---*/
  bool US_System = (config->GetSystemMeasurements() == SI);
  su2double dim_temp = T*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //hs = library->ComputePartialEnthalpy(dim_temp);
  //std::transform(hs.begin(),hs.end(),hs.begin(),[config](su2double elem){return elem*config->GetEnergy_Ref();});
  //if(US_System)
  //  std::transform(hs.begin(),hs.end(),hs.begin(),[config](su2double elem){return elem*3.28084*3.28084;});

  /*--- Extract molar fractions, their gradient and mass fractions ---*/
  std::copy(val_primvar.data() + CReactiveNSVariable::RHOS_INDEX_PRIM,
            val_primvar.data() + (CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies), Xs.begin());
  Grad_Xs = val_grad_primvar.block(RHOS_INDEX_AVGGRAD,0,nSpecies,nDim);
  //Ys = library->GetMassFromMolar(Xs);

  /*--- Compute the velocity divergence ---*/
  div_vel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    div_vel += val_grad_primvar(VX_INDEX_AVGGRAD + iDim,iDim);

  /*--- Compute the viscous stress tensor ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    std::fill(tau[iDim],tau[iDim] + nDim,0.0);

  for(iDim = 0 ; iDim < nDim; ++iDim) {
    for(jDim = 0 ; jDim < nDim; ++jDim)
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

    Flux_Tensor[CReactiveNSVariable::RHO_INDEX_SOL][iDim] = rho*Grad_sigma_iDim/sigma_alpha;

    /*--- Shear stress related terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] = 0.0;
    for(jDim = 0; jDim < nDim; ++jDim) {
      Flux_Tensor[CReactiveNSVariable::RHOVX_INDEX_SOL + jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += tau[iDim][jDim]*val_primvar[CReactiveNSVariable::VX_INDEX_PRIM + jDim];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] += ktr*val_grad_primvar(T_INDEX_AVGGRAD,iDim);

    /*--- Heat flux due to species diffusion term ---*/
    auto Jd_iDim = Solve_SM(rho, alpha, val_Dij, Xs, Grad_Xs_iDim, Ys);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim] = Jd_iDim[iSpecies];
      Flux_Tensor[CReactiveNSVariable::RHOE_INDEX_SOL][iDim] +=
      Flux_Tensor[CReactiveNSVariable::RHOS_INDEX_SOL + iSpecies][iDim]*hs[iSpecies];
    }
  }

  for(iVar = 0; iVar < nVar; ++iVar) {
    for(iDim = 0; iDim < nDim; ++iDim)
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
  }
}


//
//
/*!
 * \brief Compute Jacobian of projected viscous flux
 */
//
//
void CAvgGradReactive_Flow::GetViscousProjJacs(const Vec& val_Mean_PrimVar, const RealMatrix& val_Mean_GradPrimVar,
                                               const RealVec& val_diffusion_coeff, const su2double val_laminar_viscosity,
                                               const su2double val_thermal_conductivity, const su2double val_dist_ij, SmartArr val_normal,
                                               const su2double val_dS, su2double* val_Proj_Visc_Flux, su2double** val_Proj_Jac_Tensor_i,
                                               su2double** val_Proj_Jac_Tensor_j, CConfig* config) {

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
                                            su2double** val_Jacobian_j, CConfig* config) {
  SU2_Assert(V_i != NULL, "The array for the primitive variables at node i has not been allocated");
  SU2_Assert(V_j != NULL, "The array for the primitive variables at node j has not been allocated");
  SU2_Assert(PrimVar_Lim_i != NULL, "The array for the primitive variables at node i has not been allocated");
  SU2_Assert(PrimVar_Lim_j != NULL, "The array for the primitive variables at node j has not been allocated");

  SU2_Assert(Mean_GradPrimVar.rows() == nPrimVarAvgGrad, "The number of rows in the mean gradient is not correct");
  SU2_Assert(Mean_GradPrimVar.cols() == nDim, "The number of columns in the mean gradient is not correct");

  unsigned short iVar, iDim, jDim, iSpecies; /*!< \brief Indexes for iterations. */

  su2double Mean_Laminar_Viscosity; /*!< \brief Mean value of laminar viscosity. */
  su2double Mean_Thermal_Conductivity;  /*!< \brief Mean value of thermal conductivity. */

  /*--- Mean transport coefficients ---*/
  Mean_Laminar_Viscosity = 2.0/(1.0/Laminar_Viscosity_i + 1.0/Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 2.0/(1.0/Thermal_Conductivity_i + 1.0/Thermal_Conductivity_j);
  Mean_Dij = 2.0/(Dij_i.cwiseInverse() + Dij_j.cwiseInverse()).array();

  /*--- Copy primitive varaibles ---*/
  std::copy(V_i, V_i + nPrimVar, PrimVar_i.data());
  std::copy(V_j, V_j + nPrimVar, PrimVar_j.data());

  /*-- Use Molar fractions instead of mass fractions ---*/
  //Xs_i = library->GetMolarFromMass(RealVec(PrimVar_i.data() + CReactiveNSVariable::RHOS_INDEX_PRIM,
  //                                         PrimVar_i.data() + (CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies));
  //Xs_j = library->GetMolarFromMass(RealVec(PrimVar_j.data() + CReactiveNSVariable::RHOS_INDEX_PRIM,
  //                                         PrimVar_i.data() + (CReactiveNSVariable::RHOS_INDEX_PRIM + nSpecies));
  std::copy(Xs_i.cbegin(), Xs_i.cend(), PrimVar_i.data() + CReactiveNSVariable::RHOS_INDEX_PRIM);
  std::copy(Xs_j.cbegin(), Xs_j.cend(), PrimVar_j.data() + CReactiveNSVariable::RHOS_INDEX_PRIM);

  /*--- Compute the mean ---*/
  Mean_PrimVar = 0.5*(PrimVar_i + PrimVar_j);

  /*
  if(grid_movement)
    for(iDim = 0; iDim < nDim; ++iDim)
      Mean_PrimVar[CReactiveNSVariable::VX_INDEX_PRIM + nDim] -= 0.5*(GridVel_i[iDim] + GridVel_j[iDim]);
  */

  /*--- Compute the vector from i to j ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    Edge_Vector[iDim] = Coord_j[iDim] - Coord_i[iDim];

  /*--- Mean gradient approximation ---*/
  if(!limiter) {
    for(iDim = 0; iDim < nDim; ++iDim) {
      /*--- Temperature ---*/
      Mean_GradPrimVar(T_INDEX_AVGGRAD,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::T_INDEX_GRAD][iDim] +
                                                    PrimVar_Grad_j[CReactiveNSVariable::T_INDEX_GRAD][iDim]);

      /*--- Velocities ---*/
      for(jDim = 0; jDim < nDim; ++jDim)
        Mean_GradPrimVar(VX_INDEX_AVGGRAD + jDim,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::VX_INDEX_GRAD + jDim][iDim] +
                                                              PrimVar_Grad_i[CReactiveNSVariable::VX_INDEX_GRAD + jDim][iDim]);

      /*--- Molar Fractions ---*/
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Mean_GradPrimVar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies][iDim] +
                                                                    PrimVar_Grad_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies][iDim]);
    }
  }
  else {
    for(iDim = 0; iDim < nDim; ++iDim) {
      /*--- Temperature ---*/
      Mean_GradPrimVar(T_INDEX_AVGGRAD,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::T_INDEX_GRAD][iDim]*
                                                    PrimVar_Lim_i[CReactiveNSVariable::T_INDEX_LIM] +
                                                    PrimVar_Grad_j[CReactiveNSVariable::T_INDEX_GRAD][iDim]*
                                                    PrimVar_Lim_j[CReactiveNSVariable::T_INDEX_LIM]);
      /*--- Velocities ---*/
      for(jDim = 0; jDim < nDim; ++jDim)
        Mean_GradPrimVar(VX_INDEX_AVGGRAD + jDim,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::VX_INDEX_GRAD + jDim][iDim]*
                                                              PrimVar_Lim_i[CReactiveNSVariable::VX_INDEX_LIM + jDim] +
                                                              PrimVar_Grad_i[CReactiveNSVariable::VX_INDEX_GRAD + jDim][iDim]*
                                                              PrimVar_Lim_j[CReactiveNSVariable::VX_INDEX_LIM + jDim]);
      /*--- Molar Fractions ---*/
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        Mean_GradPrimVar(RHOS_INDEX_AVGGRAD + iSpecies,iDim) = 0.5*(PrimVar_Grad_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies][iDim] +
                                                                    PrimVar_Grad_i[CReactiveNSVariable::RHOS_INDEX_GRAD + iSpecies][iDim]);
    }
  }

  Proj_Mean_GradPrimVar_Edge = Mean_GradPrimVar*Edge_Vector;
  su2double dist_ij_2 = std::inner_product(Edge_Vector.data(), Edge_Vector.data() + Edge_Vector.size(), Edge_Vector.data(), 0.0);
  if(dist_ij_2 > EPS) {
    Diff_PrimVar[T_INDEX_AVGGRAD] = V_j[CReactiveNSVariable::T_INDEX_PRIM] - V_i[CReactiveNSVariable::T_INDEX_PRIM];
    for(iDim = 0; iDim < nDim; ++iDim)
      Diff_PrimVar[VX_INDEX_AVGGRAD + iDim] = V_j[CReactiveNSVariable::VX_INDEX_PRIM + iDim] -
                                              V_i[CReactiveNSVariable::VX_INDEX_PRIM + iDim];
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Diff_PrimVar[RHOS_INDEX_AVGGRAD + iSpecies] = V_j[CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies] -
                                                    V_i[CReactiveNSVariable::RHOS_INDEX_PRIM + iSpecies];

    Mean_GradPrimVar -= (Proj_Mean_GradPrimVar_Edge - Diff_PrimVar)*Edge_Vector.transpose()/dist_ij_2;
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Common::wrap_in_unique(Normal), Mean_Laminar_Viscosity,
                     Mean_Thermal_Conductivity, Mean_Dij, config);

	/*--- Update viscous residual with the species projected flux ---*/
  std::copy(Proj_Flux_Tensor, Proj_Flux_Tensor + nVar, val_residual);

	/*--- Implicit part ---*/
	if(implicit)
    throw Common::NotImplemented("Implicit method still not implemented. Setting explicit");
}

//
//
/*!
 * \brief Constructor for the class CSourceChemistry
 */
//
//
CSourceReactive::CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
                 CNumerics(val_nDim,val_nVar,config),library(CReactiveEulerVariable::GetLibrary()),nSpecies(library->GetNSpecies()) {
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Ys.resize(nSpecies);
}

//
//
/*!
 * \brief Compute residual for chemistry source term
 */
//
//
void CSourceReactive::ComputeChemistry(su2double* val_residual, su2double** val_Jacobian_i, CConfig* config) {
  SU2_Assert(val_residual != NULL,"The array for residuals has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive varaibles has not been allocated");

  su2double rho, temp, dim_temp, dim_rho;

  /*--- Initialize residual array ---*/
  std::fill(val_residual, val_residual + nVar, 0.0);

  /*--- Nonequilibrium chemistry source term from library ---*/
  std::copy(V_i + CReactiveEulerVariable::RHOS_INDEX_PRIM, V_i + (CReactiveEulerVariable::RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  rho = V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];
  temp = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  dim_temp = temp*config->GetTemperature_Ref();
  dim_rho = rho*config->GetDensity_Ref();

  //omega = library->GetMassProductionTerm(dim_temp,dim_rho,Ys);

  /*--- Assign to the residual. NOTE: We need to invert the sign since it is a residual that will be ADDED to the total one ---*/
  std::copy(omega.cbegin(),omega.cend(),val_residual);
  for(unsigned short iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] *= -Volume/(config->GetDensity_Ref()*config->GetTime_Ref());

  /*--- Implicit computation ---*/
  if(implicit)
    throw Common::NotImplemented("Implicit method for source chemistry residual not yet implemented. Setting explicit");
}
