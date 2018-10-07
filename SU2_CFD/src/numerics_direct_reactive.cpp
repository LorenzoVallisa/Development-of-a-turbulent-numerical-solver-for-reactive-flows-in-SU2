#include "../include/numerics_reactive.hpp"

#include <algorithm>

namespace {

  /*!
   * \brief Compute unit normal for the current cell interface
   */
  void Compute_Outward_UnitNormal(const unsigned short nDim,su2double* Normal,su2double*& UnitNormal) {

    /*--- Face area (norm or the normal vector) ---*/
    su2double Area = 0.0;
    for (auto iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);

    /*-- Unit Normal ---*/
    for (auto iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;
  }

} /*-- End of unnamed namespace ---*/


/*!
 * \brief Constructor of the class CUpwReactiveAUSM
 */
CUpwReactiveAUSM::CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig>& config):
    CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
    nSpecies(library->GetNSpecies()) {

  nVar = val_nVar;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  ProjFlux_i.resize(nVar);
  ProjFlux_j.resize(nVar);

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

  unsigned short iDim, iVar, jVar, iSpecies; /*!< \brief Indexes for iterations. */
  su2double sq_vel_i,sq_vel_j,  /*!< \brief squared velocity. */
            Energy_i, Energy_j, /*!< \brief Energy at node i and at node j. */
            Temperature_i, Temperature_j, /*!< \brief Temperature at node i and at node j. */
	          ProjVelocity, ProjVelocity_i, ProjVelocity_j; /*!< \brief Projected velocities at node i and at node j. */
  RealVec   Velocity_i(nDim),Velocity_j(nDim); /*!< \brief Velocity at node i and at node j. */

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  /*--- Point i: compute energy,pressure,sound speed and enthalpy  ---*/
  sq_vel_i = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[CReactiveEulerVariable::VX_INDEX_PRIM+iDim];
		sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
  }

  Density_i = U_i[CReactiveEulerVariable::RHO_INDEX_SOL];
  Pressure_i = V_i[CReactiveEulerVariable::P_INDEX_PRIM];
  Energy_i = U_i[CReactiveEulerVariable::RHOE_INDEX_SOL]/Density_i;
  Enthalpy_i = Energy_i + Pressure_i/Density_i + 0.5*sq_vel_i;
  Temperature_i = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  library->Gamma_SoundSpeed(Temperature_i,Pressure_i,Density_i,Gamma,SoundSpeed_i);

  /*--- Point j: compute squared velocity,energy,pressure,sound speed and ethalpy  ---*/
  sq_vel_j = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = V_j[CReactiveEulerVariable::VX_INDEX_PRIM+nDim];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }

  Density_j = U_j[CReactiveEulerVariable::RHO_INDEX_SOL];
  Pressure_j = V_j[CReactiveEulerVariable::P_INDEX_PRIM];
  Energy_j = U_j[CReactiveEulerVariable::RHOE_INDEX_SOL]/Density_j;
  Enthalpy_j = Energy_j + Pressure_j/Density_j + 0.5*sq_vel_j;
  Temperature_i = V_i[CReactiveEulerVariable::T_INDEX_PRIM];
  library->Gamma_SoundSpeed(Temperature_j,Pressure_j,Density_j,Gamma,SoundSpeed_j);

  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0;
  ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  /*--- Calculate L/R Mach numbers ---*/
  su2double mL = std::sqrt(sq_vel_i)/SoundSpeed_i;
  su2double mR = std::sqrt(sq_vel_j)/SoundSpeed_j;
  /*--- Calculate Mean sound speed ---*/
  su2double MeanSoundSpeed = 0.5*(SoundSpeed_i+SoundSpeed_j);
  /*--- Calculate adjusted polynomial function AUSM +-Up ---*/
  su2double mRef,mF;
  if(mL>=1.0 or mR >= 1.0)
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
    pLP = 0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0);
  }
  else {
    mLP = 0.5*(mL+abs(mL));
    pLP = 0.5*(1.0+std::abs(mL)/mL);
  }

  if (std::abs(mR) < 1.0) {
    mRM = -0.25*(mR-1.0)*(mR-1.0)-0.125*(mR*mR-1.0)*(mR*mR-1.0);
    pRM = 0.25*(mR-1.0)*(mR-1.0)*(mR+2.0)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0);
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

  su2double pF = pLP + pRM;
  su2double pLF = pLP*Pressure_i + pRM*Pressure_j;
  mF = MeanSoundSpeed*(mLF*Density_i + mRF*Density_j);
  //auto Phi = abs(mF);

  RealVec Ys_i,Ys_j;
  Ys_i.reserve(nVar);
  Ys_j.reserve(nVar);
  Ys_i.push_back(Density_i);
  Ys_j.push_back(Density_j);
  for(iDim = 0; iDim<nDim; ++iDim) {
    Ys_i.push_back(Density_i*Velocity_i[iDim]);
    Ys_j.push_back(Density_j*Velocity_j[iDim]);
  }
  Ys_i.push_back(Density_i*Enthalpy_i);
  Ys_j.push_back(Density_j*Enthalpy_j);
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Ys_i.push_back(V_i[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]);
    Ys_j.push_back(V_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]);
  }
  /*--- Calculate the numerical flux ---*/
  for (iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] = 0.5*(mF*(Ys_i[iVar]+Ys_j[iVar]) - std::abs(mF)*(Ys_i[iVar]-Ys_j[iVar]))*Area;
  const su2double Ku = 1.0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    val_residual[CReactiveEulerVariable::VX_INDEX_PRIM+iDim] +=
    (pLF*UnitNormal[iDim] - Ku*pLP*pRM*(Density_i+Density_j)*mF*(ProjVelocity_i)*UnitNormal[iDim])*Area;
  }

  if(implicit) {

    /*--- Initialize the Jacobians ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        val_Jacobian_i[iVar][jVar] = 0.0;
        val_Jacobian_j[iVar][jVar] = 0.0;
      }
    }
    //GetInviscidProjFlux(&Density_i, Velocity_i.data(), &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i.data());
	  //GetInviscidProjFlux(&Density_j, Velocity_j.data(), &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j.data());

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

CAvgGradReactive_Flow::CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig>& config):
  CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  nVar = val_nVar;
  nPrimVarGrad = nDim + 3;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Mean_PrimVar.resize(nSpecies+nDim+2);
  PrimVar_i.resize(nSpecies+nDim+2);
  PrimVar_j.resize(nSpecies+nDim+2);
  Mean_GradPrimVar.resize(nSpecies+nDim+2,RealVec(nDim));
  Proj_Mean_GradPrimVar_Edge.resize(nSpecies+nDim+2);

}

//
//
/*!
 * \brief Compute projection of viscous fluxes
 */
//
//
void CAvgGradReactive_Flow::GetViscousProjFlux(RealVec& val_primvar, RealMatrix& val_gradprimvar,
                                               su2double* val_normal, RealVec& val_diffusioncoeff,
                                               su2double val_viscosity,su2double val_therm_conductivity,
                                               CConfig* config) {

  // Requires a non-standard primitive vector with molar masses instead of mass fractions
 	unsigned short iSpecies, iVar, iDim, jDim;
  RealVec Ds, V;
  RealMatrix GV;
  su2double mu, ktr, div_vel;
  su2double Ru;
  su2double rho, T, P;

  /*--- Initialize ---*/
  for (iVar = 0; iVar < nVar; ++iVar) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Flux_Tensor[iVar][iDim] = 0.0;
  }

  /*--- Rename for convenience ---*/
  Ds  = val_diffusioncoeff;
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  rho = val_primvar[CReactiveEulerVariable::RHO_INDEX_PRIM];
  T   = val_primvar[CReactiveEulerVariable::T_INDEX_PRIM];
  P   = val_primvar[CReactiveEulerVariable::P_INDEX_PRIM];
  Ru  = (*library).R_ungas;
  V   = val_primvar;
  GV  = val_gradprimvar;
  RealVec hs(nSpecies);
  library->GetSpeciesTotEnthalpies(T,P,hs);

  /*--- Calculate the velocity divergence ---*/
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; ++iDim)
    div_vel += GV[CReactiveEulerVariable::VX_INDEX_PRIM+iDim][iDim];

  /*--- Pre-compute mixture quantities ---*/
  RealVec Vec(nDim);
  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Vec[iDim] += rho*Ds[iSpecies]*GV[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies][iDim];
  }

  /*--- Compute the viscous stress tensor ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      tau[iDim][jDim] = 0.0;

  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++)
      tau[iDim][jDim] += mu * (GV[CReactiveEulerVariable::VX_INDEX_PRIM+jDim][iDim] + GV[CReactiveEulerVariable::VX_INDEX_PRIM+iDim][jDim]);
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux vector ---*/
 	for (iDim = 0; iDim < nDim; iDim++) {

    /*--- Species diffusion velocity ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Flux_Tensor[iSpecies][iDim] = rho*Ds[iSpecies]*GV[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies][iDim]
                                    - V[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]*Vec[iDim];
    }

    /*--- Shear stress related terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      Flux_Tensor[nSpecies+jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[nSpecies+nDim][iDim] += tau[iDim][jDim]*V[CReactiveEulerVariable::VX_INDEX_PRIM+jDim];
    }

    /*--- Diffusion terms ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Flux_Tensor[nSpecies+nDim][iDim]   += Flux_Tensor[iSpecies][iDim] * hs[iSpecies];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim]   += ktr*GV[CReactiveEulerVariable::T_INDEX_PRIM][iDim];
  }

  for (iVar = 0; iVar < nVar; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
    }
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
                                               su2double val_dist_ij, su2double* val_normal, su2double val_dS,
                                               su2double* val_Proj_Visc_Flux,su2double**& val_Proj_Jac_Tensor_i,su2double**& val_Proj_Jac_Tensor_j,
                                               CConfig* config) {

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

  unsigned short iDim,iVar,iSpecies,nVar_species; /*!< \brief Indexes for iterations. */
  RealVec Diffusion_Coefficient_Species_i(nSpecies), /*!< \brief Diffusion coefficients species at node i. */
          Diffusion_Coefficient_Species_j(nSpecies),  /*!< \brief Diffusion coefficients species at node j. */
          Mean_Diffusion_Coefficient(nSpecies); /*!< \brief Mean value of diffusion coefficient. */
  su2double Mean_Laminar_Viscosity; /*!< \brief Mean value of laminar viscosity. */
  su2double Mean_Thermal_Conductivity;  /*!< \brief Mean value of thermal conductivity. */

  library->GetRhoUdiff(Diffusion_Coefficient_Species_i);
  library->GetRhoUdiff(Diffusion_Coefficient_Species_j);

//  RealVec Edge_Vector(nDim);

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  /*--- Compute vector going from iPoint to jPoint ---*/
/*	dist_ij_2 = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
	}
*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coefficient[iSpecies] = 1.0/(1.0/Diffusion_Coefficient_Species_i[iSpecies] + 1.0/Diffusion_Coefficient_Species_j[iSpecies]);

  Mean_Laminar_Viscosity = 1.0/(1.0/Laminar_Viscosity_i + 1.0/Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 1.0/(1.0/Thermal_Conductivity_i + 1.0/Thermal_Conductivity_j);


   /*--- Set the primitive variables and compute the mean ---*/
	 for (iVar = 0; iVar < nVar; iVar++) {
     PrimVar_i[iVar] = V_i[iVar];
     PrimVar_j[iVar] = V_j[iVar];
     Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Mean gradient approximation ---*/
  // Mass fraction
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    PrimVar_i[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] = V_i[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]/V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];
    PrimVar_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] = V_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]/V_j[CReactiveEulerVariable::RHO_INDEX_PRIM];
    Mean_PrimVar[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] = 0.5*(PrimVar_i[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies] + PrimVar_j[CReactiveEulerVariable::RHOS_INDEX_PRIM+iSpecies]);
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iSpecies][iDim] = 0.5*(1.0/V_i[0] * (PrimVar_Grad_i[iSpecies][iDim] - PrimVar_i[iSpecies] * PrimVar_Grad_i[0][iDim]) +
                                              1.0/V_j[0] * (PrimVar_Grad_j[iSpecies][iDim] - PrimVar_j[iSpecies] * PrimVar_Grad_j[0][iDim]));
    }
  }
  /*for (iVar = nSpecies; iVar < PrimVar_i.size(); iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }*/
  for (iVar = nSpecies; iVar < nPrimVarGrad; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Diffusion_Coefficient, Mean_Laminar_Viscosity,Mean_Thermal_Conductivity, config);


		/*--- Update viscous residual with the species projected flux ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			val_residual[iVar] = Proj_Flux_Tensor[iVar];

		/*--- Implicit part ---*/
		if (implicit) {

        RealVec Edge_Vector(nDim);

        /*--- Compute vector going from iPoint to jPoint ---*/
        su2double	dist_ij_2 = 0.0;
      	for (iDim = 0; iDim < nDim; iDim++) {
      		Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
      		dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
      	}

			  GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar,Mean_Diffusion_Coefficient, Mean_Laminar_Viscosity,Mean_Thermal_Conductivity,
                           std::sqrt(dist_ij_2), UnitNormal, Area, Proj_Flux_Tensor, val_Jacobian_i, val_Jacobian_j,config);
    }

}

//
//
/*!
 * \brief Constructor for the class CSourceChemistry
 */
//
//

CSourceReactive::CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig>& config):
  CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
  nSpecies(library->GetNSpecies()) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

}
