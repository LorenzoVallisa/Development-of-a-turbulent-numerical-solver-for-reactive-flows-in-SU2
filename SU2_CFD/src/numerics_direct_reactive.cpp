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
CUpwReactiveAUSM::CUpwReactiveAUSM(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig> config):
    CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Velocity_i.resize(nDim);
  Velocity_j.resize(nDim);
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


  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+2); AD::SetPreaccIn(V_j, nDim+2);

  unsigned short iDim, iVar, jVar, iSpecies; /*!< \brief Indexes for iterations. */
  su2double sq_vel_i,sq_vel_j, /*!< \brief squared velocity. */
            Energy_i, Energy_j,/*!< \brief Energy at node i and at node j. */
	          ProjVelocity, ProjVelocity_i, ProjVelocity_j; /*!< \brief Projected velocities at node i and at node j. */

  ::Compute_Outward_UnitNormal(nDim,Normal,UnitNormal);

  /*--- Point i: compute energy,pressure,sound speed and enthalpy  ---*/
  sq_vel_i = 0.0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = V_i[VX_INDEX_PRIM+iDim];
		sq_vel_i += Velocity_i[iDim]*Velocity_i[iDim];
  }

  using CReactiveEulerVariable::RHO_INDEX_SOL = RHO_INDEX_SOL;
  using CReactiveEulerVariable::VX_INDEX_PRIM = VX_INDEX_PRIM;
  using CReactiveEulerVariable::RHOS_INDEX_PRIM = RHOS_INDEX_PRIM;
  using CReactiveEulerVariable::P_INDEX_PRIM = P_INDEX_PRIM;
  using CReactiveEulerVariable::RHOE_INDEX_SOL = RHOE_INDEX_SOL;

  Density_i = U_i[RHO_INDEX_SOL];
  Pressure_i = V_i[P_INDEX_PRIM];
  Energy_i = U_i[RHOE_INDEX_SOL]/Density_i;
  Enthalpy_i = Energy_i + Pressure_i/Density_i + 0.5*sq_vel_i;
  SoundSpeed_i = std::sqrt(Pressure_i*Gamma/Density_i);

  /*--- Point j: compute squared velocity,energy,pressure,sound speed and ethalpy  ---*/
  sq_vel_j = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = V_j[VX_INDEX_PRIM+nDim];
    sq_vel_j += Velocity_j[iDim]*Velocity_j[iDim];
  }

  Density_j = U_j[RHO_INDEX_SOL];
  Pressure_j = V_j[P_INDEX_PRIM];
  Energy_j = U_j[RHOE_INDEX_SOL]/Density_j;
  Enthalpy_j = Energy_j + Pressure_j/Density_j + 0.5*sq_vel_j;
  SoundSpeed_j = std::sqrt(Pressure_j*Gamma/Density_j);

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
    Ys_i.push_back(V_i[RHOS_INDEX_PRIM+iSpecies]);
    Ys_j.push_back(V_j[RHOS_INDEX_PRIM+iSpecies];
  }
  /*--- Calculate the numerical flux ---*/
  for (iVar = 0; iVar < nVar; ++iVar)
    val_residual[iVar] = 0.5*(mF*(Ys_i[iVar]+Ys_j[iVar]) - std::abs(mF)*(Ys_i[iVar]-Ys_j[iVar]))
  const su2double Ku = 1.0;
  for (iDim = 0; iDim < nDim; ++iDim)
    val_residual[VX_INDEX_PRIM+iDim] +=
    pLF*UnitNormal[iDim] - Ku*pLP*pRM*(Density_i+Density_j)*mF*(ProjVelocity_i)*UnitNormal[iDim]

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

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

//
//
/*!
 * \brief Constructor of the class CAvgGradReactive_Flow
 */
//
//

CAvgGradReactive_Flow::CAvgGradReactive_Flow(unsigned short val_nDim, unsigned short val_nVar, std::unique_ptr<CConfig> config):
  CNumerics(val_nDim,val_nVar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

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
 * \brief Compute residual for viscous term
 */
//
//

void CAvgGradReactive_Flow::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                                            su2double** val_Jacobian_j, CConfig* config) {

  unsigned short iDim,iVar,iSpecies,nVar_species; /*!< \brief Indexes for iterations. */
  RealVec Diffusion_Coefficient_Species_i(nSpecies), /*!< \brief Diffusion coefficients species at node i. */
          Diffusion_Coefficient_Species_j(nSpecies),  /*!< \brief Diffusion coefficients species at node j. */
          Mean_Diffusion_Coefficient(nSpecies); /*!< \brief Mean value of diffusion coefficient. */

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
  using CReactiveEulerVariable::RHO_INDEX_PRIM = RHO_INDEX_PRIM;
  using CReactiveEulerVariable::RHOS_INDEX_PRIM = RHOS_INDEX_PRIM;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coefficient_Viscosity[iSpecies] = 0.5*(Diffusion_Coefficient_Species_i[iSpecies] + Diffusion_Coefficient_Species_j[iSpecies]);

  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i + Thermal_Conductivity_j);


   /*--- Set the primitive variables and compute the mean ---*/
	 for (iVar = 0; iVar < nVar; iVar++) {
     PrimVar_i[iVar] = V_i[iVar];
     PrimVar_j[iVar] = V_j[iVar];
     Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }

  /*--- Mean gradient approximation ---*/
  // Mass fraction
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    PrimVar_i[RHOS_INDEX_PRIM+iSpecies] = V_i[RHOS_INDEX_PRIM+iSpecies]/V_i[RHO_INDEX_PRIM];
    PrimVar_j[RHOS_INDEX_PRIM+iSpecies] = V_j[RHOS_INDEX_PRIM+iSpecies]/V_j[RHO_INDEX_PRIM];
    Mean_PrimVar[RHOS_INDEX_PRIM+iSpecies] = 0.5*(PrimVar_i[RHOS_INDEX_PRIM+iSpecies] + PrimVar_j[RHOS_INDEX_PRIM+iSpecies]);
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
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Normal, Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,Mean_Thermal_Conductivity, config);


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

			  GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar,Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,Mean_Thermal_Conductivity,
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

CSourceReactive::CSourceReactive(unsigned short val_nDim, unsigned short val_nVar, CConfig std::unique_ptr<CConfig>config):
  CNumerics(val_nDim,val_nVar,config),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),nSpecies(library->GetNSpecies()) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

}
