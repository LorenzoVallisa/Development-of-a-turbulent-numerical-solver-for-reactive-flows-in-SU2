#include "../include/variable_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"

//
//
/*--- Class default constructor ---*/
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
/*--- Class constructor ---*/
//
//

CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim,unsigned short val_nvar, std::uniuqe_ptr<CConfig> config):
                        CVariable(val_nDim,val_nvar,config.get()),library(new Framework::ReactingModelLibrary(config->GetLibraryName())),
                        nSpecies(library->GetNSpecies()) {

  nPrimVar = nSpecies + nDim + 5;
  nPrimVarGrad = nSpecies + nDim + 1;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;

  Primitive.resize(nPrimVar);
  Gradient_Primitive.resize(nPrimVarGrad,RealVec(nDim));
  Limiter_Primitive.resize(nPrimVarGrad);

}

//
//
/*--- Get gradient primitive variables ---*/
//
//

su2double** CReactiveEulerVariable::GetGradient_Primitive(void) {
  std::vector<su2double*> tmp;
  for(auto& elem:Gradient_Primitive)
    tmp.push_back(elem.data());
  return tmp.data();
}

//
//
/*--- Set density ---*/
//
//

bool CReactiveEulerVariable::SetDensity(void) {

  unsigned short iSpecies;
  su2double Density = 0.0;

  for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    Primitive.at(RHOS_INDEX_PRIM+iSpecies) = Solution[RHOS_INDEX_SOL+iSpecies];
    Density += Solution[RHOS_INDEX_SOL+iSpecies];
  }
  Primitive.at(RHO_INDEX) = Density;

  return false;
}

//
//
/*--- Set Pressure ---*/
//
//

void CReactiveEulerVariable::SetPressure(void) {

  unsigned short iSpecies;
  su2double Pressure = 0.0;
  RealVec Ms(nSpecies);

  /*--- Read gas mixture properties from library ---*/
  library->GetMolarMasses(Ms);

  su2double Ru = (*library).R_ungas;

  /*--- Solve for mixture pressure using ideal gas law & Dalton's law ---*/
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Pressure += Solution[RHOS_INDEX_SOL+iSpecies] * Ru/Ms[iSpecies] * Primitive[T_INDEX_PRIM];

  /*--- Store computed values and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if (Primitive[P] < 0.0)
    throw std::runtime_error("Non Physical solution");
}

//
//
/*--- Set sound speed ---*/
//
//

bool CReactiveEulerVariable::SetSoundSpeed(void) {

  su2double Gamma = 0.0;
  su2double Sound_Speed = 0.0;
  library->Gamma_SoundSpeed(Primitive.at(T_INDEX_PRIM),Primitive.at(P_INDEX_PRIM),Primitive.at(RHO_INDEX_PRIM),Gamma,Sound_Speed);
  Primitive.at(A_INDEX_PRIM) = Sound_Speed;
  if(Sound_Speed < 0.0)
    return false;
  return true;

}

//
//
/*--- Compute enthalpy for desired species ---*/
//
//

/*su2double CReactiveEulerVariable::CalcHs(su2double val_T, unsigned short val_Species) {

  assert(val_Species < nSpecies);

  RealVec Ms(nSpecies), hf(nSpecies);
  su2double hs = 0.0;*/

  /*--- Read from library ---*/
  /*library->GetMolarMasses(Ms);
  library->GetFormationEnthalpies(hf);

  su2double Ru = (*library).R_ungas;

  hs = 5.0/2.0*Ru/Ms[val_Species]*val_T + hf[val_Species];

  return hs;
}*/

//
//
/*--- Compute norm projected velocity ---*/
//
//

su2double CReactiveEulerVariable::GetProjVel(su2double* val_vector) {

  su2double ProjVel = 0.0;
  //su2double density = 0.0;
	unsigned short iDim, iSpecies;

  //for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //  density += Solution[RHOS_INDEX_SOL+iSpecies];
	for (iDim = 0; iDim < nDim; iDim++)
		ProjVel += Solution[RHOVX_INDEX_SOL+iDim]*val_vector[iDim]/Solution[RHO_INDEX];

	return ProjVel;
}

//
//
/*--- Class constructor ---*/
//
//

CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim,unsigned short val_nvar, std::unique_ptr<CConfig> config):
                     CReactiveEulerVariable(val_nDim,val_nvar,config.get()) {

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
