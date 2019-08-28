#include "../include/Framework/reacting_model_library.hpp"
#include "../include/Framework/not_setup_exception.hpp"
#include "../include/Tools/spline.hpp"
#include "../include/Tools/utility.hpp"

#include <experimental/filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

namespace Framework {
  //
  //
  /*--- Setting gas constant for each species ---*/
  void ReactingModelLibrary::SetRiGas(void) {
    Ri.clear();
    Ri.reserve(nSpecies);
    for(const auto& mass: mMasses)
      Ri.emplace_back(R_ungas/mass);
  }

  //
  //
  /*--- Setting gas constant for mixture ---*/
  inline void ReactingModelLibrary::SetRgas(const RealVec& ys) {
    /*--- Check the correct size of vectors ---*/
    SU2_Assert(ys.size() == nSpecies, "The dimension of vector ys doesn't match nSpecies");

    SetMassFractions(ys);
    Rgas = std::inner_product(Ys.cbegin(), Ys.cend(), Ri.cbegin(), 0.0);
  }

  //
  //
  /*--- Computing gas constant for mixture ---*/
  inline double ReactingModelLibrary::ComputeRgas(const RealVec& ys) {
    SetRgas(ys);
    return Rgas;
  }

  //
  //
  /*--- Setting molar fraction for each species ---*/
  void ReactingModelLibrary::SetMolarFractions(const RealVec& xs) {
    /*--- Check the correct size of vectors ---*/
    SU2_Assert(xs.size() == nSpecies, "The dimension of vector xs doesn't match nSpecies");
    Xs = xs;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      /*--- Adjust vanishing species mole fractions ---*/
      if(Xs[iSpecies] < 0.0)
        Xs[iSpecies] = 1.0e-30;

      /*--- Check physical value for mole fractions. NOTE: We add an extra tolerance in case of numerical perturbation for Jacobian
            or take into account some small instabilities ---*/
      SU2_Assert(Xs[iSpecies] <= (1.0 + 1.0e-1), std::string("The molar fraction of species number " +
                                                              std::to_string(iSpecies) + "is too greater than 1"));
    }
  }

  //
  //
  /*--- Setting mass fraction for each species ---*/
  void ReactingModelLibrary::SetMassFractions(const RealVec& ys) {
    /*--- Check the correct size of vectors ---*/
    SU2_Assert(ys.size() == nSpecies, "The dimension of vector ys doesn't match nSpecies");
    Ys = ys;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      /*--- Adjust vanishing species mass fractions ---*/
      if(Ys[iSpecies] < 0.0)
        Ys[iSpecies] = 1.0e-30;

      /*--- Check physical value for mass fractions ---*/
      SU2_Assert(Ys[iSpecies] <= (1.0 + 1.0e-1), std::string("The mass fraction of species number " +
                                                              std::to_string(iSpecies) + " is too greater than 1"));
    }
  }

  //
  //
  /*--- This function sets the molar fractions from mass fractions. ---*/
  void ReactingModelLibrary::SetMolarFromMass(const RealVec& ys) {
    SetMassFractions(ys);

    unsigned short iSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Xs[iSpecies] = Ys[iSpecies]/mMasses[iSpecies];

    double massTot = std::accumulate(Ys.cbegin(), Ys.cend(), 0.0)/std::accumulate(Xs.cbegin(), Xs.cend(), 0.0);
    std::transform(Xs.begin(), Xs.end(), Xs.begin(), std::bind1st(std::multiplies<double>(),massTot));
  }

  //
  //
  /*--- This function returns the molar fractions from mass fractions. ---*/
  inline RealVec ReactingModelLibrary::GetMolarFromMass(const RealVec& ys) {
    SetMolarFromMass(ys);
    return Xs;
  }

  //
  //
  /*--- This function sets the mass fractions from molar fractions. ---*/
  void ReactingModelLibrary::SetMassFromMolar(const RealVec& xs) {
    SetMolarFractions(xs);

    unsigned short iSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Ys[iSpecies] = Xs[iSpecies]*mMasses[iSpecies];

    double massTot = std::accumulate(Xs.cbegin(), Xs.cend(), 0.0)/std::accumulate(Ys.cbegin(), Ys.cend(), 0.0);
    std::transform(Ys.begin(), Ys.end(), Ys.begin(), std::bind1st(std::multiplies<double>(),massTot));
  }

  //
  //
  /*--- This function returns the mass fractions from molar fractions. ---*/
  inline RealVec ReactingModelLibrary::GetMassFromMolar(const RealVec& xs) {
    SetMassFromMolar(xs);
    return Ys;
  }

  //
  //
  /*--- This function computes gamma and the frozen sound speed. ---*/
  inline void ReactingModelLibrary::Gamma_FrozenSoundSpeed(const double temp, const RealVec& ys, double& gamma, double& sound_speed) {
    double Cp = ComputeCP(temp,ys);
    double Rgas = std::inner_product(Ys.cbegin(), Ys.cend(), Ri.cbegin(), 0.0);
    double Cv = Cp - Rgas;
    gamma = Cp/Cv;
    sound_speed = std::sqrt(gamma*Rgas*temp);
  }

  //
  //
  /*--- This function computes the frozen gamma from temperature. ---*/
  inline double ReactingModelLibrary::ComputeFrozenGamma(const double temp, const RealVec& ys) {
    double Cp = ComputeCP(temp,ys);
    double Rgas = std::inner_product(Ys.cbegin(), Ys.cend(), Ri.cbegin(), 0.0);
    double Cv = Cp - Rgas;
    return Cp/Cv;
  }

  //
  //
  /*--- This function computes the frozen sound speed from temeprature. ---*/
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed(const double temp, const RealVec& ys) {
    double gamma = ComputeFrozenGamma(temp,ys);
    return std::sqrt(gamma*Rgas*temp);
  }

  //
  //
  /*--- This function computes the frozen sound speed when we already computed gamma. ---*/
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed_FromGamma(const double temp, const double gamma, const RealVec& ys) {
    SetRgas(ys);
    return std::sqrt(gamma*Rgas*temp);
  }

  //
  //
  /*--- This function computes the frozen sound speed when we already computed gamma. ---*/
  inline double ReactingModelLibrary::ComputeFrozenGamma_FromSoundSpeed(const double temp, const double sound_speed, const RealVec& ys) {
    SetRgas(ys);
    return sound_speed*sound_speed/(Rgas*temp);
  }

  //
  //
  /*--- This function computes the frozen sound speed with pressure and density. ---*/
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed(const double temp,  const RealVec& ys,
                                                              const double press, const double rho) {
    double gamma = ComputeFrozenGamma(temp,ys);
    return std::sqrt(gamma*press/rho);
  }

  //
  //
  /*--- This function computes the frozen sound speed when we already computed gamma. ---*/
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed_FromGamma(const double gamma, const RealVec& ys,
                                                                        const double press, const double rho) {
    return std::sqrt(gamma*press/rho);
  }

  //
  //
  /*--- This function computes pressure at given temperature and density. ---*/
  inline double ReactingModelLibrary::ComputePressure(const double temp, const double rho, const RealVec& ys) {
    SetRgas(ys);
    return rho*temp*Rgas;
  }

  //
  //
  /*--- This function computes density at given temperature and pressure. ---*/
  inline double ReactingModelLibrary::ComputeDensity(const double temp, const double pressure, const RealVec& ys) {
    SetRgas(ys);
    return pressure/(temp*Rgas);
  }

  //
  //
  /*--- This function computes temperature at given density and pressure. ---*/
  inline double ReactingModelLibrary::ComputeTemperature(const double pressure, const double rho, const RealVec& ys) {
    SetRgas(ys);
    return pressure/(rho*Rgas);
  }

  //
  //
  /*--- This function computes density at given temperature and pressure. ---*/
  inline double ReactingModelLibrary::ComputeTemperature_FromGamma(const double sound_speed2, const double gamma, const RealVec& ys) {
    SetRgas(ys);
    return sound_speed2/(gamma*Rgas);
  }

  //
  //
  /*--- This function computes the mixture internal energy at given temperature and pressure. ---*/
  inline double ReactingModelLibrary::ComputeEnergy(const double temp, const RealVec& ys) {
    double enthalpy = ComputeEnthalpy(temp, ys);
    SetRgas(ys);
    return enthalpy - temp*Rgas;
  }

  //
  //
  /*--- This function computes density,internal energy and static enthalpy at given temperature and pressure. ---*/
  inline void ReactingModelLibrary::Density_Enthalpy_Energy(const double temp, const double pressure, const RealVec& ys, RealVec& dhe) {
    dhe.resize(3);
    dhe[0] = ComputeDensity(temp, pressure, ys);
    dhe[1] = ComputeEnthalpy(temp, ys);
    dhe[2] = dhe[1] - pressure/dhe[0];
  }

  //
  //
  /*--- This function sets the static enthalpy for each species. ---*/
  void ReactingModelLibrary::SetPartialEnthalpy(const double temp) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Enthalpies[iSpecies] = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                                  std::get<X_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                                  std::get<Y_DATA_SPLINE>(Enth_Spline[iSpecies]),temp)/mMasses[iSpecies];
  }


  //
  //
  /*--- This function returns the computed static enthalpy for each species. ---*/
  RealVec ReactingModelLibrary::ComputePartialEnthalpy(const double temp) {
    SetPartialEnthalpy(temp);
    return Enthalpies;
  }

  //
  //
  /*--- This function computes the static enthalpy of the mixture. ---*/
  double ReactingModelLibrary::ComputeEnthalpy(const double temp, const RealVec& ys) {
    SetPartialEnthalpy(temp);
    SetMassFractions(ys);
    return std::inner_product(Ys.cbegin(), Ys.cend(), Enthalpies.cbegin(), 0.0);
  }

  //
  //
  /*--- This function sets the internal energy for each species. ---*/
  void ReactingModelLibrary::SetPartialEnergy(const double temp) {
    SetPartialEnthalpy(temp);
    /*--- Set internal energies ---*/
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Internal_Energies[iSpecies] = Enthalpies[iSpecies] - Ri[iSpecies]*temp;
  }

  //
  //
  /*--- This function computes the internal energy for each species. ---*/
  RealVec ReactingModelLibrary::ComputePartialEnergy(const double temp) {
    SetPartialEnergy(temp);
    return Internal_Energies;
  }

  //
  //
  /* This function computes the pressure derivative w.r.t. partial densities. ---*/
  RealVec ReactingModelLibrary::ComputedP_dYs(const double temp, const double gamma) {
    SetPartialEnergy(temp);
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      dPdYs[iSpecies] = Ri[iSpecies]*temp - (gamma - 1.0)*Internal_Energies[iSpecies];

    return dPdYs;
  }

  //
  //
  /*--- This function sets the static enthalpy for each species. ---*/
  void ReactingModelLibrary::SetPartialEnthalpy(const double temp, unsigned short iSpecies) {
    Enthalpies.at(iSpecies) = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                                   std::get<X_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                                   std::get<Y_DATA_SPLINE>(Enth_Spline[iSpecies]),temp)/mMasses[iSpecies];
  }

  //
  //
  /*--- This function returns the computed static enthalpy for each species. ---*/
  double ReactingModelLibrary::ComputePartialEnthalpy(const double temp, unsigned short iSpecies) {
    SetPartialEnthalpy(temp, iSpecies);
    return Enthalpies[iSpecies];
  }

  //
  //
  /*--- This function sets the internal energy for each species. ---*/
  void ReactingModelLibrary::SetPartialEnergy(const double temp, unsigned short iSpecies) {
    SetPartialEnthalpy(temp, iSpecies);
    /*--- Set internal energies ---*/
    Internal_Energies.at(iSpecies) = Enthalpies.at(iSpecies) - Ri.at(iSpecies)*temp;
  }

  //
  //
  /*--- This function computes the internal energy for each species. ---*/
  double ReactingModelLibrary::ComputePartialEnergy(const double temp, unsigned short iSpecies) {
    SetPartialEnergy(temp, iSpecies);
    return Internal_Energies[iSpecies];
  }

  //
  //
  /* This function computes the pressure derivative w.r.t. partial densities. ---*/
  double ReactingModelLibrary::ComputedP_dYs(const double temp, const double gamma, unsigned short iSpecies) {
    SetPartialEnergy(temp, iSpecies);
    dPdYs.at(iSpecies) = Ri[iSpecies]*temp - (gamma - 1.0)*Internal_Energies[iSpecies];

    return dPdYs[iSpecies];
  }

  //
  //
  /*--- This function computes the specific heat at constant pressure for each species. ---*/
  RealVec ReactingModelLibrary::ComputeCps(const double temp) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      CPs[iSpecies] = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Cp_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Cp_Spline[iSpecies]),
                                           std::get<Y_DATA_SPLINE>(Cp_Spline[iSpecies]),temp)/mMasses[iSpecies];

    return CPs;
  }

  //
  //
  /*--- This function computes the specific heat at constant pressure. ---*/
  double ReactingModelLibrary::ComputeCP(const double temp, const RealVec& ys) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      CPs[iSpecies] = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Cp_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Cp_Spline[iSpecies]),
                                           std::get<Y_DATA_SPLINE>(Cp_Spline[iSpecies]),temp)/mMasses[iSpecies];

    SetMassFractions(ys);
    return std::inner_product(Ys.cbegin(), Ys.cend(), CPs.cbegin(), 0.0);
  }

  //
  //
  /*--- Computing molecular viscosity of each species. ---*/
  void ReactingModelLibrary::ComputeViscosities(const double temp) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Viscosities[iSpecies] = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Mu_Spline[iSpecies]),
                                                   std::get<X_DATA_SPLINE>(Mu_Spline[iSpecies]),
                                                   std::get<Y_DATA_SPLINE>(Mu_Spline[iSpecies]),temp);
  }

  //
  //
  /*--- Computing viscosity of the mixture. ---*/
  double ReactingModelLibrary::ComputeEta(const double temp, const RealVec& ys) {
    ComputeViscosities(temp);

    /*--- Setting mass fractions and computing Ys/mMasses ---*/
    SetMassFractions(ys);
    unsigned short iSpecies, jSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      ys_over_mm[iSpecies] = Ys[iSpecies]/mMasses[iSpecies];

    /*--- Mixture viscosity calculation as sum weighted over PHI ---*/
    double phi;
    double Viscosity_Mixture = 0.0;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)  {
      phi = 0.0;
      for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
        phi += ys_over_mm[jSpecies]/std::sqrt(8.0*(1.0 + mMasses[iSpecies]/mMasses[jSpecies]))*
               (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25))*
               (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25));

      Viscosity_Mixture += Viscosities[iSpecies]*ys_over_mm[iSpecies]/phi;
    }
    return Viscosity_Mixture;
  }

  //
  //
  /*--- Computing thermal conductivity of each species ---*/
  void ReactingModelLibrary::ComputeConductivities(const double temp) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Thermal_Conductivities[iSpecies] = MathTools::GetSpline(std::get<T_DATA_SPLINE>(Kappa_Spline[iSpecies]),
                                                              std::get<X_DATA_SPLINE>(Kappa_Spline[iSpecies]),
                                                              std::get<Y_DATA_SPLINE>(Kappa_Spline[iSpecies]),temp);
  }

  //
  //
  /*--- Computing thermal conductivity of the mixtures---*/
  double ReactingModelLibrary::ComputeLambda(const double temp, const RealVec& ys) {
    ComputeConductivities(temp);
    ComputeViscosities(temp);

    /*--- Setting mass fractions and computing Ys/mMasses ---*/
    SetMassFractions(ys);
    unsigned short iSpecies,jSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      ys_over_mm[iSpecies] = ys[iSpecies]/mMasses[iSpecies];

    /*--- Mixture thermal condictivity calculation as sum weighted over PHI ---*/
    double phi;
    double Thermal_Conductivity_Mixture = 0.0;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      phi = 0.0;
	    for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
        if(jSpecies != iSpecies)
          phi += 1.065*ys_over_mm[jSpecies]/std::sqrt(8.0*(1.0 + mMasses[iSpecies]/mMasses[jSpecies]))*
                      (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25))*
                      (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25));

      phi += ys_over_mm[iSpecies];
      Thermal_Conductivity_Mixture += Thermal_Conductivities[iSpecies]*ys_over_mm[iSpecies]/phi;
    }
	  return Thermal_Conductivity_Mixture;
  }

  //
  //
  /* This function computes the species concetration. */
  inline void ReactingModelLibrary::SetConcentration(const double rho, const RealVec& ys) {
    SetMassFractions(ys);
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Cs[iSpecies] = 1.0e3*rho*Ys[iSpecies]/mMasses[iSpecies];
  }

  //
  //
  /* This function computes the mixture concetration. */
  inline double ReactingModelLibrary::ComputeConcentration(const double rho, const RealVec& ys) {
    SetConcentration(rho,ys);
    return std::accumulate(Cs.cbegin(), Cs.cend(), 0.0);
  }

  //
  //
  /* This function computes the species diffusion in case of a constant Lewis number */
  RealVec ReactingModelLibrary::GetRhoUdiff(const double temp, const double rho, const RealVec& ys) {
    double kappa = ComputeLambda(temp,ys);
    double Cp = ComputeCP(temp,ys);
    std::fill(rhoUdiff.begin(), rhoUdiff.end(), kappa/(rho*Cp*Le));
    return rhoUdiff;
  }

  //
  //
  /* This function computes the effective diffusion coefficients for each species. */
  RealVec ReactingModelLibrary::GetDiffCoeffs(const double temp, const double pressure, const RealVec& ys) {
    /*--- Compute binary diffusion coefficients and mole fractions ---*/
    Dij = GetDij_SM(temp,pressure);
    SetMolarFromMass(ys);

    /*--- Compute mean effective diffusion coefficients ---*/
    double tmp;
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      tmp = 0.0;
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        if(iSpecies != jSpecies)
          tmp += Xs[jSpecies]/Dij(iSpecies,jSpecies);
      }
      Dm_coeffs[iSpecies] = ((1.0 - Xs[iSpecies])/tmp);
    }
    return Dm_coeffs;
  }

  //
  //
  /* This function computes the binary diffusion coefficients for Stefan-Maxwell diffusion with an empirical formula */
  RealMatrix ReactingModelLibrary::GetDij_SM(const double temp, const double pressure) {
    /*--- Local variables ---*/
    double Mij, diff_vol_i, molar_mass_i, diff_vol_j;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      molar_mass_i = mMasses[iSpecies];
      diff_vol_i = std::cbrt(Diff_Volumes[iSpecies]);
      for(unsigned short jSpecies = iSpecies; jSpecies < nSpecies; ++jSpecies) {
        Mij = std::sqrt((molar_mass_i*mMasses[jSpecies])/(molar_mass_i + mMasses[jSpecies]));
        diff_vol_j = std::cbrt(Diff_Volumes[jSpecies]);
        Dij(iSpecies,jSpecies) = 1.0e-3*std::pow(temp,1.75)/(pressure*Mij*(diff_vol_i + diff_vol_j)*(diff_vol_i + diff_vol_j));
        Dij(jSpecies,iSpecies) = Dij(iSpecies,jSpecies);
      }
    }
    return Dij;
  }

  //
  //
  /* This function computes the original matrix for Stefan-Maxwell equations */
  RealMatrix ReactingModelLibrary::GetGamma(const double rho, const RealVec& val_xs, const RealVec& val_ys,
                                                                  const RealMatrix& val_Dij) {
    /*--- Check the correct size of vectors ---*/
    SU2_Assert(val_ys.size() == nSpecies, "The dimension of vector val_ys doesn't match nSpecies");
    SU2_Assert(val_xs.size() == nSpecies, "The dimension of vector val_xs doesn't match nSpecies");

    const double sigma = std::accumulate(val_ys.cbegin(), val_ys.cend(), 0.0);
    double massTot = 0.0;
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      massTot += val_ys[iSpecies]/mMasses[iSpecies];
    massTot = 1.0/massTot;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        if(iSpecies != jSpecies)
          Gamma(iSpecies,jSpecies) = -sigma*massTot*val_xs[iSpecies]/(rho*mMasses[jSpecies]*val_Dij(iSpecies,jSpecies));
        else {
          double tmp = 0.0;
          for(unsigned short kSpecies = 0; kSpecies < nSpecies; ++kSpecies) {
            if(kSpecies != iSpecies)
              tmp += val_xs[kSpecies]/val_Dij(iSpecies,kSpecies);
          }
          Gamma(iSpecies,iSpecies) = sigma*massTot*tmp/(rho*mMasses[iSpecies]);
        }
      }
    }
    return Gamma;
  }

  //
  //
  /* This function computes the reaction rates constants for a specific reaction. */
  std::pair<double,double> ReactingModelLibrary::ComputeKeq(const double temp, unsigned short iReac) {
    /*--- Check correct index passing ---*/
    SU2_Assert(iReac < nReactions, "The index of reaction exceeds the number of reactions detected");

    /*--- Computing reaction equilibrium constant ---*/
    double dG = 0.0;
    double dnu = 0.0;
    unsigned short iSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      double dcoeff = Stoich_Coeffs_Products(iSpecies,iReac) - Stoich_Coeffs_Reactants(iSpecies,iReac);
      if(dcoeff != 0.0) {
        dG += dcoeff*
              (MathTools::GetSpline(std::get<T_DATA_SPLINE>(Enth_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                    std::get<Y_DATA_SPLINE>(Enth_Spline[iSpecies]),temp) - temp*
               MathTools::GetSpline(std::get<T_DATA_SPLINE>(Entr_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Entr_Spline[iSpecies]),
                                    std::get<Y_DATA_SPLINE>(Entr_Spline[iSpecies]),temp));
        dnu += dcoeff;
      }
    }
    const double RT = R_ungas*temp;
    const double lnKp = -dG/RT;
    const double lnKc = lnKp - dnu*std::log(R_ungas_atm*temp);
    Kc[iReac] = std::exp(lnKc);

    /*--- Return the equilibrium constants ---*/
    return std::make_pair(Kc[iReac], std::exp(lnKp));
  }


  //
  //
  /* This function computes the reaction rates constants for a specific reaction. */
  std::pair<double,double> ReactingModelLibrary::ComputeRateConstants(const double temp, unsigned short iReac) {
    /*--- Local variables ---*/
    double kf,kb;

    /*--- Forward rate constant ---*/
    kf = As[iReac]*std::pow(temp,Betas[iReac])*std::exp(-Temps_Activation[iReac]/temp);

    /*--- If the reaction is irreversible or backward data were not available use Gibbs to compute equilibrium constant ---*/
    if(Available_Backward_Rate.count(iReac) == 0) {
      /*--- Compute equilibrium constants (concentration and pressure) ---*/
      auto Keqs = ComputeKeq(temp, iReac);

      /*--- Check whether the reaction is complete or not ---*/
      bool is_complete = (Keqs.second > 1.0e10);
      if(!Reversible_Reactions[iReac]) {
        if(!is_complete)
          std::cerr<<"The equilibrium constant of reaction " + std::to_string(iReac) +
                     " is not so big even if the reaction is considered 'irreversible'"<<std::endl;
        kb = 0.0;
      }
      else if(is_complete)
        kb = 0.0;
      else
        kb = kf/Keqs.first;
    }
    /*--- If available use backward data previously read ---*/
    else {
      kb = As_back[iReac]*std::pow(temp,Betas_back[iReac])*std::exp(-Temps_Activation_back[iReac]/temp);
      Kc[iReac] = kf/kb;
    }

    return std::make_pair(kf,kb);
  }

  //
  //
  /* This function sets the forward and backward reaction rates. */
  void ReactingModelLibrary::SetReactionRates(const double temp, const double rho, const RealVec& ys) {
    SetConcentration(rho, ys);

    double for_rate, back_rate;
    /*--- Map concetnration to EIGEN for useful next computations ---*/
    auto cs = Eigen::Map<Eigen::ArrayXd>(Cs.data(), Cs.size());

    for(unsigned short iReac = 0; iReac < nReactions; ++iReac) {
      auto Rate_Const = ComputeRateConstants(temp, iReac);

      /*--- Compute forward rate ---*/
      for_rate = 0.0;
      /*--- Check if we need to compute forward reaction rate in order to avoid useless and dangerous computations ---*/
      bool found_zero_mass_negative_exp = false;
      for(auto it = Reactant_Species_Negative_Exponent[iReac].cbegin(); it != Reactant_Species_Negative_Exponent[iReac].cend(); ++it) {
        if(Ys[*it] < 1.0e-15) {
           found_zero_mass_negative_exp = true;
           break;
        }
      }
      if(!found_zero_mass_negative_exp) {
        Eigen::ArrayXd stoich_reac_exp = Stoich_Coeffs_Reactants_Exp.row(iReac);
        auto cs_exp = cs.pow(stoich_reac_exp);
        for_rate = cs_exp.prod();
        for_rate *= Rate_Const.first;
      }

      /*--- Compute backward rate ---*/
      back_rate = 0.0;
      /*--- Check if we need to compute reverse reaction rate in order to avoid useless and dangerous computations ---*/
      found_zero_mass_negative_exp = false;
      for(auto it = Product_Species_Negative_Exponent[iReac].cbegin(); it != Product_Species_Negative_Exponent[iReac].cend(); ++it) {
        if(Ys[*it] < 1.0e-15) {
           found_zero_mass_negative_exp = true;
           break;
        }
      }
      if(!found_zero_mass_negative_exp) {
        Eigen::ArrayXd stoich_prod_exp = Stoich_Coeffs_Products_Exp.row(iReac);
        auto cs_prod = cs.pow(stoich_prod_exp);
        back_rate = cs_prod.prod();
        back_rate *= Rate_Const.second;
      }

      /*--- Save forward and backward rates ---*/
      Forward_Rates[iReac] = for_rate;
      Backward_Rates[iReac] = back_rate;
    }
  }

  //
  //
  /* This function computes the omega term. */
  RealVec ReactingModelLibrary::GetMassProductionTerm(const double temp, const double rho, const RealVec& ys) {
    /*--- Initialize to zero ---*/
    std::fill(Omega.begin(), Omega.end(), 0.0);

    /*--- Set forward and backward reaction rates ---*/
    SetReactionRates(temp, rho, ys);

    /*--- Compute mass production term ---*/
    for(unsigned short iReac = 0; iReac < nReactions; ++iReac) {
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        Omega[iSpecies] += (Stoich_Coeffs_Products(iSpecies,iReac) - Stoich_Coeffs_Reactants(iSpecies,iReac))*
                           (Forward_Rates[iReac] - Backward_Rates[iReac]);
      }
    }
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Omega[iSpecies] *= 1.0e-3*mMasses[iSpecies];

    return Omega;
  }

  //
  //
  /* This function computes the source chemistry Jacobian. NOTE: It requires SetReactionRates call*/
  RealMatrix ReactingModelLibrary::GetSourceJacobian(const double temp, const double rho) {
    /*--- Initialize to zero ---*/
    Source_Jacobian.setZero();

    /*--- Set Kc derivatives ---*/
    const double epsilon = 1.0e-6;
    double temp_pert = temp + epsilon*temp;
    const double RT = R_ungas*temp_pert;
    const double lnRT = std::log(R_ungas_atm*temp_pert);
    unsigned short iReac, iSpecies, jSpecies;
    double Kc_pert;

    for(iReac = 0; iReac < nReactions; ++iReac) {
      if(Available_Backward_Rate.count(iReac) == 0) {
        if(Backward_Rates[iReac] > 0.0) {
          double dG = 0.0;
          double dnu = 0.0;
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            double dcoeff = Stoich_Coeffs_Products(iSpecies,iReac) - Stoich_Coeffs_Reactants(iSpecies,iReac);
            if(dcoeff != 0.0) {
              dG += dcoeff*
                    (MathTools::GetSpline(std::get<T_DATA_SPLINE>(Enth_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Enth_Spline[iSpecies]),
                                          std::get<Y_DATA_SPLINE>(Enth_Spline[iSpecies]),temp_pert) - temp_pert*
                     MathTools::GetSpline(std::get<T_DATA_SPLINE>(Entr_Spline[iSpecies]),std::get<X_DATA_SPLINE>(Entr_Spline[iSpecies]),
                                          std::get<Y_DATA_SPLINE>(Entr_Spline[iSpecies]),temp_pert));
              dnu += dcoeff;
            }
          }
          double lnKc_pert = -dG/RT - dnu*lnRT;
          Kc_pert = std::exp(lnKc_pert);
        }
        else
          Kc_pert = Kc[iReac];
      }
      else {
        double kf_pert = As[iReac]*std::pow(temp_pert,Betas[iReac])*std::exp(-Temps_Activation[iReac]/temp_pert);
        double kb_pert = As_back[iReac]*std::pow(temp_pert,Betas_back[iReac])*std::exp(-Temps_Activation_back[iReac]/temp_pert);
        Kc_pert = kf_pert/kb_pert;
      }
      Kc_Derivatives[iReac] = (Kc_pert - Kc[iReac])/(temp_pert - temp);
    }

    /*--- Compute Jacobian ---*/
    for(iReac = 0; iReac < nReactions; ++iReac) {
      double tmp = (Betas[iReac] + Temps_Activation[iReac]/temp)/temp;
      double for_contr = Forward_Rates[iReac]*tmp;
      double back_contr;
      if(Available_Backward_Rate.count(iReac) == 0)
        back_contr = Backward_Rates[iReac]*(tmp - Kc_Derivatives[iReac]/Kc[iReac]);
      else
        back_contr = Backward_Rates[iReac]*(Betas_back[iReac] + Temps_Activation_back[iReac]/temp)/temp;

      /*--- Compute effectively the source Jacobian ---*/
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        /*--- Derivatives with respect to density ---*/
        double fixed_contr =  1.0e-3*mMasses[iSpecies]*
                              (Stoich_Coeffs_Products(iSpecies,iReac) - Stoich_Coeffs_Reactants(iSpecies,iReac));
        Source_Jacobian(iSpecies,0) += fixed_contr*(for_contr - back_contr);

        /*--- Derivatives with respect to partial densitiy ---*/
        for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
          if(Ys[jSpecies] > 1.0e-10)
            Source_Jacobian(iSpecies,jSpecies + 1) += fixed_contr*
                                                      (Forward_Rates[iReac]*Stoich_Coeffs_Reactants_Exp(iReac,jSpecies) -
                                                       Backward_Rates[iReac]*Stoich_Coeffs_Products_Exp(iReac,jSpecies))/(rho*Ys[jSpecies]);
        }
      }
    }

    return Source_Jacobian;
  }

  //
  //
  /* This function reads the species and sets their order. */
  void ReactingModelLibrary::ReadDataMixture(const std::string& f_name) {
    /*--- Local variables ---*/
    std::string line;
    std::string curr_species;
    unsigned short n_line = 0;

    std::ifstream mixfile(Lib_Path + "/" + f_name);
    if(mixfile.is_open()) {
      /*--- Clear vectors for safety ---*/
      Species_Names.clear();
      mMasses.clear();
      Formation_Enthalpies.clear();
      Diff_Volumes.clear();
      /*---- Read lines ---*/
      while(mixfile.good() && !mixfile.eof()) {
        std::getline(mixfile,line);
        /*--- Check if we encounter the termination character ---*/
        if(line == "STOP")
          break;
        /*--- We avoid clearly reading comments and empty lines in the file ---*/
        if(!line.empty() && !std::ispunct(line.at(0))) {
          if(n_line == 0) {
            std::istringstream curr_line(line);
            curr_line>>nSpecies;
            SU2_Assert(!curr_line.fail(), "You have to specify the number of species before proceding");
            SU2_Assert(nSpecies > 0, "The number of species must be greater than zero: you can't simulate anything otherwise");
            mMasses.reserve(nSpecies);
            Formation_Enthalpies.reserve(nSpecies);
            Diff_Volumes.reserve(nSpecies);
            n_line++;
          }
          else {
            SU2_Assert(std::isalpha(line.at(0)),std::string("Empty species field at line " + std::to_string(n_line)));
            std::string curr_species;
            double curr_mass,curr_enth,curr_vol;
            std::istringstream curr_line(line);
            curr_line>>curr_species;
            SU2_Assert(!curr_line.fail(), std::string("Empty species field at line " + std::to_string(n_line)));
            /*--- Read molar mass ---*/
            curr_line>>curr_mass;
            SU2_Assert(!curr_line.fail(), std::string("The molar mass of species " + curr_species + " is missing"));
            mMasses.push_back(curr_mass);
            /*--- Read formation enthalpy ---*/
            curr_line>>curr_enth;
            SU2_Assert(!curr_line.fail(), std::string("The formation enthalpy of species " + curr_species + " is missing"));
            Formation_Enthalpies.push_back(curr_enth);
            /*--- Read formation enthalpy ---*/
            curr_line>>curr_vol;
            SU2_Assert(!curr_line.fail(), std::string("The diffusion of species " + curr_species + " is missing"));
            Diff_Volumes.push_back(curr_vol);
            /*--- Insert in the map ---*/
            auto res = Species_Names.insert(std::make_pair(curr_species, n_line - 1));
            SU2_Assert(res.second == true, std::string("The species " + curr_species + " has already been declared"));
            n_line++;
          }
        }
      }
      mixfile.close();
      SU2_Assert(nSpecies > 0, "The number of species read is equal to zero: impossible");
      SU2_Assert(Species_Names.size() == nSpecies, "The number of species detected doesn't match nSpecies");

      /*--- Resize vector that wil be often used ---*/
      Ys.resize(nSpecies);
      Xs.resize(nSpecies);
      Cs.resize(nSpecies);

      Mu_Spline.resize(nSpecies);
      Kappa_Spline.resize(nSpecies);
      Entr_Spline.resize(nSpecies);
      Cp_Spline.resize(nSpecies);
      Enth_Spline.resize(nSpecies);

      Viscosities.resize(nSpecies);
      Enthalpies.resize(nSpecies);
      Internal_Energies.resize(nSpecies);
      CPs.resize(nSpecies);
      Thermal_Conductivities.resize(nSpecies);

      dPdYs.resize(nSpecies);
      ys_over_mm.resize(nSpecies);
      rhoUdiff.resize(nSpecies);
      Dm_coeffs.resize(nSpecies);
      Omega.resize(nSpecies, 0.0);
      Dij.resize(nSpecies,nSpecies);
      Dij.setZero();
      Gamma.resize(nSpecies,nSpecies);
      Gamma.setZero();
      Source_Jacobian.resize(nSpecies, nSpecies + 1);
      Source_Jacobian.setZero();
    }
    else {
      std::cerr<<"Unable to open the mixture file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Reading data about chemistry. ---*/
  void ReactingModelLibrary::ReadDataChem(const std::string& f_name) {
    /*--- Local variables ---*/
    std::string line;
    unsigned n_line = 0;
    unsigned n_reac_read = 0;

    std::ifstream chemfile(Lib_Path + "/" + f_name);
    if(chemfile.is_open()) {
      /*--- Clear for safety ---*/
      Stoich_Coeffs_Reactants.resize(0,0);
      Stoich_Coeffs_Reactants.resize(0,0);
      Stoich_Coeffs_Products_Exp.resize(0,0);
      Stoich_Coeffs_Reactants_Exp.resize(0,0);
      Reversible_Reactions.clear();
      As.clear();
      Betas.clear();
      Temps_Activation.clear();
      while(chemfile.good() && !chemfile.eof()) {
        std::getline(chemfile,line);

        /*--- Check if we encounter the termination character ---*/
        if(line == "STOP")
          break;
        /*--- We avoid clearly reading empty lines and comments in the file ---*/
        if(!line.empty() && !std::ispunct(line.at(0))) {
          if(n_line == 0) {
            std::istringstream curr_line(line);
            curr_line>>nReactions;
            SU2_Assert(!curr_line.fail(), "You have to specify the number of reactions before proceding");

            /*--- Resize and reserve space for vectors ---*/
            Stoich_Coeffs_Reactants.resize(nSpecies,nReactions);
            Stoich_Coeffs_Reactants.setZero();
            Stoich_Coeffs_Products.resize(nSpecies,nReactions);
            Stoich_Coeffs_Products.setZero();
            Stoich_Coeffs_Reactants_Exp.resize(nReactions,nSpecies);
            Stoich_Coeffs_Reactants_Exp.setZero();
            Stoich_Coeffs_Products_Exp.resize(nReactions,nSpecies);
            Stoich_Coeffs_Products_Exp.setZero();

            Forward_Rates.resize(nReactions);
            Backward_Rates.resize(nReactions);
            Kc.resize(nReactions);
            Kc_Derivatives.resize(nReactions);

            Reversible_Reactions.reserve(nReactions);
            As.reserve(nReactions);
            Betas.reserve(nReactions);
            Temps_Activation.reserve(nReactions);

            n_line++;
          }
          else if(n_line == 1) {
            std::istringstream curr_line(line);
            std::string kind_units;
            curr_line>>kind_units;
            if(kind_units == "CGS")
              CGS_Units = true;
            else if(kind_units == "SI")
              CGS_Units = false;
            else
              throw std::out_of_range("Unknown option for the type of units measure");

            n_line++;
          }
          else {
            bool is_rev;
            if((n_line % 2 == 0) && (n_line < 2*nReactions + 1)) {
              is_rev = (line.find('<') != std::string::npos);
              Reversible_Reactions.push_back(is_rev);
              n_reac_read++;
              ReadReactSpecies(line,is_rev,n_reac_read);
            }
            else if((n_line % 2 == 1) && (n_line < 2*nReactions + 2)) {
              ReadChemCoefs(line);
            }
            else {
              ReadExtraData_Rates(line);
              ReadExtraData_ForwardExponent(line);
              ReadExtraData_BackwardExponent(line);
            }
            n_line++;
          }
        }
      }
      SU2_Assert(n_reac_read == nReactions, "The number of reactions detected doesn't match nReactions");
      chemfile.close();
      unsigned short iReac, iSpecies;

      /*--- Try automatic computations of exponents of products in case backward data were not already available ---*/
      for(iReac = 0; iReac < nReactions; ++iReac) {
        if(Reversible_Reactions[iReac] && Available_Backward_Rate.count(iReac) == 0) {
          for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
            Stoich_Coeffs_Products_Exp(iReac,iSpecies) = Stoich_Coeffs_Reactants_Exp(iReac,iSpecies) +
                                                         Stoich_Coeffs_Products(iSpecies,iReac) - Stoich_Coeffs_Reactants(iSpecies,iReac);
        }
      }

      /*--- Update to SI units the Arrhenius constants if needed ---*/
      if(CGS_Units) {
        for(iReac = 0; iReac < nReactions; ++ iReac) {
          double sum_forward_exp = Stoich_Coeffs_Reactants_Exp.row(iReac).sum();
          As[iReac] *= std::pow(10.0, 6.0*(1.0 - sum_forward_exp));
          if(Available_Backward_Rate.count(iReac) == 1) {
            double sum_backward_exp = Stoich_Coeffs_Products_Exp.row(iReac).sum();
            As_back[iReac] *= std::pow(10.0, 6.0*(1.0 - sum_backward_exp));
          }
        }
      }

      /*--- Save species with negative exponents for forward  rates---*/
      Reactant_Species_Negative_Exponent.resize(nReactions);
      for(iReac = 0; iReac < nReactions; ++iReac) {
         for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            if(Stoich_Coeffs_Reactants_Exp(iReac,iSpecies) < 0.0)
               Reactant_Species_Negative_Exponent[iReac].push_back(iSpecies);
         }
      }

      /*--- Save species with negative exponents for backward rates---*/
      Product_Species_Negative_Exponent.resize(nReactions);
      for(iReac = 0; iReac < nReactions; ++iReac) {
         for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            if(Stoich_Coeffs_Products_Exp(iReac,iSpecies) < 0.0)
               Product_Species_Negative_Exponent[iReac].push_back(iSpecies);
         }
      }

    }
    else {
      std::cerr<<"Unable to open the chemical file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- This function reads chemical reactions and store stoichiometric coefficients. ---*/
  void ReactingModelLibrary::ReadReactSpecies(const std::string& line, bool is_rev, unsigned n_reac) {
    /*--- Check correct reaction format ---*/
    auto minor_pos = line.find('<');
    auto major_pos = line.find('>');
    SU2_Assert(major_pos != std::string::npos, "No reaction in this line");
    SU2_Assert(line.find('>',major_pos+1), "Already detected > symbol for reactions");

    std::string reactants_side, products_side;

    if(is_rev) {
      SU2_Assert(minor_pos + 2 == major_pos, "Incorrect symbol to detect reactions");
      SU2_Assert(line.find('<',minor_pos + 1), "Already detected < symbol for reactions");
      reactants_side = line.substr(0,minor_pos);
    }
    else {
      auto bar_pos = line.find('=');
      SU2_Assert(bar_pos != std::string::npos, "No reaction in this line");
      SU2_Assert(bar_pos + 1 == major_pos, "Incorrect symbol to detect reactions");
      reactants_side = line.substr(0,bar_pos);
    }
    products_side = line.substr(major_pos + 1);

    Utility::Parse_Terms(reactants_side, n_reac, is_rev, true, Species_Names, Stoich_Coeffs_Reactants,
                         Stoich_Coeffs_Reactants_Exp, Stoich_Coeffs_Products_Exp);
    Utility::Parse_Terms(products_side, n_reac, is_rev, false, Species_Names, Stoich_Coeffs_Products,
                         Stoich_Coeffs_Reactants_Exp, Stoich_Coeffs_Products_Exp);
  }

  //
  //
  /*--- This function reads coefficients for equilibrium constants. ---*/
  void ReactingModelLibrary::ReadChemCoefs(const std::string& line) {
    /*--- Local variables ---*/
    double A, beta, Ta;

    std::istringstream curr_line(line);
    /*--- Reading pre-exponential factor ---*/
    curr_line>>A;
    SU2_Assert(!curr_line.fail(), "No exponential prefactor after a reaction");
    As.push_back(A);
    /*--- Reading temperature exponent ---*/
    curr_line>>beta;
    SU2_Assert(!curr_line.fail(), "No temperature exponent after a reaction");
    Betas.push_back(beta);
    /*--- Reading temperature activation ---*/
    curr_line>>Ta;
    SU2_Assert(!curr_line.fail(), "No activation temperature after a reaction");
    if(CGS_Units)
      Temps_Activation.push_back(Ta/R_ungas_scal);
    else
      Temps_Activation.push_back(Ta);
  }

  //
  //
  /*--- This function reads data of backward rates in case they are already available ---*/
  void ReactingModelLibrary::ReadExtraData_Rates(std::string line) {
    /*--- Find if the chosen string to detect data is present ---*/
    auto position = line.find("Available Backward Rate reaction",0);

    if(position != std::string::npos) {
      /*--- Start reading data after the detecting string ---*/
      line.erase(0,32);

      std::istringstream curr_line(line);
      unsigned short iReac;
      double A_back, beta_back, Ta_back;
      curr_line>>iReac;
      if(Reversible_Reactions[iReac - 1]) {
        if(Available_Backward_Rate.count(iReac - 1) == 0)
          Available_Backward_Rate[iReac - 1] = true;
        else
          throw std::runtime_error("Manual backward rate for reaction " + std::to_string(iReac) + " already read");
      }
      else
        throw std::runtime_error("You are trying to read the coefficients to compute the backward rate of a irreversible reaction");
      line.erase(0,3);
      curr_line.clear();
      curr_line.str(line);

      /*--- Reading Arrhenius constant ---*/
      curr_line>>A_back;
      SU2_Assert(!curr_line.fail(), "No exponential prefactor of the manual backward rate for reaction " + std::to_string(iReac));
      As_back[iReac - 1] = A_back;

      /*--- Reading temperature exponent ---*/
      curr_line>>beta_back;
      SU2_Assert(!curr_line.fail(), "No temperature exponent of the manual backward rate for reaction " + std::to_string(iReac));
      Betas_back[iReac - 1] = beta_back;

      /*--- Reading temperature activation ---*/
      curr_line>>Ta_back;
      SU2_Assert(!curr_line.fail(), "No activation temperature of the manual backward rate for reaction " + std::to_string(iReac));
      if(CGS_Units)
        Temps_Activation_back[iReac - 1] = Ta_back/R_ungas_scal;
      else
        Temps_Activation_back[iReac - 1] = Ta_back;
    }
  }

  //
  //
  /*--- This function reads data of extra species to include in the computation of forward rates. ---*/
  void ReactingModelLibrary::ReadExtraData_ForwardExponent(std::string line) {
    /*--- Find if the chosen string to detect data is present ---*/
    auto position = line.find("Extra Forward terms reaction",0);

    if(position != std::string::npos) {
      /*--- Start reading data after the detecting string ---*/
      line.erase(0,28);

      std::istringstream curr_line(line);
      unsigned short iReac;
      curr_line>>iReac;
      line.erase(0,3);
      curr_line.clear();
      curr_line.str(line);

      /*--- Read species and exponents involved ---*/
      Utility::Parse_ExtraTerms(line, iReac - 1, Species_Names, Stoich_Coeffs_Reactants_Exp);
    }
  }

  //
  //
  /*--- This function reads data of extra species to include in the computation of backward rates. ---*/
  void ReactingModelLibrary::ReadExtraData_BackwardExponent(std::string line) {
    /*--- Find if the chosen string to detect data is present ---*/
    auto position = line.find("Extra Backward terms reaction",0);

    if(position != std::string::npos) {
      /*--- Start reading data after the detecting string ---*/
      line.erase(0,29);

      std::istringstream curr_line(line);
      unsigned short iReac;
      curr_line>>iReac;
      line.erase(0,3);
      curr_line.clear();
      curr_line.str(line);

      /*--- Read species and exponents involved ---*/
      Utility::Parse_ExtraTerms(line, iReac - 1, Species_Names, Stoich_Coeffs_Products_Exp);
    }
  }

  //
  //
  /*--- Reading data about transport properties. ---*/
  void ReactingModelLibrary::ReadDataTransp(const std::string& f_name) {
    /*--- Local variables ---*/
    std::string line;
    unsigned n_line = 0;
    std::string curr_species;
    unsigned short iSpecies;
    double curr_temp, curr_visc, curr_cond;
    RealVec temp_data, mu_data, kappa_data;

    std::ifstream transpfile(Lib_Path + "/" + f_name);

    if(transpfile.is_open()) {
      while(transpfile.good() && !transpfile.eof()) {
        std::getline(transpfile,line);
        /*--- We avoid clearly reading empty lines and comments ---*/
        if(!line.empty() && !std::ispunct(line.at(0))) {
          if(n_line == 0) {
            SU2_Assert(std::isalpha(line.at(0)), "You have to specify the species");
            auto it = Species_Names.find(line);
            SU2_Assert(it != Species_Names.end(), "The species is not present in the mixture");
            curr_species = it->first;
            iSpecies = it->second;
            n_line++;
          }
          else {
            std::istringstream curr_line(line);

            /*--- Reading temperature ---*/
            curr_line>>curr_temp;
            SU2_Assert(!curr_line.fail(), std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            temp_data.push_back(curr_temp);

            /*--- Reading viscosity ---*/
            curr_line>>curr_visc;
            SU2_Assert(!curr_line.fail(), std::string("Empty Mu field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            mu_data.push_back(curr_visc);

            /*--- Reading conductivity ---*/
            curr_line>>curr_cond;
            SU2_Assert(!curr_line.fail(), std::string("Empty Kappa field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            kappa_data.push_back(curr_cond);

            n_line++;
          } /*--- End reading data --*/
        }
      } /*--- End reading file ---*/

      RealVec y2_mu, y2_kappa;

      MathTools::SetSpline(temp_data, mu_data, 0.0, 0.0, y2_mu);
      Mu_Spline[iSpecies] = std::make_tuple(temp_data, std::move_if_noexcept(mu_data), std::move_if_noexcept(y2_mu));

      MathTools::SetSpline(temp_data, kappa_data, 0.0, 0.0, y2_kappa);
      Kappa_Spline[iSpecies] = std::make_tuple(temp_data, std::move_if_noexcept(kappa_data), std::move_if_noexcept(y2_kappa));

      transpfile.close();
    }
    else {
      std::cerr<<"Unable to open the species file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Reading data around thermodynamic properties. ---*/
  void ReactingModelLibrary::ReadDataThermo(const std::string& f_name) {
    /*--- Local variables ---*/
    std::string line;
    unsigned n_line = 0;
    std::string curr_species;
    unsigned short iSpecies;
    double curr_temp, curr_enth, curr_Cp, curr_entr;
    RealVec temp_data, cp_data, enth_data, entr_data;

    std::ifstream thermofile(Lib_Path + "/" + f_name);

    if(thermofile.is_open()) {
      while(thermofile.good() && !thermofile.eof()) {
        std::getline(thermofile,line);
        /*--- We clearly avoid reading and empty lines ---*/
        if(!line.empty() && !std::ispunct(line.at(0))) {
          if(n_line == 0) {
            SU2_Assert(std::isalpha(line.at(0)), "You have to specify the species");
            auto it = Species_Names.find(line);
            SU2_Assert(it != Species_Names.end(), "The species is not present in the mixture");
            curr_species = it->first;
            iSpecies = it->second;
            n_line++;
          }
          else {
            std::istringstream curr_line(line);

            /*--- Reading temperature ---*/
            curr_line>>curr_temp;
            SU2_Assert(!curr_line.fail(), std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            temp_data.push_back(curr_temp);

            /*--- Reading Cp ---*/
            curr_line>>curr_Cp;
            SU2_Assert(!curr_line.fail(), std::string("Empty Cp field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            cp_data.push_back(curr_Cp);

            /*--- Reading Enthalpy ---*/
            curr_line>>curr_enth;
            SU2_Assert(!curr_line.fail(), std::string("Empty Enthalpy field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            enth_data.push_back(curr_enth);

            /*--- Reading Entropy ---*/
            curr_line>>curr_entr;
            SU2_Assert(!curr_line.fail(), std::string("Empty Entropy field at line " + std::to_string(n_line + 1) +
                                                      " for species " + curr_species));
            entr_data.push_back(curr_entr);

            n_line++;
          } /*--- End reading data ---*/
        }
      } /*--- End reading file ---*/

      RealVec y2_cp, y2_enth, y2_entr;

      MathTools::SetSpline(temp_data, cp_data, 0.0, 0.0, y2_cp);
      Cp_Spline[iSpecies] = std::make_tuple(temp_data, std::move_if_noexcept(cp_data), std::move_if_noexcept(y2_cp));

      MathTools::SetSpline(temp_data, enth_data, 0.0, 0.0, y2_enth);
      Enth_Spline[iSpecies] = std::make_tuple(temp_data, std::move_if_noexcept(enth_data), std::move_if_noexcept(y2_enth));

      MathTools::SetSpline(temp_data, entr_data, 0.0, 0.0, y2_entr);
      Entr_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(entr_data),std::move_if_noexcept(y2_entr));

      thermofile.close();
    }
    else {
      std::cerr<<"Unable to open the thermo file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /* This function reads the physical data of the fuel in case of regression boundary condition. */
  void ReactingModelLibrary::ReadDataFuel(const std::string& f_name) {
    /*--- Local variables ---*/
    std::string line;
    unsigned short n_line = 0;
    std::array<bool,5> set_data;
    std::fill(set_data.begin(), set_data.end(), false);
    std::array<std::size_t,5> found_pos;
    unsigned short n_found, elem_found;

    std::ifstream fuel_file(Lib_Path + "/" + f_name);
    if(fuel_file.is_open()) {
      /*---- Read lines ---*/
      while(fuel_file.good() && !fuel_file.eof()) {
        std::getline(fuel_file,line);
        /*--- Check if we encounter the termination character ---*/
        if(line == "STOP")
          break;
        /*--- We avoid clearly reading comments and empty lines in the file ---*/
        if(!line.empty() && !std::ispunct(line.at(0))) {
          found_pos[A1_INDEX] = line.find("A1   = ",0);
          found_pos[A2_INDEX] = line.find("A2   = ",0);
          found_pos[EA1_INDEX] = line.find("EA1  = ",0);
          found_pos[EA2_INDEX] = line.find("EA2  = ",0);
          found_pos[TBAR_INDEX] = line.find("Tbar = ",0);

          n_found = 0;
          for(unsigned short iData = 0; iData < 5; ++iData) {
            if(found_pos[iData] != std::string::npos) {
              n_found++;
              elem_found = iData;
            }
          }
          SU2_Assert(n_found == 1, "Multiple data read on single line");
          SU2_Assert(set_data[elem_found] == false, "Detected data has already been read. Check the file content");
          line.erase(0,7);
          Fuel_Data[elem_found] = std::stod(line);
          set_data[elem_found] = true;
          n_line++;
        }
      }
      fuel_file.close();
      auto it = std::find(set_data.begin(), set_data.end(), false);
      SU2_Assert(it == set_data.end(), "Some fuel properties has not been setted. Check file content");
    }
    else {
      std::cerr<<"Unable to open the fuel data file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Compute fuel regression rate at specified temperature ---*/
  inline double ReactingModelLibrary::ComputeRegressionRate(const double temp) {
    if(temp < Fuel_Data[TBAR_INDEX])
      return Fuel_Data[A2_INDEX]*std::exp(Fuel_Data[EA2_INDEX]/(R_ungas_scal*temp));
    return Fuel_Data[A1_INDEX]*std::exp(Fuel_Data[EA1_INDEX]/(R_ungas_scal*temp));
  }

  //
  //
  /*--- Setup library ---*/
  void ReactingModelLibrary::Setup(void) {
    if(!Lib_Setup) {
      Le = 1.0;
      /*--- If nobody has configured the library path, we try to do it here with a default value ---*/
      if(Lib_Path == "") {
        std::cout<<"Library path set to default"<<std::endl;
        auto base_dir = std::experimental::filesystem::current_path().string();
        Lib_Path = base_dir;
      }

      std::vector<std::string> list_file;
      std::ifstream config_file(Config_File);
      if(config_file.is_open()) {
        while(config_file.good() && !config_file.eof()) {
          std::string curr_line;
          std::getline(config_file,curr_line);
          if(!curr_line.empty() && !std::ispunct(curr_line.at(0)))
            list_file.push_back(curr_line);
        }
      }
      else {
        std::cerr<<"Unable to open the specified file with all the file names for setting library."<<std::endl;
        std::exit(1);
      }

      /*--- Read mixture file: it needs to be the first to check exactness of chemical reactions and properties ---*/
      std::string file_mix = list_file.at(0);
      ReadDataMixture(file_mix);
      std::cout<<"Mixture Data read"<<std::endl;

      /*--- Check we have the right number of files ---*/
      using size_type = std::vector<std::string>::size_type;
      size_type max_n_file = 2*nSpecies + 2;
      size_type n_file = list_file.size();
      SU2_Assert((n_file == max_n_file) || (n_file == max_n_file - 1), "The number of files present in the configuration file is wrong");

      /*--- Set the specific gas constants ---*/
      SetRiGas();

      /*--- Read chemistry file (if present) ---*/
      int buffer_chemistry = 1;
      nReactions = 0;
      if(n_file == max_n_file) {
        /*--- We assume that chemistry file is the second in the list if its present ---*/
        std::string file_chem = list_file[1];
        ReadDataChem(file_chem);
        std::cout<<"Chemical Reactions read"<<std::endl;
        buffer_chemistry = 0;
      }

      /*--- We assume that data correspond to the species declared at the beginning of the file
            and a transport file is followed by a thermodynamic file
            (I can't check the content so it seems reasonable) ---*/
      std::string file_transp, file_thermo;
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        file_transp = list_file[iSpecies*2 + 2 - buffer_chemistry];
        ReadDataTransp(file_transp);
        file_thermo = list_file[iSpecies*2 + 3 - buffer_chemistry];
        ReadDataThermo(file_thermo);
      }
      Lib_Setup = true;
      std::cout<<"Library set."<<std::endl;
      std::cout<<std::endl;
   }
   else
      throw Common::NotSetup("Trying to setup again without calling unsetup first.");
  }

  //
  //
  /*--- Unsetup library ---*/
  void ReactingModelLibrary::Unsetup(void) {
    if(Lib_Setup) {
      Species_Names.clear();

      Ri.clear();
      mMasses.clear();
      Diff_Volumes.clear();
      Ys.clear();
      Xs.clear();
      Cs.clear();
      Viscosities.clear();
      Enthalpies.clear();
      Internal_Energies.clear();
      CPs.clear();
      Thermal_Conductivities.clear();
      Formation_Enthalpies.clear();

      Stoich_Coeffs_Products.resize(0,0);
      Stoich_Coeffs_Products_Exp.resize(0,0);
      Stoich_Coeffs_Reactants.resize(0,0);
      Stoich_Coeffs_Reactants_Exp.resize(0,0);
      Forward_Rates.clear();
      Backward_Rates.clear();
      Kc.clear();
      Kc_Derivatives.clear();
      As.clear();
      Betas.clear();
      Temps_Activation.clear();
      Reversible_Reactions.clear();

      Enth_Spline.clear();
      Cp_Spline.clear();
      Mu_Spline.clear();
      Kappa_Spline.clear();
      Entr_Spline.clear();

      ys_over_mm.clear();
      rhoUdiff.clear();
      Dm_coeffs.clear();
      Omega.clear();
      Dij.resize(0,0);
      Gamma.resize(0,0);

      Available_Backward_Rate.clear();
      As_back.clear();
      Betas_back.clear();
      Temps_Activation_back.clear();
      Reactant_Species_Negative_Exponent.clear();
      Product_Species_Negative_Exponent.clear();

      Lib_Setup = false;
      std::cout<<"Library unset."<<std::endl;
      std::cout<<std::endl;
    }
    else
      throw Common::NotSetup("Trying to unsetup without calling setup first.");
  }

} /*--- End of namespace Common ---*/
