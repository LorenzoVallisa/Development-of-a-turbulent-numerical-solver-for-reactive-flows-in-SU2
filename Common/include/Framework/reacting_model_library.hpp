#ifndef SU2_REACTING_MODEL_LIBRARY
#define SU2_REACTING_MODEL_LIBRARY

#include "physical_chemical_library.hpp"

#include "su2_assert.hpp"
#include <map>
#include <set>
#include <tuple>

namespace Framework {

  /*!
   * \brief Provides a particular library definition to compute the physical and chemical properties.
   */
  class ReactingModelLibrary: public Framework::PhysicalChemicalLibrary {

  public:
    typedef std::tuple<RealVec,RealVec,RealVec> MyTuple;
    using Vec = Eigen::VectorXd;

  public:
    /*!
     * \brief Constructor with the name of the configuration file of the library and the library path.
     */
    ReactingModelLibrary(const std::string& conf_name, const std::string& lib_path_name):
    PhysicalChemicalLibrary(conf_name, lib_path_name), Rgas(), Le() {}

    /*!
     * \brief Default destructor.
     */
    virtual ~ReactingModelLibrary() {}

    /*!
     * \brief Check if library is setup.
     */
    inline bool IsSetup(void) const {
      return PhysicalChemicalLibrary::Lib_Setup;
    }

    /*!
     * \brief Setups the data of the library.
     */
    void Setup(void) override;

    /*!
     * \brief Unsetups the data of the library.
     */
    void Unsetup(void) override;

    /*!
     * \brief Get the index of a species.
     * \param[in] name_species - Name of the desired species
     */
    inline unsigned short GetIndexSpecies(const std::string& name_species) const {
      auto it = Species_Names.find(name_species);
      SU2_Assert(it != Species_Names.end(), "The species " + name_species + " which you want the index is not present in the mixture");
      return it->second;
    }

    /*!
     * \brief Set the gas constant for each species [J/(Kg*K)]
     */
    void SetRiGas(void) override;

    /*!
     * \brief Set the gas constant for each species [J/(Kg*K)]
     */
    RealVec GetRiGas(void) const override {
      return Ri;
    }

    /*!
     * Get the gas constant for a desired species [J/(Kg*K)]
     * \param[in] iSpecies - index of the desired species;
    */
    inline double GetRiGas(unsigned short iSpecies) const {
      return Ri.at(iSpecies);
    }

    /*!
     * \brief Get the gas constant for the mixture [J/(Kg*K)]
     */
    inline double GetRgas(void) const override {
      return Rgas;
    }

    /*!
     * \brief Set the gas constant for the mixture[J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    void SetRgas(const RealVec& ys) override;

    /*!
     * Compute the gas constant for the mixture [J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    double ComputeRgas(const RealVec& ys) override;

    /*!
     * \brief Get the molar masses of the species
     */
    inline RealVec GetMolarMasses(void) const override {
      SU2_Assert(mMasses.size() == nSpecies,"The number of elements in the vector of molar masses doesn't match nSpecies");
      return mMasses;
    }

    /*!
     * \brief Get the molar fractions from the mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    RealVec GetMolarFromMass(const RealVec& ys) override;

    /*!
     * Set the molar fractions of elements Xs. This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] xs - The vector of the molar fractions of elements
    */
    void SetMolarFractions(const RealVec& xs) override;

    /*!
     * \brief Set the molar fractions from mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    void SetMolarFromMass(const RealVec& ys) override;

    /*!
     * \brief Get the mass fractions from the molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
     */
    RealVec GetMassFromMolar(const RealVec& xs) override;

    /*!
     * \brief Set the mass fractions. This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] ys The vector of the mass fractions of species
    */
    void SetMassFractions(const RealVec& ys) override;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
    */
    void SetMassFromMolar(const RealVec& xs) override;

    /*!
     * \brief Computes the frozen specific heat ratio and the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[out] gamma - specific heat ratio (output)
     * \param[out] sound_speed - speed of sound (output)
    */
    void Gamma_FrozenSoundSpeed(const double temp, const RealVec& ys, double& gamma, double& sound_speed) override;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     */
    double ComputeFrozenGamma(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
    */
    double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] gamma - specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
    */
    double ComputeFrozenSoundSpeed_FromGamma(const double temp, const double gamma, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] sound_speed - frozen speed of sound
     * \param[in] ys - The vector of the mass fractions of species
    */
    double ComputeFrozenGamma_FromSoundSpeed(const double temp, const double sound_speed, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys, const double press, const double rho) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] gamma - specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    double ComputeFrozenSoundSpeed_FromGamma(const double gamma, const RealVec& ys, const double press, const double rho) override;

    /*!
     * \brief Computes the density, the enthalpy and the internal energy
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
     */
    void Density_Enthalpy_Energy(const double temp, const double pressure, const RealVec& ys, RealVec& dhe) override;

    /*!
     * \brief Computes the pressure  at given temperature and density.
     * \param[in] temp - temperature
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    double ComputePressure(const double temp, const double rho, const RealVec& ys) override;

    /*!
     * \brief Computes the density at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
    */
    double ComputeDensity(const double temp, const double pressure, const RealVec& ys) override;

    /*!
     * \brief Compute the temperature given density and pressure.
     * \param[in] pressure - pressure
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    double ComputeTemperature(const double pressure, const double rho, const RealVec& ys) override;

    /*!
     * \brief Compute the density given sound speed and gamma.
     * \param[in] sound_speed2 - squared sound speed
     * \param[in] gamma - specific heats ratio
     * \param[in] ys - mass fractions
    */
    double ComputeTemperature_FromGamma(const double sound_speed2, const double gamma, const RealVec& ys) override;

    /*!
     * \brief Computes the internal energy per unit of mass at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     */
    double ComputeEnergy(const double temp, const RealVec& ys) override;

    /*!
     * \brief Return the formation enthalpies per unit mass of species
    */
    inline RealVec GetFormationEnthalpies(void) const override {
      return Formation_Enthalpies;
    }

    /*!
     * \brief Return the static enthalpy per unit of mass
     * \param[in] temp - the mixture temperature
     * \param[in] ys - mass fractions
     * \return Mixture static enthalpy
    */
    double ComputeEnthalpy(const double temp, const RealVec& ys) override;

    /*!
     * \brief Set the static enthalpy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    void SetPartialEnthalpy(const double temp) override;

    /*!
     * \brief Set the static enthalpy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    void SetPartialEnthalpy(const double temp, unsigned short iSpecies) override;

    /*!
     * \brief Return the static enthalpy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    RealVec ComputePartialEnthalpy(const double temp) override;

    /*!
     * \brief Compute the static enthalpy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    double ComputePartialEnthalpy(const double temp, unsigned short iSpecies) override;

    /*!
     * \brief Set the internal energy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    void SetPartialEnergy(const double temp) override;

    /*!
     * \brief Set the internal energy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    void SetPartialEnergy(const double temp, unsigned short iSpecies) override;

    /*!
     * \brief Return the internal energy of each species.
     * \param[in] temp - the mixture temperature
    */
    RealVec ComputePartialEnergy(const double temp) override;

    /*!
     * \brief Compute the internal energy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    double ComputePartialEnergy(const double temp, unsigned short iSpecies) override;

    /*!
     * \brief Computes the pressure derivative w.r.t. partial densities.
     * \param[in] temp - temperature
     * \param[in] gamma - frozen specific heat ratio
    */
    RealVec ComputedP_dYs(const double temp, const double gamma) override;

    /*!
     * \brief Compute the pressure derivative w.r.t. partial density for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] gamma - frozen specific heat ratio
     * \param[in] iSpecies - index of desired species
    */
    double ComputedP_dYs(const double temp, const double gamma, unsigned short iSpecies) override;

    /*!
     * \brief Set the actual concetration for each species.
     * \param[in] rho - the mixture density
     * \param[in] ys - the actual mass fractions
    */
    void SetConcentration(const double rho, const RealVec& ys) override;

    /*!
     * \brief Computes the mixture total concentration
     * \param[in] rho - density
     */
    double ComputeConcentration(const double rho, const RealVec& ys) override;

    /*!
     * \brief Return the specific heat at constant pressure per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    RealVec ComputeCps(const double temp) override;

    /*!
     * \brief Computes the specific heat at constant pressure
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Specific heat at constant pressure
    */
    double ComputeCP(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Specific heat at constant volume
    */
    double ComputeCV(const double temp, const RealVec& ys) override {
      return ComputeCP(temp,ys) - ComputeRgas(ys);
    }

    /*!
     * \brief Compute the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] cp - specific heat at constant pressure
     * \return Cv - specific heat at constant volume
    */
    inline double ComputeCV_FromCP(const double cp, const RealVec& ys) override {
      return cp - ComputeRgas(ys);
    }

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] cp - specific heat at constant pressure
     * \param[in] ys - The vector of the mass fractions of species
     */
    inline double ComputeFrozenGamma_FromCP(const double cp, const RealVec& ys) override {
      return cp/(cp - ComputeRgas(ys));
    }

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - mxiture temperature;
     * \param[in] cp - specific heat at constant pressure
     * \param[in] ys - The vector of the mass fractions of species
     */
    inline double ComputeFrozenSoundSpeed_FromCP(const double temp, const double cp, const RealVec& ys) override {
      double gamma = ComputeFrozenGamma_FromCP(cp, ys);
      return std::sqrt(gamma*Rgas*temp);
    }

    /*!
     * \brief Computes the specific heat at constant pressure.
     * \param[in] temp - temperature
     * \param[in] sound_speed - frozen speed of sound
     * \param[in] ys - The vector of the mass fractions of species
    */
    inline double ComputeCP_FromSoundSpeed(const double temp, const double sound_speed, const RealVec& ys) override {
      SetRgas(ys);
      return (sound_speed*sound_speed*Rgas)/(sound_speed*sound_speed - Rgas*temp);
    }

    /*!
     * \brief Computes the specific heat at constant pressure.
     * \param[in] gamma - frozen specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
    */
    inline double ComputeCP_FromGamma(const double gamma, const RealVec& ys) override {
      return gamma*ComputeRgas(ys)/(gamma - 1.0);
    }

    /*!
     * \brief Computes the thermal conductivity for each species at given temperature
     * \param[in] temp - temperature
     * \return Thermal conductivity for each species
    */
    void ComputeConductivities(const double temp);

    /*!
     * \brief Computes the thermal conductivity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Mixture thermal conductivity
    */
    double ComputeLambda(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the dynamic viscosity for each species at given temperature
     * \param[in] temp - temperature
     * \return Molecular viscosity for each species
    */
    void ComputeViscosities(const double temp);

    /*!
     * \brief Computes the dynamic viscosity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Molecular viscosity for the mixture in thermal non equilibrium
    */
    double ComputeEta(const double temp, const RealVec& ys) override;

   /*!
    * Return the diffusion velocities of species multiplied by the species
    * densities for nonequilibrium computations
    * \param[in] temp - the mixture temperature
    * \param[in] rho  - the mixture density
    * \param[in] ys - the species mass fractions
    * \return Diffusion coefficient with constant Lewis number for each species
    */
    RealVec GetRhoUdiff(const double temp, const double rho, const RealVec& ys) override;

   /*!
    * Return the mass production/destruction terms [kg m^-3 s^-1] in chemical
    * non-equilibrium based on Arrhenius's formula.
    * \param[in] temp - the mixture temperature
    * \param[in] rho - the mixture density
    * \param[in] ys - the species mass fractions
    * \return Mass production terms
    */
    RealVec GetMassProductionTerm(const double temp, const double rho, const RealVec& ys) override;

    /*!
     * Compute the Jacobian of source chemistry. NOTE: It requires SetReactionRates call
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \return Contribution to source derivatives with respect to mixture density and partial densitiies
     */
    RealMatrix GetSourceJacobian(const double temp, const double rho) override;

    /*!
     * \brief Return the effective diffusion coefficients to solve Stefan-Maxwell equation
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[in] ys - mass fractions in the mixture
     */
    RealVec GetDiffCoeffs(const double temp, const double pressure, const RealVec& ys) override;

   /*!
    * \brief Return the binary diffusion coefficients
    * \param[in] temp - the mixture temperature
    * \param[in] pressure - the mixture pressure
    */
    RealMatrix GetDij_SM(const double temp, const double pressure) override;

    /*!
     * \brief Return the matrix of Stefan-Maxwell equations
     * \param[in] rho - the mixture density
     * \param[in] xs - current molar fractions
     * \param[in] ys - current mass fractions
     * \param[in] val_Dij - current binary diffusion coefficients
     */
    RealMatrix GetGamma(const double rho, const RealVec& xs, const RealVec& ys, const RealMatrix& val_Dij) override;

    /*!
     * \brief Compute the regression rate at specified temperature with an empirical law
     * \param[in] temp - the mixture temperature
     */
    double ComputeRegressionRate(const double temp) override;

    /*!
     * \brief Read all physical data about fuel
     * \param[in] f_name - File with the fuel data to be read
     */
    void ReadDataFuel(const std::string& f_name) override;

  private:
    /*!
     * \brief Return the equilibrium constants (concentration and pressure).
     * \param[in] temp - the mixture temperature
     * \param[in] iReac - index of the desired reaction
    */
    std::pair<double,double> ComputeKeq(const double temp, unsigned short iReac);

    /*!
     * \brief Return the forward and backward reaction rate constants.
     * \param[in] temp - the mixture temperature
     * \param[in] iReac - index of the desired reaction
    */
    std::pair<double,double> ComputeRateConstants(const double temp, unsigned short iReac);

    /*!
     * \brief Set forward and backward reaction rates for each reaction.
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
    */
    void SetReactionRates(const double temp, const double rho, const RealVec& ys);

    /*!
      * \brief Read transport data
      * \param[in] f_name - name of the file with the properties
    */
    void ReadDataTransp(const std::string& f_name);

    /*!
      * \brief Read thermodynamical data
      * \param[in] f_name - name of the file with the properties
    */
    void ReadDataThermo(const std::string& f_name);

    /*!
      * \brief Read mixture data
      * \param[in] f_name - name of the file with the involved species
      */
    void ReadDataMixture(const std::string& f_name);

    /*!
    * \brief Read chemistry data
    * \param[in] f_name - name of the file with the involved reactions
    */
    void ReadDataChem(const std::string& f_name);

    /*!
      * \brief Read species involved in a chemical reaction from line
      * \param[in] line - line to be read
      * \param[in] is_elem - flag whether the reaction is elementary or not
      * \parm[in]  n_reac - the cuurent number of reactions detected
      * \param[out] Species_Reactions - names of species form reactions
    */
    void ReadReactSpecies(const std::string& line, bool is_elem, unsigned n_reac);

    /*!
      * \brief Read coefficients to compute reaction rates from line
      * \param[in] line - line to be read
    */
    void ReadChemCoefs(const std::string& line);

    /*!
      * \brief Read extra data in case backward rates are already available for a reversible reaction.
      * \param[in] line - line to be read
    */
    void ReadExtraData_Rates(std::string line);

    /*!
      * \brief Read extra exponent for forward rate to compute coorectly the rates.
      * \param[in] line - line to be read
    */
    void ReadExtraData_ForwardExponent(std::string line);

    /*!
      * \brief Read extra exponent for backward rate to compute coorectly the rates.
      * \param[in] line - line to be read
    */
    void ReadExtraData_BackwardExponent(std::string line);

  protected:

    std::map<std::string,unsigned short> Species_Names;  /*!< \brief Names of species in the mixture. */

    double Rgas; /*!< \brief Gas constant of the mixture. */

    double Le; /*!< \brief Lewis number. */

    RealVec mMasses; /*!< \brief Molar mass for each species. */

    RealVec Diff_Volumes; /*!< \brief Molecular diffusion volume for each species. */

    RealVec Ri; /*!< \brief Specific gas constant for each species. */

    RealVec Ys;    /*!<  \brief Mass fraction for each species. */

    RealVec Xs;    /*!<  \brief Molar fraction for each species. */

    RealVec Cs;    /*!<  \brief Actual concentration for each species. */

    RealVec Viscosities; /*!< \brief Viscosity for each species. */

    RealVec Thermal_Conductivities; /*!< \brief Thermal conductivity for each species. */

    std::vector<MyTuple> Mu_Spline; /*!< \brief Spline interpolation coefficient for viscosity computation. */

    std::vector<MyTuple> Kappa_Spline; /*!< \brief Spline interpolation coefficient for conductivity computation. */

    std::vector<MyTuple> Entr_Spline; /*!< \brief Spline interpolation coefficient for entropy computation. */

    std::vector<MyTuple> Cp_Spline; /*!< \brief Spline interpolation coefficient for specific heat at constant pressure computation. */

    std::vector<MyTuple> Enth_Spline; /*!< \brief Spline interpolation coefficient for enthalpy computation. */

    RealVec Enthalpies; /*!< \brief Enthalpy for each species. */

    RealVec Internal_Energies; /*!< \brief Internal Energy for each species. */

    RealVec CPs; /*!< \brief Specific heat at constant pressure for each species (Cp). */

    RealVec dPdYs; /*!< \brief Auxiliary vector for pressure derivatives w.r.t partial densities. */

    RealVec Formation_Enthalpies; /*!< \brief Formation enthalpy for each species. */

    RealMatrix Stoich_Coeffs_Reactants; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Products; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Reactants_Exp; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Products_Exp; /*!< \brief Stochiometric coefficents vector. */

    std::vector<bool> Reversible_Reactions; /*!< \brief Vector to check if a reaction is elementary. */

    RealVec As;  /*!< \brief Vector with exponential pre factor. */

    RealVec Betas;  /*!< \brief Vector with temperature exponent. */

    RealVec Temps_Activation;  /*!< \brief Vector with activation temperatures to estimate reaction rates for each reaction. */

    RealVec ys_over_mm; /*!< \brief Auxiliary vector to compute viscosity and thermal conductivity. */

    RealVec rhoUdiff; /*!< \brief Auxiliary vector for mass production term in case of constant Lewis number. */

    RealVec Dm_coeffs; /*!< \brief Auxiliary vector for effective diffusion coefficients. */

    RealVec Omega; /*!< \brief Auxiliary vector for mass production term. */

    RealMatrix Dij; /*!< \brief Auxiliary matrix for diffusion binary coefficients. */

    RealMatrix Gamma; /*!< \brief Auxiliary matrix for Stefan-Maxwell equations. */

    RealVec Forward_Rates; /*!< \brief Auxiliary vector for forward rate of each reaction. */

    RealVec Backward_Rates; /*!< \brief Auxiliary vector for backward rate of each reaction. */

    RealVec Kc; /*!< \brief Auxiliary vector for equilibrium constants. */

    RealVec Kc_Derivatives; /*!< \brief Auxiliary vector for equilibrium constants derivative. */

    RealMatrix Source_Jacobian; /*!< \brief Auxiliary matrix for source chemistry Jacobian. */

    std::array<double,5> Fuel_Data; /*!< \brief Auxiliary array for fuel data in case of regression boundary condition. */

  private:

    enum {T_DATA_SPLINE = 0, X_DATA_SPLINE = 1, Y_DATA_SPLINE = 2}; /*!< \brief Enumerator for spline indexes. */

    enum {A1_INDEX = 1, A2_INDEX = 2, EA1_INDEX = 3, EA2_INDEX = 4, TBAR_INDEX = 5}; /*!< \brief Enumerator for regression rate data indexes. */

    bool CGS_Units;             /*!< \brief Bool to check wheter we are reading chemical data in CGS or SI units. */

    std::map<unsigned short, bool> Available_Backward_Rate; /*!< \brief Auxiliary map to check wheter for a reversible reaction
                                                                       we already know the backward rate data. */

    std::map<unsigned short, double> As_back; /*!< \brief Auxiliary map for Arrhenius constant for backward rate. */

    std::map<unsigned short, double> Betas_back; /*!< \brief Auxiliary map for temperature exponent for backward rate. */

    std::map<unsigned short, double> Temps_Activation_back; /*!< \brief Auxiliary map for activation temperatures for backward rate. */

    std::vector<std::vector<unsigned short>> Product_Species_Negative_Exponent; /*!< \brief Auxiliary vector of vector to detect negative exponents
                                                                                            for backward rates. */

    std::vector<std::vector<unsigned short>> Reactant_Species_Negative_Exponent; /*!< \brief Auxiliary vector of vector to detect negative exponents
                                                                                            for backward rates. */
  }; /*-- End of class ReactingModelLibrary ---*/

} /*-- End of Namespace Framework ---*/

#endif
