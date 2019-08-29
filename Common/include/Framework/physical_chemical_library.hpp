#ifndef SU2_PHYSICAL_CHEMICAL_LIBRARY
#define SU2_PHYSICAL_CHEMICAL_LIBRARY

#include "not_copyable.hpp"

#include <vector>
#include "../../../externals/Eigen/Dense"

namespace Framework  {

  typedef std::vector<double> RealVec;
  typedef Eigen::MatrixXd RealMatrix;

  /*!
   * /brief Provides an abstract interface for libraries that compute the physical properties.
  */
  class PhysicalChemicalLibrary: public Common::NotCopyable<PhysicalChemicalLibrary> {

  public:
    typedef const std::string& Arg1;
    typedef const std::string& Arg2;

  public:
    /*!
     * \brief Class constructor(two arguments).
     */
    PhysicalChemicalLibrary(Arg1 config_name, Arg2 lib_path_name): Config_File(config_name), Lib_Path(lib_path_name),
                                                                   Lib_Setup(false), nSpecies(), nReactions() {}

    /*!
     * \brief Default destructor.
     */
    virtual ~PhysicalChemicalLibrary() {}

    /*!
     * \brief Get the name of the base class of library.
    */
    inline static std::string GetBaseName(void) {
      return "PhysicalChemicalLibrary";
    }

    /*!
     *\brief Setup the path library name.
     *\param[in] lib_path_name - Path of the library
    */
    inline void SetLibPathName(const std::string& lib_path_name) {
      Lib_Path = lib_path_name;
    }

    /*!
     * Get the number of species in the mixture
     */
    inline unsigned short GetnSpecies(void) const {
      return nSpecies;
    }

    /*!
     * Set the number of species in the mixture
    */
    inline void SetnSpecies(const unsigned short ns) {
      nSpecies = ns;
    }

    /*!
     * Get the number of reactions in the mixture
    */
    inline unsigned short GetnReactions(void) const {
      return nReactions;
    }

    /*!
     * Set the number of reactions in the mixture
    */
    inline void SetNReactions(const unsigned short nr) {
      nReactions = nr;
    }

    /*!
     * \brief Setups the data of the library.
     */
    virtual void Setup(void) = 0;

    /*!
     * \brief Unsetups the data of the library.
     */
    virtual void Unsetup(void) = 0;

    /*!
     * \brief Get the index of a species.
     * \param[in] name_species - Name of the desired species
     */
    virtual unsigned short GetIndexSpecies(const std::string& name_species) const = 0;

    /*!
     * Set the gas constant for each species [J/(Kg*K)]
    */
    virtual void SetRiGas(void) = 0;

    /*!
     * Get the gas constant for each species [J/(Kg*K)]
    */
    virtual RealVec GetRiGas(void) const = 0;

    /*!
     * Get the gas constant for a desired species [J/(Kg*K)]
     * \param[in] iSpecies - index of the desired species;
    */
    virtual double GetRiGas(unsigned short iSpecies) const = 0;

    /*!
     * \brief Get the gas constant for the mixture [J/(Kg*K)]
    */
    virtual double GetRgas(void) const = 0;

    /*!
     * \brief Set the gas constant for the mixture [J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    virtual void SetRgas(const RealVec& ys) = 0;

    /*!
     * \brief Compute the gas constant for the mixture [J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    virtual double ComputeRgas(const RealVec& ys) = 0;

    /*!
     * \brief Get the molar masses of the species
    */
    virtual RealVec GetMolarMasses(void) const = 0;

    /*!
     * Set the molar fractions of elements Xn.This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] xs - The vector of the mass fractions of elements
    */
    virtual void SetMolarFractions(const RealVec& xs) = 0;

    /*!
     * Set the mass fractions. This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual void SetMassFractions(const RealVec& ys) = 0;

    /*!
     * \brief Set the molar fractions from mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    virtual void SetMolarFromMass(const RealVec& ys) = 0;

    /*!
     * \brief Get the molar fractions from mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    virtual RealVec GetMolarFromMass(const RealVec& ys) = 0;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
    */
    virtual void SetMassFromMolar(const RealVec& xs) = 0;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
    */
    virtual RealVec GetMassFromMolar(const RealVec& xs) = 0;

    /*!
     * \brief Compute the specific heat ratio and the speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[out] gamma - specific heat ratio
     * \param[out] sound_speed - speed of sound
    */
    virtual void Gamma_FrozenSoundSpeed(const double temp, const RealVec& ys, double& gamma, double& sound_speed) = 0;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     */
    virtual double ComputeFrozenGamma(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] gamma - specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual double ComputeFrozenSoundSpeed_FromGamma(const double temp, const double gamma, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] sound_speed - frozen speed of sound
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual double ComputeFrozenGamma_FromSoundSpeed(const double temp, const double sound_speed, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    virtual double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys, const double press, const double rho) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] gamma - specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    virtual double ComputeFrozenSoundSpeed_FromGamma(const double gamma, const RealVec& ys, const double press, const double rho) = 0;

    /*!
     * \brief Compute the density, the enthalpy and the internal energy
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     * \param[out] dhe - vector with density, enthalpy, energy (output) for thermal equilibrium
    */
    virtual void Density_Enthalpy_Energy(const double temp, const double pressure, const RealVec& ys, RealVec& dhe) = 0;

    /*!
     * \brief Compute the pressure given temperature and density.
     * \param[in] temp - temperature
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    virtual double ComputePressure(const double temp, const double rho, const RealVec& ys) = 0;

    /*!
     * \brief Compute the density given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
    */
    virtual double ComputeDensity(const double temp, const double pressure, const RealVec& ys) = 0;

    /*!
     * \brief Compute the temperature given density and pressure.
     * \param[in] pressure - pressure
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    virtual double ComputeTemperature(const double pressure, const double rho, const RealVec& ys) = 0;

    /*!
     * \brief Compute the density given sound speed and gamma.
     * \param[in] sound_speed2 - squared sound speed
     * \param[in] gamma - specific heats ratio
     * \param[in] ys - mass fractions
    */
    virtual double ComputeTemperature_FromGamma(const double sound_speed2, const double gamma, const RealVec& ys) = 0;

    /*!
     * \brief Compute the internal energy per unit of mass at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     */
    virtual double ComputeEnergy(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Return the formation enthalpies per unit mass of species
    */
    virtual RealVec GetFormationEnthalpies(void) const = 0;

    /*!
     * \brief Return the static enthalpy per unit of mass
     * \param[in] temp - the mixture temperature
     * \param[in] ys - mass fractions
     * \return Mixture static enthalpy (output)
    */
    virtual double ComputeEnthalpy(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Set the static enthalpy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    virtual void SetPartialEnthalpy(const double temp) = 0;

    /*!
     * \brief Set the static enthalpy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    virtual void SetPartialEnthalpy(const double temp, unsigned short iSpecies) = 0;

    /*!
     * \brief Return the static enthalpy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    virtual RealVec ComputePartialEnthalpy(const double temp) = 0;

    /*!
     * \brief Compute the static enthalpy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    virtual double ComputePartialEnthalpy(const double temp, unsigned short iSpecies) = 0;

    /*!
     * \brief Set the internal energy per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    virtual void SetPartialEnergy(const double temp) = 0;

    /*!
     * \brief Set the internal energy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    virtual void SetPartialEnergy(const double temp, unsigned short iSpecies) = 0;

    /*!
     * \brief Return the internal energy of each species.
     * \param[in] temp - the mixture temperature
    */
    virtual RealVec ComputePartialEnergy(const double temp) = 0;

    /*!
     * \brief Compute the internal energy per unit of mass for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] iSpecies - index of desired species
    */
    virtual double ComputePartialEnergy(const double temp, unsigned short iSpecies) = 0;

    /*!
     * \brief Computes the pressure derivative w.r.t. partial densities.
     * \param[in] temp - temperature
     * \param[in] gamma - frozen specific heat ratio
    */
    virtual RealVec ComputedP_dYs(const double temp, const double gamma) = 0;

    /*!
     * \brief Compute the pressure derivative w.r.t. partial density for a desired species.
     * \param[in] temp - the mixture temperature
     * \param[in] gamma - frozen specific heat ratio
     * \param[in] iSpecies - index of desired species
    */
    virtual double ComputedP_dYs(const double temp, const double gamma, unsigned short iSpecies) = 0;

    /*!
     * \brief Set the actual concetration for each species
     * \param[in] rho - the mixture density
     * \param[in] ys - the actual mass fractions
    */
    virtual void SetConcentration(const double rho, const RealVec& ys) = 0;

    /*!
     * \brief Compute the mixture total concentration
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    virtual double ComputeConcentration(const double rho, const RealVec& ys) = 0;

    /*!
     * \brief Return the specific heat at constant pressure per unit of mass of each species.
     * \param[in] temp - the mixture temperature
    */
    virtual RealVec ComputeCps(const double temp) = 0;

    /*!
     * \brief Compute the specific heat at constant pressure
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Cp - specific heat at constant pressure
     */
    virtual double ComputeCP(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Compute the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Cv - specific heat at constant volume
    */
    virtual double ComputeCV(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Compute the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] cp - specific heat at constant pressure
     * \return Cv - specific heat at constant volume
    */
    virtual double ComputeCV_FromCP(const double cp, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] cp - specific heat at constant pressure
     * \param[in] ys - The vector of the mass fractions of species
     */
    virtual double ComputeFrozenGamma_FromCP(const double cp, const RealVec& ys) = 0;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - mxiture temperature;
     * \param[in] cp - specific heat at constant pressure
     * \param[in] ys - The vector of the mass fractions of species
     */
    virtual double ComputeFrozenSoundSpeed_FromCP(const double temp, const double cp, const RealVec& ys) = 0;

    /*!
     * \brief Computes the specific heat at constant pressure.
     * \param[in] temp - temperature
     * \param[in] sound_speed - frozen speed of sound
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual double ComputeCP_FromSoundSpeed(const double temp, const double sound_speed, const RealVec& ys) = 0;

    /*!
     * \brief Computes the specific heat at constant pressure.
     * \param[in] gamma - frozen specific heat ratio
     * \param[in] ys - The vector of the mass fractions of species
    */
    virtual double ComputeCP_FromGamma(const double gamma, const RealVec& ys) = 0;

    /*!
     * \brief Compute the thermal conductivity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Total thermal conductivity for thermal non equilibrium
     */
    virtual double ComputeLambda(const double temp, const RealVec& ys) = 0;

    /*!
     * \brief Compute the dynamic viscosity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Total molecular viscosity for thermal non equilibrium
    */
    virtual double ComputeEta(const double temp, const RealVec& ys) = 0;

    /*!
     * Return the diffusion velocities of species multiplied by the species
     * densities for nonequilibrium computations
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
     * \return Diffusion coefficient with constant Lewis number for each species
    */
    virtual RealVec GetRhoUdiff(const double temp, const double rho, const RealVec& ys) = 0;

    /*!
     * Return the mass production/destruction terms [kg m^-3 s^-1] in chemical
     * non-equilibrium based on Arrhenius's formula.
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
     * \return Mass production terms
    */
    virtual RealVec GetMassProductionTerm(const double temp, const double rho, const RealVec& ys) = 0;

    /*!
     * Compute the Jacobian of source chemistry. NOTE: It requires SetReactionRates call
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \return Source Jacobian
     */
     virtual RealMatrix GetSourceJacobian(const double temp, const double rho) = 0;

    /*!
     * \brief Return the effective diffusion coefficients to solve Stefan-Maxwell equation using Sutton algorithm
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[in] ys - mass fractions in the mixture
     */
    virtual RealVec GetDiffCoeffs(const double temp, const double pressure, const RealVec& ys) = 0;

    /*!
     * \brief Return the binary diffusion coefficients
     * densities for nonequilibrium computations
     * \param[in] pressure - the mixture pressure
     * \param[in] temp - the mixture temperature
    */
    virtual RealMatrix GetDij_SM(const double pressure, const double temp) = 0;

    /*!
     * \brief Return the matrix of Stefan-Maxwell equations
     * \param[in] rho - the mixture density
     * \param[in] xs - current molar fractions
     * \param[in] ys - current mass fractions
     * \param[in] val_Dij - current binary diffusion coefficients
     */
    virtual RealMatrix GetGamma(const double rho, const RealVec& xs, const RealVec& ys, const RealMatrix& val_Dij) = 0;

    /*!
     * \brief Compute the regression rate at specified temperature with an empirical law
     * \param[in] temp - the mixture temperature
     */
    virtual double ComputeRegressionRate(const double temp) = 0;

    /*!
     * \brief Read all physical data about fuel
     * \param[in] f_name - File with the fuel data to be read
     */
    virtual void ReadDataFuel(const std::string& f_nmae) = 0;

  protected:

    std::string Config_File; /*!\brief File with all the configuration options */

    std::string Lib_Path; /*!\brief Path to some input library data */

    bool Lib_Setup; /*!\brief Bool to check if library is setup */

    unsigned short nSpecies; /*!< \brief Number of species. */

    unsigned short nReactions; /*!< \brief Number of reactions. */

  public:

    static constexpr double NA = 6.02214129*1.0e23; /*!< \brief Avogadro number. */

    static constexpr double KB = 1.3806488*1.0e-23; /*!< \brief Boltzmann constant. */

    static constexpr double R_ungas = NA*KB*1.0e3; /*!< \brief Universal constant of perfect gas. (J/kmol*K) */

    static constexpr double R_ungas_scal = 1.9858775; /*!< \brief Universal gas constant in cal/(mol K). */

    static constexpr double R_ungas_atm = 1.0e-3*0.082057338; /*!< \brief Universal gas constant in m3*atm/(mol K). */

  }; /*-- End of class PhysicalChemicalLibrary ---*/

} /*-- End of namespace Framework ---*/

#endif
