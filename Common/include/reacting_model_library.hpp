#ifndef SU2_REACTING_MODEL_LIBRARY
#define SU2_REACTING_MODEL_LIBRARY

#include "physical_property_library.hpp"

#include <vector>
#include <cassert>

#ifndef NDEBUG
#define NDEBUG

namespace Framework {

  /*!
   * \brief Provides a particular library definition to compute the physical and chemical properties.
   */

  class ReactingModelLibrary: public  Framework::PhysicalPropertyLibrary {

  public:

    typedef Framework::PhysicalPropertyLibrary::RealVec RealVec;
    typedef Framework::PhysicalPropertyLibrary::RealMatrix RealMatrix;


  public: // functions

    /*!
      * \brief Constructor with the name of the library.
      */
    explicit ReactingModelLibrary(const std::string& name):PhysicalPropertyLibrary(),Lib_Name(name) {};


    /*!
      * \brief Default destructor.
      */
    ~ReactingModelLibrary() = default;

    /*!
      * \brief Setups the library name.
      */
    inline void SetLibName(const std::string& library_name) {
      Lib_Name = library_name;
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
      * Set the constant of gases for each species [J/(Kg*K)]
      */
    inline void SetRiGas(RealVec& Ri) override {
      assert(Ri.size()==nSpecies);
      for(auto i=0;i<nSpecies;++i)
        Ri[i]=R_ungas/mMasses[i];
    }

    /*!
      * Get the constant of perfect gases [J/(Kg*K)]
      */
    inline su2double GetRgas(void) const override {
      return Rgas;
    }

    /*!
      * Set the constant of perfect gases [J/(Kg*K)]
      */
    void SetRgas(const RealVec& Ys,const RealVec& Ri) override;

    /*!
      * \brief Get the molar masses of the species
      * \pre mm.size() = number of species
      */
    inline void GetMolarMasses(RealVec& mm) override {
      assert(mm.size()==nSpecies);
      mm=mMasses;
    }

    /*!
      * \brief Set the IDs of the molecules in the mixture
      */
    void SetMoleculesIDs(std::vector<unsigned>& Ids) override;

    /*!
      * \brief Calculates the thermal conductivity given temperature and pressure
      * \param[in] temp - temperature
      * \param[in] pressure - pressure
      * \return Total thermal conductivity for thermal equilibrium
      */
    su2double Lambda(su2double& temp, su2double& pressure) override;

    /*!
      * \brief Calculates the dynamic viscosity given temperature and pressure
      * \param[in] temp - temperature
      * \param[in] pressure - pressure
      */
    su2double Eta(su2double& temp,su2double& pressure) override;

    /*!
      * \brief Calculates the specific heat ratio
      */
    su2double Gamma(void) override;

    /*!
      * \brief Calculates the specific heat ratio and the speed of sound in thermal equilibrium.
      * \param[in] temp - temperature
      * \param[in] pressure - pressure
      * \param[in] rho - density
      * \return gamma - specific heat ratio (output)
      * \return sound_speed - speed of sound (output)
    */
    void Gamma_SoundSpeed(su2double& temp,su2double& pressure,su2double& rho,su2double&gamam,su2double& sound_speed) override;

    /*!
      * \brief Calculates the density, the enthalpy and the internal energy
      * \param[in] temp - temperature
      * \param[in] pressure - pressure
      * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
    */
    void Density_Enthalpy_Energy(su2double& temp,su2double& pressure,RealVec& dhe) override;

    /*!
      * \brief Calculates the density given temperature and pressure.
      * \param[in] temp - temperature
      * \param[in] pressure - pressure
    */
    su2double Density(const su2double& temp,const su2double& pressure) override;

    /*!
      * \brief Calculates the internal energy at given temperature and pressure.
      * \param[in] temp temperature
      * \param[in] pressure pressure
    */
    su2double Energy(su2double& temp,su2double& pressure) override;

    /*!
      * Calculates the enthalpy in LTE conditions
      * at given temperature and pressure.
      * \param[in] temp      temperature
      * \param[in] pressure  pressure
    */
    su2double Enthalpy(su2double& temp,su2double& pressure) override;

    /*!
      * \brief Gets the molar fractions.
      * \param[in] ys The vector of the mass fractions of species (input)
      * \param[out] xs The vector of the molar fractions of species (output)
    */
    void GetMolarFractions(const RealVec& ys, RealVec& xs) override;

    /*!
     * \brief Sets the molar fractions of elements Xn.This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] xn - The vector of the mass fractions of elements
    */
    void SetMolarFractions(const RealVec& xs) override;

    /*!
      * \brief Sets the molar fractions of elements Xn.This function should be called before getting
      * \brief thermodynamic quantities or transport properties.
    */
    void SetMolarFractions(void) override;

    /*!
     * \brief Gets the mass fractions.
     * \param[in] xs The vector of the molar fractions of species (input)
     * \param[out] ys The vector of the mass fractions of species (output)
     */
    void GetMassFractions(const RealVec& xs, RealVec& ys) override;

    /*!
     * \brief Sets the mass fractions. This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] ys The vector of the mass fractions of species
    */
    void SetMassFractions(const RealVec& ys) override;

    /*!
     * \brief Sets the mass fractions. This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] ys The vector of the mass fractions of species
    */
    void SetMassFractions(void) override;

    /*!
     * \brief Sets the mole fractions of elements Xn starting from the given species mass fractions Yn
     * \param[in] ys - The vector of the mass fractions of species
    */
    void SetMolarFromMass(const RealVec& ys) override;

    /*!
     * \brief Returns the formation enthalpies per unit mass of species
     * \param[out] hs - species formation enthalpies (output)
    */
    inline void GetFormationEnthalpies(RealVec& hsTot) override;

    /*!
     * \brief Returns the total enthalpies per unit mass of species
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[out] hsTot - species total enthalpy (output)
    */
    void GetSpeciesTotEnthalpies(su2double& temp,su2double& pressure,RealVec& hsTot) override;

    /*!
     * \brief Returns the mass production/destruction terms [kg m^-3 s^-1] in CHEMICAL
     *  brief NONEQUILIBRIUM based on Arrhenius's formula.
     * \param[in] pressure the mixture pressure
     * \param[in] temp the mixture temperature
     * \param[in] ys the species mass fractions
     * \param[in] omega the mass production terms
     * \param[in] jacobian the Jacobian matrix of the mass production terms
    */
    void GetMassProductionTerm(su2double& temp,su2double& pressure,su2double& rho,
                               RealVec& ys,
                               RealVec& omega,
                               RealMatrix& jacobian) override;

    /*!
     * Returns the source terms species continuity and energy equations
     * \param[in] temp the mixture temperature
     * \param[in] pressure the mixture pressure
     * \param[in] ys the species mass fractions
     * \param[in] rho the mixture density
     * \param[in] omega the mass producrtion term
     * \param[in] omegav the source term
    */
    void GetSource(su2double& temp,su2double& pressure,su2double& rho,
                   RealVec& ys,
                   RealVec& omega,RealVec& omegav) override;

   /*!
    * Returns the diffusion velocities of species multiplied by the species
    * densities for nonequilibrium computations
    * \param[in] temp the mixture temperature
    * \param[in] pressure the mixture pressure
    * \param[in] normConcGradients the cell normal gradients of species mass fractions
    * \param[in] rhoUdiff   (1) Constant Lewis number
    *                   (2) Stefan-Maxwell model: rho*species diffusive velocity
    *                   (3) Fick, Fick+Ramshaw  : rowwise-ordered matrix of coefficients
    */
    void GetRhoUdiff(su2double& temp,su2double& pressure, RealVec& normConcGradients,
                     RealVec& rhoUdiff) override;

   /*!
    * \brief Returns the reference temperature
    */
    inline su2double GetTemperature_Ref(void) override {
      return T_ref;
    }

    /*!
     * \brief Returns the reference viscosity
     */
    inline su2double GetViscosity_Ref(void) override {
      return Mu_ref;
    }

    /*!
     * \brief Returns the free stream viscosity
     */
    inline su2double GetViscosity_FreeStream(void) override {
      return Viscosity_Mixture;
    }

    /*!
     * \brief Returns the laminar Prandtl number
     */
    inline su2double GetPrandtl_Lam(void) override {
      return Lam_Pr;
    }

    /*!
     * \brief Returns the free stream density
     */
    inline su2double GetDensity_FreeStream(void) override {
      return Rho_inf;
    }

    /*!
     * \brief Returns the free stream pressure
     */
    inline su2double GetPressure_FreeStream(void) override {
      return P_inf;
    }

    /*!
     * \brief Returns the free stream temperature
     */
    inline su2double GetTemperature_FreeStream(void) override {
      return T_inf;
    }

    /*!
     * \brief Returns the free stream temperature
     */
    inline su2double GetMach(void) override {
      return Mach;
    }

  private:

    /*!
    * \brief Returns the reaction contribution to source term of species continuity equations.
    * \param[in] temp - the mixture temperature
    * \param[in] pressure - the mixture pressure
    * \param[in] ys - the species mass fractions
    * \param[in] mMasses -  molar masses
    * \param[in] rho the mixture density
    * \param[in] omega the mass production term
    */
    void OmegaContribution(su2double& temp,su2double& pressure,su2double& rho,
                         const RealVec& ys,const RealVec& mMasses,RealVec& omega);

   /*!
    * \brief Returns the forward reaction rate coefficient.
    * \param[in] temp - the corrected temperature Tq
    */
    su2double Kf(su2double& temp);

    /*!
     * \brief Returns the backward reaction rate coefficient.
     * \param[in] temp the mixture temperature
    */
    su2double Kb(su2double& temp);

    /*!
      * \brief Read mixture data
      */
    void ReadDataMixture(const std::string& f_name);

    /*!
    * \brief Read chemistry data
    * \param[in] f_name - name of the file with the nvolved reactions
    */
    void ReadDataChem(const std::string& f_name);

    /*!
      * \brief Read transport data
      * \param[in] f_name - name of the file with the properties
      * \param[in] iSpecies - index of the desired species
    */
    void ReadDataTransp(const std::string& f_name, const unsigned short iSpecies);

    /*!
      * \brief Read thermodynamical data
      * \param[in] f_name - name of the file with the properties
      * \param[in] iSpecies - index of the desired species
    */
    void ReadDataThermo(const std::string& f_name, const unsigned short iSpecies);

    /*!
      * \brief Read species involved in a chemical reaction from line
      * \param[in] idx index of position in "line" separating reaction definition and coefficients
      * \param[in] StochVec Vector containing Stoichiometric coefficients
      * \param[in] Line line to be read
    */
    void ReadReactSpecies(RealVec& stoch_vec, std::string& line);

    /*!
      * Read coefficients to compute reaction rates from line
      * \param[in] idx index of position in "line" separating reaction definition and coefficients
      * \param[in] ChemCoefs Coefficients to compute reaction rates
      * \param[in] Line line to be read
    */
    void ReadChemCoefs(RealVec& chem_coeffs,std::string& Line);

  protected:

    std::string Lib_Name;  /*!< \brief Name of the library. */

    std::string File_Names; /*!< \brief Name of the file for reading thermodynamical and transport properties. */
                           // where to add it?

    std::vector<std::string> Species_Names;  /*!< \brief Names of species in the mixture. */

    std::vector<unsigned short> Atoms; /*!< \brief Constitutive atoms of species. */

    su2double Rgas; /*!< \brief Gas cosntant of the mixture. */

    su2double Le; /*!< \brief Lewis number. */

    su2double T_ref; /*!< \brief Reference tempeature. */

    su2double Mu_ref; /*!< \brief Reference viscosity. */

    su2double Lam_Pr; /*!< \brief Laminar Prandtl number. */

    su2double Mach; /*!< \brief Mach number. */

    su2double T_inf;  /*!< \brief Free stream tempeature. */

    su2double P_inf;  /*!< \brief Free stream pressure. */

    su2double Rho_inf;  /*!< \brief Free stream density. */

    su2double Viscosity_Mixture;  /*!< \brief Viscosity of the mixture. */

    su2double Thermal_Conductivity_Mixture; /*!< \brief Laminar Prandtl number. */

    su2double Cp_Mixture;

    su2double AF; /*!< \brief Acentric factor. */

    RealVec mMasses; /*!< \brief Molar mass for each species. */

    RealVec Ri; /*!< \brief Specific gas constant for each species. */

    RealVec Ys;    /*!<  \brief Mass fraction for each species. */

    RealVec Xs;    /*!<  \brief Molar fractionfor each species. */

    RealVec Viscosities; /*!< \brief Viscosity for each species. */

    RealVec Internal_Energies; /*!< \brief Internal energy for each species. */

    RealVec Enthalpies; /*!< \brief Enthalpy for each species. */

    RealVec Heat_Capacities; /*!< \brief Heat Capacity for each species. */

    RealVec CPs; /*!< \brief Specific heat at constant pressure for each species (Cp). */

    RealVec CVs; /*!< \brief Specific heat at constant volume for each species (Cv). */

    RealVec Thermal_Conductivities; /*!< \brief Thermal conductivity for each species. */

    RealVec Sound_Speeds; /*!< \brief Sound speed for each species. */

    RealVec Formation_Enthalpies; /*!< \brief Formation enthalpy for each species. */

    RealVec Stoich_Coeffs; /*!< \brief Stochiometric coefficents vector. */

    RealVec Chem_Coeffs;     /*!< \brief Vector with coefficients to estimate reaction rates. */

  }; /*-- End of class ReactingModelLibrary ---*/

} /*-- End of Namespace Frameworkn ---*/

#endif

#endif
