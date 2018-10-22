#ifndef SU2_PHYSICAL_PROPERTY_LIBRARY
#define SU2_PHYSICAL_PROPERTY_LIBRARY

#include "concrete_provider.hpp"
#include "builder_provider.hpp"
#include "datatype_structure.hpp"

namespace Framework  {

  /*!
    * /brief Provides an abstract interface for libraries that compute the physical properties.
    */
    class PhysicalPropertyLibrary {

    public:
      typedef Common::ConcreteProvider<PhysicalPropertyLibrary> Provider;
      typedef const std::string& Arg1;

      typedef std::vector<su2double> RealVec;
      typedef std::vector<RealVec> RealMatrix;

    public:

      /*!
       * \brief Default constructor.
       */
       PhysicalPropertyLibrary():Lib_Setup(false),nSpecies(0),nReactions(0) {}

      /*!
       * \brief Default destructor.
       */
       virtual ~PhysicalPropertyLibrary() = default;

       /*!
        *\brief Setups the path library name.
       */
       inline void SetLibPathName(const std::string& lib_path_name)  {
         Lib_Path = lib_path_name;
       }

       /*!
       * \brief Get the name of the library.
       */
       inline static std::string GetBaseName(void) {
          return "PhysicalPropertyLibrary";
       }

       inline bool IsSetup(void) const {
         return Lib_Setup;
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
         * Get the number of species in the mixture
         */
       inline const unsigned short GetNSpecies(void) const {
         return nSpecies;
       }

       /*!
         * Set the number of species in the mixture
         */
       inline void SetNSpecies(const unsigned short ns) {
         nSpecies = ns;
       }

       /*!
         * Get the number of reactions in the mixture
         */
       inline const unsigned GetNReactions(void) const {
         return nReactions;
       }

       /*!
         * Set the number of reactions in the mixture
         */
       inline void GetNReactions(const unsigned nr) {
         nReactions = nr;
       }

       /*!
         * Set the constant of gases for each species [J/(Kg*K)]
         */
       virtual void SetRiGas(RealVec& Ri) = 0;

       /*!
         * Get the constant of perfect gases [J/(Kg*K)]
         */
       virtual su2double GetRgas(void) const = 0;

       /*!
         * Set the constant of perfect gases [J/(Kg*K)]
         */
       virtual void SetRgas(const RealVec& Ys,const RealVec& Ri) = 0;

       /*!
         * \brief Get the molar masses of the species
         * \pre mm.size() = number of species
         */
       virtual void GetMolarMasses(RealVec& mm) = 0;

       /*!
         * \brief Set the IDs of the molecules in the mixture
         */
       virtual void SetMoleculesIDs(std::vector<unsigned>& Ids) = 0;

       /*!
         * \brief Calculates the thermal conductivity given temperature and pressure
         * \param[in] temp - temperature
         * \param[in] pressure - pressure
         * \return Total thermal conductivity for thermal equilibrium
         */
       virtual su2double Lambda(su2double& temp, su2double& pressure) = 0;

       /*!
         * \brief Calculates the dynamic viscosity given temperature and pressure
         * \param[in] temp - temperature
         * \param[in] pressure - pressure
         */
       virtual su2double Eta(su2double& temp,su2double& pressure) = 0;

       /*!
         * \brief Calculates the specific heat ratio
         */
       virtual su2double Gamma(void) = 0;

       /*!
         * \brief Calculates the specific heat ratio and the speed of sound in thermal equilibrium.
         * \param[in] temp - temperature
         * \param[in] pressure - pressure
         * \param[in] rho - density
         * \return gamma - specific heat ratio (output)
         * \return sound_speed - speed of sound (output)
       */
       virtual void Gamma_FrozenSoundSpeed(su2double& temp,su2double& pressure,su2double& rho,su2double&gamma,su2double& sound_speed) = 0;

       /*!
         * \brief Calculates the density, the enthalpy and the internal energy
         * \param[in] temp - temperature
         * \param[in] pressure - pressure
         * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
       */
       virtual void Density_Enthalpy_Energy(su2double& temp,su2double& pressure,RealVec& dhe) = 0;

       /*!
         * \brief Calculates the density given temperature and pressure.
         * \param[in] temp - temperature
         * \param[in] pressure - pressure
       */
       virtual su2double Density(const su2double& temp,const su2double& pressure) = 0;

       /*!
         * \brief Calculates the internal energy at given temperature and pressure.
         * \param[in] temp temperature
         * \param[in] pressure pressure
       */
       virtual su2double Energy(su2double& temp,su2double& pressure) = 0;

       /*!
         * Calculates the enthalpy in LTE conditions
         * at given temperature and pressure.
         * \param[in] temp      temperature
         * \param[in] pressure  pressure
       */
       virtual su2double Enthalpy(su2double& temp,su2double& pressure) = 0;

       /*!
         * \brief Gets the molar fractions.
         * \param[in] ys The vector of the mass fractions of species (input)
         * \param[out] xs The vector of the molar fractions of species (output)
       */
       virtual void GetMolarFractions(const RealVec& ys, RealVec& xs) = 0;

       /*!
        * \brief Sets the molar fractions of elements Xn.This function should be called before getting
        * \brief thermodynamic quantities or transport properties.
        * \param[in] xn - The vector of the mass fractions of elements
       */
       virtual void SetMolarFractions(const RealVec& xs) = 0;

       /*!
         * \brief Sets the molar fractions of elements Xn.This function should be called before getting
         * \brief thermodynamic quantities or transport properties.
       */
       virtual void SetMolarFractions(void) = 0;

       /*!
        * \brief Gets the mass fractions.
        * \param[in] xs The vector of the molar fractions of species (input)
        * \param[out] ys The vector of the mass fractions of species (output)
        */
       virtual void GetMassFractions(const RealVec& xs, RealVec& ys) = 0;

       /*!
        * \brief Sets the mass fractions. This function should be called before getting
        * \brief thermodynamic quantities or transport properties.
        * \param[in] ys The vector of the mass fractions of species
       */
       virtual void SetMassFractions(const RealVec& ys) = 0;

       /*!
        * \brief Sets the mass fractions. This function should be called before getting
        * \brief thermodynamic quantities or transport properties.
        * \param[in] ys The vector of the mass fractions of species
       */
       virtual void SetMassFractions(void) = 0;

       /*!
        * \brief Sets the mole fractions of elements Xn starting from the given species mass fractions Yn
        * \param[in] ys - The vector of the mass fractions of species
       */
       virtual void SetMolarFromMass(const RealVec& ys) = 0;

       /*!
        * \brief Returns the formation enthalpies per unit mass of species
        * \param[out] hs - species formation enthalpies (output)
       */
       virtual void GetFormationEnthalpies(RealVec& hs) = 0;

       /*!
        * \brief Returns the total enthalpies per unit mass of species
        * \param[in] temp - the mixture temperature
        * \param[in] pressure - the mixture pressure
        * \param[out] hsTot - species total enthalpy (output)
       */
       virtual void GetSpeciesTotEnthalpies(su2double& temp,su2double& pressure,RealVec& hsTot) = 0;

       /*!
        * \brief Returns the mass production/destruction terms [kg m^-3 s^-1] in chemical
        *  brief non-equilibrium based on Arrhenius's formula.
        * \param[in] pressure the mixture pressure
        * \param[in] temp the mixture temperature
        * \param[in] ys the species mass fractions
        * \param[in] omega the mass production terms
        * \param[in] jacobian the Jacobian matrix of the mass production terms
       */
       virtual void GetMassProductionTerm(su2double& temp,su2double& pressure,su2double& rho,
                                          RealVec& ys,
                                          RealVec& omega,
                                          std::vector<RealVec>& jacobian) = 0;

       /*!
        * Returns the source terms species continuity and energy equations
        * \param[in] temp the mixture temperature
        * \param[in] pressure the mixture pressure
        * \param[in] ys the species mass fractions
        * \param[in] rho the mixture density
        * \param[in] omega the mass producrtion term
        * \param[in] omegav the source term
       */
       virtual void GetSource(su2double& temp,su2double& pressure,su2double& rho,
                              RealVec& ys,
                              RealVec& omega,RealVec& omegav) = 0;

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
       virtual void GetRhoUdiff(su2double& temp,su2double& pressure, RealVec& normConcGradients,
                                RealVec& rhoUdiff) = 0;

      /*!
       * \brief Returns the laminar Prandtl number
       */
      virtual su2double GetPrandtl_Lam(void) = 0;

    protected:

       bool Lib_Setup; /*!\brief Path to some input library data */

       std::string Lib_Path; /*!\brief Path to some input library data */

       unsigned short nSpecies; /*!< \brief Number of species. */

       unsigned nReactions; /*!< \brief Number of reactions. */

     public:

       static constexpr su2double NA = 6.02214129*1e23; /*!< \brief Avogadro number. */

       static constexpr su2double KB = 1.3806488*1e-23; /*!< \brief Boltzmann constant. */

       static constexpr su2double R_ungas = 8.31446215; /*!< \brief Universal constant of perfect gas. */

  }; /*-- End of class PhysicalPropertyLibrary ---*/

} /*-- End of namespace Framework ---*/
#endif
