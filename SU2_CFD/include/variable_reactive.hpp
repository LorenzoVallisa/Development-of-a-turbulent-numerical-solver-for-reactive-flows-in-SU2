#ifndef SU2_VARIABLE_REACTIVE
#define SU2_VARIABLE_REACTIVE

#include "variable_structure.hpp"
#include "../../Common/include/physical_property_library.hpp"
#include "../../Common/include/su2_assert.hpp"

#include <Eigen/Dense>
#include <numeric>

/*! \class CReactiveEulerVariable
 *  \brief Main class for defining a variable for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerVariable:public CVariable {
public:
  typedef std::vector<su2double> RealVec;
  typedef su2double** SU2Matrix;
  typedef std::shared_ptr<Framework::PhysicalPropertyLibrary> LibraryPtr;

protected:
  static LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  static unsigned short nSpecies; /*!< \brief Number of species in the mixture. */
  unsigned short nPrimVarLim; /*!< \brief Number of primitive variables to limit in the problem. */

  /*--- Primitive variable definition ---*/
  RealVec    Primitive; /*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,a,Y1,...YNs) in compressible flows. */
  SU2Matrix  Gradient_Primitive; /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    dPdU;                 /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  RealVec    dTdU;                /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */

public:

  /**
   * Mapping between the primitive variable name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_PRIM = 0;
  static constexpr unsigned short VX_INDEX_PRIM = 1;
  static const unsigned short P_INDEX_PRIM;
  static const unsigned short RHO_INDEX_PRIM;
  static const unsigned short H_INDEX_PRIM;
  static const unsigned short A_INDEX_PRIM;
  static const unsigned short RHOS_INDEX_PRIM;

  /**
   * Mapping between the solution variable name and its position in the physical data
   */
  static constexpr unsigned short RHO_INDEX_SOL = 0;
  static constexpr unsigned short RHOVX_INDEX_SOL = 1;
  static const unsigned short RHOE_INDEX_SOL;
  static const unsigned short RHOS_INDEX_SOL;

  /**
   * Mapping between the primitive variable gradient name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_GRAD = 0;
  static constexpr unsigned short VX_INDEX_GRAD = 1;
  static const unsigned short P_INDEX_GRAD;

  /**
   * Mapping between the primitivelimited variable name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_LIM = 0;
  static constexpr unsigned short VX_INDEX_LIM = 1;
  static const unsigned short P_INDEX_LIM;

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerVariable();

  /*!
   * \overload Class constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config);

  /*!
	 * \overload Class constructor
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the flow energy (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity, const su2double val_temperature,
                         unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config);

	/*!
	 * \overload Class constructor
	 * \param[in] val_solution - Vector with the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                         unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerVariable();

  /*!
   * \brief Get the number of dimensions of the problem.
   * \return Number of dimensions of the problem.
   */
  inline static const unsigned short GetnDim(void) {
    return nDim;
  }

  /*!
   * \brief Get the loaded library
   * \return Pointer to loaded library
   */
  inline static LibraryPtr GetLibrary(void) {
    return library;
  }

  /*!
   * \brief Get the number of species in the mixture.
   * \return Number of species in the problem.
   */
  inline static const unsigned short GetnSpecies(void) {
    return nSpecies;
  }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Value of the primitive variable vector.
   */
  inline su2double* GetPrimitive(void) override {
    return Primitive.data();
  }

  /*!
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetPrimitive(unsigned short val_var) override {
    return Primitive.at(val_var);
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_prim - Value of the selected primitive variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitive(unsigned short val_var, su2double val_prim) override {
    Primitive.at(val_var) = val_prim;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitive(su2double* val_prim) override {
    SU2_Assert(val_prim != NULL, "The array of primitive variables has not been allocated");
    std::copy(val_prim,val_prim + nPrimVar,Primitive.begin());
  }

  /*!
   * \brief Get the value of the limiter for gradient of variables gradient.
   * \return Value of the limiter for the gradient of primitive variables.
   */
  inline su2double* GetLimiter_Primitive(void) override {
    return Limiter_Primitive.data();
  }

  /*!
   * \brief Get the value of the primitive variables limiter.
   * \param[in] val_var - Index of the variable.
   * \return Value of the selected primitive variable limiter.
   */
  inline su2double GetLimiter_Primitive(unsigned short val_var) override {
    return Limiter_Primitive.at(val_var);
  }

  /*!
   * \brief Set the value of the limiter.
   */
  inline void SetLimiter_Primitive(unsigned short val_var, su2double val_value) override {
    Limiter_Primitive.at(val_var) = val_value;
  }

  /*!
	 * \brief Set to zero the gradient of the primitive variables.
	 */
	inline void SetGradient_PrimitiveZero(unsigned short val_primvar) override {
    SU2_Assert(Gradient_Primitive != NULL, "The matrix for gradient primitive has not been allocated");
    for(unsigned short iVar = 0; iVar < val_primvar; ++iVar) {
      SU2_Assert(Gradient_Primitive[iVar] != NULL, std::string("The row " + std::to_string(iVar) + " of gradient primitive has not been allocated"));
      std::fill(Gradient_Primitive[iVar], Gradient_Primitive[iVar] + nDim, 0.0);
    }
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double** GetGradient_Primitive(void) override {
    return Gradient_Primitive;
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
   */
  su2double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    return Gradient_Primitive[val_var][val_dim];
  }

  /*!
	 * \brief Add val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the selected primitive variable.
	 */
	inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    Gradient_Primitive[val_var][val_dim] += val_value;
  }

  /*!
	 * \brief Subtract val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the selected primitive variable.
	 */
	inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    Gradient_Primitive[val_var][val_dim] -= val_value;
  }

  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    Gradient_Primitive[val_var][val_dim] = val_value;
  }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   * \param[in] config - Configuration of the particular problem.
   */
  bool SetPrimVar(CConfig* config) override;

  /*!
   * \brief Set all the primitive variables form conserved variables.
   * \param[in] config - Configuration of the particular problem.
   * \param[in] U - Storage of conservative variables.
   * \param[in] V - Storage of primitive variables.
   */
  bool Cons2PrimVar(CConfig* config, su2double* U, su2double* V);

  /*!
   * \brief Set all the conserved variables from primitive variables.
   * \param[in] config - Configuration of the particular problem.
   * \param[in] U - Storage of conservative variables.
   * \param[in] V - Storage of primitive variables.
   */
  void Prim2ConsVar(CConfig* config, su2double* V, su2double* U) override;

  /*!
	 * \brief Set the value of the mixture density.
	 */
	bool SetDensity(void) override;

	/*!
	 * \brief Set the value of the pressure.  Requires T calculation.
   * \param[in] config - Configuration of the particular problem.
	 */
	bool SetPressure(CConfig* config) override;

	/*!
	 * \brief Set the value of the speed of the sound.
	 */
	bool SetSoundSpeed(CConfig* config) override;

	/*!
	 * \brief Set the value of the total enthalpy (need a call of SetPressure()).
	 */
	inline void SetEnthalpy(void) override {
    SU2_Assert(Solution != NULL,"The array of solution variables has not been allocated");
    Primitive.at(H_INDEX_PRIM) = (Solution[RHOE_INDEX_SOL] + Primitive.at(P_INDEX_PRIM)) / Solution[RHO_INDEX_SOL];
  }

  /*!
   * \brief Calculates partial derivative of pressure w.r.t. conserved variables \f$\frac{\partial P}{\partial U}\f$
   * \param[in] config - Configuration settings
   * \param[in] dPdU - Passed-by-reference array to assign the derivatives
   */
  void CalcdPdU(su2double* V, CConfig* config, su2double* dPdU) override;

  /*!
   * \brief Calculates partial derivative of temperature w.r.t. conserved variables \f$\frac{\partial P}{\partial U}\f$
   * \param[in] config - Configuration settings
   * \param[in] dTdU - Passed-by-reference array to assign the derivatives
   */
  void CalcdTdU(su2double* V, CConfig* config, su2double* dTdU) override;

  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  inline su2double* GetdPdU(void) override {
    return dPdU.data();
  }

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
   */
  inline su2double* GetdTdU(void) override {
    return dTdU.data();
  }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(void) override {
    return Primitive.at(P_INDEX_PRIM);
  }

  /*!
   * \brief Get the speed of the sound.
   * \return Value of speed of the sound.
   */
  inline su2double GetSoundSpeed(void) override {
    return Primitive.at(A_INDEX_PRIM);
  }

  /*!
   * \brief Get the total enthalpy of the flow.
   * \return Value of the enthalpy of the flow.
   */
  inline su2double GetEnthalpy(void) override {
    return Primitive.at(H_INDEX_PRIM);
  }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(void) override {
    return Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Get the mass fraction \f$\rho_s / \rho \f$ of species s.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline su2double GetMassFraction(unsigned short val_Species) override {
    return Primitive.at(RHOS_INDEX_PRIM + val_Species);
  }

  /*!
   * \brief Get the mass fraction \f$\rho_s / \rho \f$ of species s.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline RealVec GetMassFractions(void) const {
    return RealVec(Primitive.cbegin() + RHOS_INDEX_PRIM, Primitive.cbegin() + (RHOS_INDEX_PRIM + nSpecies));
  }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(void) override {
    SU2_Assert(Solution != NULL,"The array with the solution has not been allocated");
    return Solution[RHOE_INDEX_SOL]/Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(void) override {
    return Primitive.at(T_INDEX_PRIM);
  }

  /*!
   * \brief Set the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline bool SetTemperature(su2double val_T) override {
    Primitive.at(T_INDEX_PRIM) = val_T;
    if(val_T < EPS)
      return true;

    return false;
  }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline su2double GetVelocity(unsigned short val_dim) override {
    return Primitive.at(VX_INDEX_PRIM + val_dim);
  }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(su2double* val_vector) override {
    SU2_Assert(val_vector != NULL, "The vector where to project velocity has not been allocated");
    return std::inner_product(Primitive.cbegin() + VX_INDEX_PRIM, Primitive.cbegin() + (VX_INDEX_PRIM + nDim), val_vector, 0.0);
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  void SetVelocity_Old(su2double* val_velocity) override;
};
const unsigned short CReactiveEulerVariable::P_INDEX_PRIM = CReactiveEulerVariable::VX_INDEX_PRIM + CReactiveEulerVariable::nDim;
const unsigned short CReactiveEulerVariable::RHO_INDEX_PRIM = CReactiveEulerVariable::P_INDEX_PRIM + 1;
const unsigned short CReactiveEulerVariable::H_INDEX_PRIM = CReactiveEulerVariable::RHO_INDEX_PRIM + 1;
const unsigned short CReactiveEulerVariable::A_INDEX_PRIM = CReactiveEulerVariable::H_INDEX_PRIM + 1;
const unsigned short CReactiveEulerVariable::RHOS_INDEX_PRIM = CReactiveEulerVariable::A_INDEX_PRIM + 1;

const unsigned short CReactiveEulerVariable::RHOE_INDEX_SOL = CReactiveEulerVariable::RHOVX_INDEX_SOL + CReactiveEulerVariable::nDim;
const unsigned short CReactiveEulerVariable::RHOS_INDEX_SOL = CReactiveEulerVariable::RHOE_INDEX_SOL + 1;

const unsigned short CReactiveEulerVariable::P_INDEX_GRAD = CReactiveEulerVariable::VX_INDEX_GRAD + CReactiveEulerVariable::nDim;

const unsigned short CReactiveEulerVariable::P_INDEX_LIM = CReactiveEulerVariable::VX_INDEX_LIM + CReactiveEulerVariable::nDim;

/*! \class CReactiveNSVariable
 *  \brief Main class for defining a variable for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSVariable:public CReactiveEulerVariable {
public:
  typedef Eigen::Matrix<su2double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RealMatrix;

protected:
	RealMatrix Diffusion_Coeffs;    /*!< \brief Binary diffusion coefficients of the mixture. */
  su2double  Laminar_Viscosity;	/*!< \brief Viscosity of the fluid. */
  su2double  Thermal_Conductivity;   /*!< \brief Thermal conductivity of the gas mixture. */

public:
  static const unsigned short RHOS_INDEX_GRAD;

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSVariable():CReactiveEulerVariable(),Laminar_Viscosity(),Thermal_Conductivity() {}

  /*!
   * \overloaded Class constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config);

  /*!
	 * \overload Class constructor
	 * \param[in] val_pressure - Value of the flow pressure (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the flow temperature (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,const su2double val_temperature,
                      unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad,unsigned short val_nprimvarlim, CConfig* config);

  /*!
	 * \overload Class constructor
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                      unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, CConfig* config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSVariable() {};

  /*!
   * \brief Set all primitive variables and transport properties for compressible flows.
   * \param[in] config - Configuration of the particular problem.
   */
  bool SetPrimVar(CConfig* config) override;

  /*!
   * \brief Get the laminar viscosity of the mixture.
   * \return Laminar viscoisty of the mixture
   */
  inline su2double GetLaminarViscosity(void) override {
    return Laminar_Viscosity;
  }

  /*!
   * \brief Get the thermal conductivity of the mixture.
   * \return Laminar viscoisty of the mixture
   */
  inline su2double GetThermalConductivity(void) override {
    return Thermal_Conductivity;
  }

  /*!
	 * \brief Get the species diffusion coefficient.
	 * \return Value of the species diffusion coefficient.
	 */
  inline RealMatrix GetBinaryDiffusionCoeff(void) {
    return Diffusion_Coeffs;
  }

  /*!
   * \brief Set the laminar viscosity of the mixture.
   * \param[in] laminarViscosity - value of laminar viscosity to set
   * \return Laminar viscoisty of the mixture
   */
  inline void SetLaminarViscosity(su2double laminarViscosity) override {
    Laminar_Viscosity = laminarViscosity;
  }

  /*!
   * \brief Set the thermal conductivity of the mixture.
   * \param[in] thermalConductivity - value of thermal conductivity to set
   * \return Laminar viscoisty of the mixture
   */
  inline void SetThermalConductivity(su2double thermalConductivity) override {
    Thermal_Conductivity = thermalConductivity;
  }

  /*!
   * \brief Set the temperature at the wall
   * \param[in] temperature_wall - value of the wall temperature to set
   */
	inline void SetWallTemperature(su2double temperature_wall) override {
    Primitive.at(T_INDEX_PRIM) = temperature_wall;
  }
};
const unsigned short CReactiveNSVariable::RHOS_INDEX_GRAD = CReactiveNSVariable::P_INDEX_GRAD + 1;

#endif
