#ifndef SU2_VARIABLE_REACTIVE
#define SU2_VARIABLE_REACTIVE

#include "variable_structure.hpp"
#include "../../Common/include/Framework/factory.hpp"
#include "../../Common/include/Framework/su2_assert.hpp"

#include "../../externals/Eigen/Dense"
#include <numeric>
#include <memory>

/*! \class CReactiveEulerVariable
 *  \brief Main class for defining a variable for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerVariable: public CVariable {
public:
  typedef std::vector<su2double> RealVec;
  typedef su2double** SU2Matrix;
  typedef std::shared_ptr<Framework::PhysicalChemicalLibrary> LibraryPtr;

protected:
  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  unsigned short nSpecies; /*!< \brief Number of species in the mixture. */
  unsigned short nPrimVarLim; /*!< \brief Number of primitive variables to limit in the problem. */

  bool US_System;             /*!< \brief Flag for US units. */

  su2double Cp;                 /*!< \brief Specific heat at constant pressure. */

  /*--- Primitive variable definition ---*/
  RealVec    Primitive; /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, a, Y1,...YNs) in compressible flows. */
  SU2Matrix  Gradient_Primitive; /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    dPdU;                 /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  RealVec    dTdU;                /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */

  RealVec Ys;               /*!< \brief Auxiliary vector to store mass fractions separately. */
  RealVec dTdYs,            /*!< \brief Auxiliary vector for temperature derivatives w.r.t partial densities. */
          dPdYs;            /*!< \brief Auxiliary vector for pressure derivatives w.r.t partial densities. */

  /**
   * Mapping between the primitive variable name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_PRIM = 0;
  static constexpr unsigned short VX_INDEX_PRIM = 1;
  static unsigned short P_INDEX_PRIM;
  static unsigned short RHO_INDEX_PRIM;
  static unsigned short H_INDEX_PRIM;
  static unsigned short A_INDEX_PRIM;
  static unsigned short RHOS_INDEX_PRIM;

  /**
   * Mapping between the solution variable name and its position in the physical data
   */
  static constexpr unsigned short RHO_INDEX_SOL = 0;
  static constexpr unsigned short RHOVX_INDEX_SOL = 1;
  static unsigned short RHOE_INDEX_SOL;
  static unsigned short RHOS_INDEX_SOL;

  /**
   * Mapping between the primitive variable gradient name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_GRAD = 0;
  static constexpr unsigned short VX_INDEX_GRAD = 1;
  static unsigned short P_INDEX_GRAD;

  /**
   * Mapping between the primitivelimited variable name and its position in the physical data
   */
  static constexpr unsigned short T_INDEX_LIM = 0;
  static constexpr unsigned short VX_INDEX_LIM = 1;
  static unsigned short P_INDEX_LIM;

public:

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerVariable();

  /*!
   * \overload Class constructor to initialize dimensions of the problem.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \overload Class constructor with pressure, temperature, mass fractions and velocity.
	 * \param[in] val_pressure - Value of the flow pressure (initialization value).
   * \param[in] val_massfrac - Value of mass fractions (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the temperature (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                         const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

	/*!
	 * \overload Class constructor with a complete initial state.
	 * \param[in] val_solution - Vector with the flow values (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \overload Class constructor with a complete initial state.
	 * \param[in] val_solution - Array with the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(su2double* val_solution, unsigned short val_nDim, unsigned short val_nvar,
                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                         unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

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
   * \brief Get the number of species in the mixture.
   * \return Number of species in the problem.
   */
  inline unsigned short GetnSpecies(void) {
    return nSpecies;
  }

  /*!
   * \brief Get the index of temperature in primitive variables array
   * \return Index of temperature in primitive variables array
   */
  inline static const unsigned short GetT_INDEX_PRIM(void) {
    return T_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of velocity along x in primitive variables array
   * \return Index of x-velocity component in primitive variables array
   */
  inline static const unsigned short GetVX_INDEX_PRIM(void) {
    return VX_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of pressure in primitive variables array
   * \return Index of pressure in primitive variables array
   */
  inline static const unsigned short GetP_INDEX_PRIM(void) {
    return P_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of density in primitive variables array
   * \return Index of density in primitive variables array
   */
  inline static const unsigned short GetRHO_INDEX_PRIM(void) {
    return RHO_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of total enthalpy in primitive variables array
   * \return Index of total enthalpy in primitive variables array
   */
  inline static const unsigned short GetH_INDEX_PRIM(void) {
    return H_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of speed of sound in primitive variables array
   * \return Index of speed of sound in primitive variables array
   */
  inline static const unsigned short GetA_INDEX_PRIM(void) {
    return A_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of mass fractions in primitive variables array
   * \return Index of mass fractions in primitive variables array
   */
  inline static const unsigned short GetRHOS_INDEX_PRIM(void) {
    return RHOS_INDEX_PRIM;
  }

  /*!
   * \brief Get the index of density in conserved variables array
   * \return Index of density in conserved variables array
   */
  inline static const unsigned short GetRHO_INDEX_SOL(void) {
    return RHO_INDEX_SOL;
  }

  /*!
   * \brief Get the index of momentum along x in conserved variables array
   * \return Index of x-momentum component in conserved variables array
   */
  inline static const unsigned short GetRHOVX_INDEX_SOL(void) {
    return RHOVX_INDEX_SOL;
  }

  /*!
   * \brief Get the index of density times total energy in conserved variables array
   * \return Index of density times totla energy in conserved variables array
   */
  inline static const unsigned short GetRHOE_INDEX_SOL(void) {
    return RHOE_INDEX_SOL;
  }

  /*!
   * \brief Get the index of partial densities in conserved variables array
   * \return Index of partial densities in conserved variables array
   */
  inline static const unsigned short GetRHOS_INDEX_SOL(void) {
    return RHOS_INDEX_SOL;
  }

  /*!
   * \brief Get the index of temperature in primitive gradient
   * \return Index of temperature in primitive gradient
   */
  inline static const unsigned short GetT_INDEX_GRAD(void) {
    return T_INDEX_GRAD;
  }

  /*!
   * \brief Get the index of velocity along x in primitive gradient
   * \return Index of velocity along x in primitive gradient
   */
  inline static const unsigned short GetVX_INDEX_GRAD(void) {
    return VX_INDEX_GRAD;
  }

  /*!
   * \brief Get the index of pressure in primitive gradient
   * \return Index of pressure in primitive gradient
   */
  inline static const unsigned short GetP_INDEX_GRAD(void) {
    return P_INDEX_GRAD;
  }

  /*!
   * \brief Get the index of temperature for limited variables
   * \return Index of temperature for limited variables
   */
  inline static const unsigned short GetT_INDEX_LIM(void) {
    return T_INDEX_LIM;
  }

  /*!
   * \brief Get the index of velocity along x for limited variables
   * \return Index of velocity along x for limited variables
   */
  inline static const unsigned short GetVX_INDEX_LIM(void) {
    return VX_INDEX_LIM;
  }

  /*!
   * \brief Get the index of pressure for limited variables
   * \return Index of pressure for limited variables
   */
  inline static const unsigned short GetP_INDEX_LIM(void) {
    return P_INDEX_LIM;
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
    std::copy(val_prim, val_prim + nPrimVar, Primitive.begin());
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
      SU2_Assert(Gradient_Primitive[iVar] != NULL,
                 std::string("The row " + std::to_string(iVar) + " of gradient primitive has not been allocated"));
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
    SU2_Assert(Gradient_Primitive[val_var] != NULL,
               std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    return Gradient_Primitive[val_var][val_dim];
  }

  /*!
	 * \brief Add val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the selected primitive variable.
	 */
	inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,
               std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    Gradient_Primitive[val_var][val_dim] += val_value;
  }

  /*!
	 * \brief Subtract val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the selected primitive variable.
	 */
	inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,
               std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
    Gradient_Primitive[val_var][val_dim] -= val_value;
  }

  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    SU2_Assert(Gradient_Primitive[val_var] != NULL,
               std::string("The row " + std::to_string(val_var) + " of gradient primitive hasn't been allocated"));
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
    Primitive.at(H_INDEX_PRIM) = (Solution[RHOE_INDEX_SOL] + Primitive.at(P_INDEX_PRIM))/Solution[RHO_INDEX_SOL];
  }

  /*!
   * \brief Compute partial derivative of pressure w.r.t. conserved variables \f$\frac{\partial P}{\partial U}\f$
   * \param[in] V - Actual state
   * \param[in] config - Configuration settings
   * \param[in] dPdU - Array to assign the derivatives
   */
  void CalcdPdU(su2double* V, CConfig* config, su2double* dPdU) override;

  /*!
   * \brief Compute partial derivative of temperature w.r.t. conserved variables \f$\frac{\partial T}{\partial U}\f$
   * \param[in] V - Actual state
   * \param[in] config - Configuration settings
   * \param[in] dTdU - Srray to assign the derivatives
   */
  void CalcdTdU(su2double* V, CConfig* config, su2double* dTdU) override;

  /*!
   * \brief Get partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  inline su2double* GetdPdU(void) override {
    return dPdU.data();
  }

  /*!
   * \brief Get partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
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
   * \brief Get the total energy of the flow.
   * \return Value of the total energy of the flow.
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
   * \brief Get the projected velocity in a unitary vector direction.
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(su2double* val_vector) override {
    SU2_Assert(val_vector != NULL, "The vector where to project velocity has not been allocated");
    return std::inner_product(Primitive.cbegin() + VX_INDEX_PRIM, Primitive.cbegin() + (VX_INDEX_PRIM + nDim), val_vector, 0.0);
  }

  /*!
   * \brief Get the squared velocity.
   * \return Value of the squared velocity.
   */
  inline su2double GetVelocity2(void) override {
    return std::inner_product(Primitive.cbegin() + VX_INDEX_PRIM, Primitive.cbegin() + (VX_INDEX_PRIM + nDim),
                              Primitive.cbegin() + VX_INDEX_PRIM, 0.0);
  }

  /*!
   * \brief Get the specific heat at constant pressure.
   * \return Value of the specific heat at constant pressure.
   */
  inline su2double GetSpecificHeatCp(void) override {
    return Cp;
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  void SetVelocity_Old(su2double* val_velocity) override;
};

/*! \class CReactiveNSVariable
 *  \brief Main class for defining a variable for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSVariable: public CReactiveEulerVariable {
public:
  using RealMatrix = Eigen::MatrixXd;

protected:
  su2double  Laminar_Viscosity;	      /*!< \brief Laminar viscosity of the fluid. */
  su2double  Thermal_Conductivity;   /*!< \brief Thermal conductivity of the gas mixture. */
  RealMatrix Diffusion_Coeffs;      /*!< \brief Binary diffusion coefficients of the mixture. */

public:
  static unsigned short RHOS_INDEX_GRAD; /*!< \brief Index for position of mole fractions in primitives gradient. */

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSVariable(): CReactiveEulerVariable(), Laminar_Viscosity(), Thermal_Conductivity() {}

  /*!
   * \overloaded Class constructor to initialize dimension of the problem.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad, unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \overload Class constructor with pressure, temperature, mass fractions and velocity.
	 * \param[in] val_pressure - Value of the flow pressure (initialization value).
   * \param[in] val_massfrac - Value of the mass fractions (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_viscosity - Value of the flow viscosity (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                      const su2double val_temperature, const su2double val_viscosity, unsigned short val_nDim, unsigned short val_nvar,
                      unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                      unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \overload Class constructor with a complete initial state of the flow.
	 * \param[in] val_solution - Vector with the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                      unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                      unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \overload Class constructor with a complete initial state.
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nSpecies - Number of species in the mixture
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvarlim - Number of primitive variables to limit in the problem.
   * \param[in] lib_ptr - Pointer to the external library for physical-chemical properties
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(su2double* val_solution, unsigned short val_nDim, unsigned short val_nvar,
                      unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                      unsigned short val_nprimvarlim, LibraryPtr lib_ptr, CConfig* config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSVariable() {}

  /*!
   * \brief Get the index of mole fractions in primitive gradient
   * \return Index of mole fractions in primitive gradient
   */
  inline static const unsigned short GetRHOS_INDEX_GRAD(void) {
    return RHOS_INDEX_GRAD;
  }

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
  inline su2double* GetDiffusionCoeff(void) override {
    return Diffusion_Coeffs.data();
  }

  /*!
   * \brief Set the laminar viscosity of the mixture.
   * \param[in] laminarViscosity - value of laminar viscosity to set
   */
  inline void SetLaminarViscosity(su2double laminarViscosity) override {
    Laminar_Viscosity = laminarViscosity;
  }

  /*!
   * \brief Set the thermal conductivity of the mixture.
   * \param[in] thermalConductivity - value of thermal conductivity to set
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

#endif
