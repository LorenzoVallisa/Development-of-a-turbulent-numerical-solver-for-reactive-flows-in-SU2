#ifndef SU2_VARIABLE_REACTIVE
#define SU2_VARIABLE_REACTIVE

#include <memory>

#include "variable_structure.hpp"
#include "../../Common/include/physical_property_library.hpp"
#include "../../Common/include/datatypes/vectorT.hpp"
#include "../../Common/include/datatypes/matrixT.hpp"

/*! \class CReactiveEulerVariable
 *  \brief Main class for defining a variable for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */

class CReactiveEulerVariable:public CVariable {
public:

  using RealVec = Common::RealVec;
  using RealMatrix = Common::RealMatrix;
  typedef std::shared_ptr<Framework::PhysicalPropertyLibrary> LibraryPtr;
  typedef std::unique_ptr<su2double[]> SmartArr;

protected:

  LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  const unsigned short nSpecies; /*!< \brief Number of species in the mixture. */

  /*--- Primitive variable definition ---*/
  RealVec    Primitive; /*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,a,rho1,...rhoNs) in compressible flows. */
  RealMatrix Gradient_Primitive; /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */
  RealVec    dPdU;                 /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  RealVec    dTdU;                /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */

public:

  /**
   * Mapping between the primitive variable name and its position in the physical data
   */
  static constexpr unsigned T_INDEX_PRIM = 0;
  static constexpr unsigned VX_INDEX_PRIM = 1;
  static const unsigned P_INDEX_PRIM;
  static const unsigned RHO_INDEX_PRIM;
  static const unsigned H_INDEX_PRIM;
  static const unsigned A_INDEX_PRIM;
  static const unsigned RHOS_INDEX_PRIM;

  //enum {T_INDEX_PRIM=0, VX_INDEX_PRIM=1, VY_INDEX_PRIM=2, VZ_INDEX_PRIM=3,P_INDEX_PRIM=VZ_INDEX_PRIM+1, RHO_INDEX_PRIM=P_INDEX_PRIM+1,
  //      H_INDEX_PRIM=RHO_INDEX_PRIM+1, A_INDEX_PRIM=H_INDEX_PRIM+1, RHOS_INDEX_PRIM=A_INDEX_PRIM+1};

  /**
   * Mapping between the solution variable name and its position in the physical data
   */
  static constexpr unsigned RHO_INDEX_SOL = 0;
  static constexpr unsigned RHOVX_INDEX_SOL = 1;
  static const unsigned RHOE_INDEX_SOL;
  static const unsigned RHOS_INDEX_SOL;
  //enum {RHO_INDEX_SOL=0, RHOVX_INDEX_SOL=1, RHOVY_INDEX_SOL=2,RHOVZ_INDEX_SOL=3,
  //      RHOE_INDEX_SOL=RHOVZ_INDEX_SOL+1, RHOS_INDEX_SOL=RHOE_INDEX_SOL+1};

  /**
   * Mapping between the primitive variable gradient name and its position in the physical data
   */
  static constexpr unsigned T_INDEX_GRAD = 0;
  static constexpr unsigned VX_INDEX_GRAD = 1;
  static const unsigned P_INDEX_GRAD;
  //enum {T_INDEX_GRAD=0, VX_INDEX_GRAD=1, VY_INDEX_GRAD=2,VZ_INDEX_GRAD=3,P_INDEX_GRAD = 4};

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerVariable();

  /*!
   * \overload
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad,std::shared_ptr<CConfig> config);

  /*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the flow energy (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(su2double val_pressure, RealVec& val_massfrac, RealVec& val_velocity, su2double val_temperature,
                         unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad, std::shared_ptr<CConfig> config);

	/*!
	 * \overload
	 * \param[in] val_solution - Vector with the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerVariable(RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                         unsigned short val_nprimvargrad, std::shared_ptr<CConfig> config);


  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerVariable() {}

  /*!
   * \brief Get the number of species in the mixture.
   * \return Number of species in the problem.
   */
  inline const unsigned short GetnSpecies(void) const {
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
    RealVec tmp(val_prim,val_prim + nPrimVar);
    Primitive = tmp;
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
    for (unsigned short iVar = 0; iVar < val_primvar; ++iVar) {
		  for (unsigned short iDim = 0; iDim < nDim; ++iDim)
			   Gradient_Primitive(iVar,iDim) = 0.0;
    }
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  su2double** GetGradient_Primitive(void) override;

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  //inline RealMatrix GetGradient_Primitive(void) const {
  //  return Gradient_Primitive;
  //}

  /*!
	 * \brief Add val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the selected primitive variable.
	 */
	inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    Gradient_Primitive.at(val_var,val_dim) += val_value;
  }

  /*!
	 * \brief Subtract val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the selected primitive variable.
	 */
	inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    Gradient_Primitive.at(val_var,val_dim) -= val_value;
  }

  /*!
	 * \brief Set the gradient of the primitive variables.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value of the gradient.
	 */
	inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    Gradient_Primitive.at(val_var,val_dim) = val_value;
  }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   * \param[in] config - Configuration of the particular problem.
   */
  bool SetPrimVar(CConfig* config) override;

  /*!
   * \brief Set all the primitive variables form conserved variables.
   * \param[in] U - Storage of conservative variables.
   * \param[in] V - Storage of primitive variables.
   */
  bool Cons2PrimVar(su2double* U, su2double* V);

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
	bool SetSoundSpeed(void) override;

	/*!
	 * \brief Set the value of the enthalpy.
	 */
	inline void SetEnthalpy(void) override {
    Primitive.at(H_INDEX_PRIM) = (Solution[RHOE_INDEX_SOL] + Primitive.at(P_INDEX_PRIM)) / Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Calculates enthalpy per mass for input species (not including KE)
   * \param[in] val_T - Temperature value
   * \param[in] val_Species - Index of desired species
   */
  su2double CalcHs(su2double val_T, unsigned short iSpecies);

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
  };

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
   */
  inline su2double* GetdTdU(void) override {
    return dTdU.data();
  };

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
    return Primitive.at(RHOS_INDEX_PRIM + val_Species)/Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(void) override {
    return Solution[RHOE_INDEX_SOL]/Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(void) override {
    return Primitive.at(T_INDEX_PRIM);
  };

  /*!
   * \brief Sets the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline bool SetTemperature(su2double val_T) override {
    Primitive.at(T_INDEX_PRIM) = val_T;
    return false;
  }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline su2double GetVelocity(unsigned short val_dim) override {
    return Primitive.at(VX_INDEX_PRIM+val_dim);
  }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  su2double GetProjVel(su2double* val_vector) override;

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  void SetVelocity_Old(su2double* val_velocity) override;
  /*!
   * \brief Set the temperature.
   * \param[in] config - Configuration parameters.
   */
  bool SetTemperature(CConfig* config) override;

};
const unsigned CReactiveEulerVariable::P_INDEX_PRIM = CReactiveEulerVariable::VX_INDEX_PRIM + CReactiveEulerVariable::nDim;
const unsigned CReactiveEulerVariable::RHO_INDEX_PRIM = CReactiveEulerVariable::P_INDEX_PRIM + 1;
const unsigned CReactiveEulerVariable::H_INDEX_PRIM = CReactiveEulerVariable::RHO_INDEX_PRIM + 1;
const unsigned CReactiveEulerVariable::A_INDEX_PRIM = CReactiveEulerVariable::H_INDEX_PRIM + 1;
const unsigned CReactiveEulerVariable::RHOS_INDEX_PRIM = CReactiveEulerVariable::A_INDEX_PRIM + 1;

const unsigned CReactiveEulerVariable::RHOE_INDEX_SOL = CReactiveEulerVariable::RHOVX_INDEX_SOL + CReactiveEulerVariable::nDim;
const unsigned CReactiveEulerVariable::RHOS_INDEX_SOL = CReactiveEulerVariable::RHOE_INDEX_SOL + 1;

const unsigned CReactiveEulerVariable::P_INDEX_GRAD = CReactiveEulerVariable::VX_INDEX_GRAD + CReactiveEulerVariable::nDim;


/*! \class CReactiveNSVariable
 *  \brief Main class for defining a variable for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSVariable:public CReactiveEulerVariable {
protected:
	RealVec    Diffusion_Coeffs;    /*!< \brief Diffusion coefficients of the mixture. */
  su2double  Laminar_Viscosity;	/*!< \brief Viscosity of the fluid. */
  su2double  Thermal_Conductivity;   /*!< \brief Thermal conductivity of the gas mixture. */

  unsigned short nPrimVarAvgGrad;
  RealMatrix AvgGradient_Primitive; /*!< \brief Gradient of the primitive for average computation (T, vx, vy, vz, P, rho,rho1...rhoNs). */

public:

  /**
   * Enumerator defining the mapping between the primitive variable gradient name
   * for average computations and its position in the physical data
   */
  static constexpr unsigned T_INDEX_AVGGRAD = 0;
  static constexpr unsigned VX_INDEX_AVGGRAD = 1;
  static const unsigned RHO_INDEX_AVGGRAD;
  static const unsigned RHOS_INDEX_AVGGRAD;
  //enum {T_INDEX_AVGGRAD=0,VX_INDEX_AVGGRAD=1,VY_INDEX_AVGGRAD=2,VZ_INDEX_AVGGRAD=3,RHO_INDEX_AVGGRAD=4,RHOS_INDEX_AVGGRAD=5};

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSVariable(): CReactiveEulerVariable(),Laminar_Viscosity(),Thermal_Conductivity(),nPrimVarAvgGrad() {}

  /*!
   * \overloaded Constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad, unsigned short val_nprimvar_avggrad, std::shared_ptr<CConfig> config);

  /*!
	 * \overload
	 * \param[in] val_density - Value of the flow density (initialization value).
	 * \param[in] val_velocity - Value of the flow velocity (initialization value).
	 * \param[in] val_temperature - Value of the flow temperature (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nprimvar - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(su2double val_density, RealVec& val_massfrac, RealVec& val_velocity,su2double val_temperature,
                      unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad, unsigned short val_nprimvar_avggrad, std::shared_ptr<CConfig> config);

  /*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_nDim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of conserved variables.
   * \param[in] val_nprimvar - Number of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] val_nprimvargrad - Number of gradient of primitive variables of the problem.
   * \param[in] config - Definition of the particular problem.
	 */
	CReactiveNSVariable(RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nprimvar,
                      unsigned short val_nprimvargrad, unsigned short val_nprimvar_avggrad, std::shared_ptr<CConfig> config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSVariable() {};

  /*!
	 *\ Set value of gradient for average computations
	 * \param[in] iVar - Index of the desired variable
	 * \param[in] iDim - Index of dimension.
	 * \param[in] value - Value to set.
   */
  inline void SetAvgGradient(unsigned short iVar,unsigned short iDim,su2double value) {
    AvgGradient_Primitive.at(iVar,iDim) = value;
  }

  /*!
	 *\ Set the whole matrix for average gradient computations to zero
	 */
  inline void SetAvgGradient_Zero(void) {
    std::fill(AvgGradient_Primitive.begin(),AvgGradient_Primitive.end(),0.0);
  }

  /*!
	 *\ Get value of gradient for average computations
	 * \param[in] iVar - Index of the desired variable
	 * \param[in] iDim - Index of dimension.
   */
  inline su2double GetAvgGradient(unsigned short iVar,unsigned short iDim) {
    return AvgGradient_Primitive.at(iVar,iDim);
  }

  /*!
	 *\ Get value of gradient for average computations
	 */
  inline RealMatrix GetAvgGradient() {
    return AvgGradient_Primitive;
  }


  /*!
	 *\ Add value of gradient for average computations
	 * \param[in] iVar - Index of the desired variable
	 * \param[in] iDim - Index of dimension.
	 * \param[in] value - Value to add.
   */
  inline void AddAvgGradient(unsigned short iVar,unsigned short iDim,su2double value) {
    AvgGradient_Primitive.at(iVar,iDim) += value;
  }

  /*!
	 *\ Subtract value of gradient for average computations
	 * \param[in] iVar - Index of the desired variable
	 * \param[in] iDim - Index of dimension.
	 * \param[in] value - Value to add.
   */
  inline void SubtractAvgGradient(unsigned short iVar,unsigned short iDim,su2double value) {
    AvgGradient_Primitive.at(iVar,iDim) -= value;
  }

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
const unsigned CReactiveNSVariable::RHO_INDEX_AVGGRAD = CReactiveNSVariable::VX_INDEX_AVGGRAD + CReactiveNSVariable::nDim;
const unsigned CReactiveNSVariable::RHOS_INDEX_AVGGRAD = CReactiveNSVariable::VX_INDEX_AVGGRAD + 1;

#endif
