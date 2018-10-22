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

  typedef std::vector<su2double> RealVec;
  typedef std::vector<RealVec>  RealMatrix;
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

  RealVec    Ri; /*!< \brief constant gas for each species in the mixture. */
  SmartArr   xi,  /*!< \brief Rotation modes for each species in the mixture. */
             Ms, /*!< \brief Molar masses for each species in the mixture. */
             hf, /*!< \brief Formation enthalpies for each species in the mixture. */
             Tref; /*!< \brief Reference temperatures for each species in the mixture. */

public:

  /**
   * Enumerator defining the mapping between the primitive variable name
   * and its position in the physical data
   */
  enum {T_INDEX_PRIM=0, VX_INDEX_PRIM=1, VY_INDEX_PRIM=2, VZ_INDEX_PRIM=3,P_INDEX_PRIM=VZ_INDEX_PRIM+1, RHO_INDEX_PRIM=P_INDEX_PRIM+1,
        H_INDEX_PRIM=RHO_INDEX_PRIM+1, A_INDEX_PRIM=H_INDEX_PRIM+1, RHOS_INDEX_PRIM=A_INDEX_PRIM+1};

  /**
   * Enumerator defining the mapping between the solution variable name
   * and its position in the physical data
   */
  enum {RHO_INDEX_SOL=0, RHOVX_INDEX_SOL=1, RHOVY_INDEX_SOL=2,RHOVZ_INDEX_SOL=3,
        RHOE_INDEX_SOL=RHOVZ_INDEX_SOL+1, RHOS_INDEX_SOL=RHOE_INDEX_SOL+1};


  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerVariable();

  /*!
   * \overloaded Constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, std::shared_ptr<CConfig> config);

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
    std::vector<su2double> tmp(val_prim,val_prim + nPrimVar);
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
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  su2double** GetGradient_Primitive(void) override;

  /*!
	 * \brief Add val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to add to the gradient of the selected primitive variable.
	 */
	inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    Gradient_Primitive.at(val_var).at(val_dim) += val_value;
  }

  /*!
	 * \brief Subtract val_value to the gradient of the val_var primitive variable.
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_dim - Index of the dimension.
	 * \param[in] val_value - Value to subtract to the gradient of the selected primitive variable.
	 */
	inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) override {
    Gradient_Primitive.at(val_var).at(val_dim) -= val_value;
  }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(CConfig* config) override;
  //bool SetPrimVar_Compressible(CConfig *config) ovveride;

  /*!
   * \brief Set gradient primitive variables for compressible flows.
   */
  //void SetPrimVar_Gradient(CConfig* config);

  /*!
   * \brief Set all the primitive variables form conserved variables.
   */
  bool Cons2PrimVar(CConfig* config, su2double* U, su2double* V, su2double* dPdU, su2double* dTdU);

  /*!
   * \brief Set Gradient of the primitive variables from conserved variables
   */
  //bool GradCons2GradPrimVar(std::unique_ptr<CConfig> config, RealVec& U, RealVec& V,
  //                          RealMatrix& GradU, RealMatrix& GradV);

  /*!
   * \brief Set all the conserved variables from primitive variables.
   */
  void Prim2ConsVar(CConfig* config, su2double* V, su2double* U) override;

  /*!
	 * \brief Set the value of the mixture density.
	 */
	bool SetDensity(void) override;

	/*!
	 * \brief Set the value of the pressure.  Requires T calculation.
	 */
	void SetPressure(void) override;

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
   * \brief Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  //inline su2double GetRhoCv_tr(void) override {
  //  return Primitive.at(RHOCV_INDEX_PRIM);
  //}

  /*!
   * \brief Set the squared velocity
   */
  //void SetVelocity2(void) override;

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
  inline void SetVelocity_Old(su2double* val_velocity) override {
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      Solution_Old[RHOVX_INDEX_SOL+iDim] = val_velocity[iDim]*Primitive.at(RHO_INDEX_PRIM);
  }

  /*!
   * \brief Set the temperature.
   * \param[in] config - Configuration parameters.
   */
  bool SetTemperature(CConfig *config) override;

};






















































































/*! \class CReactiveEulerVariable
 *  \brief Main class for defining a variable for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveNSVariable:public CReactiveEulerVariable {
protected:
	su2double  Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
	su2double  Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
	su2double  Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
  RealVec    Diffusion_Coeffs;    /*!< \brief Diffusion coefficients of the mixture. */
  RealMatrix Dij;             /*!< \brief Binary diffusion coefficients. */
	su2double  Laminar_Viscosity;	/*!< \brief Viscosity of the fluid. */
  su2double  Thermal_Conductivity;       /*!< \brief Thermal conductivity of the gas mixture. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSVariable() : CReactiveEulerVariable(),Temperature_Ref(),Viscosity_Ref(),Viscosity_Inf(),
                          Laminar_Viscosity(),Thermal_Conductivity() {}

  /*!
   * \overloaded Constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, std::unique_ptr<CConfig>& config);

  /*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSVariable() {};

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
  	 * \brief Set the temperature at the wall
  	 */
	inline void SetWallTemperature(su2double temperature_wall) override {
    Primitive.at(T_INDEX_PRIM) = temperature_wall;
  }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  //bool SetPrimVar_Compressible(CConfig *config) ovveride;

  /*!
   * \brief Set gradient primitive variables for compressible flows.
   */
  //void SetPrimVar_Gradient(CConfig *config);

};

#endif
